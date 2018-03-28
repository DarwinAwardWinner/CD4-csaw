#!/usr/bin/env Rscript

library(getopt)
library(optparse)
library(assertthat)
library(rex)
library(sitools)

## Deparse and then concatenate into a single string
deparse_onestring <- function(...) {
    deparse(...) %>% str_c(collapse="\n")
}

## Extension of match.arg with automatic detection of the argument
## name for use in error messages.
match.arg <- function (arg, choices, several.ok = FALSE, argname=substitute(arg), ignore.case=FALSE) {
    if (missing(choices)) {
        formal.args <- formals(sys.function(sys.parent()))
        choices <- eval(formal.args[[as.character(substitute(arg))]])
    }
    if (is.null(arg))
        return(choices[1L])
    else if (!is.character(arg))
        stop(glue("{dQuote(argname)} must be NULL or a character vector"))
    if (!several.ok) {
        if (identical(arg, choices))
            return(arg[1L])
        if (length(arg) > 1L)
            stop(glue("{dQuote(argname)} must be of length 1"))
    }
    else if (length(arg) == 0L)
        stop(glue("{dQuote(argname)} must be of length >= 1"))
    fold_case <- identity
    if (ignore.case) {
        fold_case <- tolower
    }
    i <- pmatch(fold_case(arg), fold_case(choices), nomatch = 0L, duplicates.ok = TRUE)
    if (all(i == 0L))
        stop(gettextf("%s should be one of %s", dQuote(argname), paste(dQuote(choices),
            collapse = ", ")), domain = NA)
    i <- i[i > 0L]
    if (!several.ok && length(i) > 1)
        stop("there is more than one match in 'match.arg'")
    choices[i]
}

get.options <- function(opts) {

    ## Do argument parsing early so the script exits quickly if
    ## arguments are invalid
    optlist <- list(
        make_option(c("-q", "--transcript-quant"), metavar="SEXP.RDS", type="character",
                    help="File name of an R data file containing a RangedSummarizedExperiment of transcript abundances. These will be used to select the highest-expressed TSS for each gene."),
        ## So far this script only supports TxDb objects because
        ## figuring out the first exon and TSS from other
        ## less-structured formats is a pain.
        make_option(c("-t", "--annotation-txdb"), metavar="TXDBNAME", type="character",
                    help="Name of TxDb package, or the name of a database file, to use for gene annotation"),
        make_option(c("-o", "--output-file"), metavar="FILENAME.RDS", type="character",
                    help="Output file name. The GRanges object containing the promoter regions will be saved here using saveRDS, so it should end in '.RDS'."),
        make_option(c("-a", "--additional-gene-info"), metavar="FILENAME", type="character",
                    help="RDS/RData/xlsx/csv file containing a table of gene metadata. Row names (or the first column of the file if there are no row names) should be gene/feature IDs that match the ones used in the main annotation, and these should be unique. This can also be a GFF3 file where the metadata is in the attributes of elements of type 'gene', where the 'ID' attribute specifies the gene ID."))
    progname <- na.omit(c(get_Rscript_filename(), "rnaseq-count.R"))[1]
    parser <- OptionParser(
        usage="Usage: %prog [options] -q SEXP.RDS -t TXDB -o OUTPUT.RDS",
        description="Select the most abundant TSS for each gene.

For each gene, transcripts are grouped by TSS, and their average abundances are added up. The TSS with the largest sum of average transcript abundances is selected as the representative TSS for that gene. These are all stored in a GRanges object in the output file. The resulting GRanges object will be annotated with a GeneID column. For transcripts with no associaated Gene ID, the GeneID column will be identical to the TxID column. Since a single TSS is being chosen for each gene, the GeneID column should not contain any duplicates.",
option_list = optlist,
add_help_option = TRUE,
prog = progname)

    cmdopts <- parse_args(parser, opts)
    ## Ensure that all required arguments were provided
    required.opts <- c("annotation-txdb", "output-file", "transcript-quant")
    missing.opts <- setdiff(required.opts, names(cmdopts))
    if (length(missing.opts) > 0) {
        stop(str_c("Missing required arguments: ", deparse_onestring(missing.opts)))
    }

    ## Replace dashes with underscores so that all options can easily
    ## be accessed by "$"
    cmdopts %>% setNames(chartr("-", "_", names(.)))
}

## Do this early to handle "--help" before wasting time loading
## pacakges & stuff
invisible(get.options(commandArgs(TRUE)))

library(assertthat)
library(dplyr)
library(magrittr)
library(stringr)
library(glue)
library(future)
library(GenomicRanges)
library(edgeR)
library(SummarizedExperiment)

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

## Read a single R object from an RDA file. If run on an RDA
## file containing more than one object, throws an error.
read.single.object.from.rda <- function(filename) {
    objects <- within(list(), suppressWarnings(load(filename)))
    if (length(objects) != 1) {
        stop("RDA file should contain exactly one object")
    }
    return(objects[[1]])
}

## Read a single object from RDS or RDA file
read.RDS.or.RDA <- function(filename, expected.class="ANY") {
    object <- suppressWarnings(tryCatch({
        readRDS(filename)
    }, error=function(...) {
        read.single.object.from.rda(filename)
    }))
    if (!is(object, expected.class)) {
        object <- as(object, expected.class)
    }
    return(object)
}

save.RDS.or.RDA <-
    function(object, file, ascii = FALSE, version = NULL, compress = TRUE,
             savetype=ifelse(str_detect(file, regex("\\.rda(ta)?", ignore_case = TRUE)),
                             "rda", "rds")) {
    if (savetype == "rda") {
        save(list="object", file=file, ascii=ascii, version=version, compress=compress)
    } else{
        saveRDS(object=object, file=file, ascii=ascii, version=version, compress=compress)
    }
}

is.empty <- function(x) {
    x %>% unlist %>% na.omit %>% length %>% equals(0)
}

## Get column names that are always the same for all elements of a
## gene. Used for extracting only the gene metadata from exon
## metadata.
get.gene.common.colnames <- function(df, geneids, blacklist=c("type", "Parent")) {
    if (nrow(df) < 1) {
        return(character(0))
    }
    if (any(is.na(geneids))) {
        stop("Gene IDs cannot be undefined")
    }
    if (any(lengths(geneids) > 1)) {
        stop("Gene IDs must not be a list")
    }
    if (!anyDuplicated(geneids)) {
        return(names(df))
    }
    ## Forget blacklisted columns
    df <- df[setdiff(names(df), blacklist)]
    ## Forget list columns
    df <- df[sapply(df, . %>% lengths %>% max) == 1]
    ## Forget empty columns
    df <- df[!sapply(df, is.empty)]
    if (ncol(df) < 1) {
        return(character(0))
    }
    ## Convert to Rle
    df <- DataFrame(lapply(df, . %>% unlist %>% Rle))
    geneids %<>% Rle
    genecols <- sapply(df, . %>% split(geneids) %>% runLength %>% lengths %>% max %>% is_weakly_less_than(1))
    names(which(genecols))
}

## Given a GRangesList whose underlying ranges have mcols, find mcols
## of the ranges that are constant within each gene and promote them
## to mcols of the GRangesList. For example, if exons are annotated with
promote.common.mcols <- function(grl, delete.from.source=FALSE, ...) {
    colnames.to.promote <- get.gene.common.colnames(mcols(unlist(grl)), rep(names(grl), lengths(grl)), ...)
    promoted.df <- mcols(unlist(grl))[cumsum(lengths(grl)),colnames.to.promote, drop=FALSE]
    if (delete.from.source) {
        mcols(grl@unlistData) %<>% .[setdiff(names(.), colnames.to.promote)]
    }
    mcols(grl) %<>% cbind(promoted.df)
    grl
}

get.txdb <- function(txdbname) {
    tryCatch({
        library(txdbname, character.only=TRUE)
        pos <- str_c("package:", txdbname)
        get(txdbname, pos)
    }, error=function(...) {
        library(GenomicFeatures)
        loadDb(txdbname)
    })
}

read.table.general <- function(filename, read.table.args=NULL, read.xlsx.args=NULL,
                               dataframe.class="data.frame") {
    suppressWarnings({
        read.table.args %<>% as.list
        read.table.args$file <- filename
        read.table.args$header <- TRUE
        read.xlsx.args %<>% as.list
        read.xlsx.args$xlsxFile <- filename
        lazy.results <- list(
            rdata=future(read.RDS.or.RDA(filename, dataframe.class), lazy=TRUE),
            table=future(do.call(read.table, read.table.args), lazy=TRUE),
            csv=future(do.call(read.csv, read.table.args), lazy=TRUE),
            xlsx=future(do.call(read.xlsx, read.xlsx.args), lazy=TRUE))
        for (lzresult in lazy.results) {
            result <- tryCatch({
                x <- as(value(lzresult), dataframe.class)
                assert_that(is(x, dataframe.class))
                x
            }, error=function(...) NULL)
            if (!is.null(result)) {
                return(result)
            }
        }
        stop(glue("Could not read a data frame from {deparse_onestring(filename)} as R data, csv, or xlsx"))
    })
}

read.additional.gene.info <- function(filename, gff_format="GFF3", geneFeatureType="gene", ...) {
    df <- tryCatch({
        gff <- tryCatch({
            read.RDS.or.RDA(filename, "GRanges")
        }, error=function(...) {
            import(filename, format=gff_format)
        })
        assert_that(is(gff, "GRanges"))
        gff %>% .[.$type %in% geneFeatureType] %>%
            mcols %>% cleanup.mcols(mcols_df=.)
    }, error=function(...) {
        tab <- read.table.general(filename, ..., dataframe.class="DataFrame")
        ## Nonexistent or automatic row names
        if (.row_names_info(tab) <= 0) {
            row.names(tab) <- tab[[1]]
        }
        tab
    })
    df %<>% DataFrame
    assert_that(is(df, "DataFrame"))
    return(df)
}

print.var.vector <- function(v) {
    for (i in names(v)) {
        cat(i, ": ", deparse_onestring(v[[i]]), "\n", sep="")
    }
    invisible(v)
}

## Convert strand to -1, 0, or 1
strand.sign <- function(x, allow.unstranded=FALSE) {
    s <- strand(x)
    ss <- (s == "+") - (s == "-")
    if (allow.unstranded) {
        ss[ss == 0] <- NA
    } else if (any(unlist(ss == 0))) {
        stop("Strand must be either '+' or '-'")
    }
    ss
}

{
    cmdopts <- get.options(commandArgs(TRUE))
    cmdopts$help <- NULL

    ## ## For testing only
    ## cmdopts <- list(
    ##     transcript_quant="saved_data/SummarizedExperiment_rnaseq_transcript_shoal_hg38.analysisSet_knownGene.RDS",
    ##     annotation_txdb="TxDb.Hsapiens.UCSC.hg38.knownGene",
    ##     additional_gene_info="/home/ryan/references/hg38/genemeta.org.Hs.eg.db.RDS",
    ##     output_file="test.rds")

    tsmsg("Args:")
    print.var.vector(cmdopts)

    ## Delete the output file if it exists
    suppressWarnings(file.remove(cmdopts$output_file))
    assert_that(!file.exists(cmdopts$output_file))

    ## Only chr1-chr22,chrX,chrY
    std.chr <- extractSeqlevels("Homo sapiens", "UCSC") %>% setdiff("chrM")

    tsmsg("Reading quantification data")
    sexp <- readRDS(cmdopts$transcript_quant)
    sexp %<>% keepSeqlevels(std.chr, pruning.mode="coarse")

    tsmsg("Reading annotation data")
    txdb <- get.txdb(cmdopts$annotation_txdb)

    tsmsg("Computing average transcript abundances")
    tx <- rowRanges(sexp)
    tx$abundance <- sexp %>% assay("abundance") %>% rowMeans
    tx$GeneID <- mapIds(txdb, names(tx),  keytype="TXNAME", column="GENEID", multiVals="first")

    tsmsg("Grouping transcripts by TSS and gene ID")
    tss_table <- tx %>% promoters(upstream=0, downstream=1) %>% as("data.frame") %>%
        filter(!is.na(GeneID)) %>%
        group_by(GeneID, seqnames, start, end, strand) %>%
        summarize(transcript=str_c(transcript, collapse=","),
                  abundance=sum(abundance))

    tsmsg("Selecting most abundant TSS for each gene")
    abundant_tss <- tss_table %>%
        arrange(desc(abundance)) %>% filter(!duplicated(GeneID)) %>%
        as("GRanges") %>% setNames(.$GeneID) %>%
        .[unique(tss_table$GeneID)]

    if ("additional_gene_info" %in% names(cmdopts)) {
        tsmsg("Reading additional gene annotation metadata")
        additional_gene_info <- read.additional.gene.info(cmdopts$additional_gene_info)
        ## Generate empty rows for genes that don't have additional
        ## info
        genes_without_info <- setdiff(abundant_tss$GeneID, rownames(additional_gene_info))
        if (length(genes_without_info) > 0) {
            empty_row <- list(character(0)) %>% rep(ncol(additional_gene_info)) %>% setNames(colnames(additional_gene_info))
            single.val.cols <- sapply(additional_gene_info, function(x) all(lengths(x) == 1))
            for (i in seq_along(empty_row)) {
                if (single.val.cols[i]) {
                    empty_row[[i]] <- NA
                } else {
                    empty_row[[i]] <- list(logical(0)) %>% as(class(additional_gene_info[[i]]))
                }
            }
            empty_row %<>% DataFrame
            empty_gene_table <- empty_row[rep(1, length(genes_without_info)),] %>%
                set_rownames(genes_without_info)
            additional_gene_info %<>% rbind(empty_gene_table)
        }
        assert_that(all(abundant_tss$GeneID %in% rownames(additional_gene_info)))
        mcols(abundant_tss)[colnames(additional_gene_info)] <- additional_gene_info[abundant_tss$GeneID,]
        metadata(abundant_tss) %<>% c(metadata(additional_gene_info))
    }

    tsmsg("Saving output file")
    save.RDS.or.RDA(abundant_tss, cmdopts$output_file)
    invisible(NULL)
}
