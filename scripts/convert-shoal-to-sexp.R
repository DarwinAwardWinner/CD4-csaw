#!/usr/bin/env Rscript

library(getopt)
library(optparse)
library(magrittr)
library(assertthat)
library(glue)

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
        stop(glue("{deparse(argname)} must be NULL or a character vector"))
    if (!several.ok) {
        if (identical(arg, choices))
            return(arg[1L])
        if (length(arg) > 1L)
            stop(glue("{deparse(argname)} must be of length 1"))
    }
    else if (length(arg) == 0L)
        stop(glue("{deparse(argname)} must be of length >= 1"))
    fold_case <- identity
    if (ignore.case) {
        fold_case <- tolower
    }
    i <- pmatch(fold_case(arg), fold_case(choices), nomatch = 0L, duplicates.ok = TRUE)
    if (all(i == 0L))
        stop(gettextf("%s should be one of %s", deparse(argname), paste(dQuote(choices),
            collapse = ", ")), domain = NA)
    i <- i[i > 0L]
    if (!several.ok && length(i) > 1)
        stop("there is more than one match in 'match.arg'")
    choices[i]
}

get.options <- function(opts) {
    optlist <- list(
        make_option(c("-s", "--samplemeta-file"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) RDS/RData/xlsx/csv file containing a table of sample metadata. Any existing rownames will be replaced with the values in the sample ID  column (see below)."),
        make_option(c("-c", "--sample-id-column"), type="character", default="Sample",
                    help="Sample metadata column name that holds the sample IDs. These will be used to determine the abundance file names."),
        make_option(c("-i", "--shoal-dir"), metavar="PATTERN", type="character",
                    help="(REQUIRED) Directory containing shoal output."),
        make_option(c("-l", "--aggregate-level"), metavar="LEVEL", type="character", default="auto",
                    help="Whether to save aggregated gene counts or transcript counts in the output file. By default, aggregated gene counts are saved if a gene annotation is provided, and transcript counts are saved otherwise. You can force one or the other by specifying 'gene' or 'transcript' for this option."),
        make_option(c("-o", "--output-file"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) Output file name. The SummarizedExperiment object containing the counts will be saved here using saveRDS, so it should end in '.RDS'."),
        make_option(c("-m", "--genemap-file"), metavar="FILENAME", type="character",
                    help="Genemap file in the Salmon simple gene map format (see 'salmon quant --help-reads')"),
        make_option(c("-d", "--annotation-txdb"), metavar="PACKAGE_OR_FILE_NAME", type="character",
                    help="Name of TxDb package, or the name of a database file, to use for gene annotation"),
        make_option(c("-g", "--gene-info"), metavar="FILENAME", type="character",
                    help="RDS/RData/xlsx/csv file containing a table of gene metadata. Row names (or the first column of the file if there are no row names) should be gene/feature IDs that match the ones used in the main annotation, and these should be unique. This option is ignored when not aggregating counts to the gene level."),
        make_option(c("--transcript-info"), metavar="FILENAME", type="character",
                    help="RDS/RData/xlsx/csv file containing a table of transcript metadata. Row names (or the first column of the file if there are no row names) should be transcript IDs that match the ones used in the quantification files, and these should be unique. This option is ignored when aggregating counts to the gene level."))
    progname <- na.omit(c(get_Rscript_filename(), "convert-quant-to-sexp.R"))[1]
    parser <- OptionParser(
        usage="Usage: %prog [ -d TXDB | -m GENEMAP ] [ -g GENEINFO | -t TXINFO ] -s SAMPLEMETA.RDS -a PATTERN -l (gene|transcript) -o SUMEXP.RDS",
        description="Collect RNA-seq quantification results into a SummarizedExperiment object.

TODO UPDATE Counts are stored along with the sample and gene metadata in a SummarizedExperiment object. Note that the '-s', '-a', '-t', and '-o' arguments are all required, since they specify the essential input and output files and formats.",
option_list = optlist,
add_help_option = TRUE,
prog = progname,
epilogue = "")

    cmdopts <- parse_args(parser, opts)
    ## Ensure that all required arguments were provided
    required.opts <- c("samplemeta-file", "shoal-dir", "output-file")
    missing.opts <- setdiff(required.opts, names(cmdopts))
    if (length(missing.opts) > 0) {
        stop(str_c("Missing required arguments: ", deparse(missing.opts)))
    }

    ## Ensure that no more than one annotation was provided, and that
    ## exactly one annotation was provided if gene-level aggregation
    ## was requested.
    annot.opts <- c("annotation-txdb", "genemap-file")
    provided.annot.opts <- intersect(annot.opts, names(cmdopts))
    if (length(provided.annot.opts) > 1) {
        stop("Multiple gene annotations were provided. Please provide only one.")
    }
    quant.level.options <- c("auto", "gene", "transcript", "tx")
    cmdopts[['aggregate-level']] %<>% tolower %>% match.arg(choices=quant.level.options, argname="--aggregate-level", ignore.case=TRUE)
    if (cmdopts[['aggregate-level']] == "auto") {
        cmdopts[['aggregate-level']] = ifelse(length(provided.annot.opts) == 1, "gene", "transcript")
    }
    ## "tx" is an undocumented shortcut for "transcript", for
    ## consistency with tximport, TxDb, etc.
    if (cmdopts[['aggregate-level']] == "tx") {
        cmdopts[['aggregate-level']] <- "transcript"
    }
    assert_that(cmdopts[['aggregate-level']] %in% c("gene", "transcript"))
    if (cmdopts[['aggregate-level']] == "gene" && length(provided.annot.opts) < 1) {
        stop("Gene-level quantification was requested but no gene annotations were provided")
    }

    ## Replace dashes with underscores so that all options can easily
    ## be accessed by "$"
    cmdopts %>% setNames(chartr("-", "_", names(.)))
}

## Do argument parsing early so the script exits quickly if arguments are invalid
get.options(commandArgs(TRUE))

library(assertthat)
library(dplyr)
library(future)
library(magrittr)
library(openxlsx)
library(stringr)
library(future)

library(annotate)
library(GenomicRanges)
library(rtracklayer)
library(S4Vectors)
library(SummarizedExperiment)
library(tximport)

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

## TODO: Move to utils package

## Read a table from a R data file, csv, or xlsx file. Returns a data
## frame or thorws an error.
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
            table=futureCall(read.table, read.table.args, lazy=TRUE),
            csv=futureCall(read.csv, read.table.args, lazy=TRUE),
            xlsx=futureCall(read.xlsx, read.xlsx.args, lazy=TRUE))
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
        stop(glue("Could not read a data frame from {deparse{filename}} as R data, csv, or xlsx"))
    })
}

cleanup.mcols <- function(object, mcols_df=mcols(object)) {
    nonempty <- !sapply(mcols_df, is.empty)
    mcols_df %<>% .[nonempty]
    if (!missing(object)) {
        mcols(object) <- mcols_df
        return(object)
    } else {
        return(mcols_df)
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

## This merges exons into genes (GRanges to GRangesList)
gff.to.grl <- function(gr, exonFeatureType="exon", geneIdAttr="gene_id", geneFeatureType="gene") {
    exon.gr <- gr[gr$type %in% exonFeatureType]
    exon.gr %<>% cleanup.mcols
    grl <- split(exon.gr, as.character(mcols(exon.gr)[[geneIdAttr]])) %>%
        promote.common.mcols
    if (!is.null(geneFeatureType)) {
        gene.meta <- gr[gr$type %in% geneFeatureType] %>%
            mcols %>% cleanup.mcols(mcols_df=.) %>% .[match(names(grl), .[[geneIdAttr]]),]
        for (i in names(gene.meta)) {
            if (i %in% names(mcols(grl))) {
                value <- ifelse(is.na(gene.meta[[i]]), mcols(grl)[[i]], gene.meta[[i]])
            } else {
                value <- gene.meta[[i]]
            }
            mcols(grl)[[i]] <- value
        }
    }
    return(grl)
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

get.tx2gene.from.txdb <- function(txdb) {
    k <- keys(txdb, keytype = "GENEID")
    suppressMessages(AnnotationDbi::select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")) %>%
        .[c("TXNAME", "GENEID")]
}

read.tx2gene.from.genemap <- function(fname) {
    df <- read.table.general(fname)
    df %<>% .[1:2]
    df[] %<>% lapply(as.character)
    names(df) <- c("TXNAME", "GENEID")
    df
}

read.annotation.from.gff <- function(filename, format="GFF3", ...) {
    gff <- NULL
    ## Allow the file to be an RDS file containing the GRanges
    ## resulting from import()
    gff <- tryCatch({
        read.RDS.or.RDA(filename, "GRanges")
    }, error=function(...) {
        import(filename, format=format)
    })
    assert_that(is(gff, "GRanges"))
    grl <- gff.to.grl(gff, ...)
    return(grl)
}

read.annotation.from.saf <- function(filename, ...) {
    saf <- read.table.general(filename, ...)
    assert_that("GeneID" %in% names(saf))
    gr <- as(saf, "GRanges")
    grl <- split(gr, gr$GeneID) %>% promote.common.mcols
    return(grl)
}

read.annotation.from.rdata <- function(filename) {
    read.RDS.or.RDA(filename, "GRangesList")
}

read.additional.gene.info <- function(filename, gff_format="GFF3", geneFeatureType="gene", ...) {
    df <- tryCatch({
        gff <- tryCatch({
            read.RDS.or.RDA(filename, "GRanges")
        }, error=function(...) {
            import(filename, format=gff_format)
        })
        assert_that(is(gff, "GRanges"))
        if (!is.null(geneFeatureType)) {
            gff %<>% .[.$type %in% geneFeatureType]
        }
        gff %<>% .[!is.na(.$ID) & !duplicated(.$ID)]
        gff %>% mcols %>% cleanup.mcols(mcols_df=.)
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

## This converts a GRangesList into the SAF ("Simplified annotation
## format")
grl.to.saf <- function(grl) {
    gr <- unlist(grl)
    data.frame(Chr=as.vector(seqnames(gr)),
               Start=start(gr),
               End=end(gr),
               Strand=as.vector(strand(gr)),
               GeneID=rep(names(grl), lengths(grl)))
}

print.var.vector <- function(v) {
    for (i in names(v)) {
        cat(i, ": ", deparse(v[[i]]), "\n", sep="")
    }
    invisible(v)
}

{

    cmdopts <- get.options(commandArgs(TRUE))
    ## myargs <- c('--samplemeta-file', 'saved_data/samplemeta-RNASeq.RDS',
    ##             '--sample-id-column', 'SRA_run',
    ##             '--shoal-dir', 'shoal_quant/hg38.analysisSet_ensembl.85',
    ##             '--output-file', 'saved_data/SummarizedExperiment_rnaseq_shoal_hg38.analysisSet_ensembl.85.RDS',
    ##             '--aggregate-level', 'gene',
    ##             '--annotation-txdb', '/home/ryan/references/hg38/TxDb.Hsapiens.ensembl.hg38.v85.sqlite3',
    ##             '--gene-info', '/home/ryan/references/hg38/genemeta.ensembl.85.RDS')
    ## cmdopts <- get.options(myargs)
    cmdopts$help <- NULL

    tsmsg("Args:")
    print.var.vector(cmdopts)

    ## Delete the output file if it exists
    suppressWarnings(file.remove(cmdopts$output_file))
    assert_that(!file.exists(cmdopts$output_file))

    tsmsg("Loading sample metadata")
    samplemeta <- read.table.general(cmdopts$samplemeta_file)

    tsmsg("Got metadata for ", nrow(samplemeta), " samples")

    assert_that(cmdopts$sample_id_column %in% colnames(samplemeta))
    assert_that(!anyDuplicated(samplemeta[[cmdopts$sample_id_column]]))

    rownames(samplemeta) <- samplemeta$sample <- samplemeta[[cmdopts$sample_id_column]]

    samplemeta$path <- file.path(cmdopts$shoal_dir, glue("{samplemeta[[cmdopts$sample_id_column]]}_adapt.sf"))

    assert_that(all(file.exists(samplemeta$path)))

    annot <- NULL
    annot_ranges <- NULL
    tx2gene <- NULL
    if (cmdopts$aggregate_level == "gene") {
        if ("annotation_txdb" %in% names(cmdopts)) {
            tsmsg("Reading transcript ranges from TxDb")
            txdb <- get.txdb(cmdopts$annotation_txdb)
            annot_ranges <- transcriptsBy(txdb, "gene")
            tsmsg("Reading transcript gene mappings from TxDb")
            tx2gene <- get.tx2gene.from.txdb(txdb)
        } else if ("genemap_file" %in% names(cmdopts)) {
            tsmsg("Reading transcript gene mappings from genemap file")
            tx2gene <- read.tx2gene.from.genemap(cmdopts$genemap_file)
        } else {
            stop("Need a gene annotation to aggregate at the gene level.")
        }
        if ("gene_info" %in% names(cmdopts)) {
            tsmsg("Reading gene annotations")
            annot <- read.table.general(cmdopts$gene_info, dataframe.class="DataFrame")
            ## Nonexistent or automatic row names
            if (.row_names_info(annot) <= 0) {
                row.names(annot) <- annot[[1]]
            }
        }
    } else {
        if ("annotation_txdb" %in% names(cmdopts)) {
            tsmsg("Reading transcript ranges from TxDb")
            txdb <- get.txdb(cmdopts$annotation_txdb)
            annot_ranges <- transcripts(txdb)
            names(annot_ranges) <- annot_ranges$tx_name
        }
        if ("transcript_info" %in% names(cmdopts)) {
            tsmsg("Reading transcript annotations")
            annot <- read.table.general(cmdopts$transcript_info, dataframe.class="DataFrame")
            ## Nonexistent or automatic row names
            if (.row_names_info(annot) <= 0) {
                row.names(annot) <- annot[[1]]
            }
        }
    }

    tsmsg("Reading quantification files")
    txi <- tximport(samplemeta$path, type="salmon", txOut=TRUE)
    if (cmdopts$aggregate_level == "gene") {
        txi %<>% summarizeToGene(tx2gene)
    }

    tsmsg("Matching annotation to quantification tables")
    txi_assayNames <- c("counts", "abundance", "length")
    txi_featureNames <- rownames(txi[[txi_assayNames[1]]])
    if (is.null(annot)) {
        annot <- DataFrame(row.names=txi_featureNames)
        annot[[cmdopts$aggregate_level]] <- txi_featureNames
    } else {
        annot %<>% .[txi_featureNames,] %>% set_rownames(txi_featureNames)
    }

    ## If we have ranges (from a TxDb), we make a
    ## RangedSummarizedExperiment, otherwise we make a regular one
    if (!is.null(annot_ranges)) {
        tsmsg("Constructing the RangedSummarizedExperiment object")
        annot_ranges %<>% .[txi_featureNames]
        mcols(annot_ranges) <- as(annot, "DataFrame")
        sexp <- SummarizedExperiment(
            assays=List(txi[txi_assayNames]),
            colData=as(samplemeta, "DataFrame"),
            rowRanges=annot_ranges,
            ## Put non-assay elements of txi into the metadata
            metadata=SimpleList(txi[!names(txi) %in% txi_assayNames]))
    } else {
        tsmsg("Constructing the SummarizedExperiment object")
        sexp <- SummarizedExperiment(
            assays=List(txi[txi_assayNames]),
            colData=as(samplemeta, "DataFrame"),
            rowData=as(annot, "DataFrame"),
            ## Put non-assay elements of txi into the metadata
            metadata=SimpleList(txi[!names(txi) %in% txi_assayNames]))
    }

    tsmsg("Saving SummarizedExperiment")
    save.RDS.or.RDA(sexp, cmdopts$output_file)
    invisible(NULL)
}
