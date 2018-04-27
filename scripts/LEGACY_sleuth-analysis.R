#!/usr/bin/env Rscript

library(getopt)
library(optparse)

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

known.gene.id.types <- c(entrez="ENTREZID", ensembl="ENSEMBL", unigene="UNIGENE", symbol="SYMBOL")
gene.id.type.metavar <- sprintf("(%s)", str_c(names(known.gene.id.types), collapse="|"))

get.options <- function(opts) {

    ## Do argument parsing early so the script exits quickly if arguments are invalid
    optlist <- list(
        make_option(c("-s", "--samplemeta-file"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) RDS/RData/xlsx/csv file containing a table of sample metadata. Any existing rownames will be replaced with the values in the sample ID  column (see below)."),
        make_option(c("-c", "--sample-id-column"), type="character", default="Sample",
                    help="Sample metadata column name that holds the sample IDs. These will be substituted into '--abundance-file-pattern' to determine the abundance file names."),
        make_option(c("-a", "--abundance-file-pattern"), metavar="PATTERN", type="character",
                    help="(REQUIRED) Format string to convert sample IDs into file paths to the abundance.h5 file for each sample. This should contain a '%s' wherever the sample ID should be substituted ('%s' can occur multiple times),. Example: 'kallisto_quant/Sample_%s/abundance.h5"),
        make_option(c("-o", "--output-file"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) Output file name. The SummarizedExperiment object containing the counts will be saved here using saveRDS, so it should end in '.RDS'."),
        make_option(c("-b", "--expected-abundance-files"), metavar="FILE1.h5,FILE2.h5,...", type="character",
                    help="Comma-separated list of file names expected to be used as input. This argument is optional, but if it is provided, it will be checked against the list of files determined from '--samplemeta-file' and '--abundance-file-pattern', and an error will be raised if they don't match exactly."),
        ## TODO: Allow different annotations, via txdb, or gff file
        ## make_option(c("-g", "--annotation-gff"), metavar="FILENAME", type="character",
        ##             help="File Name of GFF3 file to use for gene annotation."),
        ## make_option(c("-f", "--gff-transcript-featuretype"), metavar="FEATURETYPE", type="character", default="exon",
        ##             help="GFF feature type to use"),
        ## make_option(c("-i", "--gff-geneid-attr"), metavar="ATTRNAME", type="character", default="gene_id",
        ##             help="GFF feature attribute to use as a feature's Gene ID."),
        ## make_option(c("-e", "--gff-gene-featuretype"), metavar="FEATURETYPE", type="character", default="gene",
        ##             help="GFF feature type from which gene metadata should be extracted."),
        ## make_option(c("-r", "--annotation-rds"), metavar="FILENAME", type="character",
        ##             help="File Name of RDS or RData file to use for gene annotation. It should contain a single GRanges or GRangesList object (or something that can be coereced into one), with each element representing a gene/feature to be counted. Any metadata columns on the object will be carried through to the output SummarizedExperiment."),
        ## make_option(c("--annotation-saf"), metavar="FILENAME", type="character",
        ##             help="File Name of RDS/RData/xlsx/csv file containing gene annotations in SAF format (i.e. 5 columns named GeneID, Chr, Start, End, Strand). Additional columns beyond the first 5 will be retained as metadata columns on the genes/exons."),
        ## make_option(c("-a", "--additional-gene-info"), metavar="FILENAME", type="character",
        ##             help="RDS/RData/xlsx/csv file containing a table of gene metadata. Row names (or the first column of the file if there are no row names) should be gene/feature IDs that match the ones used in the main annotation, and these should be unique. This can also be a GFF3 file where the metadata is in the attributes of elements of type specified by '--gff-gene-featuretype' ('gene' by default), where the 'ID' attribute specifies the gene ID.")
        )
    progname <- na.omit(c(get_Rscript_filename(), "sleuth-analysis.R"))[1]
    parser <- OptionParser(
        usage="Usage: %prog [options] -s SAMPLEMETA.RDS -p PATTERN -o SUMEXP.RDS [ -t TXDB.PACKAGE.NAME | -g ANNOTATION.GFF3 | -r ANNOTATION.RDS ]",
        description="Count reads in genes using Rsubread::featureCounts.

Counts are performed for stranded, non-stranded, and reverse-stranded modes, and are stored along with the sample and gene metadata in a SummarizedExperiment object. Note that the '-s', '-p', and '-o' arguments are all required, since they specify the input and output files. Also, exactly one of '-t', '-g', or '-r' is required to specify the annotation.",
option_list = optlist,
add_help_option = TRUE,
prog = progname,
epilogue = "")

    cmdopts <- parse_args(parser, opts)
    ## Ensure that all required arguments were provided
    required.opts <- c("samplemeta-file", "abundance-file-pattern", "output-file")
    missing.opts <- setdiff(required.opts, names(cmdopts))
    if (length(missing.opts) > 0) {
        stop(str_c("Missing required arguments: ", deparse(missing.opts)))
    }

    ## ## Ensure that exactly one annotation was provided
    ## annot.opts <- c("annotation-txdb", "annotation-gff", "annotation-rds", "annotation-saf")
    ## provided.annot.opts <- intersect(annot.opts, names(cmdopts))
    ## if (length(provided.annot.opts) < 1) {
    ##     stop("No annotation provided")
    ## } else if (length(provided.annot.opts) > 1) {
    ##     stop("Multiple annotations provided. Please provide only one annotation.")
    ## }

    ## if ("txdb-geneid-type" %in% cmdopts) {
    ##     cmdopts["txdb-geneid-type"] %<>%
    ##         tolower %>%
    ##         match.arg(c(known.gene.id.types, "auto"),
    ##                   argname="txdb-geneid-type")
    ## }
    cmdopts %>% setNames(str_replace_all(names(.), "-", "_"))
}

## Do this early to handle "--help" before wasting time loading
## pacakges & stuff
## get.options(commandArgs(TRUE))

library(assertthat)
library(dplyr)
library(future)
library(magrittr)
library(openxlsx)
library(stringr)

library(annotate)
library(GenomicRanges)
library(Rsubread)
library(sleuth)
library(SummarizedExperiment)

library(org.Hs.eg.db)

library(sleuth)
library(limma)

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
            rdata=lazy(read.RDS.or.RDA(filename, dataframe.class)),
            table=lazy(do.call(read.table, read.table.args)),
            csv=lazy(do.call(read.csv, read.table.args)),
            xlsx=lazy(do.call(read.xlsx, read.xlsx.args)))
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
        stop(sprintf("Could not read a data frame from %s as R data, csv, or xlsx", deparse(filename)))
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

make.lazy <- function(func, ...) {
    lazymaker <- function(expr)
        lazy(expr, ...)
    function(...) {
        lazymaker(func(...))
    }
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

## Like sprintf, but inserts the same value into every placeholder
sprintf.single.value <- function(fmt, value) {
    ## Max function arguments is 100
    arglist = c(list(fmt=fmt), rep(list(value), 99))
    do.call(sprintf, arglist)
}

## Guess type of ID
identify.ids <- function(ids, db="org.Hs.eg.db", idtypes=c("ENTREZID", "ENSEMBL", "UNIGENE"), threshold=0.5) {
    if (is.character(db)) {
        library(db, character.only=TRUE)
        pos <- str_c("package:", db)
        db <- get(db, pos)
    }
    assert_that(is(db, "AnnotationDb"))
    idtypes %<>% intersect(keytypes(db))
    assert_that(length(idtypes) > 0)
    idcounts <- sapply(idtypes, function(idtype) {
        sum(ids %in% keys(db, idtype))
    }) %>% setNames(idtypes)
    idcounts %<>% sort(decreasing = TRUE)
    result <- names(idcounts)[1]
    if (idcounts[result] / length(ids) < threshold) {
        stop(sprintf("Could not identify more than %.0f%% of given IDs as any of %s", threshold * 100, deparse(idtypes)))
    }
    tsmsg("Detected gene IDs as ", result)
    attr(result, "counts") <- idcounts
    return(result)
}

{
    ## cmdopts <- get.options(commandArgs(TRUE))
    myargs <- c("-s", "./saved_data/samplemeta-RNASeq.RDS", "-c", "SRA_run", "-a", "salmon_quant/hg38.analysisSet_ensembl.85/%s/abundance.h5", "-o", "temp.rds", "-j", "4")
    cmdopts <- get.options(myargs)
    cmdopts$help <- NULL

    ## Expand expected_bam_files into vector
    if ("expected_abundance_files" %in% names(cmdopts)) {
        cmdopts$expected_abundance_files %<>% str_split(",") %>% unlist
    }

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

    samplemeta$path <- sprintf.single.value(cmdopts$abundance_file_pattern, samplemeta[[cmdopts$sample_id_column]])

    if ("expected_abundance_files" %in% names(cmdopts)) {
        tryCatch({
            assert_that(setequal(samplemeta$path, cmdopts$expected_abundance_files))
            tsmsg("Got all expected abundance files")
        }, error=function(...) {
            unexpected_existing <- setdiff(samplemeta$path, cmdopts$expected_abundance_files)
            expected_but_missing <- setdiff(cmdopts$expected_abundance_files, samplemeta$path)
            if (length(unexpected_existing) > 0) {
                tsmsg(sprintf("Got unexpected abundance files: %s", deparse(unexpected_existing)))
            }
            if (length(expected_but_missing) > 0) {
                tsmsg(sprintf("Didn't find expected abundance files: %s", deparse(expected_but_missing)))
            }
            stop("Abundance file list was not as expected")
        })
    }

    assert_that(all(file.exists(samplemeta$path)))

    samplemeta %<>%
        mutate(days_after_activation=factor(days_after_activation),
               Group=interaction(cell_type, "D", days_after_activation, sep=""))

    full_design <- model.matrix(~0 + Group + donor_id, samplemeta)
    colnames(full_design) %<>%
        str_replace("Group", "") %>%
        str_replace("donor_idD", "Dn")

    ct <- list(
        NaiveD0v1="NaiveD1 - NaiveD0",
        NaiveD0v5="NaiveD5 - NaiveD0",
        NaiveD0v14="NaiveD14 - NaiveD0",
        NaiveAll=c(
            NaiveD0v1="NaiveD1 - NaiveD0",
            NaiveD0v5="NaiveD5 - NaiveD0",
            NaiveD0v14="NaiveD14 - NaiveD0"),
        MemoryD0v1="MemoryD1 - MemoryD0",
        MemoryD0v5="MemoryD5 - MemoryD0",
        MemoryD0v14="MemoryD14 - MemoryD0",
        MemoryAll=c(
            MemoryD0v1="MemoryD1 - MemoryD0",
            MemoryD0v5="MemoryD5 - MemoryD0",
            MemoryD0v14="MemoryD14 - MemoryD0"),
        NvM.D0="MemoryD0 - NaiveD0",
        NvM.D1="MemoryD1 - NaiveD1",
        NvM.D5="MemoryD5 - NaiveD5",
        NvM.D14="MemoryD14 - NaiveD14",
        Nvm.All=c(
            NvM.D0="MemoryD0 - NaiveD0",
            NvM.D1="MemoryD1 - NaiveD1",
            NvM.D5="MemoryD5 - NaiveD5",
            NvM.D14="MemoryD14 - NaiveD14"),
        DonorVar=c(Donor4659vs5053="Dn5053",
                   Donor4659vs5131="Dn5131",
                   Donor4659vs5291="Dn5291"))

    reduced.designs <- lapply(
        ct, . %>% makeContrasts(contrasts=., levels=full_design) %>%
            contrastAsCoef(full_design, contrast=.) %>%
            {.$design[,-.$coef]})

    all.designs <- c(list(full=full_design, intercept=~1), reduced.designs)

    so <- sleuth_prep(samplemeta, full_design, max_bootstrap=30)
    so$sample_to_covariates %<>% mutate(sample=SampleName)

    ## so <- sleuth_fit(so, full_design, "full")
    ## so <- sleuth_fit(so, ~1, 'intercept')

    for (i in names(all.designs)) {
        tsmsg("Fitting model for ", i)
        so <- sleuth_fit(so, all.designs[[i]], fit_name=i)
    }

    for (i in names(reduced.designs)) {
        tsmsg("Doing LRT for ", i)
        so <- sleuth_lrt(so, i, "full")
    }

    coefs <- lapply(so$fits$full$models, . %$% ols_fit %$% coefficients) %>%
        do.call(what=rbind)

    results <- lapply(names(reduced.designs), function(i) {
        res <- sleuth_results(so, str_c(i, ":full"), test_type="lrt")
        ctmat <- makeContrasts(contrasts=ct[[i]], levels=full_design)
        if (ncol(ctmat) > 1) {
            if (!is.null(names(ct[[i]])))
                colnames(ctmat) <- str_c("logFC_", names(ct[[i]]))
        } else if (ncol(ctmat) == 1) {
            colnames(ctmat) <- "logFC"
        }
        logFCmat <- coefs %*% ctmat
        logFCtable <- as_data_frame(logFCmat) %>% mutate(target_id=rownames(logFCmat))
        left_join(res, logFCtable, by="target_id")
    }) %>% setNames(names(reduced.designs))

    invisible(NULL)
}
