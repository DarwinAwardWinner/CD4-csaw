#!/usr/bin/env Rscript

library(getopt)
library(optparse)
library(stringr)
library(magrittr)
library(assertthat)

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

## Kilo, mega, giga, tera
size.prefixes <- c(k=1e3, m=1e6, g=1e9, t=1e12)

parse.bp <- function(size) {
    m <- str_match(size, "^(.*?)(?:([kmgt]?)bp?)?\\s*$")
    ## m <- str_match(size, "^\\s*(\\d+(?:\\.\\d+))\\s*(?:([kmgt]?)bp?)?\\s*$")
    base <- suppressWarnings(as.numeric(m[,2]))
    if (any(is.na(base))) {
        invalid <- size[is.na(base)]
        stop(sprintf("Invalid base pair size specification: %s", deparse(head(invalid))))
    }
    multiplier <- unname(size.prefixes[m[,3]])
    multiplier[is.na(multiplier)] <- 1
    x <- base * multiplier
    assert_that(!any(is.na(x)))
    x
}

get.options <- function(opts) {
    optlist <- list(
        make_option(c("-s", "--samplemeta-file"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) RDS/RData/xlsx/csv file containing a table of sample metadata. Any existing rownames will be replaced with the values in the sample ID  column (see below)."),
        make_option(c("-c", "--sample-id-column"), type="character", default="Sample",
                    help="Sample metadata column name that holds the sample IDs. These will be substituted into '--bam-file-pattern' to determine the BAM file names."),
        make_option(c("-p", "--bam-file-pattern"), metavar="PATTERN", type="character",
                    help="(REQUIRED) Format string to convert sample IDs into BAM file paths. This should contain a '%s' wherever the sample ID should be substituted ('%s' can occur multiple times),. Example: 'bam_files/Sample_%s/Aligned.bam"),
        make_option(c("-r", "--regions"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) File specifying regions in which reads should be counted. This can be a BED file, GFF file, R data file containing a GRanges object, or csv file that can be converted to a GRanges object. If the regions have associated annotations, then a GRanges in an R data file is the recommended format."),
        make_option(c("-o", "--output-file"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) Output file name. The SummarizedExperiment object containing the counts will be saved here using saveRDS, so it should end in '.RDS'."),
        make_option(c("-b", "--expected-bam-files"), metavar="BAMFILE1,BAMFILE2,...", type="character",
                    help="Comma-separated list of bam file names expected to be used as input. This argument is optional, but if it is provided, it will be checked against the list of files determined from '--samplemeta-file' and '--bam-file-pattern', and an error will be raised if they don't match exactly."),
        make_option(c("-e", "--read-extension"), type="character", default="100bp",
                    help="Assumed fragment length of reads. Each read will be assumed to represent a DNA fragment extending this far from its 5 prime end, regardless of the actual read length."),
        make_option(c("-x", "--blacklist"), metavar="FILENAME.bed", type="character",
                    help="File describing blacklist regions to be excluded from the analysis. Reads that overlap these regions will be discarded without counting them toward any region. This can be a BED file, GFF file, R data file containing a GRanges object, or csv file that can be converted to a GRanges object."),
        make_option(c("-j", "--threads"), metavar="N", type="integer", default=1,
                    help="Number of threads to use"))
    progname <- na.omit(c(get_Rscript_filename(), "csaw-count-windows.R"))[1]
    parser <- OptionParser(
        usage="Usage: %prog [options] -s SAMPLEMETA.RDS -p PATTERN -w WSIZE -e READEXT [ -s WSPACE ] -o SUMEXP.RDS",
        description="Count ChIP-seq reads a set of specified regions",
        option_list = optlist,
        add_help_option = TRUE,
        prog = progname,
        epilogue = "Note that all base pair sizes (window width/spacing and read extension) may have an suffix of 'bp', 'kbp', 'mbp', or 'tbp'. For example, 10kb = 10000.")

    cmdopts <- parse_args(parser, opts)
    ## Ensure that all required arguments were provided
    required.opts <- c("samplemeta-file", "output-file", "bam-file-pattern", "regions")
    missing.opts <- setdiff(required.opts, names(cmdopts))
    if (length(missing.opts) > 0) {
        stop(str_c("Missing required arguments: ", deparse(missing.opts)))
    }

    ## Convert bp args to numbers
    for (i in c("read-extension")) {
        cmdopts[[i]] %<>% parse.bp
    }

    cmdopts$threads %<>% round
    assert_that(cmdopts$threads >= 1)

    cmdopts$help <- NULL

    ## Replace dashes with underscores so that all options can easily
    ## be accessed by "$"
    cmdopts %>% setNames(str_replace_all(names(.), "-", "_"))
}

## Terminate early on argument-processing errors
invisible(get.options(commandArgs(TRUE)))

library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(dplyr)
library(csaw)
library(Matrix)
library(future)

format.bp <- function(x) {
    x %<>% round
    sapply(x, function(xx) {
        ## Anything more than this, just use scientific notation
        if ((xx) >= 1e15) {
            prefix <= ""
            multiplier <- 1
        } else {
            prefix <- size.prefixes %>% .[. <= xx] %>% which.max %>% names
            if (is.null(prefix)) {
                prefix <- ""
                multiplier <- 1
            } else {
                multiplier <- size.prefixes[prefix]
            }
        }
        digits <- max(1, ceiling(log10(xx+1)))
        numfmt <- str_c("%.", digits, "g%sbp")
        sprintf(numfmt, xx / multiplier, prefix)
    })
}

regionCountsParallel <- function(bam.files, ..., BPPARAM=bpparam()) {
    reslist <- bplapply(X=bam.files, FUN=regionCounts, ..., BPPARAM=BPPARAM)
    assert_that(all(sapply(reslist, is, "SummarizedExperiment")))
    res <- do.call(cbind, reslist)
    rm(reslist)
    res
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
    if (!any(sapply(expected.class, is, object=object))) {
        object <- as(object, expected.class)
    }
    return(object)
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

read.saf <- function(filename, ...) {
    saf <- read.table.general(filename, ...)
    assert_that("GeneID" %in% names(saf))
    gr <- as(saf, "GRanges")
    grl <- split(gr, gr$GeneID) %>% promote.common.mcols
    return(grl)
}

read.regions <- function(filename) {
    suppressWarnings({
        lazy.results <- list(
            rdata=lazy(read.RDS.or.RDA(filename)),
            bed=lazy(import(filename, format="bed")),
            gff=lazy(import(filename, format="gff")),
            saf=lazy(read.saf(filename)),
            table=lazy(read.table.general(filename)))
        for (lzresult in lazy.results) {
            result <- tryCatch({
                x <- value(lzresult)
                if (is(x, "List")) {
                    x <- unlist(x)
                }
                x <- as(x, "GRanges")
                x
            }, error=function(...) NULL)
            if (!is.null(result)) {
                return(result)
            }
        }
        stop(sprintf("Could not read genomic regions from %s as R data, bed, gff, SAF, or csv", deparse(filename)))
    })
}

print.var.vector <- function(v) {
    for (i in names(v)) {
        cat(i, ": ", deparse(v[[i]]), "\n", sep="")
    }
    invisible(v)
}

{
    cmdopts <- get.options(commandArgs(TRUE))
    tryCatch(setwd(file.path(dirname(na.omit(get_Rscript_filename())), "..")),
             error=function(...) tsmsg("WARNING: Could not determine script path. Ensure that you are already in the correct directory."))

    tsmsg("Args:")
    print.var.vector(cmdopts)

    if (cmdopts$threads > 1) tryCatch({
        library(doParallel)
        library(BiocParallel)
        registerDoParallel(cores=cmdopts$threads)
        register(DoparParam())
    }, error=function(...){
        tsmsg("Could not initialize parallel backend. Falling back to single-core mode.")
        cmdopts$threads <- 1
    })
    tsmsg("Using ", cmdopts$threads, " cores.")

    tsmsg("Assuming a fragment size of ", format.bp(cmdopts$read_extension))

    tsmsg("Loading sample data")

    sample.table <- readRDS(cmdopts$samplemeta_file) %>%
        ## Compute full path to BAM file
        mutate(bam_file=sprintf(cmdopts$bam_file_pattern, .[[cmdopts$sample_id_column]])) %>%
        ## Ensure that days_after_activation is a factor and can't be
        ## interpreted as a numeric
        mutate(days_after_activation=days_after_activation %>%
                   factor %>% `levels<-`(str_c("Day", levels(.)))) %>%
        rename(time_point=days_after_activation)

    assert_that(all(file.exists(sample.table$bam_file)))

    tsmsg("Loading regions")
    target.regions <- read.regions(cmdopts$regions)
    assert_that(is(target.regions, "GRanges"))
    ## Analysis is not stranded
    strand(target.regions) <- "*"

    blacklist.regions <- GRanges()
    if (!is.null(cmdopts$blacklist)) {
        tsmsg("Loading blacklist regions")
        blacklist.regions <- read.regions(cmdopts$blacklist)
        assert_that(is(blacklist.regions, "GRanges"))
        ## Blacklist applies to both strands
        strand(blacklist.regions) <- "*"
    }
    rparam <- readParam(discard=blacklist.regions)

    tsmsg(sprintf("Counting reads in %i regions in %i samples.",
                  length(target.regions), nrow(sample.table)))
    if (cmdopts$threads > 1) {
        rCountsFun <- regionCountsParallel
        options(mc.cores=cmdopts$threads, mc.preschedule=FALSE)
        registerDoParallel(cores=cmdopts$threads)
        register(DoparParam())
    } else {
        rCountsFun <- regionCounts
    }
    wcounts <- rCountsFun(
        sample.table$bam_file, regions=target.regions,
        ext=cmdopts$read_extension, param=rparam)
    colData(wcounts) %<>% {cbind(sample.table, .[c("totals", "ext")])}
    saveRDS(wcounts, cmdopts$output_file)
}
