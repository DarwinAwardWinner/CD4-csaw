#!/usr/bin/env Rscript

library(getopt)
library(optparse)
library(stringr)
library(magrittr)
library(GenomicRanges)
library(rtracklayer)
library(SummarizedExperiment)
library(dplyr)
library(csaw)
library(Matrix)
library(assertthat)

library(doParallel)
library(BiocParallel)

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

get.options <- function(opts) {
    optlist <- list(
        make_option(c("-s", "--samplemeta-file"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) RDS/RData/xlsx/csv file containing a table of sample metadata. Any existing rownames will be replaced with the values in the sample ID  column (see below)."),
        make_option(c("-c", "--sample-id-column"), type="character", default="Sample",
                    help="Sample metadata column name that holds the sample IDs. These will be substituted into '--bam-file-pattern' to determine the BAM file names."),
        make_option(c("-p", "--bam-file-pattern"), metavar="PATTERN", type="character",
                    help="(REQUIRED) Format string to convert sample IDs into BAM file paths. This should contain a '%s' wherever the sample ID should be substituted ('%s' can occur multiple times),. Example: 'bam_files/Sample_%s/Aligned.bam"),
        make_option(c("-o", "--output-file"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) Output file name. The SummarizedExperiment object containing the counts will be saved here using saveRDS, so it should end in '.RDS'."),
        make_option(c("-b", "--expected-bam-files"), metavar="BAMFILE1,BAMFILE2,...", type="character",
                    help="Comma-separated list of bam file names expected to be used as input. This argument is optional, but if it is provided, it will be checked against the list of files determined from '--samplemeta-file' and '--bam-file-pattern', and an error will be raised if they don't match exactly."),
        make_option(c("-w", "--window-width"), type="character", default="150bp",
                    help="Width of windows in which to count."),
        make_option(c("-s", "--window-spacing"), metavar="BP", type="character",
                    help="Spacing between the start points of consecutive windows. By default, this is identical to the window width, so that the windows exactly tile the genome. Changing this results in either gapped windows (spacing > width) or overlapping windows (spacing < width)."),
        make_option(c("-e", "--read-extension"), type="character", default="100bp",
                    help="Assumed fragment length of reads. Each read will be assumed to represent a DNA fragment extending this far from its 5 prime end, regardless of the actual read length."),
        make_option(c("--bin"), action="store_true", default=FALSE,
                    help="Run in bin mode, where each read is counted into exactly one bin."),
        make_option(c("-j", "--threads"), metavar="N", type="integer", default=1,
                    help="Number of threads to use"))
    progname <- na.omit(c(get_Rscript_filename(), "csaw-count-windows.R"))[1]
    parser <- OptionParser(
        usage="Usage: %prog [options] -s SAMPLEMETA.RDS -p PATTERN -w WSIZE -e READEXT [ -s WSPACE ] -o SUMEXP.RDS",
        description="Do window counting across the genome for ChIP-Seq data",
option_list = optlist,
add_help_option = TRUE,
prog = progname,
epilogue = "Note that all base pair sizes (window width/spacing and read extension) may have an suffix of 'bp', 'kbp', 'mbp', or 'tbp'. For example, 10kb = 10000")

    cmdopts <- parse_args(parser, opts)
    ## Ensure that all required arguments were provided
    required.opts <- c("samplemeta-file", "output-file", "bam-file-pattern")
    missing.opts <- setdiff(required.opts, names(cmdopts))
    if (length(missing.opts) > 0) {
        stop(str_c("Missing required arguments: ", deparse(missing.opts)))
    }

    if (! "window-spacing" %in% names(cmdopts)) {
        cmdopts[["window-spacing"]] <- cmdopts[["window-width"]]
    }
    ## Convert bp args to numbers
    for (i in c("window-width", "window-spacing", "read-extension")) {
        cmdopts[[i]] %<>% parse.bp
    }

    cmdopts$threads %<>% round
    assert_that(cmdopts$threads >= 1)

    cmdopts$help <- NULL

    ## Replace dashes with underscores so that all options can easily
    ## be accessed by "$"
    cmdopts %>% setNames(str_replace_all(names(.), "-", "_"))
}

windowCountsParallel <- function(bam.files, ..., filter=10, BPPARAM=bpparam()) {
    reslist <- bplapply(X=bam.files, FUN=windowCounts, ..., filter=0, BPPARAM=BPPARAM)
    assert_that(all(sapply(reslist, is, "SummarizedExperiment")))
    res <- do.call(cbind, reslist)
    rm(reslist)
    keep <- rowSums(assay(res)) >= filter
    res[keep,]
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

    tsmsg("Using ", cmdopts$threads, " cores.")
    registerDoParallel(cores=cmdopts$threads)
    register(DoparParam())

    if (cmdopts$window_width == cmdopts$window_spacing) {
        tsmsg("Using a window size and spacing of ", format.bp(cmdopts$window_width), ".")
    } else {
        tsmsg("Using a window size of ", format.bp(cmdopts$window_width),
              " and a spacing of ", format.bp(cmdopts$window_spacing), ".")
    }
    ## Fragment size is not used for bins
    if (!cmdopts$bin) {
        tsmsg("Assuming a fragment size of ", format.bp(cmdopts$read_extension))
    }

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

    tsmsg("Loading blacklist regions")
    blacklist <- import("saved_data/ChIPSeq-merged-blacklist.bed", format="bed")

    ## Standard nuclear chromosomes only. (chrM is excluded because it is
    ## not located in the nucleus and is thus not subject to histone
    ## modification. The unplaced scaffolds are mostly not large enough to
    ## contain even a single typically-sized peak, so little is lost by
    ## excluding them for this analysis.)
    std.chr <- extractSeqlevels("Homo sapiens", "UCSC") %>% setdiff("chrM")
    rparam <- readParam(restrict=std.chr, discard=blacklist)

    tsmsg(sprintf("Counting reads in %s %s in %i samples.",
                  format.bp(cmdopts$window_width), ifelse(cmdopts$bin, "bins", "windows"),
                  nrow(sample.table)))
    if (cmdopts$threads > 1) {
        wCountsFun <- windowCountsParallel
        options(mc.cores=cmdopts$threads, mc.preschedule=FALSE)
        registerDoParallel(cores=cmdopts$threads)
        register(DoparParam())
    } else {
        wCountsFun <- windowCounts
    }
    wcounts <- wCountsFun(
        sample.table$bam_file, spacing=cmdopts$window_spacing,
        width=cmdopts$window_width, ext=cmdopts$read_extension,
        param=rparam, bin=cmdopts$bin)
    colData(wcounts) %<>% {cbind(sample.table, .[c("totals", "ext")])}
    saveRDS(wcounts, cmdopts$output_file)
}
