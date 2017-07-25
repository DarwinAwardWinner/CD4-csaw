#!/usr/bin/env Rscript

# This is a script to reproduce the plots of
# https://github.com/nboley/idr using ggplot2.

library(getopt)
library(optparse)

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

get.options <- function(opts) {

    ## Do argument parsing early so the script exits quickly if arguments are invalid
    optlist <- list(
        make_option(c("-p", "--peak-file"), metavar="FILE.narrowPeak", type="character",
                    help="(REQUIRED) NarrowPeak file containing peaks called from the combination of all samples."),
        make_option(c("-o", "--output-file"), metavar="FILE.narrowPeak", type="character",
                    help="(REQUIRED) Output file in which to save filtered peaks. Output peaks will always be sorted."),
        make_option(c("-i", "--idr-files"), metavar="FILE1[,FILE2,...]", type="character",
                    help="(REQUIRED) Comma-separated list of output files from idr script containing computed IDR values. Typically these represent every possible pairing of individual samples."),
        make_option(c("-t", "--idr-threshold"), type="numeric", default=1,
                    help="IDR threshold at which to filter peaks. A threshold of 1 will perform no filtering (useful with '-r')."),
        make_option(c("-r", "--replace-qValue-with-min-IDR-theshold"), action="store_true", default=FALSE,
                    help="In the output, replace the qValue column from the input with the minimum IDR threshold at which the peak would be kept. In keeping with the convention of the narrowPeak file format, the negative base 10 logarithm of the threshold value will be stored in the output file."),
        make_option(c("-s", "--sort-by"), type="character", default="pValue,score,signalValue",
                    help="Comma-separated list of columns to sort by. Later columns will be used to resolve any ties in earlier columns. Available columns are: pValue, qValue, score, signalValue."))
    progname <- na.omit(c(get_Rscript_filename(), "filter-by-idr.R"))[1]
    parser <- OptionParser(
        usage="Usage: %prog [options] -p combined_peaks.narrowPeak -o combined_filtered_peaks.narrowPeak -i idrfile1[,idrfile2,...]",
        description="Filter a peak list by IDR and/or annotate the peaks with computed IDR thresholds.",
        option_list = optlist,
        add_help_option = TRUE,
        prog = progname,
        epilogue = "Note that at least one of '-t' (with a value less than 1) or '-r' is required, since otherwise this script would not perform any useful action.")

    cmdopts <- parse_args(parser, opts)
    ## Ensure that all required arguments were provided
    required.opts <- c("peak-file", "output-file", "idr-files")
    missing.opts <- setdiff(required.opts, names(cmdopts))
    if (length(missing.opts) > 0) {
        stop(str_c("Missing required arguments: ", deparse(missing.opts)))
    }
    if (cmdopts[['idr-threshold']] >= 1 && !cmdopts[['replace-qValue-with-min-IDR-theshold']]) {
        stop("You must either use the '-r' option or use the '-t' option with a value less than 1.")
    }
    cmdopts %>% setNames(str_replace_all(names(.), "-", "_"))
}

## Do this early to handle "--help" before wasting time loading
## pacakges & stuff
invisible(get.options(commandArgs(TRUE)))

library(magrittr)
library(dplyr)
library(ggplot2)
library(scales)
library(ks)
library(reshape2)
library(stringr)
library(stringi)
library(rex)
library(Biobase)

print.var.vector <- function(v) {
    for (i in names(v)) {
        cat(i, ": ", deparse(v[[i]]), "\n", sep="")
    }
    invisible(v)
}

read.idr.table <- function(file) {
    idrcols <- c("chr", "start", "end", "name", "score", "strand",
                 "LocalIDR", "GlobalIDR", "startA", "endA", "scoreA", "startB", "endB", "scoreB")
    read.table(file, header=FALSE, sep="\t", col.names=idrcols) %>%
        mutate(LocalIDR=10^-LocalIDR, GlobalIDR=10^-GlobalIDR)
}

read.narrowPeak <- function(file, ...) {
    peaks.df <- read.table(file, sep="\t", row.names=NULL, ...)
    names(peaks.df) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")
    peaks.df$name <- as.character(peaks.df$name)
    ## havenames <- !any(peaks.df$name == ".")
    ## res <- data.frame2GRanges(peaks.df, keepColumns=TRUE, startOffset=1, endOffset=0)
    ## ## Eliminate the dummy row names from the data
    ## if (havenames)
    ##     names(res) <- res$name
    ## else
    ##     names(res) <- NULL
    ## res
    peaks.df
}

write.narrowPeak <- function(x, file, ...) {
    x <- as(x, "data.frame")
    if("seqnames" %in% names(x))
        names(x)[names(x) == "seqnames"] <- "chr"
    x <- x[c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")]
    write.table(x, file, sep="\t", row.names=FALSE, col.names=FALSE, ...)
}

{
    cmdopts <- get.options(commandArgs(TRUE))
    ## myargs <- c("-p", "peak_calls/epic_hg38.analysisSet/H3K4me3_condition.ALL_donor.ALL/peaks_noBL.narrowPeak", "-i", "idr_analysis/epic_hg38.analysisSet/H3K4me3_condition.ALL_D4659vsD5053/idrValues.txt,idr_analysis/epic_hg38.analysisSet/H3K4me3_condition.ALL_D4659vsD5131/idrValues.txt,idr_analysis/epic_hg38.analysisSet/H3K4me3_condition.ALL_D4659vsD5291/idrValues.txt,idr_analysis/epic_hg38.analysisSet/H3K4me3_condition.ALL_D5053vsD5131/idrValues.txt,idr_analysis/epic_hg38.analysisSet/H3K4me3_condition.ALL_D5053vsD5291/idrValues.txt,idr_analysis/epic_hg38.analysisSet/H3K4me3_condition.ALL_D5131vsD5291/idrValues.txt",
    ##             "-o", "peak_calls/epic_hg38.analysisSet/H3K4me3_condition.ALL_donor.ALL/peaks_noBL_IDR.narrowPeak",
    ##             "-t", "0.05", "-r")
    ## cmdopts <- get.options(myargs)
    cmdopts$help <- NULL

    tsmsg("Args:")
    print.var.vector(cmdopts)

    tsmsg("Reading peaks")
    peaks <- read.narrowPeak(cmdopts$peak_file)

    tsmsg("Sorting peaks")
    sort.cols <- str_split(cmdopts$sort_by, ",")[[1]]
    sort.args <- sapply(sort.cols, . %>% as.name %>% interp(quote(desc(x)), x=.))
    peaks %<>% arrange_(.dots=sort.args)

    tsmsg("Reading IDR files")
    idr.files <- str_split(cmdopts$idr_files, ",")[[1]]
    idr.tables <- lapply(idr.files, read.idr.table)

    tsmsg("computing minimum IDR thresholds")
    idr.thresh <- rep(1, nrow(peaks))
    idr.mat <- idr.tables %>% lapply(. %$% GlobalIDR %>% sort %>% rbind) %>%
        do.call(what=rbind.fill.matrix) %>% t
    idr.mat[is.na(idr.mat)] <- 1
    n <- min(length(idr.thresh), nrow(idr.mat))
    idr.thresh[seq_len(n)] <- rowMin(idr.mat)

    if (cmdopts$replace_qValue_with_min_IDR_theshold) {
        tsmsg("Replacing qValue column with minimum IDR threshold")
        peaks$qValue <- -log10(idr.thresh)
    }

    if (cmdopts$idr_threshold < 1) {
        tsmsg("Filtering peaks at an IDR threshold of ", cmdopts$idr_threshold)
        nfilter <- sum(idr.thresh <= cmdopts$idr_threshold)
        tsmsg("Selecting the top ", nfilter, " peaks out of ", nrow(peaks), ".")
        peaks <- peaks[seq_len(nfilter),]
        tsmsg("Saving filtered peaks")
    } else if (cmdopts$replace_qValue_with_min_IDR_theshold) {
        tsmsg("Saving peaks")
    } else {
        stop("No action was requested.")
    }

    write.narrowPeak(peaks, cmdopts$output_file)
}
