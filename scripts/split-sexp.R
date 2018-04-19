#!/usr/bin/env Rscript

library(getopt)
library(glue)
library(optparse)
library(stringr)
library(assertthat)
library(magrittr)
library(SummarizedExperiment)
library(dplyr)

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

get.options <- function(opts) {
    optlist <- list(
        make_option(c("-i", "--input-file"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) Input file name. This should be an RDS file containing a SummarizedExperiment object, whose containing the counts will be saved here using saveRDS, so it should end in '.RDS'."),
        make_option(c("-o", "--output-file-pattern"), metavar="TEMPLATE.RDS", type="character",
                    help="(REQUIRED) Output file name pattern. This should contain one or more column names from the sample metadata enclosed in curly braces. For example: 'csaw-counts-{chip_antibody}.RDS'. These will be filled in for each sample based on that sample's metadata, and the samples will be split into those files accordingly. Each SummarizedExperiment object will be saved using saveRDS, so it should end in '.RDS'."))
    progname <- na.omit(c(get_Rscript_filename(), "sexp-split.R"))[1]
    parser <- OptionParser(
        usage="Usage: %prog -i INFILE.RDS -o OUTTEMPLATE.RDS",
        description="Split SummarizedExperiment file by metadata",
        option_list = optlist,
        add_help_option = TRUE,
        prog = progname)

    cmdopts <- parse_args(parser, opts)
    cmdopts$help <- NULL
    ## Ensure that all required arguments were provided
    required.opts <- c("input-file", "output-file-pattern")
    missing.opts <- setdiff(required.opts, names(cmdopts))
    if (length(missing.opts) > 0) {
        stop(str_c("Missing required arguments: ", deparse(missing.opts)))
    }
    assert_that(str_detect(cmdopts[['output-file-pattern']], c("\\{.*\\}")))
    ## Replace dashes with underscores so that all options can easily
    ## be accessed by "$"
    cmdopts %>% setNames(chartr("-", "_", names(.)))
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

    tsmsg("Reading SummarizedExperiment file")
    sexp <- readRDS(cmdopts$input_file)
    assert_that(is(sexp, "SummarizedExperiment"))

    output_filenames = glue_data(as.list(colData(sexp)), cmdopts$output_file_pattern)
    output_groups <- split(seq_len(ncol(sexp)), output_filenames)
    output_sexps <- lapply(output_groups, . %>% sexp[,.])
    for (fname in names(output_sexps)) {
        tsmsg("Writing ", fname)
        saveRDS(output_sexps[[fname]], fname)
    }
}
