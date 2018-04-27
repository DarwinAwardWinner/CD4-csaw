#!/usr/bin/env Rscript

suppressMessages({
    library(getopt)
    library(optparse)
    library(stringr)
    library(glue)
    library(rex)
    library(sitools)
    library(magrittr)
    library(GenomicRanges)
    library(rtracklayer)
    library(SummarizedExperiment)
    library(dplyr)
    library(csaw)
    library(Matrix)
    library(assertthat)
    library(rctutils)

    library(future)
    library(doParallel)
    library(BiocParallel)
})

get_options <- function(opts) {
    optlist <- list(
        make_option(c("-s", "--samplemeta-file"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) RDS/RData/xlsx/csv file containing a table of sample metadata. Any existing rownames will be replaced with the values in the sample ID  column (see below)."),
        make_option(c("-c", "--sample-id-column"), type="character", default="Sample",
                    help="Sample metadata column name that holds the sample IDs. These will be substituted into '--bam-file-pattern' to determine the BAM file names."),
        make_option(c("-f", "--filter-sample-ids"), type="character",
                    help="Comma-separated list of sample IDs. If this options is provided, only the specified sample IDs will be used."),
        make_option(c("-p", "--bam-file-pattern"), metavar="PATTERN", type="character",
                    help="(REQUIRED) Format string to convert sample IDs into BAM file paths. This should contain the string '{SAMPLE}' wherever the sample ID should be substituted (this can occur multiple times),. Example: 'bam_files/Sample_{SAMPLE}/Aligned.bam"),
        make_option(c("-o", "--output-file"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) Output file name. The SummarizedExperiment object containing the counts will be saved here using saveRDS, so it should end in '.RDS'."),
        make_option(c("-b", "--expected-bam-files"), metavar="BAMFILE1,BAMFILE2,...", type="character",
                    help="Comma-separated list of bam file names expected to be used as input. This argument is optional, but if it is provided, it will be checked against the list of files determined from '--samplemeta-file' and '--bam-file-pattern', and an error will be raised if they don't match exactly."),
        make_option(c("-w", "--window-width"), type="character", default="150bp",
                    help="Width of windows in which to count."),
        make_option(c("--window-spacing"), metavar="BP", type="character",
                    help="Spacing between the start points of consecutive windows. By default, this is identical to the window width, so that the windows exactly tile the genome. Changing this results in either gapped windows (spacing > width) or overlapping windows (spacing < width)."),
        make_option(c("-e", "--read-extension"), type="character", default="100bp",
                    help="Assumed fragment length of reads. Each read will be assumed to represent a DNA fragment extending this far from its 5 prime end, regardless of the actual read length."),
        make_option(c("-x", "--blacklist"), metavar="FILENAME.bed", type="character",
                    help="File describing blacklist regions to be excluded from the analysis. Reads that overlap these regions will be discarded without counting them toward any window. This can be a BED file, GFF file, R data file containing a GRanges object, or csv file that can be converted to a GRanges object."),
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
        epilogue = "Note that all base pair sizes (window width/spacing and read extension) may have a suffix of 'bp', 'kbp', 'mbp', or 'tbp'. For example, 10kbp = 10000")

    cmdopts <- parse_args(parser, opts)
    ## Ensure that all required arguments were provided
    required.opts <- c("samplemeta-file", "output-file", "bam-file-pattern")
    missing.opts <- setdiff(required.opts, names(cmdopts))
    if (length(missing.opts) > 0) {
        stop(str_c("Missing required arguments: ", deparse(missing.opts)))
    }

    ## Split list arguments
    for (i in c("filter-sample-ids", "expected-bam-files")) {
        if (i %in% names(cmdopts)) {
            cmdopts[[i]] %<>% str_split(",") %>% unlist
        }
    }

    if (! "window-spacing" %in% names(cmdopts)) {
        cmdopts[["window-spacing"]] <- cmdopts[["window-width"]]
    }
    ## Convert bp args to numbers
    for (i in c("window-width", "window-spacing", "read-extension")) {
        cmdopts[[i]] %<>% parse_bp
    }

    cmdopts$threads %<>% round
    assert_that(cmdopts$threads >= 1)

    cmdopts$help <- NULL

    ## Replace dashes with underscores so that all options can easily
    ## be accessed by "$"
    cmdopts %>% setNames(chartr("-", "_", names(.)))
}

## Terminate early on argument-processing errors
invisible(get_options(commandArgs(TRUE)))

{
    cmdopts <- get_options(commandArgs(TRUE))
    ## TODO: Eliminate all setwd
    tryCatch(setwd(file.path(dirname(na.omit(get_Rscript_filename())), "..")),
             error=function(...) tsmsg("WARNING: Could not determine script path. Ensure that you are already in the correct directory."))

    tsmsg("Args:")
    print_var_vector(cmdopts)

    if (cmdopts$threads > 1) {
        setup_multicore()
    } else {
        registerDoSEQ()
        register(SerialParam())
    }

    if (cmdopts$window_width == cmdopts$window_spacing) {
        tsmsg("Using a window size and spacing of ", format_bp(cmdopts$window_width), ".")
    } else {
        tsmsg("Using a window size of ", format_bp(cmdopts$window_width),
              " and a spacing of ", format_bp(cmdopts$window_spacing), ".")
    }
    ## Fragment size is not used for bins
    if (!cmdopts$bin) {
        tsmsg("Assuming a fragment size of ", format_bp(cmdopts$read_extension))
    }

    tsmsg("Loading sample data")

    sample.table <- readRDS(cmdopts$samplemeta_file) %>%
        ## Compute full path to BAM file
        mutate(bam_file=glue(cmdopts$bam_file_pattern, SAMPLE=.[[cmdopts$sample_id_column]], .envir=emptyenv())) %>%
        ## Ensure that days_after_activation is a factor and can't be
        ## interpreted as a numeric
        mutate(days_after_activation=days_after_activation %>%
                   factor %>% `levels<-`(str_c("Day", levels(.)))) %>%
        rename(time_point=days_after_activation)

    if (!is.null(cmdopts$filter_sample_ids)) {
        tsmsg("Selecting only ", length(cmdopts$filter_sample_ids), " specified samples.")
        assert_that(all(cmdopts$filter_sample_ids %in% sample.table[[cmdopts$sample_id_column]]))
        sample.table %<>% .[.[[cmdopts$sample_id_column]] %in% cmdopts$filter_sample_ids,]
    }

    assert_that(all(file.exists(sample.table$bam_file)))

    if ("expected_bam_files" %in% names(cmdopts)) {
        tryCatch({
            assert_that(setequal(samplemeta$bam_file, cmdopts$expected_bam_files))
            tsmsg("Sample metadata contains all expected bam files")
        }, error=function(...) {
            unexpected_existing <- setdiff(samplemeta$bam_file, cmdopts$expected_bam_files)
            expected_but_missing <- setdiff(cmdopts$expected_bam_files, samplemeta$bam_file)
            if (length(unexpected_existing) > 0) {
                tsmsg(glue("Got unexpected bam files: {deparse(unexpected_existing)}"))
            }
            if (length(expected_but_missing) > 0) {
                tsmsg(glue("Didn't find expected bam files: {deparse(expected_but_missing)}"))
            }
            stop("Bam file list was not as expected")
        })
    }

    blacklist_regions <- GRanges()
    if (!is.null(cmdopts$blacklist)) {
        tsmsg("Loading blacklist regions")
        blacklist_regions <- read_regions(cmdopts$blacklist)
        assert_that(is(blacklist_regions, "GRanges"))
        ## Blacklist applies to both strands
        strand(blacklist_regions) <- "*"
    }

    ## Standard nuclear chromosomes only. (chrM is excluded because it is
    ## not located in the nucleus and is thus not subject to histone
    ## modification. The unplaced scaffolds are mostly not large enough to
    ## contain even a single typically-sized peak, so little is lost by
    ## excluding them for this analysis.)
    std.chr <- extractSeqlevels("Homo sapiens", "UCSC") %>% setdiff("chrM")
    rparam <- readParam(discard=blacklist)

    tsmsg(glue("Counting reads in {width} {type} in {scount} samples.",
        width=format_bp(cmdopts$window_width),
        type=ifelse(cmdopts$bin, "bins", "windows"),
        scount=nrow(sample.table)))
    if (cmdopts$threads > 1) {
        wCountsFun <- windowCountsParallel
    } else {
        wCountsFun <- windowCounts
    }
    wcounts <- wCountsFun(
        sample.table$bam_file, spacing=cmdopts$window_spacing,
        width=cmdopts$window_width, ext=cmdopts$read_extension,
        param=rparam, bin=cmdopts$bin)
    colData(wcounts) %<>% {cbind(sample.table, .[c("totals", "ext")])}
    ## Save command and options in the metadata
    metadata(sexp)$cmd_name <- na.omit(c(get_Rscript_filename(), "csaw-count-windows.R"))[1]
    metadata(sexp)$cmd_opts <- cmdopts

    tsmsg("Saving output file")
    saveRDS(wcounts, cmdopts$output_file)
    tsmsg("Finished.")
}
