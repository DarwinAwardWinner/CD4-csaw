#!/usr/bin/env Rscript

library(getopt)
library(optparse)
library(stringr)
library(magrittr)
library(assertthat)
library(rctutils)

get_options <- function(opts) {
    optlist <- list(
        make_option(c("-s", "--samplemeta-file"), metavar = "FILENAME.RDS", type = "character",
                    help = "(REQUIRED) RDS/RData/xlsx/csv file containing a table of sample metadata. Any existing rownames will be replaced with the values in the sample ID  column (see below)."),
        make_option(c("-c", "--sample-id-column"), type = "character", default = "Sample",
                    help = "Sample metadata column name that holds the sample IDs. These will be substituted into '--bam-file-pattern' to determine the BAM file names."),
        make_option(c("-f", "--filter-sample-ids"), type = "character",
                    help = "Comma-separated list of sample IDs. If this options is provided, only the specified sample IDs will be used."),
        make_option(c("-p", "--bam-file-pattern"), metavar = "PATTERN", type = "character",
                    help = "(REQUIRED) Format string to convert sample IDs into BAM file paths. This should contain the string '{SAMPLE}' wherever the sample ID should be substituted (this can occur multiple times),. Example: 'bam_files/Sample_{SAMPLE}/Aligned.bam"),
        make_option(c("-r", "--regions"), metavar = "FILENAME.RDS", type = "character",
                    help = "(REQUIRED) File specifying regions in which reads should be counted. This can be a BED file, GFF file, narrowPeak file, R data file containing a GRanges object, or csv file that can be converted to a GRanges object. If the regions have associated annotations, then a GRanges in an R data file is the recommended format."),
        make_option(c("-o", "--output-file"), metavar = "FILENAME.RDS", type = "character",
                    help = "(REQUIRED) Output file name. The SummarizedExperiment object containing the counts will be saved here using saveRDS, so it should end in '.RDS'."),
        make_option(c("-b", "--expected-bam-files"), metavar = "BAMFILE1,BAMFILE2,...", type = "character",
                    help = "Comma-separated list of bam file names expected to be used as input. This argument is optional, but if it is provided, it will be checked against the list of files determined from '--samplemeta-file' and '--bam-file-pattern', and an error will be raised if they don't match exactly."),
        make_option(c("-e", "--read-extension"), type = "character", default = "100bp",
                    help = "Assumed fragment length of reads. Each read will be assumed to represent a DNA fragment extending this far from its 5 prime end, regardless of the actual read length."),
        make_option(c("-x", "--blacklist"), metavar = "FILENAME.bed", type = "character",
                    help = "File describing blacklist regions to be excluded from the analysis. Reads that overlap these regions will be discarded without counting them toward any region. This can be a BED file, GFF file, R data file containing a GRanges object, or csv file that can be converted to a GRanges object."),
        make_option(c("-j", "--threads"), metavar = "N", type = "integer", default = 1,
                    help = "Number of threads to use"))
    progname <- na.omit(c(get_Rscript_filename(), "chipseq-count-windows.R"))[1]
    parser <- OptionParser(
        usage = "Usage: %prog [options] -s SAMPLEMETA.RDS -p PATTERN -r REGIONS.RDS -o SUMEXP.RDS",
        description = "Count ChIP-seq reads a set of specified regions",
        option_list = optlist,
        add_help_option = TRUE,
        prog = progname,
        epilogue = "Note that all base pair sizes (window width/spacing and read extension) may have a suffix of 'bp', 'kbp', 'mbp', or 'tbp'. For example, 10kbp = 10000.")

    cmdopts <- parse_args(parser, opts)
    ## Ensure that all required arguments were provided
    required.opts <- c("samplemeta-file", "output-file", "bam-file-pattern", "regions")
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

    ## Convert bp args to numbers
    for (i in c("read-extension")) {
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
invisible(get.options(commandArgs(TRUE)))

library(dplyr)
library(glue)
library(future)
library(GenomicRanges)
library(SummarizedExperiment)
library(csaw)
library(forcats)

## cmdopts <- list(
##     samplemeta_file = "saved_data/samplemeta-ChIPSeq.RDS",
##     sample_id_column = "SRA_run",
##     bam_file_pattern = "aligned/chipseq_bowtie2_hg38.analysisSet/{SAMPLE}/Aligned.bam",
##     regions = "saved_data/promoter-regions_hg38.analysisSet_knownGene_2.5kbp.RDS",
##     output_file = "saved_data/promoter-counts_hg38.analysisSet_knownGene_2.5kbp-radius_147bp-reads.RDS",
##     read_extension = 147,
##     blacklist = "saved_data/ChIPSeq-merged-blacklist.bed",
##     threads = 4)

{
    cmdopts <- get.options(commandArgs(TRUE))
    tryCatch(setwd(file.path(dirname(na.omit(get_Rscript_filename())), "..")),
             error = function(...) tsmsg("WARNING: Could not determine script path. Ensure that you are already in the correct directory."))

    tsmsg("Args:")
    print_var_vector(cmdopts)

    if (cmdopts$threads > 1) {
        use_futures("multicore", workers = cmdopts$threads, quiet = TRUE)
    } else {
        use_futures("sequential", quiet = TRUE)
    }
    tsmsg(glue("Using {cmdopts$threads} cores."))

    tsmsg(glue("Assuming a fragment size of {format_bp(cmdopts$read_extension)} for unpaired reads."))

    tsmsg("Loading sample data")

    sample_table <- readRDS(cmdopts$samplemeta_file) %>%
        ## Compute full path to BAM file
        mutate(bam_file = glue(cmdopts$bam_file_pattern, SAMPLE = .[[cmdopts$sample_id_column]])) %>%
        ## Ensure that days_after_activation is a factor and can't be
        ## interpreted as a numeric
        mutate(days_after_activation = days_after_activation %>%
                   factor %>% fct_relabel(~str_c("Day", .))) %>%
        rename(time_point = days_after_activation)

    if (!is.null(cmdopts$filter_sample_ids)) {
        tsmsg("Selecting only ", length(cmdopts$filter_sample_ids), " specified samples.")
        assert_that(all(cmdopts$filter_sample_ids %in% sample_table[[cmdopts$sample_id_column]]))
        sample_table %<>% .[.[[cmdopts$sample_id_column]] %in% cmdopts$filter_sample_ids,]
    }

    assert_that(all(file.exists(sample_table$bam_file)))

    if ("expected_bam_files" %in% names(cmdopts)) {
        tryCatch({
            assert_that(setequal(samplemeta$bam_file, cmdopts$expected_bam_files))
            tsmsg("Sample metadata contains all expected bam files")
        }, error = function(...) {
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

    tsmsg("Loading regions")
    target_regions <- read_regions(cmdopts$regions)
    assert_that(is(target_regions, "GRanges"))
    ## Analysis is not stranded
    strand(target_regions) <- "*"

    blacklist_regions <- GRanges()
    if (!is.null(cmdopts$blacklist)) {
        tsmsg("Loading blacklist regions")
        blacklist_regions <- read_regions(cmdopts$blacklist)
        assert_that(is(blacklist_regions, "GRanges"))
        ## Blacklist applies to both strands
        strand(blacklist_regions) <- "*"
    }
    rparam <- readParam(discard = blacklist_regions)

    tsmsg(glue("Counting reads in {length(target_regions)} regions in {nrow(sample_table)} samples."))
    if (cmdopts$threads > 1) {
        rCountsFun <- regionCountsParallel
    } else {
        rCountsFun <- regionCounts
    }
    rcounts <- rCountsFun(
        sample_table$bam_file, regions=target_regions,
        ext=cmdopts$read_extension, param=rparam)

    ## Add sample metadata to colData in front of mapping stats
    colData(rcounts) %<>% cbind(sample_table, .)
    colnames(rcounts) <- sample_table[[cmdopts$sample_id_column]]

    ## Save command and options in the metadata
    metadata(sexp)$cmd_name <- na.omit(c(get_Rscript_filename(), "chipseq-count-regions.R"))[1]
    metadata(sexp)$cmd_opts <- cmdopts

    tsmsg("Saving output file")
    saveRDS(rcounts, cmdopts$output_file)
    tsmsg("Finished.")
}
