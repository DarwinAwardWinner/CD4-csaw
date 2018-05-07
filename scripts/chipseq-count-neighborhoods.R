#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(getopt)
    library(optparse)
    library(stringr)
    library(magrittr)
    library(assertthat)
    library(rctutils)
})

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
        make_option(c("-t", "--targets"), metavar = "FILENAME.RDS", type = "character",
                    help = "(REQUIRED) File specifying target genomic positions around which reads should be counted. This can be a BED file, GFF file, narrowPeak file, R data file containing a GRanges object, or csv file that can be converted to a GRanges object. If the regions have associated annotations, then a GRanges in an R data file is the recommended format. Generally the ranges specified should each be only a single base pair, which will be used as the center of the neighborhood. If any ranges are longer than 1bp, the neighborhoods will be formed around the 5-prime ends (or, for ranges with no strand information, the end closest to the beginning of the chromosome)."),
        make_option(c("-u", "--upstream-neighborhood"), type = "character", default = "5kbp", metavar = "5kbp",
                    help = "How far upstream to extend the neighborhood around each target."),
        make_option(c("-d", "--downstream-neighborhood"), type = "character", default = "5kbp", metavar = "5kbp",
                    help = "How far downstream to extend the neighborhood around each target."),
        make_option(c("-w", "--window-width"), type = "character", default = "1kbp", metavar = "1kbp",
                    help = "Width of windows in which to count."),
        make_option(c("--window-spacing"), metavar = "BP", type = "character",
                    help = "Spacing between the start points of consecutive windows. By default, this is identical to the window width, so that the windows exactly tile the neighborhood. Changing this results in either gapped windows (spacing > width) or overlapping windows (spacing < width)."),
        make_option(c("--initial-window-offset"), type = "character", default = "0bp", metavar = "0bp",
                    help = "Offset of each neighborhood's central window relative to each target. The default of 0 means that the central window will be centered directly around the target itself. Negative values place the center of the window further upstream (in the 5-prime direction), while positive values place it downstream (in the 3-prime direction). If you want a border between adjacent windows to fall on the target, set this to half the window width."),
        make_option(c("-o", "--output-file"), metavar = "FILENAME.RDS", type = "character",
                    help = "(REQUIRED) Output file name. The SummarizedExperiment object containing the counts will be saved here using saveRDS, so it should end in '.RDS'."),
        make_option(c("-b", "--expected-bam-files"), metavar = "BAMFILE1,BAMFILE2,...", type = "character",
                    help = "Comma-separated list of bam file names expected to be used as input. This argument is optional, but if it is provided, it will be checked against the list of files determined from '--samplemeta-file' and '--bam-file-pattern', and an error will be raised if they don't match exactly."),
        make_option(c("-e", "--read-extension"), type = "character", default = "100bp", metavar = "100bp",
                    help = "Assumed fragment length of unpaired reads. Each single read will be assumed to represent a DNA fragment extending this far from its 5 prime end, regardless of the actual read length. (Mated reads pairs already define a fragment length and are unaffected by this option.) Note that each read is counted into the window that contains the midpoint of the represented fragment."),
        make_option(c("-x", "--blacklist"), metavar = "FILENAME.bed", type = "character",
                    help = "File describing blacklist regions to be excluded from the analysis. Windows that overlap these regions will have their counts set to NA. This can be a BED file, GFF file, R data file containing a GRanges object, or csv file that can be converted to a GRanges object."),
        make_option(c("-a", "--blacklist-action"), metavar = "ACTION", type = "character", default = "mark",
                    help = "What action to take on windows that overlap blacklisted regions. Options are 'mark', 'setNA', and 'discard'. The default, 'mark', adds an additional logical column to the rowData of the output named 'blacklist' that is TRUE for windows overlapping the blacklist and FALSE for the rest. 'setNA' additionally sets the read count for blacklisted windows to NA. 'discard' throws away any blacklisted windows, so that they will not be present at all in the output."),
        make_option(c("-j", "--threads"), metavar = "N", type = "integer", default = 1,
                    help = "Number of threads to use"))
    progname <- na.omit(c(get_Rscript_filename(), "csaw-count-neighborhoods.R"))[1]
    parser <- OptionParser(
        usage = "Usage: %prog [options] -s SAMPLEMETA.RDS -p PATTERN -t TARGETS.RDS -o SUMEXP.RDS",
        description = "Count ChIP-seq reads in neighborhoods around a set of specified genomic positions.",
        option_list = optlist,
        add_help_option = TRUE,
        prog = progname,
        epilogue = "A \"neighborhood\" consists of a set of windows at regular offsets relative to a genomic position. For example, windows every 200 bp from 5kbp upstream to 2kb downstream. Each window is annotated with an 'offset' column that indicates the distance from the center of that window to the specified genomic position. Note that the windows at the edges of the neighborhood will likely extend past the specified distances upstream and downstream, since these distances refer to the center of the window, not the outer edge. In addition, note that not every neighborhood is guaranteed to contain a complete \"set\" of windows, since it may extend off the edge of a chromosome,

Note that all base pair sizes (window width/spacing and read extension) may have a suffix of 'bp', 'kbp', 'mbp', or 'tbp'. For example, 10kbp = 10000.")

    cmdopts <- parse_args(parser, opts)
    ## Ensure that all required arguments were provided
    required.opts <- c("samplemeta-file", "output-file", "bam-file-pattern", "targets")
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
    for (i in c("upstream-neighborhood", "downstream-neighborhood", "window-width", "window-spacing", "initial-window-offset", "read-extension")) {
        if (i %in% names(cmdopts)) {
            cmdopts[[i]] %<>% parse_bp
        }
    }

    if (is.null(cmdopts[["window-spacing"]])) {
        cmdopts[["window-spacing"]] <- cmdopts[["window-width"]]
    }

    ## Validate blacklist action
    if ("blacklist-action" %in% names(cmdopts)) {
        cmdopts[["blacklist-action"]] %<>%
            match_arg(c("mark", "setNA", "discard"), arg_name = "--blacklist-action", ignore.case = TRUE)
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

suppressPackageStartupMessages({
    library(dplyr)
    library(glue)
    library(future)
    library(GenomicRanges)
    library(SummarizedExperiment)
    library(GenomicAlignments)
    library(rlang)
    library(forcats)
    library(csaw)
})

options(future.globals.maxSize = 4 * 1024^3)

## cmdopts <- list(
##     samplemeta_file = "saved_data/samplemeta-ChIPSeq.RDS",
##     sample_id_column = "SRA_run",
##     bam_file_pattern = "aligned/chipseq_bowtie2_hg38.analysisSet/{SAMPLE}/Aligned.bam",
##     targets = "saved_data/tss_shoal_hg38.analysisSet_ensembl.85.RDS",
##     upstream_neighborhood = 5000,
##     downstream_neighborhood = 5000,
##     window_width = 500,
##     initial_window_offset = 0,
##     output_file = "saved_data/tss-neighborhood-counts_hg38.analysisSet_ensembl.85_5kbp-radius_500bp-windows_147bp-reads.RDS",
##     read_extension = 147,
##     blacklist = "saved_data/ChIPSeq-merged-blacklist.bed",
##     blacklist_action = "mark",
##     threads = 4,
##     window_spacing = 500)

{
    cmdopts <- get_options(commandArgs(TRUE))
    ## TODO: Don't use setwd in this or any other script
    tryCatch(setwd(file.path(dirname(na.omit(get_Rscript_filename())), "..")),
             error = function(...) tsmsg("WARNING: Could not determine script path. Ensure that you are already in the correct directory."))

    tsmsg("Args:")
    print_var_vector(cmdopts)

    if (cmdopts$threads > 1) {
        use_futures("multicore", workers = cmdopts$threads, quiet = TRUE)
    } else {
        use_futures("sequential", quiet = TRUE)
        register(SerialParam())
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

    tsmsg("Loading target positions")
    targets <- read_regions(cmdopts$targets)
    assert_that(is(targets, "GRanges"))
    if (any(strand(targets) == "*")) {
        warning("Some targets have no strand information, and will be treated as being on the plus strand.")
    }
    ## Reduce targets to 5-prime end only
    targets %<>% resize(width = 1, fix = "start")

    blacklist_regions <- GRanges()
    if (!is.null(cmdopts$blacklist)) {
        tsmsg("Loading blacklist regions")
        blacklist_regions <- read_regions(cmdopts$blacklist)
        assert_that(is(blacklist_regions, "GRanges"))
        ## Blacklist applies to both strands
        strand(blacklist_regions) <- "*"
    }

    nhood_offsets <- cmdopts %$% c(
        rev(seq(from = initial_window_offset,
                to = -upstream_neighborhood,
                by = -window_spacing)),
        seq(from = initial_window_offset + window_spacing,
            to = downstream_neighborhood,
            by = window_spacing))

    assert_that(is_named(targets))
    nhood_windows <- rep(targets, each = length(nhood_offsets))
    nhood_windows$offset <- rep(nhood_offsets, length.out = length(nhood_windows))
    nhood_windows %<>%
        shift(.$offset * ifelse(strand(.) == "-", -1, 1)) %>%
        resize(width = cmdopts$window_width, fix = "center")
    ## Add offset to window names
    names(nhood_windows) %<>% str_c(sprintf("%+ibp", nhood_windows$offset))

    if (length(blacklist_regions) > 0) {
        blacklisted <- overlapsAny(nhood_windows, blacklist_regions, ignore.strand = TRUE)
        nhood_windows$blacklist <- blacklisted
        if (cmdopts$blacklist_action == "discard") {
            tsmsg("Discarding blacklisted windows")
            nhood_windows %<>% .[!.$blacklisted]
        } else if (cmdopts$blacklist_action %in% c("mark", "setNA")) {
            ## Makring is handled above, and setNA will be handled
            ## later
            NULL
        } else {
            stop(glue("Unknown blacklist action '{cmdopts$blacklist_action}'"))
        }
    } else {
        ## We always add the blacklist column, but without a blacklist
        ## it is simply false for everything
        nhood_windows$blacklist <- FALSE
    }

    tsmsg(glue("Counting reads in neighborhoods around {length(targets)} regions in {nrow(sample_table)} samples."))
    tsmsg(glue("Neighborhoods consist of {length(nhood_offsets)} windows of width {format_bp(cmdopts$window_width)} tiled from {format_bp(cmdopts$upstream_neighborhood)} upstream to {format_bp(cmdopts$downstream_neighborhood)} downstream."))

    param <- readParam(BPPARAM=bpparam())
    rcounts <- regionCounts(
        sample_table$bam_file, regions=nhood_windows,
        # See ?windowCounts for explanation of "ext=list(...)"
        ext=list(rep(cmdopts$read_extension, nrow(sample_table)), 1),
        param = param)

    ## Add sample metadata to colData in front of mapping stats
    colData(rcounts) %<>% cbind(sample_table, .)
    colnames(rcounts) <- sample_table[[cmdopts$sample_id_column]]

    ## Set blacklisted window counts to NA, if requested
    if (cmdopts$blacklist_action == "setNA") {
        bl <- rowRanges(rcounts)$blacklist
        assay(rcounts, "counts")[bl,] <- NA
    }

    ## Save command and options in the metadata
    metadata(rcounts)$cmd.name <- na.omit(c(get_Rscript_filename(), "csaw-count-neighborhoods.R"))[1]
    metadata(rcounts)$cmd.opts <- cmdopts

    tsmsg("Saving output file")
    saveRDS(rcounts, cmdopts$output_file)
    tsmsg("Finished.")
}
