#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(getopt)
    library(optparse)
    library(stringr)
    library(rex)
    library(sitools)
    library(magrittr)
    library(assertthat)
})

tsmsg <- function(...) {
    message(date(), ": ", ...)
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

# Inverse of sitools::f2si
si2f <- function(string, unit="") {
    if (length(string) == 0) {
        return(numeric(0))
    }
    sifactor <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-09, 1e-06,
                  0.001, 1, 1000, 1e+06, 1e+09, 1e+12, 1e+15, 1e+18, 1e+21,
                  1e+24)
    pre <- c("y", "z", "a", "f", "p", "n", "u", "m",
             "", "k", "M", "G", "T", "P", "E", "Z", "Y")

    rx <- rex(
        ## Leading whitespace
        start,
        zero_or_more(space),

        ## Capture a floating point number
        capture(
            ## Sign
            maybe(one_of("+", "-")),
            ## Integer part
            zero_or_more(digit),
            ## Decimal point
            maybe("."),
            ## Fractional part (or integer part when decimal is not
            ## present)
            one_or_more(digit),
            ## Exponential notation
            maybe(
                one_of("e", "E"),
                maybe(one_of("+", "-")),
                one_or_more(digit)
            )
        ),

        ## Space between number and unit
        zero_or_more(space),

        ## Capture SI prefix
        capture(maybe(one_of(pre))),

        unit,

        ## Trailing whitespace
        zero_or_more(space),
        end
    )

    m <- str_match(string, rx)
    base <- as.numeric(m[,2])
    p <- m[,3]
    fac <- sifactor[match(p, pre)]
    base * fac
}

parse.bp <- function(size) {
    suppressWarnings({
        result <- si2f(size, "bp")
        ## Fall back to just parsing a number without the "bp" suffix
        result[is.na(result)] <- si2f(size[is.na(result)])
    })
    assert_that(!any(is.na(result)))
    result
}

format.bp <- function(x, use_si=TRUE, always_signed=FALSE) {
    result <- round(x)
    if (use_si) {
        result %<>% f2si("bp") %>% str_replace_all(rex(one_or_more(space)), "")
    } else {
        result %<>% str_c("bp")
    }
    if (always_signed) {
        result %<>% str_c(ifelse(str_detect(., "^-"), "", "+"), .)
    }
    result
}

get.options <- function(opts) {
    optlist <- list(
        make_option(c("-s", "--samplemeta-file"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) RDS/RData/xlsx/csv file containing a table of sample metadata. Any existing rownames will be replaced with the values in the sample ID  column (see below)."),
        make_option(c("-c", "--sample-id-column"), type="character", default="Sample",
                    help="Sample metadata column name that holds the sample IDs. These will be substituted into '--bam-file-pattern' to determine the BAM file names."),
        make_option(c("-f", "--filter-sample-ids"), type="character",
                    help="Comma-separated list of sample IDs. If this options is provided, only the specified sample IDs will be used."),
        make_option(c("-p", "--bam-file-pattern"), metavar="PATTERN", type="character",
                    help="(REQUIRED) Format string to convert sample IDs into BAM file paths. This should contain the string '{SAMPLE}' wherever the sample ID should be substituted (this can occur multiple times),. Example: 'bam_files/Sample_{SAMPLE}/Aligned.bam"),
        make_option(c("-t", "--targets"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) File specifying target genomic positions around which reads should be counted. This can be a BED file, GFF file, narrowPeak file, R data file containing a GRanges object, or csv file that can be converted to a GRanges object. If the regions have associated annotations, then a GRanges in an R data file is the recommended format. Generally the ranges specified should each be only a single base pair, which will be used as the center of the neighborhood. If any ranges are longer than 1bp, the neighborhoods will be formed around the 5-prime ends (or, for ranges with no strand information, the end closest to the beginning of the chromosome)."),
        make_option(c("-u", "--upstream-neighborhood"), type="character", default="5kbp", metavar="5kbp",
                    help="How far upstream to extend the neighborhood around each target."),
        make_option(c("-d", "--downstream-neighborhood"), type="character", default="5kbp", metavar="5kbp",
                    help="How far downstream to extend the neighborhood around each target."),
        make_option(c("-w", "--window-width"), type="character", default="1kbp", metavar="1kbp",
                    help="Width of windows in which to count."),
        make_option(c("--window-spacing"), metavar="BP", type="character",
                    help="Spacing between the start points of consecutive windows. By default, this is identical to the window width, so that the windows exactly tile the neighborhood. Changing this results in either gapped windows (spacing > width) or overlapping windows (spacing < width)."),
        make_option(c("--initial-window-offset"), type="character", default="0bp", metavar="0bp",
                    help="Offset of each neighborhood's central window relative to each target. The default of 0 means that the central window will be centered directly around the target itself. Negative values place the center of the window further upstream (in the 5-prime direction), while positive values place it downstream (in the 3-prime direction). If you want a border between adjacent windows to fall on the target, set this to half the window width."),
        make_option(c("-o", "--output-file"), metavar="FILENAME.RDS", type="character",
                    help="(REQUIRED) Output file name. The SummarizedExperiment object containing the counts will be saved here using saveRDS, so it should end in '.RDS'."),
        make_option(c("-b", "--expected-bam-files"), metavar="BAMFILE1,BAMFILE2,...", type="character",
                    help="Comma-separated list of bam file names expected to be used as input. This argument is optional, but if it is provided, it will be checked against the list of files determined from '--samplemeta-file' and '--bam-file-pattern', and an error will be raised if they don't match exactly."),
        make_option(c("-e", "--read-extension"), type="character", default="100bp", metavar="100bp",
                    help="Assumed fragment length of unpaired reads. Each single read will be assumed to represent a DNA fragment extending this far from its 5 prime end, regardless of the actual read length. (Mated reads pairs already define a fragment length and are unaffected by this option.) Note that each read is counted into the window that contains the midpoint of the represented fragment."),
        make_option(c("-x", "--blacklist"), metavar="FILENAME.bed", type="character",
                    help="File describing blacklist regions to be excluded from the analysis. Windows that overlap these regions will have their counts set to NA. This can be a BED file, GFF file, R data file containing a GRanges object, or csv file that can be converted to a GRanges object."),
        make_option(c("-a", "--blacklist-action"), metavar="ACTION", type="character", default="mark",
                    help="What action to take on windows that overlap blacklisted regions. Options are 'mark', 'setNA', and 'discard'. The default, 'mark', adds an additional logical column to the rowData of the output named 'blacklist' that is TRUE for windows overlapping the blacklist and FALSE for the rest. 'setNA' additionally sets the read count for blacklisted windows to NA. 'discard' throws away any blacklisted windows, so that they will not be present at all in the output."),
        make_option(c("-j", "--threads"), metavar="N", type="integer", default=1,
                    help="Number of threads to use"))
    progname <- na.omit(c(get_Rscript_filename(), "csaw-count-neighborhoods.R"))[1]
    parser <- OptionParser(
        usage="Usage: %prog [options] -s SAMPLEMETA.RDS -p PATTERN -t TARGETS.RDS -o SUMEXP.RDS",
        description="Count ChIP-seq reads in neighborhoods around a set of specified genomic positions.",
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
            cmdopts[[i]] %<>% parse.bp
        }
    }

    if (is.null(cmdopts[["window-spacing"]])) {
        cmdopts[["window-spacing"]] <- cmdopts[["window-width"]]
    }

    ## Validate blacklist action
    if ("blacklist-action" %in% names(cmdopts)) {
        cmdopts[["blacklist-action"]] %<>%
            match.arg(c("mark", "setNA", "discard"), argname="--blacklist-action", ignore.case=TRUE)
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

suppressPackageStartupMessages({
    library(dplyr)
    library(glue)
    library(Matrix)
    library(future)
    library(GenomicRanges)
    library(rtracklayer)
    library(SummarizedExperiment)
    library(GenomicAlignments)
    library(doParallel)
    library(BiocParallel)
    library(rlang)
})

## Should really be an S4 method, but writing S4 methods is a pain
readsToFragmentMidpoints <- function(reads, fraglen) {
    if (is(reads, "GAlignmentsList")) {
        isSingle <- lengths(reads) == 1
        mated.frags <- granges(reads[!isSingle], ignore.strand=TRUE)
        ## Extend smaller reads to fraglen, but don't shrink longer
        ## reads
        single.frags <- granges(reads[isSingle]) %>% resize(width=pmax(fraglen, width(.)), fix="start")
        frags <- c(mated.frags, single.frags)
    } else if (is(reads, "GAlignmentPairs")) {
        frags <- granges(reads)
    } else if (is(reads, "GAlignments")) {
        ## Extend smaller reads to fraglen, but don't shrink longer
        ## reads
        frags <- reads %>% granges %>% resize(width=pmax(fraglen, width(.)), fix="start")
    } else {
        warning(glue("Unknown reads type: {class(reads)}. Attempting to coerce to GRanges."))
        reads %<>% as("GRanges")
        ## If all ranges are the same width, assume they represent
        ## single-end reads and resize them to fraglen
        if (all(width(reads) == width(reads[1])) && width(reads[1]) < fraglen) {
            warning(glue("All reads from unknown class have the same length, {width(reads[1])}, and are therefore assumed to be single-end reads, which will be resized to {fraglen}."))
            frags <- resize(reads, width=fraglen, fix="start")
        } else {
            frags <- reads
        }
    }
    ## Finally, get the center of each fragment
    resize(frags, width=1, fix="center")
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

## TODO: Move to utilities.R

## Read a table from a R data file, csv, or xlsx file. Returns a data
## frame or throws an error.
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
        stop(glue("Could not read a data frame from {deparse{filename}} as R data, csv, or xlsx"))
    })
}

read.saf <- function(filename, ...) {
    saf <- read.table.general(filename, ...)
    assert_that("GeneID" %in% names(saf))
    gr <- as(saf, "GRanges")
    grl <- split(gr, gr$GeneID) %>% promote.common.mcols
    return(grl)
}

# Functions for reading and writing narrowPeak files
read.narrowPeak <- function(file, ...) {
    peaks.df <- read.table(file, sep="\t", row.names=NULL, ...)
    names(peaks.df) <- c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")
    peaks.df$name <- as.character(peaks.df$name)
    peaks.df
}

write.narrowPeak <- function(x, file, ...) {
    x <- as(x, "data.frame")
    if("seqnames" %in% names(x))
        names(x)[names(x) == "seqnames"] <- "chr"
    x <- x[c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "summit")]
    write.table(x, file, sep="\t", row.names=FALSE, col.names=FALSE, ...)
}

read.regions <- function(filename) {
    suppressWarnings({
        lazy.results <- list(
            rdata=future(read.RDS.or.RDA(filename), lazy=TRUE),
            narrowPeak=future(read.narrowPeak(filename), lazy=TRUE),
            bed=future(import(filename, format="bed"), lazy=TRUE),
            gff=future(import(filename, format="gff"), lazy=TRUE),
            saf=future(read.saf(filename), lazy=TRUE),
            table=future(read.table.general(filename), lazy=TRUE))
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
        stop(glue("Could not read genomic regions from {deparse(filename)} as R data, narrowPeak, bed, gff, SAF, or csv"))
    })
}

print.var.vector <- function(v) {
    for (i in names(v)) {
        cat(i, ": ", deparse(v[[i]]), "\n", sep="")
    }
    invisible(v)
}

## cmdopts <- list(
##     samplemeta_file="saved_data/samplemeta-ChIPSeq.RDS",
##     sample_id_column="SRA_run",
##     bam_file_pattern="aligned/chipseq_bowtie2_hg38.analysisSet/{SAMPLE}/Aligned.bam",
##     targets="saved_data/tss_shoal_hg38.analysisSet_ensembl.85.RDS",
##     upstream_neighborhood=5000,
##     downstream_neighborhood=5000,
##     window_width=500,
##     initial_window_offset=0,
##     output_file="saved_data/tss-neighborhood-counts_hg38.analysisSet_ensembl.85_5kbp-radius_500bp-windows_147bp-reads.RDS",
##     read_extension=147,
##     blacklist="saved_data/ChIPSeq-merged-blacklist.bed",
##     blacklist_action="mark",
##     threads=4,
##     window_spacing=500)

{
    cmdopts <- get.options(commandArgs(TRUE))
    ## TODO: Don't use setwd in this or any other script
    tryCatch(setwd(file.path(dirname(na.omit(get_Rscript_filename())), "..")),
             error=function(...) tsmsg("WARNING: Could not determine script path. Ensure that you are already in the correct directory."))

    tsmsg("Args:")
    print.var.vector(cmdopts)

    if (cmdopts$threads > 1) tryCatch({
        suppressPackageStartupMessages({
        })
        registerDoParallel(cores=cmdopts$threads)
        register(DoparParam())
    }, error=function(...){
        tsmsg("Could not initialize parallel backend. Falling back to single-core mode.")
        cmdopts$threads <- 1
    })

    if (cmdopts$thread <= 1) {
        registerDoSEQ()
        register(SerialParam())
    }

    tsmsg(glue("Using {cmdopts$threads} cores."))

    tsmsg(glue("Assuming a fragment size of {format.bp(cmdopts$read_extension)} for unpaired reads."))

    tsmsg("Loading sample data")

    sample.table <- readRDS(cmdopts$samplemeta_file) %>%
        ## Compute full path to BAM file
        mutate(bam_file=glue(cmdopts$bam_file_pattern, SAMPLE=.[[cmdopts$sample_id_column]])) %>%
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

    tsmsg("Loading target positions")
    targets <- read.regions(cmdopts$targets)
    assert_that(is(targets, "GRanges"))
    if (any(strand(targets) == "*")) {
        warning("Some targets have no strand information, and will be treated as being on the plus strand.")
    }
    ## Reduce targets to 5-prime end only
    targets %<>% resize(width=1, fix="start")

    blacklist_regions <- GRanges()
    if (!is.null(cmdopts$blacklist)) {
        tsmsg("Loading blacklist regions")
        blacklist_regions <- read.regions(cmdopts$blacklist)
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
    nhood_windows <- rep(targets, each=length(nhood_offsets))
    nhood_windows$offset <- rep(nhood_offsets, length.out=length(nhood_windows))
    nhood_windows %<>%
        shift(.$offset * ifelse(strand(.) == "-", -1, 1)) %>%
        resize(width=cmdopts$window_width, fix="center")
    ## Add offset to window names
    names(nhood_windows) %<>% str_c(format.bp(nhood_windows$offset, use_si=FALSE, always_signed=TRUE))

    if (length(blacklist_regions) > 0) {
        blacklisted <- overlapsAny(nhood_windows, blacklist_regions, ignore.strand=TRUE)
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

    tsmsg(glue("Counting reads in neighborhoods around {length(targets)} regions in {nrow(sample.table)} samples."))
    tsmsg(glue("Neighborhoods consist of {length(nhood_offsets)} windows of width {format.bp(cmdopts$window_width)} tiled from {format.bp(cmdopts$upstream_neighborhood)} upstream to {format.bp(cmdopts$downstream_neighborhood)} downstream."))

    bfl <- BamFileList(sample.table$bam_file)
    names(bfl) <- sample.table[[cmdopts$sample_id_column]]

    sexp <- summarizeOverlaps(
        features=nhood_windows, reads=bfl,
        ## Get mapping stats in colData
        count.mapped.reads=TRUE,
        ## Read pairs as fragments, treat all fragments as unstranded
        fragments=TRUE, singleEnd=FALSE, ignore.strand=TRUE,
        ## Reduce every fragment to only the single base pair at its
        ## midpoint
        preprocess.reads=. %>% readsToFragmentMidpoints(fraglen=cmdopts$read_extension),
        ## If a fragment (midpoint) overlaps multiple features, count it
        ## once for each one.
        inter.feature=FALSE,
        ## This doesn't matter since every fragment is reduced to a single
        ## base pair.
        mode="IntersectionStrict")

    ## Set blacklisted window counts to NA, if requested
    if (cmdopts$blacklist_action == "setNA") {
        bl <- rowRanges(sexp)$blacklist
        assay(sexp, "counts")[bl,] <- NA
    }

    ## Add sample metadata to colData in front of mapping stats
    colData(sexp) %<>% cbind(sample.table, .)
    ## Save command and options in the metadata
    metadata(sexp)$cmd.name <- na.omit(c(get_Rscript_filename(), "csaw-count-neighborhoods.R"))[1]
    metadata(sexp)$cmd.opts <- cmdopts

    tsmsg("Saving output file")
    saveRDS(sexp, cmdopts$output_file)
    tsmsg("Finished.")
}
