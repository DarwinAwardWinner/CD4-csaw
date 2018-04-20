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
    library(future)

    library(doParallel)
    library(BiocParallel)
})

tsmsg <- function(...) {
    message(date(), ": ", ...)
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

format.bp <- function(x) {
    x %>% round %>% f2si("bp") %>% str_replace_all(rex(one_or_more(space)), "")
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
        cmdopts[[i]] %<>% parse.bp
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

windowCountsParallel <- function(bam.files, ..., filter=10, BPPARAM=bpparam()) {
    reslist <- bplapply(X=bam.files, FUN=windowCounts, ..., filter=0, BPPARAM=BPPARAM)
    assert_that(all(sapply(reslist, is, "SummarizedExperiment")))
    res <- do.call(cbind, reslist)
    rm(reslist)
    keep <- rowSums(assay(res)) >= filter
    res[keep,]
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
        stop(glue("Could not read genomic regions from {deparse(filename)} as R data, bed, gff, SAF, or csv"))
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

    blacklist.regions <- GRanges()
    if (!is.null(cmdopts$blacklist)) {
        tsmsg("Loading blacklist regions")
        blacklist.regions <- read.regions(cmdopts$blacklist)
        assert_that(is(blacklist.regions, "GRanges"))
        ## Blacklist applies to both strands
        strand(blacklist.regions) <- "*"
    }

    ## Standard nuclear chromosomes only. (chrM is excluded because it is
    ## not located in the nucleus and is thus not subject to histone
    ## modification. The unplaced scaffolds are mostly not large enough to
    ## contain even a single typically-sized peak, so little is lost by
    ## excluding them for this analysis.)
    std.chr <- extractSeqlevels("Homo sapiens", "UCSC") %>% setdiff("chrM")
    rparam <- readParam(discard=blacklist)

    tsmsg(glue("Counting reads in {width} {type} in {scount} samples.",
        width=format.bp(cmdopts$window_width),
        type=ifelse(cmdopts$bin, "bins", "windows"),
        scount=nrow(sample.table)))
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
    ## Save command and options in the metadata
    metadata(sexp)$cmd.name <- na.omit(c(get_Rscript_filename(), "csaw-count-windows.R"))[1]
    metadata(sexp)$cmd.opts <- cmdopts

    tsmsg("Saving output file")
    saveRDS(wcounts, cmdopts$output_file)
    tsmsg("Finished.")
}
