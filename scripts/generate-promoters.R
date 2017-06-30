#!/usr/bin/env Rscript

library(getopt)
library(optparse)
library(assertthat)
num.cores <- 1
## Don't default to more than 4 cores
try({library(parallel); num.cores <- detectCores(); }, silent=TRUE)

match.arg <- function (arg, choices, several.ok = FALSE, argname=substitute(arg))
{
    if (missing(choices)) {
        formal.args <- formals(sys.function(sys.parent()))
        choices <- eval(formal.args[[as.character(substitute(arg))]])
    }
    if (is.null(arg))
        return(choices[1L])
    else if (!is.character(arg))
        stop(sprintf("%s must be NULL or a character vector", deparse(argname)))
    if (!several.ok) {
        if (identical(arg, choices))
            return(arg[1L])
        if (length(arg) > 1L)
            stop(sprintf("%s must be of length 1", deparse(argname)))
    }
    else if (length(arg) == 0L)
        stop(sprintf("%s must be of length >= 1", deparse(argname)))
    i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
    if (all(i == 0L))
        stop(gettextf("%s should be one of %s", deparse(argname), paste(dQuote(choices),
            collapse = ", ")), domain = NA)
    i <- i[i > 0L]
    if (!several.ok && length(i) > 1)
        stop("there is more than one match in 'match.arg'")
    choices[i]
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

    ## Do argument parsing early so the script exits quickly if arguments are invalid
    optlist <- list(
        ## So far this script only supports TxDb objects because
        ## figuring out the first exon and TSS from other
        ## less-structured formats is a pain.
        make_option(c("-t", "--annotation-txdb"), metavar="TXDBNAME", type="character",
                    help="Name of TxDb package, or the name of a database file, to use for gene annotation"),
        make_option(c("-r", "--promoter-radius"), metavar="RADIUS", type="character",
                    help="Maximum distance from a gene's transcription start site that is considered part of the promoter."),
        make_option(c("-o", "--output-file"), metavar="FILENAME.RDS", type="character",
                    help="Output file name. The GRanges object containing the promoter regions will be saved here using saveRDS, so it should end in '.RDS'."),
        make_option(c("-j", "--threads"), metavar="N", type="integer", default=num.cores,
                    help="Number of threads to use while counting reads"))
    progname <- na.omit(c(get_Rscript_filename(), "rnaseq-count.R"))[1]
    parser <- OptionParser(
        usage="Usage: %prog [options] -t TXDB -r RADIUS -o OUTPUT.RDS",
        description="Generate promoter regions for a gene annotation.

From each transcription start site, a region is extended to the specified radius in both directions to define the promoter region for that TSS. Then, any overlapping promoters that share a gene ID are merged. Note that the number of non-overlapping promoters necessarily depends on the chosen promoter radius, and different radii will give different numbers of promoters. The resulting GRanges object will have several annotation columns added: 'PromoterID', which is an arbitrary unique key for each promoter; 'GeneID', which may be NA, may be equal to the TxID for transcripts with no Gene ID, and may be non-unique if a single gene has multiple non-overlapping promoters; TxID, a CharacterList of all transcript IDs whose TSS are included in the promoter.",
option_list = optlist,
add_help_option = TRUE,
prog = progname,
epilogue = "Note that all base pair sizes (for the promoter radius) may have an suffix of 'bp', 'kbp', 'mbp', or 'tbp' (and the 'p' is optional). For example, 10kb = 10000")

    cmdopts <- parse_args(parser, opts)
    ## Ensure that all required arguments were provided
    required.opts <- c("annotation-txdb", "output-file", "promoter-radius")
    missing.opts <- setdiff(required.opts, names(cmdopts))
    if (length(missing.opts) > 0) {
        stop(str_c("Missing required arguments: ", deparse(missing.opts)))
    }

    ## Convert bp args to numbers
    for (i in c("promoter-radius")) {
        cmdopts[[i]] %<>% parse.bp
    }

    ## Replace dashes with underscores so that all options can easily
    ## be accessed by "$"
    cmdopts %>% setNames(str_replace_all(names(.), "-", "_"))
}

## Do this early to handle "--help" before wasting time loading
## pacakges & stuff
invisible(get.options(commandArgs(TRUE)))

library(assertthat)
library(dplyr)
library(magrittr)
library(stringr)
library(GenomicRanges)
library(BiocParallel)
library(doParallel)
options(mc.preschedule=FALSE)
register(DoparParam())

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

is.empty <- function(x) {
    x %>% unlist %>% na.omit %>% length %>% equals(0)
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

print.var.vector <- function(v) {
    for (i in names(v)) {
        cat(i, ": ", deparse(v[[i]]), "\n", sep="")
    }
    invisible(v)
}

{
    cmdopts <- get.options(commandArgs(TRUE))
    cmdopts$help <- NULL

    tsmsg("Args:")
    print.var.vector(cmdopts)

    cmdopts$threads %<>% round %>% max(1)
    tsmsg("Running with ", cmdopts$threads, " threads")
    registerDoParallel(cores=cmdopts$threads)

    ## Delete the output file if it exists
    suppressWarnings(file.remove(cmdopts$output_file))
    assert_that(!file.exists(cmdopts$output_file))

    tsmsg("Reading annotation data")
    txdb <- get.txdb(cmdopts$annotation_txdb)
    tsmsg(sprintf("Getting %s-radius promoters", format.bp(cmdopts$promoter_radius)))
    all.promoters <- suppressWarnings(promoters(txdb, upstream=cmdopts$promoter_radius, downstream = cmdopts$promoter_radius)) %>%
        trim
    tsmsg("Annotating promoters")
    mcols(all.promoters) %<>%
        ## Not using dplyr because it's a BioC DataFrame
        transform(TxID=tx_name,
                  tx_name=NULL,
                  tx_id=NULL)
    all.promoters$GeneID <- suppressMessages(mapIds(txdb, all.promoters$TxID, keytype="TXNAME", column="GENEID", multiVals="first")) %>%
        ifelse(is.na(.), all.promoters$TxID, .)
    tsmsg("Splitting promoters by gene ID")
    gene.promoters <- split(all.promoters, all.promoters$GeneID)
    assert_that(all(lengths(gene.promoters) >= 1))
    tsmsg("Merging overlapping promoters from the same gene")
    merged.promoters <- bplapply(gene.promoters, function(gp) {
        gp.reduced <- reduce(gp)
        gp.reduced$GeneID <- gp$GeneID[1]
        names(gp.reduced) <- gp.reduced$PromoterID <- sprintf("%s-P%s", gp.reduced$GeneID, seq_along(gp.reduced))
        pgroup <- gp.reduced$PromoterID[nearest(gp, gp.reduced)]
        gp.reduced$TxID <- CharacterList(split(gp$TxID, pgroup))[gp.reduced$PromoterID]
        gp.reduced
    }) %>% unname %>% GRangesList %>% unlist

    tsmsg("Saving output file")
    save.RDS.or.RDA(merged.promoters, cmdopts$output_file)
    invisible(NULL)
}
