#!/usr/bin/env Rscript

library(getopt)
library(optparse)

## Do argument parsing early so the script exits quickly if arguments are invalid
optlist <- list(
    make_option(c("-s", "--samplemeta-file"), metavar="FILENAME.RDS", type="character",
                help="(REQUIRED) RDS file containing a data frame of sample metadata."),
    make_option(c("-c", "--sample-id-column"), type="character", default="Sample",
                help="Sample metadata column name that holds the sample IDs. These will be substituted into '--bam-file-pattern' to determine the BAM file names."),
    make_option(c("-p", "--bam-file-pattern"), metavar="PATTERN", type="character",
                help="(REQUIRED) Format string to convert sample IDs into BAM file paths. This should contain a '%s' wherever the sample ID should be substituted ('%s' can occur multiple times),. Example: 'bam_files/Sample_%s/Aligned.bam"),
    make_option(c("-o", "--output-file"), metavar="FILENAME.RDS", type="character",
                help="(REQUIRED) Output file name. The SummarizedExperiment object containing the counts will be saved here using saveRDS, so it should end in '.RDS'."),
    make_option(c("-b", "--expected-bam-files"), metavar="BAMFILE1,BAMFILE2,...", type="character",
                help="Comma-separated list of bam file names expected to be used as input. This argument is optional, but if it is provided, it will be checked against the list of files determined from '--samplemeta-file' and '--bam-file-pattern', and an error will be raised if they don't match exactly."),
    make_option(c("-j", "--threads"), metavar="N", type="integer", default=1,
                help="Number of threads to use while counting reads"),
    ## TODO: Allow different annotations, via txdb, or gff file
    make_option(c("-t", "--annotation-txdb"), metavar="PACKAGENAME", type="character",
                help="Name of TxDb package to use for gene annotation"),
    make_option(c("-g", "--annotation-gff"), metavar="FILENAME", type="character",
                help="File Name of GFF3 file to use for gene annotation."),
    make_option(c("-f", "--gff-featuretype"), metavar="FEATURETYPE", type="character", default="exon",
                help="GFF feature type to use"),
    make_option(c("-i", "--gff-geneid-attr"), metavar="ATTRNAME", type="character", default="gene_id",
                help="GFF feature attribute to use as a feature's Gene ID."),
    make_option(c("-r", "--annotation-rds"), metavar="FILENAME", type="character",
                help="File Name of RDS or RData file to use for gene annotation. It should contain a single GRanges or GRangesList object, with each element representing a gene/feature to be counted. Metadata columns will be carried through to the output SummarizedExperiment."),
    make_option(c("-a", "--additional-gene-info"), metavar="FILENAME", type="character",
                help="RDS/RData/xlsx/csv file containing a table of gene metadata. Row names should be gene IDs."))
progname <- na.omit(c(get_Rscript_filename(), "rnaseq-count.R"))[1]
parser <- OptionParser(
    usage="Usage: %prog [options] -s SAMPLEMETA.RDS -p PATTERN -o SUMEXP.RDS [ -t TXDB.PACKAGE.NAME | -g ANNOTATION.GFF3 | -r ANNOTATION.RDS ]",
    description="Count reads in genes using Rsubread::featureCounts.

Counts are performed for stranded, non-stranded, and reverse-stranded modes, and are stored along with the sample and gene metadata in a SummarizedExperiment object. Note that the '-s', '-p', and '-o' arguments are all required, since they specify the input and output files. Also, exactly one of '-t', '-g', or '-r' is required to specify the annotation.",
    option_list = optlist,
    add_help_option = TRUE,
    prog = progname,
    epilogue = "")
cmdopts <- parse_args(parser, commandArgs(TRUE))
required.args <- c("samplemeta-file", "bam-file-pattern", "output-file")

## myargs <- c("-s", "samplemeta.RDS", "-p", "bam_files/Sample_%s/Aligned.bam", "-o", "sexp.rds", "-j", "2")
## parse_args(parser, "-h")
## cmdargs <- parse_args(parser, myargs)

missing.required.args <- setdiff(required.args, names(cmdargs))
if (length(missing.required.args) > 0) {
    stop(str_c("Missing required arguments: ", repr(missing.required.args)))
}

stop("Unimplemented")

library(stringr)
library(magrittr)
library(dplyr)
library(assertthat)

library(GenomicRanges)
library(Rsubread)
library(annotate)
library(SummarizedExperiment)

library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

## This merges exons into genes (GRanges to GRangesList)
gr.to.grl <- function(gr, featureType="exon", attrType="gene_id") {
    gr <- gr[gr$type %in% featureType]
    split(gr, as.character(mcols(gr)[[attrType]]))
}

## This converts a GRangesList into the SAF ("Simplified annotation
## format")
grl.to.saf <- function(grl) {
    gr <- unlist(grl)
    data.frame(Chr=as.vector(seqnames(gr)),
               Start=start(gr),
               End=end(gr),
               Strand=as.vector(strand(gr)),
               GeneID=rep(names(grl), lengths(grl)))
}

parse.vars.from.args <- function(...) {
    args <- unlist(list(...))
    if (length(args) == 0)
        return(character(0))
    m <- str_match(args, "^(.+?)=(.*)$")
    invalid <- is.na(m[,1])
    if (any(invalid)) {
        invalid.args <- args[invalid]
        stop("The following arguments could not be parsed as variable specifications: ",
             deparse(invalid.args))
    }
    ## Now guaranteed none are missing
    setNames(m[,3], nm=m[,2])
}

print.var.vector <- function(v) {
    for (i in names(v)) {
        cat(i, ": ", deparse(v[[i]]), "\n", sep="")
    }
    v
}

## Like sprintf, but inserts the same value into every placeholder
sprintf.single.value <- function(fmt, value) {
    ## Max function arguments is 100
    arglist = c(list(fmt=fmt), rep(list(value), 99))
    do.call(sprintf, arglist)
}

## NULL = required arg; anything else is a default, and will be
## coerced to character
argspec <- list(
    SAMPLEMETA_FILE=NULL,
    BAM_FILE_PATTERN=NULL,
    SUMEXP_OUTPUT_FILE=NULL,
    BAM_FILES="",
    THREADS=getOption("mc.cores", 2))

parse.args <- function(argspec, args=commandArgs(TRUE)) {
    assert_that(max(lengths(argspec)) <= 1)

    arg.vars <- parse.vars.from.args(args)
    wanted <- names(arg.vars) %in% names(argspec)
    extra.arg.vars <- arg.vars[!wanted]
    if (length(extra.arg.vars)) {
        tsmsg("Warning: got extra args: ")
        print.var.vector(extra.arg.vars)
    }
    arg.vars <- arg.vars[wanted]
    required.argnames <- names(argspec)[lengths(argspec) == 0]
    missing.argnames <- setdiff(required.argnames, names(arg.vars))
    if (length(missing.argnames) > 0) {
        stop(str_c("Missing arguments: ", str_c(missing.argnames, sep=", ")))
    }
    assert_that(all(required.argnames %in% names(arg.vars)))

    optional.arg.defaults <- unlist(argspec[lengths(argspec) == 1])
    defaults.to.use <- optional.arg.defaults %>%
        .[! names(.) %in% names(arg.vars)]

    arg.vars %<>% c(defaults.to.use)

    return(lapply(arg.vars, as.character))
}

{
    args <- parse.args(argspec)

    tsmsg("Args:")
    print.var.vector(args)

    num.threads <- as.numeric(args$THREADS)
    if (str_length(args$BAM_FILES) == 0) {
        bam.files <- NULL
    } else {
        bam.files <- str_split(args$BAM_FILES, ",")[[1]]
    }
    saveRDS(NULL, args$SUMEXP_OUTPUT_FILE)
    file.remove(args$SUMEXP_OUTPUT_FILE)

    tsmsg("Loading sample metadata")
    samplemeta <- readRDS(args$SAMPLEMETA_FILE) %>%
        mutate(bam_file=file.path("bam_files", str_c(SRA_run, ".bam")))

    if (!is.null(bam.files)) {
        assert_that(setequal(samplemeta$bam_file, bam.files))
    }

    tsmsg("Reading annotation data")

    gene.exons <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
    ## Add gene annotations as metadata on the exons-by-gene list
    metacols <- c("SYMBOL", "GENENAME", "UNIGENE", "ENSEMBL", "REFSEQ", "UCSCKG", "UNIPROT")
    gene.annot <-
        DataFrame(row.names=names(gene.exons),
                  ENTREZ=names(gene.exons),
                  lapply(setNames(nm=metacols),
                         function(i) {
                             CharacterList(lookUp(names(gene.exons),
                                                  data="org.Hs.eg",
                                                  what=i))
                         }))
    mcols(gene.exons) <- gene.annot
    saf <- grl.to.saf(gene.exons)

    tsmsg("Computing sense counts")
    sense.fc <- featureCounts(
        samplemeta$bam_file, annot.ext=saf,
        strandSpecific=1,
        nthreads=num.threads)
    tsmsg("Computing antisense counts")
    antisense.fc <- featureCounts(
        samplemeta$bam_file, annot.ext=saf,
        strandSpecific=2,
        nthreads=num.threads)
    tsmsg("Computing unstranded counts")
    unstranded.fc <- featureCounts(
        samplemeta$bam_file, annot.ext=saf,
        strandSpecific=0,
        nthreads=num.threads)

    tsmsg("Saving SummarizedExperiment")
    sexp <- SummarizedExperiment(
        assays=List(
            counts=unstranded.fc$counts,
            sense.counts=sense.fc$counts,
            antisense.counts=antisense.fc$counts),
        colData=as(samplemeta, "DataFrame"),
        rowRanges=gene.exons,
        metadata=List(
            stat=List(
                counts=unstranded.fc$stat,
                sense.counts=sense.fc$stat,
                antisense.counts=antisense.fc$stat) %>%
                endoapply(. %>% { set_colnames(t(.[-1]), .[[1]])} %>% DataFrame)))
    colnames(sexp) <- colData(sexp)$title
    rownames(sexp) <- mcols(sexp)$ENTREZ

    saveRDS(sexp, args$SUMEXP_OUTPUT_FILE)
    invisible(NULL)
}
