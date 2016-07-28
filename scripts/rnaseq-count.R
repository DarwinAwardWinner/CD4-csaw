#!/usr/bin/env Rscript

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

## NULL = required arg; anything else is a default, and will be
## coerced to character
argspec <- list(
    SAMPLEMETA_FILE=NULL,
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
