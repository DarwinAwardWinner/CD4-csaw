#!/usr/bin/env Rscript

getScriptPath <- function() {
    argv <-commandArgs()
    dir <- na.omit(stringr::str_match(argv, "^--file=(.*)$")[,2])[1]
    if (!is.na(dir) && !is.null(dir))
        return(dir)
}
setwd(file.path(dirname(getScriptPath()), ".."))

library(stringr)
library(magrittr)
library(GenomicRanges)
library(Rsubread)
library(openxlsx)
library(annotate)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(doParallel)
options(mc.cores=parallel::detectCores())
registerDoParallel(cores=parallel::detectCores())
library(GenomicRanges)
library(SummarizedExperiment)
library(plyr)
library(dplyr)

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

samplemeta.file <- "saved_data/samplemeta-RNASeq.RDS"
sumexp.outfile <- "saved_data/SummarizedExperiment-RNASeq.RDS"

tsmsg("Loading sample metadata")
samplemeta <- readRDS(samplemeta.file) %>%
    mutate(bam_file=file.path("bam_files", str_c(SRA_run, ".bam")))

tsmsg("Reading annotation data")
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

gene.exons <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, "gene")
## Add gene annotations as metadata on the exons-by-gene list
metacols <- c("SYMBOL", "GENENAME", "UNIGENE", "ENSEMBL", "REFSEQ", "UCSCKG", "UNIPROT")
gene.annot <-
  DataFrame(row.names=names(gene.exons),
            ENTREZ=names(gene.exons),
            llply(setNames(nm=metacols),
                  function(i) {
                    CharacterList(lookUp(names(gene.exons),
                                         data="org.Hs.eg",
                                         what=i))
                  }, .parallel=FALSE))
mcols(gene.exons) <- gene.annot
saf <- grl.to.saf(gene.exons)

tsmsg("Computing sense counts")
sense.fc <- featureCounts(
    samplemeta$bam_file, annot.ext=saf,
    strandSpecific=1,
    nthreads=getOption("mc.cores", 2))
tsmsg("Computing antisense counts")
antisense.fc <- featureCounts(
    samplemeta$bam_file, annot.ext=saf,
    strandSpecific=2,
    nthreads=getOption("mc.cores", 2))
tsmsg("Computing unstranded counts")
unstranded.fc <- featureCounts(
    samplemeta$bam_file, annot.ext=saf,
    strandSpecific=0,
    nthreads=getOption("mc.cores", 2))

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

saveRDS(sexp, sumexp.outfile)
