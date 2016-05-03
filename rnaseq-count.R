#!/usr/bin/env Rscript

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

setwd("~/Projects/CD4-csaw")

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

tsmsg("Loading sample data")

sample.table <- read.xlsx("data_files/RNA-Seq/sample-tables.xlsx", "Samples") %>%
    rename(
        Celltype=`characteristics:.celltype`,
        Activated=`characteristics:.activated`,
        Day=`characteristics:.days.after.activation`,
        Donor=`characteristics:.donor.ID`,
        Batch=`characteristics:.technical.batch`,
        file.name=raw.file) %>%
    set_colnames(make.names(colnames(.))) %>%
    ## Drop some columns we don't care about
    dplyr::select(
        -Sample.name,
        -source.name,
        -organism,
        -molecule,
        -description,
        -processed.data.file
    ) %>%
    ## Make sure no factorial variables can be accidentally
    ## numericized
    mutate(
        Batch=sprintf("B%i", Batch),
        Donor=sprintf("Dn%s", Donor),
        Day=sprintf("D%i", Day)
    ) %>%
    ## Compute full path to bam file
    mutate(
        bampath=file.path("data_files/RNA-Seq", file.name)
    )

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
               GeneID=rep(names(grl), elementLengths(grl)))
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

sense.fc <- featureCounts(
    sample.table$bampath, annot.ext=saf,
    strandSpecific=1,
    nthreads=getOption("mc.cores", 2))
antisense.fc <- featureCounts(
    sample.table$bampath, annot.ext=saf,
    strandSpecific=2,
    nthreads=getOption("mc.cores", 2))
unstranded.fc <- featureCounts(
    sample.table$bampath, annot.ext=saf,
    strandSpecific=0,
    nthreads=getOption("mc.cores", 2))

sexp <- SummarizedExperiment(
    assays=List(
        counts=unstranded.fc$counts,
        sense.counts=sense.fc$counts,
        antisense.counts=antisense.fc$counts),
    colData=as(sample.table, "DataFrame"),
    rowRanges=gene.exons,
    metadata=List(
        stat=List(
            counts=unstranded.fc$stat,
            sense.counts=sense.fc$stat,
            antisense.counts=antisense.fc$stat) %>%
            endoapply(. %>% { set_colnames(t(.[-1]), .[[1]])} %>% DataFrame)))
colnames(sexp) <- colData(sexp)$title
rownames(sexp) <- mcols(sexp)$ENTREZ

saveRDS(sexp, "saved_data/RNASeq-SummarizedExperiment.RDS")

save.image("saved_data/rnaseq-count.rda")
