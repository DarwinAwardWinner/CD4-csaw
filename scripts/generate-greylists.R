#!/usr/bin/env Rscript

getScriptPath <- function() {
    argv <-commandArgs()
    dir <- na.omit(stringr::str_match(argv, "^--file=(.*)$")[,2])[1]
    if (!is.na(dir) && !is.null(dir))
        return(dir)
}
setwd(file.path(dirname(getScriptPath()), ".."))

library(foreach)
library(doParallel)
options(mc.cores=parallel::detectCores())
registerDoParallel(cores=parallel::detectCores())
library(doRNG)
library(magrittr)
library(dplyr)
library(openxlsx)
library(GreyListChIP)
library(BSgenome.Hsapiens.UCSC.hg19)

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

tsmsg("Loading sample data")

sample.table <- read.xlsx("data_files/ChIP-Seq/sample-tables.xlsx", "Samples") %>%
    set_colnames(make.unique(colnames(.))) %>%
    ## Select/compute/rename desired columns
    transmute(
        Sample=title,
        ChIP=`characteristics:.sampletype`,
        Celltype=`characteristics:.celltype`,
        Activated=`characteristics:.activated`,
        Day=`characteristics:.days.after.activation`,
        Donor=`characteristics:.donor.ID`,
        file.name=raw.file) %>%
    ## Make sure no factorial variables can be accidentally
    ## numericized by prefixing them with letters
    mutate(
        Donor=sprintf("Dn%s", Donor),
        Day=sprintf("D%i", Day)
    ) %>%
    ## Compute full path to bam file
    mutate(
        bampath=file.path("data_files/ChIP-Seq", file.name)
    ) %>%
    ## Factor variables with proper levels (not ASCIIbetical order)
    mutate(
        ChIP=factor(ChIP),
        Celltype=factor(Celltype, levels=c("Naive", "Memory")),
        Day=factor(Day, levels=sprintf("D%s", c(0,1,5,14))),
        TreatmentGroup=interaction(Celltype, Day, sep=""),
        Group=interaction(ChIP, TreatmentGroup, sep="."),
        Donor=factor(Donor)
    ) %>% set_rownames(.$Sample)

stopifnot(all(file.exists(sample.table$bampath)))

input.sample.table <- sample.table %>% filter(ChIP == "input")

std.chr <- extractSeqlevels("Homo sapiens", "UCSC") %>% setdiff("chrM")
kt <- BSgenome.Hsapiens.UCSC.hg19 %>%
    seqinfo %>% keepSeqlevels(std.chr) %>%
    { data.frame(Chrom=seqnames(.), Length=seqlengths(.))}
## Write the karyotype to a file because there's no other way to
## provide it to GreyListChIP.
ktfile <- "saved_data/GreyListChIP-karyotype-stdChrom.txt"
write.table(kt, ktfile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

set.seed(1986)
glists <- foreach(bamfile = input.sample.table$bampath) %dorng% {
    tsmsg("Generating Grey List for ", bamfile)
    gl <- new("GreyList", karyoFile=ktfile)
    gl %<>% countReads(bamfile)
    gl %<>% calcThreshold(reps=100, p=0.99)
    gl %<>% makeGreyList(maxGap=2048)
}
saveRDS(glists, "saved_data/GreyLists-maxgap2048bp.RDS")
greylist.regions <- glists %>% lapply(slot, "regions") %>% GRangesList

## TODO: Re-implement greylist functionality on top of csaw
