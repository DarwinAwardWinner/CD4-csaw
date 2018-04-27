#!/usr/bin/env Rscript

library(here)
library(stringr)
library(glue)
library(magrittr)
library(GenomicRanges)
library(Rsubread)
library(openxlsx)
library(annotate)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(SummarizedExperiment)
library(dplyr)
library(ChIPQC)
library(rctutils)

## Don't run in parallel because of memory issues
library(doParallel)
options(mc.cores=1)
registerDoSEQ()
library(BiocParallel)
register(SerialParam())

tsmsg("Loading sample data")

## TODO: Rewrite into the new pipeline

sample.table <- read.xlsx(here("data_files/ChIP-Seq/sample-tables.xlsx"), "Samples") %>%
    set_colnames(make.unique(colnames(.))) %>%
    ## Select/compute/rename desired columns
    transmute(
        Sample=title,
        ChIP=`characteristics:.sampletype`,
        Celltype=`characteristics:.celltype`,
        Activated=`characteristics:.activated`,
        Day=`characteristics:.days.after.activation`,
        Donor=`characteristics:.donor.ID`,
        bam.file.name=raw.file,
        peak.file.name=processed.data.file.3 %>% str_replace("bb$", "bed")
    ) %>%
    ## Make sure no factorial variables can be accidentally
    ## numericized by prefixing them with letters
    mutate(
        Donor=glue("Dn{Donor}"),
        Day=glue("D{Day}")
    ) %>%
    ## Compute full path to files
    mutate(
        bampath=here("data_files/ChIP-Seq", bam.file.name),
        peakpath=here("data_files/ChIP-Seq", peak.file.name)
    ) %>%
    ## Reorder levels on factors (not ASCIIbetical order)
    mutate(
        ChIP=factor(ChIP),
        Celltype=factor(Celltype, levels=c("Naive", "Memory")),
        Day=factor(Day, levels=glue("D{day}", day=c(0,1,5,14))),
        TreatmentGroup=interaction(Celltype, Day, sep=""),
        Group=interaction(ChIP, TreatmentGroup, sep="."),
        Donor=factor(Donor)
    )

## > colnames(sample.table)
##  [1] "Sample"         "ChIP"           "Celltype"       "Activated"
##  [5] "Day"            "Donor"          "bam.file.name"      "bampath"
##  [9] "TreatmentGroup" "Group"

inputs <- sample.table %>% filter(ChIP == "input") %>%
    transmute(Celltype, Activated, Day, Donor, input.bampath=bampath)
chipqc.exp.df <- sample.table %>% filter(ChIP != "input") %>%
    droplevels %>% merge(inputs) %>%
    transmute(
        SampleID=Sample,
        Tissue=Celltype,
        Factor=ChIP,
        Treatment=Day,
        Condition=TreatmentGroup,
        Replicate=Donor,
        bamReads=bampath,
        bamControl=input.bampath,
        Peaks=peakpath
    )

cqc <- ChIPQC(
    experiment = chipqc.exp.df,
    annotation = "hg19",
    chromosomes = extractSeqlevels("Homo sapiens", "UCSC") %>%
        setdiff("chrM"),
    profileWin = 600)

## Error: https://support.bioconductor.org/p/68215/#80421
## temp <- readGAlignments(chipqc.exp.df$bamReads[1])
## Sample_GIT <- GNCList(GRanges(seqnames = seqnames(temp),
##                               ranges = ranges(temp), strand = strand(temp),
##                               elementMetadata(temp)))

saveRDS(cqc, here("saved_data/ChIPQC.RDS"))
save.image(here("saved_data/ChIPQC.rda"))
