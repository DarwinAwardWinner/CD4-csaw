#!/usr/bin/env Rscript

tsmsg <- function(...) {
    message(date(), ": ", ...)
}

getScriptPath <- function() {
    argv <-commandArgs()
    dir <- na.omit(stringr::str_match(argv, "^--file=(.*)$")[,2])[1]
    if (!is.na(dir) && !is.null(dir))
        return(dir)
}
tryCatch(setwd(file.path(dirname(getScriptPath()), "..")),
         error=function(...) tsmsg("WARNING: Could not determine script path. Ensure that you are already in the correct directory."))

library(stringr)
library(glue)
library(magrittr)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(reshape2)
library(purrr)
library(csaw)
library(Matrix)
library(assertthat)

library(doParallel)
options(mc.preschedule=FALSE)
ncores <- getOption("mc.cores", default=1)
registerDoParallel(cores=ncores)
library(BiocParallel)
register(DoparParam())

tsmsg("Loading sample data")

sample.table <- readRDS("saved_data/samplemeta-ChIPSeq.RDS") %>%
    ## Compute full path to BAM file
    mutate(bam_file=glue("aligned/chipseq_bowtie2_hg38.analysisSet/{SRA_run}/Aligned.bam")) %>%
    ## Ensure that days_after_activation is a factor and can't be
    ## interpreted as a numeric
    mutate(days_after_activation=days_after_activation %>%
               factor %>% `levels<-`(str_c("Day", levels(.)))) %>%
    rename(time_point=days_after_activation) %>%
    mutate(treatment_group=interaction(cell_type, time_point, drop=TRUE, lex.order = TRUE))

assert_that(all(file.exists(sample.table$bam_file)))

tsmsg("Loading blacklist regions")
blacklist <- import("saved_data/ChIPSeq-merged-blacklist.bed", format="bed")

## Standard nuclear chromosomes only. (chrM is excluded because it is
## not located in the nucleus and is thus not subject to histone
## modification. The unplaced scaffolds are mostly not large enough to
## contain even a single typically-sized peak, so little is lost by
## excluding them for this analysis.)
std.chr <- extractSeqlevels("Homo sapiens", "UCSC") %>% setdiff("chrM")
param <- readParam(restrict=std.chr, discard=blacklist)
param.dedup.on <- reform(param, dedup=TRUE)

## Determining window size using findMaxima & profileSites, see csaw
## UG 2.5

tsmsg("Finding maxima")
sample.maxima <- bplapply(sample.table$bam_file, function(bam) {
    windowed <- windowCounts(bam, spacing=50, width=50, ext=147, param=param.dedup.on, filter=20)
    rwsms <- rowSums(assay(windowed))
    maxed <- findMaxima(rowRanges(windowed), range=5000, metric=rwsms)
    maxwindowed <- windowed[maxed,]
    maxranges <- rowRanges(maxwindowed)
    maxranges$Count <- rwsms[maxed]
    maxranges
})
names(sample.maxima) <- sample.table$SampleName
saveRDS(sample.maxima, "saved_data/chipseq-sample-maxima.RDS")

weights <- lapply(sample.maxima, . %$% {1/Count})

tsmsg("Profiling maxima")
sample.mean.profiles <- bpmapply(
    profileSites,
    bam.files=sample.table$bam_file,
    regions=sample.maxima,
    weight=weights,
    MoreArgs = list(
        average = TRUE,
        ext=147,
        range=10000,
        param=param.dedup.on))
colnames(sample.mean.profiles) <- sample.table$SampleName
saveRDS(sample.mean.profiles, "saved_data/chipseq-siteprof.RDS")

profile.table <- sample.mean.profiles %>%
    melt(varnames=c("Distance", "SampleName"), value.name="RelativeCoverage") %>%
    mutate(Distance=Distance-25) %>%
    inner_join(sample.table, by="SampleName") %>%
    group_by(SampleName)

{
    baseplot <- ggplot(profile.table) +
        facet_wrap(~chip_antibody, scales="free") +
        aes(x=Distance, y=RelativeCoverage, group=SampleName, color=treatment_group) +
        ylim(0,NA)
    p <- list(
        Raw10kb=baseplot +
            geom_line(size=0.25) +
            coord_cartesian(xlim=c(-10000, 10000)) +
            ggtitle("Relative Coverage Around Maxima, 10kb Radius"),
        Raw1kb=baseplot +
            geom_line(size=0.25) +
            coord_cartesian(xlim=c(-1000, 1000)) +
            ggtitle("Relative Coverage Around Maxima, 1kb Radius"),
        Raw300bp=baseplot +
            geom_line(size=0.25) +
            coord_cartesian(xlim=c(-300, 300)) +
            ggtitle("Relative Coverage Around Maxima, 300bp Radius"))
    pdf("plots/csaw/site-profile-plots.pdf", width=16, height=16)
    print(p)
    dev.off()
}
