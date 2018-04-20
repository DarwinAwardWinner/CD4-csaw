#!/usr/bin/env Rscript

library(magrittr)
library(dplyr)
library(assertthat)
library(ggplot2)
library(ggforce)

ccf <- readRDS("saved_data/chipseq-ccf.RDS")
ccf.nbl <- readRDS("saved_data/chipseq-ccf-noBL.RDS")

sample.table <- readRDS("saved_data/samplemeta-ChIPSeq.RDS") %>%
    ## Ensure that days_after_activation is a factor and can't be
    ## interpreted as a numeric
    mutate(days_after_activation=days_after_activation %>%
               factor %>% `levels<-`(str_c("Day", levels(.)))) %>%
    rename(Sample=SampleName,
           ChIP=chip_antibody,
           TimePoint=days_after_activation,
           CellType=cell_type,
           Donor=donor_id) %>%
    mutate(TreatmentGroup=interaction(CellType, TimePoint, sep="."))

ccftable <- lapply(names(ccf), function(i) {
    ccfValues <- ccf[[i]]
    ccf.nblValues <- ccf.nbl[[i]]
    assert_that(length(ccfValues) == length(ccf.nblValues))
    data_frame(
        Sample=i, Delay=seq(from=0, length.out=length(ccfValues)),
        CCF=ccfValues,
        RelCCF=CCF/max(CCF),
        CCF.noBL=ccf.nblValues,
        RelCCF.noBL=CCF.noBL/max(CCF.noBL))
}) %>% do.call(what=rbind) %>% inner_join(sample.table, ., by="Sample")

refline.table <- data.frame(
    Reference=c("Read Length (100bp)", "Nucleosome Footprint (147bp)"),
    Xintercept=c(100, 147))

{
    baseplot <- ggplot(ccftable) +
        facet_wrap(~ChIP, scales="free") +
        aes(x=Delay, y=CCF, group=Sample, color=TreatmentGroup, linetype=NA) +
        ylim(0,NA) +
        geom_vline(data=refline.table,
                   aes(xintercept=Xintercept, linetype=Reference, color=NA),
                   color="black", alpha=0.5) +
        ## geom_rug(data=ccfmaxtable, sides="b") +
        scale_color_hue(name="Group") +
        scale_linetype(name="Reference") +
        theme(legend.position="bottom")
    p <- list(
        Raw=baseplot +
            geom_line(size=0.25, linetype="solid") +
            ggtitle("Cross-Correlation Function, Raw"),
        loess_span0.05=baseplot +
            geom_smooth(fill=NA, method="loess", span=0.05, n=500, size=0.25, linetype="solid") +
            ggtitle("Cross-Correlation Function, Loess-Smoothed (span = 0.05)"),
        loess_span0.075=baseplot +
            geom_smooth(fill=NA, method="loess", span=0.075, n=500, size=0.25, linetype="solid") +
            ggtitle("Cross-Correlation Function, Loess-Smoothed (span = 0.075)"),
        loess_span0.1=baseplot +
            geom_smooth(fill=NA, method="loess", span=0.1, n=500, size=0.25, linetype="solid") +
            ggtitle("Cross-Correlation Function, Loess-Smoothed (span = 0.1)"))
    pdf("plots/csaw/CCF-plots.pdf", width=12, height=8)
    print(p)
    dev.off()
    pdf("plots/csaw/CCF-plots-relative.pdf", width=12, height=8)
    print(lapply(p, . %>% add(aes(y=RelCCF))))
    dev.off()
    pdf("plots/csaw/CCF-plots-noBL.pdf", width=12, height=8)
    print(lapply(p, . %>% add(aes(y=CCF.noBL))))
    dev.off()
    pdf("plots/csaw/CCF-plots-relative-noBL.pdf", width=12, height=8)
    print(lapply(p, . %>% add(aes(y=RelCCF.noBL))))
    dev.off()
}


ccfmaxtable <- ccftable %>%
    filter(Delay>=50) %>%
    group_by(Sample, CellType, TimePoint, Donor, ChIP, TreatmentGroup) %>%
    summarize(Delay=Delay[which.max(CCF)],
              CCF=max(CCF),
              Delay.noBL=.$Delay[which.max(CCF.noBL)],
              CCF.noBL=max(CCF.noBL)) %>%
    ungroup

ccfmaxtable.noBL <- ccftable %>%
    filter(Delay>=50) %>%
    group_by(Sample) %>%
    do(.[which.max(.$CCF.noBL),]) %>%
    ungroup

lims <- range(c(ccfmaxtable$Delay, ccfmaxtable$Delay.noBL))

p <- ggplot(ccfmaxtable) +
    facet_wrap(~ChIP) +
    coord_fixed(xlim=lims, ylim=lims) +
    ## Reference guide lines with circles at intersection
    geom_hline(data=refline.table,
               aes(yintercept=Xintercept, linetype=Reference),
               color="grey40") +
    geom_vline(data=refline.table,
               aes(xintercept=Xintercept, linetype=Reference),
               color="grey40") +
    geom_circle(data=refline.table,
                aes(x0=Xintercept, y0=Xintercept, linetype=Reference, r=5),
                color="grey40", show_guide=FALSE) +
    ## Plot the actual data
    geom_point(aes(x=Delay, y=Delay.noBL, color=TreatmentGroup)) +
    scale_shape_manual(values=c(`Read Length (100bp)`=1, `Nucleosome Footprint (147bp)`=1)) +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    theme(legend.position="bottom") +
    xlab("Delay of Maximum Cross-Correlation (With Blacklist)") +
    ylab("Delay of Maximum Cross-Correlation (No Blacklist)") +
    ggtitle("Delay of Maximum Cross-Correlation With and Without Blacklist")

pdf("plots/csaw/CCF-max-plot.pdf", width=10, height=10)
print(p)
dev.off()
