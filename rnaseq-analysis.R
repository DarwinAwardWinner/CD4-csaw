#!/usr/bin/env Rscript

library(stringr)
library(magrittr)
library(openxlsx)
library(doParallel)
options(mc.cores=parallel::detectCores())
registerDoParallel(cores=parallel::detectCores())
library(SummarizedExperiment)
library(dplyr)
library(edgeR)
library(limma)
library(GGally)

suppressPlot <- function(arg) {
    dev.new()
    result <- arg
    dev.off()
    result
}

add.numbered.colnames <- function(x, prefix="C") {
    x %>% set_colnames(sprintf("%s%i", prefix, seq(from=1, length.out=ncol(x))))
}

## Version of cpm that uses an offset matrix instead of lib sizes
cpmWithOffset <- function(dge, offset=expandAsMatrix(getOffset(dge), dim(dge)),
                          log = FALSE, prior.count = 0.25, ...) {
    x <- dge$counts
    effective.lib.size <- exp(offset)
    if (log) {
        prior.count.scaled <- effective.lib.size/mean(effective.lib.size) * prior.count
        effective.lib.size <- effective.lib.size + 2 * prior.count.scaled
    }
    effective.lib.size <- 1e-06 * effective.lib.size
    if (log)
        log2((x + prior.count.scaled) / effective.lib.size)
    else x / effective.lib.size
}

## Like limma's removeBatchEffect, but for DGEList. Modifies the
## offsets instead of the data.
offsetBatchEffect <- function (dge, batch = NULL, batch2 = NULL, covariates = NULL,
                               design = matrix(1, ncol(dge), 1), ...) {
    if (is.null(batch) && is.null(batch2) && is.null(covariates))
        return(dge)
    if (!is.null(batch)) {
        batch <- as.factor(batch)
        contrasts(batch) <- contr.sum(levels(batch))
        batch <- model.matrix(~batch)[, -1, drop = FALSE]
    }
    if (!is.null(batch2)) {
        batch2 <- as.factor(batch2)
        contrasts(batch2) <- contr.sum(levels(batch2))
        batch2 <- model.matrix(~batch2)[, -1, drop = FALSE]
    }
    if (!is.null(covariates))
        covariates <- as.matrix(covariates)
    X.batch <- cbind(batch, batch2, covariates)
    fit <- glmFit(dge, cbind(design, X.batch), ..., priot.count = 0)
    beta <- fit$coefficients[, -(1:ncol(design)), drop = FALSE]
    beta[is.na(beta)] <- 0
    dge$offset <- fit$offset + beta %*% t(X.batch)
    dge
}

## Version of voom that uses an offset matrix instead of lib sizes
voomWithOffset <- function (dge, design = NULL, offset=expandAsMatrix(getOffset(dge), dim(dge)),
                            normalize.method = "none", plot = FALSE, span = 0.5, ...)
{
    out <- list()
    out$genes <- dge$genes
    out$targets <- dge$samples
    if (is.null(design) && diff(range(as.numeric(counts$sample$group))) >
        0)
        design <- model.matrix(~group, data = counts$samples)
    counts <- dge$counts
    if (is.null(design)) {
        design <- matrix(1, ncol(counts), 1)
        rownames(design) <- colnames(counts)
        colnames(design) <- "GrandMean"
    }

    effective.lib.size <- exp(offset)

    y <- log2((counts + 0.5)/(effective.lib.size + 1) * 1e+06)
    y <- normalizeBetweenArrays(y, method = normalize.method)
    fit <- lmFit(y, design, ...)
    if (is.null(fit$Amean))
        fit$Amean <- rowMeans(y, na.rm = TRUE)
    sx <- fit$Amean + mean(log2(effective.lib.size + 1)) - log2(1e+06)
    sy <- sqrt(fit$sigma)
    allzero <- rowSums(counts) == 0
    if (any(allzero)) {
        sx <- sx[!allzero]
        sy <- sy[!allzero]
    }
    l <- lowess(sx, sy, f = span)
    if (plot) {
        plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )",
            pch = 16, cex = 0.25)
        title("voom: Mean-variance trend")
        lines(l, col = "red")
    }
    f <- approxfun(l, rule = 2)
    if (fit$rank < ncol(design)) {
        j <- fit$pivot[1:fit$rank]
        fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[,
            j, drop = FALSE])
    }
    else {
        fitted.values <- fit$coef %*% t(fit$design)
    }
    fitted.cpm <- 2^fitted.values
    ## fitted.count <- 1e-06 * t(t(fitted.cpm) * (lib.size + 1))
    fitted.count <- 1e-06 * fitted.cpm * (effective.lib.size + 1)
    fitted.logcount <- log2(fitted.count)
    w <- 1/f(fitted.logcount)^4
    dim(w) <- dim(fitted.logcount)
    out$E <- y
    out$weights <- w
    out$design <- design
    out$effective.lib.size <- effective.lib.size
    if (is.null(out$targets))
        out$targets <- data.frame(lib.size = exp(colMeans(offset)))
    else out$targets$lib.size <- exp(colMeans(offset))
    new("EList", out)
}


## Version of voom that uses an offset matrix instead of lib sizes
voomWithQualityWeightsAndOffset <-
    function (dge, design = NULL,
              offset=expandAsMatrix(getOffset(dge), dim(dge)),
              normalize.method = "none",
              plot = FALSE, span = 0.5, var.design = NULL, method = "genebygene",
              maxiter = 50, tol = 1e-10, trace = FALSE, replace.weights = TRUE,
              col = NULL, ...)
{
    counts <- dge$counts
    if (plot) {
        oldpar <- par(mfrow = c(1, 2))
        on.exit(par(oldpar))
    }
    v <- voomWithOffset(dge, design = design, offset = offset, normalize.method = normalize.method,
        plot = FALSE, span = span, ...)
    aw <- arrayWeights(v, design = design, method = method, maxiter = maxiter,
        tol = tol, var.design = var.design)
    v <- voomWithOffset(dge, design = design, weights = aw, offset = offset,
        normalize.method = normalize.method, plot = plot, span = span,
        ...)
    aw <- arrayWeights(v, design = design, method = method, maxiter = maxiter,
        tol = tol, trace = trace, var.design = var.design)
    wts <- asMatrixWeights(aw, dim(v)) * v$weights
    attr(wts, "arrayweights") <- NULL
    if (plot) {
        barplot(aw, names = 1:length(aw), main = "Sample-specific weights",
            ylab = "Weight", xlab = "Sample", col = col)
        abline(h = 1, col = 2, lty = 2)
    }
    if (replace.weights) {
        v$weights <- wts
        v$sample.weights <- aw
        return(v)
    }
    else {
        return(wts)
    }
}

estimateDispByGroup <- function(dge, group=as.factor(dge$samples$group), batch, ...) {
    stopifnot(nlevels(group) > 1)
    stopifnot(length(group) == ncol(dge))
    if (!is.list(batch)) {
        batch <- list(batch=batch)
    }
    batch <- as.data.frame(batch)
    stopifnot(nrow(batch) == ncol(dge))
    colnames(batch) %>% make.names(unique=TRUE)
    igroup <- seq_len(ncol(dge)) %>% split(group)
    lapply(igroup, function(i) {
        group.dge <- dge[,i]
        group.batch <- droplevels(batch[i,, drop=FALSE])
        group.batch <- group.batch[sapply(group.batch, . %>% unique %>% length %>% is_greater_than(1))]
        group.vars <- names(group.batch)
        if (length(group.vars) == 0)
            group.vars <- "1"
        group.batch.formula <- as.formula(str_c("~", str_c(group.vars, collapse="+")))
        des <- model.matrix(group.batch.formula, group.batch)
        estimateDisp(group.dge, des, ...)
    })
}

dir.create("results/RNA-seq/", FALSE, TRUE)

sexp <- readRDS("saved_data/RNASeq-SummarizedExperiment.RDS")
colnames(sexp) <- colData(sexp)$title
rownames(sexp) <- mcols(sexp)$ENTREZ

sample.table <- colData(sexp) %>% as.data.frame %>%
    mutate(Celltype=factor(Celltype, levels=c("Naive", "Memory")),
           Day=factor(Day, levels=sprintf("D%s", c(0,1,5,14))),
           Group=interaction(Celltype, Day, sep=""),
           Donor=factor(Donor),
           Batch=factor(Batch),
           CountType=factor(ifelse(Batch == "B1", "sense.counts", "antisense.counts")))
contrasts(sample.table$Donor) <- contr.sum(nlevels(sample.table$Donor))

design <- model.matrix(~0 + Group + Donor, sample.table)
colnames(design)[1:nlevels(sample.table$Group)] = levels(sample.table$Group)

assay(sexp, "correct.counts") <- lapply(seq_len(nrow(sample.table)), function(i) {
    assay(sexp, sample.table[i,]$CountType %>% as.character)[,i]
}) %>% do.call(what=cbind)

## Sanity check on selection of correct strand sense
total.counts <- sexp %>% assays %>% sapply(colSums) %>% data.frame
total.counts %$% stopifnot(all(correct.counts == pmax(sense.counts, antisense.counts)))

## Extract gene metadata and colapse lists
all.gene.meta <- mcols(sexp) %>% as.data.frame
all.gene.meta[] %<>% lapply(function(x) if (is.list(x)) sapply(x, str_c, collapse=",") else x)

## Fit models, correcting for prep batch effect
dge <- DGEList(counts=assay(sexp, "correct.counts"))
dge$genes <- all.gene.meta
dge %<>% calcNormFactors %>% .[aveLogCPM(.) > 0,] %>% estimateDisp(design, robust=TRUE)
dge2 <- offsetBatchEffect(dge, batch=sample.table$Batch,
                          design=model.matrix(~Donor + Celltype, sample.table))

## Dispersions with/without eBayes, with/without robust
dge.with.eBayes <- dge %>% estimateDisp(design, robust=FALSE)
dge.with.robust.eBayes <- dge %>% estimateDisp(design, robust=TRUE)
dge.without.eBayes <- dge %>% estimateDisp(design, prior.df=0)

disptable <- data.frame(
    logCPM=dge.without.eBayes$AveLogCPM,
    TrendBCV=dge.with.eBayes$trended.dispersion %>% sqrt,
    GeneWiseBCV=dge.without.eBayes$tagwise.dispersion %>% sqrt,
    eBayesBCV=dge.with.eBayes$tagwise.dispersion %>% sqrt,
    RobustBCV=dge.with.robust.eBayes$tagwise.dispersion %>% sqrt) +
    ylab("BCV")

raw.disp.plot <- ggplot(disptable) +
    aes(x=logCPM) +
    geom_point(aes(y=GeneWiseBCV), size=1.5, color="black") +
    geom_line(aes(y=TrendBCV), color="red") +
    ylab("BCV")

eBayes.disp.plot <- ggplot(disptable) +
    aes(x=logCPM) +
    geom_point(aes(y=GeneWiseBCV), size=1.5, color="gray") +
    geom_point(aes(y=eBayesBCV), size=1.5, color="darkblue") +
    geom_line(aes(y=TrendBCV), color="red") +
    ylab("BCV")

robust.eBayes.disp.plot <- ggplot(disptable) +
    aes(x=logCPM) +
    geom_point(aes(y=GeneWiseBCV), size=1.5, color="gray") +
    geom_point(aes(y=eBayesBCV), size=1.5, color="dodgerblue") +
    geom_point(aes(y=RobustBCV), size=1.5, color="darkgreen") +
    geom_line(aes(y=TrendBCV), color="red") +
    ylab("BCV")

png("results/RNA-seq/disp-plot-raw.png", res=300,
    height=8, width=8, units="in")
print(raw.disp.plot)
dev.off()
png("results/RNA-seq/disp-plot-eBayes.png", res=300,
    height=8, width=8, units="in")
print(eBayes.disp.plot)
dev.off()
png("results/RNA-seq/disp-plot-eBayes-robust.png", res=300,
    height=8, width=8, units="in")
print(robust.eBayes.disp.plot)
dev.off()

## Quality weight analysis
dbg <- estimateDispByGroup(dge2, sample.table$Group, sample.table$Batch)
elist.w <- voomWithQualityWeightsAndOffset(dge2, design, var.design=design)
covars <- sample.table %>% select(Group, Day, Donor, Batch, Celltype)
qcmetrics <- data.frame(Weights=elist.w$sample.weights,
                        GroupBCV=sapply(dbg, `[[`, "common.dispersion")[as.character(covars$Group)])
p <- ggpairs(cbind(covars, qcmetrics))
pdf("results/RNA-seq/qc-weights.pdf", width = 12, height=12)
print(p)
dev.off()

## Specific quality weight plots
awdf <- data.frame(covars, qcmetrics)
aw.plot.base <- ggplot(awdf) +
    aes(y=Weights) +
    geom_boxplot()

pdf("results/RNA-seq/weights-vs-covars.pdf")
aw.plot.base + aes(x=Group) + ggtitle("Array weights by Group")
aw.plot.base + aes(x=Celltype) + ggtitle("Array weights by Cell Type")
aw.plot.base + aes(x=Day) + ggtitle("Array weights by Time Point")
aw.plot.base + aes(x=Donor) + ggtitle("Array weights by Donor")
aw.plot.base + aes(x=Batch) + ggtitle("Array weights by Batch")
dev.off()

## Model fitting
gfit <- glmFit(dge2, design)
qfit <- glmQLFit(dge2, design, robust=TRUE)

elist <- voom(dge, design)
elist.bc <- voomWithOffset(dge2, design)

## Just for comparing voom %>% removeBatchEffect vs offsetBatchEffect %>% voom
elist.bc2 <- elist
elist.bc2$E %<>%
    removeBatchEffect(batch=sample.table$Batch,
                      design=model.matrix(~Donor + Celltype, sample.table))

fit <- lmFit(elist.bc, design)
fit.w <- lmFit(elist.w, design)

dmat <- suppressPlot(plotMDS(elist)$distance.matrix) %>% as.dist
mds <- cmdscale(dmat, k=attr(dmat, "Size") - 1, eig=TRUE)
mds$points %<>% add.numbered.colnames("Dim") %>% data.frame(sample.table, .)
dmat.bc <- suppressPlot(plotMDS(elist.bc)$distance.matrix) %>% as.dist
mds.bc <- cmdscale(dmat.bc, k=attr(dmat.bc, "Size") - 1, eig=TRUE)
mds.bc$points %<>% add.numbered.colnames("Dim") %>% data.frame(sample.table, .)

edger.dmat <- suppressPlot(plotMDS(cpm(dge, log=TRUE))$distance.matrix) %>% as.dist
edger.mds <- cmdscale(edger.dmat, k=attr(dmat, "Size") - 1, eig=TRUE)
edger.mds$points %<>% add.numbered.colnames("Dim") %>% data.frame(sample.table, .)
edger.dmat.bc <- suppressPlot(plotMDS(cpmWithOffset(dge2, log=TRUE))$distance.matrix) %>% as.dist
edger.mds.bc <- cmdscale(edger.dmat.bc, k=attr(dmat.bc, "Size") - 1, eig=TRUE)
edger.mds.bc$points %<>% add.numbered.colnames("Dim") %>% data.frame(sample.table, .)

pdf("results/RNA-seq/rnaseq-MDSPlots.pdf", width=12, height=12)
print(ggplot(mds$points) +
      aes(x=Dim1, y=Dim2, color=Batch, label=title) +
      geom_text() +
      scale_x_continuous(expand=c(0.15, 0)) +
      coord_equal() +
      ggtitle("limma voom PC1 & 2 of data (no batch correction)"))
print(ggplot(mds.bc$points) +
      aes(x=Dim1, y=Dim2, color=Batch, label=title) +
      geom_text() +
      scale_x_continuous(expand=c(0.15, 0)) +
      coord_equal() +
      ggtitle("limma voom PC1 & 2 of data (after batch correction)"))
print(ggplot(mds$points) +
      aes(x=Dim2, y=Dim3, color=Batch, label=title) +
      geom_text() +
      scale_x_continuous(expand=c(0.15, 0)) +
      coord_equal() +
      ggtitle("limma voom PC2 & 3 of data (no batch correction)"))
print(ggplot(edger.mds$points) +
      aes(x=Dim1, y=Dim2, color=Batch, label=title) +
      geom_text() +
      scale_x_continuous(expand=c(0.15, 0)) +
      coord_equal() +
      ggtitle("edgeR PC1 & 2 of data (no batch correction)"))
print(ggplot(edger.mds.bc$points) +
      aes(x=Dim1, y=Dim2, color=Batch, label=title) +
      geom_text() +
      scale_x_continuous(expand=c(0.15, 0)) +
      coord_equal() +
      ggtitle("edgeR PC1 & 2 of data (after batch correction)"))
print(ggplot(edger.mds$points) +
      aes(x=Dim2, y=Dim3, color=Batch, label=title) +
      geom_text() +
      scale_x_continuous(expand=c(0.15, 0)) +
      coord_equal() +
      ggtitle("edgeR PC2 & 3 of data (no batch correction)"))
dev.off()

## Test definitions
celltypes <- unique(sample.table$Celltype)
all.timepoints <- unique(sample.table$Day)
nonzero.timepoints <- setdiff(all.timepoints, "D0")

timepoint.anova.tests <- setNames(llply(celltypes, function(ct) {
    setNames(sprintf("%s%s - %sD0", ct, nonzero.timepoints, ct),
             sprintf("%s.D0v%s", ct, nonzero.timepoints))
}), nm=str_c(celltypes, ".AllT"))
timepoint.single.tests <- as.list(unlist(timepoint.anova.tests))
celltype.singlet.tests <-
    as.list(setNames(sprintf("Memory%s - Naive%s", all.timepoints, all.timepoints),
                     sprintf("NvM.%s", all.timepoints)))
celltype.allt.test <- list(NvM.AllT=unlist(celltype.singlet.tests))
factorial.singlet.tests <-
    as.list(setNames(sprintf("(Memory%s - MemoryD0) - (Naive%s - NaiveD0)",
                             nonzero.timepoints, nonzero.timepoints),
                     sprintf("Fac.%s", nonzero.timepoints)))
factorial.allt.test <- list(Fac.AllT=unlist(factorial.singlet.tests))
mi.vs.nf.test <- list(MD0vND14="MemoryD0 - NaiveD14")
donor.var.test <- list(InterDonor=sprintf("Donor%s", 1:3))
alltests <- c(timepoint.anova.tests, timepoint.single.tests,
              celltype.allt.test, celltype.singlet.tests,
              factorial.allt.test, factorial.singlet.tests,
              mi.vs.nf.test, donor.var.test)

limma.results <- lapply(
    alltests,
    . %>% { set_colnames(makeContrasts(contrasts=., levels=design),
                         str_c("logFC.", if (is.null(names(.))) seq_len(length(.)) else names(.))) } %>%
        contrasts.fit(fit, contrasts=.) %>%
        eBayes(robust=TRUE) %>%
        topTable( n=Inf, sort.by = "none") %>%
        .[!colnames(.) == "B"] %>%
        dplyr::rename(PValue=P.Value, FDR=adj.P.Val, logCPM=AveExpr) %>%
        merge(x=., y=all.gene.meta, all=TRUE) %>%
        select(ENTREZ, SYMBOL, GENENAME, PValue, FDR,
               logCPM, starts_with("logFC"),
               matches("^(F|t)$"), everything()))

limma.results.w <- lapply(
    alltests,
    . %>% { set_colnames(makeContrasts(contrasts=., levels=design),
                         str_c("logFC.", if (is.null(names(.))) seq_len(length(.)) else names(.))) } %>%
        contrasts.fit(fit.w, contrasts=.) %>%
        eBayes(robust=TRUE) %>%
        topTable( n=Inf, sort.by = "none") %>%
        .[!colnames(.) == "B"] %>%
        dplyr::rename(PValue=P.Value, FDR=adj.P.Val, logCPM=AveExpr) %>%
        merge(x=., y=all.gene.meta, all=TRUE) %>%
        select(ENTREZ, SYMBOL, GENENAME, PValue, FDR,
               logCPM, starts_with("logFC"),
               matches("^(F|t)$"), everything()))

edgeR.results <- lapply(
    alltests,
    . %>% { set_colnames(makeContrasts(contrasts=., levels=design), names(.)) } %>%
        glmQLFTest(glmfit=qfit, contrast=.) %>%
        topTags(n=Inf, sort.by = "none") %>%
        as.data.frame %>%
        merge(x=., y=all.gene.meta, all=TRUE) %>%
        select(ENTREZ, SYMBOL, GENENAME, PValue, FDR,
               logCPM, starts_with("logFC"),
               matches("^(F|t)$"), everything()))

limma.destats <- data_frame(
    Test=names(alltests),
    Percent_DE=limma.results %>%
        sapply(. %$% PValue %>% na.omit %>% propTrueNull) %>%
        subtract(1, .) %>% multiply_by(100),
    Percent_Sig=limma.results %>%
        sapply(. %$% FDR %>% na.omit %>% is_weakly_less_than(0.1) %>% {sum(.) / length(.)}) %>%
        multiply_by(100))

limma.destats.w <- data_frame(
    Test=names(alltests),
    Percent_DE=limma.results.w %>%
        sapply(. %$% PValue %>% na.omit %>% propTrueNull) %>%
        subtract(1, .) %>% multiply_by(100),
    Percent_Sig=limma.results.w %>%
        sapply(. %$% FDR %>% na.omit %>% is_weakly_less_than(0.1) %>% {sum(.) / length(.)}) %>%
        multiply_by(100))

edgeR.destats <- data_frame(
    Test=names(alltests),
    Percent_DE=edgeR.results %>%
        sapply(. %$% PValue %>% na.omit %>% propTrueNull) %>%
        subtract(1, .) %>% multiply_by(100),
    Percent_Sig=edgeR.results %>%
        sapply(. %$% FDR %>% na.omit %>% is_weakly_less_than(0.1) %>% {sum(.) / length(.)}) %>%
        multiply_by(100))

plot.pvals <- function(pvals, nullprop=propTrueNull(pvals)) {
    x <- data.frame(PValue=pvals)
    hlines <- c(Uniform=1, `Null Fraction`=nullprop) %>%
        data.frame(Intercept=.,
                   Line=names(.) %>% factor(levels=unique(.)))
    ggplot(x) +
        aes(x=PValue, y=..density..) + geom_histogram(binwidth=0.02) +
        geom_hline(aes(yintercept=Intercept, group=Line, color=Line, linetype=Line),
                   data=hlines, alpha=0.5, show_guide=TRUE) +
        theme(legend.title=element_blank(),
              legend.justification=c(1,1),
              legend.position=c(0.95,0.95),
              legend.background=element_rect(fill="gray90")) +
        xlim(0,1) +
        xlab("p-value") + ylab("Relative Density")
}

{
    p <- lapply(names(limma.results), function(i) {
        x <- limma.results[[i]]
        percent.sig <- limma.destats %>% filter(Test == i) %$% Percent_Sig
        percent.de <- limma.destats %>% filter(Test == i) %$% Percent_DE
        title <- sprintf(
            "(Unweighted) limma-voom P-value histogram for %s\n(%.3g%% FDR <= 0.1; Est. %.3g%% DE)",
            i, percent.sig, percent.de)
        plot.pvals(x$PValue, 1 - percent.de / 100) +
            ggtitle(title)
    })
    pdf("results/RNA-seq/rnaseq-limma-unweighted-pvalue-histograms.pdf")
    print(p)
    dev.off()
}

{
    p <- lapply(names(limma.results.w), function(i) {
        x <- limma.results.w[[i]]
        percent.sig <- limma.destats.w %>% filter(Test == i) %$% Percent_Sig
        percent.de <- limma.destats.w %>% filter(Test == i) %$% Percent_DE
        title <- sprintf(
            "limma-voom P-value histogram for %s\n(%.3g%% FDR <= 0.1; Est. %.3g%% DE)",
            i, percent.sig, percent.de)
        plot.pvals(x$PValue, 1 - percent.de / 100) +
            ggtitle(title)
    })
    pdf("results/RNA-seq/rnaseq-limma-pvalue-histograms.pdf")
    print(p)
    dev.off()
}

{
    p <- lapply(names(edgeR.results), function(i) {
        x <- edgeR.results[[i]]
        percent.sig <- edgeR.destats %>% filter(Test == i) %$% Percent_Sig
        percent.de <- edgeR.destats %>% filter(Test == i) %$% Percent_DE
        title <- sprintf(
            "edgeR P-value histogram for %s\n(%.3g%% FDR <= 0.1; Est. %.3g%% DE)",
            i, percent.sig, percent.de)
        plot.pvals(x$PValue, 1 - percent.de / 100) +
            ggtitle(title)
    })
    pdf("results/RNA-seq/rnaseq-edgeR-pvalue-histograms.pdf")
    print(p)
    dev.off()
}

{
    p <- lapply(names(alltests), function(i) {
        title <- sprintf("edgeR vs limma-voom P-values for %s", i)
        df <- data_frame(limma=limma.results.w[[i]]$PValue,
                         edgeR=edgeR.results[[i]]$PValue)
        ggplot(df) +
            aes(x=limma,
                y=edgeR,) +
            geom_point(size=0.5) +
            coord_equal() +
            geom_abline(slope=1, intercept=0, color="red", alpha=0.5, linetype=2) +
            scale_x_log10(name="limma voom P-value") +
            scale_y_log10(name="edgeR QLFTest P-value") +
            ggtitle(title)
    })
    pdf("results/RNA-seq/rnaseq-edgeR-vs-limma.pdf")
    print(p)
    dev.off()
}

{
    p <- lapply(names(alltests), function(i) {
        title <- sprintf("limma-voom weighted vs unweighted P-values for %s", i)
        df <- data_frame(Unweighted=limma.results[[i]]$PValue,
                         Weighted=limma.results.w[[i]]$PValue)
        ggplot(df) +
            aes(x=Unweighted,
                y=Weighted,) +
            geom_point(size=0.5) +
            coord_equal() +
            geom_abline(slope=1, intercept=0, color="red", alpha=0.5, linetype=2) +
            scale_x_log10(name="Unweighted limma voom P-value") +
            scale_y_log10(name="Weighted limma voom P-value") +
            ggtitle(title)
    })
    pdf("results/RNA-seq/rnaseq-limma-weighted-vs-uw.pdf")
    print(p)
    dev.off()
}

saveRDS(limma.results, "saved_data/rnaseq-limma-results-unweighted.RDS")
saveRDS(limma.results.w, "saved_data/rnaseq-limma-results.RDS")
saveRDS(edgeR.results, "saved_data/rnaseq-edgeR-results.RDS")

write.xlsx(limma.results %>%
           lapply(. %>% filter(!is.na(FDR) & FDR <= 0.1) %>% arrange(PValue)) %>%
           c(list(Stats=limma.destats), .),
           "results/RNA-seq/rnaseq-limma-results-unweighted.xlsx")
write.xlsx(limma.results.w %>%
           lapply(. %>% filter(!is.na(FDR) & FDR <= 0.1) %>% arrange(PValue)) %>%
           c(list(Stats=limma.destats.w), .),
           "results/RNA-seq/rnaseq-limma-results.xlsx")
write.xlsx(edgeR.results %>%
           lapply(. %>% filter(!is.na(FDR) & FDR <= 0.1) %>% arrange(PValue)) %>%
           c(list(Stats=edgeR.destats), .),
           "results/RNA-seq/rnaseq-edgeR-results.xlsx")

save.image("saved_data/rnaseq-analysis.rda")


## Minor breaks at every 1 logFC
maplot.minor.break.fun <- function(limits) {
    seq(from=ceiling(min(limits)),
        to=floor(max(limits)),
        by=1)
}

do.maplot <- function(tab) {
    tab %<>% arrange(PValue, logCPM)
    ggplot(tab) + aes(x=logCPM, y=logFC) +
        geom_hline(yintercept=0, alpha=0.5) +
        geom_point(aes(size=-log10(PValue), color=FDR <= 0.1)) +
        geom_density2d(show_guide=FALSE) +
        geom_smooth(color="red") +
        xlab("log2(CPM)") + ylab("log2(FC)") +
        scale_size(
            range=c(.1,1.5),
            trans="sqrt",
            name="-log10(PValue)") +
        scale_x_continuous(
            expand=c(0.05, 0),
            ## labels=function(x) parse(text=sprintf("2^%s", x)),
            minor_breaks=maplot.minor.break.fun) +
        scale_y_continuous(
            expand=c(0.05, 0),
            ## labels=function(x) parse(text=sprintf("2^%s", x)),
            minor_breaks=maplot.minor.break.fun,
            limits=c(-1, 1) * max(abs(tab$logFC))) +
        theme_bw() +
        theme(
            legend.position=c(1,1),
            legend.justification=c(1.1, 1.1))
}

single.logFC.tests <- alltests %>% {names(.)[elementLengths(.) == 1]}

{
    p <- lapply(single.logFC.tests, function(i) {
            do.maplot(limma.results[[i]]) +
                ggtitle(sprintf("MA plot for %s", i))
        })
    pdf("results/RNA-seq/rnaseq-maplots-limma-unweighted.pdf", width=8, height=8)
    suppressMessages(print(p))
    dev.off()
}

{
    p <- lapply(single.logFC.tests, function(i) {
        do.maplot(limma.results.w[[i]]) +
            ggtitle(sprintf("MA plot for %s", i))
    })
    pdf("results/RNA-seq/rnaseq-maplots-limma-arrayweights.pdf", width=8, height=8)
    suppressMessages(print(p))
    dev.off()
}

{
    p <- lapply(single.logFC.tests, function(i) {
        do.maplot(edgeR.results[[i]]) +
            ggtitle(sprintf("MA plot for %s", i))
    })
    pdf("results/RNA-seq/rnaseq-maplots-edger.pdf", width=8, height=8)
    suppressMessages(print(p))
    dev.off()
}
