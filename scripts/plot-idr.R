#!/usr/bin/env Rscript

# This is a script to reproduce the plots of
# https://github.com/nboley/idr using ggplot2.

library(getopt)
library(optparse)

get_options <- function(opts) {

    ## Do argument parsing early so the script exits quickly if arguments are invalid
    optlist <- list(
        make_option(c("-i", "--idr-file"), metavar = "FILE", type = "character",
                    help = "(REQUIRED) Output file from idr script containing computed IDR values."),
        make_option(c("-o", "--output-file"), metavar = "PLOTS.PDF", type = "character",
                    help = "(REQUIRED) PDF file in which to create plots"),
        make_option(c("-A", "--sample-A-name"), type = "character", default = "Sample A",
                    help = "Name of the first sample. This will be used to identify the sample in plot titles and axis labels."),
        make_option(c("-B", "--sample-B-name"), type = "character", default = "Sample B",
                    help = "Name of the second sample. This will be used to identify the sample in plot titles and axis labels."),
        make_option(c("-P", "--sample-name-common-prefix"), type = "character",
                    help = "Common prefix to both sample names. Used to reduced redundancy in plot titles. Will be stripped from given sample names if present."),
        make_option(c("-s", "--prefix-separator"), type = "character", default = "._-/:| ",
                    help = "Characters that will be stripped from the split point between prefix and sample name."))
    progname <- na.omit(c(get_Rscript_filename(), "plot-idr.R"))[1]
    parser <- OptionParser(
        usage = "Usage: %prog [options] -i FILE -o PLOTS.PDF [ -A \"Sample A\" -B \"Sample B\" ]",
        description = "Generate QC plots for an IDR analysis of two samples.",
        option_list = optlist,
        add_help_option = TRUE,
        prog = progname,
        epilogue = "")

    cmdopts <- parse_args(parser, opts)
    ## Ensure that all required arguments were provided
    required.opts <- c("idr-file", "output-file")
    missing.opts <- setdiff(required.opts, names(cmdopts))
    if (length(missing.opts) > 0) {
        stop(str_c("Missing required arguments: ", deparse(missing.opts)))
    }
    cmdopts %>% setNames(str_replace_all(names(.), "-", "_"))
}

## Do this early to handle "--help" before wasting time loading
## pacakges & stuff
get_options(commandArgs(TRUE))

library(magrittr)
library(dplyr)
library(ggplot2)
library(scales)
library(ks)
library(reshape2)
library(stringr)
library(rex)
library(glue)
library(rctutils)

cutIDR <- function(x, thresholds = c(0.01, 0.05, 0.1)) {
    fullbreaks <- c(0, thresholds, 1)
    labels <- c(glue(" <= {thresholds}"), glue(" > {tail(thresholds, 1)}"))
    cut(x, breaks = fullbreaks, labels = labels) %>%
        factor(levels = rev(levels(.)))
}

{
    cmdopts <- get_options(commandArgs(TRUE))
    ## myargs <- c("-i", "idr_analysis/epic_hg38.analysisSet/H3K4me3_condition.ALL_D4659vsD5053/idrValues.txt",
    ##             "-o", "idr_analysis/epic_hg38.analysisSet/H3K4me3_condition.ALL_D4659vsD5053/idrplots.pdf",
    ##             "-A", "H3K4me3_ALL_D4659", "-B", "H3K4me3_ALL_D5053",
    ##             "-P", "H3K4me3_ALL")
    ## cmdopts <- get_options(myargs)
    cmdopts$help <- NULL

    tsmsg("Args:")
    print_var_vector(cmdopts)

    if (!is.null(cmdopts$sample_name_common_prefix)) {
        prefix.rx <- rex(start, cmdopts$sample_name_common_prefix)
        cmdopts$sample_A_name %<>% str_replace(prefix.rx, "")
        cmdopts$sample_B_name %<>% str_replace(prefix.rx, "")

        if (!is.null(cmdopts$prefix_separator) &&
            str_length(cmdopts$prefix_separator) > 0) {
            post.sep.rx <- rex(some_of(cmdopts$prefix_separator), end)
            pre.sep.rx <- rex(start, some_of(cmdopts$prefix_separator))

            cmdopts$sample_name_common_prefix %<>%
                str_replace(post.sep.rx, "")
            cmdopts$sample_A_name %<>% str_replace(pre.sep.rx, "")
            cmdopts$sample_B_name %<>% str_replace(pre.sep.rx, "")
        }
    }

    idrtab <- read_idr_table(cmdopts$idr_file)

    idrtab %<>%
        mutate(GlobalIDR.Cut = cutIDR(GlobalIDR),
               LocalIDR.Cut = cutIDR(LocalIDR),
               rankA = rank(scoreA),
               rankB = rank(scoreB),
               rankBinA = ceiling(rankA / length(rankA) * 20) %>% factor,
               rankBinB = ceiling(rankB / length(rankB) * 20) %>% factor)

    title_samples <- str_interp("${sample_A_name} vs ${sample_B_name}", cmdopts)
    title_sampleA <- str_interp("${sample_A_name}", cmdopts)
    title_sampleB <- str_interp("${sample_B_name}", cmdopts)
    if (!is.null(cmdopts$sample_name_common_prefix)) {
        title_samples <- str_c(cmdopts$sample_name_common_prefix, ", ", title_samples)
        title_sampleA <- str_c(cmdopts$sample_name_common_prefix, " ", title_sampleA)
        title_sampleB <- str_c(cmdopts$sample_name_common_prefix, " ", title_sampleB)
    }

    plotlist <- list(
        RankCons = ggplot(idrtab %>% arrange(desc(GlobalIDR))) +
            aes(x = rankA, y = rankB,
                color = GlobalIDR.Cut) +
            geom_point() +
            scale_color_manual(name = "IDR", values = discrete_gradient(nlevels(idrtab$GlobalIDR.Cut))) +
            coord_fixed() +
            theme(legend.position = "bottom") +
            xlab(str_interp("Peak Rank in ${cmdopts$sample_A_name}")) +
            ylab(str_interp("Peak Rank in ${cmdopts$sample_B_name}")) +
            ggtitle(str_interp("Rank consistency plot for ${title_samples}")),
        ScoreCons = ggplot(idrtab %>% arrange(desc(GlobalIDR))) +
            aes(x = scoreA, y = scoreB,
                color = cutIDR(GlobalIDR)) +
            geom_point() +
            scale_x_log10() + scale_y_log10() +
            scale_color_manual(name = "IDR", values = discrete_gradient(nlevels(idrtab$GlobalIDR.Cut))) +
            coord_fixed() +
            theme(legend.position = "bottom") +
            xlab(str_interp("Peak Score in ${cmdopts$sample_A_name}")) +
            ylab(str_interp("Peak Score in ${cmdopts$sample_B_name}")) +
            ggtitle(str_interp("Score consistency plot for ${title_samples}")))

    neglog10_trans <- neglog_trans(10)

    ## Sample B rank vs IDR
    pointdata <- idrtab %>% transmute(x = rankA, y = -log10(GlobalIDR))
    H <- pointdata %>% Hbcv.diag(binned = TRUE)
    k <- pointdata %>%
        as.matrix %>%
        kde(gridsize = 1024, bgridsize = rep(1024, 2), verbose = TRUE,
            H = H/8, binned = TRUE)
    ## Sometimes the estimate goes a bit negative, which is no good
    densdata <- melt(k$estimate) %>%
        transmute(
            x = k$eval.points[[1]][Var1],
            y = k$eval.points[[2]][Var2],
            Density = value %>% pmax(0),
            ## Part of a hack to make the alpha look less bad
            AlphaDens = value %>% pmax(1e-15))

    plotlist$DensA <-
        ggplot(densdata) +
        aes(x = x, y = neglog10_trans$inverse(y), alpha = Density) +
        geom_raster(fill = muted("blue",c = 90), interpolate = TRUE) +
        scale_alpha(limits = c(0, max(densdata$Density)/3), range = c(0,1), guide = FALSE) +
        theme(legend.position = "bottom") +
        coord_cartesian(expand = FALSE) +
        xlab(str_interp("Peak Rank in ${cmdopts$sample_A_name}")) +
        scale_y_continuous(name = "IDR", trans = neglog10_trans) +
        ggtitle(str_interp("IDR vs Peak Rank for ${title_sampleA}"))

    ## Sample B rank vs IDR
    pointdata <- idrtab %>% transmute(x = rankB, y = -log10(GlobalIDR))
    H <- pointdata %>% Hbcv.diag(binned = TRUE)
    k <- pointdata %>%
        as.matrix %>%
        kde(gridsize = 1024, bgridsize = rep(1024, 2), verbose = TRUE,
            H = H/8, binned = TRUE)
    ## Sometimes the estimate goes a bit negative, which is no good
    densdata <- melt(k$estimate) %>%
        transmute(
            x = k$eval.points[[1]][Var1],
            y = k$eval.points[[2]][Var2],
            Density = value %>% pmax(0),
            ## Part of a hack to make the alpha look less bad
            AlphaDens = value %>% pmax(1e-15))

    plotlist$DensB <-
        ggplot(densdata) +
        aes(x = x, y = neglog10_trans$inverse(y), alpha = Density) +
        geom_raster(fill = muted("blue",c = 90), interpolate = TRUE) +
        scale_alpha(limits = c(0, max(densdata$Density)/3), range = c(0,1), guide = FALSE) +
        theme_bw() +
        coord_cartesian(expand = FALSE) +
        xlab(str_interp("Peak Rank in ${cmdopts$sample_B_name}")) +
        scale_y_continuous(name = "IDR", trans = neglog10_trans) +
        ggtitle(str_interp("IDR vs Peak Rank for ${title_sampleB}"))

    pdf(cmdopts$output_file)
    print(plotlist)
    dev.off()
}
