#!/usr/bin/env Rscript

library(getopt)
library(optparse)
library(assertthat)
library(rex)
library(rctutils)

get_options <- function(opts) {

    ## Do argument parsing early so the script exits quickly if
    ## arguments are invalid
    optlist <- list(
        make_option(c("-q", "--transcript-quant"), metavar = "SEXP.RDS", type = "character",
                    help = "File name of an R data file containing a RangedSummarizedExperiment of transcript abundances. These will be used to select the highest-expressed TSS for each gene."),
        ## So far this script only supports TxDb objects because
        ## figuring out the first exon and TSS from other
        ## less-structured formats is a pain.
        make_option(c("-t", "--annotation-txdb"), metavar = "TXDBNAME", type = "character",
                    help = "Name of TxDb package, or the name of a database file, to use for gene annotation"),
        make_option(c("-o", "--output-file"), metavar = "FILENAME.RDS", type = "character",
                    help = "Output file name. The GRanges object containing the promoter regions will be saved here using saveRDS, so it should end in '.RDS'."),
        make_option(c("-a", "--additional-gene-info"), metavar = "FILENAME", type = "character",
                    help = "RDS/RData/xlsx/csv file containing a table of gene metadata. Row names (or the first column of the file if there are no row names) should be gene/feature IDs that match the ones used in the main annotation, and these should be unique. This can also be a GFF3 file where the metadata is in the attributes of elements of type 'gene', where the 'ID' attribute specifies the gene ID."))
    progname <- na.omit(c(get_Rscript_filename(), "rnaseq-count.R"))[1]
    parser <- OptionParser(
        usage = "Usage: %prog [options] -q SEXP.RDS -t TXDB -o OUTPUT.RDS",
        description = "Select the most abundant TSS for each gene.

For each gene, transcripts are grouped by TSS, and their average abundances are added up. The TSS with the largest sum of average transcript abundances is selected as the representative TSS for that gene. These are all stored in a GRanges object in the output file. The resulting GRanges object will be annotated with a GeneID column. For transcripts with no associaated Gene ID, the GeneID column will be identical to the TxID column. Since a single TSS is being chosen for each gene, the GeneID column should not contain any duplicates.",
option_list = optlist,
add_help_option = TRUE,
prog = progname)

    cmdopts <- parse_args(parser, opts)
    ## Ensure that all required arguments were provided
    required.opts <- c("annotation-txdb", "output-file", "transcript-quant")
    missing.opts <- setdiff(required.opts, names(cmdopts))
    if (length(missing.opts) > 0) {
        stop(str_c("Missing required arguments: ", deparse_onestring(missing.opts)))
    }

    ## Replace dashes with underscores so that all options can easily
    ## be accessed by "$"
    cmdopts %>% setNames(chartr("-", "_", names(.)))
}

## Do this early to handle "--help" before wasting time loading
## pacakges & stuff
invisible(get_options(commandArgs(TRUE)))

library(assertthat)
library(dplyr)
library(magrittr)
library(stringr)
library(glue)
library(future)
library(SummarizedExperiment)

{
    cmdopts <- get_options(commandArgs(TRUE))
    cmdopts$help <- NULL

    ## ## For testing only
    ## cmdopts <- list(
    ##     transcript_quant = "saved_data/SummarizedExperiment_rnaseq_transcript_shoal_hg38.analysisSet_knownGene.RDS",
    ##     annotation_txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    ##     additional_gene_info = "/home/ryan/references/hg38/genemeta.org.Hs.eg.db.RDS",
    ##     output_file = "test.rds")

    tsmsg("Args:")
    print_var_vector(cmdopts)

    ## Delete the output file if it exists
    suppressWarnings(file.remove(cmdopts$output_file))
    assert_that(!file.exists(cmdopts$output_file))

    ## Only chr1-chr22,chrX,chrY
    std.chr <- extractSeqlevels("Homo sapiens", "UCSC") %>% setdiff("chrM")

    tsmsg("Reading quantification data")
    sexp <- readRDS(cmdopts$transcript_quant)
    sexp %<>% keepSeqlevels(std.chr, pruning.mode = "coarse")

    tsmsg("Reading annotation data")
    txdb <- get_txdb(cmdopts$annotation_txdb)

    tsmsg("Computing average transcript abundances")
    tx <- rowRanges(sexp)
    tx$abundance <- sexp %>% assay("abundance") %>% rowMeans
    tx$GeneID <- mapIds(txdb, names(tx),  keytype = "TXNAME", column = "GENEID", multiVals = "first")

    tsmsg("Grouping transcripts by TSS and gene ID")
    tss_table <- tx %>% promoters(upstream = 0, downstream = 1) %>% as("data.frame") %>%
        filter(!is.na(GeneID)) %>%
        group_by(GeneID, seqnames, start, end, strand) %>%
        summarize(transcript = str_c(transcript, collapse = ","),
                  abundance = sum(abundance))

    tsmsg("Selecting most abundant TSS for each gene")
    abundant_tss <- tss_table %>%
        arrange(desc(abundance)) %>% filter(!duplicated(GeneID)) %>%
        as("GRanges") %>% setNames(.$GeneID) %>%
        .[unique(tss_table$GeneID)]

    if ("additional_gene_info" %in% names(cmdopts)) {
        tsmsg("Reading additional gene annotation metadata")
        additional_gene_info <- read_additional_gene_info(cmdopts$additional_gene_info)
        ## Generate empty rows for genes that don't have additional
        ## info
        genes_without_info <- setdiff(abundant_tss$GeneID, rownames(additional_gene_info))
        if (length(genes_without_info) > 0) {
            empty_row <- list(character(0)) %>% rep(ncol(additional_gene_info)) %>% setNames(colnames(additional_gene_info))
            single.val.cols <- sapply(additional_gene_info, function(x) all(lengths(x) == 1))
            for (i in seq_along(empty_row)) {
                if (single.val.cols[i]) {
                    empty_row[[i]] <- NA
                } else {
                    empty_row[[i]] <- list(logical(0)) %>% as(class(additional_gene_info[[i]]))
                }
            }
            empty_row %<>% DataFrame
            empty_gene_table <- empty_row[rep(1, length(genes_without_info)),] %>%
                set_rownames(genes_without_info)
            additional_gene_info %<>% rbind(empty_gene_table)
        }
        assert_that(all(abundant_tss$GeneID %in% rownames(additional_gene_info)))
        mcols(abundant_tss)[colnames(additional_gene_info)] <- additional_gene_info[abundant_tss$GeneID,]
        metadata(abundant_tss) %<>% c(metadata(additional_gene_info))
    }

    tsmsg("Saving output file")
    save.RDS.or.RDA(abundant_tss, cmdopts$output_file)
    invisible(NULL)
}
