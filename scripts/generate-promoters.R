#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(getopt)
    library(optparse)
    library(stringr)
    library(glue)
    library(magrittr)
    library(assertthat)
    library(GenomicRanges)
    library(rctutils)
})

get_options <- function(opts) {

    ## Do argument parsing early so the script exits quickly if arguments are invalid
    optlist <- list(
        ## So far this script only supports TxDb objects because
        ## figuring out the first exon and TSS from other
        ## less-structured formats is a pain.
        make_option(c("-t", "--annotation-txdb"), metavar = "TXDBNAME", type = "character",
                    help = "Name of TxDb package, or the name of a database file, to use for gene annotation"),
        make_option(c("-r", "--promoter-radius"), metavar = "RADIUS", type = "character",
                    help = "Maximum distance from a gene's transcription start site that is considered part of the promoter."),
        make_option(c("-o", "--output-file"), metavar = "FILENAME.RDS", type = "character",
                    help = "Output file name. The GRanges object containing the promoter regions will be saved here using saveRDS, so it should end in '.RDS'."),
        make_option(c("-a", "--additional-gene-info"), metavar = "FILENAME", type = "character",
                    help = "RDS/RData/xlsx/csv file containing a table of gene metadata. Row names (or the first column of the file if there are no row names) should be gene/feature IDs that match the ones used in the main annotation, and these should be unique. This can also be a GFF3 file where the metadata is in the attributes of elements of type 'gene', where the 'ID' attribute specifies the gene ID."))
    progname <- na.omit(c(get_Rscript_filename(), "generate-promoter.R"))[1]
    parser <- OptionParser(
        usage = "Usage: %prog [options] -t TXDB -r RADIUS -o OUTPUT.RDS",
        description = "Generate promoter regions for a gene annotation.

From each transcription start site, a region is extended to the specified radius in both directions to define the promoter region for that TSS. Then, any overlapping promoters that share a gene ID are merged. Note that the number of non-overlapping promoters necessarily depends on the chosen promoter radius, and different radii will give different numbers of promoters. The resulting GRanges object will have several annotation columns added: 'PromoterID', which is an arbitrary unique key for each promoter; 'GeneID', which may be NA, may be equal to the TxID for transcripts with no Gene ID, and may be non-unique if a single gene has multiple non-overlapping promoters; TxID, a CharacterList of all transcript IDs whose TSS are included in the promoter.",
option_list = optlist,
add_help_option = TRUE,
prog = progname,
epilogue = "Note that all base pair sizes (for the promoter radius) may have a suffix of 'bp', 'kbp', 'mbp', or 'tbp' (and the 'p' is optional). For example, 10kbp = 10000")

    cmdopts <- parse_args(parser, opts)
    ## Ensure that all required arguments were provided
    required.opts <- c("annotation-txdb", "output-file", "promoter-radius")
    missing.opts <- setdiff(required.opts, names(cmdopts))
    if (length(missing.opts) > 0) {
        stop(str_c("Missing required arguments: ", deparse(missing.opts)))
    }

    ## Convert bp args to numbers
    for (i in c("promoter-radius")) {
        cmdopts[[i]] %<>% parse_bp
    }

    ## Replace dashes with underscores so that all options can easily
    ## be accessed by "$"
    cmdopts %>% setNames(chartr("-", "_", names(.)))
}

## Do this early to handle "--help" before wasting time loading
## pacakges & stuff
invisible(get_options(commandArgs(TRUE)))

{
    cmdopts <- get_options(commandArgs(TRUE))
    cmdopts$help <- NULL

    tsmsg("Args:")
    print_var_vector(cmdopts)

    ## Delete the output file if it exists
    suppressWarnings(file.remove(cmdopts$output_file))
    assert_that(!file.exists(cmdopts$output_file))

    ## Only chr1-chr22,chrX,chrY
    std.chr <- extractSeqlevels("Homo sapiens", "UCSC") %>% setdiff("chrM")

    tsmsg("Reading annotation data")
    txdb <- get_txdb(cmdopts$annotation_txdb)
    tsmsg(glue("Getting {format_bp(cmdopts$promoter_radius)}-radius promoters"))
    all_promoters <- suppressWarnings(promoters(txdb, upstream = cmdopts$promoter_radius, downstream = cmdopts$promoter_radius)) %>%
        trim %>% keepSeqlevels(std.chr, pruning.mode = "coarse")

    tsmsg("Annotating promoters")
    mcols(all_promoters) %<>%
        ## Not using dplyr because it's a BioC DataFrame
        transform(TxID = tx_name,
                  tx_name = NULL,
                  tx_id = NULL)
    all_promoters$GeneID <- suppressMessages(mapIds(txdb, all_promoters$TxID, keytype = "TXNAME", column = "GENEID", multiVals = "first")) %>%
        ifelse(is.na(.), all_promoters$TxID, .)
    tsmsg("Splitting promoters by gene ID")
    gene_promoters <- split(all_promoters, all_promoters$GeneID)
    assert_that(all(lengths(gene_promoters) >= 1))
    tsmsg("Merging overlapping promoters from the same gene")
    merged_promoters <- bplapply(gene_promoters, function(gp) {
        gp_reduced <- reduce(gp)
        gp_reduced$GeneID <- gp$GeneID[1]
        names(gp_reduced) <- gp_reduced$PromoterID <-
            glue("{gene}-P{pnum}", gene = gp_reduced$GeneID, pnum = seq_along(gp_reduced))
        pgroup <- gp_reduced$PromoterID[nearest(gp, gp_reduced)]
        gp_reduced$TxID <- CharacterList(split(gp$TxID, pgroup))[gp_reduced$PromoterID]
        gp_reduced
    }) %>% unname %>% GRangesList %>% unlist

    if ("additional_gene_info" %in% names(cmdopts)) {
        tsmsg("Reading additional gene annotation metadata")
        additional_gene_info <- read_additional_gene_info(cmdopts$additional_gene_info)
        genes_without_info <- setdiff(names(gene_promoters), rownames(additional_gene_info))
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
        assert_that(all(merged_promoters$GeneID %in% rownames(additional_gene_info)))
        mcols(merged_promoters)[colnames(additional_gene_info)] <- additional_gene_info[merged_promoters$GeneID,]
        metadata(merged_promoters) %<>% c(metadata(additional_gene_info))
    }

    tsmsg("Saving output file")
    save_RDS_or_RDA(merged_promoters, cmdopts$output_file)
    invisible(NULL)
}
