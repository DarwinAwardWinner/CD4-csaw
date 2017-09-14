library(here)
source(here("scripts", "utilities.R"))

liftOver_motifMap <- function(infile, chainfile, outfile, allow.gap=2) {
    gr <- read.motifmap(infile, parse_name=TRUE)
    chain <- import.chain(chainfile)
    gr2 <- liftOverLax(gr, chain, allow.gap=allow.gap)
    gr2 <- unlist(gr2[lengths(gr2) == 1])
    write.motifmap(gr2, outfile)
}

liftOver_motifMap(infile=snakemake@input[["bed"]], chainfile=snakemake@input[["chain"]],
                  outfile=snakemake@output[["bed"]], allow.gap=2)
