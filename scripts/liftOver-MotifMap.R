library(rctutils)

liftOver_motifMap(infile=snakemake@input[["bed"]], chainfile=snakemake@input[["chain"]],
                  outfile=snakemake@output[["bed"]], allow.gap=2)
