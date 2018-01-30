suppressMessages({
    library(rtracklayer)
    library(stringr)
    library(assertthat)
    library(dplyr)
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(here)
    source(here("scripts", "utilities.R"))
})

{
    chainfile <- snakemake@input[["chain"]]
    chain <- import.chain(chainfile)

    outfile <- snakemake@output[[1]]
    assert_that(is.character(outfile))

    mySession <- browserSession()
    genome(mySession) <- "hg19"
    sites.table <- getTable(ucscTableQuery(mySession, track="tfbsConsSites", table="tfbsConsSites")) %>%
        mutate_if(is.factor, as.character)
    names.table <- getTable(ucscTableQuery(mySession, track="tfbsConsSites", table="tfbsConsFactors")) %>%
        mutate_if(is.factor, as.character)
    ## Keep only human entries, interpret "N" as NA
    names.table$id[names.table$id == "N"] <- NA
    names.table %<>% filter(species=="human") %>% select(-species) %>% droplevels %>%
        group_by(name) %>% summarize_all(. %>% str_c(collapse=","))
    assert_that(!anyDuplicated(names.table$name))
    full.table <- sites.table %>% inner_join(names.table, "name")

    gr <- makeGRangesFromDataFrame(full.table, start.field="chromStart", end.field="chromEnd",
                                   starts.in.df.are.0based=TRUE, keep.extra.columns=TRUE,
                                   seqinfo=seqinfo(BSgenome.Hsapiens.UCSC.hg19))
    gr.lifted <- liftOverLax(gr, chain, allow.gap=2) %>% .[lengths(.) == 1] %>% unlist
    saveRDS(gr.lifted, outfile)
}
