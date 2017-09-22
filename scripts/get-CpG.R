suppressMessages({
    library(rtracklayer)
    library(stringr)
    library(assertthat)
    library(BSgenome.Hsapiens.UCSC.hg38)
})

{
    outfile <- snakemake@output[[1]]
    assert_that(is.character(outfile))

    mySession <- browserSession()
    genome(mySession) <- "hg38"
    tab <- getTable(ucscTableQuery(mySession, "cpgIslandExtUnmasked"))
    gr <- makeGRangesFromDataFrame(tab, start.field="chromStart", end.field="chromEnd",
                                   starts.in.df.are.0based=TRUE, keep.extra.columns=TRUE,
                                   seqinfo=seqinfo(BSgenome.Hsapiens.UCSC.hg38))
    ## GRanges already knows the length of each feature, so this field is
    ## redundant.
    assert_that(all(width(gr) == gr$length))
    mcols(gr)$length <- NULL
    seqinfo(gr) <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

    saveRDS(gr, outfile)
}
