library(GSEABase)
library(BiocParallel)
library(org.Hs.eg.db)
library(rtracklayer)
library(glue)
library(tidyverse)
library(assertthat)

library(doParallel)
ncores <- getOption("mc.cores", default=parallel::detectCores(logical = FALSE))
registerDoParallel(cores=ncores)
library(BiocParallel)
register(DoparParam())

quotemeta <- function (string) {
  str_replace_all(string, "(\\W)", "\\\\\\1")
}

split.by.category.and.subcategory <- function(x, cat, subcat) {
    cat <- droplevels(as.factor(cat))
    subcat <- droplevels(as.factor(subcat))
    subcat.by.cat <- lapply(split(subcat, cat), droplevels)
    res <- split(x, cat)
    for (ci in levels(cat)) {
        cat.subcat <- subcat.by.cat[[ci]]
        for (sci in levels(cat.subcat)) {
            scname <- str_c(ci, ".", sci)
            ## This regexp stuff makes sure that "CP:KEGG" etc. are
            ## included in the "CP" subcategory
            scregexp <- regex(str_c("^", quotemeta(sci), "(\\:|$)"))
            res[[scname]] <- res[[ci]][!is.na(cat.subcat) & str_detect(cat.subcat, scregexp)]
        }
    }
    SimpleList(res)
}

extract.bset.metadata <- function(x) {
    ## Extract relevant metadata about each gene set and save it in a data frame
    bset.urls <- CharacterList(lapply(x, function(a) grep("^file:", urls(a), perl=TRUE, value=TRUE, invert=TRUE)))
    bset.urls <- rtracklayer:::pasteCollapse(bset.urls)
    bset.urls[bset.urls == ""] <- NA
    assert_that(all(lengths(bset.urls) == 1))
    bsets.meta <- data.frame(row.names=names(x), setName=names(x),
                             category=categories, subCategory=subcategories,
                             setIdentifier=unlist(lapply(x, setIdentifier)),
                             contributor=unlist(lapply(x, contributor)),
                             description=unlist(lapply(x, description)),
                             url=bset.urls)
    bsets.meta
}

## This just loads all the msigdb gene sets into R and saves them as an
## RDS file. Or it reads that RDS file if it already exists.
{
    ## Load MSigDB
    if (file.exists("saved_data/msigdb_v6.1.RDS")) {
        bsets_symbol <- readRDS("saved_data/msigdb_v6.1.RDS")
    } else {
        bsets_symbol <- getBroadSets("saved_data/msigdb_v6.1.xml")
        saveRDS(bsets_symbol, "saved_data/msigdb_v6.1.RDS")
    }
    id_converters <- list(symbol=identity,
                          entrez = . %>% mapIdentifiers(EntrezIdentifier("org.Hs.eg")),
                          ensembl = . %>% mapIdentifiers(ENSEMBLIdentifier("org.Hs.eg")))
    bsets <- bplapply(id_converters, function(converter) converter(bsets_symbol))
    sids <- unlist(lapply(as.list(bsets_symbol), setIdentifier))
    # Need to tell the converted gene sets that their gene IDs have
    # been converted
    for (idtype in setdiff(names(bsets), "symbol")) {
        message("Fixing idtype ", idtype)
        bsetlist <- as.list(bsets[[idtype]])
        bsetlist <- bplapply(seq_along(bsetlist), function(i) {
            x <- bsetlist[[i]]
            setIdentifier(x) <- sids[i]
            x
        })
        bsets[[idtype]] <- GeneSetCollection(bsetlist)
    }
    ## Split into categories
    ctypes <- lapply(bsets_symbol, collectionType)
    categories <- factor(unlist(lapply(ctypes, bcCategory)))
    subcategories <- factor(unlist(lapply(ctypes, bcSubCategory)))
    all.collections <- lapply(bsets, split.by.category.and.subcategory,
                             cat=categories, subcat=subcategories)
    for (i in names(all.collections)) {
        fname <- glue("saved_data/msigdb-{i}.RDS")
        saveRDS(all.collections[[i]], fname)
    }

    bsets.meta <- extract.bset.metadata(bsets_symbol)
    ## Munge some overly verbose descriptions (Note: These may no
    ## longer be applicable, but running them anyway won't hurt.)
    bsets.meta$description %<>%
        str_replace(
            pattern="Genes with promoter regions \\[-2kb,2kb\\] around transcription start site containing motif (\\w+)\\. Motif does not match any known transcription factor",
            replacement="Motif \\1 (no known TF) in gene promoter (2kb radius)") %>%
        str_replace(
            pattern="Genes with promoter regions \\[-2kb,2kb\\] around transcription start site containing the motif (\\w+) which matches annotation for (.+)",
            replacement="Motif \\1 (matches \\2) in gene promoter (2kb radius)") %>%
        str_replace(
            pattern="^Genes involved in ",
            replacement="") %>%
        str_replace(
            pattern="^Genes annotated by the GO term ",
            replacement="")
    saveRDS(bsets.meta, "saved_data/msigdb-meta.RDS")
}
