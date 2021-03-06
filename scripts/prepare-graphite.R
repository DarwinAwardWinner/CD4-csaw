#!/usr/bin/Rscript
library(graphite)
library(here)
library(glue)
library(BiocParallel)
library(rctutils)

use_futures("multicore")

target.species <- "hsapiens"
dbnames <- pathwayDatabases() %>%
    filter(species == target.species) %$% database %>%
    as.character

graphite.species <- "hsapiens"
dbs.raw <- dbnames %>% set_names %>%
    lapply(pathways, species = graphite.species)
idtypes <- c(entrez = "entrez", ensembl = "ENSEMBL", symbol = "symbol")
dbs <- suppressMessages(lapply(idtypes, . %>% bplapply(dbs.raw, convertIdentifiers, to = .)))

for (i in names(dbs)) {
    fname <- here("saved_data", glue("graphite-{i}.RDS"))
    saveRDS(dbs[[i]], fname)
}

spia.outdir <- here("saved_data", "SPIA")
dir.create(spia.outdir, recursive = TRUE, showWarnings = FALSE)

## Prepare SPIA databases
for (idtype in names(dbs)) {
    db <- dbs[[idtype]]
    bpmapply(prepareSPIA,
             db = db,
             pathwaySetName = file.path(spia.outdir, glue("graphite-{idtype}-{names(db)}Ex")))
}
