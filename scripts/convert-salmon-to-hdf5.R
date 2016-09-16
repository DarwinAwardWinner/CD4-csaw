#!/usr/bin/env Rscript

## Copy function, giving the copy a new empty environment whose parent
## is the original function's environment
copy_func <- function(func) {
    newfunc <- function () {}
    formals(newfunc) <- formals(func)
    body(newfunc) <- body(func)
    environment(newfunc) <- new.env(parent=environment(func))
    newfunc
}

## Fixed copy of fish_to_hdf5. See:
## https://github.com/COMBINE-lab/wasabi/issues/4
fish_to_hdf5 <- copy_func(wasabi::fish_to_hdf5)
body(fish_to_hdf5) <- expression({
    h5file <- file.path(fish_dir, "abundance.h5")
    if (!force && file.exists(h5file)) {
        print(paste("Skipping conversion: abundance.h5 already in ",
                    fish_dir))
        return()
    }
    if (file.exists(h5file)) {
        file.remove(h5file)
    }
    quant <- fread(file.path(fish_dir, "quant.sf"))
    setnames(quant, c("target_id", "length", "eff_length", "tpm",
                      "est_counts"))
    if (!is.character(quant$target_id)) {
        quant$target_id <- as.character(quant$target_id)
    }
    aux_dir = "aux"
    cmd_info <- rjson::fromJSON(file = file.path(fish_dir, "cmd_info.json"))
    if (is.element("auxDir", names(cmd_info))) {
        aux_dir = cmd_info$auxDir
    }
    auxPath <- file.path(fish_dir, aux_dir)
    minfo <- rjson::fromJSON(file = file.path(auxPath, "meta_info.json"))
    numBoot <- minfo$num_bootstraps
    if (numBoot > 0) {
        bootCon <- gzcon(file(file.path(auxPath, "bootstrap",
                                        "bootstraps.gz"), "rb"))
        dtype <- switch(minfo$samp_type, gibbs="integer", bootstrap="double",
                        stop(sprintf("Unknown sample type: %s", deparse(minfo$samp_type))))
        boots <- readBin(bootCon, dtype, n = minfo$num_targets *
                                             minfo$num_bootstraps)

        close(bootCon)
        boots <- as.numeric(boots)
        dim(boots) <- c(minfo$num_targets, minfo$num_bootstraps)
    }
    numProcessed <- minfo$num_processed
    if (is.element("mapping_type", names(minfo))) {
        if (minfo$mapping_type == "alignment") {
            if (fallback_num_reads == -1) {
                msg <- paste("Since salmon was run in alignment mode, it is recommended you provide",
                             "the total number of input reads via the fallback_num_reads argument to",
                             "prepare_fish_for_sleuth()")
                message(msg)
            }
            else {
                msg <- paste("Fish dir =", fish_dir, ":: fallback # reads =",
                             fallback_num_reads)
                message(msg)
                numProcessed <- fallback_num_reads
            }
        }
        else if (minfo$mapping_type == "mapping" && fallback_num_reads !=
                 -1) {
            msg <- paste("It is generally NOT recommended to provide a fallback # of reads",
                         "when salmon was run in mapping-based mode, as it should already",
                         "accurately record this information. Please make sure you really",
                         "mean to do this!")
            message(msg)
            msg <- paste("Fish dir =", fish_dir, ":: fallback # reads =",
                         fallback_num_reads)
            message(msg)
            numProcessed <- fallback_num_reads
        }
    }
    rhdf5::h5createFile(h5file)
    rhdf5::h5write(quant$est_counts, h5file, "/est_counts")
    rhdf5::h5createGroup(h5file, "aux")
    rhdf5::h5write(numProcessed, h5file, "aux/num_processed")
    rhdf5::h5write(numBoot, h5file, "aux/num_bootstrap")
    rhdf5::h5write(quant$length, h5file, "aux/lengths")
    rhdf5::h5write(quant$eff_length, h5file, "aux/eff_lengths")
    rhdf5::h5write(quant$target_id, h5file, "aux/ids")
    rhdf5::h5write("10", h5file, "aux/index_version")
    rhdf5::h5write("sailfish", h5file, "aux/kallisto_version")
    rhdf5::h5write(timestamp(prefix = "", suffix = ""), h5file,
                   "aux/start_time")
    if (numBoot > 0) {
        rhdf5::h5createGroup(h5file, "bootstrap")
        sapply(0:(numBoot - 1), function(i) {
            bootid <- paste("bs", i, sep = "")
            rhdf5::h5write(unlist(boots[, i + 1]), h5file, paste("bootstrap",
                                                                 bootid, sep = "/"))
        })
    }
    bootCon <- gzcon(file(file.path(auxPath, "fld.gz"), "rb"))
    fld <- readBin(bootCon, "int", n = minfo$frag_dist_length)
    close(bootCon)
    rhdf5::h5write(fld, h5file, "aux/fld")
    bObsCon <- gzcon(file(file.path(auxPath, "observed_bias.gz"),
                          "rb"))
    bObs <- readBin(bObsCon, "int", n = minfo$num_bias_bins)
    close(bObsCon)
    rhdf5::h5write(bObs, h5file, "aux/bias_observed")
    bExpCon <- gzcon(file(file.path(auxPath, "expected_bias.gz"),
                          "rb"))
    bExp <- readBin(bObsCon, "double", n = minfo$num_bias_bins)
    close(bExpCon)
    rhdf5::h5write(bExp, h5file, "aux/bias_normalized")
    rhdf5::H5close()
    print(paste("Successfully converted sailfish / salmon results in",
                fish_dir, "to kallisto HDF5 format"))
})
prepare_fish_for_sleuth <- copy_func(wasabi::prepare_fish_for_sleuth)
environment(prepare_fish_for_sleuth)$fish_to_hdf5 <- fish_to_hdf5

prepare_fish_for_sleuth(fish_dirs=commandArgs(TRUE), force=TRUE)
