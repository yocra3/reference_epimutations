#'##############################################################################
# Combine reference datasets
#'##############################################################################

## Docker command
docker run -dit -v "$PWD":"$PWD" -w "$PWD" --memory '250Gb' yocra3/reference_epimutations:1.1  /bin/bash

## Load libraries ####
library(minfi)
library(epimutacions)

## Run epimutations
dir.create("results/epimutations/")

runEpimutations <- function(geoid){
  print(geoid)
  load(paste0("results/preprocess/", geoid, "/", geoid, ".autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))
  gse_epimutations <- mclapply(seq_len(ncol(gset)), function(i) {
    epimutations(gset[, i], gset[, -i], method = "quantile")
  }, mc.cores = 10)
  gse_epimutations <- do.call(rbind, gse_epimutations)
  save(gse_epimutations, file = paste0("results/epimutations/GSE", geoid, "_independent.Rdata"))
  gse_epimutations
}

gseDatasets <- paste0("GSE", c(105018, 42861, 51057, 115278, 121633, 125105, 145361, 163970, 73103, 85210, 87571))

names(gseDatasets) <- gseDatasets
all_epimutations <- lapply(gseDatasets, runEpimutations)


runEpimutations2 <- function(geoid){
  print(geoid)
  load(paste0("results/preprocess/", geoid, "/", geoid, ".autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))
  gse_epimutations <- mclapply(seq_len(ncol(gset)), function(i) {
    epimutations(gset[, i], gset[, -i], method = "quantile")
  }, mc.cores = 10)
  # gse_epimutations <- do.call(rbind, gse_epimutations)
  # save(gse_epimutations, file = paste0("results/epimutations/GSE", geoid, "_independent.Rdata"))
  gse_epimutations
}


cot <- runEpimutations2("GSE105018")
failed <- which(sapply(cot, class) == "NULL")

load("results/preprocess/GSE105018/GSE105018.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
epi_test <- mclapply(1:15, function(i) {
    epimutations(gset[, i], gset[, -i], method = "quantile")
  }, mc.cores = 10)

cot <- epimutations(gset[, 1], gset[, -1], method = "quantile")
