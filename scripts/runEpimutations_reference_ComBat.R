#'##############################################################################
# Run Epimutations using GSE105018 as reference
#'##############################################################################

## Docker command
docker run -dit -v "$PWD":"$PWD" -w "$PWD" --memory '250Gb' yocra3/reference_epimutations:1.2  /bin/bash

## Load libraries
library(minfi)
library(sva)
library(epimutacions)

## Load reference dataset
load("results/preprocess/GSE105018/GSE105018.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")
dir.create("results/epimutations/combat_reference/")

applyComBat <- function(combGset){
    
    m <- getM(combGset)
    m[m == -Inf] <- -10
    m[m == Inf] <- 10

    combat_M <- ComBat(dat = m, batch = combGset$Batch, par.prior=TRUE, prior.plots = FALSE, ref.batch = "Reference") 
    beta <- ilogit2(combat_M)
    assay(combGset) <- beta
    combGset
}

correctCombat <- function(newgset, refgset){

    newgset$Batch <- "New"
    refgset$Batch <- "Reference"
    combGset <- combineArrays(refgset, newgset)
    combat_set <- applyComBat(combGset)
    combat_set[, combat_set$Batch == "New"]
}

runEpimutations_combatReference <- function(geoid, ref){

    message(paste("Loading dataset", geoid))
    load(paste0("results/preprocess/", geoid, "/", geoid, ".autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))
  
    if ("Disease" %in% colnames(colData(gset)) & geoid != "GSE51057"){
        gset <- gset[, gset$Disease %in% c("Normal", "control", "Control")]
    }

    message("Combine dataset with reference")
    combatGset <- correctCombat(gset, ref)
    save(combatGset, file = paste0("results/epimutations/combat_reference/", geoid, "_mergedGset.Rdata"))

    combat_epimutations <- epimutations(combatGset, ref, method = "quantile")
    save(combat_epimutations, file = paste0("results/epimutations/combat_reference/", geoid, "_combat_reference.Rdata"))
    numeric(0)
}
gseDatasets <- paste0("GSE", c(42861, 51057, 115278, 121633, 125105, 145361, 163970, 73103, 85210, 87571))

names(gseDatasets) <- gseDatasets
all_epimutations <- lapply(gseDatasets, runEpimutations_combatReference, ref = gset)