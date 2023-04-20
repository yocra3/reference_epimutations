#'##############################################################################
# Run Epimutations using GSE105018 as reference
#'##############################################################################

## Docker command
docker run -dit -v "$PWD":"$PWD" -w "$PWD" --memory '250Gb' yocra3/reference_epimutations:1.2  /bin/bash

## Install devel epimutations
devtools::install_github("isglobal-brge/epimutacions@quantile_reference", upgrade = "never")

## Load libraries
library(minfi)
library(epimutacions)

## Load reference dataset
load("results/preprocess/GSE105018/GSE105018.allCpGs.withNA.GenomicRatioSet.Rdata")
dir.create("results/epimutations/quantile_reference/")

ref_quantiles <-  rowQuantiles(getBeta(gset), probs =  c(0.005, 0.995, seq(0.05, 0.95, 0.05)), na.rm = TRUE)

runEpimutations_quantileReference <- function(geoid, ref_quantile){

    message(paste("Loading dataset", geoid))
    load(paste0("results/preprocess/", geoid, "/", geoid, ".autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))
  
    # if ("Disease" %in% colnames(colData(gset)) & geoid != "GSE51057"){
    #     gset <- gset[, gset$Disease %in% c("Normal", "control", "Control")]
    # }
    com_cpgs <- intersect(rownames(gset), rownames(ref_quantiles))
    message("Correcting betas")
    betas_corrected <- epimutacions:::integrateReferenceQuantile(getBeta(gset[com_cpgs, ]), ref_quantiles[com_cpgs, ])
    save(betas_corrected, file = paste0("results/epimutations/quantile_reference/", geoid, "_quantile_reference_betas.Rdata"))

    quantile_epimutations <- epimutations(gset, gset, quantile_reference = ref_quantile, method = "quantile_reference")
    save(quantile_epimutations, file = paste0("results/epimutations/quantile_reference/", geoid, "_quantile_reference.Rdata"))
    numeric(0)
}
gseDatasets <- paste0("GSE", c(42861, 51057, 115278, 121633, 125105, 145361, 163970, 73103, 85210, 87571))

names(gseDatasets) <- gseDatasets
all_epimutations <- lapply(gseDatasets, runEpimutations_quantileReference, ref_quantile = ref_quantiles)


