#'##############################################################################
# Compare quantile transformation with ComBat
#'##############################################################################

## Docker command
docker run -it -v "$PWD":"$PWD" -w "$PWD" --memory '250Gb' yocra3/reference_epimutations:1.3  /bin/bash

## Load libraries
library(minfi)
library(epimutacions)
library(matrixStats)
library(tidyverse)
library(sva)
library(meffil)
library(ExperimentHub)

eh <- ExperimentHub()
## Load data
load("results/preprocess/GSE197676/GSE197676.allCpGs.withNA.GenomicRatioSet.Rdata")

load( "results/simulations/simulated_GRS_Epic.Rdata")
load( "results/simulations/simulated_GRS_450K.Rdata")


# sim.newGRS$Batch <- "New"
# sim.refGRS$Batch <- "Reference"

# ## Check batch effect in simulated data
# comb_raw <- combineArrays(sim.newGRS, sim.refGRS)

# pc_raw <- meffil.methylation.pcs(getBeta(comb_raw), probe.range = 40000, full.obj = TRUE)

# pc_raw.vars <- pc_raw$sdev^2/sum(pc_raw$sdev^2)
# pc_raw_df <- data.frame(pc_raw$x[, 1:10], Batch = comb_raw$Batch)

# png("figures/simulations.combRaw_PC.png")
# ggplot(pc_raw_df, aes(x = PC1, y = PC2, color = Batch)) +
#   geom_point() +
#   scale_x_continuous(name = paste0("PC1 (", round(pc_raw.vars[1]*100, 1), "%)")) +
#   scale_y_continuous(name = paste0("PC2 (", round(pc_raw.vars[2]*100, 1), "%)")) +
#   theme_bw()
# dev.off()

# ## Compute epimutations with LOO
# epi_loo <- epimutations_one_leave_out(sim.newGRS, method = "quantile")
# save(epi_loo, file = "results/simulations/epimutations_LOO.Rdata")

## Infer epimutations without correcting for batch
com_cpgsEpic <- intersect(rownames(sim.EpicGRS), rownames(gset))
com_cpgs450 <- intersect(rownames(sim.450GRS), rownames(gset))

annotation(sim.450GRS) <- annotation(gset)
epi_raw_Epic <- epimutations(sim.EpicGRS[com_cpgsEpic, ], gset[com_cpgsEpic, ], method = "quantile")
epi_raw_450k <- epimutations(sim.450GRS[com_cpgs450, ], gset[com_cpgs450, ], method = "quantile")

## Infer epimutations with quantile reference
ref_quantiles <- rowQuantiles(getBeta(gset), probs =  c(0.005, 0.995, seq(0.05, 0.95, 0.05)), na.rm = TRUE)

epi_reference_Epic <- epimutations(sim.EpicGRS[com_cpgsEpic, ], sim.EpicGRS[com_cpgsEpic, ], 
                      quantile_reference = ref_quantiles[com_cpgsEpic, ], method = "quantile_reference")

epi_reference_450 <- epimutations(sim.450GRS[com_cpgs450, ], sim.450GRS[com_cpgs450, ], 
                      quantile_reference = ref_quantiles[com_cpgs450, ], method = "quantile_reference")

## Correct batch with ComBat
applyComBat <- function(combGset){
    
    m <- getM(combGset)
    m[m == -Inf] <- -10
    m[m == Inf] <- 10

    combat_M <- ComBat(dat = m, batch = combGset$Batch, par.prior=TRUE, prior.plots = FALSE, ref.batch = "Reference") 
    beta <- ilogit2(combat_M)
    assay(combGset) <- beta
    combGset
}

mergeDataset <- function(newgset, refgset){

    combGset <- combineArrays(refgset, newgset)
    combGset$Batch <- rep(c("Reference", "New"), c(ncol(refgset), ncol(newgset)))
    combat_set <- applyComBat(combGset)
    combat_set
}
combat_Epic <- mergeDataset(sim.EpicGRS, gset)
combat_450 <- mergeDataset(sim.450GRS, gset)
save(combat_Epic, combat_450, file = "results/simulations/combat_datasets.Rdata")
epi_combatEpic <- epimutations(combat_Epic[, combat_Epic$Batch == "New"], combat_Epic[, combat_Epic$Batch == "Reference"], method = "quantile")
epi_combat450 <- epimutations(combat_450[, combat_450$Batch == "New"], combat_450[, combat_450$Batch == "Reference"], method = "quantile")

save(epi_combatEpic, epi_combat450, epi_raw_450k, epi_raw_Epic, epi_reference_450, epi_reference_Epic, file = "results/simulations/epimutations_results.Rdata")



getValues <- function(res, unique_regs, nonunique_regs, all_ranges){
  
  epi_regs <- c(unique_regs, nonunique_regs)
  all_neg <- all_ranges[all_ranges %outside% epi_regs]
  
  if(length(res) == 0){
    return(  c(TPu = 0,
               TPn = 0, 
               FP = 0, 
               TN = length(all_neg),
               FNu = length(unique_regs), 
               FNn = length(nonunique_reg)
            
    ) )
  } else {
  c(TPu = sum(sapply(unique(res$sample), function(x) sum(subset(res, sample == x) %over% subset(unique_regs, sample == x)))), 
    TPn = sum(sapply(unique(res$sample), function(x) sum(subset(res, sample == x) %over% subset(nonunique_regs, sample == x)))), 
    FP = sum(res %outside% epi_regs),
    TN = sum(all_neg %outside%  res),
    FNu = sum(sapply(unique(res$sample), function(x) sum(subset(unique_regs, sample == x) %outside% subset(res, sample == x)))),
    FNn = sum(sapply(unique(res$sample), function(x) sum(subset(nonunique_regs, sample == x) %outside% subset(res, sample == x))))

  ) 
  }
}
# library(ExperimentHub)
# candRegsGR <- epimutacions:::get_candRegsGR()

epimut_res450 <- list(quant_ref = epi_reference_450, raw = epi_raw_450k, ComBat = epi_combat450)
epimut_resEpic <- list(quant_ref = epi_reference_Epic, raw = epi_raw_Epic, ComBat = epi_combatEpic)

epimut_resGR450 <- lapply(epimut_res450, function(x) {
  GR <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  GR[seqnames(GR) != 0]
})
epimut_resGREpic <- lapply(epimut_resEpic, function(x) {
  GR <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  GR[seqnames(GR) != 0]
})
summary_tab450 <- sapply(epimut_resGR450, 
                        getValues, unique_regs = amrs.unique450,
                        nonunique_reg = nonunique450, 
                        all_ranges = all.ranges450)

summary_tabEpic <- sapply(epimut_resGREpic, 
                        getValues, unique_regs = amrs.uniqueEpic,
                        nonunique_reg = nonuniqueEpic, 
                        all_ranges = all.rangesEpic)



summary_tab450 <- summary_tab450 %>%
  t() %>%
  as.tibble(summary_tab450) %>%
  mutate(Methods = colnames(summary_tab450))

summary_tab450 %>%
  mutate(TPRu = TPu/(TPu + FNu),
         TPRn = TPn/(TPn + FNn),
         FDRu = FP / (TPu + FP), 
         FDRn = FP / (TPn + FP)) %>%
  dplyr::select(Methods, starts_with(c("TPR", "FDR")))


summary_tabEpic <- summary_tabEpic %>%
  t() %>%
  as.tibble(summary_tabEpic) %>%
  mutate(Methods = colnames(summary_tabEpic))

summary_tabEpic %>%
  mutate(TPRu = TPu/(TPu + FNu),
         TPRn = TPn/(TPn + FNn),
         FDRu = FP / (TPu + FP), 
         FDRn = FP / (TPn + FP)) %>%
  dplyr::select(Methods, starts_with(c("TPR", "FDR")))

