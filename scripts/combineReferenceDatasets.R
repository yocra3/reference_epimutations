#'##############################################################################
# Combine reference datasets
#'##############################################################################

## Docker command
docker run -it --cpus 5.0 --memory 58Gb -v "$PWD":"$PWD" -w "$PWD" reference_epimutations:1.0  R

## Load libraries ####
library(minfi)
library(meffil)
library(tidyverse)
library(sva)

## Load files
loadDataset <- function(geoid){
  load(paste0("results/preprocess/", geoid, "/", geoid, ".autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))
  gset$GEOID <- geoid
  gset
}

gseDatasets <- paste0("GSE", c(42861, 51057, 105018, 115278, 121633, 125105, 145361, 163970, 73103, 85210, 87571))
names(gseDatasets) <- gseDatasets
allGRs <- lapply(gseDatasets, loadDataset)

## Remove non-control samples ####
allGRs$GSE42861 <- allGRs$GSE42861[, allGRs$GSE42861$Disease == "Normal"]
allGRs$GSE125105 <- allGRs$GSE125105[, allGRs$GSE125105$Disease == "control"]
allGRs$GSE145361 <- allGRs$GSE145361[, allGRs$GSE145361$Disease == "Control"]
allGRs$GSE163970 <- allGRs$GSE163970[, allGRs$GSE163970$Disease == "Control"]

## Combine datasets ####
### Do in different steps due to memory limitations
allGRs1 <- allGRs[1:5]
allGRs2 <- allGRs[6:11]

combGset1 <- Reduce(combineArrays, allGRs1)
rm(allGRs1, allGRs)
gc()

combGset2 <- Reduce(combineArrays, allGRs2)
rm(allGRs2)
gc()

combGset <- combineArrays(combGset1, combGset2)
rm(combGset1, combGset2)
gc()

dir.create("results/preprocess/combined_reference/")
save(combGset, file = "results/preprocess/combined_reference/combined_reference.autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata")

## Explore global methylation patterns ####
pc_comb <- meffil.methylation.pcs(getBeta(combGset), probe.range = 40000, full.obj = TRUE)

pc_comb.vars <- pc_comb$sdev^2/sum(pc_comb$sdev^2)
pc_comb_df <- data.frame(pc_comb$x[, 1:10], GEO = combGset$GEOID)

png("figures/combGset_PC_raw.png")
ggplot(pc_comb_df, aes(x = PC1, y = PC2, color = GEO)) +
  geom_point() +
  scale_x_continuous(name = paste0("PC1 (", round(pc_comb.vars[1]*100, 1), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(pc_comb.vars[2]*100, 1), "%)")) +
  theme_bw()
dev.off()

### Check differences effect ####
mod <- model.matrix(~ GEOID, colData(combGset))
comb_fit <- lmFit(getBeta(combGset), mod) %>% eBayes()
comb_res <- topTable(comb_fit, coef = 2:11, n = Inf)

### Plot top 9 CpGs with largest differences
top_cpgs <- rownames(comb_res)[1:9]

png("figures/comb_topCpgs_raw.png", width = 3000, height = 3000, res = 300)
getBeta(combGset[top_cpgs, ]) %>%
  t() %>%
  data.frame() %>%
  mutate(GEO = combGset$GEOID,
         sample = colnames(combGset)) %>%
  gather(CpG, Methylation, 1:9) %>%
  ggplot(aes(x = GEO, y = Methylation)) +
    geom_boxplot() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  facet_wrap(~ CpG)
dev.off()

### Apply ComBat ####
m <- getM(combGset)
m[m == -Inf] <- -10
m[m == Inf] <- 10

combat_M <- ComBat(dat = m, batch = combGset$GEOID, par.prior=TRUE, prior.plots = FALSE)

beta <- ilogit2(combat_M)
assay(gset) <- beta
save(gset, file = "results/preprocess/combined_reference/combined_reference.autosomic.filterAnnotatedProbes.withNA.normalizedComBat.GenomicRatioSet.Rdata")



