#'##############################################################################
# Simulate epimutations with ramr
#'##############################################################################

## Docker command
docker run -it -v "$PWD":"$PWD" -w "$PWD" --memory '250Gb' yocra3/epimutations_simulations  /bin/bash

## Load libraries
library(minfi)
library(ramr)

## Load data
load("results/preprocess/GSE105018/GSE105018.allCpGs.withNA.GenomicRatioSet.Rdata")
gset450 <- gset

load("results/preprocess/GSE131433/GSE131433.allCpGs.withNA.GenomicRatioSet.Rdata")
gsetEpic <- gset[, gset$Age == "adulthood"]


## Convert GRSet to ramr format
ranges450 <- granges(gset450)
mcols(ranges450) <- getBeta(gset450)

rangesEpic <- granges(gsetEpic)
mcols(rangesEpic) <- getBeta(gsetEpic)

## Simulate data ####
set.seed(27)

### 450K
amrs.unique450 <-
  simulateAMR(ranges450, nsamples = 100, regions.per.sample = 10,
              min.cpgs = 3, merge.window = 1000, dbeta = 0.40)

nonunique450 <-
   simulateAMR(ranges450, nsamples = 100, 
              exclude.ranges = amrs.unique450,
              regions.per.sample = 10, samples.per.region = 3, min.cpgs = 3,
              merge.window = 1000, dbeta = 0.4)
noise.450 <-
  simulateAMR(ranges450, nsamples = 100, regions.per.sample = 10,
              exclude.ranges = c(amrs.unique450, nonunique450),
              min.cpgs = 1, max.cpgs = 1, merge.window = 1, dbeta = 0.5)

sim.450 <-
  simulateData(ranges450, nsamples = 100,
               amr.ranges = c(amrs.unique450, nonunique450, noise.450), cores = 5)

sim.450GRS <- makeGenomicRatioSetFromMatrix(data.matrix(mcols(sim.450)))
all.ranges450 <- getUniverse(sim.450, min.cpgs = 3, merge.window = 1000)
save(sim.450GRS, amrs.unique450, nonunique450, all.ranges450, file = "results/simulations/simulated_GRS_450K.Rdata")

#### EPIC
set.seed(27)
amrs.uniqueEpic <-
  simulateAMR(rangesEpic, nsamples = 100, regions.per.sample = 10,
              min.cpgs = 3, merge.window = 1000, dbeta = 0.40)

nonuniqueEpic <-
   simulateAMR(rangesEpic, nsamples = 100, 
              exclude.ranges = amrs.uniqueEpic,
              regions.per.sample = 10, samples.per.region = 3, min.cpgs = 3,
              merge.window = 1000, dbeta = 0.4)
noise.Epic <-
  simulateAMR(rangesEpic, nsamples = 100, regions.per.sample = 10,
              exclude.ranges = c(amrs.uniqueEpic, nonuniqueEpic),
              min.cpgs = 1, max.cpgs = 1, merge.window = 1, dbeta = 0.5)

sim.Epic <-
  simulateData(rangesEpic, nsamples = 100,
               amr.ranges = c(amrs.uniqueEpic, nonuniqueEpic, noise.Epic), cores = 5)

sim.EpicGRS <- makeGenomicRatioSetFromMatrix(data.matrix(mcols(sim.Epic)),
                                      array = "IlluminaHumanMethylationEPIC",
                                      annotation = "ilm10b2.hg19")
all.rangesEpic <- getUniverse(sim.Epic, min.cpgs = 3, merge.window = 1000)
save(sim.EpicGRS, amrs.uniqueEpic, nonuniqueEpic, all.rangesEpic, file = "results/simulations/simulated_GRS_Epic.Rdata")




comb <- combineArrays(gset, sim.newGRS)
comb$Dataset <- ifelse(grepl("samp", colnames(comb)), "Simulation", "Original")
pc_comb <- meffil.methylation.pcs(getBeta(comb), probe.range = 40000, full.obj = TRUE)

pc_comb.vars <- pc_comb$sdev^2/sum(pc_comb$sdev^2)
pc_comb_df <- data.frame(pc_comb$x[, 1:10], Dataset = comb$Dataset)

png("figures/simulation_comb450K_PC.png")
ggplot(pc_comb_df, aes(x = PC1, y = PC2, color = Dataset)) +
  geom_point() +
  scale_x_continuous(name = paste0("PC1 (", round(pc_comb.vars[1]*100, 1), "%)")) +
  scale_y_continuous(name = paste0("PC2 (", round(pc_comb.vars[2]*100, 1), "%)")) +
  theme_bw()
dev.off()
