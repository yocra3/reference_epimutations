#'##############################################################################
# Preprocess GSE42861
#'##############################################################################

## Docker command
docker run -it --cpus 5.0 --memory 58Gb -v "$PWD":"$PWD" -w "$PWD" reference_epimutations:1.0  R

## Load libraries  ####
library(GEOquery)
library(minfi)
library(meffil)
library(tidyverse)

options(timeout = 10000)
options(mc.cores = 5)

## Load phenotypes ####
### Download series to data
#### Shorten the file with
###  gzip -cd GSE42861_series_matrix.txt.gz | head -n 1000 > GSE42861_series_matrix.txt
#### Add !series_matrix_table_end  as final line
gse42861.full <- getGEO(filename = "data/GSE42861_series_matrix.txt", getGPL = FALSE)

pheno  <- data.frame(Sample_Name = gse42861.full$geo_accession ,
                     Disease = gse42861.full$`disease state:ch1`,
                     Sex =  gse42861.full$`gender:ch1`, 
                     Smoking =  gse42861.full$`smoking status:ch1`, 
                     Age = as.numeric(gse42861.full$`age:ch1`))
rownames(pheno) <- pheno$Sample_Name

pheno$Sex <- toupper(pheno$Sex)

## Create samplesheet
samplesheet <- meffil.create.samplesheet(path = "data/GSE42861/")
samplesheet$Sex <- pheno[samplesheet$Sample_Name, "Sex"]
samplesheet$Disease <- pheno[samplesheet$Sample_Name, "Disease"]
samplesheet$Age <- pheno[samplesheet$Sample_Name, "Age"]
samplesheet$Smoking <- pheno[samplesheet$Sample_Name, "Smoking"]


## Run QC ####
qc.objects <- meffil.qc(samplesheet, verbose = TRUE,  cell.type.reference="blood gse35069 complete")
qc.parameters <- meffil.qc.parameters(
  beadnum.samples.threshold             = 0.1,
  detectionp.samples.threshold          = 0.05,
  detectionp.cpgs.threshold             = 0.05,
  beadnum.cpgs.threshold                = 0.1,
  sex.outlier.sd                        = 5,
)
qc.summary <- meffil.qc.summary(
  qc.objects,
  parameters = qc.parameters,
)

## Define function to select outliers
filterOutliers <- function(outlier){
  subset(outlier, issue %in% c("Control probe (dye.bias)", 
                               "Methylated vs Unmethylated",
                               "X-Y ratio outlier",
                               "Low bead numbers",
                               "Detection p-value",
                               "Sex mismatch",
                               "Genotype mismatch",
                               "Control probe (bisulfite1)",
                               "Control probe (bisulfite2)")
  )
}

## Remove bad samples based on QC report and rerun QC
outlier <- qc.summary$bad.samples
outlier <- filterOutliers(outlier)
round <- 1

dir.create("results/preprocess/GSE42861/")
  
outPrefix <- "results/preprocess/GSE42861/GSE42861"
save(qc.objects, file = paste0(outPrefix, ".qc.objects.round", round, ".Rdata"))
save(qc.summary, file = paste0(outPrefix, ".qcsummary.round", round, ".Rdata"))
meffil.qc.report(qc.summary, output.file = paste0(outPrefix, ".methylationQC.raw.html"))


outs <- c()
while (nrow(outlier)> 0){
  outs <- rbind(outs, outlier)
  round <- round + 1
  qc.objects <- meffil.remove.samples(qc.objects, outlier$sample.name)
  save(qc.objects, file = paste0(outPrefix,".qc.objects.round", round, ".Rdata"))
  
  qc.summary <- meffil.qc.summary(qc.objects, parameters = qc.parameters)
  save(qc.summary, file = paste0(outPrefix, ".qcsummary.round", round, ".Rdata"))
  outlier <- qc.summary$bad.samples
  outlier <- filterOutliers(outlier)
}
save(qc.objects, file = paste0(outPrefix, ".qc.objects.clean.Rdata"))
save(qc.summary, file = paste0(outPrefix, ".qcsummary.clean.Rdata"))
meffil.qc.report(qc.summary, output.file = paste0(outPrefix, ".methylationQC.clean.html"))


## Report filtered samples and probes
write.table(outs, file = paste0(outPrefix, ".removed.samples.txt"), quote = FALSE, row.names = FALSE,
            sep = "\t")
write.table(qc.summary$bad.cpgs, file = paste0(outPrefix, ".removed.probes.txt"), quote = FALSE, row.names = FALSE,
            sep = "\t")


## Run functional normalization ####
## To be changed in other projects
### Select number PCs (run, see plot and adapt pcs number)
y <- meffil.plot.pc.fit(qc.objects)
ggsave(y$plot, filename = paste0(outPrefix, ".pc.fit.pdf"), height = 6, width = 6)


pcs <- 13
norm.objects <- meffil.normalize.quantiles(qc.objects, number.pcs = pcs)

## Add predicted sex as sample sheet variable
for (i in seq_len(length(norm.objects))){
  norm.objects[[i]]$samplesheet$pred.sex <- norm.objects[[i]]$predicted.sex
}
save(norm.objects, file = paste0(outPrefix, ".norm.obj.pc.Rdata"))

norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove = qc.summary$bad.cpgs$name,
                                      verbose = TRUE)
beta.pcs <- meffil.methylation.pcs(norm.beta, probe.range = 40000)

batch_var <- c("Slide", "sentrix_col",  "sentrix_row", "Sex", "Disease", "Age", "Smoking")

norm.parameters <- meffil.normalization.parameters(
  norm.objects,
  variables = batch_var,
  control.pcs = seq_len(8),
  batch.pcs = seq_len(8),
  batch.threshold = 0.01
)
norm.summary <- meffil.normalization.summary(norm.objects, pcs = beta.pcs, parameters = norm.parameters)
save(norm.summary, file = paste0(outPrefix, ".norm.summary.Rdata"))
meffil.normalization.report(norm.summary, output.file = paste0(outPrefix, ".methylationQC.normalization.html"))

rownames(samplesheet) <- samplesheet$Sample_Name
samplesheet.final <- samplesheet[colnames(norm.beta), ]
## Add predicted sex
samplesheet.final$pred.sex <- vapply(norm.objects, function(x) x$predicted.sex, character(1))

## Add cell counts
cc <- t(sapply(qc.objects, function(obj) obj$cell.counts$counts))
cc <- data.frame(IID=row.names(cc),cc)
samplesheet.final <- cbind(samplesheet.final, cc[rownames(samplesheet.final), ])

## Save genomicratioset
gset <- makeGenomicRatioSetFromMatrix(norm.beta, pData = samplesheet.final,
                                      array = "IlluminaHumanMethylation450k",
                                      annotation = "ilmn12.hg19")
save(gset, file = paste0(outPrefix, ".GenomicRatioSet.Rdata"))

## Create GRSets  ####

ori <- gset

grAnnot <- readRDS("data/HM450.hg19.manifest.rds")

### Probes not measuring methylation
gset <- dropMethylationLoci(gset)
save(gset, file = paste0(outPrefix, ".allCpGs.GenomicRatioSet.Rdata"))

### Remove crosshibridizing and probes with SNPs
gset <- gset[!grAnnot[rownames(gset)]$MASK_general, ]
save(gset, file =  paste0(outPrefix, ".filterAnnotatedProbes.GenomicRatioSet.Rdata"))

### Remove probes in sexual chromosomes
gset <- gset[!seqnames(rowRanges(gset)) %in% c("chrX", "chrY"), ]
save(gset, file = paste0(outPrefix, ".autosomic.filterAnnotatedProbes.GenomicRatioSet.Rdata"))

## Create Initial and final dataset with missings
final <- gset

dp <- meffil.load.detection.pvalues(qc.objects)
gset <- ori
dp.f <- dp[rownames(gset), colnames(gset)]

beta <- assay(gset)
beta[dp.f > 2e-16] <- NA
assay(gset) <- beta
save(gset, file = paste0(outPrefix, ".allCpGs.withNA.GenomicRatioSet.Rdata"))

gset <- gset[rownames(final), colnames(final)]
save(gset, file = paste0(outPrefix, ".autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))

## Se ven bastantes diferencias entre casos y controles y un grupo de muestras separadas en el PC2