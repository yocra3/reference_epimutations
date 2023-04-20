#'##############################################################################
# Compare quantile transformation with ComBat
#'##############################################################################

## Docker command
docker run -it -v "$PWD":"$PWD" -w "$PWD" --memory '250Gb' yocra3/reference_epimutations:1.3  /bin/bash

## Install devel epimutations
devtools::install_github("isglobal-brge/epimutacions@quantile_reference")

## Load libraries
library(minfi)
library(epimutacions)
library(matrixStats)
library(tidyverse)
library(sva)

gseDatasets <- paste0("GSE", c(42861, 51057, 115278, 121633, 125105, 145361, 163970, 73103, 85210, 87571))
geoid <- gseDatasets[1]
geoid <- "GSE51057"
load(paste0("results/preprocess/", geoid, "/", geoid, ".autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))
gsetTest <- gset

load("results/preprocess/GSE105018/GSE105018.allCpGs.withNA.GenomicRatioSet.Rdata")
refgset <- gset
rm(gset)

## Correct batch with quantile
ref_quantiles <- rowQuantiles(getBeta(refgset), probs =  c(0.005, 0.995, seq(0.05, 0.95, 0.05)), na.rm = TRUE)
com_cpgs <- intersect(rownames(refgset), rownames(gsetTest))
betas_quant <- epimutacions:::integrateReferenceQuantile(getBeta(gset42861[com_cpgs, ]), ref_quantiles[com_cpgs, ])

integrateReferenceQuantile2 <- function(betas, quantiles){
  
  joint <- cbind(quantiles[, -c(1:2)], betas)
  corrected <- apply(joint, 1, function(x) normalizeQuantile2(x[-c(1:19)], x[1:19]))
  corrected <- t(corrected)
  dimnames(corrected) <- dimnames(betas)
  corrected
}


normalizeQuantile2 <- function(vals, quantiles){
  ### Select values from quantiles 5-95%
  vals <- logit2(vals)
  quantiles <- logit2(quantiles)
  center <- vals > quantile(vals, 0.05, na.rm = TRUE) & vals < quantile(vals, 0.95, na.rm = TRUE)
  valsf <- vals[center & !is.na(vals)]
  valsQ <- preprocessCore::normalize.quantiles.use.target(matrix(valsf, ncol = 1), quantiles)
  moddf <- list(Q = valsQ, f = valsf)
  
  # diff <- valsf - valsQ
  mod <- .lm.fit(cbind(matrix(valsf), 1), valsQ)
  
  vals_out <- vals*mod$coefficients[1] + mod$coefficients[2]
  # bpars <- getBetaParams(matrix(valsf, ncol = 1))
  # ps <- pbeta(vals, bpars[1], bpars[2])
  # quants <- pmax(ps, 1 - ps)
  # vals_out <- (1 - quants)/0.5*(vals*mod$coefficients[1] + mod$coefficients[2]) + (quants - 0.5)/0.5*(vals - mean(diff))
  # 
  vals_out <- ilogit2(vals_out)
  vals_out[vals_out < 1e-3] <- 1e-3
  vals_out[vals_out > 1-1e-3] <- 1-1e-3
  vals_out
  
}
betas_quant2 <- integrateReferenceQuantile2(getBeta(gset42861[com_cpgs, ]), ref_quantiles[com_cpgs, ])


integrateReferenceQuantile3 <- function(betas, quantiles){
  
  joint <- cbind(quantiles[, -c(1:2)], betas)
  corrected <- apply(joint, 1, function(x) normalizeQuantile3(x[-c(1:19)], x[1:19]))
  corrected <- t(corrected)
  dimnames(corrected) <- dimnames(betas)
  corrected
}

betas_quant3 <- integrateReferenceQuantile3(getBeta(gsetTest[com_cpgs, ]), ref_quantiles[com_cpgs, ])

integrateReferenceQuantile4 <- function(betas, quantiles){
  
  joint <- cbind(quantiles[, -c(1:2)], betas)
  corrected <- apply(joint, 1, function(x) normalizeQuantile4(x[-c(1:19)], x[1:19]))
  corrected <- t(corrected)
  dimnames(corrected) <- dimnames(betas)
  corrected
}

betas_quant4 <- integrateReferenceQuantile4(getBeta(gset42861[com_cpgs, ]), ref_quantiles[com_cpgs, ])


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

    newgset$Batch <- "New"
    refgset$Batch <- "Reference"
    combGset <- combineArrays(refgset, newgset)
    combat_set <- applyComBat(combGset)
    combat_set[, combat_set$Batch == "New"]
}
combatGset <- mergeDataset(gsetTest, refgset)

diff <-  betas_quant - getBeta(combatGset)
diff2 <-  betas_quant2 - getBeta(combatGset)
diff3 <-  betas_quant3 - getBeta(combatGset)
diff4 <-  betas_quant4 - getBeta(combatGset)

summary(colMeans(abs(diff), na.rm = TRUE))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0008519 0.0009767 0.0010711 0.0011817 0.0012580 0.0057530 

summary(colMeans(abs(diff2), na.rm = TRUE))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0007572 0.0008918 0.0009838 0.0010738 0.0011453 0.0044005

summary(colMeans(abs(diff3), na.rm = TRUE))
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.002249 0.003052 0.003621 0.004038 0.004549 0.019810 

summary(colMeans(abs(diff4), na.rm = TRUE))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0008759 0.0010266 0.0011557 0.0012955 0.0014024 0.0089223 

getOutliers <- function(vec){
  med <- median(vec, na.rm = TRUE)
  top_val <- quantile(vec, 0.995, na.rm = TRUE)
  top <- vec > top_val & vec > med + 0.15
  
  low_val <- quantile(vec, 0.005, na.rm = TRUE)
  low <- vec < low_val & vec < med - 0.15
  top | low
}
outliers <- apply(getBeta(gset42861[com_cpgs, ]), 1, getOutliers)
outliers[is.na(outliers)] <- FALSE
outliers <- t(outliers)
summary(abs(diff[outliers]), na.rm = TRUE)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.00523 0.01153 0.01853 0.02215 0.74718 
summary(abs(diff2[outliers]), na.rm = TRUE)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.003536 0.008091 0.013598 0.016010 0.769546 
summary(abs(diff3[outliers]), na.rm = TRUE)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.0000001 0.0238693 0.0466110 0.0566784 0.0758673 0.8980004 
summary(abs(diff4[outliers]), na.rm = TRUE)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000000 0.009464 0.018377 0.025950 0.032375 0.889611 

mean(abs(diff[outliers]) < 0.01, na.rm = TRUE)
# [1] 0.4462118
mean(abs(diff2[outliers]) < 0.01, na.rm = TRUE)
# [1] 0.5780384
mean(abs(diff3[outliers]) < 0.01, na.rm = TRUE)
# [1] 0.1025128
mean(abs(diff4[outliers]) < 0.01, na.rm = TRUE)
# [1] 0.2658716

png("figures/combat_vs_quantile_GSE42861.png")
data.frame(diff = c(colMeans(abs(diff), na.rm = TRUE), diff[outliers]),
            group = rep(c("Sample Mean", "Outlier"), c(ncol(diff), sum(outliers, na.rm = TRUE)))) %>%
            ggplot(aes(x = diff, color = group)) + 
            geom_density(alpha = 0.1) +
            theme_bw()
dev.off()

sds <- rowSds(getBeta(gset42861), na.rm = TRUE)
medians <- rowMedians(getBeta(gset42861), na.rm = TRUE)

out_mat <- which(outliers, arr.ind = TRUE)

out_df <- data.frame(beta = diff[outliers], 
                    M = diff2[outliers],
                    comb = diff3[outliers],
                    comb2 = diff4[outliers],
                    sd = sds[which(outliers, arr.ind = TRUE)[, 1]],
                    median = medians[which(outliers, arr.ind = TRUE)[, 1]]
) %>%
  gather(Method, Difference, 1:4) %>%
  gather(Estimate, Value, 1:2)

png("figures/combat_vs_quantile_GSE42861_diffs.png")
ggplot(out_df, aes(x = Value, y = Difference)) + 
            geom_point() +
            theme_bw() +
            facet_grid(Estimate ~ Method, scales = "free_y")
dev.off()

out_df <- data.frame(out_mat, beta = abs(diff[outliers]), M = abs(diff2[outliers]), comb = abs(diff3[outliers]), comb2 = abs(diff4[outliers]),
                    sd = sds[which(outliers, arr.ind = TRUE)[, 1]],
                    median = medians[which(outliers, arr.ind = TRUE)[, 1]])
out_df_sel <- subset(out_df, beta > 0.1 | M > 0.1  | comb  > 0.1 | comb  > 0.1)

badcpgs <- rownames(outliers)[subset(out_df, abs(diff) > 0.2)$row]
badcpgs2 <- badcpgs[!badcpgs %in% rownames(multiDistr$proberesults)]


normalizeQuantile3 <- function(vals, quantiles){
  ### Select values from quantiles 5-95%
  center <- vals > quantile(vals, 0.05, na.rm = TRUE) & vals < quantile(vals, 0.95, na.rm = TRUE)
  valsf <- vals[center & !is.na(vals)]
  valsQ <- preprocessCore::normalize.quantiles.use.target(matrix(valsf, ncol = 1), quantiles)
  moddf <- list(Q = valsQ, f = valsf)
  
  diff <- valsf - valsQ
  mod <- .lm.fit(cbind(matrix(valsf), 1), valsQ)
  
  bpars <- epimutacions:::getBetaParams(matrix(valsf, ncol = 1))
  ps <- pbeta(vals, bpars[1], bpars[2])
  quants <- pmax(ps, 1 - ps)
  vals_out <- (1 - quants)/0.5*(vals*mod$coefficients[1] + mod$coefficients[2]) + (quants - 0.5)/0.5*(vals - mean(diff))
   
  vals_out[vals_out < 1e-3] <- 1e-3
  vals_out[vals_out > 1-1e-3] <- 1-1e-3
  vals_out
  
}


normalizeQuantile4 <- function(vals, quantiles){
  ### Select values from quantiles 5-95%
  center <- vals > quantile(vals, 0.05, na.rm = TRUE) & vals < quantile(vals, 0.95, na.rm = TRUE)
  valsf <- vals[center & !is.na(vals)]
  valsQ <- preprocessCore::normalize.quantiles.use.target(matrix(valsf, ncol = 1), quantiles)
  moddf <- list(Q = valsQ, f = valsf)
  
  diff <- valsf - valsQ
  mod <- .lm.fit(cbind(matrix(valsf), 1), valsQ)
  
  # bpars <- epimutacions:::getBetaParams(matrix(valsf, ncol = 1))
  # ps <- pbeta(vals, bpars[1], bpars[2])
  # quants <- pmax(ps, 1 - ps)
  quants <- abs(vals - median(vals, na.rm = TRUE))
  vals_out <- (0.5 - quants)/0.5*(vals*mod$coefficients[1] + mod$coefficients[2]) + quants/0.5 *(vals - median(diff))
   
  vals_out[vals_out < 1e-3] <- 1e-3
  vals_out[vals_out > 1-1e-3] <- 1-1e-3
  vals_out
  
}

cpg_sel <- "cg27506082"
cot <- normalizeQuantile3(getBeta(gset42861[cpg_sel, ]), ref_quantiles[cpg_sel, ])
cot2 <- normalizeQuantile4(getBeta(gset42861[cpg_sel, ]), ref_quantiles[cpg_sel, ])

png("cot.png", width = 1000, height = 1000)
par(mfrow = c(2, 2))
plot(getBeta(gset42861[cpg_sel, ]), getBeta(combatGset[cpg_sel, ]), main = "ComBat")
abline(a = 0, b = 1)
# plot(getBeta(gset42861[cpg_sel, ]), betas_quant["cg19876388", ], main = "Quantile - Beta")
plot(getBeta(gset42861[cpg_sel, ]), cot2, main = "Quantile -  merged")
abline(a = 0, b = 1)
plot(getBeta(combatGset[cpg_sel, ]), cot2, main = "Quantile -  merged")
abline(a = 0, b = 1)
plot(getBeta(combatGset[cpg_sel, ]), betas_quant[cpg_sel, ], main = "Quantile - M")
abline(a = 0, b = 1)

dev.off()

summary(as.vector(abs(getBeta(combatGset[cpg_sel, ]) - betas_quant[cpg_sel, ])))
summary(as.vector(abs(getBeta(combatGset[cpg_sel, ]) - betas_quant2[cpg_sel, ])))
summary(as.vector(abs(getBeta(combatGset[cpg_sel, ]) - cot)))
summary(as.vector(abs(getBeta(combatGset[cpg_sel, ]) - cot2)))

badcpgs2[1:3]
multiDistr <- gaphunter(gset42861, threshold=0.05, verbose=TRUE)




png("figures/combat_vs_quantile_GSE42861_diffMedian.png")
ggplot(out_df, aes(x = median, y = diff)) + 
            geom_point() +
            theme_bw()
dev.off()

epi_reference <- epimutations2(gsetTest, gsetTest, quantile_reference = ref_quantiles, method = "quantile_reference")
epi_combat <- epimutations(combatGset, refgset, method = "quantile")
save(epi_combat, file = "results/epimutations/GSE51057_combRef.Rdata")
load("results/epimutations/GSEGSE51057_independent.Rdata")

epi_reference_loc <- mutate(epi_reference, name = paste(sample, epi_region_id)) %>%
  subset(chromosome != 0 & sample != "GSM1051807") %>%
  pull(name) 
gse_epimutations_loc <- mutate(gse_epimutations, name = paste(sample, epi_region_id)) %>%
  subset(chromosome != 0 & sample != "GSM1051807") %>%
  pull(name) 
epi_combat_loc <- mutate(epi_combat, name = paste(sample, epi_region_id)) %>%
  subset(chromosome != 0 & sample != "GSM1051807") %>%
  pull(name) 

mean(epi_reference_loc %in% epi_combat_loc)
sum(epi_reference_loc %in% epi_combat_loc)
sum(!epi_reference_loc %in% epi_combat_loc)
sum(!epi_combat_loc %in% epi_reference_loc)

mean(epi_reference_loc %in% gse_epimutations_loc)
sum(epi_reference_loc %in% gse_epimutations_loc)
sum(!epi_reference_loc %in% gse_epimutations_loc)
sum(!gse_epimutations_loc %in% epi_reference_loc)


mean(gse_epimutations_loc %in% epi_combat_loc)
sum(gse_epimutations_loc %in% epi_combat_loc)
sum(!gse_epimutations_loc %in% epi_combat_loc)
sum(!epi_combat_loc %in% gse_epimutations_loc)

BiocManager::install("VennDetail")
library(VennDetail)

ven <- venndetail(list(ComBat = epi_combat_loc, Reference = epi_reference_loc,
                    LOO = gse_epimutations_loc))

png("pac.png")
plot(ven)
dev.off()


bad_regs <- setdiff(setdiff(epi_reference_loc, epi_combat_loc), gse_epimutations_loc)

epi_reference_bad <- subset(epi_reference, paste(sample, epi_region_id) %in% bad_regs)

reg <- epi_reference_bad[3, ]
samp <- reg$sample
cpgs <- strsplit(reg$cpg_ids, ",")[[1]]

png("reg3.png", width = 2000, height = 1300)
cpg_plot <- cbind(getBeta(gsetTest[cpgs, ]), betas_quant3[cpgs, , drop = FALSE], getBeta(combatGset[cpgs, ]),  getBeta(refgset[cpgs, ])) %>%
  t() %>%
  data.frame() %>%
  mutate(Dataset = rep(c("Raw", "Quantile", "ComBat", "Reference"), c(ncol(gsetTest), ncol(betas_quant3), ncol(combatGset), ncol(refgset))),
        Dataset = factor(Dataset, levels = c("Raw", "Quantile", "ComBat", "Reference"))) %>%
  gather(CpG, Methylation, seq_len(length(cpgs))) %>%
  ggplot(aes(x = Dataset, y = Methylation)) +
    geom_boxplot() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  facet_wrap(~ CpG)

plot_grid(
  plot_epimutations(reg, combineArrays(makeGenomicRatioSetFromMatrix(betas_quant3[cpgs, samp, drop = FALSE]), refgset[cpgs, ])) + ggtitle("Reference"),
  plot_epimutations(reg, combineArrays(combatGset[cpgs, samp], refgset[cpgs, ])) + ggtitle("Combat"),
  plot_epimutations(reg, gsetTest) + ggtitle("original"),
  cpg_plot,
  ncol = 2)
dev.off()

data.frame(diffMedian = rowMedians(getBeta(refgset[cpgs, ])) - rowMedians(getBeta(gsetTest[cpgs, ])) ,
           diffQuantile = as.vector(betas_quant3[cpgs, samp, drop = FALSE] - getBeta(gsetTest[cpgs, samp])),
          diffComBat = as.vector(getBeta(combatGset[cpgs, samp]) - getBeta(gsetTest[cpgs, samp])))

epi_mini <- epimutations(gsetTest[cpgs, samp], gsetTest[cpgs, colnames(gsetTest) != samp], method = "quantile")


epimutations2 <- function(case_samples, control_panel,
                        method = "manova", 
                        chr = NULL, start = NULL, end = NULL, 
                        epi_params = epi_parameters(), 
                        maxGap = 1000, bump_cutoff =  0.1, 
                        min_cpg = 3, verbose = TRUE, quantile_reference = NULL)
{
    ### Feature annotation
    fd <- as.data.frame(GenomicRanges::granges(case_samples))
    rownames(fd) <- rownames(case_samples)
    ### Betas
    betas_case <- minfi::getBeta(case_samples)
    betas_case <- betas_case[rownames(fd), , drop = FALSE]
    betas_control <- minfi::getBeta(control_panel)
    betas_control <- betas_control[rownames(fd), ]
   
    if (!is.null(chr)) {
        if (!is.null(start) & !is.null(end)) {
            fd <- fd[fd$seqnames %in% chr & fd$start>=start & fd$end <= end, ]
        } else{
            fd <- fd[fd$seqnames %in% chr, ]
        }
        if (nrow(fd) == 0) {
            stop("No CpG was found in the specified region")
        }
        
        betas_case <- betas_case[rownames(fd), , drop = FALSE]
        betas_control <- betas_control[rownames(fd), ]
    }
    
    ### Identify case and control samples
    cas_sam <- colnames(betas_case)
    ctr_sam <- colnames(betas_control)
    # 2. Epimutations definition (using different methods)
    ##Methods that need bumphunter
    ##("manova", "mlm", "mahdist" and "iForest")
    if (method == "quantile_reference") {
      # Compute reference statistics
      if (verbose)
        message("Using quantiles reference")
      
      com_cpgs <- intersect(rownames(betas_case), rownames(quantile_reference))
      quant_ref_com <- quantile_reference[com_cpgs, ]
      
      bctr_pmin <- quant_ref_com[, "0.5%"]
      bctr_pmax <- quant_ref_com[, "99.5%"]
     
      if (verbose)
        message("Integrating cases into reference")
      
      betas <- cbind(betas_case, betas_control)
      betas_corrected <- integrateReferenceQuantile3(betas[com_cpgs, ], quant_ref_com)
      
      if (verbose)
        message("Computing epimutations")
      
      # Run region detection
      rst <- do.call(rbind, lapply(cas_sam, function(case) {
        x <- epi_quantile( betas_corrected[, case, drop = FALSE],
                           fd, bctr_pmin, bctr_pmax,
                           quant_ref_com[, "50%"], epi_params$quantile$window_sz,
                           min_cpg, epi_params$quantile$offset_abs )
        if (is.null(x) || nrow(x) == 0) {
          x <- data.frame( chromosome = 0, start = 0, end = 0, sz = NA,
                           cpg_n = NA, cpg_ids = NA, outlier_score = NA,
                           outlier_direction = NA, pvalue = NA,
                           adj_pvalue = NA, delta_beta = NA, 
                           sample = case )
        } else {
          x$sample <- case
        }
        x
      }))
    } 
    # 3. Prepare the output and addition of CREs
    ## Prepare the output
    rst$epi_id <- vapply(seq_len(nrow(rst)), function(ii)
                            paste0("epi_", method, "_", ii), character(1))
    rownames(rst) <- seq_len(nrow(rst))
    rst <- rst[, c(13, 12, seq_len(11))]
    
    ## Add CREs and epi_region_id
    rst$CRE_type <- rst$CRE <- rst$epi_region_id <- NA
    rst_c <- rst
    rst_c <- tryCatch({
        rstGR <- GenomicRanges::makeGRangesFromDataFrame(rst)
        ## Ensure chromosomes have the same format
        seqlevelsStyle(rstGR) <- "UCSC"
        #Get candidate regions
        candRegsGR <- get_candRegsGR()
        over <- GenomicRanges::findOverlaps(rstGR, candRegsGR)
        #variables (avoid long code)
        ids <- names(candRegsGR[S4Vectors::to(over)])
        cre <- candRegsGR[S4Vectors::to(over)]$CRE
        cre_type <- candRegsGR[S4Vectors::to(over)]$CRE_type
        
        rst$epi_region_id[S4Vectors::from(over)] <- ids
        rst$CRE[S4Vectors::from(over)] <- cre
        rst$CRE_type[S4Vectors::from(over)] <- cre_type
        rm(c("ids", "cre", "cre_type"))
        rst
    }, error = function(e) {
        rst
    })
    
    ## Convert rst into a tibble class
    if (requireNamespace("tibble", quietly = TRUE)) {
        rst <- tibble::as_tibble(rst_c)
    } else {
        stop("'tibble' package not avaibale")
    }
    return(rst)
}


pdf("cot.pdf")
plot_epimutations(epis_quantile[5, ], combineArrays(combatGset[, epis_quantile[5, ]$sample], refgset)) + ggtitle("Combat")
plot_epimutations(epis_quantile[5, ], gset42861) + ggtitle("original")
dev.off()

sel_cpgs <- strsplit(quantile_epimutations[1, ]$cpg_ids, ",")[[1]]

ref_quantiles[sel_cpgs,]
rowQuantiles(getBeta(gset42861[sel_cpgs,]), probs =  c(0.005, 0.995, seq(0.05, 0.95, 0.05)), na.rm = TRUE)

betas_reg <- epimutacions:::integrateReferenceQuantile(getBeta(gset42861[sel_cpgs, ]), ref_quantiles[sel_cpgs, ])
betas_reg[, "GSM1051533"]
getBeta(gset42861[sel_cpgs, "GSM1051533"])

ep <- epimutations(gset42861[sel_cpgs, "GSM1051533"], gset42861[sel_cpgs, -1], method = "quantile")

combGset$Study <- ifelse(colnames(combGset) %in% colnames(gset42861), "Case-Combat", "Reference-Combat")

gset$Study <- "Reference"
gset42861$Study <- "Case"

mergeGset <- combineArrays(gset[sel_cpgs, ], gset42861[sel_cpgs, ])

combGset2 <- mergeDataset(gset42861, gset)

largeGset <- combineArrays(mergeGset, combGset2[sel_cpgs, ])
largeGset$Study <- c(mergeGset$Study, paste(combGset2$Study, "-ComBat"))


integrateReferenceQuantile2 <- function(betas, quantiles){
  
  joint <- cbind(quantiles[, -c(1:2)], betas)
  corrected <- apply(joint, 1, function(x) normalizeQuantile2(x[-c(1:19)], x[1:19]))
  corrected <- t(corrected)
  dimnames(corrected) <- dimnames(betas)
  corrected
}


normalizeQuantile2 <- function(vals, quantiles){
  ### Select values from quantiles 5-95%
  center <- vals > quantile(vals, 0.05, na.rm = TRUE) & vals < quantile(vals, 0.95, na.rm = TRUE)
  valsf <- vals[center & !is.na(vals)]
  valsQ <- preprocessCore::normalize.quantiles.use.target(matrix(valsf, ncol = 1), quantiles)
  moddf <- list(Q = valsQ, f = valsf)
  
  diff <- valsf - valsQ
  mod <- .lm.fit(cbind(matrix(valsf), 1), valsQ)
  bpars <- epimutacions:::getBetaParams(matrix(valsf, ncol = 1))
  ps <- pbeta(vals, bpars[1], bpars[2])
  quants <- pmax(ps, 1 - ps)
  vals_out <- vals*mod$coefficients[1] + mod$coefficients[2]
#   vals_out <- (1 - quants)/0.5*(vals*mod$coefficients[1] + mod$coefficients[2]) + (quants - 0.5)/0.5*(vals - mean(diff))
  
  vals_out[vals_out < 1e-3] <- 1e-3
  vals_out[vals_out > 1-1e-3] <- 1-1e-3
  vals_out
  
}

betas_reg2 <- integrateReferenceQuantile2(getBeta(gset42861[sel_cpgs, ]), ref_quantiles[sel_cpgs, ])


pdf("pac.pdf")
  cbind(getBeta(largeGset), betas_reg, betas_reg2) %>%
  t() %>%
  data.frame() %>%
  mutate(GEO = c(largeGset$Study, rep(c("mod_ref1", "mod_ref2"), c(ncol(betas_reg), ncol(betas_reg2))))) %>%
  gather(CpG, Methylation, 1:3) %>%
  ggplot(aes(x = GEO, y = Methylation)) +
    geom_boxplot() +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5)) +
  facet_wrap(~ CpG)
dev.off()


## Try new methods
## Install devel epimutations
devtools::install_github("isglobal-brge/epimutacions@quantile_reference")

## Load libraries
library(minfi)
library(epimutacions)
library(matrixStats)
library(tidyverse)
library(sva)

## Load reference dataset
load("results/preprocess/GSE105018/GSE105018.allCpGs.withNA.GenomicRatioSet.Rdata")
ref_quantiles <-  rowQuantiles(getBeta(gset), probs =  c(0.005, 0.995, seq(0.05, 0.95, 0.05)), na.rm = TRUE)

geoid <- "GSE42861"
load(paste0("results/preprocess/", geoid, "/", geoid, ".autosomic.filterAnnotatedProbes.withNA.GenomicRatioSet.Rdata"))

if ("Disease" %in% colnames(colData(gset)) & geoid != "GSE51057"){
    gset <- gset[, gset$Disease %in% c("Normal", "control", "Control")]
}

quantile_epimutations_new <- epimutations(gset, gset, quantile_reference = ref_quantiles, method = "quantile_reference")

gset42861 <- gset
load("results/preprocess/GSE105018/GSE105018.allCpGs.withNA.GenomicRatioSet.Rdata")

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

    newgset$Batch <- "New"
    refgset$Batch <- "Reference"
    combGset <- combineArrays(refgset, newgset)
    combat_set <- applyComBat(combGset)
    combat_set
}
combGset <- mergeDataset(gset42861, gset)

combat_epimutations_new <- epimutations(combGset[, combGset$Batch == "New"], combGset[, combGset$Batch == "Reference"], method = "quantile")
