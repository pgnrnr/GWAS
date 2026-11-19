# ACE1 forking GWAS. Run from the folder containing data.rds.

library(asreml)
library(ASRgenomics)
library(ASRgwas)
library(GAPIT)
library(ggplot2)
library(qvalue)

set.seed(2026)

# 1. Read the phenotype, genotype, and marker map.
input <- readRDS("data.rds")
pheno_ace1 <- input$pheno
geno_ace1 <- as.matrix(input$geno)
map <- input$map
colnames(map) <- c("marker", "chrom", "pos")
map$marker <- as.character(map$marker)
map$chrom <- as.numeric(sub("^LG", "", as.character(map$chrom)))
map$pos <- as.numeric(as.character(map$pos))
stopifnot(!anyDuplicated(map$marker), !anyDuplicated(rownames(geno_ace1)),
          !anyDuplicated(colnames(geno_ace1)))

pheno_ace1[c("id", "test", "row", "col")] <-
  lapply(pheno_ace1[c("id", "test", "row", "col")], factor)

# 2. Estimate clone means and transform logits to probabilities.
blue_model <- asreml(
  fixed = forkram ~ id + test,
  random = ~ row + col,
  family = asr_binomial(link = "logit"),
  data = pheno_ace1,
  workspace = "8gb"
)
stopifnot(isTRUE(blue_model$converge))

blue_values <- predict(blue_model, classify = "id", pworkspace = "8gb")$pvals
pheno_for_gapit <- data.frame(
  Taxa = as.character(blue_values$id),
  Trait = plogis(blue_values$predicted.value)
)
pheno_for_gapit <- pheno_for_gapit[is.finite(pheno_for_gapit$Trait), ]

# 3. Match clone IDs and marker names before either GWAS.
ids <- intersect(rownames(geno_ace1), pheno_for_gapit$Taxa)
pheno_for_gapit <- pheno_for_gapit[match(ids, pheno_for_gapit$Taxa), ]
map <- map[map$marker %in% colnames(geno_ace1) &
             is.finite(map$chrom) & map$chrom > 0 & is.finite(map$pos), ]
geno_ace1_filtered <- geno_ace1[ids, map$marker, drop = FALSE]
rownames(map) <- map$marker
stopifnot(length(ids) > 5, nrow(map) > 3,
          identical(rownames(geno_ace1_filtered), pheno_for_gapit$Taxa),
          identical(colnames(geno_ace1_filtered), map$marker))

dir.create("results", showWarnings = FALSE)
saveRDS(pheno_for_gapit, "results/clone_probabilities.rds")

# 4. Run the ASRgwas analysis using the original ACE1 settings.
gwas_data_ace1 <- pre.gwas(
  pheno.data = pheno_for_gapit,
  indiv = "Taxa", resp = "Trait",
  geno.data = geno_ace1_filtered, map.data = map,
  Q.method = "geno.data",
  maf = 0.02, marker.callrate = 0.10, pruning.thr = 0.95,
  by.chrom = TRUE, clean.diagonal = TRUE, impute = TRUE
)

gwas_results_ace1 <- gwas.asreml(
  pheno.data = gwas_data_ace1$pheno.data,
  resp = "Trait", gen = "Taxa",
  Kinv = gwas_data_ace1$Kinv, Q = gwas_data_ace1$Q, npc = 3,
  geno.data = gwas_data_ace1$geno.data,
  map.data = gwas_data_ace1$map.data,
  pvalue.thr = 0.05, bonferroni = TRUE, family = "gaussian",
  workspace = "8gb", threads = 2
)

# 5. Calculate significance using raw P values.
assoc <- gwas_results_ace1$gwas.all
valid <- is.finite(assoc$p.value)
stopifnot(any(valid), all(assoc$p.value[valid] >= 0 & assoc$p.value[valid] <= 1))
assoc$p.value[!valid] <- NA_real_
n_markers <- sum(valid)
bonferroni_cutoff <- 0.05 / n_markers
assoc$fdr <- p.adjust(assoc$p.value, method = "BH", n = n_markers)
assoc$significant <- valid & assoc$p.value <= bonferroni_cutoff

qval <- tryCatch(qvalue(assoc$p.value[valid]), error = function(error) {
  warning(conditionMessage(error))
  list(qvalues = rep(NA_real_, n_markers), pi0 = NA_real_)
})
assoc$q.value <- NA_real_
assoc$q.value[valid] <- qval$qvalues
gwas_results_ace1$gwas.all <- assoc

write.csv(assoc, "results/ASRgwas_all.csv", row.names = FALSE)
write.csv(assoc[assoc$significant, ], "results/ASRgwas_significant.csv", row.names = FALSE)
saveRDS(gwas_results_ace1, "results/ASRgwas_results.rds")
print(gwas_results_ace1$heritability)
print(gwas_results_ace1$aov)
print(qval$pi0)

# 6. Plot the ASRgwas results using consensus-map positions.
plot_data <- assoc[valid, ]
map_rows <- match(plot_data$marker, map$marker)
stopifnot(!anyNA(map_rows))
plot_data$chrom <- map$chrom[map_rows]
plot_data$pos <- map$pos[map_rows]
plot_data$logP <- -log10(pmax(plot_data$p.value, .Machine$double.xmin))

manhattan <- ggplot(plot_data, aes(x = pos, y = logP)) +
  geom_point(size = 1, colour = "steelblue4", alpha = 0.7) +
  geom_hline(yintercept = -log10(bonferroni_cutoff), linetype = "dashed") +
  facet_grid(. ~ chrom, scales = "free_x", space = "free_x") +
  labs(x = "Linkage group and position (cM)", y = expression(-log[10](P))) +
  theme_classic()
ggsave("results/ASRgwas_manhattan.pdf", manhattan, width = 12, height = 4)

pdf("results/ASRgwas_QQ.pdf", width = 5, height = 5)
qq.plot(gwas.table = assoc[valid, ])
dev.off()

# 7. Run GAPIT on the original mapped genotypes after QC and imputation.
stopifnot(!anyNA(geno_ace1_filtered))
myGD <- data.frame(Taxa = ids, geno_ace1_filtered, check.names = FALSE)
myGM <- map
colnames(myGM) <- c("SNP", "Chromosome", "Position")

project_dir <- getwd()
dir.create("results/GAPIT", showWarnings = FALSE)
setwd("results/GAPIT")
gwas_gapit <- tryCatch(
  GAPIT(
    Y = pheno_for_gapit, GD = myGD, GM = myGM,
    model = c("MLMM", "FarmCPU", "BLINK"),
    PCA.total = 3, FDRcut = TRUE, cutOff = 0.10,
    file.output = TRUE, Multiple_analysis = FALSE,
    PCA.3d = FALSE, Inter.Plot = FALSE, Geno.View.output = FALSE
  ),
  finally = setwd(project_dir)
)
writeLines(capture.output(sessionInfo()), "results/session-info.txt")
