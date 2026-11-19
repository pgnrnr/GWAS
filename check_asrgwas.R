args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
source("R/forking.R")
source("R/asrgwas.R")
source("R/output.R")
source("R/asrgwas_output.R")

data <- readRDS(args[1L])
data$pheno$Trait[data$pheno$Taxa == data$geno$Taxa[1L]] <- NA_real_
data$geno[2L, 2L] <- NA_real_
inputs <- ace1_asrgwas_inputs(data$pheno, data$geno, data$map)
stopifnot(identical(inputs$pheno$Taxa, rownames(inputs$geno)), nrow(inputs$geno) == 219L,
          identical(inputs$map$marker, colnames(inputs$geno)), identical(inputs$map$marker, data$map$marker),
          is.na(inputs$geno[data$geno$Taxa[2L], names(data$geno)[2L]]))

# Verify postprocessing using an ASRgwas-shaped result with 200 tested SNPs.
p <- c(0, 5e-5, 7e-4, 0.004, seq(0.025, 0.997, length.out = 196L), NA_real_)
marker <- rev(inputs$map$marker)[seq_along(p)]
fit <- list(gwas.all = data.frame(marker = marker, p.value = p), heritability = 0.3)
out <- tempfile("ace1-results-")
in_results_dir(out, function() write_ace1_asrgwas(fit, inputs$map))
results <- read.csv(file.path(out, "gwas_results_all_ace1.csv"))
summary <- readLines(file.path(out, "gwas_summary_ace1.txt"))
bh_line <- grep("^BH raw P cutoff ", summary, value = TRUE)
stopifnot(identical(which(results$significant_bonferroni), 1:2),
          identical(which(results$suggestive), 1:4), identical(which(results$significant_fdr), 1:3),
          all(is.finite(results$q.value[seq_len(200L)])), is.na(tail(results$q.value, 1L)),
          all(results$pos == inputs$map$pos[match(marker, inputs$map$marker)]),
          identical(as.numeric(sub("^BH raw P cutoff ", "", bh_line)), 7e-4),
          all(file.exists(file.path(out, c("manhattan_ace1.pdf", "QQ_plot.pdf")))))
message("ASRgwas input alignment, raw P thresholds, BH and Storey q values verified.")
