#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript scripts/make_demo.R <data.rds>", call. = FALSE)
if (file.exists(args[1L])) stop("Input already exists: ", args[1L], call. = FALSE)

set.seed(824)
n <- 220L
markers_per_group <- 40L
group <- rep(seq_len(12L), each = markers_per_group)
m <- length(group)
snp <- sprintf("demo_LG%02d_%03d", group, rep(seq_len(markers_per_group), 12L))
lead <- match("demo_LG03_021", snp)
freq <- runif(m, 0.15, 0.45)
dosage <- matrix(rbinom(n * m, 2, prob = rep(freq, each = n)), nrow = n, dimnames = list(NULL, snp))
ids <- sprintf("demo_%03d", seq_len(n))
trait <- plogis(-1.2 + 0.65 * dosage[, lead] + rnorm(n, sd = 0.80))
pheno <- data.frame(Taxa = ids, Trait = trait)
geno <- data.frame(Taxa = ids, dosage, check.names = FALSE)
map <- data.frame(marker = snp, chrom = sprintf("LG%02d", group), pos = rep(seq(0, 156, by = 4), 12L))

dir.create(dirname(args[1L]), recursive = TRUE, showWarnings = FALSE)
saveRDS(list(pheno = pheno[sample(n), ], geno = geno, map = map[sample(m), ]), args[1L])
message("Synthetic clone-level input written to ", args[1L])
