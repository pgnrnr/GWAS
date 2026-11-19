#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript scripts/run_dwarfism.R <gwas_dataset.rds> <results_directory>", call. = FALSE)
}
source("R/output.R")
source("R/asrgwas.R")

gwas_dataset <- readRDS(args[1L])
if (!is.list(gwas_dataset) || !all(c("genotype", "phenotype", "map") %in% names(gwas_dataset))) {
  stop("gwas_dataset.rds must contain genotype, phenotype, and map.", call. = FALSE)
}
geno <- gwas_dataset$genotype
pheno <- gwas_dataset$phenotype
map <- gwas_dataset$map
if (!is.data.frame(pheno) || !all(c("unique_id", "height") %in% names(pheno))) {
  stop("The phenotype table must contain unique_id and height.", call. = FALSE)
}
prep_args <- list(
  pheno.data = pheno, indiv = "unique_id", resp = "height", geno.data = geno,
  Q.method = "K", map.data = map, maf = 0.01, method = "VanRaden",
  marker.callrate = 0.90, ind.callrate = 0.95, heterozygosity = 0.95,
  pruning.thr = 0.80, Fis = 0.98, by.chrom = TRUE, clean.diagonal = TRUE, impute = TRUE
)
in_results_dir(args[2L], function() {
  set.seed(2026)
  gwas_data_dwarf <- asr_pre_gwas(prep_args)
  gwas <- ASRgwas::gwas.asreml(
    pheno.data = gwas_data_dwarf$pheno.data,
    resp = "height", gen = "unique_id",
    Kinv = gwas_data_dwarf$Kinv, Q = gwas_data_dwarf$Q, npc = 5,
    geno.data = gwas_data_dwarf$geno.data, map.data = gwas_data_dwarf$map.data,
    pvalue.thr = 0.05, bonferroni = TRUE, family = "gaussian",
    workspace = "4gb", threads = 2
  )
  saveRDS(gwas, "dwarfism_gwas.rds")
})
