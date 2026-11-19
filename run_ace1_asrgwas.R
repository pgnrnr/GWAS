#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || !args[3L] %in% c("ramet", "adjusted")) {
  stop("Usage: Rscript scripts/run_ace1_asrgwas.R <data.rds> <results_directory> <ramet|adjusted>", call. = FALSE)
}
source("R/forking.R")
source("R/output.R")
source("R/asrgwas.R")
source("R/asrgwas_output.R")
run_ace1_asrgwas(readRDS(args[1L]), args[3L], args[2L])
