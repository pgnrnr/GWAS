#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L || !args[3L] %in% c("ramet", "adjusted")) {
  stop("Usage: Rscript scripts/run_forking.R <data.rds> <results_directory> <ramet|adjusted>", call. = FALSE)
}
source("R/output.R")
source("R/forking.R")
Sys.setenv(RGL_USE_NULL = Sys.getenv("RGL_USE_NULL", "TRUE"))
run_forking(readRDS(args[1L]), args[3L], args[2L])
