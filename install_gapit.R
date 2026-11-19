# GAPIT3.5 release commit. The package DESCRIPTION at this commit says 3.4.0.
gapit_commit <- "8d6651c719484c9f6c844144783dca1e4ef85b3e"
options(repos = c(CRAN = "https://cloud.r-project.org"))
Sys.setenv(RGL_USE_NULL = "TRUE")

dependencies <- c("BiocManager", "remotes", "ggplot2", "patchwork")
missing <- dependencies[!vapply(dependencies, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing)
if (!requireNamespace("multtest", quietly = TRUE)) {
  BiocManager::install("multtest", ask = FALSE, update = FALSE)
}
remotes::install_github(paste0("jiabowang/GAPIT@", gapit_commit), dependencies = NA, upgrade = "never")
stopifnot(is.function(GAPIT::GAPIT))
