# ASRgwas preprocessing moved between package builds. Resolve the installed export.
asr_pre_gwas <- function(args) {
  for (package in c("asreml", "ASRgenomics", "ASRgwas")) {
    if (!requireNamespace(package, quietly = TRUE)) {
      stop("Install ", package, " for the ASRgwas workflow. See the README.", call. = FALSE)
    }
  }
  package <- if ("pre.gwas" %in% getNamespaceExports("ASRgenomics")) "ASRgenomics" else "ASRgwas"
  if (!"pre.gwas" %in% getNamespaceExports(package)) stop("Install a package build that exports pre.gwas().")
  pre_gwas <- getExportedValue(package, "pre.gwas")
  missing <- setdiff(names(args), names(formals(pre_gwas)))
  if (length(missing) && !"..." %in% names(formals(pre_gwas))) {
    stop(package, "::pre.gwas() does not accept: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  do.call(pre_gwas, args)
}

ace1_asrgwas_inputs <- function(pheno, geno, map) {
  aligned <- match_forking_inputs(pheno, geno, map)
  names(aligned$pheno)[2L] <- "Trait"
  if (any(aligned$pheno$Trait < 0 | aligned$pheno$Trait > 1)) {
    stop("Adjusted ACE1 phenotypes must be predicted probabilities between 0 and 1.")
  }
  order <- order(aligned$map_rows)
  geno <- as.matrix(data.frame(aligned$dosages[order], check.names = FALSE))
  rownames(geno) <- aligned$pheno$Taxa
  map <- aligned$map[order, , drop = FALSE]
  rownames(map) <- map$marker
  stopifnot(identical(colnames(geno), map$marker))
  message(sprintf("Matched %d clones and %d mapped SNPs for pre.gwas().", nrow(geno), ncol(geno)))
  list(pheno = aligned$pheno, geno = geno, map = map)
}

run_ace1_asrgwas <- function(data, mode, out) {
  if (!is.list(data) || !all(c("pheno", "geno", "map") %in% names(data))) {
    stop("data.rds must contain a named list with pheno, geno, and map.")
  }
  for (package in c("ggplot2", "qvalue")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Install ", package, " before running this workflow.")
  }
  workspace <- Sys.getenv("ASREML_WORKSPACE", "8gb")
  pheno <- data$pheno
  if (mode == "ramet") {
    adjusted <- clone_probabilities(pheno, spatial = ~ row + col, workspace = workspace)
    pheno <- adjusted[c("Taxa", "Trait")]
  }
  inputs <- ace1_asrgwas_inputs(pheno, data$geno, data$map)

  in_results_dir(out, function() {
    if (mode == "ramet") {
      saveRDS(pheno, "pheno_blue_prob.rds")
      saveRDS(adjusted, "pheno_blue_full.rds")
    }
    set.seed(2026)
    prepared <- asr_pre_gwas(list(
      pheno.data = inputs$pheno, indiv = "Taxa", resp = "Trait",
      geno.data = inputs$geno, map.data = inputs$map,
      Q.method = "geno.data", maf = 0.02, marker.callrate = 0.10,
      pruning.thr = 0.95, by.chrom = TRUE, clean.diagonal = TRUE, impute = TRUE
    ))
    fit <- ASRgwas::gwas.asreml(
      pheno.data = prepared$pheno.data, resp = "Trait", gen = "Taxa",
      Kinv = prepared$Kinv, Q = prepared$Q, npc = 3,
      geno.data = prepared$geno.data, map.data = prepared$map.data,
      pvalue.thr = 0.05, bonferroni = TRUE, family = "gaussian",
      workspace = workspace, threads = 2
    )
    write_ace1_asrgwas(fit, inputs$map)
  })
}
