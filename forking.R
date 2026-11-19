# The two-stage loblolly pine forking workflow.

check_ids <- function(ids, label) {
  ids <- as.character(ids)
  if (anyNA(ids) || any(!nzchar(ids)) || any(trimws(ids) != ids) || anyDuplicated(ids)) {
    stop(label, " must have unique, nonempty IDs with no leading or trailing spaces.")
  }
  ids
}

as_number <- function(x, label) {
  if (is.numeric(x)) {
    if (any(is.infinite(x))) stop(label, " contains a nonfinite value.")
    return(as.numeric(x))
  }
  x <- as.character(x)
  missing <- is.na(x) | x %in% c("", "NA", "NaN", ".")
  values <- suppressWarnings(as.numeric(x))
  if (any(!missing & !is.finite(values))) stop(label, " contains a nonnumeric or nonfinite value.")
  values[missing] <- NA_real_
  values
}

clone_probabilities <- function(pheno, spatial = ~ test:row + test:col,
                                workspace = Sys.getenv("ASREML_WORKSPACE", "8gb")) {
  if (!requireNamespace("asreml", quietly = TRUE)) {
    stop("Ramet-level analysis requires a licensed ASReml-R installation. Use adjusted mode with saved clone means.")
  }
  columns <- c("forkram", "test", "id", "row", "col")
  if (!is.data.frame(pheno) || !all(columns %in% names(pheno))) {
    stop("Ramet phenotypes need forkram, test, id, row, and col columns.")
  }
  pheno <- pheno[columns]
  pheno$forkram <- as_number(pheno$forkram, "forkram")
  pheno <- pheno[!is.na(pheno$forkram), , drop = FALSE]
  if (!setequal(unique(pheno$forkram), c(0, 1))) stop("forkram must contain observed 0 and 1 values.")
  for (name in c("test", "id", "row", "col")) {
    ids <- as.character(pheno[[name]])
    if (anyNA(ids) || any(!nzchar(ids)) || any(trimws(ids) != ids)) {
      stop("Observed ramets must have valid ", name, " values.")
    }
    pheno[[name]] <- factor(ids)
  }

  fit <- asreml::asreml(
    fixed = forkram ~ test + id,
    random = spatial,
    family = asreml::asr_binomial(link = "logit"),
    data = pheno,
    na.action = asreml::na.method(x = "fail", y = "fail"),
    workspace = workspace
  )
  if (!isTRUE(fit$converge)) stop("The ASReml-R clone-mean model did not converge.")

  predicted <- stats::predict(fit, classify = "id", aliased = TRUE, pworkspace = workspace)$pvals
  estimable <- is.finite(predicted$predicted.value)
  if ("status" %in% names(predicted)) {
    estimable <- estimable & tolower(as.character(predicted$status)) %in% c("estimable", "estimated")
  }
  if (!any(estimable)) stop("ASReml-R did not return any estimable clone means.")
  data.frame(
    Taxa = as.character(predicted$id),
    eta = ifelse(estimable, predicted$predicted.value, NA_real_),
    se_eta = ifelse(estimable, predicted$standard.error, NA_real_),
    Trait = ifelse(estimable, stats::plogis(predicted$predicted.value), NA_real_)
  )
}

match_forking_inputs <- function(pheno, geno, map) {
  if (!is.data.frame(pheno) || ncol(pheno) != 2L) {
    stop("Adjusted phenotypes must have two columns, Taxa and one clone-level trait.")
  }
  ids <- check_ids(pheno[[1L]], "Phenotypes")
  trait <- names(pheno)[2L]
  if (!nzchar(trait) || grepl("[[:cntrl:]/\\\\:*?\"<>|]", trait) || identical(trait, "Taxa")) {
    stop("Use a nonempty trait name that is safe in output filenames.")
  }
  values <- as_number(pheno[[2L]], trait)

  if (!is.matrix(geno) && !is.data.frame(geno)) stop("geno must be a matrix or data frame.")
  geno <- as.data.frame(geno, stringsAsFactors = FALSE, optional = TRUE)
  if (ncol(geno) && tolower(names(geno)[1L]) %in% c("taxa", "id", "sample_id")) {
    geno_ids <- check_ids(geno[[1L]], "Genotypes")
    geno <- geno[-1L]
  } else {
    if (.row_names_info(geno, type = 1L) < 0L) stop("geno needs Taxa as its first column or named sample rows.")
    geno_ids <- check_ids(rownames(geno), "Genotypes")
  }
  markers <- check_ids(names(geno), "Genotype marker columns")

  if (!is.data.frame(map) || ncol(map) != 3L) {
    stop("map must have three columns in order: marker ID, linkage group, position.")
  }
  map_ids <- check_ids(map[[1L]], "Map markers")
  linkage_group <- sub("^(LG|CHR)[ _-]*", "", toupper(trimws(as.character(map[[2L]]))))
  unplaced <- is.na(linkage_group) | linkage_group %in% c("", "NA", "UNMAPPED", "UNGROUPED", ".")
  if (any(!unplaced & !grepl("^[0-9]+$", linkage_group))) {
    stop("Map linkage groups must be integers or labels such as LG03.")
  }
  groups <- suppressWarnings(as.integer(linkage_group))
  if (any(!unplaced & is.na(groups))) stop("Map linkage group number is out of range.")
  position <- as_number(map[[3L]], "Map position")
  if (any(position < 0, na.rm = TRUE)) stop("Map positions must be nonnegative.")

  rows <- match(geno_ids, ids)
  keep_rows <- !is.na(rows) & !is.na(values[rows])
  if (sum(keep_rows) < 6L || length(unique(values[rows[keep_rows]])) < 3L) {
    stop("At least six matched clones and a varying clone-level response are required.")
  }
  matched_map <- match(markers, map_ids)
  anchored <- !is.na(groups[matched_map]) & groups[matched_map] > 0L & !is.na(position[matched_map])
  dosages <- lapply(geno[keep_rows, anchored, drop = FALSE], as_number, label = "Genotype dosage")
  if (!length(dosages)) stop("No genotype markers have mapped linkage-group positions.")
  if (any(vapply(dosages, function(x) any(x < 0 | x > 2, na.rm = TRUE), logical(1)))) {
    stop("Observed genotype dosages must be between 0 and 2.")
  }
  y <- data.frame(Taxa = geno_ids[keep_rows], value = values[rows[keep_rows]], check.names = FALSE)
  names(y)[2L] <- trait
  map_rows <- matched_map[anchored]
  gm <- data.frame(marker = markers[anchored], chrom = groups[map_rows], pos = position[map_rows])
  list(pheno = y, dosages = dosages, map = gm, map_rows = map_rows,
       n_genotyped = length(geno_ids), n_markers = length(markers))
}

gapit_inputs <- function(pheno, geno, map) {
  aligned <- match_forking_inputs(pheno, geno, map)
  if (any(vapply(aligned$dosages, anyNA, logical(1)))) {
    stop("GAPIT expects complete genotypes for this workflow. Apply genotype QC and imputation first.")
  }
  varying <- vapply(aligned$dosages, function(x) min(x) < max(x), logical(1))
  if (sum(varying) < 4L) stop("At least four variable, mapped SNPs are needed for three PCs and association.")
  y <- aligned$pheno
  gd <- data.frame(Taxa = y$Taxa, aligned$dosages[varying], check.names = FALSE)
  gm <- aligned$map[varying, , drop = FALSE]
  names(gm) <- c("SNP", "Chromosome", "Position")
  stopifnot(identical(names(gd)[-1L], gm$SNP), identical(y$Taxa, gd$Taxa))
  message(sprintf("Matched %d of %d genotyped clones and kept %d of %d SNPs with map positions and variation.",
                  nrow(y), aligned$n_genotyped, nrow(gm), aligned$n_markers))
  list(Y = y, GD = gd, GM = gm)
}

run_forking <- function(data, mode, out) {
  if (!is.list(data) || !all(c("pheno", "geno", "map") %in% names(data))) {
    stop("data.rds must contain a named list with pheno, geno, and map.")
  }
  for (package in c("GAPIT", "ggplot2", "patchwork")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Install ", package, " for the GAPIT workflow. See the README.")
  }
  adjusted <- NULL
  pheno <- data$pheno
  if (mode == "ramet") {
    adjusted <- clone_probabilities(pheno)
    pheno <- adjusted[c("Taxa", "Trait")]
  }
  inputs <- gapit_inputs(pheno, data$geno, data$map)

  in_results_dir(out, function() {
    if (!is.null(adjusted)) {
      saveRDS(pheno, "pheno_blue_prob.rds")
      saveRDS(adjusted, "pheno_blue_full.rds")
    }
    models <- c("MLMM", "FarmCPU", "BLINK")
    set.seed(2026)
    GAPIT::GAPIT(
      Y = inputs$Y, GD = inputs$GD, GM = inputs$GM,
      model = models,
      PCA.total = 3,
      Multiple_analysis = FALSE,
      FDRcut = TRUE,
      cutOff = 0.1,
      file.output = TRUE,
      PCA.3d = FALSE,
      Inter.Plot = FALSE,
      Geno.View.output = FALSE
    )
    collect_gapit(models, names(inputs$Y)[2L], inputs$GM)
  })
}
