# Keep a failed run from leaving a directory that looks complete.
in_results_dir <- function(out, run) {
  if (file.exists(out) || dir.exists(out)) stop("Output already exists: ", out)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  parent <- normalizePath(dirname(out), mustWork = TRUE)
  target <- file.path(parent, basename(out))
  temp <- tempfile(".gwas-", tmpdir = parent)
  if (!dir.create(temp)) stop("Cannot create a temporary results directory in ", parent)
  original <- getwd()
  on.exit({
    setwd(original)
    unlink(temp, recursive = TRUE)
  }, add = TRUE)
  setwd(temp)
  run()
  writeLines(capture.output(sessionInfo()), "session-info.txt")
  setwd(original)
  if (file.exists(target) || !file.rename(temp, target)) stop("Could not publish results to ", target)
  message("Results written to ", target)
}

collect_gapit <- function(models, trait, map) {
  results <- lapply(models, function(model) {
    files <- c(
      paste0("GAPIT.Association.GWAS_Results.", model, ".", trait, ".csv"),
      paste0("GAPIT.Association.GWAS_Results.", model, ".Trait(", trait, ").csv")
    )
    found <- files[file.exists(files)]
    if (length(found) != 1L) stop("Expected one GAPIT association table for ", model, " and trait ", trait, ".")
    tab <- utils::read.csv(found, check.names = FALSE, colClasses = c(SNP = "character"))
    columns <- c("SNP", "P.value", "MAF", "nobs", "Effect")
    if (!all(columns %in% names(tab))) stop("Unexpected columns in ", found)
    idx <- match(tab$SNP, map$SNP)
    if (anyNA(idx) || anyDuplicated(tab$SNP)) stop("SNP IDs do not match the map uniquely in ", found)
    p <- tab$P.value
    tested <- is.finite(p) & p >= 0 & p <= 1
    if (!any(tested)) stop("No valid P values in ", found)
    p[!tested] <- NA_real_
    data.frame(
      model = model, trait = trait,
      SNP = tab$SNP, Chr = map$Chromosome[idx], Pos = map$Position[idx],
      P.value = p, MAF = tab$MAF, nobs = tab$nobs, Effect = tab$Effect,
      P.Bonferroni = stats::p.adjust(p, method = "bonferroni", n = sum(tested)),
      FDR.BH = stats::p.adjust(p, method = "BH", n = sum(tested))
    )
  })
  assoc <- do.call(rbind, results)
  assoc <- assoc[order(match(assoc$model, models), assoc$P.value), ]
  utils::write.csv(assoc, "associations.csv", row.names = FALSE)
  significant <- !is.na(assoc$P.value) & (assoc$P.Bonferroni <= 0.05 | assoc$FDR.BH <= 0.10)
  utils::write.csv(assoc[significant, ], "significant.csv", row.names = FALSE)
  plot_gapit(assoc[!is.na(assoc$P.value), ], models)
}

manhattan_axis <- function(chrom, pos) {
  groups <- sort(unique(chrom))
  groups_min <- tapply(pos, factor(chrom, levels = groups), min)
  span <- tapply(pos, factor(chrom, levels = groups), max) - groups_min
  offset <- c(0, head(cumsum(span + max(1, 0.03 * max(span))), -1L))
  group <- match(chrom, groups)
  list(x = pos - groups_min[group] + offset[group], band = factor(group %% 2L),
       centers = offset + span / 2, groups = groups)
}

plot_gapit <- function(assoc, models) {
  library(ggplot2)
  axis <- manhattan_axis(assoc$Chr, assoc$Pos)
  assoc$x <- axis$x
  assoc$band <- axis$band
  assoc$model <- factor(assoc$model, levels = models)
  assoc$logp <- -log10(pmax(assoc$P.value, .Machine$double.xmin))
  cutoffs <- data.frame(
    model = factor(models, levels = models),
    bonf = -log10(0.05 / as.numeric(table(assoc$model))),
    fdr = vapply(models, function(model) {
      p <- assoc$P.value[assoc$model == model & assoc$FDR.BH <= 0.10]
      if (length(p)) -log10(max(p, .Machine$double.xmin)) else NA_real_
    }, numeric(1))
  )
  manhattan_panel <- function(data, thresholds) {
    ggplot(data, aes(x = x, y = logp, colour = band)) +
      geom_point(size = 1.15, alpha = 0.80) +
      geom_hline(data = thresholds, aes(yintercept = bonf), linetype = "dashed", colour = "#8B3C28") +
      geom_hline(data = thresholds[is.finite(thresholds$fdr), ], aes(yintercept = fdr), colour = "#344955") +
      scale_colour_manual(values = c("#174A6E", "#2C8577"), guide = "none") +
      scale_x_continuous(breaks = axis$centers, labels = axis$groups, limits = range(assoc$x)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      labs(x = "Linkage group", y = expression(-log[10](italic(P)))) +
      theme_bw(base_size = 12) +
      theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
  }
  manhattan <- manhattan_panel(assoc, cutoffs) + facet_grid(model ~ .)
  ggsave("manhattan.png", manhattan, width = 8, height = 7.5, dpi = 160)

  qq <- do.call(rbind, lapply(split(assoc, assoc$model), function(group) {
    data.frame(model = group$model[1L], expected = -log10(stats::ppoints(nrow(group), a = 0.5)),
               observed = sort(group$logp, decreasing = TRUE))
  }))
  qq_plot <- ggplot(qq, aes(x = expected, y = observed, colour = model)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 1.2, alpha = 0.75) +
    scale_colour_manual(values = c(MLMM = "#174A6E", FarmCPU = "#A7572F", BLINK = "#2C8577")) +
    labs(x = expression(Expected ~ -log[10](italic(P))),
         y = expression(Observed ~ -log[10](italic(P))), colour = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  ggsave("qq.png", qq_plot, width = 5.5, height = 5, dpi = 160)

  farmcpu <- assoc[assoc$model == "FarmCPU", , drop = FALSE]
  farmcpu_cutoffs <- cutoffs[cutoffs$model == "FarmCPU", , drop = FALSE]
  figure2 <- patchwork::wrap_plots(manhattan_panel(farmcpu, farmcpu_cutoffs), qq_plot, widths = c(1, 0.65))
  ggsave("figure2_GWAS.pdf", figure2, width = 10, height = 5)
  saveRDS(list(manhattan_data = farmcpu, qq_plot_data = qq,
               bonf_threshold = farmcpu_cutoffs$bonf, fdr_threshold = farmcpu_cutoffs$fdr,
               n_markers = nrow(farmcpu)), "Figure2_GWAS_plot_input.rds")
}
