write_ace1_asrgwas <- function(fit, map) {
  assoc <- fit$gwas.all
  if (!is.data.frame(assoc) || !all(c("marker", "p.value") %in% names(assoc))) {
    stop("ASRgwas gwas.all must contain marker and p.value columns.")
  }
  marker <- check_ids(assoc$marker, "ASRgwas results")
  idx <- match(marker, map$marker)
  if (anyNA(idx)) stop("Some ASRgwas marker IDs are missing from the supplied map.")
  assoc$chrom <- map$chrom[idx]
  assoc$pos <- map$pos[idx]

  p <- as_number(assoc$p.value, "ASRgwas P values")
  tested <- !is.na(p)
  n <- sum(tested)
  if (!n || any(p[tested] < 0 | p[tested] > 1)) stop("Association results require valid P values in [0, 1].")
  assoc$p.value <- p
  assoc$fdr <- stats::p.adjust(p, method = "BH", n = n)
  assoc$significant_bonferroni <- tested & p <= 0.05 / n
  assoc$suggestive <- tested & p <= 1 / n
  assoc$significant_fdr <- tested & assoc$fdr <= 0.05

  q <- tryCatch(qvalue::qvalue(p[tested], lfdr.out = FALSE), error = function(error) {
    warning("Storey q values unavailable: ", conditionMessage(error), call. = FALSE)
    NULL
  })
  assoc$q.value <- NA_real_
  pi0 <- NA_real_
  if (!is.null(q)) {
    assoc$q.value[tested] <- q$qvalues
    pi0 <- q$pi0
  }
  assoc$significant_q <- assoc$q.value <= 0.05
  fit$gwas.all <- assoc

  utils::write.csv(assoc, "gwas_results_all_ace1.csv", row.names = FALSE)
  utils::write.csv(assoc[assoc$significant_bonferroni, ], "gwas_results_significant_ace1.csv", row.names = FALSE)
  saveRDS(fit, "gwas_results_ace1.rds")
  thresholds <- c(bonferroni = 0.05 / n, suggestive = 1 / n, BH = NA_real_)
  if (any(assoc$significant_fdr)) thresholds["BH"] <- max(p[assoc$significant_fdr])
  writeLines(c(
    paste("Tested SNPs", n), paste("Estimated pi0", pi0),
    paste(names(thresholds), "raw P cutoff", format(thresholds, digits = 8)),
    "ASRgwas heritability", capture.output(print(fit$heritability)),
    "ASRgwas ANOVA", capture.output(print(fit$aov))
  ), "gwas_summary_ace1.txt")
  plot_ace1_asrgwas(assoc[tested, , drop = FALSE], thresholds)
}

plot_ace1_asrgwas <- function(assoc, thresholds) {
  library(ggplot2)
  axis <- manhattan_axis(assoc$chrom, assoc$pos)
  assoc$x <- axis$x
  assoc$band <- axis$band
  assoc$logp <- -log10(pmax(assoc$p.value, .Machine$double.xmin))
  lines <- data.frame(method = names(thresholds), raw = unname(thresholds))
  lines <- lines[!is.na(lines$raw), ]
  lines$logp <- -log10(pmax(lines$raw, .Machine$double.xmin))

  manhattan <- ggplot(assoc, aes(x = x, y = logp, colour = band)) +
    geom_point(size = 1.15, alpha = 0.8) +
    geom_hline(data = lines, aes(yintercept = logp, linetype = method), colour = "#8B3C28") +
    scale_colour_manual(values = c("#174A6E", "#2C8577"), guide = "none") +
    scale_linetype_manual(values = c(BH = "solid", bonferroni = "dashed", suggestive = "dotted"),
                          labels = c(BH = "BH FDR (5%)", bonferroni = "Bonferroni (5%)", suggestive = "Suggestive")) +
    scale_x_continuous(breaks = axis$centers, labels = axis$groups) +
    labs(x = "Linkage group", y = expression(-log[10](italic(P))), linetype = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  ggsave("manhattan_ace1.pdf", manhattan, width = 9, height = 4.5)

  qq <- data.frame(expected = -log10(stats::ppoints(nrow(assoc), a = 0.5)),
                   observed = sort(assoc$logp, decreasing = TRUE))
  qq_plot <- ggplot(qq, aes(expected, observed)) +
    geom_abline(colour = "grey50", linetype = "dashed") +
    geom_point(size = 1.2, alpha = 0.75, colour = "#174A6E") +
    labs(x = expression(Expected ~ -log[10](italic(P))),
         y = expression(Observed ~ -log[10](italic(P)))) +
    theme_bw(base_size = 12)
  ggsave("QQ_plot.pdf", qq_plot, width = 5, height = 5)
}
