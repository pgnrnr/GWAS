args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[1L]
assoc <- read.csv(file.path(out, "associations.csv"))
models <- c("MLMM", "FarmCPU", "BLINK")
stopifnot(setequal(assoc$model, models))
for (model in models) {
  tab <- assoc[assoc$model == model & is.finite(assoc$P.value), ]
  lead <- tab[which.min(tab$P.value), ]
  stopifnot(nrow(tab) == 480L, lead$SNP == "demo_LG03_021", lead$Chr == 3L,
            lead$Pos == 80, lead$P.Bonferroni < 0.05)
  stopifnot(isTRUE(all.equal(tab$FDR.BH, p.adjust(tab$P.value, "BH"))))
  message(model, " recovered ", lead$SNP, " (P = ", format(lead$P.value, digits = 3), ").")
}
figure <- readRDS(file.path(out, "Figure2_GWAS_plot_input.rds"))
stopifnot(all(figure$manhattan_data$model == "FarmCPU"), figure$n_markers == 480L,
          setequal(figure$qq_plot_data$model, models), nrow(figure$qq_plot_data) == 1440L)
stopifnot(all(file.exists(file.path(out, c("manhattan.png", "qq.png", "figure2_GWAS.pdf",
                                         "significant.csv", "session-info.txt")))))
