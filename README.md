# Loblolly pine GWAS workflows

R scripts for genome-wide association studies of stem forking and seedling height in loblolly pine (*Pinus taeda*). The workflows use the R packages and analysis settings from Pinar Guner's forking and dwarfism analyses.

| Workflow | Phenotype | Association methods |
| --- | --- | --- |
| Stem forking with GAPIT3 | ASReml-R binomial GLMM with fixed test and clone effects and random rows and columns within tests. The inverse-logit of each clone mean is the GWAS response. | MLMM, FarmCPU, and BLINK with three PCs. Bonferroni at 0.05 and Benjamini–Hochberg at 0.10. |
| ACE1 stem forking with ASRgwas | ASReml-R binomial GLMM with fixed clone and test effects and random row and column effects. The GWAS uses the inverse-logit clone means. | `pre.gwas()` and `gwas.asreml()` with three PCs, Bonferroni at 0.05, Benjamini–Hochberg at 0.05, and Storey q values. |
| Dwarfism | Measured seedling height at 20 weeks, analyzed as a continuous trait. | `pre.gwas()` and `ASRgwas::gwas.asreml()`. Gaussian mixed model with VanRaden genomic relationships, five PCs, and Bonferroni at 0.05. |

## Run the GAPIT example

Run the commands from the repository root. The example uses synthetic clone-level probabilities with one planted linkage-group 3 association. It exercises the GAPIT stage and figures using R 4.5.0 and the GAPIT3.5 release commit. This upstream commit reports package version `3.4.0` in its `DESCRIPTION` file. [GAPIT3.5 release](https://github.com/jiabowang/GAPIT/releases/tag/GAPIT3.5)

```bash
export RGL_USE_NULL=TRUE
Rscript scripts/install_gapit.R
Rscript scripts/make_demo.R work/demo.rds
Rscript scripts/run_forking.R work/demo.rds results/demo adjusted
Rscript tests/check_demo.R results/demo
```

All three models recovered `demo_LG03_021` as the top marker among 480 SNPs in 220 synthetic clones.

![Synthetic linkage-group Manhattan plots for MLMM, FarmCPU, and BLINK](docs/demo-manhattan.png)

The GitHub Actions workflow runs this synthetic example. The local installation script installs the pinned GAPIT commit and resolves its R dependencies using CRAN and Bioconductor. On Linux, an R toolchain and the system libraries for these packages may be needed.

## Stem forking

Create a single RDS with the existing phenotype, genotype, and consensus-map objects. Store study inputs under `data/`, which Git ignores.

```r
dir.create("data", showWarnings = FALSE)
saveRDS(list(pheno = pheno, geno = geno, map = map), "data/data.rds")
```

| Object | Expected format |
| --- | --- |
| `pheno` | In `ramet` mode, a data frame with `forkram` (0/1), `test`, `id`, `row`, and `col`. In `adjusted` mode, a two-column data frame with clone IDs first and one numeric trait second, such as `Taxa` and `Trait` from `pheno_ace1_blues_prob.rds`. |
| `geno` | A marker matrix with clone IDs as row names or a data frame with `Taxa` in the first column. Observed dosages use the 0–2 scale. GAPIT expects complete genotypes after QC and imputation. The ASRgwas branch passes missing calls to `pre.gwas()`. |
| `map` | A three-column data frame in marker, linkage group, position order. For example `marker`, `chrom`, `pos` or `SNP`, `Chromosome`, `Position`. Linkage groups can be `3` or `LG03`. Positions can be decimal cM coordinates. |

Choose the entry point for the phenotype object you have.

```bash
# Start with individual ramet scores. Requires licensed ASReml-R.
Rscript scripts/run_forking.R data/data.rds results/forking ramet

# Start with saved ASReml-R clone probabilities.
Rscript scripts/run_forking.R data/data.rds results/forking_adjusted adjusted
```

For ramet records, the script saves `pheno_blue_full.rds` with the link-scale estimate, its standard error, and the probability, along with the two-column `pheno_blue_prob.rds` passed to GAPIT. Clone identity is fitted as a fixed effect. GWAS uses the probability-scale clone means as a continuous response. The script matches by clone and marker ID, uses mapped SNPs, and removes markers with zero dosage variance in the matched clones.

The forking study used sample call rate ≥0.90, marker call rate ≥0.95, MAF ≥0.01, and Wright imputation before the GAPIT scan. The wrapper keeps the GAPIT models, three PCs, and BLINK cutoff settings. It disables the built-in genotype view and combined plots because the GAPIT3 release's circular plot can fail with linkage maps. The script produces cross-model plots with ggplot2 and patchwork from the association tables.

| Output | Contents |
| --- | --- |
| `GAPIT.Association.GWAS_Results.*.csv` | Native GAPIT association table for each model. |
| `associations.csv` | Results joined to the input map by SNP ID, with raw and adjusted P values. Corrections are applied per model across finite tests. |
| `significant.csv` | Rows passing Bonferroni at 0.05 or BH at 0.10. |
| `manhattan.png`, `qq.png` | Linkage-map Manhattan plots by model and an overlaid Q–Q plot. The dashed Manhattan line marks Bonferroni. The solid line appears when there are BH discoveries. |
| `figure2_GWAS.pdf`, `Figure2_GWAS_plot_input.rds` | FarmCPU Manhattan plot paired with the MLMM, FarmCPU, and BLINK Q–Q plot, plus the data and thresholds used to make both panels. |
| `session-info.txt` | R and loaded package versions. |

### ACE1 stem forking with ASRgwas

Use the same `data.rds` structure and choose the command that matches the available phenotype object. The ramet model uses `forkram ~ id + test` with `row + col` as random effects. ASReml-R predicts clone means on the logit scale and `plogis()` converts them to the `Taxa` and `Trait` input. Use a clone-level probability in `adjusted` mode.

```bash
Rscript scripts/run_ace1_asrgwas.R data/data.rds results/ace1_asrgwas ramet
Rscript scripts/run_ace1_asrgwas.R data/data.rds results/ace1_asrgwas_adjusted adjusted
```

Install licensed ASReml-R, ASRgenomics, ASRgwas, ggplot2, and the Bioconductor `qvalue` package. The ACE1 code passes `Q.method = "geno.data"`, `maf = 0.02`, `marker.callrate = 0.10`, and `pruning.thr = 0.95` to `pre.gwas()`, with pruning within chromosomes, diagonal cleaning, and imputation enabled. The ASRgenomics QC convention expresses `marker.callrate` as a missing-data cutoff. Use the `pre.gwas()` build from the original analysis to reproduce its filtering. The Gaussian scan uses `Kinv`, `Q`, three PCs, two threads, and an `8gb` workspace by default. [ASRgwas](https://vsni.co.uk/asrgwas/)

The full CSV includes BH and Storey q values and flags for Bonferroni, suggestive, BH, and Storey q-value cutoffs. The significant CSV uses Bonferroni at 0.05. The scan saves the model results as `gwas_results_ace1.rds`, a summary with the estimated null proportion, and Manhattan and Q–Q PDFs. The Manhattan plot uses the numeric consensus-map position within each linkage group. An unavailable Storey estimate is recorded as missing in the table.

## Dwarfism

Use the original `gwas_dataset.rds` structure with `genotype`, `phenotype`, and `map`. Genotype rows are named by `unique_id`. The phenotype data frame contains `unique_id` and `height`. The input map is in marker, chromosome, position order. The published analysis used the Pita50K subset from Nana, WB-10, and WB-20. [Guner et al., *Heredity*](https://doi.org/10.1038/s41437-025-00814-5)

```bash
Rscript scripts/run_dwarfism.R data/gwas_dataset.rds results/dwarfism
```

The script preserves the arguments from the original analysis, including the LD pruning setting of 0.80, imputation, five PCs, Gaussian errors, and Bonferroni correction. It looks for `pre.gwas()` in the original `ASRgenomics` namespace and then in `ASRgwas` for installations that provide it there. Install `ASRgenomics`, the vendor's `ASRgwas`, and a licensed version of ASReml-R. Use the package builds from the original study if you need to reproduce its results. The output directory contains `dwarfism_gwas.rds`, the native files produced by the package, and session information. [ASRgwas download](https://vsni.co.uk/asrgwas/)

## Verification

The synthetic example checks the three GAPIT models, SNP alignment, significance correction, and plot creation. `Rscript tests/check_asrgwas.R work/demo.rds` checks input alignment and reporting with synthetic inputs and an ASRgwas-shaped results table. Verifying the ramet model and both ASRgwas scans requires a licensed ASReml-R installation and the corresponding study inputs.

## Publish on GitHub

Create an empty repository on GitHub, then run these commands from the extracted project directory with your repository URL.

```bash
git init -b main
git add .
git commit -m "Add loblolly pine GWAS workflows"
git remote add origin https://github.com/pgnrnr/plant-gwas.git
git push -u origin main
```

## References

- Wang and Zhang (2021). GAPIT Version 3: Boosting Power and Accuracy for Genomic Association and Prediction. [Genomics, Proteomics & Bioinformatics](https://doi.org/10.1016/j.gpb.2021.08.005)
- Guner, Shalizi, Isik, and Walker (2026). Dwarfism in *Pinus taeda* originates from independent somatic mutations co-localized in a shared genomic region. [Heredity](https://doi.org/10.1038/s41437-025-00814-5)
