# Mutation annotation beyond binary status and paralog dependency in DepMap

Analysis of DepMap Public 26Q3 data in R. I look at four known paralog pairs where cancer cell lines with a damaging mutation in gene A tend to depend on its paralog, gene B:

- ARID1A → ARID1B (also analysed separately as a case study)
- SMARCA4 → SMARCA2
- STAG2 → STAG1
- CREBBP → EP300

None of these pairs is new. The question is how much detail about the mutation in gene A helps to predict dependency on gene B, compared with a simple mutated / not mutated call. The extra layers tested are:

- the DepMap mutation category (0/1/2)
- copy number of gene A
- variant allele fraction and counts of loss-of-function variants
- expression of gene B

The project is exploratory: the questions were refined after looking at the same data release.

## Results

All numbers are from the run committed in `results/final_mutation_biomarkers/`. Prediction error is out-of-fold RMSE on mutant lines, from nested cross-validation with folds grouped by patient. Intervals are 95% patient-bootstrap intervals over those predictions.

- A binary mutation call already picks up the known dependencies. For ARID1A → ARID1B, adding it lowers RMSE by 0.051 (0.026 to 0.071).
- Replacing the binary call with the 0/1/2 category gives a small further gain for ARID1A: 0.023 (0.003 to 0.043), BH q = 0.087. The gain comes from category-2 lines.
- Adjusted for lineage, category 2 is associated with a more negative gene effect of gene B in all four pairs: −0.39 for ARID1B, −0.33 for SMARCA2, −0.86 for STAG1, −0.51 for EP300. Category-1 SMARCA4 mutations show no SMARCA2 dependency (+0.04, p = 0.63).
- Copy number of gene A, exact allele fraction and variant counts add nothing detectable in any pair.
- SMARCA2 expression matters a lot for SMARCA4 → SMARCA2. SMARCA4-mutant lines with almost no SMARCA2 expression do not depend on it.
- Within mutant-only cohorts, nothing survives FDR correction (minimum q = 0.23). There are 29–139 mutant lines per pair.
- The category-2 effect changes by at most about 14% after adjusting for MSI, TP53, ploidy, whole-genome doubling, LOH, aneuploidy, sex, growth format and screen quality.

So for these pairs, "probably complete loss of gene A" plus the expression of gene B carries about as much information as the detailed variant annotation.

## Variables

- `effect`: CRISPR gene effect (Chronos) of gene B. More negative means stronger dependency: about 0 is no effect, about −1 is similar to common essential genes.
- `status`: 0/1/2 from `OmicsSomaticMutationsMatrixDamaging.csv`. 0 means no damaging (LikelyLoF) variant. DepMap describes 2 as summed allele fraction > 0.95.
  - In 26Q3, 17 lines are labelled 2 with a sum of 0.905–0.950 (see `tables/01_status_af_audit.csv`).
  - For ARID1A, 31 of 68 category-2 lines have more than one LoF variant.
  - So category 2 does not mean confirmed biallelic loss.
- `binary`: `status > 0`.
- `max_vaf`, `sum_af`, `n_lof`, `n_truncating`, `n_splice`: summaries of LikelyLoF variants in gene A, from `OmicsSomaticMutations.csv`.
- `expression_A`, `expression_B`: log2(TPM + 1).
- `copy_number_A`: relative copy number from `OmicsCNGeneWGS.csv`.
- `lineage_group`: OncotreeLineage. Lineages with fewer than 15 lines are grouped as "Other".

## Data

Raw data are not in the repository. Download DepMap Public 26Q3 from the [DepMap portal](https://depmap.org/portal/data_page/) and put the files in `data/raw/`. Don't mix files from different releases.

Required:

- `Model.csv`
- `CRISPRGeneEffect.csv`
- `OmicsSomaticMutationsMatrixDamaging.csv`
- `OmicsSomaticMutations.csv`
- `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`
- `OmicsCNGeneWGS.csv`

Optional, used only by the covariate checks in steps 09 and 16:

- `OmicsGlobalSignatures.csv`
- `CRISPRConfounders.csv`

If an optional file is missing, its variables are skipped and this is recorded in `tables/09_covariate_availability.csv`.

Only rows with `IsDefaultEntryForModel == "Yes"` are used, and tables are joined by `ModelID`. The per-pair feature tables are committed in `results/final_mutation_biomarkers/features/`.

## How to run

You need R 4.3 or newer and these packages:

```r
install.packages(c("data.table", "glmnet", "sandwich", "lmtest", "ggplot2"))
```

The committed run used R 4.6.1 and glmnet 5.2. Exact versions are in `results/final_mutation_biomarkers/session_info.txt`.

Open `ARID1B-dependency.Rproj` in RStudio, so that the working directory is the project root, and run:

```r
source("scripts/research/run_all.R")
```

Settings are in `cfg` in `scripts/research/00_setup.R`: folders, number of CV repeats and folds, bootstrap and permutation counts, and the seed.

A quick test on synthetic data, with no DepMap files needed:

```r
source("tests/test_invariants.R")
```

## Scripts

All in `scripts/research/`:

| Step | Script | What it does |
|---|---|---|
| 00 | `00_setup.R` | Settings and helper functions |
| 01 | `01_prepare_features.R` | Builds one feature table per pair; audits cohorts and status vs allele fraction |
| 02 | `02_descriptive_tests.R` | Group summaries, Welch ANOVA, Kruskal–Wallis, Holm pairwise tests, lineage-adjusted regression with HC3 standard errors |
| 03 | `03_arid1a_case_study.R` | ARID1A → ARID1B case study: lines with ARID1B mutations excluded, status 1 vs 2 compared within shared lineages |
| 04 | `04_predictive_ablation.R` | Nested ridge CV for models A (binary), B (0/1/2), C (+ copy number), D (+ variant features) |
| 05 | `05_bootstrap_and_diagnostics.R` | Bootstrap of RMSE gains; correlations between predictors |
| 06 | `06_lineage_holdout.R` | Leave-one-lineage-out test |
| 07–08 | `07_figures.R`, `08_results_summary.R` | Figures 1–9 and `RESULTS_SUMMARY.md` |
| 09 | `09_covariate_checks.R` | Covariate balance, change-in-estimate, interactions |
| 10–11 | `10_exploratory_screen.R`, `11_covariate_summary.R` | Optional exhaustive covariate screen (off by default) |
| 12 | `12_extended_model_cv.R` | 15 model variants, patient-grouped nested CV |
| 13 | `13_extended_inference.R` | Paired bootstrap and sign-flip tests for model comparisons |
| 14 | `14_mutant_heterogeneity.R` | Variation in dependency among mutant lines only |
| 15 | `15_target_expression_sensitivity.R` | The same, restricted to lines that express gene B |
| 16 | `16_specificity_and_sensitivity.R` | Decoy targets from the other pairs; excluding lines with gene B mutations; screen-quality adjustment |
| 17 | `17_expression_response.R` | Whether mutation status is associated with expression of gene A and gene B |
| 18 | `18_extended_lineage_holdout.R` | Lineage holdout for the extended models |
| 19 | `19_extended_report.R` | Figures 12–15 and `EXTENDED_RESULTS.md` |

Outputs go to `results/final_mutation_biomarkers/`: `features/`, `tables/` and `figures/`, plus `RESULTS_SUMMARY.md`, `EXTENDED_RESULTS.md` and `session_info.txt`.

## Limitations

- Only one data release, with no independent validation set. Lineage holdout tests transfer to other DepMap lineages, not to new data.
- Cross-validation repeats are correlated. Bootstrap intervals are computed over fixed out-of-fold predictions.
- Mutant groups are small, for example 29 STAG2-mutant lines.
- Mutation category, allele fraction and mRNA level are only proxies for protein loss.
- Results may partly reflect lineage, co-mutations, culture conditions or screen quality.
- Complete-case filtering changes cohort sizes; for example, the full covariate model uses 804 of 1,214 lines.
- p-values are corrected within each family of tests, not across everything that was tried.

## Known issues

- Step 17 adjusts `expression_B` for `expression_A`, which is itself affected by the mutation. This inflates the apparent STAG1 increase in STAG2-mutant lines: +1.31 with the adjustment, +0.35 without it.
- In step 16, the decoy targets are not clean controls. 38% of CREBBP-mutant lines are also ARID1A-mutant, compared with 10% of other lines.
- In step 13, the `comparison` column is labelled incorrectly. The `from` and `to` columns are correct.
- Lineage order, and with it some random seeds, depends on the system locale. Fix: `Sys.setlocale("LC_COLLATE", "C")`.
- `tests/test_invariants.R` overwrites `session_info.txt` in the results folder.

## References

- DepMap, Broad Institute. DepMap Public 26Q3 (cite the figshare record of the release). https://depmap.org/portal
- Tsherniak A et al. (2017). Defining a cancer dependency map. Cell.
- Dempster JM et al. (2021). Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. Genome Biology. doi:10.1186/s13059-021-02540-7
- Helming KC et al. (2014). ARID1B is a specific vulnerability in ARID1A-mutant cancers. Nature Medicine. doi:10.1038/nm.3480
- Kelso et al. (2017). Chromatin accessibility underlies synthetic lethality of SWI/SNF subunits in ARID1A-mutant cancers. eLife. doi:10.7554/eLife.30506
- Schick S et al. (2019). Systematic characterization of BAF mutations provides insights into intracomplex synthetic lethalities in human cancers. Nature Genetics. doi:10.1038/s41588-019-0477-9
- Hoffman GR et al. (2014). Functional epigenetics approach identifies BRM/SMARCA2 as a critical synthetic lethal target in BRG1-deficient cancers. PNAS.
- Ehrenhöfer-Wölfer K et al. (2019). SMARCA2-deficiency confers sensitivity to targeted inhibition of SMARCA4 in esophageal squamous cell carcinoma cell lines. Scientific Reports 9:11661.
- van der Lelij P et al. (2017). Synthetic lethality between the cohesin subunits STAG1 and STAG2 in diverse cancer contexts. eLife. doi:10.7554/eLife.26980
- Benedetti L et al. (2017). Synthetic lethal interaction between the tumour suppressor STAG2 and its paralog STAG1. Oncotarget. doi:10.18632/oncotarget.16838
- Ogiwara H et al. (2016). Targeting p300 addiction in CBP-deficient cancers causes synthetic lethality by apoptotic cell death due to abrogation of MYC expression. Cancer Discovery.
- De Kegel B, Ryan CJ (2019). Paralog buffering contributes to the variable essentiality of genes in cancer cell lines. PLoS Genetics.

## License

DepMap data are generally released under CC BY 4.0. No license has been chosen for the code yet.

Anna Mikhailova, University of Toronto
