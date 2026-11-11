# generate a research summary from computed tables

fmt <- function(x, digits = 4L) {
  if (is.na(x)) return("NA")
  formatC(x, digits = digits, format = "f")
}

lines <- c(
  "# Mutation annotation and paralog dependency — Results summary",
  "",
  paste0("Data source: **", cfg[["release"]], "**. Analysis type: secondary-data research."),
  "",
  "## Research question",
  "Does mutation characterization improve prediction of known paralog dependencies beyond binary mutation status?",
  "",
  "## Original A–D comparison",
  "A: binary damaging status + gene A expression + lineage; ",
  "B: A with mutation categories 0/1/2 instead of binary; ",
  "C: B + gene A copy number; ",
  "D: C + VAF and counts/types of likely loss-of-function variants.",
  "",
  "All four models use identical complete-case cohorts and identical outer CV folds per pair.",
  "",
  "## Cohort sizes",
  "| Altered gene -> target gene | All | Status 0 | Status 1 | Status 2 |",
  "|---|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(cohort_audit))) {
  row <- cohort_audit[i]
  lines <- c(lines, sprintf("| %s | %d | %d | %d | %d |",
                            row[["pair"]], row[["n"]], row[["n_0"]],
                            row[["n_1"]], row[["n_2"]]))
}

lines <- c(lines, "", "## Prediction on mutant cell lines",
           "| Pair | RMSE A | RMSE B | RMSE C | RMSE D | Gain A-D | Conditional 95% bootstrap interval |",
           "|---|---:|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(pairs))) {
  pair_name <- pairs[["pair"]][i]
  tab <- cv_summary[pair == pair_name & group == "mutated"]
  vals <- setNames(tab[["RMSE_mean"]], tab[["model"]])
  gain <- ablation_steps[pair == pair_name & group == "mutated", D_minus_A]
  boot <- bootstrap_results[pair == pair_name & group == "mutated" &
                              from == "A_binary" & to == "D_full"]
  ci <- paste0("[", fmt(boot[["bootstrap_low_2.5"]]), ", ",
               fmt(boot[["bootstrap_high_97.5"]]), "]")
  lines <- c(lines, paste0("| ", pair_name, " | ",
    paste(vapply(c("A_binary", "B_category", "C_plus_CN", "D_full"),
                 function(nm) fmt(vals[[nm]]), character(1)), collapse = " | "),
    " | ", fmt(gain), " | ", ci, " |"))
}
lines <- c(lines,
  "",
  "*RMSE gain is positive when the model with detailed annotations predicts better.*",
  "The displayed bootstrap intervals resample cell-line IDs in existing out-of-fold predictions;",
  "they do **not** account for feature-selection, model-selection or complete model-refitting uncertainty.",
  "Repeated CV folds are correlated and are not independent biological replicates.",
  "",
  "## Incremental benefit on mutant lines",
  "| Pair | A→B | B→C | C→D | A→D |",
  "|---|---:|---:|---:|---:|"
)
for (i in seq_len(nrow(pairs))) {
  p <- pairs[["pair"]][i]
  d <- ablation_steps[pair == p & group == "mutated"]
  lines <- c(lines, sprintf("| %s | %s | %s | %s | %s |",
    p, fmt(d[["B_minus_A"]]), fmt(d[["C_minus_B"]]),
    fmt(d[["D_minus_C"]]), fmt(d[["D_minus_A"]])))
}


lines <- c(lines, "", "## Conditional bootstrap by incremental feature group",
           "| Pair | Step | n mutant | RMSE gain | Conditional 95% interval |",
           "|---|---|---:|---:|---:|")
for (i in seq_len(nrow(pairs))) {
  p <- pairs[["pair"]][i]
  for (from_name in c("A_binary", "B_category", "C_plus_CN")) {
    d <- bootstrap_results[pair == p & group == "mutated" & from == from_name]
    if (nrow(d) != 1L) next
    lines <- c(lines, sprintf("| %s | %s → %s | %d | %s | [%s, %s] |",
      p, d[["from"]], d[["to"]], d[["n_lines"]],
      fmt(d[["rmse_gain"]]), fmt(d[["bootstrap_low_2.5"]]),
      fmt(d[["bootstrap_high_97.5"]])))
  }
}

if (exists("holdout_metrics")) {
  lines <- c(lines, "", "## Held-out cancer-lineage validation",
    "| Pair | Subset | Number tested | RMSE A | RMSE D | Gain A-D |",
    "|---|---|---:|---:|---:|---:|")
  for (i in seq_len(nrow(holdout_metrics))) {
    d <- holdout_metrics[i]
    lines <- c(lines, sprintf("| %s | %s | %d | %s | %s | %s |",
      d[["pair"]], d[["group"]], d[["n"]],
      fmt(d[["RMSE_A"]]), fmt(d[["RMSE_D"]]), fmt(d[["gain_A_to_D"]])))
  }
  lines <- c(lines, "Held-out lineage tests evaluate transfer to observed, well-represented lineages;",
             "they do not validate on an independent biological dataset.")
}

if (!is.null(case_data)) {
  case_reg <- robust_coefficients(case_models[["both_expressions"]])
  main_terms <- c("status_factor2", "ARID1B_expression_c",
                  "ARID1A_expression_c:status_factor2")
  lines <- c(lines, "", "## ARID1A–ARID1B detailed case study",
    paste0("Eligible lines (ARID1B damaging mutation excluded): **", nrow(case_data), "**."),
    paste0("Shared-lineage ARID1A-mutant lines: **", nrow(case_shared), "**."),
    "| Coefficient, joint model | Estimate | HC3 SE | HC3 p-value |",
    "|---|---:|---:|---:|")
  # use a loop variable distinct from the term column
  for (term_name in main_terms) {
    d <- case_reg[case_reg[["term"]] == term_name, ]
    if (nrow(d) == 1L) {
      lines <- c(lines, paste0("| ", term_name, " | ", fmt(d[["estimate"]]), " | ",
                               fmt(d[["robust_SE"]]), " | ",
                               format(d[["p_value"]], scientific = TRUE, digits = 3L), " |"))
    }
  }
  lines <- c(lines, "Case-study inference concerns a different selected cohort than the",
             "four-pair prediction study, so coefficient values and sample sizes are not directly interchangeable.")
}

lines <- c(lines, "", "## Figures",
  "![Mutation status comparisons](figures/Figure_01_mutation_groups.png)",
  "![Incremental feature benefit](figures/Figure_02_ablation_steps.png)",
  "![Mutant line prediction error](figures/Figure_03_mutant_RMSE.png)",
  "![Lineage-adjusted estimates](figures/Figure_04_adjusted_effects.png)",
  "![Out-of-fold predictions](figures/Figure_05_OOF_mutant_predictions.png)",
  "",
  "## Interpretation and limitations",
  "- These are associations and prediction benchmarks, not proof of synthetic lethality.",
  "- The 0/1/2 mutation categories are not direct measurements of monoallelic/biallelic loss.",
  "- Most lines in each pair are mutation-negative; report mutant-only metrics alongside overall RMSE.",
  "- The same data informed the hypotheses and candidate-pair choices; discovery claims need independent validation.",
  "- Complete-case selection may affect representativeness and should be disclosed.",
  "- Copy-number measurement units depend on the selected DepMap source file; do not interpret raw CN values as absolute allele counts.",
  "- LOOCV/repeated CV, fold variance, and conditional bootstrap do not replace an independent holdout dataset.",
  "- Mutation associations may reflect tissue lineage, co-mutations, culture conditions or assay characteristics.",
  "",
  "## Main files",
  "- `tables/01_cohort_audit.csv` — included line counts and input provenance.",
  "- `tables/02_global_tests.csv` — Welch and Kruskal-Wallis results.",
  "- `tables/02_lineage_adjusted_HC3.csv` — descriptive adjusted associations.",
  "- `tables/04_oof_predictions.csv` — all repeated held-out predictions.",
  "- `tables/04_cv_summary.csv` and `tables/04_ablation_steps.csv` — central predictive comparisons.",
  "- `tables/05_cluster_bootstrap_oof_gain.csv` — conditional bootstrap intervals.",
  "- `session_info.txt` — R and package versions.",
  "",
  "## Source literature",
  "- Helming et al. (2014), Nature Medicine, DOI: 10.1038/nm.3480.",
  "- Kelso et al. (2017), eLife, DOI: 10.7554/eLife.30506.",
  "- Schick et al. (2019), Nature Genetics, DOI: 10.1038/s41588-019-0477-9.",
  "- DepMap Public 26Q3 documentation and release notes (cite the exact downloaded files).",
  "",
  "Report generated from the present run, not a statement of new biological discovery."
)

writeLines(lines, file.path(cfg[["out_dir"]], "RESULTS_SUMMARY.md"), useBytes = TRUE)
message("Research summary: ", file.path(cfg[["out_dir"]], "RESULTS_SUMMARY.md"))
