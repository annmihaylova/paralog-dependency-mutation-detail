# create summary figures and a reproducible exploratory report

plot_rmse <- copy(extended_summary[group == "mutated" & model != "SELECTED_INNER_CV"])
plot_rmse[, model := factor(model, levels = names(extended_model_formulas))]
p <- ggplot2::ggplot(plot_rmse, ggplot2::aes(x = model, y = RMSE_mean)) +
  ggplot2::geom_point(size = 2.1) +
  ggplot2::facet_wrap(~pair, scales = "free_y", ncol = 2) +
  ggplot2::labs(title = "Mutation-positive prediction across 15 feature models",
                subtitle = "Patient-grouped nested CV; lower RMSE is better",
                x = NULL, y = "Mean held-out RMSE") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 65, hjust = 1))
ggplot2::ggsave(file.path(cfg[["figure_dir"]], "Figure_12_extended_models.png"),
                p, width = 15, height = 10, dpi = 220, bg = "white")

# display mutant heterogeneity without dichotomizing the outcome for inference
scatter_data <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  d <- mutant_data[[pairs[["pair"]][i]]][
    , .(ModelID, status, expression_B, effect)]
  d[, pair := pairs[["pair"]][i]]
  d
}))
p <- ggplot2::ggplot(scatter_data,
  ggplot2::aes(x = expression_B, y = effect, color = factor(status))) +
  ggplot2::geom_point(alpha = 0.7, size = 1.9) +
  ggplot2::geom_smooth(method = "lm", se = TRUE, formula = y ~ x) +
  ggplot2::facet_wrap(~pair, scales = "free", ncol = 2) +
  ggplot2::labs(title = "Dependency heterogeneity among damaging-mutant models",
                subtitle = "Lines are unadjusted descriptive fits",
                x = "Target paralog expression (log2[TPM + 1])",
                y = "CRISPR gene-effect", color = "Mutation status") +
  ggplot2::theme_minimal(base_size = 11)
ggplot2::ggsave(file.path(cfg[["figure_dir"]], "Figure_13_mutant_expression.png"),
                p, width = 11, height = 8, dpi = 220, bg = "white")

# show uncertainty for mutant-only predictor comparisons
important_gains <- mutant_gain[
  (from == "U0_status" & to == "U1_target_expression") |
  (from == "U1_target_expression" & to == "U3_target_expression_AF") |
  (from == "U3_target_expression_AF" & to == "U4_both_expressions_AF")]
important_gains[, contrast := paste(from, to, sep = " -> ")]
p <- ggplot2::ggplot(important_gains,
  ggplot2::aes(x = contrast, y = rmse_gain, ymin = ci_low, ymax = ci_high)) +
  ggplot2::geom_pointrange() +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::facet_wrap(~pair, scales = "free_y", ncol = 2) +
  ggplot2::labs(title = "Conditional gains within mutant-only cohorts",
                subtitle = "Patient bootstrap on fixed out-of-fold predictions",
                x = NULL, y = "RMSE reduction") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 28, hjust = 1))
ggplot2::ggsave(file.path(cfg[["figure_dir"]], "Figure_14_mutant_predictor_gains.png"),
                p, width = 13, height = 8, dpi = 220, bg = "white")

if (nrow(expression_threshold_gains)) {
  threshold_plot <- copy(expression_threshold_gains[
    from == "T1_status_expression" & to == "T2_plus_max_AF"])
  p <- ggplot2::ggplot(threshold_plot,
    ggplot2::aes(x = cutoff_label, y = rmse_gain, ymin = ci_low, ymax = ci_high)) +
    ggplot2::geom_pointrange() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::facet_wrap(~pair, scales = "free_y", ncol = 2) +
    ggplot2::labs(title = "Additional value of maximum allele fraction",
                  subtitle = "Expression thresholds are exploratory, not validated biomarkers",
                  x = "Target-expression subset", y = "RMSE reduction") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
  ggplot2::ggsave(file.path(cfg[["figure_dir"]],
                             "Figure_15_expression_thresholds.png"),
                  p, width = 12, height = 7, dpi = 220, bg = "white")
}

# assemble findings from this run rather than copying earlier numbers
format_number <- function(x, digits = 4L) {
  if (!length(x) || is.na(x) || !is.finite(x)) return("NA")
  formatC(x, digits = digits, format = "f")
}
report <- c(
  "# Paralog dependency heterogeneity — exploratory analysis",
  "",
  paste0("Source: ", cfg[["release"]], " | produced ", format(Sys.Date(), "%Y-%m-%d")),
  "",
  "## Research question",
  "",
  "How much do variant-level mutation annotations and target-paralog expression explain variation in CRISPR dependency among mutation-positive cancer cell models, and do mutation features improve prediction after accounting for target expression and cancer lineage?",
  "",
  "## Model definitions",
  "",
  "| Model | Formula |", "|---|---|"
)
for (nm in names(extended_model_formulas)) {
  report <- c(report, paste0("| ", nm, " | `",
    paste(deparse(extended_model_formulas[[nm]]), collapse = " "), "` |"))
}
report <- c(report, "", "## Common-cohort audit", "",
            "| Pair | Lines | Mutants | Patients | Excluded after feature checks |",
            "|---|---:|---:|---:|---:|")
for (i in seq_len(nrow(extended_cohort_audit))) {
  d <- extended_cohort_audit[i]
  report <- c(report, sprintf("| %s | %d | %d | %d | %d |",
    d[["pair"]], d[["n"]], d[["n_mutated"]], d[["n_patients"]],
    d[["excluded_vs_step01"]]))
}
report <- c(report, "", "## Held-out RMSE on mutant lines", "",
            "| Pair | Model | RMSE | MAE |", "|---|---|---:|---:|")
for (i in seq_len(nrow(extended_summary))) {
  d <- extended_summary[i]
  if (d[["group"]] != "mutated") next
  report <- c(report, sprintf("| %s | %s | %s | %s |",
    d[["pair"]], d[["model"]], format_number(d[["RMSE_mean"]]),
    format_number(d[["MAE_mean"]])))
}
report <- c(report, "", "## Mutation-positive variability", "",
            "| Pair | Status | n | Weak > -0.25 | Strong <= -0.5 |",
            "|---|---:|---:|---:|---:|")
for (i in seq_len(nrow(mutant_descriptives))) {
  d <- mutant_descriptives[i]
  report <- c(report, sprintf("| %s | %d | %d | %d | %d |",
    d[["pair"]], d[["status"]], d[["n"]],
    d[["weak_above_minus_025"]], d[["strong_at_most_minus_05"]]))
}
report <- c(report, "", "## Central mutant-only comparisons", "",
            "| Pair | From | To | RMSE gain | Conditional 95% interval | Exploratory q |",
            "|---|---|---|---:|---|---:|")
for (i in seq_len(nrow(mutant_gain))) {
  d <- mutant_gain[i]
  report <- c(report, paste0("| ", d[["pair"]], " | ", d[["from"]], " | ",
    d[["to"]], " | ", format_number(d[["rmse_gain"]]), " | [",
    format_number(d[["ci_low"]]), ", ", format_number(d[["ci_high"]]),
    "] | ", format_number(d[["q_BH_exploratory"]]), " |"))
}
report <- c(report, "", "## Target-expression subgroup analyses", "",
            "Thresholds are exploratory cutoffs in log2(TPM + 1) and were not chosen on an independent dataset",
            "")
for (i in seq_len(nrow(expression_threshold_audit))) {
  d <- expression_threshold_audit[i]
  report <- c(report, paste0("- ", d[["pair"]], " / ", d[["cutoff_label"]],
    ": n=", d[["n"]], ", tested=", d[["tested"]],
    if (d[["tested"]]) "" else paste0(", reason=", d[["skip_reason"]])))
}
report <- c(report, "", "## Sensitivity and specificity checks", "",
            "- `16_cross_target_decoy_interactions.csv` compares different target genes; these are exploratory decoys, not validated null targets",
            "- `16_mutant_sensitivity_regressions.csv` tests target-gene mutation exclusions and CRISPR screen quality",
            "- `17_mutation_expression_associations.csv` tests transcript-level associations rather than causal compensation",
            "- `18_lineage_cluster_intervals.csv` estimates transfer to excluded cancer lineages",
            "", "## Interpretation limits", "",
            "- this is exploratory secondary-data analysis; feature and hypothesis choices were informed by prior examination of the same DepMap release",
            "- no result establishes causality, biallelic inactivation, or a new synthetic-lethal mechanism",
            "- mutation status 2 is not direct evidence of loss of both alleles",
            "- allele fraction and transcript expression are imperfect proxies for protein activity",
            "- the selected model is chosen on inner CV scores, whereas uncertainty intervals condition on existing out-of-fold predictions",
            "- patient-grouped CV and lineage holdout are not independent external datasets",
            "- multiple testing is controlled within displayed exploratory families, not across every analysis performed during hypothesis development",
            "- the expression threshold of 2 and weak-effect thresholds are descriptive post hoc choices",
            "", "## Figures", "",
            "![extended models](figures/Figure_12_extended_models.png)",
            "![mutant expression](figures/Figure_13_mutant_expression.png)",
            "![mutant prediction gains](figures/Figure_14_mutant_predictor_gains.png)",
            if (nrow(expression_threshold_gains))
              "![expression thresholds](figures/Figure_15_expression_thresholds.png)"
            else "expression-threshold figure omitted because no subset qualified"
)
writeLines(report, file.path(cfg[["out_dir"]], "EXTENDED_RESULTS.md"), useBytes = TRUE)
message("Extended report written: ", file.path(cfg[["out_dir"]], "EXTENDED_RESULTS.md"))
