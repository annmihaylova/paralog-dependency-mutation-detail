# draw figures for the original analysis

save_plot <- function(p, name, w = 10, h = 5.8) {
  ggplot2::ggsave(file.path(cfg[["figure_dir"]], name), p,
                  width = w, height = h, dpi = 220, bg = "white")
}

# figure 1 shows dependency by mutation status
plot_data <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  d <- feature_data[[pairs[["stem"]][i]]][, .(status, effect)]
  d[, pair := pairs[["pair"]][i]]
  d
}))
plot_data[, status := factor(status, levels = 0:2)]
p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = status, y = effect)) +
  ggplot2::geom_boxplot(width = 0.60, outlier.alpha = 0.26) +
  ggplot2::facet_wrap(~pair, scales = "free_y", ncol = 2) +
  ggplot2::labs(title = "Paralog dependency by damaging-mutation category",
                subtitle = "Exploratory comparison; status 2 does not prove biallelic loss",
                x = "Mutation status of gene A", y = "CRISPR gene-effect score") +
  ggplot2::theme_minimal(base_size = 12)
save_plot(p1, "Figure_01_mutation_groups.png", 11, 7)

# figure 2 shows the incremental contribution to RMSE
steps_plot <- copy(ablation_steps[group == "all" | group == "mutated"])
step_long <- melt(steps_plot,
                  id.vars = c("pair", "group"),
                  measure.vars = c("B_minus_A", "C_minus_B", "D_minus_C"),
                  variable.name = "step", value.name = "rmse_gain")
step_long[, step := factor(step, levels = c("B_minus_A", "C_minus_B", "D_minus_C"),
                          labels = c("Mutation category", "+ Copy number", "+ Detailed variants"))]
step_long[, group := factor(group, levels = c("all", "mutated"),
                            labels = c("All cell lines", "Mutant lines only"))]
p2 <- ggplot2::ggplot(step_long,
  ggplot2::aes(x = step, y = rmse_gain, fill = group)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.72), width = 0.67) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  ggplot2::facet_wrap(~pair, scales = "free_y", ncol = 2) +
  ggplot2::labs(title = "Incremental predictive value of mutation annotations",
                subtitle = "Positive values indicate reduced RMSE",
                x = NULL, y = "RMSE reduction", fill = "Evaluation cohort") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))
save_plot(p2, "Figure_02_ablation_steps.png", 12, 7)

# figure 3 compares mutant-only prediction errors
model_plot <- copy(cv_summary[group == "mutated"])
model_plot[, model := factor(model, levels = names(model_formulas),
                            labels = c("A: binary", "B: 0/1/2", "C: + CN", "D: + variants"))]
p3 <- ggplot2::ggplot(model_plot,
    ggplot2::aes(x = model, y = RMSE_mean, group = pair)) +
  ggplot2::geom_line(linewidth = 0.65) +
  ggplot2::geom_point(size = 2.3) +
  ggplot2::facet_wrap(~pair, scales = "free_y", ncol = 2) +
  ggplot2::labs(title = "Cross-validated prediction error in mutant cell lines",
                subtitle = "Repeated nested five-fold CV; same cell lines in all models",
                x = NULL, y = "Mean RMSE over five repeats") +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 18, hjust = 1))
save_plot(p3, "Figure_03_mutant_RMSE.png", 11, 7)

# figure 4 displays adjusted status associations
forest <- copy(adjusted_coefficients[
  term == "status_factor1" | term == "status_factor2"])
forest[, category := factor(term, levels = c("status_factor1", "status_factor2"),
                            labels = c("Status 1 vs 0", "Status 2 vs 0"))]
p4 <- ggplot2::ggplot(forest,
    ggplot2::aes(x = pair, y = estimate, ymin = CI_low, ymax = CI_high,
                 shape = category)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
  ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = 0.50)) +
  ggplot2::coord_flip() +
  ggplot2::labs(title = "Lineage- and expression-adjusted associations",
                subtitle = "HC3 robust 95% intervals; relative to mutation status 0",
                x = NULL, y = "Difference in gene-effect score", shape = "Comparison") +
  ggplot2::theme_minimal(base_size = 12)
save_plot(p4, "Figure_04_adjusted_effects.png", 11, 6)

# figure 5 compares observed and predicted gene effects
scatter <- copy(all_predictions[repeat_id == 1 & (status == 1 | status == 2)])
p5 <- ggplot2::ggplot(scatter,
    ggplot2::aes(x = observed, y = pred_D_full, color = factor(status))) +
  ggplot2::geom_point(alpha = 0.68, size = 1.8) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggplot2::facet_wrap(~pair, scales = "free", ncol = 2) +
  ggplot2::labs(title = "Out-of-fold predictions for mutant cell lines",
                subtitle = "Full model D; first CV repeat only",
                x = "Observed gene-effect", y = "Predicted gene-effect",
                color = "Mutation status") +
  ggplot2::theme_minimal(base_size = 12)
save_plot(p5, "Figure_05_OOF_mutant_predictions.png", 10, 8)

# figure 6 displays the ARID1A ARID1B case study
if (!is.null(case_data)) {
  p6 <- ggplot2::ggplot(case_data,
      ggplot2::aes(x = factor(status), y = ARID1B_effect)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.45, width = 0.62) +
    ggplot2::labs(title = "ARID1B dependency by ARID1A mutation status",
                  subtitle = "ARID1B damaging-mutation models excluded",
                  x = "ARID1A mutation category", y = "ARID1B gene-effect score") +
    ggplot2::theme_minimal(base_size = 12)
  save_plot(p6, "Figure_06_ARID1A_ARID1B_case.png", 8.5, 5.5)

  # figure 7 compares ARID1A expression with dependency
  p7 <- ggplot2::ggplot(case_shared,
      ggplot2::aes(x = ARID1A_expression, y = ARID1B_effect)) +
    ggplot2::geom_point(alpha = 0.66) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
    ggplot2::facet_wrap(~status_factor) +
    ggplot2::labs(title = "ARID1A expression vs ARID1B dependency",
                  subtitle = "Mutant categories; shared cancer lineages only",
                  x = "ARID1A expression (log2[TPM + 1])",
                  y = "ARID1B gene-effect score") +
    ggplot2::theme_minimal(base_size = 12)
  save_plot(p7, "Figure_07_ARID1A_expression.png", 10, 5.5)

  # figure 8 checks residuals of the case-study model
  diag <- data.table(fitted = fitted(case_models[["both_expressions"]]),
                     residual = residuals(case_models[["both_expressions"]]))
  p8 <- ggplot2::ggplot(diag, ggplot2::aes(x = fitted, y = residual)) +
    ggplot2::geom_point(alpha = 0.65) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_smooth(method = "loess", se = FALSE) +
    ggplot2::labs(title = "Residual diagnostics: ARID1A–ARID1B case study",
                  x = "Fitted gene-effect", y = "Residual") +
    ggplot2::theme_minimal(base_size = 12)
  save_plot(p8, "Figure_08_case_residuals.png", 8.5, 5)
}

# figure 9 shows generalization to held-out lineages
if (exists("holdout_metrics")) {
  p9 <- ggplot2::ggplot(holdout_metrics,
      ggplot2::aes(x = pair, y = gain_A_to_D, fill = group)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Generalization to held-out cancer lineages",
                  subtitle = "Positive RMSE gain favors full model D",
                  x = NULL, y = "RMSE(A) - RMSE(D)", fill = "Test subset") +
    ggplot2::theme_minimal(base_size = 12)
  save_plot(p9, "Figure_09_lineage_holdout.png", 10, 5.5)
}
message("Figures saved to: ", cfg[["figure_dir"]])
