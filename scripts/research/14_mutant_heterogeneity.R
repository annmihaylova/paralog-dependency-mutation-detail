# study dependency variation among mutation-positive cell lines

mutant_model_formulas <- list(
  U0_status = ~ status_mutant + lineage_group,
  U1_target_expression = ~ status_mutant + expression_B + lineage_group,
  U2_max_AF = ~ status_mutant + max_vaf + lineage_group,
  U3_target_expression_AF = ~ status_mutant + expression_B + max_vaf + lineage_group,
  U4_both_expressions_AF = ~ status_mutant + expression_A + expression_B +
    max_vaf + lineage_group,
  U5_full_variants = ~ status_mutant + expression_A + expression_B + max_vaf +
    copy_number_A + n_lof + n_truncating + n_splice + lineage_group,
  U6_target_expression_sum_AF = ~ status_mutant + expression_B + sum_af + lineage_group,
  U7_status_expression_interaction = ~ status_mutant * expression_B +
    max_vaf + lineage_group
)

mutant_data <- list()
mutant_cohort_summary <- list()
mutant_regressions <- list()
for (i in seq_len(nrow(pairs))) {
  pname <- pairs[["pair"]][i]
  d <- copy(feature_data[[pairs[["stem"]][i]]][status > 0])
  d[, status_mutant := factor(status, levels = c(1, 2))]
  d <- prepare_extended_data(d, mutant_model_formulas, min_rows = 15L)
  mutant_data[[pname]] <- d
  mutant_cohort_summary[[i]] <- d[, .(
    n = .N, n_patients = uniqueN(PatientID),
    mean_effect = mean(effect), median_effect = median(effect),
    sd_effect = sd(effect), weak_above_minus_025 = sum(effect > -0.25),
    strong_at_most_minus_05 = sum(effect <= -0.5)
  ), by = status][, pair := pname]

  # estimate adjusted associations as descriptive effect sizes
  fit <- tryCatch(lm(effect ~ status_mutant + expression_A + expression_B +
                       max_vaf + lineage_group, data = d),
                  error = function(e) NULL)
  if (is.null(fit)) next
  vc <- tryCatch(sandwich::vcovCL(fit, cluster = d[["PatientID"]], type = "HC1"),
                 error = function(e) NULL)
  if (is.null(vc)) next
  test <- tryCatch(lmtest::coeftest(fit, vcov. = vc,
                                    df = uniqueN(d[["PatientID"]]) - 1L),
                   error = function(e) NULL)
  if (is.null(test)) next
  mutant_regressions[[i]] <- data.table(
    pair = pname, n = nrow(d), n_patients = uniqueN(d[["PatientID"]]),
    term = rownames(test), estimate = unname(test[, "Estimate"]),
    robust_SE = unname(test[, "Std. Error"]),
    p_value = unname(test[, ncol(test)])
  )
}
mutant_descriptives <- rbindlist(mutant_cohort_summary)
save_table(mutant_descriptives, "14_mutant_status_heterogeneity.csv")
mutant_adjusted_associations <- rbindlist(mutant_regressions, fill = TRUE)
if (nrow(mutant_adjusted_associations)) {
  mutant_adjusted_associations[term %in% c("expression_A", "expression_B", "max_vaf") &
                                 is.finite(p_value), q_BH_exploratory :=
                                   p.adjust(p_value, "BH")]
}
save_table(mutant_adjusted_associations, "14_mutant_adjusted_associations.csv")

# check predictive performance inside mutant-only cohorts
mutant_cv <- lapply(seq_len(nrow(pairs)), function(i) {
  pname <- pairs[["pair"]][i]
  d <- mutant_data[[pname]]
  run_grouped_cv(d, mutant_model_formulas, pname,
                 cfg[["extended_cv_repeats"]], cfg[["extended_outer_folds"]],
                 cfg[["extended_inner_folds"]],
                 cfg[["seed"]] + 39000L + 100L * i, select_best = TRUE)
})
mutant_predictions <- rbindlist(lapply(mutant_cv, `[[`, "predictions"))
mutant_metrics <- score_extended(mutant_predictions)
mutant_summary <- mutant_metrics[, .(
  n = unique(n)[1L], RMSE_mean = mean(RMSE), MAE_mean = mean(MAE),
  RMSE_SD_across_repeats = sd(RMSE)
), by = .(pair, model, group)]
mutant_comparisons <- data.table(
  from = c("U0_status", "U0_status", "U1_target_expression",
           "U1_target_expression", "U3_target_expression_AF",
           "U3_target_expression_AF", "U6_target_expression_sum_AF",
           "U3_target_expression_AF"),
  to = c("U1_target_expression", "U2_max_AF", "U3_target_expression_AF",
         "U6_target_expression_sum_AF", "U4_both_expressions_AF",
         "U5_full_variants", "U3_target_expression_AF",
         "U7_status_expression_interaction")
)
mutant_gain <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  pname <- pairs[["pair"]][i]
  pred <- mutant_predictions[pair == pname]
  rbindlist(lapply(seq_len(nrow(mutant_comparisons)), function(j) {
    cmp <- mutant_comparisons[j]
    res <- patient_gain(pred, cmp[["from"]], cmp[["to"]],
                        cfg[["extended_bootstrap_reps"]],
                        cfg[["extended_permutation_reps"]],
                        cfg[["seed"]] + 49000L + i * 100L + j)
    res[, `:=`(pair = pname, from = cmp[["from"]], to = cmp[["to"]])]
    res
  }))
}))
mutant_gain[is.finite(p_swap), q_BH_exploratory := p.adjust(p_swap, "BH")]

# flag unusually weak dependence using held-out model residuals
mutant_outliers <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  pname <- pairs[["pair"]][i]
  pred <- mutant_predictions[pair == pname & model == "U4_both_expressions_AF"]
  rows <- pred[, .(observed = observed[1L], predicted = mean(predicted)),
               by = ModelID]
  d <- mutant_data[[pname]]
  rows <- merge(rows, d[, .(ModelID, PatientID, status, expression_A,
                             expression_B, max_vaf, sum_af, status_B,
                             copy_number_A, lineage_group)], by = "ModelID")
  rows[, `:=`(pair = pname, residual = observed - predicted,
               weak_target_effect = observed > -0.25,
               target_expressing_example = expression_B > 2)]
  setorder(rows, -residual)
  rows
}))
save_table(mutant_predictions, "14_mutant_only_oof_predictions.csv")
save_table(rbindlist(lapply(mutant_cv, `[[`, "choices")), "14_mutant_inner_model_choices.csv")
save_table(mutant_summary, "14_mutant_only_cv_summary.csv")
save_table(mutant_gain, "14_mutant_only_paired_gains.csv")
save_table(mutant_outliers, "14_mutant_residual_review.csv")
