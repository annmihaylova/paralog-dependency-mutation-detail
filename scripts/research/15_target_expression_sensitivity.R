# test whether mutation annotations help when the target paralog is expressed

expressing_formulas <- list(
  T0_expression_only = ~ expression_B + lineage_group,
  T1_status_expression = ~ status_mutant + expression_B + lineage_group,
  T2_plus_max_AF = ~ status_mutant + expression_B + max_vaf + lineage_group,
  T3_plus_sum_AF = ~ status_mutant + expression_B + sum_af + lineage_group,
  T4_plus_expression_A = ~ status_mutant + expression_A + expression_B +
    max_vaf + lineage_group
)
expressing_comparisons <- data.table(
  from = c("T0_expression_only", "T1_status_expression",
           "T1_status_expression", "T2_plus_max_AF"),
  to = c("T1_status_expression", "T2_plus_max_AF",
         "T3_plus_sum_AF", "T4_plus_expression_A")
)

expression_threshold_audit <- list()
expression_threshold_predictions <- list()
expression_threshold_gains <- list()
iter <- 1L
for (i in seq_len(nrow(pairs))) {
  pname <- pairs[["pair"]][i]
  initial <- mutant_data[[pname]]
  for (cutoff in c(-Inf, cfg[["expression_cutoffs"]])) {
    cutoff_label <- if (is.infinite(cutoff)) "all_mutants" else paste0("exprB_gt_", cutoff)
    d <- copy(initial[expression_B > cutoff])
    n1 <- sum(d[["status"]] == 1L)
    n2 <- sum(d[["status"]] == 2L)
    n_patients <- uniqueN(d[["PatientID"]])
    reason <- ""
    if (nrow(d) < cfg[["min_expressing_mutants"]]) reason <- "too few mutant lines"
    if (n1 < 4L || n2 < 4L) reason <- "too few lines in a mutation category"
    if (n_patients < cfg[["extended_outer_folds"]] * 2L) reason <- "too few patients"
    if (uniqueN(d[["lineage_group"]]) < 2L) reason <- "only one lineage group"
    expression_threshold_audit[[iter]] <- data.table(
      pair = pname, cutoff_label = cutoff_label, cutoff = cutoff,
      n = nrow(d), n_status1 = n1, n_status2 = n2,
      n_patients = n_patients, tested = reason == "", skip_reason = reason
    )
    if (reason != "") {
      iter <- iter + 1L
      next
    }
    cv <- run_grouped_cv(d, expressing_formulas, pname,
                         cfg[["extended_cv_repeats"]],
                         cfg[["extended_outer_folds"]],
                         cfg[["extended_inner_folds"]],
                         cfg[["seed"]] + 57000L + i * 700L + iter,
                         select_best = FALSE)
    pred <- cv[["predictions"]]
    pred[, cutoff_label := cutoff_label]
    expression_threshold_predictions[[length(expression_threshold_predictions) + 1L]] <- pred
    tests <- rbindlist(lapply(seq_len(nrow(expressing_comparisons)), function(j) {
      cmp <- expressing_comparisons[j]
      ans <- patient_gain(pred, cmp[["from"]], cmp[["to"]],
                          cfg[["extended_bootstrap_reps"]],
                          cfg[["extended_permutation_reps"]],
                          cfg[["seed"]] + 65000L + i * 700L + iter * 9L + j)
      ans[, `:=`(pair = pname, cutoff_label = cutoff_label,
                 from = cmp[["from"]], to = cmp[["to"]])]
      ans
    }))
    expression_threshold_gains[[length(expression_threshold_gains) + 1L]] <- tests
    iter <- iter + 1L
  }
}

expression_threshold_audit <- rbindlist(expression_threshold_audit)
expression_threshold_predictions <- rbindlist(expression_threshold_predictions)
expression_threshold_gains <- rbindlist(expression_threshold_gains)
if (nrow(expression_threshold_gains)) {
  expression_threshold_gains[is.finite(p_swap),
    q_BH_exploratory := p.adjust(p_swap, "BH")]
}
if (nrow(expression_threshold_predictions)) {
  threshold_cv_metrics <- expression_threshold_predictions[, .(
    n = uniqueN(ModelID), RMSE = sqrt(mean((observed - predicted)^2)),
    MAE = mean(abs(observed - predicted))
  ), by = .(pair, cutoff_label, repeat_id, model)]
  threshold_cv_summary <- threshold_cv_metrics[, .(
    n = unique(n)[1L], RMSE_mean = mean(RMSE), MAE_mean = mean(MAE)
  ), by = .(pair, cutoff_label, model)]
} else {
  threshold_cv_summary <- data.table()
}
save_table(expression_threshold_audit, "15_expression_threshold_audit.csv")
save_table(expression_threshold_predictions, "15_expression_threshold_oof.csv")
save_table(threshold_cv_summary, "15_expression_threshold_cv_summary.csv")
save_table(expression_threshold_gains, "15_expression_threshold_gains.csv")
