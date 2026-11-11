# test selected feature blocks on cancer lineages absent from training

lineage_model_names <- c(
  "B_category", "C_plus_CN", "D_full", "E_target_expression",
  "F_mutation_interaction", "G_target_expression_CN", "H_mutation_interaction_CN"
)
new_lineage_results <- list()
if (isTRUE(cfg[["run_lineage_holdout"]])) {
  for (i in seq_len(nrow(pairs))) {
    pname <- pairs[["pair"]][i]
    d <- extended_data[[i]]
    d[, lineage_group := factor(lineage_group)]
    lineage_contrasts <- list(
      lineage_group = contrasts(d[["lineage_group"]], contrasts = FALSE)
    )
    matrices <- lapply(extended_model_formulas[lineage_model_names], function(f) {
      model.matrix(f, d, contrasts.arg = lineage_contrasts)[, -1L, drop = FALSE]
    })
    lineages <- names(which(table(d[["lineage_group"]]) >= 15L))
    lineages <- setdiff(lineages, "Other")
    for (j in seq_along(lineages)) {
      held <- lineages[[j]]
      test <- which(as.character(d[["lineage_group"]]) == held)
      test_patients <- unique(d[["PatientID"]][test])
      train <- which(as.character(d[["lineage_group"]]) != held &
                       !d[["PatientID"]] %in% test_patients)
      if (length(train) < 100L || length(test) < 15L) next
      if (any(table(factor(d[["status"]][train], levels = 0:2)) < 5L)) next
      inner_k <- min(cfg[["extended_inner_folds"]],
                     uniqueN(d[["PatientID"]][train]))
      if (inner_k < 3L) next
      inner <- patient_folds(d[["status"]][train], d[["PatientID"]][train],
                             inner_k, cfg[["seed"]] + 80000L + i * 200L + j)
      base <- data.table(
        pair = pname, lineage = held, ModelID = d[["ModelID"]][test],
        PatientID = d[["PatientID"]][test], status = d[["status"]][test],
        observed = d[["effect"]][test]
      )
      for (nm in lineage_model_names) {
        x <- matrices[[nm]]
        fit <- glmnet::cv.glmnet(x[train, , drop = FALSE],
                                 d[["effect"]][train], alpha = 0,
                                 family = "gaussian", foldid = inner,
                                 standardize = TRUE, nlambda = 60L,
                                 type.measure = "mse")
        tab <- copy(base)
        tab[, `:=`(model = nm,
                   predicted = as.numeric(predict(fit,
                     newx = x[test, , drop = FALSE], s = "lambda.min")))]
        new_lineage_results[[length(new_lineage_results) + 1L]] <- tab
      }
    }
  }
}
new_lineage_predictions <- rbindlist(new_lineage_results)
if (nrow(new_lineage_predictions)) {
  new_lineage_metrics <- score_extended(new_lineage_predictions[
    , repeat_id := 1L])
  new_lineage_summary <- new_lineage_metrics[, .(
    n = unique(n)[1L], RMSE = mean(RMSE), MAE = mean(MAE)
  ), by = .(pair, group, model)]
} else {
  new_lineage_summary <- data.table()
}
save_table(new_lineage_predictions, "18_heldout_lineage_oof.csv")
save_table(new_lineage_summary, "18_heldout_lineage_summary.csv")

# resample held-out lineages rather than correlated cell lines
lineage_gain_rows <- list()
if (nrow(new_lineage_predictions)) {
  lineage_comparisons <- data.table(
    from = c("B_category", "C_plus_CN", "E_target_expression",
             "G_target_expression_CN"),
    to = c("E_target_expression", "D_full", "F_mutation_interaction",
           "H_mutation_interaction_CN")
  )
  for (pname in pairs[["pair"]]) {
    for (group_name in c("all", "mutated")) {
      pred <- new_lineage_predictions[pair == pname]
      if (group_name == "mutated") pred <- pred[status > 0]
      for (j in seq_len(nrow(lineage_comparisons))) {
        cmp <- lineage_comparisons[j]
        x <- pred[model %in% c(cmp[["from"]], cmp[["to"]])]
        if (uniqueN(x[["model"]]) != 2L) next
        wide <- dcast(x, lineage + ModelID + observed ~ model, value.var = "predicted")
        sums <- wide[, .(n = .N,
          ss_a = sum((observed - get(cmp[["from"]]))^2),
          ss_b = sum((observed - get(cmp[["to"]]))^2)), by = lineage]
        if (nrow(sums) < 3L) next
        gain <- sqrt(sum(sums[["ss_a"]]) / sum(sums[["n"]])) -
                sqrt(sum(sums[["ss_b"]]) / sum(sums[["n"]]))
        set.seed(cfg[["seed"]] + 90000L + j + which(pairs[["pair"]] == pname) * 40L)
        draws <- replicate(cfg[["extended_bootstrap_reps"]], {
          ii <- sample.int(nrow(sums), nrow(sums), replace = TRUE)
          sqrt(sum(sums[["ss_a"]][ii]) / sum(sums[["n"]][ii])) -
            sqrt(sum(sums[["ss_b"]][ii]) / sum(sums[["n"]][ii]))
        })
        lineage_gain_rows[[length(lineage_gain_rows) + 1L]] <- data.table(
          pair = pname, group = group_name, from = cmp[["from"]],
          to = cmp[["to"]], n_lineages = nrow(sums),
          n_lines = sum(sums[["n"]]), rmse_gain = gain,
          ci_low = unname(quantile(draws, 0.025)),
          ci_high = unname(quantile(draws, 0.975))
        )
      }
    }
  }
}
new_lineage_gains <- rbindlist(lineage_gain_rows)
save_table(new_lineage_gains, "18_lineage_cluster_intervals.csv")
