# test prediction on cancer lineages absent from training

lineage_predictions <- list()
if (exists("holdout_metrics")) rm(holdout_metrics)  # clear previous results
if (isTRUE(cfg[["run_lineage_holdout"]])) {
  pos <- 1L
  for (i in seq_len(nrow(pairs))) {
    d <- copy(feature_data[[pairs[["stem"]][i]]])
    pair_name <- pairs[["pair"]][i]
    d[, binary := as.integer(status == 1 | status == 2)]
    d[, status_factor := factor(status, levels = 0:2)]
    d[, lineage_group := factor(lineage_group)]
    # encode every lineage as its own dummy column
    # avoid mapping the held-out lineage to the reference category
    lineage_contrast <- list(lineage_group = contrasts(d[["lineage_group"]], contrasts = FALSE))
    x_basic <- model.matrix(model_formulas[["A_binary"]], d,
                            contrasts.arg = lineage_contrast)[, -1L, drop = FALSE]
    x_full <- model.matrix(model_formulas[["D_full"]], d,
                           contrasts.arg = lineage_contrast)[, -1L, drop = FALSE]
    groups <- names(which(table(d[["lineage_group"]]) >= 15L))
    groups <- groups[groups != "Other"]

    for (j in seq_along(groups)) {
      lineage <- groups[[j]]
      train <- which(as.character(d[["lineage_group"]]) != lineage)
      test <- which(as.character(d[["lineage_group"]]) == lineage)
      if (length(test) == 0L || length(train) < 100L) next
      if (any(table(factor(d[["status"]][train], levels = 0:2)) < 5L)) next
      set.seed(cfg[["seed"]] + 1000L * i + j)
      inner <- stratified_folds(d[["status"]][train], cfg[["cv_inner_folds"]])
      y <- d[["effect"]]
      preds <- list()
      for (nm in c("A_binary", "D_full")) {
        x <- if (nm == "A_binary") x_basic else x_full
        fit <- glmnet::cv.glmnet(
          x[train, , drop = FALSE], y[train], alpha = 0,
          family = "gaussian", foldid = inner, standardize = TRUE,
          type.measure = "mse", nlambda = 60L
        )
        preds[[nm]] <- as.numeric(predict(fit, newx = x[test, , drop = FALSE],
                                         s = "lambda.min"))
      }
      lineage_predictions[[pos]] <- data.table(
        pair = pair_name, lineage = lineage,
        ModelID = d[["ModelID"]][test], status = d[["status"]][test],
        observed = y[test], pred_A = preds[["A_binary"]],
        pred_D = preds[["D_full"]]
      )
      pos <- pos + 1L
    }
    message(pair_name, ": held-out lineage checks done")
  }
}

if (length(lineage_predictions) > 0L) {
  holdout_pred <- rbindlist(lineage_predictions)
  save_table(holdout_pred, "06_heldout_lineage_predictions.csv")

  # report error for all and mutant-only test observations
  holdout_metrics <- rbindlist(list(
    holdout_pred[, .(
      n = .N,
      RMSE_A = sqrt(mean((observed - pred_A)^2)),
      RMSE_D = sqrt(mean((observed - pred_D)^2))
    ), by = pair][, group := "all"],
    holdout_pred[status == 1 | status == 2, .(
      n = .N,
      RMSE_A = sqrt(mean((observed - pred_A)^2)),
      RMSE_D = sqrt(mean((observed - pred_D)^2))
    ), by = pair][, group := "mutated"]
  ))
  holdout_metrics[, gain_A_to_D := RMSE_A - RMSE_D]
  save_table(holdout_metrics, "06_heldout_lineage_summary.csv")
  save_table(holdout_pred[, .(
    n = .N,
    n_mutated = sum(status == 1 | status == 2)
  ), by = .(pair, lineage)], "06_heldout_lineage_sizes.csv")
} else {
  message("Leave-one-lineage-out check skipped (cfg$run_lineage_holdout = FALSE).")
}
