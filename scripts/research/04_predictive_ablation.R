# compare incremental mutation features for CRISPR dependency prediction

model_formulas <- list(
  A_binary = ~ binary + expression_A + lineage_group,
  B_category = ~ status_factor + expression_A + lineage_group,
  C_plus_CN = ~ status_factor + expression_A + copy_number_A + lineage_group,
  D_full = ~ status_factor + expression_A + copy_number_A +
    max_vaf + n_lof + n_truncating + n_splice + lineage_group
)

# tune ridge only inside training folds of nested CV
cv_pair <- function(d, pair_name) {
  d <- copy(d)
  d[, binary := as.integer(status == 1 | status == 2)]
  d[, status_factor := factor(status, levels = 0:2)]
  d[, lineage_group := factor(lineage_group)]
  xmat <- lapply(model_formulas, function(f) {
    x <- model.matrix(f, data = d)[, -1L, drop = FALSE]
    if (any(!is.finite(x))) stop("Invalid model matrix: ", pair_name)
    x
  })
  y <- d[["effect"]]
  statuses <- d[["status"]]
  out <- vector("list", cfg[["cv_repeats"]] * cfg[["cv_outer_folds"]])
  pos <- 1L

  for (r in seq_len(cfg[["cv_repeats"]])) {
    set.seed(cfg[["seed"]] + r)
    outer <- stratified_folds(statuses, cfg[["cv_outer_folds"]])
    for (f in seq_len(cfg[["cv_outer_folds"]])) {
      train <- which(outer != f)
      test <- which(outer == f)
      set.seed(cfg[["seed"]] + 1000L + 100L * r + f)
      inner <- stratified_folds(statuses[train], cfg[["cv_inner_folds"]])
      pred <- matrix(NA_real_, nrow = length(test), ncol = length(xmat),
                     dimnames = list(NULL, names(xmat)))
      for (nm in names(xmat)) {
        x <- xmat[[nm]]
        fit <- glmnet::cv.glmnet(
          x[train, , drop = FALSE], y[train],
          alpha = 0, family = "gaussian", foldid = inner,
          standardize = TRUE, type.measure = "mse", nlambda = 60L
        )
        pred[, nm] <- as.numeric(predict(fit, newx = x[test, , drop = FALSE],
                                        s = "lambda.min"))
      }
      tab <- data.table(
        pair = pair_name, repeat_id = r, fold = f,
        ModelID = d[["ModelID"]][test], status = statuses[test],
        observed = y[test]
      )
      for (nm in colnames(pred)) tab[, (paste0("pred_", nm)) := pred[, nm]]
      out[[pos]] <- tab
      pos <- pos + 1L
    }
    message(pair_name, ": CV repeat ", r, "/", cfg[["cv_repeats"]])
  }
  result <- rbindlist(out)
  if (nrow(result) != nrow(d) * cfg[["cv_repeats"]]) stop("Incorrect OOF rows")
  if (anyNA(result)) stop("Missing OOF prediction")
  result
}

# generate out-of-fold predictions for four pairs
all_predictions <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  cv_pair(feature_data[[pairs[["stem"]][i]]], pairs[["pair"]][i])
}))
save_table(all_predictions, "04_oof_predictions.csv")

pred_cols <- grep("^pred_", names(all_predictions), value = TRUE)
cv_long <- melt(all_predictions,
                id.vars = c("pair", "repeat_id", "fold", "ModelID", "status", "observed"),
                measure.vars = pred_cols, variable.name = "model", value.name = "predicted")
cv_long[, model := sub("^pred_", "", as.character(model))]

# evaluate prediction error for all and mutation-positive lines
score_group <- function(x, group_label) {
  x[, .(
    n = uniqueN(ModelID),
    RMSE = sqrt(mean((observed - predicted)^2)),
    MAE = mean(abs(observed - predicted))
  ), by = .(pair, repeat_id, model)][, group := group_label]
}
cv_metrics <- rbindlist(list(
  score_group(cv_long, "all"),
  score_group(cv_long[status == 1 | status == 2], "mutated"),
  score_group(cv_long[status == 0], "status_0"),
  score_group(cv_long[status == 1], "status_1"),
  score_group(cv_long[status == 2], "status_2")
))
setcolorder(cv_metrics, c("pair", "repeat_id", "model", "group", "n", "RMSE", "MAE"))

# report variation across repeats without interpreting it as a confidence interval
cv_summary <- cv_metrics[, .(
  n = unique(n)[1L], RMSE_mean = mean(RMSE),
  RMSE_SD_across_repeats = sd(RMSE), MAE_mean = mean(MAE),
  MAE_SD_across_repeats = sd(MAE)
), by = .(pair, group, model)]

# compare sequential feature blocks on identical observations
cv_wide <- dcast(cv_metrics, pair + group + repeat_id ~ model, value.var = "RMSE")
cv_wide[, `:=`(
  B_minus_A = A_binary - B_category,
  C_minus_B = B_category - C_plus_CN,
  D_minus_C = C_plus_CN - D_full,
  D_minus_A = A_binary - D_full,
  relative_D_gain_pct = 100 * (A_binary - D_full) / A_binary
)]
ablation_steps <- cv_wide[, .(
  B_minus_A = mean(B_minus_A), C_minus_B = mean(C_minus_B),
  D_minus_C = mean(D_minus_C), D_minus_A = mean(D_minus_A),
  relative_D_gain_pct = mean(relative_D_gain_pct),
  D_minus_A_SD_across_repeats = sd(D_minus_A)
), by = .(pair, group)]

save_table(cv_metrics, "04_cv_metrics_each_repeat.csv")
save_table(cv_summary, "04_cv_summary.csv")
save_table(ablation_steps, "04_ablation_steps.csv")
message("Predictive ablation complete: ", nrow(all_predictions), " OOF rows.")
