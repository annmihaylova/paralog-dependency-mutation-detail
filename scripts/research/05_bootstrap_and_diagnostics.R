# assess the stability of the original A to D comparisons

comparisons <- data.table(
  from = c("A_binary", "B_category", "C_plus_CN", "A_binary"),
  to = c("B_category", "C_plus_CN", "D_full", "D_full")
)

# resample cell-line identities rather than repeated CV folds
bootstrap_gain <- function(pred, subset_name, from, to, n_boot, seed) {
  if (subset_name == "mutated") pred <- pred[status == 1 | status == 2]
  from_col <- paste0("pred_", from)
  to_col <- paste0("pred_", to)
  pred[, `:=`(
    loss_from = (observed - get(from_col))^2,
    loss_to = (observed - get(to_col))^2
  )]
  by_cell <- pred[, .(
    mse_from = mean(loss_from), mse_to = mean(loss_to)
  ), by = ModelID]
  n <- nrow(by_cell)
  gain <- sqrt(mean(by_cell[["mse_from"]])) - sqrt(mean(by_cell[["mse_to"]]))
  set.seed(seed)
  draws <- replicate(n_boot, {
    idx <- sample.int(n, size = n, replace = TRUE)
    sqrt(mean(by_cell[["mse_from"]][idx])) -
      sqrt(mean(by_cell[["mse_to"]][idx]))
  })
  data.table(
    group = subset_name, n_lines = n, from = from, to = to,
    rmse_gain = gain,
    bootstrap_low_2.5 = unname(quantile(draws, 0.025)),
    bootstrap_high_97.5 = unname(quantile(draws, 0.975)),
    fraction_bootstrap_positive = mean(draws > 0)
  )
}

# compare models using averaged errors per cell line
bootstrap_results <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  pair_name <- pairs[["pair"]][i]
  d <- all_predictions[pair == pair_name]
  tab <- rbindlist(lapply(seq_len(nrow(comparisons)), function(j) {
    rbindlist(lapply(seq_along(c("all", "mutated")), function(g) {
      group_name <- c("all", "mutated")[[g]]
      bootstrap_gain(copy(d), group_name,
                     comparisons[["from"]][j], comparisons[["to"]][j],
                     cfg[["bootstrap_reps"]], cfg[["seed"]] + 100L * i + 10L * j + g)
    }))
  }))
  tab[, pair := pair_name]
  tab
}))
save_table(bootstrap_results, "05_cluster_bootstrap_oof_gain.csv")

# inspect correlations between original predictors
feature_correlations <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  d <- feature_data[[pairs[["stem"]][i]]]
  data.table(
    pair = pairs[["pair"]][i], n = nrow(d),
    expression_vs_CN_pearson = cor(d[["expression_A"]], d[["copy_number_A"]]),
    vaf_vs_n_lof_spearman = cor(d[["max_vaf"]], d[["n_lof"]],
                               method = "spearman")
  )
}))
save_table(feature_correlations, "05_predictor_correlations.csv")
message("Conditional cell-line bootstrap comparisons complete.")
