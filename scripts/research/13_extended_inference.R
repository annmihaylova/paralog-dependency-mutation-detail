# estimate paired prediction gains using patient-level resampling

# check that the setup and model definitions have been loaded
if (!exists("cfg", inherits = TRUE) || !exists("pairs", inherits = TRUE) ||
    !exists("save_table", inherits = TRUE)) {
  stop("Run 00_setup.R and 01_prepare_features.R before 13_extended_inference.R", call. = FALSE)
}

# load saved predictions when the model script was run in a previous session
if (!exists("extended_predictions", inherits = TRUE)) {
  prediction_path <- file.path(cfg[["table_dir"]], "12_extended_oof_predictions.csv")
  if (!file.exists(prediction_path)) {
    stop(paste0(
      "Object 'extended_predictions' is missing and no saved prediction file was found. ",
      "Run 12_extended_model_cv.R successfully before running 13_extended_inference.R."
    ), call. = FALSE)
  }
  extended_predictions <- data.table::fread(prediction_path, showProgress = FALSE)
  message("Loaded saved out-of-fold predictions from: ", prediction_path)
}

data.table::setDT(extended_predictions)
required_prediction_columns <- c(
  "pair", "repeat_id", "fold", "ModelID", "PatientID",
  "status", "observed", "model", "predicted"
)
missing_prediction_columns <- setdiff(required_prediction_columns, names(extended_predictions))
if (length(missing_prediction_columns)) {
  stop("Saved predictions are missing columns: ",
       paste(missing_prediction_columns, collapse = ", "), call. = FALSE)
}

# define paired comparisons between biologically motivated models
extended_comparisons <- data.table::data.table(
  from = c("K_no_mutation", "K_no_mutation", "A_binary", "B_category",
           "C_plus_CN", "B_category", "C_plus_CN", "L_expression_only",
           "E_target_expression", "G_target_expression_CN",
           "E_target_expression", "G_target_expression_CN",
           "A_binary", "B_category", "G_target_expression_CN",
           "E_target_expression"),
  to = c("A_binary", "B_category", "B_category", "C_plus_CN",
         "D_full", "E_target_expression", "G_target_expression_CN",
         "E_target_expression", "F_mutation_interaction",
         "H_mutation_interaction_CN", "N_expression_interaction",
         "O_expression_interaction_CN", "M_continuous_AF",
         "M_continuous_AF", "I_compact_variant", "J_expression_ratio")
)
extended_comparisons[, comparison := paste(from, "to", sep = " -> ")]

# verify that all requested models are present for every pair
required_models <- unique(c(extended_comparisons[["from"]], extended_comparisons[["to"]]))
model_coverage <- extended_predictions[, .(models = list(unique(as.character(model)))), by = pair]
coverage_problems <- model_coverage[, .(
  missing_models = paste(setdiff(required_models, models[[1L]]), collapse = "; ")
), by = pair][nzchar(missing_models)]
if (nrow(coverage_problems)) {
  stop("Some required models are absent from the saved predictions: ",
       paste(paste0(coverage_problems[["pair"]], " [", coverage_problems[["missing_models"]], "]"),
             collapse = "; "),
       ". Re-run 12_extended_model_cv.R successfully.", call. = FALSE)
}

# calculate conditional intervals from existing held-out predictions
patient_gain <- function(pred, from, to, n_boot, n_perm, seed) {
  x <- pred[model %in% c(from, to)]
  if (uniqueN(x[["model"]]) != 2L) stop("Missing comparison models: ", from, " and/or ", to)
  wide <- data.table::dcast(
    x,
    pair + repeat_id + fold + ModelID + PatientID + status + observed ~ model,
    value.var = "predicted"
  )
  if (anyNA(wide[, c(from, to), with = FALSE])) stop("Missing paired predictions")
  by_line <- wide[, .(
    err_from = mean((observed - get(from))^2),
    err_to = mean((observed - get(to))^2)
  ), by = .(PatientID, ModelID)]
  by_patient <- by_line[, .(
    n_lines = .N, ss_from = sum(err_from), ss_to = sum(err_to)
  ), by = PatientID]
  n <- nrow(by_patient)
  if (n < 5L) return(data.table::data.table(
    n_patients = n, n_lines = nrow(by_line),
    rmse_from = NA_real_, rmse_to = NA_real_, rmse_gain = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_, p_swap = NA_real_
  ))
  ss_a <- by_patient[["ss_from"]]
  ss_b <- by_patient[["ss_to"]]
  weights <- by_patient[["n_lines"]]
  calc_gain <- function(a, b, w) sqrt(sum(a) / sum(w)) - sqrt(sum(b) / sum(w))
  gain <- calc_gain(ss_a, ss_b, weights)
  set.seed(seed)
  bootstrap <- replicate(n_boot, {
    idx <- sample.int(n, n, replace = TRUE)
    calc_gain(ss_a[idx], ss_b[idx], weights[idx])
  })
  swapped <- replicate(n_perm, {
    mask <- sample(c(FALSE, TRUE), n, replace = TRUE)
    calc_gain(ifelse(mask, ss_b, ss_a), ifelse(mask, ss_a, ss_b), weights)
  })
  data.table::data.table(
    n_patients = n, n_lines = nrow(by_line),
    rmse_from = sqrt(sum(ss_a) / sum(weights)),
    rmse_to = sqrt(sum(ss_b) / sum(weights)),
    rmse_gain = gain,
    ci_low = unname(stats::quantile(bootstrap, 0.025)),
    ci_high = unname(stats::quantile(bootstrap, 0.975)),
    p_swap = (1 + sum(abs(swapped) >= abs(gain))) / (n_perm + 1)
  )
}

# compare models using paired patient-level bootstrap and sign-flip tests
extended_tests <- data.table::rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  pname <- pairs[["pair"]][i]
  pred <- extended_predictions[pair == pname]
  data.table::rbindlist(lapply(seq_len(nrow(extended_comparisons)), function(j) {
    cmp <- extended_comparisons[j]
    data.table::rbindlist(lapply(c("all", "mutated"), function(group_name) {
      x <- if (group_name == "mutated") pred[status > 0] else pred
      res <- patient_gain(
        x, cmp[["from"]], cmp[["to"]],
        cfg[["extended_bootstrap_reps"]],
        cfg[["extended_permutation_reps"]],
        cfg[["seed"]] + i * 20000L + j * 13L +
          ifelse(group_name == "mutated", 1L, 0L)
      )
      res[, `:=`(pair = pname, group = group_name,
                 from = cmp[["from"]], to = cmp[["to"]],
                 comparison = cmp[["comparison"]])]
      res
    }), fill = TRUE)
  }), fill = TRUE)
}), fill = TRUE)

# adjust exploratory p-values across all pairs and comparisons within each evaluation group
if (nrow(extended_tests)) {
  extended_tests[!is.na(p_swap), q_BH_exploratory := p.adjust(p_swap, "BH"), by = group]
}
save_table(extended_comparisons, "13_comparison_definitions.csv")
save_table(extended_tests, "13_patient_bootstrap_comparisons.csv")

# use the same complete-case feature data when available
get_pair_features <- function(i) {
  stem <- pairs[["stem"]][i]
  if (exists("extended_data", inherits = TRUE) && length(extended_data) >= i) {
    return(data.table::copy(extended_data[[i]]))
  }
  if (exists("feature_data", inherits = TRUE) && stem %in% names(feature_data)) {
    return(data.table::copy(feature_data[[stem]]))
  }
  feature_path <- file.path(cfg[["feature_dir"]], paste0(stem, "_features.csv"))
  if (file.exists(feature_path)) {
    return(data.table::fread(feature_path, showProgress = FALSE))
  }
  stop("Feature data are unavailable for ", stem,
       ". Run 01_prepare_features.R or restore the saved feature table.", call. = FALSE)
}

# inspect predictor redundancy inside and outside the mutant group
extended_feature_correlations <- data.table::rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  d <- get_pair_features(i)
  corr_vars <- c("max_vaf", "sum_af", "n_lof", "expression_A", "expression_B")
  missing_vars <- setdiff(corr_vars, names(d))
  if (length(missing_vars)) {
    stop(pairs[["pair"]][i], ": missing correlation variables: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  d <- d[stats::complete.cases(d[, ..corr_vars])]
  data.table::rbindlist(lapply(c("all", "mutated"), function(g) {
    x <- if (g == "mutated") d[status > 0] else d
    safe_cor <- function(a, b) {
      ok <- is.finite(a) & is.finite(b)
      if (sum(ok) < 3L || data.table::uniqueN(a[ok]) < 2L ||
          data.table::uniqueN(b[ok]) < 2L) return(NA_real_)
      suppressWarnings(stats::cor(a[ok], b[ok], method = "spearman"))
    }
    data.table::data.table(
      pair = pairs[["pair"]][i], group = g, n = nrow(x),
      max_vaf_vs_sum_af = safe_cor(x[["max_vaf"]], x[["sum_af"]]),
      max_vaf_vs_n_lof = safe_cor(x[["max_vaf"]], x[["n_lof"]]),
      expr_A_vs_expr_B = safe_cor(x[["expression_A"]], x[["expression_B"]])
    )
  }), fill = TRUE)
}), fill = TRUE)
save_table(extended_feature_correlations, "13_feature_correlations_by_group.csv")

message("Extended inference complete: paired model comparisons and feature correlations saved")
