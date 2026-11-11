# verify grouped folds and model definitions on synthetic data

source(file.path("scripts", "research", "00_setup.R"))
old_out <- cfg[["out_dir"]]
cfg[["out_dir"]] <- tempfile("depmap_smoke_")
cfg[["table_dir"]] <- file.path(cfg[["out_dir"]], "tables")
cfg[["figure_dir"]] <- file.path(cfg[["out_dir"]], "figures")
dir.create(cfg[["table_dir"]], recursive = TRUE)
dir.create(cfg[["figure_dir"]], recursive = TRUE)
cfg[["extended_cv_repeats"]] <- 1L
cfg[["extended_outer_folds"]] <- 3L
cfg[["extended_inner_folds"]] <- 3L
cfg[["extended_bootstrap_reps"]] <- 50L
cfg[["extended_permutation_reps"]] <- 50L

set.seed(41)
n <- 120L
patient <- paste0("P", rep(seq_len(n / 2L), each = 2L))
patient_status <- sample(c(0L, 1L, 2L), n / 2L, replace = TRUE,
                         prob = c(0.6, 0.2, 0.2))
status <- rep(patient_status, each = 2L)
eA <- rnorm(n, 4, 1)
eB <- rnorm(n, 3.5, 1)
vaf <- ifelse(status == 0L, 0, pmin(1, runif(n, 0.25, 0.95)))
lof <- as.integer(status > 0L)
sum_af <- ifelse(status == 0L, 0, vaf)
synthetic <- data.table(
  ModelID = paste0("ACH-", sprintf("%06d", seq_len(n))),
  PatientID = patient, status = status, status_B = 0L,
  effect = -0.2 - 0.3 * as.integer(status == 2L) - 0.06 * eB + rnorm(n, 0, 0.25),
  expression_A = eA, expression_B = eB, ratio_AB = eA - eB,
  copy_number_A = runif(n, 0.5, 2), copy_number_B = runif(n, 0.5, 2),
  max_vaf = vaf, sum_af = sum_af, n_lof = lof,
  n_truncating = lof, n_splice = 0L, has_noncanonical_lof = 0L,
  lineage_group = rep(c("lung", "breast", "other"), length.out = n)
)

pairs <- pairs[1L]
feature_data <- setNames(list(synthetic), pairs[["stem"]])
source(file.path("scripts", "research", "12_extended_model_cv.R"))
stopifnot(length(extended_model_formulas) == 15L)
stopifnot(nrow(extended_cohort_audit) == 1L)
stopifnot(all(names(extended_model_formulas) %in% extended_summary[["model"]]))
stopifnot("SELECTED_INNER_CV" %in% extended_summary[["model"]])

selected <- extended_predictions[model == "SELECTED_INNER_CV"]
stopifnot(nrow(selected) == n)
stopifnot(all(selected[, uniqueN(fold) == 1L, by = PatientID][["V1"]]))
stopifnot(all(is.finite(selected[["predicted"]])))

source(file.path("scripts", "research", "13_extended_inference.R"))
stopifnot(nrow(extended_tests) > 0L)
stopifnot(all(is.finite(extended_tests[["rmse_gain"]])))
unlink(cfg[["out_dir"]], recursive = TRUE)
message("Synthetic checks passed")
