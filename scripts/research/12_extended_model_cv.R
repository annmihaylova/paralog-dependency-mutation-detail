# compare fifteen biologically motivated feature representations

extended_model_formulas <- list(
  K_no_mutation = ~ expression_A + lineage_group,
  A_binary = ~ binary + expression_A + lineage_group,
  B_category = ~ status_factor + expression_A + lineage_group,
  C_plus_CN = ~ status_factor + expression_A + copy_number_A + lineage_group,
  D_full = ~ status_factor + expression_A + copy_number_A +
    max_vaf + n_lof + n_truncating + n_splice + lineage_group,
  L_expression_only = ~ expression_A + expression_B + lineage_group,
  E_target_expression = ~ status_factor + expression_A + expression_B + lineage_group,
  G_target_expression_CN = ~ status_factor + expression_A + expression_B +
    copy_number_A + lineage_group,
  F_mutation_interaction = ~ status_factor + expression_A + expression_B +
    binary:expression_B + lineage_group,
  H_mutation_interaction_CN = ~ status_factor + expression_A + expression_B +
    copy_number_A + binary:expression_B + lineage_group,
  N_expression_interaction = ~ status_factor + expression_A * expression_B +
    lineage_group,
  O_expression_interaction_CN = ~ status_factor + expression_A * expression_B +
    copy_number_A + lineage_group,
  M_continuous_AF = ~ binary + sum_af + expression_A + lineage_group,
  I_compact_variant = ~ status_factor + expression_A + expression_B +
    copy_number_A + max_vaf + has_noncanonical_lof + lineage_group,
  J_expression_ratio = ~ status_factor + ratio_AB + lineage_group
)

# group models by patient before assigning folds
patient_folds <- function(status, patient, k, seed) {
  if (length(status) != length(patient)) stop("Status and patient lengths differ")
  patient <- as.character(patient)
  bad <- is.na(patient) | trimws(patient) == ""
  patient[bad] <- paste0("missing_patient_", which(bad))
  by_patient <- data.table(patient = patient, status = as.integer(status))[
    , .(stratum = max(status), n_lines = .N), by = patient]
  if (nrow(by_patient) < k) stop("Not enough independent patients for CV")
  set.seed(seed)
  assigned <- integer(nrow(by_patient))
  load_total <- numeric(k)
  for (s in sort(unique(by_patient[["stratum"]]))) {
    ix <- which(by_patient[["stratum"]] == s)
    ix <- ix[sample.int(length(ix), length(ix))]
    load_status <- numeric(k)
    for (j in ix) {
      choices <- which(load_status == min(load_status))
      target <- choices[which.min(load_total[choices])]
      assigned[j] <- target
      load_status[target] <- load_status[target] + by_patient[["n_lines"]][j]
      load_total[target] <- load_total[target] + by_patient[["n_lines"]][j]
    }
  }
  by_patient[, fold := assigned]
  by_patient[["fold"]][match(patient, by_patient[["patient"]])]
}

# prepare the same complete-case cohort for every candidate model
prepare_extended_data <- function(d, formulas = extended_model_formulas, min_rows = 30L) {
  d <- copy(d)
  d[, binary := as.integer(status > 0)]
  d[, status_factor := factor(status, levels = 0:2)]
  d[, lineage_group := factor(lineage_group)]
  if (!"PatientID" %in% names(d)) d[, PatientID := NA_character_]
  d[, PatientID := as.character(PatientID)]
  d[is.na(PatientID) | trimws(PatientID) == "", PatientID := paste0("unknown_", ModelID)]
  required <- unique(c("ModelID", "PatientID", "status", "effect", "lineage_group",
                       unlist(lapply(formulas, all.vars), use.names = FALSE)))
  absent <- setdiff(required, names(d))
  if (length(absent)) stop("Missing model inputs: ", paste(absent, collapse = ", "))
  n_before <- nrow(d)
  d <- d[complete.cases(d[, ..required])]
  numeric_inputs <- intersect(required, names(d)[vapply(d, is.numeric, logical(1))])
  for (nm in numeric_inputs) d <- d[is.finite(get(nm))]
  setorder(d, ModelID)
  d[, lineage_group := droplevels(lineage_group)]
  d[, status_factor := droplevels(status_factor)]
  attr(d, "n_excluded") <- n_before - nrow(d)
  if (anyDuplicated(d[["ModelID"]])) stop("Duplicate ModelID in extended dataset")
  if (nrow(d) < min_rows) stop("Too few rows in extended dataset")
  d
}

# check only grouping variables used by the current model family
prepare_model_formulas <- function(d, formulas, pair_name) {
  formula_vars <- unique(unlist(lapply(formulas, all.vars), use.names = FALSE))

  if ("status_factor" %in% formula_vars && uniqueN(d[["status_factor"]]) < 2L) {
    stop(pair_name, ": fewer than two observed mutation-status categories remain after complete-case filtering")
  }
  if ("status_mutant" %in% formula_vars && uniqueN(d[["status_mutant"]]) < 2L) {
    stop(pair_name, ": mutant-only models require both status 1 and status 2")
  }
  if ("binary" %in% formula_vars && uniqueN(d[["binary"]]) < 2L) {
    stop(pair_name, ": models using binary mutation status require both mutant and non-mutant groups")
  }
  if ("lineage_group" %in% formula_vars && uniqueN(d[["lineage_group"]]) < 2L) {
    message(pair_name, ": lineage_group has fewer than two observed levels and is omitted from all models")
    formulas <- lapply(formulas, function(f) update(f, . ~ . - lineage_group))
  }
  formulas
}

# build design matrices with explicit factor diagnostics
build_design_matrices <- function(d, formulas, pair_name) {
  x_list <- lapply(names(formulas), function(model_name) {
    f <- formulas[[model_name]]
    vars <- intersect(all.vars(f), names(d))
    factor_vars <- vars[vapply(vars, function(v) is.factor(d[[v]]), logical(1))]
    one_level <- factor_vars[vapply(factor_vars, function(v) {
      uniqueN(d[[v]][!is.na(d[[v]])]) < 2L
    }, logical(1))]
    if (length(one_level)) {
      details <- vapply(one_level, function(v) {
        paste0(v, "={", paste(unique(as.character(d[[v]])), collapse = ","), "}")
      }, character(1))
      stop(pair_name, " / ", model_name,
           ": factor(s) with fewer than two observed levels: ",
           paste(details, collapse = "; "),
           ". Check complete-case filtering and the cohort composition.")
    }
    x <- tryCatch(
      model.matrix(f, data = d),
      error = function(e) {
        factor_info <- if (length(factor_vars)) paste(vapply(factor_vars, function(v) {
          paste0(v, " (", uniqueN(d[[v]][!is.na(d[[v]])]), " observed levels)")
        }, character(1)), collapse = "; ") else "no factor predictors"
        stop(pair_name, " / ", model_name, ": model.matrix failed; ",
             factor_info, "; original error: ", conditionMessage(e), call. = FALSE)
      }
    )
    x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
    varying <- vapply(seq_len(ncol(x)), function(j) uniqueN(x[, j]) > 1L, logical(1))
    x <- x[, varying, drop = FALSE]
    if (!ncol(x) || any(!is.finite(x))) {
      stop(pair_name, " / ", model_name, ": no valid varying predictors remain")
    }
    x
  })
  names(x_list) <- names(formulas)
  x_list
}

# fit every model and choose regularization using training patients only
run_grouped_cv <- function(d, formulas, pair_name, repeats, outer_folds,
                           inner_folds, seed, select_best = FALSE) {
  if (anyNA(d[, .(ModelID, PatientID, status, effect)])) stop("Missing CV keys")
  formulas <- prepare_model_formulas(d, formulas, pair_name)
  x_list <- build_design_matrices(d, formulas, pair_name)
  y <- d[["effect"]]
  predictions <- vector("list", repeats * outer_folds)
  model_choices <- list()
  pos <- 1L
  for (r in seq_len(repeats)) {
    outer <- patient_folds(d[["status"]], d[["PatientID"]], outer_folds,
                           seed + r * 31L)
    for (fold in seq_len(outer_folds)) {
      train <- which(outer != fold)
      test <- which(outer == fold)
      if (length(test) == 0L) next
      inner_k <- min(inner_folds, uniqueN(d[["PatientID"]][train]))
      if (inner_k < 3L) stop("Too few training patients for inner CV")
      inner <- patient_folds(d[["status"]][train], d[["PatientID"]][train],
                             inner_k, seed + 10000L + r * 113L + fold)
      fold_results <- vector("list", length(formulas) + as.integer(select_best))
      cv_errors <- setNames(rep(Inf, length(formulas)), names(formulas))
      names(fold_results) <- c(names(formulas),
                               if (select_best) "SELECTED_INNER_CV" else character())
      for (nm in names(formulas)) {
        x <- x_list[[nm]]
        fit <- glmnet::cv.glmnet(
          x[train, , drop = FALSE], y[train], alpha = 0, family = "gaussian",
          foldid = inner, type.measure = "mse", standardize = TRUE, nlambda = 60L
        )
        cv_errors[nm] <- fit[["cvm"]][which.min(abs(fit[["lambda"]] - fit[["lambda.min"]]))]
        fold_results[[nm]] <- as.numeric(predict(
          fit, newx = x[test, , drop = FALSE], s = "lambda.min"))
      }
      chosen <- names(cv_errors)[which.min(cv_errors)]
      if (select_best) fold_results[["SELECTED_INNER_CV"]] <- fold_results[[chosen]]
      if (select_best) {
        model_choices[[length(model_choices) + 1L]] <- data.table(
          pair = pair_name, repeat_id = r, fold = fold, selected_model = chosen,
          inner_mse = unname(cv_errors[[chosen]]), n_train = length(train),
          n_test = length(test)
        )
      }
      base <- data.table(
        pair = pair_name, repeat_id = r, fold = fold,
        ModelID = d[["ModelID"]][test], PatientID = d[["PatientID"]][test],
        status = d[["status"]][test], observed = y[test]
      )
      predictions[[pos]] <- rbindlist(lapply(names(fold_results), function(nm) {
        ans <- copy(base)
        ans[, `:=`(model = nm, predicted = fold_results[[nm]])]
        ans
      }))
      pos <- pos + 1L
    }
    message(pair_name, ": extended CV repeat ", r, "/", repeats)
  }
  oof <- rbindlist(predictions, use.names = TRUE)
  if (anyNA(oof[["predicted"]])) stop("Missing extended OOF predictions")
  expected <- nrow(d) * repeats * (length(formulas) + as.integer(select_best))
  if (nrow(oof) != expected) stop("Unexpected number of extended OOF predictions")
  list(predictions = oof,
       choices = if (length(model_choices)) rbindlist(model_choices) else data.table())
}

# calculate errors separately for all and mutant-only cell lines
score_extended <- function(pred) {
  groups <- list(all = pred, mutated = pred[status > 0],
                 status_0 = pred[status == 0], status_1 = pred[status == 1],
                 status_2 = pred[status == 2])
  rbindlist(lapply(names(groups), function(g) {
    x <- groups[[g]]
    if (!nrow(x)) return(NULL)
    x[, .(n = uniqueN(ModelID),
           RMSE = sqrt(mean((observed - predicted)^2)),
           MAE = mean(abs(observed - predicted))),
      by = .(pair, repeat_id, model)][, group := g]
  }))
}

extended_audit <- list()
extended_data <- vector("list", nrow(pairs))
names(extended_data) <- pairs[["stem"]]
for (i in seq_len(nrow(pairs))) {
  d <- prepare_extended_data(feature_data[[pairs[["stem"]][i]]])
  extended_data[[i]] <- d
  extended_audit[[i]] <- data.table(
    pair = pairs[["pair"]][i], n = nrow(d),
    n_mutated = sum(d[["status"]] > 0), n_patients = uniqueN(d[["PatientID"]]),
    excluded_vs_step01 = attr(d, "n_excluded")
  )
}
extended_cohort_audit <- rbindlist(extended_audit)
save_table(extended_cohort_audit, "12_extended_cohort_audit.csv")
save_table(data.table(model = names(extended_model_formulas),
                      formula = vapply(extended_model_formulas, function(f)
                        paste(deparse(f), collapse = " "), character(1))),
           "12_model_definitions.csv")

extended_fits <- lapply(seq_len(nrow(pairs)), function(i) {
  run_grouped_cv(extended_data[[i]], extended_model_formulas,
                 pairs[["pair"]][i], cfg[["extended_cv_repeats"]],
                 cfg[["extended_outer_folds"]], cfg[["extended_inner_folds"]],
                 cfg[["seed"]] + i * 2000L, select_best = TRUE)
})
extended_predictions <- rbindlist(lapply(extended_fits, `[[`, "predictions"))
extended_model_choices <- rbindlist(lapply(extended_fits, `[[`, "choices"))
extended_metrics <- score_extended(extended_predictions)
extended_summary <- extended_metrics[, .(
  n = unique(n)[1L], RMSE_mean = mean(RMSE),
  RMSE_SD_across_repeats = sd(RMSE), MAE_mean = mean(MAE)
), by = .(pair, group, model)]
save_table(extended_predictions, "12_extended_oof_predictions.csv")
save_table(extended_model_choices, "12_inner_selected_models.csv")
save_table(extended_metrics, "12_extended_cv_metrics.csv")
save_table(extended_summary, "12_extended_cv_summary.csv")
