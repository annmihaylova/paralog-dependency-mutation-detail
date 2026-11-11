# detailed case study of ARID1A mutation and ARID1B dependency

case_available <- all(vapply(c("Model.csv", "CRISPRGeneEffect.csv",
                                "OmicsSomaticMutationsMatrixDamaging.csv"),
  function(x) !is.null(find_raw(x, FALSE)), logical(1))) &&
  !is.null(find_raw(c("OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
                      "OmicsExpressionProteinCodingGenesTPMLogp1.csv"), FALSE))

case_data <- NULL
case_shared <- NULL
case_models <- NULL

if (!case_available) {
  warning("Case study skipped: original CRISPR, mutation, expression or model CSV is absent.")
  writeLines("Case study not run: original data/raw CSVs are required.",
             file.path(cfg[["out_dir"]], "case_study_not_run.txt"))
} else {
  m <- read_models(c("ModelID", "CellLineName", "OncotreeLineage"))
  effect <- read_gene_matrix(find_raw("CRISPRGeneEffect.csv"), "ARID1B")
  mutation <- read_gene_matrix(find_raw("OmicsSomaticMutationsMatrixDamaging.csv"),
                               c("ARID1A", "ARID1B"))
  expression <- read_gene_matrix(find_raw(c(
    "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv")), c("ARID1A", "ARID1B"))
  setnames(effect, "ARID1B", "ARID1B_effect")
  setnames(mutation, c("ARID1A", "ARID1B"), c("ARID1A_mutation", "ARID1B_mutation"))
  setnames(expression, c("ARID1A", "ARID1B"),
           c("ARID1A_expression", "ARID1B_expression"))

  # merge measurements and exclude damaging mutations of ARID1B
  case_data <- Reduce(function(x, y) merge(x, y, by = "ModelID"),
                      list(m, effect, mutation, expression))
  case_data <- case_data[ARID1B_mutation == 0 &
                         (ARID1A_mutation == 0 | ARID1A_mutation == 1 |
                          ARID1A_mutation == 2) &
                         complete.cases(case_data[, .(
                           ARID1B_effect, ARID1A_expression,
                           ARID1B_expression, OncotreeLineage)])]
  case_data[, status := as.integer(ARID1A_mutation)]
  if (anyDuplicated(case_data[["ModelID"]])) stop("Duplicate case study ModelID")
  save_table(case_data[, .(ModelID, CellLineName, OncotreeLineage,
                           ARID1A_mutation, ARID1B_mutation,
                           ARID1A_expression, ARID1B_expression,
                           ARID1B_effect)], "03_case_study_cohort.csv")

  # summarize mutation groups and run global tests
  save_table(case_data[, .(
    n = .N, mean = mean(ARID1B_effect), median = median(ARID1B_effect),
    SD = sd(ARID1B_effect), IQR = IQR(ARID1B_effect)
  ), by = status][order(status)], "03_case_study_descriptives.csv")
  w <- oneway.test(ARID1B_effect ~ factor(status), data = case_data)
  k <- kruskal.test(ARID1B_effect ~ factor(status), data = case_data)
  save_table(data.table(
    test = c("Welch ANOVA", "Kruskal-Wallis"),
    statistic = c(unname(w[["statistic"]]), unname(k[["statistic"]])),
    p_value = c(w[["p.value"]], k[["p.value"]])
  ), "03_case_study_global_tests.csv")
  pw <- pairwise.t.test(case_data[["ARID1B_effect"]], factor(case_data[["status"]]),
                        pool.sd = FALSE, p.adjust.method = "holm")
  pmat <- pw[["p.value"]]
  ptab <- rbindlist(lapply(rownames(pmat), function(r) {
    rbindlist(lapply(colnames(pmat), function(c) {
      if (is.na(pmat[r, c])) return(NULL)
      data.table(group_a = as.integer(c), group_b = as.integer(r),
                 holm_p = as.numeric(pmat[r, c]))
    }))
  }))
  save_table(ptab, "03_case_study_pairwise_holm.csv")

  # compare statuses 1 and 2 within shared cancer lineages
  case_shared <- case_data[status == 1 | status == 2]
  overlap <- case_shared[, .(number_of_statuses = uniqueN(status)),
                         by = OncotreeLineage][number_of_statuses == 2, OncotreeLineage]
  case_shared <- case_shared[OncotreeLineage %chin% overlap]
  case_shared[, status_factor := factor(status, levels = c(1, 2))]
  case_shared[, ARID1A_expression_c := ARID1A_expression - mean(ARID1A_expression)]
  case_shared[, ARID1B_expression_c := ARID1B_expression - mean(ARID1B_expression)]
  save_table(case_shared[, .(ModelID, CellLineName, status, OncotreeLineage,
                             ARID1A_expression, ARID1B_expression, ARID1B_effect)],
             "03_case_study_shared_lineages.csv")

  # adjust for lineage and both gene expressions
  case_models <- list(
    intercept_only = lm(ARID1B_effect ~ 1, data = case_shared),
    lineage_mutation = lm(ARID1B_effect ~ status_factor +
                            factor(OncotreeLineage), data = case_shared),
    arid1a_expression = lm(ARID1B_effect ~
      ARID1A_expression_c * status_factor + factor(OncotreeLineage),
      data = case_shared),
    both_expressions = lm(ARID1B_effect ~
      ARID1A_expression_c * status_factor + ARID1B_expression_c +
        factor(OncotreeLineage), data = case_shared),
    both_interactions = lm(ARID1B_effect ~
      ARID1A_expression_c * status_factor +
        ARID1B_expression_c * status_factor + factor(OncotreeLineage),
      data = case_shared)
  )
  case_coef <- rbindlist(lapply(names(case_models)[-1L], function(name) {
    d <- robust_coefficients(case_models[[name]])
    d[, model := name]
    d
  }), fill = TRUE)
  save_table(case_coef, "03_case_study_robust_HC3.csv")

  # estimate held-out error using leave-one-out cross-validation
  case_cv <- rbindlist(lapply(names(case_models), function(name) {
    fit <- case_models[[name]]
    h <- hatvalues(fit)
    pred_error <- residuals(fit) / (1 - h)
    # omit leverage-one observations from leave-one-out evaluation
    pred_error[h > 1 - 1e-8] <- NA_real_
    n_bad <- sum(!is.finite(pred_error))
    if (n_bad > 0L) warning(name, ": ", n_bad, " lines with leverage 1 skipped in LOOCV")
    pred_error <- pred_error[is.finite(pred_error)]
    data.table(
      model = name, n = length(pred_error), n_skipped_leverage1 = n_bad,
      LOOCV_RMSE = sqrt(mean(pred_error^2)),
      LOOCV_MAE = mean(abs(pred_error)),
      adjusted_R2 = if (name == "intercept_only") NA_real_
                    else summary(fit)[["adj.r.squared"]]
    )
  }))
  save_table(case_cv, "03_case_study_LOOCV.csv")

  # assess sensitivity to the three most influential cell lines
  fit <- case_models[["both_expressions"]]
  cooks <- cooks.distance(fit)
  top <- order(cooks, decreasing = TRUE)[1:3]
  save_table(case_shared[top, .(ModelID, CellLineName, OncotreeLineage,
                                 status, ARID1B_effect)][,
                  cooks_distance := as.numeric(cooks[top])],
             "03_case_study_influential_cells.csv")
  fit_reduced <- update(fit, data = case_shared[-top])
  case_sensitivity <- rbindlist(list(
    robust_coefficients(fit)[, scenario := "all"],
    robust_coefficients(fit_reduced)[, scenario := "without_top3"]
  ))
  save_table(case_sensitivity, "03_case_study_sensitivity.csv")
  message("Case study: ", nrow(case_data), " eligible, ", nrow(case_shared),
          " shared-lineage mutant lines.")
}

