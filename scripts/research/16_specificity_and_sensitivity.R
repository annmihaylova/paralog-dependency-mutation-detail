# compare true paralog targets with cross-pair decoy targets

cross_target_tests <- list()
for (i in seq_len(nrow(pairs))) {
  gene_A <- pairs[["gene_A"]][i]
  true_B <- pairs[["gene_B"]][i]
  d0 <- copy(feature_data[[pairs[["stem"]][i]]])
  for (target in pairs[["gene_B"]]) {
    if (!target %in% names(effects_all) || !target %in% names(expression_all)) next
    e <- effects_all[, .(ModelID, effect_target = get(target))]
    x <- expression_all[, .(ModelID, expression_target = get(target))]
    d <- merge(merge(d0, e, by = "ModelID", all.x = TRUE),
               x, by = "ModelID", all.x = TRUE)
    cols <- c("effect_target", "expression_target", "expression_A",
              "status", "lineage_group", "PatientID")
    d <- d[complete.cases(d[, ..cols])]
    d <- d[is.finite(effect_target) & is.finite(expression_target) &
             is.finite(expression_A)]
    if (nrow(d) < 100L || sd(d[["effect_target"]]) == 0 ||
        sd(d[["expression_target"]]) == 0 || sd(d[["expression_A"]]) == 0) next
    d[, `:=`(
      binary = as.integer(status > 0),
      status_factor = factor(status, levels = 0:2),
      lineage_group = factor(lineage_group),
      effect_z = as.numeric(scale(effect_target)),
      expression_target_z = as.numeric(scale(expression_target)),
      expression_A_z = as.numeric(scale(expression_A))
    )]
    fit <- lm(effect_z ~ status_factor + expression_A_z + expression_target_z +
                binary:expression_target_z + lineage_group, data = d)
    V <- tryCatch(sandwich::vcovCL(fit, cluster = d[["PatientID"]], type = "HC1"),
                  error = function(e) NULL)
    if (is.null(V)) next
    tab <- tryCatch(lmtest::coeftest(fit, vcov. = V,
                     df = uniqueN(d[["PatientID"]]) - 1L), error = function(e) NULL)
    if (is.null(tab)) next
    term <- "expression_target_z:binary"
    if (!term %in% rownames(tab)) term <- "binary:expression_target_z"
    if (!term %in% rownames(tab)) next
    cross_target_tests[[length(cross_target_tests) + 1L]] <- data.table(
      gene_A = gene_A, actual_target = true_B, tested_target = target,
      is_true_partner = target == true_B, n = nrow(d),
      n_mutant = sum(d[["status"]] > 0),
      estimate_standardized = unname(tab[term, "Estimate"]),
      robust_SE = unname(tab[term, "Std. Error"]),
      p_value = unname(tab[term, ncol(tab)])
    )
  }
}
cross_target_results <- rbindlist(cross_target_tests, fill = TRUE)
if (nrow(cross_target_results)) {
  cross_target_results[is.finite(p_value), q_BH_exploratory := p.adjust(p_value, "BH")]
}
save_table(cross_target_results, "16_cross_target_decoy_interactions.csv")

# read technical screen quality without changing the primary cohorts
conf_path <- find_raw("CRISPRConfounders.csv", required = FALSE)
if (!is.null(conf_path)) {
  headers <- names(fread(conf_path, nrows = 0L))
  if ("ScreenNNMD" %in% headers) {
    id_name <- if ("ModelID" %in% headers) "ModelID" else headers[[1L]]
    keep <- c(id_name, "ScreenNNMD")
    if ("IsDefaultEntryForModel" %in% headers) keep <- c(keep, "IsDefaultEntryForModel")
    screen_info <- fread(conf_path, select = keep)
    if (id_name != "ModelID") setnames(screen_info, id_name, "ModelID")
    if ("IsDefaultEntryForModel" %in% names(screen_info)) {
      screen_info <- screen_info[tolower(trimws(as.character(IsDefaultEntryForModel))) == "yes"]
    }
    screen_info[, ScreenNNMD := suppressWarnings(as.numeric(ScreenNNMD))]
    screen_info <- screen_info[, .(screen_nnmd = if (all(is.na(ScreenNNMD)))
                                     NA_real_ else median(ScreenNNMD, na.rm = TRUE)),
                               by = ModelID]
  } else screen_info <- data.table(ModelID = character(), screen_nnmd = numeric())
} else screen_info <- data.table(ModelID = character(), screen_nnmd = numeric())

# test whether associations persist after removing mutated targets and adjusting for screen quality
sensitivity_results <- list()
for (i in seq_len(nrow(pairs))) {
  pname <- pairs[["pair"]][i]
  d <- copy(mutant_data[[pname]])
  d <- merge(d, screen_info, by = "ModelID", all.x = TRUE)
  scenarios <- list(all_mutants = d, target_intact = d[status_B == 0],
                    target_expressing_example = d[expression_B > 2],
                    quality_adjusted = d[!is.na(screen_nnmd)])
  for (scenario in names(scenarios)) {
    x <- scenarios[[scenario]]
    if (nrow(x) < 20L || uniqueN(x[["PatientID"]]) < 15L ||
        uniqueN(x[["lineage_group"]]) < 2L) next
    x[, status_mutant := factor(status, levels = c(1, 2))]
    x[, lineage_group := droplevels(factor(lineage_group))]
    f <- if (scenario == "quality_adjusted")
      effect ~ status_mutant + expression_A + expression_B + max_vaf +
        screen_nnmd + lineage_group
    else effect ~ status_mutant + expression_A + expression_B + max_vaf + lineage_group
    fit <- tryCatch(lm(f, data = x), error = function(e) NULL)
    if (is.null(fit)) next
    V <- tryCatch(sandwich::vcovCL(fit, cluster = x[["PatientID"]], type = "HC1"),
                  error = function(e) NULL)
    if (is.null(V)) next
    tab <- tryCatch(lmtest::coeftest(fit, vcov. = V,
                                     df = uniqueN(x[["PatientID"]]) - 1L),
                    error = function(e) NULL)
    if (is.null(tab)) next
    relevant <- intersect(c("expression_B", "max_vaf", "expression_A"), rownames(tab))
    if (!length(relevant)) next
    sensitivity_results[[length(sensitivity_results) + 1L]] <- data.table(
      pair = pname, scenario = scenario, n = nrow(x),
      n_patients = uniqueN(x[["PatientID"]]), term = relevant,
      estimate = unname(tab[relevant, "Estimate"]),
      robust_SE = unname(tab[relevant, "Std. Error"]),
      p_value = unname(tab[relevant, ncol(tab)])
    )
  }
}
sensitivity_results <- rbindlist(sensitivity_results, fill = TRUE)
if (nrow(sensitivity_results)) {
  sensitivity_results[is.finite(p_value), q_BH_exploratory := p.adjust(p_value, "BH")]
}
save_table(sensitivity_results, "16_mutant_sensitivity_regressions.csv")
