# test whether mutation status is associated with paralog transcript levels

expression_response_rows <- list()
for (i in seq_len(nrow(pairs))) {
  pname <- pairs[["pair"]][i]
  d <- copy(feature_data[[pairs[["stem"]][i]]])
  d[, status_factor := factor(status, levels = 0:2)]
  d[, lineage_group := factor(lineage_group)]
  d[is.na(PatientID) | trimws(as.character(PatientID)) == "",
    PatientID := paste0("unknown_", ModelID)]
  analyses <- list(
    expression_A = expression_A ~ status_factor + copy_number_A + lineage_group,
    expression_B = expression_B ~ status_factor + expression_A + lineage_group
  )
  for (nm in names(analyses)) {
    needed <- all.vars(analyses[[nm]])
    x <- d[complete.cases(d[, ..needed])]
    if (nrow(x) < 30L) next
    fit <- tryCatch(lm(analyses[[nm]], data = x), error = function(e) NULL)
    if (is.null(fit)) next
    V <- tryCatch(sandwich::vcovCL(fit, cluster = x[["PatientID"]], type = "HC1"),
                  error = function(e) NULL)
    if (is.null(V)) next
    tab <- tryCatch(lmtest::coeftest(fit, vcov. = V,
                                     df = uniqueN(x[["PatientID"]]) - 1L),
                    error = function(e) NULL)
    if (is.null(tab)) next
    relevant <- intersect(c("status_factor1", "status_factor2"), rownames(tab))
    if (!length(relevant)) next
    expression_response_rows[[length(expression_response_rows) + 1L]] <- data.table(
      pair = pname, outcome = nm, term = relevant,
      n = nrow(x), n_mutants = sum(x[["status"]] > 0),
      estimate = unname(tab[relevant, "Estimate"]),
      robust_SE = unname(tab[relevant, "Std. Error"]),
      p_value = unname(tab[relevant, ncol(tab)])
    )
  }
}
expression_response <- rbindlist(expression_response_rows, fill = TRUE)
if (nrow(expression_response)) {
  expression_response[is.finite(p_value), q_BH_exploratory := p.adjust(p_value, "BH")]
}
save_table(expression_response, "17_mutation_expression_associations.csv")
