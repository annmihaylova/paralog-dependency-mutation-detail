# describe mutation-status groups and test their associations with paralog dependency
all_descriptives <- list()
all_global <- list()
all_pairwise <- list()
all_adjusted <- list()

for (i in seq_len(nrow(pairs))) {
  stem <- pairs[["stem"]][i]
  pair_name <- pairs[["pair"]][i]
  d <- copy(feature_data[[stem]])

  # summarize gene-effect distributions across mutation-status categories (0, 1, 2)
  desc <- d[, .(
    n = .N,
    mean_effect = mean(effect),
    median_effect = median(effect),
    SD = sd(effect),
    IQR = IQR(effect)
  ), by = status][order(status)]
  desc[, pair := pair_name]
  all_descriptives[[i]] <- desc

  # test group differences using Welch ANOVA and the Kruskal-Wallis rank test
  welch <- oneway.test(effect ~ factor(status), data = d, var.equal = FALSE)
  kruskal <- kruskal.test(effect ~ factor(status), data = d)
  all_global[[i]] <- data.table(
    pair = pair_name,
    welch_F = unname(welch[["statistic"]]),
    welch_df1 = unname(welch[["parameter"]][1]),
    welch_df2 = unname(welch[["parameter"]][2]),
    welch_p = welch[["p.value"]],
    kruskal_H = unname(kruskal[["statistic"]]),
    kruskal_df = unname(kruskal[["parameter"]]),
    kruskal_p = kruskal[["p.value"]]
  )

  # pairwise Welch t-tests with Holm adjustment within each gene pair
  pw <- pairwise.t.test(d[["effect"]], factor(d[["status"]]),
                        pool.sd = FALSE, p.adjust.method = "holm")
  pm <- pw[["p.value"]]
  rows <- list()
  for (r in rownames(pm)) {
    for (c in colnames(pm)) {
      value <- pm[r, c]
      if (!is.na(value)) {
        rows[[length(rows) + 1L]] <- data.table(
          pair = pair_name, status_a = as.integer(c),
          status_b = as.integer(r), holm_p = as.numeric(value)
        )
      }
    }
  }
  all_pairwise[[i]] <- rbindlist(rows)

  # estimate mutation-status associations adjusted for gene A expression and cancer lineage
  d[, status_factor := factor(status, levels = 0:2)]
  d[, lineage_group := factor(lineage_group)]
  fit <- lm(effect ~ status_factor + expression_A + lineage_group, data = d)
  robust <- robust_coefficients(fit)
  robust[, `:=`(
    pair = pair_name,
    CI_low = estimate - qt(0.975, df.residual(fit)) * robust_SE,
    CI_high = estimate + qt(0.975, df.residual(fit)) * robust_SE
  )]
  all_adjusted[[i]] <- robust
}

descriptive_table <- rbindlist(all_descriptives)
global_tests <- rbindlist(all_global)
pairwise_tests <- rbindlist(all_pairwise)
adjusted_coefficients <- rbindlist(all_adjusted)

# apply Benjamini-Hochberg FDR correction across the four gene pairs within each test family
global_tests[, `:=`(
  welch_FDR_4pairs = p.adjust(welch_p, method = "BH"),
  kruskal_FDR_4pairs = p.adjust(kruskal_p, method = "BH")
)]

save_table(descriptive_table, "02_mutation_descriptives.csv")
save_table(global_tests, "02_global_tests.csv")
save_table(pairwise_tests, "02_pairwise_welch_holm.csv")
save_table(adjusted_coefficients, "02_lineage_adjusted_HC3.csv")

# exploratory associations involving target-paralog expression

extra_tests <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  d <- copy(feature_data[[pairs[["stem"]][i]]])
  required <- c(
    "PatientID", "status", "effect",
    "expression_A", "expression_B", "lineage_group"
  )
  d <- d[complete.cases(d[, ..required])]
  d[, status_factor := factor(status, levels = 0:2)]
  d[, lineage_group := factor(lineage_group)]
  d[, expression_B_c := expression_B - mean(expression_B)]
  d[, mutant_x_exprB := as.integer(status > 0) * expression_B_c]
  fit_B <- lm(
    effect ~ status_factor + expression_A +
      expression_B_c + lineage_group,
    data = d
  )
  fit_interaction <- lm(
    effect ~ status_factor + expression_A +
      expression_B_c + mutant_x_exprB + lineage_group,
    data = d
  )
  extract_term <- function(fit, term) {
    V <- sandwich::vcovCL(
      fit,
      cluster = d[["PatientID"]],
      type = "HC3"
    )
    result <- lmtest::coeftest(
      fit,
      vcov. = V,
      df = uniqueN(d[["PatientID"]]) - 1L
    )
    c(
      estimate = result[term, "Estimate"],
      p_value = result[term, "Pr(>|t|)"]
    )
  }
  b <- extract_term(fit_B, "expression_B_c")
  interaction <- extract_term(fit_interaction, "mutant_x_exprB")
  data.table(
    pair = pairs[["pair"]][i],
    test = c("expression_B", "mutation_x_expression_B"),
    n = nrow(d),
    n_mutant = sum(d[["status"]] > 0),
    estimate = c(b[["estimate"]], interaction[["estimate"]]),
    p_value = c(b[["p_value"]], interaction[["p_value"]])
  )
}))

# adjust exploratory p-values across all eight tests
extra_tests[, FDR_BH_8 := p.adjust(p_value, method = "BH")]

save_table(
  extra_tests,
  "02_expressionB_exploratory.csv"
)

message("Descriptive and inferential statistics completed for four pairs.")
