# summarize exploratory covariate checks from steps 9 and 10

f2 <- function(x, d = 3L) if (is.na(x)) "NA" else formatC(x, digits = d, format = "f")
fp <- function(p) if (is.na(p)) "NA" else if (p < 1e-3) formatC(p, digits = 1, format = "e") else formatC(p, digits = 3, format = "f")

lines <- c(
  "# Covariate checks — summary",
  "",
  paste0("Data source: **", cfg[["release"]], "**. Exploratory analysis plan: `docs/ANALYSIS_PLAN.md`."),
  "Base model (as in step 02): `effect ~ status(0/1/2) + expression_A + lineage_group`, HC3 standard errors.",
  "",
  "## Covariates analysed",
  "| Covariate | Tier | Source | Available |",
  "|---|---:|---|---|"
)
for (r in seq_len(nrow(cov_spec))) {
  lines <- c(lines, sprintf("| %s | %d | %s | %s |", cov_spec[["id"]][r], cov_spec[["tier"]][r],
                            cov_spec[["source"]][r], if (cov_spec[["available"]][r]) "yes" else "**no — file absent**"))
}

# compare status-two coefficients before and after adjustment
lines <- c(lines, "", "## Is the status effect robust to adjustment?",
           "Status-2 coefficient in the base model and after adjusting for all Tier 1–2 covariates together (same cell lines).", "",
           "| Pair | Coefficient | n | Base | Full model [95% CI] | Change |",
           "|---|---|---:|---:|---|---:|")
for (r in seq_len(nrow(full_model))) {
  x <- full_model[r]
  lines <- c(lines, sprintf("| %s | %s | %d | %s | %s [%s, %s] | %s%% |", x[["pair"]],
                            sub("status_factor", "status ", x[["coefficient"]]), x[["n"]],
                            f2(x[["beta_base"]]), f2(x[["beta_full"]]), f2(x[["CI_low_full"]]),
                            f2(x[["CI_high_full"]]), f2(x[["pct_change"]], 1)))
}
flags <- cie[confounder_flag == TRUE]
lines <- c(lines, "", sprintf("Covariates passing the exploratory change-in-estimate rule (>10%% and >0.02): **%d** of %d coefficient checks.",
                              nrow(flags), nrow(cie)))
if (nrow(flags)) {
  lines <- c(lines, "", "| Pair | Covariate | Coefficient | Base | Adjusted | Change |", "|---|---|---|---:|---:|---:|")
  for (r in seq_len(nrow(flags))) {
    x <- flags[r]
    lines <- c(lines, sprintf("| %s | %s | %s | %s | %s | %s%% |", x[["pair"]], x[["covariate"]],
                              sub("status_factor", "status ", x[["coefficient"]]), f2(x[["beta_base"]]),
                              f2(x[["beta_adjusted"]]), f2(x[["pct_change"]], 1)))
  }
}

# report strong group imbalances
strong <- balance[imbalance == "strong (>= 0.25)"]
lines <- c(lines, "", "## Strong imbalance between mutation-status groups (|SMD| >= 0.25)",
           "Imbalance alone is not confounding: it matters only if the covariate also changes the status coefficient (table above).", "",
           "| Pair | Covariate | Mean/proportion in 0 / 1 / 2 | Max abs SMD |", "|---|---|---|---:|")
for (r in seq_len(nrow(strong))) {
  x <- strong[r]
  lines <- c(lines, sprintf("| %s | %s | %s / %s / %s | %s |", x[["pair"]], x[["covariate"]],
                            f2(x[["value_0"]], 2), f2(x[["value_1"]], 2), f2(x[["value_2"]], 2), f2(x[["max_abs_SMD"]], 2)))
}

# report independent covariate associations
sig_main <- main_effects[!is.na(FDR_BH) & FDR_BH < 0.05]
lines <- c(lines, "", "## Covariates associated with the dependency itself (BH FDR < 0.05)",
           "| Pair | Covariate | Estimate [95% CI] | Per SD | FDR |", "|---|---|---|---:|---:|")
for (r in seq_len(nrow(sig_main))) {
  x <- sig_main[r]
  lines <- c(lines, sprintf("| %s | %s | %s [%s, %s] | %s | %s |", x[["pair"]], x[["covariate"]], f2(x[["estimate"]]),
                            f2(x[["CI_low"]]), f2(x[["CI_high"]]), f2(x[["estimate_per_SD"]]), fp(x[["FDR_BH"]])))
}

# report mutation-status interactions
it <- interactions[tested == TRUE]
lines <- c(lines, "", "## Interactions with mutation status",
           sprintf("Tested: %d; not tested because of too few mutant lines: %d.",
                   nrow(it), nrow(interactions[tested == FALSE])), "",
           "| Pair | Covariate | Form | Family | Wald chi2 (df) | p | FDR |", "|---|---|---|---|---|---:|---:|")
for (r in seq_len(nrow(it))) {
  x <- it[r]
  if (!(x[["family"]] == "literature_motivated" || (!is.na(x[["FDR_BH"]]) && x[["FDR_BH"]] < 0.1))) next
  lines <- c(lines, sprintf("| %s | %s | %s | %s | %s (%d) | %s | %s |", x[["pair"]], x[["covariate"]], x[["form"]],
                            x[["family"]], f2(x[["wald_chi2"]], 2), as.integer(x[["wald_df"]]), fp(x[["p_value"]]), fp(x[["FDR_BH"]])))
}
lines <- c(lines, "", "Shown: all literature-motivated tests, plus exploratory tests with FDR < 0.1.")

# exploratory screen
rep_tab <- screen[replicated == TRUE]
lines <- c(lines, "", "## Exploratory screen (200 lineage splits)",
           sprintf("Replicated in >= 50%% of splits: %d of %d tests.", nrow(rep_tab), nrow(screen)), "",
           "| Pair | Covariate | Test | Estimate (all lines) | Replication rate |", "|---|---|---|---:|---:|")
for (r in seq_len(nrow(rep_tab))) {
  x <- rep_tab[r]
  lines <- c(lines, sprintf("| %s | %s | %s | %s | %s |", x[["pair"]], x[["covariate"]], x[["test"]],
                            f2(x[["estimate_all"]]), f2(x[["replication_rate"]], 2)))
}

lines <- c(lines, "", "## Figures",
           "- `figures/Figure_10_covariate_balance.png` — SMD between status groups",
           "- `figures/Figure_11_change_in_estimate.png` — status-2 coefficient before/after adjustment",
           "", "## Cautions",
           "- All associations are observational (cell lines, DepMap); no causal claims.",
           "- Paralog expression can be a mediator of the mutation effect, not only a confounder; adjusting for it can remove part of a real effect.",
           "- Interaction tests have low power (29–139 mutant lines); a non-significant interaction is not evidence of equal effects.",
           "- Covariates whose source file is absent were not analysed; see the table at the top.")

writeLines(lines, file.path(cfg[["out_dir"]], "COVARIATE_SUMMARY.md"))
message("Covariate summary written: ", file.path(cfg[["out_dir"]], "COVARIATE_SUMMARY.md"))
