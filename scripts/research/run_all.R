# run the original analyses and exploratory extensions

script_dir <- file.path("scripts", "research")

source(file.path(script_dir, "00_setup.R"))
source(file.path(script_dir, "01_prepare_features.R"))
source(file.path(script_dir, "02_descriptive_tests.R"))
source(file.path(script_dir, "03_arid1a_case_study.R"))
source(file.path(script_dir, "04_predictive_ablation.R"))
source(file.path(script_dir, "05_bootstrap_and_diagnostics.R"))
source(file.path(script_dir, "06_lineage_holdout.R"))
source(file.path(script_dir, "07_figures.R"))
source(file.path(script_dir, "08_results_summary.R"))
source(file.path(script_dir, "09_covariate_checks.R"))

# run the exhaustive covariate screen only when requested

if (isTRUE(cfg[["run_original_covariate_screen"]])) {
  source(file.path(script_dir, "10_exploratory_screen.R"))
  source(file.path(script_dir, "11_covariate_summary.R"))
} else {
  message("original covariate screen is disabled")
}

# run the extended analyses

source(file.path(script_dir, "12_extended_model_cv.R"))
source(file.path(script_dir, "13_extended_inference.R"))
source(file.path(script_dir, "14_mutant_heterogeneity.R"))
source(file.path(script_dir, "15_target_expression_sensitivity.R"))
source(file.path(script_dir, "16_specificity_and_sensitivity.R"))
source(file.path(script_dir, "17_expression_response.R"))
source(file.path(script_dir, "18_extended_lineage_holdout.R"))
source(file.path(script_dir, "19_extended_report.R"))

message("all research analyses completed")
