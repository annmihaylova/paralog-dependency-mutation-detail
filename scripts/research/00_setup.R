suppressPackageStartupMessages(library(data.table))


cfg <- list(
  release = "DepMap Public 26Q3",
  raw_dir = "data/raw",
  out_dir = "results/final_mutation_biomarkers",
  cv_repeats = 5L,
  cv_outer_folds = 5L,
  cv_inner_folds = 5L,
  seed = 20261107L,
  bootstrap_reps = 2000L,
  run_lineage_holdout = TRUE,
  extended_cv_repeats = 3L,
  extended_outer_folds = 5L,
  extended_inner_folds = 5L,
  extended_bootstrap_reps = 2000L,
  extended_permutation_reps = 2000L,
  expression_cutoffs = c(1, 2, 3),
  min_expressing_mutants = 25L,
  run_original_covariate_screen = FALSE
)

cfg[["feature_dir"]] <- file.path(cfg[["out_dir"]], "features")
cfg[["table_dir"]] <- file.path(cfg[["out_dir"]], "tables")
cfg[["figure_dir"]] <- file.path(cfg[["out_dir"]], "figures")
for (p in c(cfg[["out_dir"]], cfg[["feature_dir"]], cfg[["table_dir"]], cfg[["figure_dir"]])) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

# analyze four predefined pairs supported by literature
pairs <- data.table(
  gene_A = c("ARID1A", "SMARCA4", "STAG2", "CREBBP"),
  gene_B = c("ARID1B", "SMARCA2", "STAG1", "EP300")
)
pairs[, pair := paste(gene_A, "->", gene_B)]
pairs[, stem := paste(gene_A, gene_B, sep = "_")]

# locate input files in the raw data directory
find_raw <- function(names, required = TRUE) {
  for (nm in names) {
    path <- file.path(cfg[["raw_dir"]], nm)
    if (file.exists(path)) return(path)
  }
  if (required) {
    stop("Required DepMap file not found: ",
         paste(names, collapse = " OR "))
  }
  NULL
}

# validate column names before loading large datasets
choose_column <- function(headers, options, required = TRUE) {
  hit <- intersect(options, headers)
  if (length(hit) > 0L) return(hit[[1L]])
  if (required) stop("Missing column. Expected: ", paste(options, collapse = ", "),
                     ". Available: ", paste(headers, collapse = ", "))
  NULL
}

# load only the required genes and ensure one record per cell-line model
read_gene_matrix <- function(path, genes) {
  headers <- names(fread(path, nrows = 0L, check.names = FALSE))
  id <- choose_column(headers, c("ModelID", "DepMap_ID", "DepMapID", "Model_ID"), FALSE)
  if (is.null(id)) id <- headers[[1L]]
  default <- choose_column(headers, "IsDefaultEntryForModel", FALSE)
  gene_cols <- vapply(genes, function(g) {
    hits <- headers[headers == g | startsWith(headers, paste0(g, " ("))]
    if (length(hits) == 0L) return(NA_character_)
    if (length(hits) > 1L) stop("Ambiguous gene column for: ", g)
    hits[[1L]]
  }, character(1))
  cols <- unique(c(id, default, unname(gene_cols[!is.na(gene_cols)])))
  d <- fread(path, select = cols, check.names = FALSE, showProgress = FALSE)
  setnames(d, id, "ModelID")
  if (!is.null(default)) {
    d <- d[tolower(trimws(as.character(get(default)))) == "yes"]
    d[, (default) := NULL]
  }
  if (anyDuplicated(d[["ModelID"]])) stop("Duplicate ModelID: ", basename(path))
  for (g in names(gene_cols)) {
    if (!is.na(gene_cols[[g]])) setnames(d, gene_cols[[g]], g)
  }
  d
}

# empty OncotreeLineage values is missing data, not a separate cancer lineage
read_models <- function(cols) {
  m <- fread(find_raw("Model.csv"), select = cols, check.names = FALSE, showProgress = FALSE)
  if (anyDuplicated(m[["ModelID"]])) stop("Duplicate ModelID in Model.csv")
  m[trimws(OncotreeLineage) == "", OncotreeLineage := NA_character_]
  m
}

# generate stratified cross-validation folds shared across competing models
stratified_folds <- function(status, k) {
  folds <- integer(length(status))
  for (s in sort(unique(status))) {
    idx <- which(status == s)
    folds[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  folds
}

# extract regression coefficients with HC3
robust_coefficients <- function(fit) {
  m <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
  data.table(
    term = rownames(m),
    estimate = as.numeric(m[, "Estimate"]),
    robust_SE = as.numeric(m[, "Std. Error"]),
    t_value = as.numeric(m[, "t value"]),
    p_value = as.numeric(m[, "Pr(>|t|)"])
  )
}

# save output tables
save_table <- function(x, name) {
  fwrite(as.data.table(x), file.path(cfg[["table_dir"]], name))
}

# record R version and package versions
writeLines(c(
  paste("DepMap release:", cfg[["release"]]),
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  capture.output(sessionInfo())
), file.path(cfg[["out_dir"]], "session_info.txt"))

message("Research setup: ", cfg[["out_dir"]])
