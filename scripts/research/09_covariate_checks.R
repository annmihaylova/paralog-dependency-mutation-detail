# test exploratory covariates under the documented analysis plan
# requires setup helpers and prepared feature data

cov_cfg <- list(
  smd_minor = 0.10,          # based on austin 2009
  smd_strong = 0.25,         # common heuristic threshold
  cie_relative = 0.10,       # based on maldonado and greenland 1993
  cie_absolute = 0.02,       # protect estimates close to zero
  corr_drop = 0.80,          # collinearity threshold for the full model
  min_level_lines = 5L,      # collapse rare lineage levels in subsets
  int_bin_mutant = 8L,       # minimum mutant lines per binary covariate level
  int_bin_cell = 5L,         # minimum lines per status by covariate cell
  int_num_mutant = 20L,      # minimum nonmissing mutant lines for numeric covariates
  int_num_status = 15L,      # minimum lines in each mutant category
  primary_disease_min = 10L, # collapse rare primary disease groups
  msi_threshold = 20         # msi threshold from the depmap readme
)

# define complex members excluding the target pair
complex_genes <- list(
  SWI_SNF = c("ARID1A", "ARID1B", "ARID2", "PBRM1", "SMARCA4", "SMARCA2", "SMARCB1",
              "SMARCC1", "SMARCC2", "SMARCD1", "SMARCE1", "ACTL6A", "DPF2", "BRD7",
              "BRD9", "BCL7A", "BCL11B", "SS18"),
  cohesin = c("STAG2", "STAG1", "RAD21", "SMC1A", "SMC3", "PDS5A", "PDS5B", "NIPBL", "WAPL")
)
pair_complex <- c(ARID1A = "SWI_SNF", SMARCA4 = "SWI_SNF", STAG2 = "cohesin", CREBBP = NA)

# define covariate types availability and exploratory priorities
cov_spec <- data.table(
  id = c("msi", "log_damaging_burden", "tp53_damaging", "ploidy", "wgd", "loh_fraction",
         "aneuploidy", "expression_B", "copy_number_B", "complex_comutation", "sex_female",
         "primary_disease", "screen_doubling_time", "cas_activity", "screen_nnmd",
         "growth_3d", "growth_suspension", "organoid", "metastatic", "pediatric", "age"),
  type = c("binary", "numeric", "binary", "numeric", "binary", "numeric",
           "numeric", "numeric", "numeric", "binary", "binary",
           "categorical", "numeric", "numeric", "numeric",
           "binary", "binary", "binary", "binary", "binary", "numeric"),
  tier = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3),
  source = c("OmicsGlobalSignatures", "DamagingMatrix", "DamagingMatrix",
             "OmicsGlobalSignatures", "OmicsGlobalSignatures", "OmicsGlobalSignatures",
             "OmicsGlobalSignatures", "Expression", "CopyNumber", "DamagingMatrix", "Model",
             "Model", "CRISPRConfounders", "CRISPRConfounders", "CRISPRConfounders",
             "Model", "Model", "Model", "Model", "Model", "Model")
)

# identify literature-motivated interaction candidates
literature_motivated_interactions <- rbind(
  data.table(gene_A = "STAG2", covariate = "sex_female"),
  data.table(gene_A = "ARID1A", covariate = c("msi", "log_damaging_burden")),
  data.table(gene_A = pairs[["gene_A"]], covariate = "expression_B"),
  data.table(gene_A = pairs[["gene_A"]], covariate = "loh_fraction")
)

# assemble covariates

# match column names without case sensitivity
pick_col <- function(headers, options) {
  hit <- headers[tolower(headers) %in% tolower(options)]
  if (length(hit)) hit[[1L]] else NULL
}

# read model-level measurements from tabular files
read_model_table <- function(path, wanted) {
  headers <- names(fread(path, nrows = 0L, check.names = FALSE))
  id <- pick_col(headers, c("ModelID", "DepMap_ID"))
  if (is.null(id)) id <- headers[[1L]]
  def <- pick_col(headers, "IsDefaultEntryForModel")
  cols <- unlist(lapply(wanted, function(w) pick_col(headers, w)))
  d <- fread(path, select = unique(c(id, def, cols)), check.names = FALSE, showProgress = FALSE)
  setnames(d, id, "ModelID")
  if (!is.null(def)) {
    d <- d[tolower(trimws(as.character(get(def)))) == "yes"]
    d[, (def) := NULL]
  }
  for (w in wanted) {
    hit <- pick_col(names(d), w)
    if (!is.null(hit) && hit != w) setnames(d, hit, w)
  }
  # reconcile multiple default rows from sequencing profiles
  # retain nonmissing values and log conflicts
  if (anyDuplicated(d[["ModelID"]])) {
    cols <- intersect(wanted, names(d))
    dup_ids <- unique(d[["ModelID"]][duplicated(d[["ModelID"]])])
    conflicts <- rbindlist(lapply(cols, function(cl) {
      x <- d[ModelID %in% dup_ids, .(n_values = uniqueN(na.omit(get(cl)))), by = ModelID]
      x[n_values > 1L, .(file = basename(path), ModelID, column = cl, n_values)]
    }))
    duplicate_log[[basename(path)]] <<- list(n_models = length(dup_ids), conflicts = conflicts)
    message(basename(path), ": ", length(dup_ids), " models have several default rows; ",
            "non-missing values combined, ", nrow(conflicts), " conflicting values (first non-missing kept).")
    d <- d[, lapply(.SD, function(v) { k <- which(!is.na(v)); if (length(k)) v[[k[1L]]] else v[[1L]] }),
           by = ModelID, .SDcols = cols]
  }
  d
}
duplicate_log <- list()

availability <- list()
note_source <- function(source, path) {
  availability[[source]] <<- data.table(
    source = source, file = if (is.null(path)) NA_character_ else basename(path),
    available = !is.null(path))
}

# load demographic and cell-model metadata
models_path <- find_raw("Model.csv")
note_source("Model", models_path)
model_meta <- fread(models_path, select = c(
  "ModelID", "OncotreeLineage", "OncotreePrimaryDisease", "Sex", "PrimaryOrMetastasis",
  "GrowthPattern", "ModelType", "AgeCategory", "Age"), showProgress = FALSE)
blank_to_na <- function(x) { x <- trimws(as.character(x)); x[x == "" | x == "Unknown"] <- NA; x }
for (nm in setdiff(names(model_meta), "ModelID")) set(model_meta, j = nm, value = blank_to_na(model_meta[[nm]]))
covariates <- model_meta[, .(
  ModelID,
  lineage_original = OncotreeLineage,
  primary_disease_raw = OncotreePrimaryDisease,
  sex_female = fifelse(Sex == "Female", 1L, fifelse(Sex == "Male", 0L, NA_integer_)),
  metastatic = fifelse(PrimaryOrMetastasis %in% c("Metastatic", "Recurrent"), 1L,
                       fifelse(PrimaryOrMetastasis == "Primary", 0L, NA_integer_)),
  growth_3d = fifelse(GrowthPattern %in% c("Dome", "Spheroid"), 1L,
                      fifelse(is.na(GrowthPattern), NA_integer_, 0L)),
  growth_suspension = fifelse(GrowthPattern == "Suspension", 1L,
                              fifelse(is.na(GrowthPattern), NA_integer_, 0L)),
  organoid = fifelse(ModelType == "Organoid", 1L,
                     fifelse(is.na(ModelType), NA_integer_, 0L)),
  pediatric = fifelse(AgeCategory == "Pediatric", 1L,
                      fifelse(AgeCategory == "Adult", 0L, NA_integer_)),
  age = suppressWarnings(as.numeric(Age))
)]

# derive mutation burden and complex comutation flags
dmg_path <- find_raw("OmicsSomaticMutationsMatrixDamaging.csv")
note_source("DamagingMatrix", dmg_path)
message("Reading damaging-mutation matrix (whole file, for mutation burden)...")
dmg <- fread(dmg_path, check.names = FALSE, showProgress = FALSE)
if (!"ModelID" %in% names(dmg)) stop("Damaging matrix has no ModelID column.")
if ("IsDefaultEntryForModel" %in% names(dmg)) {
  dmg <- dmg[tolower(trimws(as.character(IsDefaultEntryForModel))) == "yes"]
}
if (anyDuplicated(dmg[["ModelID"]])) stop("Duplicate ModelID in damaging matrix")
gene_cols <- grep(" \\([0-9]+\\)$", names(dmg), value = TRUE)
gene_symbol <- sub(" \\([0-9]+\\)$", "", gene_cols)
dmg_mat <- as.matrix(dmg[, ..gene_cols])
dmg_burden <- data.table(ModelID = dmg[["ModelID"]],
                         log_damaging_burden = log10(1 + rowSums(dmg_mat > 0, na.rm = TRUE)))
col_of <- function(g) { i <- which(gene_symbol == g); if (length(i) == 1L) gene_cols[[i]] else NA_character_ }
if (is.na(col_of("TP53"))) stop("TP53 column not found in damaging matrix")
dmg_burden[, tp53_damaging := as.integer(dmg_mat[, col_of("TP53")] > 0)]
complex_flags <- list()
for (a in names(pair_complex)) {
  cx <- pair_complex[[a]]
  if (is.na(cx)) next
  b <- pairs[gene_A == a, gene_B]
  members <- setdiff(complex_genes[[cx]], c(a, b))
  cols <- na.omit(vapply(members, col_of, character(1)))
  if (length(cols) < length(members)) {
    warning(a, ": complex genes missing in matrix: ",
            paste(setdiff(members, gene_symbol), collapse = ", "))
  }
  complex_flags[[a]] <- data.table(ModelID = dmg[["ModelID"]], gene_A = a,
                                   complex_comutation = as.integer(rowSums(dmg_mat[, cols, drop = FALSE] > 0, na.rm = TRUE) > 0),
                                   complex_genes_used = paste(sub(" \\(.*", "", cols), collapse = ";"))
}
complex_flags <- rbindlist(complex_flags)
rm(dmg, dmg_mat); invisible(gc())
covariates <- merge(covariates, dmg_burden, by = "ModelID", all.x = TRUE)

# reuse target-paralog features already prepared in step 01
has_expression_B <- all(vapply(feature_data, function(x) "expression_B" %in% names(x), logical(1)))
has_copy_number_B <- all(vapply(feature_data, function(x) "copy_number_B" %in% names(x), logical(1)))
expr_path <- find_raw(c("OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv",
                        "OmicsExpressionProteinCodingGenesTPMLogp1.csv"), FALSE)
note_source("Expression", expr_path)
expr_B <- if (!has_expression_B && !is.null(expr_path)) {
  read_gene_matrix(expr_path, pairs[["gene_B"]])
} else NULL

# avoid reading the large copy-number matrix twice
cn_path <- find_raw(c("OmicsCNGeneWGS.csv", "OmicsCNGene.csv", "OmicsAbsoluteCNGene.csv"), FALSE)
note_source("CopyNumber", cn_path)
cn_B <- if (!has_copy_number_B && !is.null(cn_path)) {
  read_gene_matrix(cn_path, pairs[["gene_B"]])
} else NULL

# load optional genomic signatures
sig_path <- find_raw(c("OmicsGlobalSignatures.csv", "OmicsSignatures.csv"), FALSE)
note_source("OmicsGlobalSignatures", sig_path)
if (!is.null(sig_path)) {
  sig <- read_model_table(sig_path, c("MSIScore", "Ploidy", "CIN", "WGD", "LoHFraction", "Aneuploidy"))
  sig_out <- data.table(ModelID = sig[["ModelID"]])
  if ("MSIScore" %in% names(sig)) sig_out[, msi := as.integer(as.numeric(sig[["MSIScore"]]) >= cov_cfg[["msi_threshold"]])]
  if ("Ploidy" %in% names(sig)) sig_out[, ploidy := as.numeric(sig[["Ploidy"]])]
  if ("WGD" %in% names(sig)) sig_out[, wgd := as.integer(as.numeric(sig[["WGD"]]) > 0)]
  if ("LoHFraction" %in% names(sig)) sig_out[, loh_fraction := as.numeric(sig[["LoHFraction"]])]
  if ("Aneuploidy" %in% names(sig)) sig_out[, aneuploidy := as.numeric(sig[["Aneuploidy"]])]
  covariates <- merge(covariates, sig_out, by = "ModelID", all.x = TRUE)
}

# load optional CRISPR screen quality measures
conf_path <- find_raw("CRISPRConfounders.csv", FALSE)
note_source("CRISPRConfounders", conf_path)
if (!is.null(conf_path)) {
  conf <- read_model_table(conf_path, c("ScreenDoublingTime", "CasActivity", "ScreenNNMD"))
  conf_out <- data.table(ModelID = conf[["ModelID"]])
  if ("ScreenDoublingTime" %in% names(conf)) conf_out[, screen_doubling_time := as.numeric(conf[["ScreenDoublingTime"]])]
  if ("CasActivity" %in% names(conf)) conf_out[, cas_activity := as.numeric(conf[["CasActivity"]])]
  if ("ScreenNNMD" %in% names(conf)) conf_out[, screen_nnmd := as.numeric(conf[["ScreenNNMD"]])]
  covariates <- merge(covariates, conf_out, by = "ModelID", all.x = TRUE)
}

availability_table <- rbindlist(availability)
# log duplicated default entries and conflicting values
availability_table[, models_with_several_default_rows := vapply(file, function(f) {
  if (is.na(f) || is.null(duplicate_log[[f]])) 0L else as.integer(duplicate_log[[f]][["n_models"]])
}, integer(1))]
save_table(rbindlist(c(list(data.table(file = character(), ModelID = character(), column = character(), n_values = integer())),
                       lapply(duplicate_log, `[[`, "conflicts")), fill = TRUE),
           "09_duplicate_row_conflicts.csv")
available_covariates <- c(
  names(covariates), "complex_comutation", "primary_disease",
  if (has_expression_B || !is.null(expr_B)) "expression_B",
  if (has_copy_number_B || !is.null(cn_B)) "copy_number_B"
)
cov_spec[, available := id %in% available_covariates]
save_table(availability_table, "09_covariate_sources.csv")
save_table(cov_spec, "09_covariate_availability.csv")
if (any(!cov_spec[["available"]])) {
  message("Skipped (source file absent): ", paste(cov_spec[available == FALSE, id], collapse = ", "))
}

# merge pair features with supplementary covariates
pair_data <- function(i) {
  a <- pairs[["gene_A"]][i]; b <- pairs[["gene_B"]][i]
  d <- merge(copy(feature_data[[pairs[["stem"]][i]]]), covariates, by = "ModelID", all.x = TRUE)
  if (!"expression_B" %in% names(d)) {
    if (is.null(expr_B)) stop("Target expression is unavailable for ", a, " -> ", b)
    d <- merge(d, expr_B[, .(ModelID, expression_B = get(b))], by = "ModelID", all.x = TRUE)
  }
  if (!"copy_number_B" %in% names(d) && !is.null(cn_B)) {
    d <- merge(d, cn_B[, .(ModelID, copy_number_B = get(b))], by = "ModelID", all.x = TRUE)
  }
  if (!is.na(pair_complex[[a]])) {
    d <- merge(d, complex_flags[gene_A == a, .(ModelID, complex_comutation)], by = "ModelID", all.x = TRUE)
  }
  # group rare primary cancer subtypes
  pd_counts <- table(d[["primary_disease_raw"]])
  d[, primary_disease := ifelse(is.na(primary_disease_raw), NA_character_,
                                ifelse(pd_counts[primary_disease_raw] >= cov_cfg[["primary_disease_min"]],
                                       primary_disease_raw, "Other"))]
  d[, status_factor := factor(status, levels = 0:2)]
  d[, mutant := as.integer(status > 0)]
  if (nrow(d) != nrow(feature_data[[pairs[["stem"]][i]]])) stop("Row count changed after merge: ", a)
  d
}
pair_tables <- lapply(seq_len(nrow(pairs)), pair_data)
names(pair_tables) <- pairs[["stem"]]

# identify available covariates for this pair
covs_for_pair <- function(d) {
  ids <- cov_spec[available == TRUE, id]
  ids[ids %in% names(d) & vapply(ids, function(x) sum(!is.na(d[[x]])) > 0, logical(1))]
}

save_table(rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  d <- pair_tables[[i]]
  cols <- c("ModelID", "status", "effect", covs_for_pair(d))
  out <- d[, ..cols]; out[, pair := pairs[["pair"]][i]]; out
}), fill = TRUE), "09_covariates_per_model.csv")

# shared model functions

# compute hc3 estimates with leverage safeguards
hc3_fit <- function(formula, data) {
  fit <- lm(formula, data = data)
  b <- coef(fit)
  ok <- !is.na(b)
  # mark leverage-one tests as unavailable
  if (max(hatvalues(fit)) > 0.99) return(list(fit = fit, table = NULL, V = NULL))
  V <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)
  if (is.null(V) || any(!is.finite(diag(V)))) return(list(fit = fit, table = NULL, V = NULL))
  b <- b[ok]
  se <- sqrt(diag(V))[names(b)]
  df <- df.residual(fit)
  tab <- data.table(term = names(b), estimate = unname(b), robust_SE = unname(se),
                    p_value = 2 * pt(abs(unname(b / se)), df, lower.tail = FALSE))
  tab[, `:=`(CI_low = estimate - qt(0.975, df) * robust_SE,
             CI_high = estimate + qt(0.975, df) * robust_SE)]
  list(fit = fit, table = tab, V = V)
}

# conduct joint hc3 wald tests
wald_hc3 <- function(m, terms) {
  if (is.null(m[["V"]])) return(c(stat = NA_real_, df = NA_real_, p = NA_real_))
  terms <- intersect(terms, rownames(m[["V"]]))
  if (!length(terms)) return(c(stat = NA_real_, df = NA_real_, p = NA_real_))
  b <- coef(m[["fit"]])[terms]
  W <- tryCatch(as.numeric(t(b) %*% solve(m[["V"]][terms, terms, drop = FALSE]) %*% b),
                error = function(e) NA_real_)
  c(stat = W, df = length(terms), p = pchisq(W, length(terms), lower.tail = FALSE))
}

# regroup rare lineages within each subset
recollapse <- function(x, min_n = cov_cfg[["min_level_lines"]]) {
  tab <- table(x)
  out <- ifelse(tab[x] >= min_n, x, "Other")
  if (sum(out == "Other") > 0L && sum(out == "Other") < min_n) {
    # merge very small other groups with the largest category
    out[out == "Other"] <- names(which.max(table(out[out != "Other"])))
  }
  factor(out)
}

# build complete-case subsets and center numeric covariates
prepare_subset <- function(d, cols, center = character()) {
  keep <- complete.cases(d[, c("effect", "status", "expression_A", "lineage_group", cols), with = FALSE])
  s <- copy(d[keep])
  s[, lineage_group := recollapse(lineage_group)]
  if ("primary_disease" %in% cols) s[, primary_disease := recollapse(primary_disease)]
  for (nm in center) set(s, j = nm, value = s[[nm]] - mean(s[[nm]]))
  s
}

status_terms <- c("status_factor1", "status_factor2")

# step 1 covariate balance

smd <- function(x1, x0, type) {
  x1 <- x1[!is.na(x1)]; x0 <- x0[!is.na(x0)]
  if (length(x1) < 2L || length(x0) < 2L) return(NA_real_)
  if (type == "binary") {
    p1 <- mean(x1); p0 <- mean(x0); s <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
    if (s == 0) return(NA_real_)
    return((p1 - p0) / s)
  }
  s <- sqrt((var(x1) + var(x0)) / 2)
  if (!is.finite(s) || s == 0) return(NA_real_)
  (mean(x1) - mean(x0)) / s
}

cramers_v <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  tab <- table(x[ok], g[ok])
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (min(dim(tab)) < 2L) return(NA_real_)
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE)[["statistic"]])
  as.numeric(sqrt(chi / (sum(tab) * (min(dim(tab)) - 1L))))
}

set.seed(cfg[["seed"]])
balance <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  d <- pair_tables[[i]]
  rbindlist(lapply(covs_for_pair(d), function(v) {
    type <- cov_spec[id == v, type]
    x <- d[[v]]; g <- d[["status"]]
    row <- data.table(pair = pairs[["pair"]][i], covariate = v, type = type,
                      tier = cov_spec[id == v, tier],
                      n_0 = sum(!is.na(x) & g == 0), n_1 = sum(!is.na(x) & g == 1),
                      n_2 = sum(!is.na(x) & g == 2))
    if (type == "categorical") {
      row[, `:=`(cramers_V = cramers_v(x, g),
                 p_value = tryCatch(suppressWarnings(chisq.test(table(x, g), simulate.p.value = TRUE, B = 5000)[["p.value"]]),
                                    error = function(e) NA_real_),
                 test = "chi-square (Monte Carlo p)")]
      return(row)
    }
    row[, `:=`(value_0 = mean(x[g == 0], na.rm = TRUE), value_1 = mean(x[g == 1], na.rm = TRUE),
               value_2 = mean(x[g == 2], na.rm = TRUE),
               SMD_1_vs_0 = smd(x[g == 1], x[g == 0], type),
               SMD_2_vs_0 = smd(x[g == 2], x[g == 0], type),
               SMD_2_vs_1 = smd(x[g == 2], x[g == 1], type))]
    if (type == "binary") {
      tab <- table(factor(x, levels = 0:1), g)
      row[, `:=`(p_value = if (all(rowSums(tab) > 0)) fisher.test(tab)[["p.value"]] else NA_real_,
                 test = "Fisher exact")]
    } else {
      row[, `:=`(p_value = kruskal.test(x ~ factor(g))[["p.value"]], test = "Kruskal-Wallis")]
    }
    row
  }), fill = TRUE)
}), fill = TRUE)
balance[, max_abs_SMD := pmax(abs(SMD_1_vs_0), abs(SMD_2_vs_0), abs(SMD_2_vs_1), na.rm = TRUE)]
balance[, imbalance := fcase(is.na(max_abs_SMD), NA_character_,
                             max_abs_SMD >= cov_cfg[["smd_strong"]], "strong (>= 0.25)",
                             max_abs_SMD >= cov_cfg[["smd_minor"]], "noticeable (>= 0.10)",
                             default = "balanced (< 0.10)")]
save_table(balance, "09_covariate_balance.csv")

# steps 2 and 3 coefficient sensitivity and main effects

base_rhs <- "status_factor + expression_A + lineage_group"

cie_rows <- list(); main_rows <- list()
for (i in seq_len(nrow(pairs))) {
  d <- pair_tables[[i]]
  for (v in covs_for_pair(d)) {
    type <- cov_spec[id == v, type]
    s <- prepare_subset(d, v)
    if (type == "binary" && length(unique(s[[v]])) < 2L) next
    m0 <- hc3_fit(as.formula(paste("effect ~", base_rhs)), s)
    m1 <- if (v == "primary_disease") {
      hc3_fit(effect ~ status_factor + expression_A + primary_disease, s)
    } else {
      hc3_fit(as.formula(paste("effect ~", base_rhs, "+", v)), s)
    }
    for (tm in status_terms) {
      r0 <- if (is.null(m0[["table"]])) NULL else m0[["table"]][term == tm]
      r1 <- if (is.null(m1[["table"]])) NULL else m1[["table"]][term == tm]
      if (is.null(r0) || is.null(r1) || !nrow(r0) || !nrow(r1)) next
      cie_rows[[length(cie_rows) + 1L]] <- data.table(
        pair = pairs[["pair"]][i], covariate = v, tier = cov_spec[id == v, tier],
        model = if (v == "primary_disease") "lineage replaced by primary disease" else "base + covariate",
        coefficient = tm, n = nrow(s), n_mutant = sum(s[["mutant"]]),
        beta_base = r0[["estimate"]], CI_low_base = r0[["CI_low"]], CI_high_base = r0[["CI_high"]],
        p_base = r0[["p_value"]],
        beta_adjusted = r1[["estimate"]], CI_low_adjusted = r1[["CI_low"]], CI_high_adjusted = r1[["CI_high"]],
        p_adjusted = r1[["p_value"]])
    }
    if (v != "primary_disease" && !is.null(m1[["table"]])) {
      rx <- m1[["table"]][term == v]
      if (nrow(rx)) {
        main_rows[[length(main_rows) + 1L]] <- data.table(
          pair = pairs[["pair"]][i], covariate = v, tier = cov_spec[id == v, tier], type = type,
          n = nrow(s), estimate = rx[["estimate"]], CI_low = rx[["CI_low"]], CI_high = rx[["CI_high"]],
          # scale numeric effects per one standard deviation
          estimate_per_SD = if (type == "numeric") rx[["estimate"]] * sd(s[[v]]) else NA_real_,
          p_value = rx[["p_value"]])
      } else {
        main_rows[[length(main_rows) + 1L]] <- data.table(
          pair = pairs[["pair"]][i], covariate = v, tier = cov_spec[id == v, tier], type = type,
          n = nrow(s), note = "aliased with lineage (no variation within lineage groups)")
      }
    }
  }
}
cie <- rbindlist(cie_rows)
cie[, `:=`(delta = beta_adjusted - beta_base,
           pct_change = 100 * (beta_adjusted - beta_base) / abs(beta_base))]
cie[, confounder_flag := abs(delta) / abs(beta_base) > cov_cfg[["cie_relative"]] &
                         abs(delta) > cov_cfg[["cie_absolute"]]]
save_table(cie, "09_change_in_estimate.csv")

main_effects <- rbindlist(main_rows, fill = TRUE)
main_effects[!is.na(p_value), FDR_BH := p.adjust(p_value, method = "BH")]
save_table(main_effects, "09_covariate_main_effects.csv")

# step 2b correlations and full model

numeric_like <- cov_spec[type != "categorical", id]
corr_rows <- list(); full_rows <- list(); dropped_rows <- list()
drop_priority <- list(c("msi", "log_damaging_burden"), c("loh_fraction", "ploidy", "wgd"),
                      c("growth_3d", "organoid"))
for (i in seq_len(nrow(pairs))) {
  d <- pair_tables[[i]]
  vars <- intersect(covs_for_pair(d), numeric_like)
  if (length(vars) >= 2L) {
    cm <- suppressWarnings(cor(as.matrix(d[, ..vars]), use = "pairwise.complete.obs", method = "spearman"))
    ut <- which(upper.tri(cm), arr.ind = TRUE)
    corr_rows[[i]] <- data.table(pair = pairs[["pair"]][i], var_1 = rownames(cm)[ut[, 1]],
                                 var_2 = colnames(cm)[ut[, 2]], spearman_rho = cm[ut])
  }
  # omit collinear lower-priority predictors in the full model
  full_vars <- intersect(covs_for_pair(d), cov_spec[tier <= 2 & id != "primary_disease", id])
  for (grp in drop_priority) {
    present <- intersect(grp, full_vars)
    if (length(present) < 2L) next
    for (j in 2:length(present)) {
      keep_v <- present[1]
      r <- suppressWarnings(cor(d[[keep_v]], d[[present[j]]], use = "pairwise.complete.obs", method = "spearman"))
      if (is.finite(r) && abs(r) > cov_cfg[["corr_drop"]]) {
        full_vars <- setdiff(full_vars, present[j])
        dropped_rows[[length(dropped_rows) + 1L]] <- data.table(
          pair = pairs[["pair"]][i], dropped = present[j], kept = keep_v, spearman_rho = r)
      }
    }
  }
  # discard predictors with no within-cohort variation
  full_vars <- full_vars[vapply(full_vars, function(v) length(unique(na.omit(d[[v]]))) > 1L, logical(1))]
  s <- prepare_subset(d, full_vars)
  full_vars <- full_vars[vapply(full_vars, function(v) length(unique(s[[v]])) > 1L, logical(1))]
  m0 <- hc3_fit(as.formula(paste("effect ~", base_rhs)), s)
  m1 <- hc3_fit(as.formula(paste("effect ~", base_rhs, if (length(full_vars)) paste("+", paste(full_vars, collapse = " + ")) else "")), s)
  for (tm in status_terms) {
    r0 <- m0[["table"]][term == tm]; r1 <- m1[["table"]][term == tm]
    if (!nrow(r0) || !nrow(r1)) next
    full_rows[[length(full_rows) + 1L]] <- data.table(
      pair = pairs[["pair"]][i], coefficient = tm, n = nrow(s), n_mutant = sum(s[["mutant"]]),
      covariates = paste(full_vars, collapse = ", "),
      beta_base = r0[["estimate"]], CI_low_base = r0[["CI_low"]], CI_high_base = r0[["CI_high"]],
      beta_full = r1[["estimate"]], CI_low_full = r1[["CI_low"]], CI_high_full = r1[["CI_high"]],
      p_full = r1[["p_value"]])
  }
}
full_model <- rbindlist(full_rows)
full_model[, `:=`(delta = beta_full - beta_base,
                  pct_change = 100 * (beta_full - beta_base) / abs(beta_base))]
save_table(rbindlist(corr_rows), "09_covariate_correlations.csv")
save_table(full_model, "09_full_model.csv")
save_table(if (length(dropped_rows)) rbindlist(dropped_rows) else
             data.table(pair = character(), dropped = character(), kept = character(), spearman_rho = numeric()),
           "09_full_model_dropped_collinear.csv")

# step 4 mutation-status interactions

# check eligibility before fitting interaction models
interaction_eligibility <- function(s, v, type) {
  mut <- s[["mutant"]] == 1
  if (type == "binary") {
    m_levels <- table(factor(s[[v]][mut], levels = 0:1))
    simple_ok <- all(m_levels >= cov_cfg[["int_bin_mutant"]])
    cells <- table(factor(s[["status"]][mut], levels = 1:2), factor(s[[v]][mut], levels = 0:1))
    detailed_ok <- simple_ok && all(cells >= cov_cfg[["int_bin_cell"]])
    reason <- if (!simple_ok) sprintf("mutant lines per level of X: %s (need >= %d)",
                                      paste(m_levels, collapse = "/"), cov_cfg[["int_bin_mutant"]]) else ""
  } else {
    n_mut <- sum(mut)
    simple_ok <- n_mut >= cov_cfg[["int_num_mutant"]]
    n_st <- table(factor(s[["status"]][mut], levels = 1:2))
    detailed_ok <- simple_ok && all(n_st >= cov_cfg[["int_num_status"]])
    reason <- if (!simple_ok) sprintf("mutant lines with X: %d (need >= %d)", n_mut, cov_cfg[["int_num_mutant"]]) else ""
  }
  list(simple = simple_ok, detailed = detailed_ok, reason = reason)
}

# compare binary and categorical interaction specifications
interaction_test <- function(s, v, type, form) {
  if (form == "mutant") {
    m <- hc3_fit(as.formula(paste("effect ~ mutant *", v, "+ expression_A + lineage_group")), s)
    int_terms <- paste0("mutant:", v)
  } else {
    m <- hc3_fit(as.formula(paste("effect ~ status_factor *", v, "+ expression_A + lineage_group")), s)
    int_terms <- paste0(status_terms, ":", v)
  }
  w <- wald_hc3(m, int_terms)
  out <- data.table(wald_chi2 = w[["stat"]], wald_df = w[["df"]], p_value = w[["p"]])
  if (form == "mutant" && !is.null(m[["table"]])) {
    ri <- m[["table"]][term == int_terms]
    rm_ <- m[["table"]][term == "mutant"]
    if (nrow(ri) && nrow(rm_)) {
      out[, `:=`(interaction_estimate = ri[["estimate"]], interaction_CI_low = ri[["CI_low"]],
                 interaction_CI_high = ri[["CI_high"]],
                 # evaluate mutation association at the covariate reference value
                 mutant_effect_at_ref = rm_[["estimate"]],
                 mutant_effect_at_X1 = if (type == "binary") rm_[["estimate"]] + ri[["estimate"]] else NA_real_)]
    }
  }
  out
}

int_rows <- list()
for (i in seq_len(nrow(pairs))) {
  d <- pair_tables[[i]]
  a <- pairs[["gene_A"]][i]
  for (v in setdiff(covs_for_pair(d), "primary_disease")) {
    type <- cov_spec[id == v, type]
    s <- prepare_subset(d, v, center = if (type == "numeric") v else character())
    pre <- nrow(literature_motivated_interactions[gene_A == a & covariate == v]) > 0L
    el <- interaction_eligibility(s, v, type)
    base_row <- data.table(pair = pairs[["pair"]][i], covariate = v, tier = cov_spec[id == v, tier],
                           type = type, literature_motivated = pre, n = nrow(s), n_mutant = sum(s[["mutant"]]))
    for (form in c("mutant", "status_0_1_2")) {
      ok <- if (form == "mutant") el[["simple"]] else el[["detailed"]]
      row <- copy(base_row)[, form := form]
      if (!ok) {
        row[, tested := FALSE]
        row[, reason_not_tested := if (form == "mutant") el[["reason"]] else
          if (!el[["simple"]]) el[["reason"]] else "too few lines in status 1 or 2 for this X"]
      } else {
        row <- cbind(row[, tested := TRUE], interaction_test(s, v, type, form))
      }
      int_rows[[length(int_rows) + 1L]] <- row
    }
  }
}
interactions <- rbindlist(int_rows, fill = TRUE)
interactions[, family := fifelse(literature_motivated, "literature_motivated", "exploratory")]
interactions[tested == TRUE & !is.na(p_value), FDR_BH := p.adjust(p_value, method = "BH"), by = .(family, form)]
setcolorder(interactions, c("pair", "covariate", "form", "family", "tested", "n", "n_mutant"))
save_table(interactions, "09_interactions.csv")

# figures

pal_div <- c(low = "#2a78d6", mid = "#f0efec", high = "#e34948")
ink <- "#2b2b29"

# figure 10 displays covariate imbalance
bal_long <- melt(balance[type != "categorical"],
                 id.vars = c("pair", "covariate"),
                 measure.vars = c("SMD_1_vs_0", "SMD_2_vs_0", "SMD_2_vs_1"),
                 variable.name = "comparison", value.name = "SMD")
bal_long[, comparison := factor(comparison, levels = c("SMD_1_vs_0", "SMD_2_vs_0", "SMD_2_vs_1"),
                                labels = c("Status 1 vs 0", "Status 2 vs 0", "Status 2 vs 1"))]
bal_long[, covariate := factor(covariate, levels = rev(intersect(cov_spec[["id"]], unique(covariate))))]
bal_long[, pair := factor(pair, levels = pairs[["pair"]])]
bal_long[, label := ifelse(!is.na(SMD) & abs(SMD) >= cov_cfg[["smd_strong"]], sprintf("%.2f", SMD), "")]
p10 <- ggplot2::ggplot(bal_long, ggplot2::aes(x = comparison, y = covariate, fill = SMD)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 2.7, colour = ink) +
  ggplot2::scale_fill_gradient2(low = pal_div[["low"]], mid = pal_div[["mid"]], high = pal_div[["high"]],
                                midpoint = 0, limits = c(-1, 1), oob = scales::squish,
                                na.value = "grey92", name = "SMD") +
  ggplot2::facet_wrap(~pair, nrow = 1) +
  ggplot2::labs(title = "Covariate balance across mutation-status groups",
                subtitle = "Standardized mean difference; labels shown where |SMD| >= 0.25; blank row = not applicable to this pair",
                x = NULL, y = NULL) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1), panel.grid = ggplot2::element_blank())
ggplot2::ggsave(file.path(cfg[["figure_dir"]], "Figure_10_covariate_balance.png"), p10,
                width = 13, height = 6.5, dpi = 220, bg = "white")

# figure 11 shows changes in the status-two coefficient
cie_plot <- cie[coefficient == "status_factor2"]
cie_plot[, covariate := factor(covariate, levels = rev(intersect(cov_spec[["id"]], unique(covariate))))]
cie_long <- rbind(
  cie_plot[, .(pair, covariate, which = "Base model (same lines)", beta = beta_base, lo = CI_low_base, hi = CI_high_base)],
  cie_plot[, .(pair, covariate, which = "Adjusted for covariate", beta = beta_adjusted, lo = CI_low_adjusted, hi = CI_high_adjusted)])
cie_long[, pair := factor(pair, levels = pairs[["pair"]])]
cie_long[, which := factor(which, levels = c("Base model (same lines)", "Adjusted for covariate"))]
p11 <- ggplot2::ggplot(cie_long, ggplot2::aes(x = beta, y = covariate, colour = which)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.5,
                          position = ggplot2::position_dodge(width = 0.6)) +
  ggplot2::geom_point(size = 1.9, position = ggplot2::position_dodge(width = 0.6)) +
  ggplot2::scale_colour_manual(values = c("#2a78d6", "#eb6834"), name = NULL) +
  ggplot2::facet_wrap(~pair, nrow = 1, scales = "free_x") +
  ggplot2::labs(title = "Status-2 coefficient before and after adjustment",
                subtitle = "HC3 95% CI; both models fitted on the same cell lines",
                x = "Coefficient for status 2 vs 0 (gene effect units)", y = NULL) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(legend.position = "top")
ggplot2::ggsave(file.path(cfg[["figure_dir"]], "Figure_11_change_in_estimate.png"), p11,
                width = 13, height = 6.5, dpi = 220, bg = "white")

message("Covariate checks done: balance, change-in-estimate, full model, interactions.")
