# load gene-level data for both members of each paralog pair

models <- read_models(
  c("ModelID", "PatientID", "OncotreeLineage")
)

a_genes <- pairs[["gene_A"]]
b_genes <- pairs[["gene_B"]]
all_genes <- unique(c(a_genes, b_genes))

damaging <- read_gene_matrix(
  find_raw("OmicsSomaticMutationsMatrixDamaging.csv"),
  all_genes
)

effects <- read_gene_matrix(
  find_raw("CRISPRGeneEffect.csv"),
  b_genes
)

expression <- read_gene_matrix(
  find_raw("OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"),
  all_genes
)

cn <- read_gene_matrix(
  find_raw("OmicsCNGeneWGS.csv"),
  all_genes
)

# keep the full pair matrices for cross-target controls
expression_all <- copy(expression)
effects_all <- copy(effects)

# extract variants annotated as likely loss-of-function (LoF)
variant_path <- find_raw("OmicsSomaticMutations.csv")
headers <- names(fread(variant_path, nrows = 0L, check.names = FALSE))
id_col <- choose_column(headers, "ModelID")
gene_col <- choose_column(headers, c("HugoSymbol", "Hugo_Symbol", "GeneSymbol"))
af_col <- choose_column(headers, c("AlleleFraction", "AF", "VariantAlleleFraction", "VAF"))
# molecular consequences are stored in VariantInfo
# variant type describes general classes such as snv deletion and insertion
type_col <- choose_column(headers, c("VariantInfo", "VariantClassification",
                                     "Variant_Classification", "MolecularConsequence"))
lof_col <- choose_column(headers, c("LikelyLoF", "IsLikelyLoF"), FALSE)
def_col <- choose_column(headers, "IsDefaultEntryForModel", FALSE)
selected <- unique(c(id_col, gene_col, af_col, type_col, lof_col, def_col))
variants <- fread(variant_path, select = selected, showProgress = FALSE)
setnames(variants, c(id_col, gene_col, af_col, type_col),
         c("ModelID", "gene", "vaf", "type"))
variants <- variants[gene %in% pairs[["gene_A"]]]
if (!is.null(def_col)) {
  variants <- variants[tolower(trimws(as.character(get(def_col)))) == "yes"]
}
variants[, vaf := suppressWarnings(as.numeric(vaf))]
variants[, type := tolower(as.character(type))]
if (any(!is.na(variants[["vaf"]]) &
        (variants[["vaf"]] < 0 | variants[["vaf"]] > 1))) {
  stop("VAF outside [0, 1]; verify column units in OmicsSomaticMutations.csv")
}
if (!is.null(lof_col)) {
  variants <- variants[tolower(as.character(get(lof_col))) == "true" |
                         tolower(as.character(get(lof_col))) == "t" |
                         tolower(as.character(get(lof_col))) == "yes" |
                         as.character(get(lof_col)) == "1"]
} else {
  warning("No LikelyLoF column; using a conservative variant-type proxy.")
  variants <- variants[grepl("frame.?shift|nonsense|stop.?gain|splice.?site|splice.?donor|splice.?acceptor",
                             type, ignore.case = TRUE)]
}
variants[, truncating := as.integer(grepl("frame.?shift|nonsense|stop.?gain", type))]
variants[, splice := as.integer(grepl("splice.?(site|donor|acceptor)", type))]
if (nrow(variants) > 0L && sum(variants[["truncating"]] + variants[["splice"]]) == 0L) {
  warning("No truncating/splice consequences found; check the '", type_col, "' column.")
}

# summarize likely-LoF allele fractions and mutation consequences
var_summary <- variants[, .(
  n_lof = .N,
  max_vaf = if (all(is.na(vaf))) NA_real_
  else max(vaf, na.rm = TRUE),
  sum_af = if (anyNA(vaf)) NA_real_
  else sum(vaf),
  n_truncating = sum(truncating),
  n_splice = sum(splice),
  has_noncanonical_lof = as.integer(
    any(truncating == 0 & splice == 0, na.rm = TRUE)
  )
), by = .(ModelID, gene)]

# merge source datasets by ModelID
exclusion_log <- list()
make_features <- function(a, b) {
  if (!(a %in% names(damaging)) || !(a %in% names(expression)) ||
      !(a %in% names(cn)) || !(b %in% names(effects)) ||
      !(b %in% names(damaging)) || !(b %in% names(expression)) ||
      !(b %in% names(cn))) {
    stop("Gene not found in a required matrix for ", a, " -> ", b)
  }
  d <- merge(models, damaging[, .(ModelID, status = get(a))], by = "ModelID")
  d <- merge(d, effects[, .(ModelID, effect = get(b))], by = "ModelID")
  d <- merge(d, expression[, .(ModelID, expression_A = get(a))], by = "ModelID")
  d <- merge(d, cn[, .(ModelID, copy_number_A = get(a))], by = "ModelID")
  # add target-paralog features without changing the original cohort
  d <- merge(d, expression[, .(ModelID, expression_B = get(b))], by = "ModelID", all.x = TRUE)
  d <- merge(d, cn[, .(ModelID, copy_number_B = get(b))], by = "ModelID", all.x = TRUE)
  d <- merge(d, damaging[, .(ModelID, status_B = get(b))], by = "ModelID", all.x = TRUE)
  d[, ratio_AB := expression_A - expression_B]
  d <- merge(d, var_summary[gene == a,
                           .(ModelID, max_vaf, sum_af, n_lof,
                             n_truncating, n_splice, has_noncanonical_lof)],
             by = "ModelID", all.x = TRUE)

  # status 0 indicates no damaging variants recorded in the mutation matrix
  for (nm in c("max_vaf", "sum_af", "n_lof", "n_truncating", "n_splice", "has_noncanonical_lof")) {
    set(d, i = which(d[["status"]] == 0 & is.na(d[[nm]])), j = nm, value = 0)
  }
  d <- d[status == 0 | status == 1 | status == 2]
  # identify and track discrepancies between mutation data sources
  conflict_0 <- d[["status"]] == 0 & d[["n_lof"]] > 0
  missing_detail <- d[["status"]] > 0 & is.na(d[["n_lof"]])
  if (any(conflict_0)) warning(a, ": ", sum(conflict_0),
                               " status-0 lines have LikelyLoF variants; excluded.")
  if (any(missing_detail)) warning(a, ": ", sum(missing_detail),
                                   " mutant lines lack variant details; excluded.")
  d <- d[!conflict_0]
  needed <- c("OncotreeLineage", "effect", "expression_A", "copy_number_A",
              "max_vaf", "n_lof", "n_truncating", "n_splice")
  d <- d[complete.cases(d[, ..needed])]
  exclusion_log[[a]] <<- c(conflict_status0 = sum(conflict_0),
                           mutant_without_details = sum(missing_detail))
  numeric_cols <- c("effect", "expression_A", "copy_number_A",
                    "max_vaf", "n_lof", "n_truncating", "n_splice")
  finite_rows <- Reduce(`&`, lapply(numeric_cols, function(nm) is.finite(d[[nm]])))
  d <- d[finite_rows]
  # group cancer lineages with fewer than 15 cell lines into "Other"
  counts <- table(d[["OncotreeLineage"]])
  d[, lineage_group := ifelse(counts[OncotreeLineage] >= 15,
                              OncotreeLineage, "Other")]
  d[, OncotreeLineage := NULL]
  setcolorder(d, c("ModelID", "PatientID", "status", "status_B", "effect",
                   "expression_A", "expression_B", "ratio_AB",
                   "copy_number_A", "copy_number_B", "max_vaf", "sum_af",
                   "n_lof", "n_truncating", "n_splice", "has_noncanonical_lof",
                   "lineage_group"))
  d
}
feature_data <- lapply(seq_len(nrow(pairs)), function(i) {
  make_features(pairs[["gene_A"]][i], pairs[["gene_B"]][i])
})
names(feature_data) <- pairs[["stem"]]
source_label <- "raw_DepMap_26Q3"

# compare mutation categories with the reported allele-fraction rule
status_af_audit <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  d <- feature_data[[pairs[["stem"]][i]]]
  d[!is.na(sum_af), .(
    pair = pairs[["pair"]][i], ModelID, status, sum_af,
    expected_status = fifelse(n_lof == 0, 0L,
                       fifelse(sum_af > 0.95, 2L, 1L)),
    agrees = status == fifelse(n_lof == 0, 0L,
                       fifelse(sum_af > 0.95, 2L, 1L))
  )]
}))
save_table(status_af_audit, "01_status_af_audit.csv")


# validate feature completeness, mutation-status groups, and cohort sizes
needed_cols <- c("ModelID", "status", "effect", "expression_A", "copy_number_A",
                 "max_vaf", "n_lof", "n_truncating", "n_splice", "lineage_group")
cohort_audit <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  stem <- pairs[["stem"]][i]
  d <- feature_data[[stem]]
  absent <- setdiff(needed_cols, names(d))
  if (length(absent)) stop(stem, ": missing ", paste(absent, collapse = ", "))
  if (anyDuplicated(d[["ModelID"]])) stop(stem, ": duplicate ModelID")
  if (anyNA(d[, ..needed_cols])) stop(stem, ": NA values in features")
  for (nm in c("status", "effect", "expression_A", "copy_number_A",
               "max_vaf", "n_lof", "n_truncating", "n_splice")) {
    if (!all(is.finite(d[[nm]]))) stop(stem, ": non-finite ", nm)
  }
  if (!all(d[["status"]] == 0 | d[["status"]] == 1 | d[["status"]] == 2)) {
    stop(stem, ": statuses must be 0/1/2")
  }
  if (any(d[["status"]] == 0 & d[["n_lof"]] > 0)) {
    stop(stem, ": contradictory status and n_lof")
  }
  n_status <- table(factor(d[["status"]], levels = 0:2))
  if (any(n_status < cfg[["cv_outer_folds"]] + 1L)) stop(stem, ": too few lines for CV")
  setorder(d, ModelID)
  fwrite(d, file.path(cfg[["feature_dir"]], paste0(stem, "_features.csv")))
  excl <- exclusion_log[[pairs[["gene_A"]][i]]]
  data.table(pair = pairs[["pair"]][i], n = nrow(d),
             n_0 = as.integer(n_status[1]), n_1 = as.integer(n_status[2]),
             n_2 = as.integer(n_status[3]),
             excluded_status0_with_LoF = if (is.null(excl)) NA_integer_ else excl[[1]],
             excluded_mutant_without_details = if (is.null(excl)) NA_integer_ else excl[[2]],
             n_truncating_lines = sum(d[["n_truncating"]] > 0),
             n_splice_lines = sum(d[["n_splice"]] > 0),
             source = source_label)
}), fill = TRUE)

save_table(cohort_audit, "01_cohort_audit.csv")
message("Cohorts ready: ", paste(cohort_audit[["pair"]], cohort_audit[["n"]], collapse = "; "))

