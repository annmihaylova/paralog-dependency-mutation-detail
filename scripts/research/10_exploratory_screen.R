# explore covariate associations across all four pairs
# check stability across repeated lineage splits

screen_cfg <- list(n_splits = 200L, seed = 20261104L, alpha = 0.05, replication_rule = 0.5)

# distribute lineages across halves with balanced cell counts
split_lineages <- function(lin_counts) {
  lin <- sample(names(lin_counts))
  half <- setNames(integer(length(lin)), lin)
  size <- c(0, 0)
  for (l in lin) {
    h <- which.min(size)
    half[[l]] <- h
    size[h] <- size[h] + lin_counts[[l]]
  }
  half
}

# test main effects or mutant interactions in one half
test_in_half <- function(s, v, type, what) {
  if (what == "main") {
    if (type == "binary" && any(table(factor(s[[v]], levels = 0:1)) < cov_cfg[["min_level_lines"]])) return(NULL)
    if (type == "numeric" && (nrow(s) < 30L || sd(s[[v]]) == 0)) return(NULL)
    m <- hc3_fit(as.formula(paste("effect ~", base_rhs, "+", v)), s)
    if (is.null(m[["table"]])) return(NULL)
    r <- m[["table"]][term == v]
  } else {
    if (!interaction_eligibility(s, v, type)[["simple"]]) return(NULL)
    m <- hc3_fit(as.formula(paste("effect ~ mutant *", v, "+ expression_A + lineage_group")), s)
    if (is.null(m[["table"]])) return(NULL)
    r <- m[["table"]][term == paste0("mutant:", v)]
  }
  if (!nrow(r) || !is.finite(r[["p_value"]])) return(NULL)
  c(estimate = r[["estimate"]], p = r[["p_value"]])
}

set.seed(screen_cfg[["seed"]])
screen_rows <- list()
for (i in seq_len(nrow(pairs))) {
  d <- pair_tables[[i]][!is.na(lineage_original)]
  lin_counts <- table(d[["lineage_original"]])
  splits <- replicate(screen_cfg[["n_splits"]], split_lineages(lin_counts), simplify = FALSE)
  vars <- setdiff(covs_for_pair(d), "primary_disease")
  message("Exploratory screen: ", pairs[["pair"]][i], " (", length(vars), " covariates)")
  for (v in vars) {
    type <- cov_spec[id == v, type]
    full <- prepare_subset(d, v, center = if (type == "numeric") v else character())
    for (what in c("main", "mutant_interaction")) {
      # record the complete-sample estimate for comparison
      all_est <- test_in_half(full, v, type, what)
      rep_ok <- 0L; valid <- 0L
      for (sp in splits) {
        h <- sp[full[["lineage_original"]]]
        res <- lapply(1:2, function(k) {
          s <- copy(full[h == k])
          if (nrow(s) < 30L) return(NULL)
          s[, lineage_group := recollapse(as.character(lineage_group))]
          if (type == "numeric") set(s, j = v, value = s[[v]] - mean(s[[v]]))
          test_in_half(s, v, type, what)
        })
        if (is.null(res[[1]]) || is.null(res[[2]])) next
        valid <- valid + 1L
        if (res[[1]][["p"]] < screen_cfg[["alpha"]] && res[[2]][["p"]] < screen_cfg[["alpha"]] &&
            sign(res[[1]][["estimate"]]) == sign(res[[2]][["estimate"]])) rep_ok <- rep_ok + 1L
      }
      screen_rows[[length(screen_rows) + 1L]] <- data.table(
        pair = pairs[["pair"]][i], covariate = v, tier = cov_spec[id == v, tier], test = what,
        estimate_all = if (is.null(all_est)) NA_real_ else all_est[["estimate"]],
        p_all = if (is.null(all_est)) NA_real_ else all_est[["p"]],
        valid_splits = valid,
        # count invalid splits as nonreplicated
        replication_rate = rep_ok / screen_cfg[["n_splits"]],
        replication_rate_among_valid = if (valid > 0L) rep_ok / valid else NA_real_)
    }
  }
}
screen <- rbindlist(screen_rows)
screen[!is.na(p_all), FDR_BH_all := p.adjust(p_all, method = "BH"), by = test]
screen[, replicated := replication_rate >= screen_cfg[["replication_rule"]]]
setorder(screen, test, -replication_rate, p_all)
save_table(screen, "10_exploratory_screen.csv")

message("Exploratory screen done: ", sum(screen[["replicated"]], na.rm = TRUE),
        " of ", nrow(screen), " tests replicated in >= 50% of lineage splits.")

