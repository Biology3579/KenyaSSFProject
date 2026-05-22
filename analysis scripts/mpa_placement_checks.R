# ============================================================
#  MPA PLACEMENT CHECKS
#  Tests whether MPA placement is non-random with respect
#  to human pressure and larval connectivity.
#
#  Tests:
#  1. Kruskal-Wallis: MPA status vs each pressure metric
#     and connectivity
#  2. Pairwise Wilcoxon with Bonferroni correction: which
#     MPA categories differ
#  3. Descriptive summaries by MPA category
#
#  All tests use scaled predictors from final_predictors.
#  Rank-based tests are scale-invariant so results are
#  identical to unscaled versions.
# ============================================================

source(here::here("predictor_preparation.R"))


# ============================================================
#  DATA PREPARATION
# ============================================================

# ── Site-level dataset ────────────────────────────────────────
mpa_placement_data <- final_predictors %>%
  mutate(
    mpa_status = factor(mpa_status,
                        levels  = c("none", "low", "medium"),
                        ordered = FALSE)
  )

cat("Sites:", nrow(mpa_placement_data), "\n")
cat("\nMPA status counts:\n")
print(table(mpa_placement_data$mpa_status))

# ── Filter to model sites only ────────────────────────────────
# lalane (Tanzania S / Mozambique N, none MPA) excluded from
# all models — only 2 transects, below minimum threshold of 3.
# Placement checks use same 54 sites as all model analyses.
mpa_placement_data <- mpa_placement_data %>%
  filter(site != "lalane")

cat("Sites after filtering:", nrow(mpa_placement_data), "\n")
cat("\nMPA status counts:\n")
print(table(mpa_placement_data$mpa_status))

# ============================================================
#  DESCRIPTIVE SUMMARIES BY MPA CATEGORY
# ============================================================

cat("\n--- Descriptive summaries by MPA status ---\n")

vars <- c("log_settlement_grav_sc",
          "log_market_gravity_sc",
          "log_settlement_pop_sc",
          "connectivity_sc")

var_labels <- c("Settlement gravity (z)",
                "Market gravity (z)",
                "Settlement population (z)",
                "Connectivity (z)")

for (i in seq_along(vars)) {
  cat(sprintf("\n%s:\n", var_labels[i]))
  mpa_placement_data %>%
    group_by(mpa_status) %>%
    summarise(
      n      = n(),
      mean   = round(mean(.data[[vars[i]]], na.rm = TRUE), 3),
      median = round(median(.data[[vars[i]]], na.rm = TRUE), 3),
      sd     = round(sd(.data[[vars[i]]], na.rm = TRUE), 3),
      .groups = "drop"
    ) %>%
    print()
}


# ============================================================
#  KRUSKAL-WALLIS TESTS
#  Tests whether predictor distributions differ across
#  MPA categories.
# ============================================================

cat("\n--- Kruskal-Wallis tests ---\n")

kw_results <- list()

for (i in seq_along(vars)) {
  kw <- kruskal.test(
    as.formula(paste(vars[i], "~ mpa_status")),
    data = mpa_placement_data
  )
  kw_results[[vars[i]]] <- kw
  cat(sprintf("\n%s:\n", var_labels[i]))
  cat(sprintf("  H = %.3f, df = %d, p = %.4f\n",
              kw$statistic, kw$parameter, kw$p.value))
}


# ============================================================
#  PAIRWISE WILCOXON TESTS
#  Bonferroni correction applied.
#  Three comparisons per variable:
#    none vs low
#    none vs medium
#    low vs medium
# ============================================================

cat("\n--- Pairwise Wilcoxon tests (Bonferroni correction) ---\n")

for (i in seq_along(vars)) {
  cat(sprintf("\n%s:\n", var_labels[i]))
  pw <- pairwise.wilcox.test(
    x           = mpa_placement_data[[vars[i]]],
    g           = mpa_placement_data$mpa_status,
    p.adjust.method = "bonferroni",
    exact       = FALSE
  )
  print(pw$p.value)
}


# ============================================================
#  RESULTS SUMMARY
# ============================================================

cat("\n--- MPA placement results summary ---\n")

cat("\nPressure metrics — key comparisons:\n")
cat("Settlement gravity:\n")
cat(sprintf("  KW: H = %.3f, p = %.4f\n",
            kw_results[["log_settlement_grav_sc"]]$statistic,
            kw_results[["log_settlement_grav_sc"]]$p.value))

cat("Market gravity:\n")
cat(sprintf("  KW: H = %.3f, p = %.4f\n",
            kw_results[["log_market_gravity_sc"]]$statistic,
            kw_results[["log_market_gravity_sc"]]$p.value))

cat("Settlement population:\n")
cat(sprintf("  KW: H = %.3f, p = %.4f\n",
            kw_results[["log_settlement_pop_sc"]]$statistic,
            kw_results[["log_settlement_pop_sc"]]$p.value))

cat("\nConnectivity:\n")
cat(sprintf("  KW: H = %.3f, p = %.4f\n",
            kw_results[["connectivity_sc"]]$statistic,
            kw_results[["connectivity_sc"]]$p.value))

#Cis
ci_df <- mpa_placement_data %>%
  group_by(mpa_status) %>%
  summarise(
    across(
      c(log_settlement_grav_sc, 
        log_settlement_pop_sc, 
        connectivity_sc),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        lwr  = ~mean(.x, na.rm = TRUE) - 1.96 * (sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x)))),
        upr  = ~mean(.x, na.rm = TRUE) + 1.96 * (sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))))
      )
    ),
    .groups = "drop"
  )

print(ci_df)

# ── End of script ─────────────────────────────────────────────

