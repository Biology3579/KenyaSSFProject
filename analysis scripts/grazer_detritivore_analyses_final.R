# ============================================================
#  DRIVERS OF GRAZER-DETRITIVORE BIOMASS
#  Functional Group Analysis: Grazers/Detritivores
# ============================================================
#
#  Scientific questions:
#
#  Q1 — Human pressure
#       Does human pressure explain variation in grazer-
#       detritivore biomass beyond local ecological context,
#       and which spatial metric best captures SSF
#       exploitation intensity for this functional group?
#
#  Q2 — Larval connectivity
#       Does larval connectivity explain additional variation
#       beyond the best Q1 model?
#       Connectivity × pressure interaction not evaluated —
#       pressure not supported in Q1.
#
#  Q3 — Formal protection
#       Does MPA status explain additional variation beyond
#       the best Q2 model?
#       Q3 is a separate question — does not update best model.
#       MPA tested last — non-randomly placed with respect
#       to pressure and connectivity (see
#       mpa_placement_checks.R).
#
#  Rationale for sequence:
#       Identical to total biomass, browsers, corallivores —
#       pressure first, connectivity second (extends Warmuth
#       et al. 2024), MPA last as governance response
#       downstream of both.
#
#  Baseline model (fixed a priori, never tested):
#       log(biomass) ~ rugosity_sc + log_chla_sc
#
#  Key difference from browsers/corallivores:
#       Grazer-detritivore biomass is zero-free at site
#       level — lm() on log_mean_biomass used throughout
#       (confirmed in family selection below).
#
#  Sensitivity analyses:
#       (a) Alternative pressure metrics — only metrics
#           within DAICc < 2 of best Q1 model evaluated
#       (b) Transect-level replication (lmer)
# ============================================================

source(here::here("predictor_preparation.R"))


# ============================================================
#  RESPONSE VARIABLE — DATA AGGREGATION AND CHECKS
# ============================================================

# ── Transect-level aggregation ────────────────────────────────
grazer_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_grazer_biomass = sum(
      ifelse(trophic_group %in% c("grazers", "detritivores",
                                  "grazer-detritivores"),
             tot_wt_g, 0),
      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(site = as.factor(site))

cat("Grazer transects:", nrow(grazer_transects), "\n")
cat("Sites:",            n_distinct(grazer_transects$site), "\n")

# ── Transect-level dataset ────────────────────────────────────
grazer_transect_data <- grazer_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(grazer_transect_data$transect_grazer_biomass == 0),
    "/", nrow(grazer_transect_data),
    sprintf("(%.3f)\n",
            mean(grazer_transect_data$transect_grazer_biomass == 0)))

# ── Site-level dataset ────────────────────────────────────────
grazer_model_data <- grazer_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    mean_biomass           = mean(transect_grazer_biomass,
                                  na.rm = TRUE),
    n_transects            = n(),
    ecoregion              = first(ecoregion),
    rugosity_sc            = first(rugosity_sc),
    log_chla_sc            = first(log_chla_sc),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_market_gravity_sc  = first(log_market_gravity_sc),
    log_settlement_pop_sc  = first(log_settlement_pop_sc),
    connectivity_sc        = first(connectivity_sc),
    mpa_status             = first(mpa_status),
    .groups = "drop"
  ) %>%
  mutate(
    site             = as.factor(site),
    ecoregion        = as.factor(ecoregion),
    mpa_status       = factor(mpa_status,
                              levels  = c("none", "low", "medium"),
                              ordered = FALSE),
    log_mean_biomass = log(mean_biomass)
  )

cat("\nGrazer model data:", nrow(grazer_model_data), "sites\n")

# ── Data checks ───────────────────────────────────────────────
grazer_model_data %>%
  dplyr::select(site, rugosity_sc, log_chla_sc,
                log_settlement_grav_sc, connectivity_sc,
                mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  { if (nrow(.) > 0) { warning("NAs in predictors:"); print(.) }
    else cat("NA check passed.\n") }

cat("\nZeros in mean_biomass:",
    sum(grazer_model_data$mean_biomass == 0), "\n")
cat("Site-level zero proportion:",
    round(mean(grazer_model_data$mean_biomass == 0), 3), "\n")
cat("-Inf in log_mean_biomass:",
    sum(is.infinite(grazer_model_data$log_mean_biomass)), "\n")
cat("\nResponse summary:\n")
print(summary(grazer_model_data$mean_biomass))
cat("\nMPA status counts:\n")
print(table(grazer_model_data$mpa_status))


# ============================================================
#  MODEL FAMILY SELECTION
#
#  Grazer-detritivore biomass is zero-free at site level
#  (confirmed above) — log transformation valid without
#  offset. Tweedie convergence failure on baseline model
#  consistent with zero-free response. Gaussian log
#  (lm on log_mean_biomass) is the correct choice.
# ============================================================

grazer_lm_baseline <- lm(
  log_mean_biomass ~ rugosity_sc + log_chla_sc,
  data = grazer_model_data)

par(mfrow = c(2, 2))
plot(grazer_lm_baseline, main = "Gaussian log — baseline")
par(mfrow = c(1, 1))

# Family selection:
#   Gaussian log: SELECTED — lm() on log_mean_biomass.
#     Zero-free response — log transformation valid.
#     Diagnostics (baseline model):
#       Residuals vs Fitted: broadly flat with minor
#         curve — acceptable. Sites 6, 21, 48 flagged.
#       Q-Q: upper tail deviation at sites 6 and 48,
#         lower tail at site 21 — consistent with
#         Shapiro-Wilk p = 0.003. Moderate departure
#         driven by extreme biomass sites.
#       Scale-Location: broadly flat — minor upward
#         trend at high fitted values.
#       Residuals vs Leverage: site 6 has high leverage
#         (leverage ~0.18) and large standardised
#         residual (~4) — approaches Cook's distance
#         threshold. Monitor throughout.
#
#   Tweedie: REJECTED — convergence failure on baseline.
#   Proceed: lm() on log_mean_biomass throughout.


# ============================================================
#  RANDOM EFFECT STRUCTURE
#
#  Tested on baseline model — same rationale as all
#  previous analyses. ML (REML = FALSE) for AICc comparison.
# ============================================================

g_re_null <- lm(
  log_mean_biomass ~ rugosity_sc + log_chla_sc,
  data = grazer_model_data)

g_re_ecoregion <- lmer(
  log_mean_biomass ~ rugosity_sc + log_chla_sc +
    (1 | ecoregion),
  data = grazer_model_data,
  REML = FALSE)

cat("\n--- Random effect structure: ecoregion ---\n")
print(make_aicc_df(list(
  "No RE"           = g_re_null,
  "(1 | ecoregion)" = g_re_ecoregion
)))

# Ecoregion RE has stronger AICc support here than for
# total biomass (DAICc = 2.93, weight = 0.187 for No RE)
# but still not pursued:
# (1) Only 4 ecoregions — insufficient group-level
#     replication (Gelman & Hill 2007).
# (2) Severely uneven group sizes (n = 2, 8, 9, 35).
# Between-ecoregion variation more prominent for this
# group than others — acknowledged as a limitation.
# Stronger ecoregion signal consistent with significant
# residual spatial autocorrelation (see below).
# All models fitted as lm() throughout.


# ============================================================
#  Q1 — HUMAN PRESSURE
#
#  Scientific question:
#  Does human pressure explain variation in grazer-
#  detritivore biomass beyond local ecological context,
#  and which spatial metric best captures SSF exploitation?
#
#  Approach: AICc comparison of baseline vs baseline + each
#  pressure metric. Best metric = highest AICc weight AND
#  outperforms baseline. Coefficient directions checked
#  across all metrics regardless of support.
# ============================================================

g_baseline <- lm(
  log_mean_biomass ~ rugosity_sc + log_chla_sc,
  data = grazer_model_data)

g_q1_settgrav <- lm(
  log_mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc,
  data = grazer_model_data)

g_q1_mktgrav <- lm(
  log_mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc,
  data = grazer_model_data)

g_q1_settpop <- lm(
  log_mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_pop_sc,
  data = grazer_model_data)

cat("\n--- Q1: Pressure metric comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                      = g_baseline,
  "Baseline + settlement gravity" = g_q1_settgrav,
  "Baseline + market gravity"     = g_q1_mktgrav,
  "Baseline + settlement pop."    = g_q1_settpop
)))

# Results:
#   Baseline:           AICc = 124.30, weight = 0.497 (BEST)
#   Market gravity:     DAICc = 1.81,  weight = 0.201
#   Settlement gravity: DAICc = 2.34,  weight = 0.154
#   Settlement pop.:    DAICc = 2.43,  weight = 0.147
#
#   Baseline best supported — no pressure metric outperforms
#   ecological context. All three metrics positive and
#   non-significant — opposite direction to exploitation
#   hypothesis. No pressure term carried forward.

# ── Best Q1 model ─────────────────────────────────────────────
# Baseline best supported — no pressure term carried forward.
# Market gravity (DAICc = 1.81) within threshold — evaluated
# in sensitivity (a). Settlement gravity (DAICc = 2.34) and
# settlement pop. (DAICc = 2.43) outside threshold — not
# evaluated in sensitivity (a).
g_best_q1 <- g_baseline


# ============================================================
#  Q2 — LARVAL CONNECTIVITY
#
#  Scientific question:
#  Does larval connectivity explain additional variation
#  in grazer-detritivore biomass beyond the ecological
#  baseline?
#
#  Connectivity × pressure interaction not evaluated —
#  pressure not supported in Q1.
#
#  Criterion: AICc weight for model support; coefficients
#  and fold differences reported if supported.
# ============================================================

g_q2_conn <- lm(
  log_mean_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc,
  data = grazer_model_data)

cat("\n--- Q2: Connectivity comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                = g_best_q1,
  "Baseline + connectivity" = g_q2_conn
)))

# ── Best Q2 model ─────────────────────────────────────────────
g_best_q2 <- g_q2_conn

# ── Best model summary — reported in results ─────────────────
# Full coefficient reporting for best model across Q1-Q2.
# Q3 is a separate question and does not update this model.

cat("\n--- Best model: full summary ---\n")
summary(g_best_q2)

cat("\n--- Best model: 95% CIs ---\n")
print(confint(g_best_q2))

cat("\n--- Best model: R² ---\n")
cat(sprintf("R²     = %.3f\n", summary(g_best_q2)$r.squared))
cat(sprintf("Adj R² = %.3f\n", summary(g_best_q2)$adj.r.squared))

# ── Effect sizes: fold differences across observed range ──────
b_rug_g  <- coef(g_best_q2)["rugosity_sc"]
b_conn_g <- coef(g_best_q2)["connectivity_sc"]

rug_span_g  <- diff(range(grazer_model_data$rugosity_sc,
                          na.rm = TRUE))
conn_span_g <- diff(range(grazer_model_data$connectivity_sc,
                          na.rm = TRUE))

cat(sprintf("\nRugosity span: %.3f SD units\n", rug_span_g))
cat(sprintf("Fold difference (low vs high rugosity): %.2fx\n",
            exp(abs(b_rug_g * rug_span_g))))

cat(sprintf("\nConnectivity span: %.3f SD units\n", conn_span_g))
cat(sprintf("Fold difference (low vs high connectivity): %.2fx\n",
            exp(abs(b_conn_g * conn_span_g))))

# ── CI-based fold difference for connectivity ─────────────────
ci_conn_g <- confint(g_best_q2)["connectivity_sc",]
cat(sprintf("Connectivity fold difference 95%% CI: %.2fx to %.2fx\n",
            exp(abs(ci_conn_g[1] * conn_span_g)),
            exp(abs(ci_conn_g[2] * conn_span_g))))

# ── Model diagnostics — best Q2 model ────────────────────────
par(mfrow = c(2, 2))
plot(g_best_q2, main = "Best Q2 model diagnostics")
par(mfrow = c(1, 1))

cat("\n--- Shapiro-Wilk: best Q2 model residuals ---\n")
shapiro.test(residuals(g_best_q2))

# ── Q2 results ────────────────────────────────────────────────
# Results:
#   Baseline + connectivity: AICc = 122.82, weight = 0.677
#   Baseline:                DAICc =  1.48, weight = 0.323
#   Connectivity marginally supported (weight = 0.677,
#   DAICc = 1.48 — genuine model selection uncertainty).
#
#   Connectivity: b = -0.196, t(50) = -1.939, p = 0.058 .
#     95% CI [-0.399, 0.007] — CI just overlaps zero
#     Marginal negative effect — same direction as
#     corallivores (b = -0.221, p = 0.015).
#     Fold difference: 1.69x across observed range
#     (span = 2.667 SD units)
#     95% CI on fold difference: 2.90x to 1.02x
#   Rugosity: b = +0.241, t(50) = 2.491, p = 0.016 *
#     95% CI [0.047, 0.435]
#     Fold difference: 3.46x across observed range
#     (span = 5.151 SD units)
#   Chla: b = -0.125, t(50) = -1.201, p = 0.235 ns
#     95% CI [-0.335, 0.084]
#   R² = 0.176, adj. R² = 0.127
#   F(3,50) = 3.56, p = 0.021
#
# Diagnostics — best Q2 model:
#   Residuals vs Fitted: broadly flat — acceptable.
#     Sites 6, 21, 48 consistently flagged.
#   Q-Q: upper tail deviation at sites 6 and 48,
#     lower tail at site 21 — pattern consistent with
#     baseline. Shapiro-Wilk p = 0.003 reflects these
#     extreme sites rather than systematic non-normality.
#   Scale-Location: broadly flat — acceptable.
#   Residuals vs Leverage: site 6 remains high leverage
#     (~0.20) with large standardised residual (~4).
#     Site 48 also elevated. No sites exceed Cook's
#     distance threshold but site 6 warrants monitoring.
#   Shapiro-Wilk: W = 0.928, p = 0.003 — significant
#     departure driven by sites 6, 21, 48 (extreme
#     biomass values). Coefficient estimates reliable
#     but p-values for marginal effects (particularly
#     connectivity, p = 0.058) should be interpreted
#     with additional caution.


# ============================================================
#  Q3 — FORMAL PROTECTION
#
#  Scientific question:
#  Does MPA status explain additional variation in grazer-
#  detritivore biomass beyond the best Q2 model?
#
#  Q3 is a separate question — does not update best model.
#  MPA tested last — non-randomly placed with respect to
#  pressure and connectivity (see mpa_placement_checks.R).
#
#  Criterion: AICc weight for model support.
# ============================================================

g_q3_mpa <- lm(
  log_mean_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc +
    mpa_status,
  data = grazer_model_data)

cat("\n--- Q3: MPA comparison ---\n")
print(make_aicc_df(list(
  "Best Q2"       = g_best_q2,
  "Best Q2 + MPA" = g_q3_mpa
)))

cat("\n--- Q3: MPA coefficients ---\n")
summary(g_q3_mpa)

cat("\n--- Q3: MPA 95% CIs ---\n")
print(confint(g_q3_mpa))

# ── Raw biomass by MPA — sanity check ────────────────────────
cat("\n--- Q3: Raw biomass by MPA status ---\n")
grazer_model_data %>%
  group_by(mpa_status) %>%
  summarise(
    n              = n(),
    mean_biomass   = round(mean(mean_biomass), 1),
    median_biomass = round(median(mean_biomass), 1),
    sd_biomass     = round(sd(mean_biomass), 1),
    .groups = "drop"
  ) %>% print()

# ── Model diagnostics — Q3 model ─────────────────────────────
par(mfrow = c(2, 2))
plot(g_q3_mpa, main = "Q3 MPA model diagnostics")
par(mfrow = c(1, 1))

cat("\n--- Shapiro-Wilk: Q3 model residuals ---\n")
shapiro.test(residuals(g_q3_mpa))

# ── Q3 results ────────────────────────────────────────────────
# Results:
#   Best Q2:       AICc = 122.82, weight = 0.925
#   Best Q2 + MPA: DAICc =  5.03, weight = 0.075
#   MPA clearly not supported — connectivity model retained.
#
#   Low MPA:    b = +0.062, t(48) = 0.183, p = 0.856 ns
#     95% CI [-0.618, 0.741]
#   Medium MPA: b = -0.053, t(48) = -0.227, p = 0.822 ns
#     95% CI [-0.519, 0.414]
#   Both near zero and non-significant.
#   Raw means: none = 3818g, low = 2187g, medium = 2750g
#     Unprotected sites have highest mean biomass —
#     no protection signal.
#   Connectivity: b = -0.199, p = 0.096 . — direction
#     consistent, marginal.
#   Rugosity: b = +0.241, p = 0.019 * — stable.
#
# Diagnostics — Q3 model:
#   Residuals vs Fitted: [update from plots]
#   Q-Q: [update from plots]
#   Scale-Location: [update from plots]
#   Residuals vs Leverage: [update from plots]
#   Shapiro-Wilk: W = 0.929, p = 0.003 — consistent
#     with Q2 model; MPA addition does not improve
#     normality. Same extreme sites driving departure.
#     
# ── Best Q3 model ─────────────────────────────────────────────
# Q3 is a separate question — best Q2 model unchanged.
g_best_q3 <- g_best_q2

# ── Predicted vs observed ─────────────────────────────────────
pred_g <- predict(g_best_q3)
cat(sprintf("\nPearson r (predicted vs observed): %.3f\n",
            cor(pred_g, grazer_model_data$log_mean_biomass)))

# ── Flagged sites ─────────────────────────────────────────────
flagged_rows_g <- c(6, 21, 48)

if (length(flagged_rows_g) > 0) {
  grazer_model_data %>%
    slice(flagged_rows_g) %>%
    dplyr::select(site, mpa_status, ecoregion,
                  mean_biomass, rugosity_sc,
                  log_settlement_grav_sc, connectivity_sc) %>%
    print()
}

# Flagged sites (row indices 6, 21, 48):
#   Row 6  — ankao_s (none, Madagascar):
#     mean_biomass = 34,482g — highest grazer-detritivore
#     biomass in the dataset. Same site flagged in total
#     biomass analysis. Remote unprotected Madagascar site
#     with high rugosity — genuine ecological outlier,
#     not a data error. Primary driver of upper tail
#     Q-Q deviation and Shapiro-Wilk failure.
#     High leverage (~0.20) in diagnostic plots.
#
#   Row 21 — malinde_kiwe (none, Tanzania S-Mozambique):
#     mean_biomass = 418g — lowest grazer-detritivore
#     biomass in the dataset. Low rugosity (z = -0.740)
#     consistent with low biomass. Drives lower tail
#     Q-Q deviation. Genuine ecological observation.
#
#   Row 48 — vamizine (medium, Tanzania S-Mozambique):
#     mean_biomass = 13,916g — high biomass medium MPA
#     site with high rugosity (z = 1.32). Contributes
#     to upper tail deviation alongside ankao_s.
#     Genuine ecological observation consistent with
#     protection and habitat complexity.
#
# None warrant exclusion. Shapiro-Wilk failure driven
# by the contrast between ankao_s (extreme high) and
# malinde_kiwe (extreme low) — reflects genuine
# ecological variation at the tails of the biomass
# distribution rather than data errors or model
# misspecification.

# ============================================================
#  SPATIAL AUTOCORRELATION CHECK
#  Residuals from best Q3 model.
# ============================================================

site_coords <- location_data %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(lon = first(longitude),
            lat = first(latitude),
            .groups = "drop")

grazer_model_data_coords <- grazer_model_data %>%
  left_join(site_coords, by = "site")

coords_mat_g <- cbind(grazer_model_data_coords$lon,
                      grazer_model_data_coords$lat)
listw5_g <- nb2listw(knn2nb(knearneigh(coords_mat_g, k = 5)),
                     style = "W")

# Warning: knearneigh — identical points found, kd_tree not
# available. 4 sub-graphs expected — discontinuous sampling
# design across four countries. See total_biomass.R.

cat("\n--- Spatial autocorrelation: grazer best model ---\n")
print(moran.test(residuals(g_best_q3), listw5_g))

# Moran's I = 0.210, p = 0.001 — significant spatial
# autocorrelation in residuals.
# Stronger than total biomass (I = 0.140, p = 0.015)
# and contrasts with browsers (I = -0.078, p = 0.811)
# and corallivores (I = -0.021, p = 0.509).
# Consistent with stronger ecoregion RE support at
# baseline (DAICc = 2.93) — between-ecoregion variation
# in grazer-detritivore communities not fully captured
# by rugosity, chla, and connectivity.
# Spatial error modelling not pursued — same reasons
# as total biomass (see total_biomass.R).
# Marginal connectivity result (b = -0.196, p = 0.058)
# should be interpreted with additional caution given
# residual spatial structure.
# No spatial error modelling required.


# ============================================================
#  SENSITIVITY ANALYSES
#
#  (a) Alternative pressure metrics
#      Only metrics within DAICc < 2 of best Q1 model
#      are evaluated. Market gravity (DAICc = 1.81)
#      meets this threshold. Settlement gravity
#      (DAICc = 2.34) and settlement pop. (DAICc = 2.43)
#      do not — not evaluated.
#
#  (b) Transect-level replication
#      4 zero-biomass transects (1.6%) filtered before
#      analysis — log transformation requires positive
#      values. lmer with (1 | site) random intercept.
# ============================================================

# ── (a) Market gravity ────────────────────────────────────────
cat("\n--- Sensitivity (a): market gravity ---\n")

g_sens_mg_conn <- lm(
  log_mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc + connectivity_sc,
  data = grazer_model_data)

cat("\nQ2 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + MG"        = g_q1_mktgrav,
  "Baseline + MG + conn" = g_sens_mg_conn
)))

cat("\nQ2 — connectivity coefficients:\n")
print(summary(g_sens_mg_conn)$coefficients)
print(confint(g_sens_mg_conn))

g_sens_mg_conn_int <- lm(
  log_mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc + connectivity_sc +
    log_market_gravity_sc:connectivity_sc,
  data = grazer_model_data)

print(make_aicc_df(list(
  "Baseline + MG"              = g_q1_mktgrav,
  "Baseline + MG + conn"       = g_sens_mg_conn,
  "Baseline + MG + conn + int" = g_sens_mg_conn_int
)))

summary(g_sens_mg_conn_int)
confint(g_sens_mg_conn_int)

g_sens_mg_mpa <- lm(
  log_mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc + connectivity_sc + mpa_status,
  data = grazer_model_data)

cat("\nQ3 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + MG + conn"       = g_sens_mg_conn,
  "Baseline + MG + conn + MPA" = g_sens_mg_mpa
)))

cat("\nQ3 — MPA coefficients:\n")
print(summary(g_sens_mg_mpa)$coefficients)
print(confint(g_sens_mg_mpa))

# ── Sensitivity (a): market gravity results ───────────────────
# Q2:
#   Baseline + MG + conn: AICc = 123.65, weight = 0.773
#   Baseline + MG:        DAICc = 2.45,  weight = 0.227
#   Connectivity: b = -0.224, t = -2.178, p = 0.034 *
#     95% CI [-0.431, -0.017] — significant with MG included
#     Stronger than primary analysis (p = 0.058) —
#     market gravity absorbs some shared variance,
#     sharpening connectivity signal.
#   Market gravity: b = +0.133, p = 0.217 ns
#
# Q3:
#   Baseline + MG + conn + MPA: DAICc = 5.35, weight = 0.064
#   Baseline + MG + conn:       weight = 0.936
#   MPA clearly not supported — consistent with primary Q3.
#   Low MPA: b = +0.044, p = 0.896 ns
#   Medium MPA: b = -0.030, p = 0.897 ns
#
# Sensitivity (a) conclusion:
#   Connectivity negative and significant with MG included
#   (p = 0.034) — strengthens primary finding (p = 0.058).
#   MPA not supported under either specification.
#   Q2 and Q3 conclusions robust to metric choice.


# ── (b) Transect-level replication ───────────────────────────
grazer_transect_data_nz <- grazer_transect_data %>%
  filter(transect_grazer_biomass > 0) %>%
  mutate(log_transect_biomass = log(transect_grazer_biomass))

cat("Transects retained:", nrow(grazer_transect_data_nz),
    "of", nrow(grazer_transect_data), "\n")

g_trans_null <- lmer(
  log_transect_biomass ~ 1 + (1 | site),
  data = grazer_transect_data_nz, REML = TRUE)

g_trans_baseline <- lmer(
  log_transect_biomass ~ rugosity_sc + log_chla_sc +
    (1 | site),
  data = grazer_transect_data_nz, REML = TRUE)

g_trans_pressure <- lmer(
  log_transect_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc + (1 | site),
  data = grazer_transect_data_nz, REML = TRUE)

g_trans_conn <- lmer(
  log_transect_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc + (1 | site),
  data = grazer_transect_data_nz, REML = TRUE)

g_trans_mpa <- lmer(
  log_transect_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc + mpa_status + (1 | site),
  data = grazer_transect_data_nz, REML = TRUE)

# ── Diagnostics on baseline model ────────────────────────────
par(mfrow = c(1, 2))
plot(fitted(g_trans_baseline), residuals(g_trans_baseline),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Grazer transect baseline: Residuals vs Fitted")
abline(h = 0, lty = 2, col = "grey60")
lines(lowess(fitted(g_trans_baseline),
             residuals(g_trans_baseline)), col = "red")
qqnorm(residuals(g_trans_baseline),
       main = "Grazer transect baseline: Q-Q")
qqline(residuals(g_trans_baseline), col = "red")
par(mfrow = c(1, 1))

# ── AICc comparison — ML refits ──────────────────────────────
g_trans_null_ml     <- update(g_trans_null,     REML = FALSE)
g_trans_baseline_ml <- update(g_trans_baseline, REML = FALSE)
g_trans_pressure_ml <- update(g_trans_pressure, REML = FALSE)
g_trans_conn_ml     <- update(g_trans_conn,     REML = FALSE)
g_trans_mpa_ml      <- update(g_trans_mpa,      REML = FALSE)

cat("\n--- Sensitivity (b): transect model comparison ---\n")
print(make_aicc_df(list(
  "Null"                = g_trans_null_ml,
  "Baseline"            = g_trans_baseline_ml,
  "Baseline + pressure" = g_trans_pressure_ml,
  "Baseline + conn"     = g_trans_conn_ml,
  "Best + MPA"          = g_trans_mpa_ml
)))

cat("\n--- Sensitivity (b): connectivity model coefficients ---\n")
summary(g_trans_conn)

# ── ICC ───────────────────────────────────────────────────────
vc_g         <- as.data.frame(VarCorr(g_trans_conn))
site_var_g   <- vc_g$vcov[vc_g$grp == "site"]
resid_var_g  <- vc_g$vcov[vc_g$grp == "Residual"]
icc_g        <- site_var_g / (site_var_g + resid_var_g)
cat(sprintf("\nICC = %.3f — %.1f%% of variance attributable
to between-site differences\n", icc_g, icc_g * 100))

# Results:
#   Baseline + conn:  AICc = 625.50, weight = 0.453
#   Baseline:         DAICc = 0.71,  weight = 0.317
#   Baseline + pres:  DAICc = 2.39,  weight = 0.137
#   Best + MPA:       DAICc = 4.14,  weight = 0.057
#   Null:             DAICc = 5.04,  weight = 0.036
#
#   Connectivity top-ranked at transect level —
#   consistent with site-level Q2 (weight = 0.677).
#   Higher model uncertainty at transect level
#   (weight = 0.453 vs 0.677) but direction consistent.
#   Pressure not supported (DAICc = 2.39) — consistent
#   with Q1. MPA not supported (DAICc = 4.14) —
#   consistent with Q3.
#
#   Connectivity: b = -0.146, t = -1.637
#     Consistent direction with site-level (b = -0.196).
#   Rugosity: b = +0.235, t = 2.729 ** — stable.
#   ICC = 0.295 — 29.5% of variance between-site;
#     higher than total biomass (0.205), consistent
#     with stronger ecoregion signal for this group.

# Diagnostics — transect baseline:
#   Residuals vs Fitted: broadly flat with minor
#     downward trend — acceptable for n = 239.
#     Slight funnel pattern at high fitted values.
#   Q-Q: upper tail deviation — a handful of high-
#     biomass transects pull above theoretical line.
#     Lower tail follows line closely. Acceptable
#     for lmer with large n.
#   Gaussian lmer structure confirmed appropriate.

# ============================================================
#  RESULTS SUMMARY
# ============================================================

cat("\n--- Grazer-detritivore results summary ---\n")
tribble(
  ~Question,   ~Result,           ~Key_finding,
  "Q1",        "Baseline best",   "weight = 0.497; all metrics positive ns — no pressure signal, wrong direction",
  "Q1 rug",    "Significant",     "b = +0.241, p = 0.016; fold diff = 3.46x (span = 5.151 SD)",
  "Q2 conn",   "Suggestive",      "weight = 0.677, DAICc = 1.48; b = -0.196, p = 0.058; fold diff = 1.69x (CI: 2.90x–1.02x) — NEGATIVE, marginal",
  "Q3 MPA",    "Not supported",   "DAICc = 5.03, weight = 0.075; connectivity model retained",
  "Moran",     "Significant SAC", "I = 0.210, p = 0.001 — strongest of all groups; additional caution for connectivity",
  "Pearson r", "0.420",           "log scale — moderate fit",
  "Shapiro",   "p = 0.003",       "significant departure — upper tail; interpret p-values with caution",
  "Sens (a)",  "Consistent",      "MG: conn b = -0.224 p = 0.034 — strengthens; MPA DAICc = 5.35 ns",
  "Sens (b)",  "Consistent",      "conn top-ranked (weight = 0.453); b = -0.146 negative — consistent; ICC = 0.295"
) %>% print()


# ── End of script ─────────────────────────────────────────────