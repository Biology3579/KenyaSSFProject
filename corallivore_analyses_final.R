# ============================================================
#  DRIVERS OF CORALLIVORE BIOMASS
#  Functional Group Analysis: Corallivores
# ============================================================
#
#  Scientific questions:
#
#  Q1 — Human pressure
#       Does human pressure explain variation in corallivore
#       biomass beyond local ecological context, and which
#       spatial metric best captures SSF exploitation
#       intensity for this functional group?
#
#  Q2 — Larval connectivity
#       Does larval connectivity explain additional variation
#       in corallivore biomass beyond the best Q1 model?
#       Connectivity × pressure interaction not evaluated —
#       pressure not supported in Q1.
#
#  Q3 — Formal protection
#       Does MPA status explain additional variation in
#       corallivore biomass beyond the best Q2 model?
#       Q3 is a separate question — does not update best model.
#       MPA tested last — non-randomly placed with respect
#       to pressure and connectivity (see
#       mpa_placement_checks.R).
#
#  Rationale for sequence:
#       Identical to total biomass and browsers — pressure
#       first, connectivity second (extends Warmuth et al.
#       2024), MPA last as governance response downstream
#       of both.
#
#  Baseline model (fixed a priori, never tested):
#       biomass ~ rugosity_sc + log_chla_sc
#
#  Key difference from total biomass:
#       Corallivore biomass may have zeros at site level
#       (check below). Tweedie GLM used throughout for
#       consistency across all functional group analyses.
#
#  Sensitivity analyses:
#       (a) Alternative pressure metrics
#       (b) Transect-level replication (Tweedie GLMM)
# ============================================================

source(here::here("predictor_preparation.R"))


# ============================================================
#  RESPONSE VARIABLE — DATA AGGREGATION AND CHECKS
# ============================================================

# ── Transect-level aggregation ────────────────────────────────
coralliv_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_coralliv_biomass = sum(
      ifelse(trophic_group == "corallivores", tot_wt_g, 0),
      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(site = as.factor(site))

cat("Corallivore transects:", nrow(coralliv_transects), "\n")
cat("Sites:",                 n_distinct(coralliv_transects$site), "\n")

# ── Transect-level dataset ────────────────────────────────────
coralliv_transect_data <- coralliv_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(coralliv_transect_data$transect_coralliv_biomass == 0),
    "/", nrow(coralliv_transect_data),
    sprintf("(%.3f)\n",
            mean(coralliv_transect_data$transect_coralliv_biomass == 0)))

# ── Site-level dataset ────────────────────────────────────────
coralliv_model_data <- coralliv_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    mean_biomass           = mean(transect_coralliv_biomass,
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
    site      = as.factor(site),
    ecoregion = as.factor(ecoregion),
    mpa_status = factor(mpa_status,
                        levels  = c("none", "low", "medium"),
                        ordered = FALSE)
  )

cat("\nCorallivore model data:", nrow(coralliv_model_data), "sites\n")

# ── Data checks ───────────────────────────────────────────────
coralliv_model_data %>%
  dplyr::select(site, rugosity_sc, log_chla_sc,
                log_settlement_grav_sc, connectivity_sc,
                mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  { if (nrow(.) > 0) { warning("NAs in predictors:"); print(.) }
    else cat("NA check passed.\n") }

cat("\nZeros in mean_biomass:",
    sum(coralliv_model_data$mean_biomass == 0), "\n")
cat("Site-level zero proportion:",
    round(mean(coralliv_model_data$mean_biomass == 0), 3), "\n")
cat("\nResponse summary:\n")
print(summary(coralliv_model_data$mean_biomass))
cat("\nMPA status counts:\n")
print(table(coralliv_model_data$mpa_status))


# ============================================================
#  MODEL FAMILY SELECTION
#
#  Corallivore biomass may have zeros at site level —
#  check above. Tweedie handles zeros natively without
#  offset and is used for consistency across all functional
#  group analyses.
#
#  Two families tested on baseline structure:
#    Gaussian log: lm() on log(mean_biomass + offset)
#    Tweedie:      glmmTMB() on raw mean_biomass, log link
# ============================================================

# ── Gaussian log ──────────────────────────────────────────────
coralliv_lm_baseline <- lm(
  log(mean_biomass + 0.01) ~ rugosity_sc + log_chla_sc,
  data = coralliv_model_data)

par(mfrow = c(2, 2))
plot(coralliv_lm_baseline, main = "Gaussian log — baseline")
par(mfrow = c(1, 1))

# ── Tweedie ───────────────────────────────────────────────────
coralliv_tw_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data)

coralliv_tw_res <- simulateResiduals(coralliv_tw_baseline, n = 1000)
plot(coralliv_tw_res)
testZeroInflation(coralliv_tw_res)
testDispersion(coralliv_tw_res)

# Family selection:
#
#   Gaussian log: REJECTED
#     Residuals vs Fitted: broadly flat but minor curve at
#       low fitted values. Sites 13 and 19 pull residuals
#       to -2 — lower tail deviation.
#     Q-Q: site 13 falls below theoretical line at lower
#       tail — moderate deviation.
#     Scale-Location: slight downward trend, sites 13 and
#       19 elevated. Acceptable but not ideal.
#     Residuals vs Leverage: site 13 approaches Cook's
#       distance threshold.
#
#   Tweedie (log link): SELECTED
#     Handles zeros natively without offset.
#     DHARMa diagnostics (n = 1000 simulations):
#       KS test:        p = 0.795 — excellent fit
#       Dispersion:     p = 0.804, ratio = 0.892 — acceptable
#       Zero inflation: p = 1.000, ratio = 0 — no zeros
#       Outlier test:   p = 1.000 — no outliers
#     Selected for fit quality and consistency across
#     all functional group analyses.
#
#   Proceed: glmmTMB(family = tweedie(link = "log")) on
#   raw mean_biomass throughout all corallivore analyses.


# ============================================================
#  RANDOM EFFECT STRUCTURE
#
#  Tested on baseline model — same rationale as total
#  biomass and browsers. See total_biomass.R for full
#  justification.
# ============================================================

coralliv_re_null <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data)

coralliv_re_ecoregion <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + (1 | ecoregion),
  family = tweedie(link = "log"),
  data   = coralliv_model_data)

cat("\n--- Random effect structure: ecoregion ---\n")
print(make_aicc_df(list(
  "No RE"           = coralliv_re_null,
  "(1 | ecoregion)" = coralliv_re_ecoregion
)))

# Ecoregion RE not supported (DAICc = 2.54, weight = 0.220)
# — consistent with total biomass (DAICc = 2.25) and
# browsers (DAICc = 2.13). All corallivore models fitted
# without RE throughout.


# ============================================================
#  Q1 — HUMAN PRESSURE
#
#  Scientific question:
#  Does human pressure explain variation in corallivore
#  biomass beyond local ecological context, and which
#  spatial metric best captures SSF exploitation intensity?
#
#  Approach: AICc comparison of baseline vs baseline + each
#  pressure metric. Best metric = highest AICc weight AND
#  outperforms baseline. Coefficient directions checked
#  across all metrics regardless of support.
# ============================================================

c_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data)

c_q1_settgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data)

c_q1_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data)

c_q1_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_pop_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data)

cat("\n--- Q1: Pressure metric comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                      = c_baseline,
  "Baseline + settlement gravity" = c_q1_settgrav,
  "Baseline + market gravity"     = c_q1_mktgrav,
  "Baseline + settlement pop."    = c_q1_settpop
)))

# Results:
#   Baseline:           AICc = 560.23, weight = 0.499 (BEST)
#   Settlement gravity: DAICc = 1.86,  weight = 0.197
#   Settlement pop.:    DAICc = 2.26,  weight = 0.161
#   Market gravity:     DAICc = 2.50,  weight = 0.143
#
#   Baseline best supported — no pressure metric outperforms
#   ecological context for corallivores. High model selection
#   uncertainty (all metrics within DAICc < 2.5).
#   No pressure term carried forward.

# ── Best Q1 model ─────────────────────────────────────────────
# Baseline best supported — no pressure term carried forward.
# Coefficient directions and robustness of Q2-Q3 conclusions
# checked in sensitivity analyses (a) below.
c_best_q1 <- c_baseline


# ============================================================
#  Q2 — LARVAL CONNECTIVITY
#
#  Scientific question:
#  Does larval connectivity explain additional variation
#  in corallivore biomass beyond the ecological baseline?
#
#  Connectivity × pressure interaction not evaluated —
#  pressure not supported in Q1. Testing an interaction
#  without a supported main effect would be
#  overparameterised and uninterpretable.
#
#  Criterion: AICc weight for model support; coefficients
#  and fold differences reported if supported.
# ============================================================

c_q2_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data)

cat("\n--- Q2: Connectivity comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                = c_best_q1,
  "Baseline + connectivity" = c_q2_conn
)))

# ── Best Q2 model ─────────────────────────────────────────────
c_best_q2 <- c_q2_conn

# ── Best model summary — reported in results ─────────────────
# Full coefficient reporting for the best model identified
# across Q1-Q2. This is the primary reported model.
# Q3 is a separate question and does not update this model.

cat("\n--- Best model: full summary ---\n")
summary(c_best_q2)

cat("\n--- Best model: 95% CIs ---\n")
print(confint(c_best_q2))

# ── Effect sizes: fold differences across observed range ──────
b_rug_c  <- fixef(c_best_q2)$cond["rugosity_sc"]
b_conn_c <- fixef(c_best_q2)$cond["connectivity_sc"]

rug_span_c  <- diff(range(coralliv_model_data$rugosity_sc,
                          na.rm = TRUE))
conn_span_c <- diff(range(coralliv_model_data$connectivity_sc,
                          na.rm = TRUE))

cat(sprintf("\nRugosity span: %.3f SD units\n", rug_span_c))
cat(sprintf("Fold difference (low vs high rugosity): %.2fx\n",
            exp(abs(b_rug_c * rug_span_c))))

cat(sprintf("\nConnectivity span: %.3f SD units\n", conn_span_c))
cat(sprintf("Fold difference (low vs high connectivity): %.2fx\n",
            exp(abs(b_conn_c * conn_span_c))))

# ── CI-based fold difference for connectivity ─────────────────
ci_conn_c <- confint(c_best_q2)["connectivity_sc",
                                c("2.5 %", "97.5 %")]
cat(sprintf("Connectivity fold difference 95%% CI: %.2fx to %.2fx\n",
            exp(abs(ci_conn_c[1] * conn_span_c)),
            exp(abs(ci_conn_c[2] * conn_span_c))))

# ── DHARMa diagnostics — best Q2 model ───────────────────────
cat("\n--- DHARMa diagnostics: best Q2 model ---\n")
c_best_q2_sim <- simulateResiduals(c_best_q2, n = 1000)
plot(c_best_q2_sim)
testOutliers(c_best_q2_sim)

# ── Q2 results ────────────────────────────────────────────────
# Results:
#   Baseline + connectivity: AICc = 557.13, weight = 0.825
#   Baseline:                DAICc = 3.10,  weight = 0.175
#   Connectivity strongly supported (weight = 0.825).
#
#   Connectivity: b = -0.221, z = -2.44, p = 0.015 *
#     95% CI [-0.399, -0.043]
#     Significant NEGATIVE effect — well-connected sites
#     have lower corallivore biomass. Contrasts with
#     browsers (b = +0.421, positive direction).
#     Fold difference: 1.80x across observed range
#     (span = 2.667 SD units)
#     95% CI on fold difference: 2.90x to 1.12x
#     Note: CI spans from larger to smaller fold difference
#     reflecting the negative direction.
#   Rugosity: b = +0.108, z = 1.24, p = 0.214 ns
#     95% CI [-0.062, 0.278]
#     Fold difference: 1.74x across observed range
#   Chla: b = +0.121, z = 1.29, p = 0.196 ns
#     95% CI [-0.062, 0.304]
#
# DHARMa diagnostics — best Q2 model:
#   KS test:        p = 0.649 — no significant deviation
#   Dispersion:     p = 0.778 — no significant deviation
#   Outlier test:   p = 1.000 — 0 outliers / 54 obs
#   Residuals vs predicted: no significant problems detected
#   Overall: good fit.

# ============================================================
#  Q3 — FORMAL PROTECTION
#
#  Scientific question:
#  Does MPA status explain additional variation in
#  corallivore biomass beyond the best Q2 model?
#
#  Q3 is a separate question — does not update best model.
#  MPA tested last — non-randomly placed with respect to
#  pressure and connectivity (see mpa_placement_checks.R).
#
#  Criterion: AICc weight for model support.
#  p-values for coefficient direction if supported.
# ============================================================

c_q3_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = coralliv_model_data)

cat("\n--- Q3: MPA comparison ---\n")
print(make_aicc_df(list(
  "Best Q2"       = c_best_q2,
  "Best Q2 + MPA" = c_q3_mpa
)))

cat("\n--- Q3: MPA coefficients ---\n")
summary(c_q3_mpa)

cat("\n--- Q3: MPA 95% CIs ---\n")
print(confint(c_q3_mpa))

# ── Raw biomass by MPA — sanity check ────────────────────────
cat("\n--- Q3: Raw biomass by MPA status ---\n")
coralliv_model_data %>%
  group_by(mpa_status) %>%
  summarise(
    n              = n(),
    mean_biomass   = round(mean(mean_biomass), 1),
    median_biomass = round(median(mean_biomass), 1),
    sd_biomass     = round(sd(mean_biomass), 1),
    .groups = "drop"
  ) %>% print()

# ── DHARMa diagnostics — Q3 model ────────────────────────────
cat("\n--- DHARMa diagnostics: Q3 model ---\n")
c_q3_sim <- simulateResiduals(c_q3_mpa, n = 1000)
plot(c_q3_sim)
testOutliers(c_q3_sim)

# ── Q3 results ────────────────────────────────────────────────
# Results:
#   Best Q2:       AICc = 557.13, weight = 0.549
#   Best Q2 + MPA: DAICc =  0.39, weight = 0.451
#   Genuine model selection uncertainty — MPA does not
#   substantially improve on connectivity baseline.
#
#   Low MPA:    b = +0.622, z = 2.22, p = 0.027 *
#     95% CI [0.072, 1.171] — significant but likely
#     artefactual. Low MPA sites in narrow high-connectivity
#     range; site characteristics elevate biomass
#     independently of protection status.
#   Medium MPA: b = +0.248, z = 1.30, p = 0.192 ns
#     95% CI [-0.125, 0.621]
#     Raw mean (68.4g) essentially identical to
#     unprotected sites (67.6g).
#   Connectivity: b = -0.337, z = -3.37, p < 0.001 ***
#     Strengthens and remains highly significant once
#     MPA included — negative connectivity signal robust
#     and not confounded by MPA placement.
#
#   Conclusion: MPA not supported (DAICc = 0.39).
#   Connectivity model retained as best model throughout.
#
# DHARMa diagnostics — Q3 model:
#   KS test:        p = 0.879 — no significant deviation
#   Dispersion:     p = 0.848 — no significant deviation
#   Outlier test:   p = 1.000 — 0 outliers / 54 obs
#   Residuals vs predicted: quantile deviations detected
#     (red curve, upper quantile elevated across fitted range).
#     KS and dispersion tests pass — mild heteroscedasticity
#     introduced by discrete MPA categories, consistent
#     with browser Q3 pattern. Acceptable overall.

# Pearson r (predicted vs observed): 0.409
# Flagged sites: kifinge (row 13) and makunga_n (row 19)
#   Both low biomass sites (~9g), different MPA status
#   (medium and none), Tanzania S-Mozambique ecoregion.
#   Low corallivore biomass consistent with habitat
#   characteristics — not data errors.

# ── Best Q3 model ─────────────────────────────────────────────
# Q3 is a separate question — best Q2 model unchanged.
c_best_q3 <- c_best_q2

# ── Predicted vs observed ─────────────────────────────────────
pred_c <- predict(c_best_q3, type = "response")
cat(sprintf("\nPearson r (predicted vs observed): %.3f\n",
            cor(pred_c, coralliv_model_data$mean_biomass)))

# ── Flagged sites ─────────────────────────────────────────────
# Identify influential sites from diagnostic plots
# Row indices from diagnostic plots — update after running
flagged_rows_c <- c(13, 19)  # update from plots

coralliv_model_data %>%
  slice(flagged_rows_c) %>%
  dplyr::select(site, mpa_status, ecoregion,
                mean_biomass, rugosity_sc,
                log_settlement_grav_sc, connectivity_sc) %>%
  print()


# ============================================================
#  SPATIAL AUTOCORRELATION CHECK
#  Pearson residuals from best Q3 model.
# ============================================================

site_coords <- location_data %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(lon = first(longitude),
            lat = first(latitude),
            .groups = "drop")

coralliv_model_data_coords <- coralliv_model_data %>%
  left_join(site_coords, by = "site")

coords_mat_c <- cbind(coralliv_model_data_coords$lon,
                      coralliv_model_data_coords$lat)
listw5_c <- nb2listw(knn2nb(knearneigh(coords_mat_c, k = 5)),
                     style = "W")

# Warning: knearneigh — identical points found, kd_tree not
# available. 4 sub-graphs expected — discontinuous sampling
# design across four countries. See total_biomass.R.

cat("\n--- Spatial autocorrelation: corallivore best model ---\n")
print(moran.test(residuals(c_best_q3, type = "pearson"),
                 listw5_c))

# Moran's I = -0.021, p = 0.509 — no significant spatial
# autocorrelation in residuals.
# Consistent with browsers (I = -0.078, p = 0.811).
# Contrasts with total biomass (I = 0.140, p = 0.015).
# Connectivity predictor adequately captures spatial
# structure in corallivore biomass.
# No spatial error modelling required.


# ============================================================
#  SENSITIVITY ANALYSES
#
#  (a) Alternative pressure metrics
#      Only metrics within DAICc < 2 of the best Q1 model
#      are evaluated. Settlement gravity (DAICc = 1.86)
#      meets this threshold. Settlement pop. (DAICc = 2.26)
#      and market gravity (DAICc = 2.50) do not — not
#      evaluated.
#
#  (b) Transect-level replication
#      Mirrors Q1-Q3 at transect level using Tweedie GLMM
#      with (1 | site) random intercept.
# ============================================================

# ── (a) Settlement gravity ────────────────────────────────────
cat("\n--- Sensitivity (a): settlement gravity ---\n")

c_sens_sg_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc + connectivity_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data)

cat("\nQ2 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + SG"        = c_q1_settgrav,
  "Baseline + SG + conn" = c_sens_sg_conn
)))

cat("\nQ2 — connectivity coefficients:\n")
print(summary(c_sens_sg_conn)$coefficients$cond)
print(confint(c_sens_sg_conn))

c_sens_sg_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc + connectivity_sc + mpa_status,
  family = tweedie(link = "log"),
  data   = coralliv_model_data)

cat("\nQ3 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + SG + conn"       = c_sens_sg_conn,
  "Baseline + SG + conn + MPA" = c_sens_sg_mpa
)))

cat("\nQ3 — MPA coefficients:\n")
print(summary(c_sens_sg_mpa)$coefficients$cond)
print(confint(c_sens_sg_mpa))

# Q3 — settlement gravity sensitivity results:
#   Baseline + SG + conn + MPA: AICc = 557.55, weight = 0.616
#   Baseline + SG + conn:       DAICc =  0.95,  weight = 0.384
#   Genuine uncertainty — consistent with primary Q3 result.
#
#   Connectivity: b = -0.359, z = -3.674, p < 0.001 ***
#     95% CI [-0.550, -0.167] — strengthens once MPA added,
#     consistent with primary analysis.
#   Low MPA:    b = +0.589, z = 2.143, p = 0.032 *
#     95% CI [0.050, 1.128] — significant but same caveat
#     as primary analysis: likely artefactual given low MPA
#     site characteristics.
#   Medium MPA: b = +0.453, z = 2.049, p = 0.040 *
#     95% CI [0.020, 0.886] — note this is significant here
#     whereas it was ns in primary analysis (b = +0.248,
#     p = 0.192). Difference driven by settlement gravity
#     absorbing some shared variance — interpret cautiously.
#     Raw medium MPA mean (68.4g) still essentially identical
#     to unprotected sites (67.6g).
#   Settlement gravity: b = +0.201, p = 0.088 ns — consistent
#     with Q1 finding that pressure does not structure
#     corallivore biomass.


# ── (b) Transect-level replication ───────────────────────────
# Mirrors Q1-Q3 at transect level using standard Tweedie
# GLMM with (1 | site) random intercept — consistent with
# site-level family selection.

c_trans_null <- glmmTMB(
  transect_coralliv_biomass ~ 1 + (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data)

c_trans_baseline <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc + log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data)

c_trans_pressure <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data)

c_trans_conn <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data)

c_trans_mpa <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc + mpa_status + (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data)

cat("\n--- Sensitivity (b): transect model comparison ---\n")
print(make_aicc_df(list(
  "Null"                = c_trans_null,
  "Baseline"            = c_trans_baseline,
  "Baseline + pressure" = c_trans_pressure,
  "Baseline + conn"     = c_trans_conn,
  "Best + MPA"          = c_trans_mpa
)))

cat("\n--- Sensitivity (b): best transect model coefficients ---\n")
summary(c_trans_conn)

# ── Site random intercept SD ──────────────────────────────────
vc_c      <- VarCorr(c_trans_conn)
site_sd_c <- sqrt(as.numeric(vc_c$cond$site))
cat(sprintf("\nSite random intercept SD = %.3f\n", site_sd_c))

# Results:
#   Best + MPA:      AICc = 2319.06, weight = 0.676
#   Baseline + conn: DAICc =  1.90,  weight = 0.261
#   Baseline:        DAICc =  6.10,  weight = 0.032
#   Null:            DAICc =  7.33,  weight = 0.017
#   Baseline + pres: DAICc =  7.71,  weight = 0.014
#
#   MPA + conn best supported at transect level —
#   stronger than site-level (0.676 vs 0.549).
#   Connectivity: b = -0.231, z = -2.60, p = 0.009 **
#     Negative direction robust at transect level —
#     consistent with site-level Q2 (b = -0.221).
#   Site random intercept SD = 0.493 — between-site
#     variance moderate; (1 | site) justified.


# ============================================================
#  RESULTS SUMMARY
# ============================================================

cat("\n--- Corallivore results summary ---\n")
tribble(
  ~Question,   ~Result,                    ~Key_finding,
  "Q1",        "Baseline best",            "weight = 0.499; no metric outperforms baseline; all metrics ns, inconsistent directions",
  "Q1 chla",   "p = 0.049 (baseline only)","b = +0.190 in baseline; weakens to b = +0.121, p = 0.196 once connectivity added",
  "Q2 conn",   "Supported",                "weight = 0.825, DAICc = 3.10; b = -0.221, p = 0.015; fold diff = 1.80x (CI: 2.90x–1.12x) — NEGATIVE",
  "Q3 MPA",    "Not supported",            "DAICc = 0.39, weight = 0.451 — genuine uncertainty; low MPA artefactual; connectivity model retained",
  "Moran",     "No SAC",                   "I = -0.021, p = 0.509 — no spatial autocorrelation",
  "Pearson r", "0.409",                    "moderate fit",
  "Sens (a)",  "Consistent",               "SG: conn b = -0.231 p = 0.010 — robust; MPA uncertain (DAICc = 0.95)",
  "Sens (b)",  "Consistent",               "Best + MPA weight = 0.676, conn b = -0.231 p = 0.009; site SD = 0.493"
) %>% print()


# ── End of script ─────────────────────────────────────────────