# ============================================================
#  DRIVERS OF BROWSER BIOMASS
#  Functional Group Analysis: Browsers
# ============================================================
#
#  Scientific questions:
#
#  Q1 — Human pressure
#       Does human pressure explain variation in browser
#       biomass beyond local ecological context, and which
#       spatial metric best captures human pressure for this functional group?
#
#  Q2 — Larval connectivity
#       Does larval connectivity explain additional variation
#       in browser biomass beyond the best-supported
#       pressure model?
#       Connectivity × pressure interaction not evaluated —
#       pressure not supported in Q1, interaction would
#       be overparameterised without a pressure main effect.
#
#  Q3 — Formal protection
#       Does MPA status explain additional variation in
#       browser biomass beyond the best-supported
#       connectivity model?
#       MPA tested last — non-randomly placed with respect
#       to pressure and connectivity (see
#       mpa_placement_checks.R). Effect confounded unless
#       both controlled first.
#
#  Rationale for sequence:
#       Identical to total biomass — pressure first,
#       connectivity second (extends Warmuth et al. 2024),
#       MPA last as governance response downstream of both.
#
#  Baseline model (fixed a priori, never tested):
#       biomass ~ rugosity_sc + log_chla_sc
#
#  Key difference from total biomass:
#       Browser biomass has ~11% zeros at site level.
#       Tweedie GLM selected over Gaussian log — handles
#       zeros natively without offset (see family
#       selection below). All models use
#       glmmTMB(family = tweedie(link = "log")) on
#       raw mean_biomass throughout.
#
#  Sensitivity analyses:
#       (a) Alternative pressure metrics — robustness of
#           Q2 and Q3 conclusions to Q1 metric uncertainty
#       (b) Transect-level replication (Tweedie GLMM)
# ============================================================

source(here::here("predictor_preparation.R"))


# ============================================================
#  RESPONSE VARIABLE — DATA AGGREGATION AND CHECKS
# ============================================================

# ── Transect-level aggregation ────────────────────────────────
browser_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_browser_biomass = sum(
      ifelse(trophic_group == "browsers", tot_wt_g, 0),
      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(site = as.factor(site))

cat("Browser transects:", nrow(browser_transects), "\n")
cat("Sites:",             n_distinct(browser_transects$site), "\n")

# ── Transect-level dataset ────────────────────────────────────
browser_transect_data <- browser_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(browser_transect_data$transect_browser_biomass == 0),
    "/", nrow(browser_transect_data),
    sprintf("(%.3f)\n",
            mean(browser_transect_data$transect_browser_biomass == 0)))

# ── Site-level dataset ────────────────────────────────────────
browser_model_data <- browser_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    mean_biomass           = mean(transect_browser_biomass,
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

cat("\nBrowser model data:", nrow(browser_model_data), "sites\n")

# ── Data checks ───────────────────────────────────────────────
browser_model_data %>%
  dplyr::select(site, rugosity_sc, log_chla_sc,
                log_settlement_grav_sc, connectivity_sc,
                mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  { if (nrow(.) > 0) { warning("NAs in predictors:"); print(.) }
    else cat("NA check passed.\n") }

cat("\nZeros in mean_biomass:",
    sum(browser_model_data$mean_biomass == 0), "\n")
cat("Site-level zero proportion:",
    round(mean(browser_model_data$mean_biomass == 0), 3), "\n")
cat("\nResponse summary:\n")
print(summary(browser_model_data$mean_biomass))
cat("\nMPA status counts:\n")
print(table(browser_model_data$mpa_status))


# ============================================================
#  MODEL FAMILY SELECTION
#
#  Browser biomass has ~11% zeros at site level.
#  Gaussian log requires an offset constant for zeros —
#  log-transformed distribution is bimodal with zero sites
#  clustering at log(offset), separated from the main
#  distribution. No offset constant resolves this.
#  Tweedie handles zeros natively without offset.
#
#  Two families tested on baseline structure:
#    Gaussian log: lm() on log(mean_biomass + offset)
#    Tweedie:      glmmTMB() on raw mean_biomass, log link
#
#  DHARMa residual diagnostics used for Tweedie.
#  Standard plot() used for Gaussian log.
# ============================================================

# ── Gaussian log ──────────────────────────────────────────────
browser_lm_baseline <- lm(
  log(mean_biomass + 0.01) ~ rugosity_sc + log_chla_sc,
  data = browser_model_data)

par(mfrow = c(2, 2))
plot(browser_lm_baseline, main = "Gaussian log — baseline")
par(mfrow = c(1, 1))

# ── Tweedie ───────────────────────────────────────────────────
browser_tw_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data)

browser_tw_res <- simulateResiduals(browser_tw_baseline, n = 1000)
plot(browser_tw_res)
testZeroInflation(browser_tw_res)
testDispersion(browser_tw_res)

# Family selection:
#
#   Gaussian log: REJECTED
#     Residuals vs Fitted: strong downward curve at low
#       fitted values — zero sites pulling residuals to
#       extreme values, no offset constant resolves this.
#     Q-Q: severe lower tail deviation.
#     Scale-Location: strong downward trend —
#       heteroscedasticity throughout fitted range.
#
#   Tweedie (log link): SELECTED
#     Handles zeros natively without offset.
#     DHARMa diagnostics (n = 1000 simulations):
#       KS test:        p = 0.762 — good fit
#       Dispersion:     p = 0.222, ratio = 1.615 — acceptable
#       Zero inflation: p = 0.902, ratio = 0.872 — not needed
#       Outlier test:   p = 0.102 — no significant outliers
#
#   Proceed: glmmTMB(family = tweedie(link = "log")) on
#   raw mean_biomass throughout all browser analyses.


# ============================================================
#  RANDOM EFFECT STRUCTURE
#
#  Tested on baseline model — same rationale as total
#  biomass. Tweedie family used throughout.
#  See total_biomass.R for full justification.
# ============================================================

browser_re_null <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data)

browser_re_ecoregion <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + (1 | ecoregion),
  family = tweedie(link = "log"),
  data   = browser_model_data)

cat("\n--- Random effect structure: ecoregion ---\n")
print(make_aicc_df(list(
  "No RE"           = browser_re_null,
  "(1 | ecoregion)" = browser_re_ecoregion
)))

# Ecoregion RE not supported — consistent with total biomass.
# All browser models fitted without RE.


# ============================================================
#  Q1 — HUMAN PRESSURE
#
#  Scientific question:
#  Does human pressure explain variation in browser biomass
#  beyond local ecological context, and which spatial metric
#  best captures SSF exploitation intensity for browsers?
#
#  Approach: AICc comparison of baseline vs baseline + each
#  pressure metric. Best metric = highest AICc weight AND
#  outperforms baseline. Coefficient direction checked
#  across all metrics regardless of support.
# ============================================================

b_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data)

b_q1_settgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data)

b_q1_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data)

b_q1_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_pop_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data)

cat("\n--- Q1: Pressure metric comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                      = b_baseline,
  "Baseline + settlement gravity" = b_q1_settgrav,
  "Baseline + market gravity"     = b_q1_mktgrav,
  "Baseline + settlement pop."    = b_q1_settpop
)))

# Results:
#   Baseline:           AICc = 835.89, weight = 0.459 (BEST)
#   Settlement gravity: DAICc = 1.61,  weight = 0.205
#   Market gravity:     DAICc = 1.65,  weight = 0.201
#   Settlement pop.:    DAICc = 2.44,  weight = 0.135
#
#   Baseline best supported — no pressure metric outperforms
#   ecological context alone for browsers.
#   Contrast with total biomass (settlement gravity
#   weight = 0.826) — browser biomass is not structured
#   by the exploitation pressure gradient in this system.
#   No pressure term carried forward into Q2 and Q3.

# ── Best Q1 model ─────────────────────────────────────────────
# Baseline best supported — no pressure term carried forward.
# All three pressure metrics within DAICc < 2.5 of baseline
# but none outperform it. Coefficient directions and
# robustness of Q2-Q3 conclusions checked in sensitivity
# analyses (a) below.
b_best_q1 <- b_baseline

# ============================================================
#  Q2 — LARVAL CONNECTIVITY
#
#  Scientific question:
#  Does larval connectivity explain additional variation
#  in browser biomass beyond the ecological baseline?
#
#  Connectivity × pressure interaction not evaluated —
#  pressure not supported in Q1. Testing an interaction
#  without a supported main effect would be
#  overparameterised and uninterpretable.
#
#  Criterion: AICc weight for model support; coefficients
#  and fold differences reported if supported.
# ============================================================

b_q2_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data)

cat("\n--- Q2: Connectivity comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                = b_best_q1,
  "Baseline + connectivity" = b_q2_conn
)))

# Results:
#   Baseline + connectivity: AICc = 831.89, weight = 0.881
#   Baseline:                DAICc = 4.00,  weight = 0.119
#   Connectivity strongly supported (weight = 0.881).

# ── Best Q2 model ─────────────────────────────────────────────
b_best_q2 <- b_q2_conn

# ── Best model summary — reported in results ─────────────────
# Full coefficient reporting for the best model identified
# across Q1-Q2. This is the primary reported model.
# Q3 is a separate question and does not update this model.

cat("\n--- Best model: full summary ---\n")
summary(b_best_q2)

# ── DHARMa diagnostics — best model ───────────────────────
cat("\n--- Q3: DHARMa diagnostics ---\n")
b_best_q2_sim <- simulateResiduals(b_best_q2, n = 1000)
plot(b_best_q2_sim)
testOutliers(b_best_q2_sim)

# DHARMa diagnostics — best Q2 model (baseline + connectivity)
# simulateResiduals n = 1000:
#   KS test:        p = 0.502 — no significant deviation
#   Dispersion:     p = 0.406 — no significant deviation
#   Outlier test:   p = 0.102 — 1 outlier / 54 obs, ns
#   Residuals vs predicted: no significant problems detected
#   Overall: good fit.

cat("\n--- Best model: 95% CIs ---\n")
print(confint(b_best_q2))

# ── Effect sizes: fold differences across observed range ──────
b_rug_b    <- fixef(b_best_q2)$cond["rugosity_sc"]
b_conn_b   <- fixef(b_best_q2)$cond["connectivity_sc"]

rug_span_b  <- diff(range(browser_model_data$rugosity_sc,
                          na.rm = TRUE))
conn_span_b <- diff(range(browser_model_data$connectivity_sc,
                          na.rm = TRUE))

cat(sprintf("\nRugosity span: %.3f SD units\n", rug_span_b))
cat(sprintf("Fold difference (low vs high rugosity): %.2fx\n",
            exp(abs(b_rug_b * rug_span_b))))

cat(sprintf("\nConnectivity span: %.3f SD units\n", conn_span_b))
cat(sprintf("Fold difference (low vs high connectivity): %.2fx\n",
            exp(abs(b_conn_b * conn_span_b))))

# ── CI-based fold difference for connectivity ─────────────────
ci_conn_b <- confint(b_best_q2)["connectivity_sc",
                                c("2.5 %", "97.5 %")]
cat(sprintf("Connectivity fold difference 95%% CI: %.2fx to %.2fx\n",
            exp(abs(ci_conn_b[1] * conn_span_b)),
            exp(abs(ci_conn_b[2] * conn_span_b))))

# Results:
#   Rugosity:     b = +0.478, p = 0.013 *
#     95% CI [update], fold difference = [update]x
#   Chla:         b = -0.192, p = 0.322 ns — retained as baseline
#   Connectivity: b = +0.421, p = 0.008 **
#     95% CI [update]
#     Span = 2.667 SD units
#     Fold difference: ~3.07x (95% CI: ~1.33x to ~7.09x)


# ============================================================
#  Q3 — FORMAL PROTECTION
#
#  Scientific question:
#  Does MPA status explain additional variation in browser
#  biomass beyond the best-supported connectivity model?
#
#  Q3 is a separate question — does not update best model.
#  MPA tested last — non-randomly placed with respect to
#  pressure and connectivity (see mpa_placement_checks.R).
#
#  Criterion: AICc weight for model support.
#  p-values for coefficient direction if supported.
# ============================================================

b_q3_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = browser_model_data)

cat("\n--- Q3: MPA comparison ---\n")
print(make_aicc_df(list(
  "Best Q2"       = b_best_q2,
  "Best Q2 + MPA" = b_q3_mpa
)))

cat("\n--- Q3: MPA coefficients ---\n")
summary(b_q3_mpa)

cat("\n--- Q3: MPA 95% CIs ---\n")
print(confint(b_q3_mpa))

# ── Medium MPA fold difference ────────────────────────────────
b_med_mpa  <- fixef(b_q3_mpa)$cond["mpa_statusmedium"]
ci_mpa     <- confint(b_q3_mpa)["mpa_statusmedium",
                                c("2.5 %", "97.5 %")]
cat(sprintf("\nMedium MPA fold difference: %.2fx\n",
            exp(b_med_mpa)))
cat(sprintf("95%% CI: %.2fx to %.2fx\n",
            exp(ci_mpa[1]), exp(ci_mpa[2])))

# ── Raw biomass by MPA — sanity check ────────────────────────
# Model fold difference should be broadly consistent with
# raw means pattern.
cat("\n--- Q3: Raw biomass by MPA status ---\n")
browser_model_data %>%
  group_by(mpa_status) %>%
  summarise(
    n              = n(),
    mean_biomass   = round(mean(mean_biomass), 1),
    median_biomass = round(median(mean_biomass), 1),
    sd_biomass     = round(sd(mean_biomass), 1),
    .groups = "drop"
  ) %>% print()

# ── DHARMa diagnostics — best Q3 model ───────────────────────
cat("\n--- Q3: DHARMa diagnostics ---\n")
b_q3_sim <- simulateResiduals(b_q3_mpa, n = 1000)
plot(b_q3_sim)
testOutliers(b_q3_sim)

# DHARMa diagnostics — best Q3 model (baseline + connectivity + MPA)
# simulateResiduals n = 1000:
#   KS test:        p = 0.894 — no significant deviation
#   Dispersion:     p = 0.406 — no significant deviation  
#   Outlier test:   p = 1.000 — 0 outliers / 54 obs
#   Residuals vs predicted: quantile deviations detected
#     (red curve, combined adjusted quantile test significant)
#     — downward trend in lower quantile across fitted range.
#     Suggests mild heteroscedasticity at low predicted values.
#     KS and dispersion tests both pass — overall fit acceptable
#     but residual structure worth noting as a limitation.
#   Outlier resolved relative to Q2 model — consistent with
#     MPA explaining the previously underfitted high-biomass site.

# Results:
#   Best Q2 + MPA: AICc = 826.43, weight = 0.939
#   Best Q2:       DAICc = 5.47,  weight = 0.061
#   MPA strongly supported.
#
#   Medium MPA: b = +1.143, p = 0.002 **
#     3.13x higher biomass (95% CI: 1.50x to 6.53x)
#   Low MPA:    b =  0.000, p = 0.999 ns — no effect
#   Rugosity:   b = +0.464, p = 0.004 ** — stable
#   Connectivity: b = +0.208, p = 0.251 — weakens once
#     MPA included; shared variance between connectivity
#     and MPA placement (see mpa_placement_checks.R)
#
#   Raw means: medium = 2461g vs none = 679g (3.6x raw)
#   Model estimate (3.1x) slightly lower — consistent
#   once connectivity context controlled.

# ── Best Q3 model ─────────────────────────────────────────────
# Q3 is a separate question — best Q2 model unchanged.
b_best_q3 <- b_q3_mpa

# ── Predicted vs observed ─────────────────────────────────────
pred_b <- predict(b_best_q3, type = "response")
cat(sprintf("\nPearson r (predicted vs observed): %.3f\n",
            cor(pred_b, browser_model_data$mean_biomass)))


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

browser_model_data_coords <- browser_model_data %>%
  left_join(site_coords, by = "site")

coords_mat_b <- cbind(browser_model_data_coords$lon,
                      browser_model_data_coords$lat)
listw5_b <- nb2listw(knn2nb(knearneigh(coords_mat_b, k = 5)),
                     style = "W")

cat("\n--- Spatial autocorrelation: browser best model ---\n")
print(moran.test(residuals(b_best_q3, type = "pearson"),
                 listw5_b))

# Moran's I = -0.078, p = 0.811 — no significant
# spatial autocorrelation in residuals.
# Contrasts with total biomass (I = 0.140, p = 0.015).
# Rugosity, connectivity, and MPA status adequately
# capture spatial variation in browser biomass.
# No spatial error modelling required.


# ============================================================
#  SENSITIVITY ANALYSES
#
#  (a) Alternative pressure metrics
#      Primary analysis uses baseline (no pressure) as Q1
#      reference. Sensitivity tests Q2 and Q3 conclusions
#      with settlement gravity and market gravity included,
#      confirming results robust to Q1 metric uncertainty.
#
#  (b) Transect-level replication
#      Mirrors Q1-Q3 at transect level using Tweedie GLMM
#      with (1 | site) random intercept.
# ============================================================

# ── (a) Settlement gravity ────────────────────────────────────
cat("\n--- Sensitivity (a): settlement gravity ---\n")

b_sens_sg_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc + connectivity_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data)

cat("\nQ2 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + SG"        = b_q1_settgrav,
  "Baseline + SG + conn" = b_sens_sg_conn
)))

cat("\nQ2 — connectivity coefficients:\n")
print(summary(b_sens_sg_conn)$coefficients$cond)
print(confint(b_sens_sg_conn))

b_sens_sg_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc + connectivity_sc + mpa_status,
  family = tweedie(link = "log"),
  data   = browser_model_data)

cat("\nQ3 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + SG + conn"       = b_sens_sg_conn,
  "Baseline + SG + conn + MPA" = b_sens_sg_mpa
)))

cat("\nQ3 — MPA coefficients:\n")
print(summary(b_sens_sg_mpa)$coefficients$cond)
print(confint(b_sens_sg_mpa))

# ── (a) Settlement gravity results ───────────────────────────
# Q2 — model comparison:
#   Baseline + SG + conn: AICc = 832.03, weight = 0.939
#   Baseline + SG:        DAICc = 5.48,  weight = 0.061
#   Connectivity strongly supported with SG included.
#
# Q2 — connectivity coefficients:
#   Connectivity: b = +0.473, z = 2.951, p = 0.003 **
#     95% CI [0.159, 0.788] — robust
#   Settlement gravity: b = -0.349, p = 0.108 ns — not significant
#
# Q3 — model comparison:
#   Baseline + SG + conn + MPA: AICc = 829.31, weight = 0.795
#   Baseline + SG + conn:       DAICc = 2.71,  weight = 0.205
#   MPA supported with SG included.
#
# Q3 — MPA coefficients:
#   Medium MPA: b = +1.134, z = 2.779, p = 0.005 **
#     95% CI [0.334, 1.934] — robust
#   Low MPA:    b = +0.003, p = 0.996 ns
#   Settlement gravity: b = -0.012, p = 0.957 ns


# ── (a) Market gravity ────────────────────────────────────────
cat("\n--- Sensitivity (a): market gravity ---\n")

b_sens_mg_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc + connectivity_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data)

cat("\nQ2 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + MG"        = b_q1_mktgrav,
  "Baseline + MG + conn" = b_sens_mg_conn
)))

cat("\nQ2 — connectivity coefficients:\n")
print(summary(b_sens_mg_conn)$coefficients$cond)
print(confint(b_sens_mg_conn))

b_sens_mg_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc + connectivity_sc + mpa_status,
  family = tweedie(link = "log"),
  data   = browser_model_data)

cat("\nQ3 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + MG + conn"       = b_sens_mg_conn,
  "Baseline + MG + conn + MPA" = b_sens_mg_mpa
)))

cat("\nQ3 — MPA coefficients:\n")
print(summary(b_sens_mg_mpa)$coefficients$cond)
print(confint(b_sens_mg_mpa))

# ── (a) Market gravity results ────────────────────────────────
# Q2 — model comparison:
#   Baseline + MG + conn: AICc = 834.54, weight = 0.818
#   Baseline + MG:        DAICc = 3.01,  weight = 0.182
#   Connectivity supported with MG included.
#
# Q2 — connectivity coefficients:
#   Connectivity: b = +0.418, z = 2.464, p = 0.014 *
#     95% CI [0.085, 0.750] — robust
#   Market gravity: b = +0.011, p = 0.956 ns — not significant
#
# Q3 — model comparison:
#   Baseline + MG + conn + MPA: AICc = 829.16, weight = 0.936
#   Baseline + MG + conn:       DAICc = 5.38,  weight = 0.064
#   MPA strongly supported with MG included.
#
# Q3 — MPA coefficients:
#   Medium MPA: b = +1.154, z = 3.085, p = 0.002 **
#     95% CI [0.421, 1.887] — robust
#   Low MPA:    b = -0.005, p = 0.993 ns
#   Market gravity: b = +0.073, p = 0.690 ns
#
# Sensitivity (a) conclusion:
#   Connectivity and medium MPA effects robust across both
#   alternative pressure metrics. Neither SG nor MG
#   significant in any model — consistent with Q1 finding
#   that pressure does not structure browser biomass.


# ── (b) Transect-level replication ───────────────────────────
# 43% zeros at transect level — Tweedie required.
# ZI Tweedie tested against standard Tweedie on baseline.

b_trans_tw_base <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc + log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data)

b_trans_tw_zi_base <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc + log_chla_sc +
    (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = browser_transect_data)

b_trans_res    <- simulateResiduals(b_trans_tw_base,    n = 1000)
b_trans_res_zi <- simulateResiduals(b_trans_tw_zi_base, n = 1000)

plot(b_trans_res);    testZeroInflation(b_trans_res)
plot(b_trans_res_zi); testZeroInflation(b_trans_res_zi)

cat("\n--- Sensitivity (b): transect family selection ---\n")
print(make_aicc_df(list(
  "Tweedie"    = b_trans_tw_base,
  "ZI Tweedie" = b_trans_tw_zi_base
)))

# ── (b) Transect family selection results ────────────────────
# ZI Tweedie: AICc = 2652.49, weight = 0.725
# Tweedie:    DAICc = 1.93,   weight = 0.276
# ZI Tweedie marginally preferred (DAICc = 1.93) but
# zero inflation test p = 0.666 (standard Tweedie) and
# p = 0.760 (ZI Tweedie) — neither significant.
# Standard Tweedie retained for consistency with
# site-level family.

# ── Transect Q1-Q3 sequence ───────────────────────────────────
b_trans_null <- glmmTMB(
  transect_browser_biomass ~ 1 + (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data)

b_trans_baseline <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc + log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data)

b_trans_pressure <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data)

b_trans_conn <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data)

b_trans_mpa <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc + mpa_status + (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data)

cat("\n--- Sensitivity (b): transect model comparison ---\n")
print(make_aicc_df(list(
  "Null"                = b_trans_null,
  "Baseline"            = b_trans_baseline,
  "Baseline + pressure" = b_trans_pressure,
  "Baseline + conn"     = b_trans_conn,
  "Best + MPA"          = b_trans_mpa
)))

# ── (b) Transect model comparison results ────────────────────
# Best + MPA:          AICc = 2653.34, weight = 0.403
# Baseline:            DAICc = 1.08,   weight = 0.235
# Baseline + conn:     DAICc = 1.17,   weight = 0.225
# Baseline + pressure: DAICc = 2.25,   weight = 0.131
# Null:                DAICc = 8.36,   weight = 0.006
#
# MPA model best supported — consistent with site-level
# Q3 (weight = 0.939). Weaker support at transect level
# (0.403 vs 0.939) — greater within-site variance reduces
# discriminatory power. Baseline and connectivity models
# competitive (DAICc < 2) — model selection less certain
# at transect level than site level, but MPA consistently
# top-ranked. Primary Q3 conclusion robust to aggregation.

cat("\n--- Sensitivity (b): best transect model coefficients ---\n")
summary(b_trans_mpa)

# ── Site random intercept SD ──────────────────────────────────
vc_b      <- VarCorr(b_trans_baseline)
site_sd_b <- sqrt(as.numeric(vc_b$cond$site))
cat(sprintf("\nSite random intercept SD = %.3f\n", site_sd_b))

# ── (b) Best transect model coefficients ─────────────────────
# Medium MPA:     b = +0.949, z = 2.330, p = 0.020 *
#   Consistent direction with site-level (b = +1.143)
#   though attenuated — expected given greater within-site
#   variance at transect level diluting the signal.
# Rugosity:       b = +0.548, z = 3.228, p = 0.001 ** — stable
# Connectivity:   b = +0.127, z = 0.629, p = 0.529 ns
#   Weaker than site level — absorbed by site random intercept
# Chla:           b = -0.353, z = -1.867, p = 0.062 . — marginal
# Low MPA:        b = +0.187, p = 0.753 ns
# Site random intercept SD = 1.069 — substantial between-site
#   variance; (1 | site) random intercept well justified.

# ============================================================
#  RESULTS SUMMARY
# ============================================================

cat("\n--- Browser results summary ---\n")
tribble(
  ~Question,   ~Result,          ~Key_finding,
  "Q1",        "Baseline best",  "weight = 0.459; no metric outperforms baseline; all metrics ns, inconsistent directions",
  "Q1 rug",    "Significant",    "b = +0.479, p = 0.006; fold diff = 11.81x (span = 5.151 SD)",
  "Q2 conn",   "Supported",      "weight = 0.881, DAICc = 4.00; b = +0.421, p = 0.008; fold diff = 3.07x (CI: 1.33x-7.09x)",
  "Q3 MPA",    "Supported",      "weight = 0.939, DAICc = 5.47; medium b = +1.143, p = 0.002; 3.13x (CI: 1.50x-6.53x); low ns",
  "Moran",     "No SAC",         "I = -0.078, p = 0.811 — no spatial autocorrelation",
  "Pearson r", "0.482",          "predicted vs observed correlation — moderate fit",
  "Sens (a)",  "Consistent",     "SG: conn b = +0.473 p = 0.003, MPA b = +1.134 p = 0.005; MG: conn b = +0.418 p = 0.014, MPA b = +1.154 p = 0.002",
  "Sens (b)",  "Consistent",     "MPA top-ranked (weight = 0.403); medium MPA b = +0.949 p = 0.020; site SD = 1.069"
) %>% print()


# ── End of script ─────────────────────────────────────────────