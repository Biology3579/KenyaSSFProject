# ============================================================
#  DRIVERS OF LARGE EXCAVATOR BIOMASS
#  Functional Group Analysis: Large Excavators/Bioeroders
# ============================================================
#
#  Scientific questions:
#
#  Q1 — Human pressure
#       Does human pressure explain variation in large
#       excavator biomass beyond local ecological context,
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
#       Identical to all previous analyses — pressure first,
#       connectivity second (extends Warmuth et al. 2024),
#       MPA last as governance response downstream of both.
#
#  Baseline model (fixed a priori, never tested):
#       biomass ~ rugosity_sc + log_chla_sc
#
#  Key difference from other groups:
#       ~13% zeros at site level — Tweedie GLM used.
#       ~58% zeros at transect level — ZI Tweedie tested
#       in sensitivity (b).
#
#  Sensitivity analyses:
#       (a) Alternative pressure metrics — only metrics
#           within DAICc < 2 of best Q1 model evaluated
#       (b) Transect-level replication (Tweedie GLMM)
# ============================================================

source(here::here("predictor_preparation.R"))


# ============================================================
#  RESPONSE VARIABLE — DATA AGGREGATION AND CHECKS
# ============================================================

# ── Transect-level aggregation ────────────────────────────────
excavator_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_excavator_biomass = sum(
      ifelse(trophic_group == "large_excavators", tot_wt_g, 0),
      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(site = as.factor(site))

cat("Excavator transects:", nrow(excavator_transects), "\n")
cat("Sites:",               n_distinct(excavator_transects$site), "\n")

# ── Transect-level dataset ────────────────────────────────────
excav_transect_data <- excavator_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(excav_transect_data$transect_excavator_biomass == 0),
    "/", nrow(excav_transect_data),
    sprintf("(%.3f)\n",
            mean(excav_transect_data$transect_excavator_biomass == 0)))

# ── Site-level dataset ────────────────────────────────────────
excavator_model_data <- excavator_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    mean_biomass           = mean(transect_excavator_biomass,
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

cat("\nExcavator model data:", nrow(excavator_model_data), "sites\n")

# ── Data checks ───────────────────────────────────────────────
excavator_model_data %>%
  dplyr::select(site, rugosity_sc, log_chla_sc,
                log_settlement_grav_sc, connectivity_sc,
                mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  { if (nrow(.) > 0) { warning("NAs in predictors:"); print(.) }
    else cat("NA check passed.\n") }

cat("\nZeros in mean_biomass:",
    sum(excavator_model_data$mean_biomass == 0), "\n")
cat("Site-level zero proportion:",
    round(mean(excavator_model_data$mean_biomass == 0), 3), "\n")
cat("\nResponse summary:\n")
print(summary(excavator_model_data$mean_biomass))
cat("\nMPA status counts:\n")
print(table(excavator_model_data$mpa_status))


# ============================================================
#  MODEL FAMILY SELECTION
#
#  ~13% zeros at site level — Tweedie handles zeros
#  natively without offset. Gaussian log rejected for
#  same reasons as browsers (zero sites drive bimodal
#  log distribution, no offset resolves this).
# ============================================================

excav_lm_baseline <- lm(
  log(mean_biomass + 0.01) ~ rugosity_sc + log_chla_sc,
  data = excavator_model_data)

par(mfrow = c(2, 2))
plot(excav_lm_baseline, main = "Gaussian log — baseline")
par(mfrow = c(1, 1))

excav_tw_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data)

excav_tw_res <- simulateResiduals(excav_tw_baseline, n = 1000)
plot(excav_tw_res)
testZeroInflation(excav_tw_res)
testDispersion(excav_tw_res)

# Family selection — Tweedie DHARMa diagnostics (baseline):
#   KS test:        p = 0.251 — no significant deviation
#   Dispersion:     p = 0.528, ratio = 1.104 — acceptable
#   Zero inflation: p = 0.942, ratio = 0.893 — not needed
#   Outlier test:   p = 0.102 — ns
#   Residuals vs predicted: quantile deviations detected
#     (red curve, lower quantile dips at mid-range predictions)
#     Combined adjusted quantile test n.s. — acceptable overall.
#   Proceed: Tweedie selected.


# ============================================================
#  RANDOM EFFECT STRUCTURE
#
#  Tested on baseline model — same rationale as all
#  previous analyses. Tweedie family throughout.
# ============================================================

excav_re_null <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data)

excav_re_ecoregion <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + (1 | ecoregion),
  family = tweedie(link = "log"),
  data   = excavator_model_data)

cat("\n--- Random effect structure: ecoregion ---\n")
print(make_aicc_df(list(
  "No RE"           = excav_re_null,
  "(1 | ecoregion)" = excav_re_ecoregion
)))

# Ecoregion RE not supported (DAICc = 2.54, weight = 0.220)
# — consistent with corallivores (DAICc = 2.54), total
# biomass (DAICc = 2.25), and browsers.
# All excavator models fitted without RE throughout.


# ============================================================
#  Q1 — HUMAN PRESSURE
#
#  Scientific question:
#  Does human pressure explain variation in large excavator
#  biomass beyond local ecological context, and which
#  spatial metric best captures SSF exploitation intensity?
#
#  Approach: AICc comparison of baseline vs baseline + each
#  pressure metric. Best metric = highest AICc weight AND
#  outperforms baseline. Coefficient directions checked
#  across all metrics regardless of support.
# ============================================================

e_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data)

e_q1_settgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data)

e_q1_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data)

e_q1_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_pop_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data)

cat("\n--- Q1: Pressure metric comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                      = e_baseline,
  "Baseline + settlement gravity" = e_q1_settgrav,
  "Baseline + market gravity"     = e_q1_mktgrav,
  "Baseline + settlement pop."    = e_q1_settpop
)))

# Results:
#   Baseline:           AICc = 686.05, weight = 0.494 (BEST)
#   Settlement pop.:    DAICc = 1.85,  weight = 0.196
#   Settlement gravity: DAICc = 2.13,  weight = 0.171
#   Market gravity:     DAICc = 2.53,  weight = 0.139
#
#   Baseline best supported — no pressure metric outperforms
#   ecological context. Inconsistent directions across
#   metrics (settlement gravity and market gravity negative,
#   settlement population positive) — all non-significant.
#   No coherent pressure signal for excavators.
#   No pressure term carried forward.
#
#   Settlement pop. (DAICc = 1.85) within threshold —
#   evaluated in sensitivity (a).
#   Settlement gravity (DAICc = 2.13) and market gravity
#   (DAICc = 2.53) outside threshold — not evaluated.

# ── Best Q1 model ─────────────────────────────────────────────
# Baseline best supported — no pressure term carried forward.
# Robustness of Q2-Q3 conclusions checked in sensitivity (a).
e_best_q1 <- e_baseline


# ============================================================
#  Q2 — LARVAL CONNECTIVITY
#
#  Scientific question:
#  Does larval connectivity explain additional variation
#  in large excavator biomass beyond the ecological baseline?
#
#  Connectivity × pressure interaction not evaluated —
#  pressure not supported in Q1.
#
#  Criterion: AICc weight for model support; coefficients
#  and fold differences reported if supported.
# ============================================================

e_q2_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data)

cat("\n--- Q2: Connectivity comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                = e_best_q1,
  "Baseline + connectivity" = e_q2_conn
)))

# ── Best Q2 model ─────────────────────────────────────────────
# Connectivity not supported — baseline retained.
e_best_q2 <- e_best_q1

# ── Best model summary — reported in results ─────────────────
# Full coefficient reporting for best model across Q1-Q2.
# Q3 is a separate question and does not update this model.

cat("\n--- Best model: full summary ---\n")
summary(e_best_q2)

cat("\n--- Best model: 95% CIs ---\n")
print(confint(e_best_q2))

# ── Effect sizes: fold differences across observed range ──────
b_rug_e <- fixef(e_best_q2)$cond["rugosity_sc"]
rug_span_e <- diff(range(excavator_model_data$rugosity_sc,
                         na.rm = TRUE))

cat(sprintf("\nRugosity span: %.3f SD units\n", rug_span_e))
cat(sprintf("Fold difference (low vs high rugosity): %.2fx\n",
            exp(abs(b_rug_e * rug_span_e))))

ci_rug_e <- confint(e_best_q2)["rugosity_sc", c("2.5 %", "97.5 %")]
cat(sprintf("Rugosity fold difference 95%% CI: %.2fx to %.2fx\n",
            exp(abs(ci_rug_e[1] * rug_span_e)),
            exp(abs(ci_rug_e[2] * rug_span_e))))

# ── DHARMa diagnostics — best Q2 model ───────────────────────
cat("\n--- DHARMa diagnostics: best Q2 model ---\n")
e_best_q2_sim <- simulateResiduals(e_best_q2, n = 1000)
plot(e_best_q2_sim)
testOutliers(e_best_q2_sim)

# ── Q2 results ────────────────────────────────────────────────
# Results:
#   Baseline:                DAICc = 0.00, weight = 0.758
#   Baseline + connectivity: DAICc = 2.28, weight = 0.242
#   Connectivity not supported.
#
#   Best model: rugosity + chla (baseline only)
#   Rugosity: b = +0.506, z = 2.625, p = 0.009 **
#     95% CI [0.128, 0.884]
#     Fold difference: 13.54x across observed range
#     (span = 5.151 SD units)
#     95% CI on fold difference: 1.93x to 94.76x
#     Wide CI reflects high within-group variance
#     (dispersion = 8.51).
#   Chla: b = -0.300, z = -1.414, p = 0.157 ns
#     95% CI [-0.715, 0.116]
#
# DHARMa diagnostics — best Q2 model (baseline):
#   KS test:        p = 0.161 — no significant deviation
#   Dispersion:     p = 0.604 — no significant deviation
#   Outlier test:   p = 1.000 — 0 outliers / 54 obs
#     (1 outlier at n = 1000 sims previously — resolved)
#   Residuals vs predicted: quantile deviations detected
#     (red curve, lower quantile below expected at high
#     fitted values). Combined adjusted quantile test
#     significant — mild heteroscedasticity at high
#     predicted values. Likely driven by near-zero and
#     zero sites at low predicted biomass end.
#     KS and dispersion pass — acceptable overall but
#     note as a limitation.
#   Newton convergence warning noted — quantile smoothing
#     only; does not affect model estimates or tests.

# ============================================================
#  Q3 — FORMAL PROTECTION
#
#  Scientific question:
#  Does MPA status explain additional variation in large
#  excavator biomass beyond the best Q2 model?
#
#  Q3 is a separate question — does not update best model.
#  MPA tested last — non-randomly placed with respect to
#  pressure and connectivity (see mpa_placement_checks.R).
#
#  Criterion: AICc weight for model support.
# ============================================================

e_q3_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = excavator_model_data)

cat("\n--- Q3: MPA comparison ---\n")
print(make_aicc_df(list(
  "Best Q2"       = e_best_q2,
  "Best Q2 + MPA" = e_q3_mpa
)))

cat("\n--- Q3: MPA coefficients ---\n")
summary(e_q3_mpa)

cat("\n--- Q3: MPA 95% CIs ---\n")
print(confint(e_q3_mpa))

# ── Raw biomass by MPA — sanity check ────────────────────────
cat("\n--- Q3: Raw biomass by MPA status ---\n")
excavator_model_data %>%
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
e_q3_sim <- simulateResiduals(e_q3_mpa, n = 1000)
plot(e_q3_sim)
testOutliers(e_q3_sim)

# Results:
#   Best Q2:       AICc = 686.05, weight = 0.648
#   Best Q2 + MPA: DAICc =  1.22, weight = 0.353
#   Genuine model selection uncertainty — best Q2 preferred
#   but MPA model competitive.
#
#   Low MPA:    b = -1.760, p = 0.029 * — significant but
#     NEGATIVE and artefactual. Low MPA sites confined to
#     narrow high-connectivity range (z = 0.66 to max, n = 7)
#     with very low mean biomass (71g vs 621g unprotected).
#     Negative coefficient reflects site characteristics,
#     not a genuine negative protection effect.
#   Medium MPA: b = -0.437, p = 0.427 ns
#     Raw mean (640g) essentially identical to unprotected
#     sites (621g) — no protection signal.
#   Rugosity: b = +0.516, p = 0.014 * — stable.
#
#   Conclusion: MPA not supported (DAICc = 1.22 — within
#   uncertainty threshold). Low MPA artefact consistent
#   with corallivores. Baseline model retained.
#   
#   DHARMa diagnostics — Q3 model (baseline + MPA):
#   KS test:        p = 0.161 — no significant deviation
#   Dispersion:     p = 0.604 — no significant deviation
#   Outlier test:   p = 1.000 — 0 outliers / 54 obs
#   Residuals vs predicted: no significant problems detected
#   Improvement over Q2 model — quantile deviation
#     resolved once MPA categories added. Consistent with
#     MPA absorbing some of the low-biomass site structure
#     (low MPA sites have very low biomass, mean = 71g).
#     However MPA not supported by AICc — this improvement
#     reflects site characteristics not genuine protection.

# ── Best Q3 model ─────────────────────────────────────────────
# Q3 is a separate question — best Q2 model unchanged.
e_best_q3 <- e_best_q2

# ── Predicted vs observed ─────────────────────────────────────
pred_e <- predict(e_best_q3, type = "response")
cat(sprintf("\nPearson r (predicted vs observed): %.3f\n",
            cor(pred_e, excavator_model_data$mean_biomass)))

# ── Flagged sites ─────────────────────────────────────────────
# Update row indices from diagnostic plots after running
flagged_rows_e <- c()  # update from plots

if (length(flagged_rows_e) > 0) {
  excavator_model_data %>%
    slice(flagged_rows_e) %>%
    dplyr::select(site, mpa_status, ecoregion,
                  mean_biomass, rugosity_sc,
                  log_settlement_grav_sc, connectivity_sc) %>%
    print()
}


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

excav_model_data_coords <- excavator_model_data %>%
  left_join(site_coords, by = "site")

coords_mat_e <- cbind(excav_model_data_coords$lon,
                      excav_model_data_coords$lat)
listw5_e <- nb2listw(knn2nb(knearneigh(coords_mat_e, k = 5)),
                     style = "W")

# Warning: knearneigh — identical points found, kd_tree not
# available. 4 sub-graphs expected — discontinuous sampling
# design across four countries. See total_biomass.R.

cat("\n--- Spatial autocorrelation: excavator best model ---\n")
print(moran.test(residuals(e_best_q3, type = "pearson"),
                 listw5_e))

# Moran's I = -0.015, p = 0.479 — no significant spatial
# autocorrelation in residuals.
# Consistent with browsers (I = -0.078, p = 0.811) and
# corallivores (I = -0.021, p = 0.509).
# Contrasts with total biomass (I = 0.140, p = 0.015)
# and grazer-detritivores (I = 0.210, p = 0.001).
# Rugosity adequately captures spatial variation in
# excavator biomass without residual geographic signal.
# No spatial error modelling required.


# ============================================================
#  SENSITIVITY ANALYSES
#
#  (a) Alternative pressure metrics
#      Only metrics within DAICc < 2 of best Q1 model
#      are evaluated. Settlement pop. (DAICc = 1.85)
#      meets this threshold. Settlement gravity
#      (DAICc = 2.13) and market gravity (DAICc = 2.53)
#      do not — not evaluated.
#
#  (b) Transect-level replication
#      ~58% zeros at transect level — ZI Tweedie tested
#      against standard Tweedie on baseline structure.
# ============================================================

# ── (a) Settlement population ─────────────────────────────────
cat("\n--- Sensitivity (a): settlement population ---\n")

e_sens_sp_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_pop_sc + connectivity_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data)

cat("\nQ2 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + SP"        = e_q1_settpop,
  "Baseline + SP + conn" = e_sens_sp_conn
)))

e_sens_sp_conn_int <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_pop_sc + connectivity_sc +
    log_settlement_pop_sc:connectivity_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data)

print(make_aicc_df(list(
  "Baseline + SP"              = e_q1_settpop,
  "Baseline + SP + conn"       = e_sens_sp_conn,
  "Baseline + SP + conn + int" = e_sens_sp_conn_int
)))

summary(e_sens_sp_conn_int)
confint(e_sens_sp_conn_int)

cat("\nQ2 — connectivity coefficients:\n")
print(summary(e_sens_sp_conn)$coefficients$cond)
print(confint(e_sens_sp_conn))

e_sens_sp_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_pop_sc + connectivity_sc + mpa_status,
  family = tweedie(link = "log"),
  data   = excavator_model_data)

cat("\nQ3 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + SP + conn"       = e_sens_sp_conn,
  "Baseline + SP + conn + MPA" = e_sens_sp_mpa
)))

cat("\nQ3 — MPA coefficients:\n")
print(summary(e_sens_sp_mpa)$coefficients$cond)
print(confint(e_sens_sp_mpa))

# ── Sensitivity (a): settlement population results ────────────
# Q2:
#   Baseline + SP:        AICc = 687.90, weight = 0.729
#   Baseline + SP + conn: DAICc = 1.98,  weight = 0.271
#   Connectivity not supported (DAICc = 1.98) — borderline
#   but consistent with primary Q2 (DAICc = 2.28).
#   Settlement pop.: b = +0.359, p = 0.283 ns
#   Connectivity:    b = -0.224, p = 0.412 ns
#
# Q3:
#   Baseline + SP + conn:       AICc = 689.87, weight = 0.724
#   Baseline + SP + conn + MPA: DAICc = 1.92,  weight = 0.277
#   MPA not clearly supported (DAICc = 1.92 — borderline).
#   Low MPA:    b = -1.908, p = 0.033 * — artefactual,
#     95% CI [-3.660, -0.156] — same pattern as primary Q3.
#   Medium MPA: b = -0.377, p = 0.603 ns
#   Settlement pop.: b = +0.328, p = 0.359 ns
#   Connectivity:    b = +0.014, p = 0.970 ns — disappears
#
# Conclusion: connectivity and MPA conclusions consistent
#   with primary analysis. Rugosity remains sole robust
#   predictor (b = +0.602–0.643, p < 0.01 across models).

# ── (b) Transect-level replication ───────────────────────────
# ~58% zeros at transect level — ZI Tweedie tested.

e_trans_tw_base <- glmmTMB(
  transect_excavator_biomass ~ rugosity_sc + log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = excav_transect_data)

e_trans_tw_zi_base <- glmmTMB(
  transect_excavator_biomass ~ rugosity_sc + log_chla_sc +
    (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = excav_transect_data)

e_trans_res    <- simulateResiduals(e_trans_tw_base,    n = 500)
e_trans_res_zi <- simulateResiduals(e_trans_tw_zi_base, n = 500)

plot(e_trans_res);    testZeroInflation(e_trans_res)
plot(e_trans_res_zi); testZeroInflation(e_trans_res_zi)

cat("\n--- Sensitivity (b): transect family selection ---\n")
print(make_aicc_df(list(
  "Tweedie"    = e_trans_tw_base,
  "ZI Tweedie" = e_trans_tw_zi_base
)))

# Transect family selection — DHARMa:
#   Standard Tweedie: zero inflation p = 0.724, ratio = 0.975
#   ZI Tweedie:       zero inflation p = 0.716, ratio = 0.975
#   Neither ZI test significant — standard Tweedie sufficient.
#   Newton convergence warning in both — quantile smoothing
#     only, does not affect model estimates.
#   AICc: standard Tweedie weight = 0.743 — selected.
#   

# ── Transect Q1-Q3 sequence ───────────────────────────────────
e_trans_null <- glmmTMB(
  transect_excavator_biomass ~ 1 + (1 | site),
  family = tweedie(link = "log"),
  data   = excav_transect_data)

e_trans_baseline <- glmmTMB(
  transect_excavator_biomass ~ rugosity_sc + log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = excav_transect_data)

e_trans_pressure <- glmmTMB(
  transect_excavator_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = excav_transect_data)

e_trans_conn <- glmmTMB(
  transect_excavator_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = excav_transect_data)

e_trans_mpa <- glmmTMB(
  transect_excavator_biomass ~ rugosity_sc + log_chla_sc +
    mpa_status + (1 | site),
  family = tweedie(link = "log"),
  data   = excav_transect_data)

cat("\n--- Sensitivity (b): transect model comparison ---\n")
print(make_aicc_df(list(
  "Null"                = e_trans_null,
  "Baseline"            = e_trans_baseline,
  "Baseline + pressure" = e_trans_pressure,
  "Baseline + conn"     = e_trans_conn,
  "Best + MPA"          = e_trans_mpa
)))

cat("\n--- Sensitivity (b): baseline model coefficients ---\n")
summary(e_trans_baseline)

# ── Site random intercept SD ──────────────────────────────────
vc_e      <- VarCorr(e_trans_baseline)
site_sd_e <- sqrt(as.numeric(vc_e$cond$site))
cat(sprintf("\nSite random intercept SD = %.3f\n", site_sd_e))

# Results:
#   Best + MPA:        AICc = 1855.63, weight = 0.308
#   Baseline:          DAICc = 0.28,   weight = 0.268
#   Null:              DAICc = 1.16,   weight = 0.173
#   Baseline + conn:   DAICc = 1.36,   weight = 0.156
#   Baseline + pres:   DAICc = 2.36,   weight = 0.095
#
#   Complete model selection uncertainty at transect level
#   — all five models within DAICc < 2.5. Null model
#   third-ranked (weight = 0.173) confirms excavator
#   biomass is the most weakly structured functional group.
#   Consistent with Q1-Q3 at site level.
#
#   Rugosity: b = +0.605, t = 1.953 . — direction
#     consistent with site-level (b = +0.506, p = 0.009)
#     but significance lost at transect level due to
#     extremely high within-site variance.
#   Site random intercept SD = 1.843 — extremely high,
#     reflecting patchy distribution of large-bodied
#     excavators within sites.
#     
#   Transect baseline diagnostics:
#   KS test:        p = 0.468 — excellent fit
#   Dispersion:     p = 0.636 — no significant deviation
#   Outlier test:   p = 1.000 — no outliers
#   Residuals vs predicted: no significant problems detected
#   Gaussian lmer structure confirmed appropriate.

#   Transect standard Tweedie diagnostics:
#   KS test:        p = 0.275 — no significant deviation
#   Dispersion:     p = 0.604 — no significant deviation
#   Outlier test:   p = 1.000 — no outliers
#   Residuals vs predicted: no significant problems detected
#   Good fit at transect level.


# ============================================================
#  RESULTS SUMMARY
# ============================================================

cat("\n--- Excavator results summary ---\n")
tribble(
  ~Question,   ~Result,          ~Key_finding,
  "Q1",        "Baseline best",  "weight = 0.494; all metrics ns, inconsistent directions",
  "Q1 rug", "Significant", "b = +0.506, p = 0.009; fold diff = 13.54x (CI: 1.93x–94.76x) — strongest rugosity signal",
  "Q2 conn",   "Not supported",  "DAICc = 2.28, weight = 0.242; b = -0.129, p = 0.614 ns",
  "Q3 MPA",    "Not supported",  "DAICc = 1.22, weight = 0.353 — within uncertainty; low MPA artefactual",
  "Moran",     "No SAC",         "I = -0.015, p = 0.479 — no spatial autocorrelation",
  "Pearson r", "0.420", "moderate fit",
  "Sens (a)", "Consistent", "SP: conn DAICc = 1.98 ns; MPA DAICc = 1.92 borderline; low MPA artefactual; rugosity sole driver",
  "Sens (b)", "Consistent", "complete uncertainty at transect level; rugosity b = +0.605 p = 0.051; site SD = 1.843",
) %>% print()


# ── End of script ─────────────────────────────────────────────
