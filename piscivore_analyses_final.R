# ============================================================
#  DRIVERS OF PISCIVORE BIOMASS
#  Functional Group Analysis: Piscivores
# ============================================================
#
#  Scientific questions:
#
#  Q1 — Human pressure
#       Does human pressure explain variation in piscivore
#       biomass beyond local ecological context, and which
#       spatial metric best captures SSF exploitation
#       intensity for this functional group?
#       Piscivores are large-bodied commercial species —
#       market gravity expected to outperform settlement
#       gravity, which captures subsistence rather than
#       market-oriented exploitation.
#
#  Q2 — Larval connectivity
#       Does larval connectivity explain additional variation
#       beyond the best Q1 model?
#       Pressure supported in Q1 — connectivity × pressure
#       interaction evaluated alongside main effect.
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
#  Key differences from other groups:
#       ~11% zeros at site level — Tweedie GLM used.
#       ~47% zeros at transect level — ZI Tweedie tested
#       in sensitivity (b).
#       Market gravity preferred in Q1 — first group besides
#       total biomass with decisive metric differentiation.
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
pisc_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_pisc_biomass = sum(
      ifelse(trophic_group == "piscivores", tot_wt_g, 0),
      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(site = as.factor(site))

cat("Piscivore transects:", nrow(pisc_transects), "\n")
cat("Sites:",               n_distinct(pisc_transects$site), "\n")

# ── Transect-level dataset ────────────────────────────────────
pisc_transect_data <- pisc_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(pisc_transect_data$transect_pisc_biomass == 0),
    "/", nrow(pisc_transect_data),
    sprintf("(%.3f)\n",
            mean(pisc_transect_data$transect_pisc_biomass == 0)))

# ── Site-level dataset ────────────────────────────────────────
pisc_model_data <- pisc_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    mean_biomass           = mean(transect_pisc_biomass,
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

cat("\nPiscivore model data:", nrow(pisc_model_data), "sites\n")

# ── Data checks ───────────────────────────────────────────────
pisc_model_data %>%
  dplyr::select(site, rugosity_sc, log_chla_sc,
                log_settlement_grav_sc, connectivity_sc,
                mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  { if (nrow(.) > 0) { warning("NAs in predictors:"); print(.) }
    else cat("NA check passed.\n") }

cat("\nZeros in mean_biomass:",
    sum(pisc_model_data$mean_biomass == 0), "\n")
cat("Site-level zero proportion:",
    round(mean(pisc_model_data$mean_biomass == 0), 3), "\n")
cat("\nResponse summary:\n")
print(summary(pisc_model_data$mean_biomass))
cat("\nMPA status counts:\n")
print(table(pisc_model_data$mpa_status))


# ============================================================
#  MODEL FAMILY SELECTION
#
#  ~11% zeros at site level — Tweedie handles zeros
#  natively without offset. Gaussian log rejected for
#  same reasons as browsers and excavators.
# ============================================================

pisc_lm_baseline <- lm(
  log(mean_biomass + 0.01) ~ rugosity_sc + log_chla_sc,
  data = pisc_model_data)

par(mfrow = c(2, 2))
plot(pisc_lm_baseline, main = "Gaussian log — baseline")
par(mfrow = c(1, 1))

pisc_tw_base <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data)

pisc_tw_res <- simulateResiduals(pisc_tw_base, n = 1000)
plot(pisc_tw_res)
testZeroInflation(pisc_tw_res)
testDispersion(pisc_tw_res)

# Family selection:
#   Gaussian log: REJECTED — zero sites pull residuals
#     to extreme lower tail values (sites 4, 41, 15 at
#     theoretical quantile ~-2.5, Q-Q plot). Scale-location
#     shows strong downward trend — heteroscedasticity
#     driven by zero/near-zero sites. Residuals vs Leverage:
#     site 47 (high leverage ~0.18, large positive residual),
#     sites 21 and 42 have large negative residuals at
#     moderate leverage. Bimodal log distribution from zero
#     sites not resolvable by offset.
#
#   Tweedie (log link): SELECTED
#   DHARMa diagnostics (n = 1000):
#       KS test:        p = 0.987 — excellent fit
#       Dispersion:     p = 0.634, ratio = 1.117 — acceptable
#       Zero inflation: p = 1.000, ratio = 0.970 — not needed
#       Outlier test:   p = 1.000 — no outliers
#     Residuals vs predicted: slight upward trend in
#       smoothing line but within confidence band.
#       No significant problems detected.
#
#   Proceed: glmmTMB(family = tweedie(link = "log")) on
#   raw mean_biomass throughout all piscivore analyses.


# ============================================================
#  RANDOM EFFECT STRUCTURE
#
#  Tested on baseline model — same rationale as all
#  previous analyses. Tweedie family throughout.
# ============================================================

pisc_re_null <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data)

pisc_re_ecoregion <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + (1 | ecoregion),
  family = tweedie(link = "log"),
  data   = pisc_model_data)

cat("\n--- Random effect structure: ecoregion ---\n")
print(make_aicc_df(list(
  "No RE"           = pisc_re_null,
  "(1 | ecoregion)" = pisc_re_ecoregion
)))

# Ecoregion RE not supported (DAICc = 2.80, weight = 0.198)
# — consistent with all previous functional groups.
# All piscivore models fitted without RE throughout.


# ============================================================
#  Q1 — HUMAN PRESSURE
#
#  Scientific question:
#  Does human pressure explain variation in piscivore
#  biomass beyond local ecological context, and which
#  spatial metric best captures SSF exploitation intensity?
#
#  Approach: AICc comparison of baseline vs baseline + each
#  pressure metric. Best metric = highest AICc weight AND
#  outperforms baseline. Coefficient directions checked
#  across all metrics regardless of support.
# ============================================================

p_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data)

p_q1_settgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data)

p_q1_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data)

p_q1_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_pop_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data)

cat("\n--- Q1: Pressure metric comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                      = p_baseline,
  "Baseline + settlement gravity" = p_q1_settgrav,
  "Baseline + market gravity"     = p_q1_mktgrav,
  "Baseline + settlement pop."    = p_q1_settpop
)))

# Results:
#   Market gravity:     AICc = 799.72, weight = 0.488 (BEST)
#   Baseline:           DAICc = 0.96,  weight = 0.302
#   Settlement pop.:    DAICc = 2.96,  weight = 0.111
#   Settlement gravity: DAICc = 3.19,  weight = 0.099
#
#   Market gravity preferred — weight = 0.488, nearly
#   twice the baseline weight. Genuine model selection
#   uncertainty (DAICc = 0.96 vs baseline) — not decisive.
#   Market gravity carried forward as ecologically motivated
#   (piscivores are commercially targeted large-bodied
#   species) and consistent with prior analysis.
#
#   Market gravity (DAICc = 0.96) within threshold —
#   only metric meeting DAICc < 2 criterion.
#   Settlement pop. (DAICc = 2.96) and settlement gravity
#   (DAICc = 3.19) outside threshold — not evaluated in
#   sensitivity (a).

# ── Best Q1 model ─────────────────────────────────────────────
p_best_q1 <- p_q1_mktgrav


# ============================================================
#  Q2 — LARVAL CONNECTIVITY
#
#  Scientific question:
#  Does larval connectivity explain additional variation
#  in piscivore biomass beyond the pressure baseline, and
#  does it modify the relationship between market gravity
#  and biomass?
#
#  Pressure supported in Q1 — connectivity × pressure
#  interaction evaluated alongside main effect.
#
#  Approach: three-model AICc comparison.
#  (1) Best Q1 — pressure only
#  (2) Best Q1 + connectivity main effect
#  (3) Best Q1 + connectivity + connectivity × pressure
#
#  Criterion: AICc weight for model support; coefficients
#  and fold differences reported if supported.
# ============================================================

p_q2_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data)

p_q2_conn_int <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc +
    connectivity_sc +
    log_market_gravity_sc:connectivity_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data)

cat("\n--- Q2: Connectivity main effect and interaction ---\n")
print(make_aicc_df(list(
  "Best Q1"                      = p_best_q1,
  "Best Q1 + conn"               = p_q2_conn,
  "Best Q1 + conn + interaction" = p_q2_conn_int
)))

# ── Best Q2 model ─────────────────────────────────────────────
p_best_q2 <- p_q2_conn_int

# ── Best model summary — reported in results ─────────────────
# Full coefficient reporting for best model across Q1-Q2.
# Q3 is a separate question and does not update this model.

cat("\n--- Best model: full summary ---\n")
summary(p_best_q2)

cat("\n--- Best model: 95% CIs ---\n")
print(confint(p_best_q2))

# ── Effect sizes ──────────────────────────────────────────────
mkt_span_p  <- diff(range(pisc_model_data$log_market_gravity_sc,
                          na.rm = TRUE))
conn_span_p <- diff(range(pisc_model_data$connectivity_sc,
                          na.rm = TRUE))

cat(sprintf("\nMarket gravity span: %.3f SD units\n", mkt_span_p))
cat(sprintf("Connectivity span:   %.3f SD units\n", conn_span_p))

# ── DHARMa diagnostics — best Q2 model ───────────────────────
cat("\n--- DHARMa diagnostics: best Q2 model ---\n")
p_best_q2_sim <- simulateResiduals(p_best_q2, n = 1000)
plot(p_best_q2_sim)
testOutliers(p_best_q2_sim)

# Results:
#   Best Q1 + conn + int: AICc = 795.27, weight = 0.836
#   Best Q1:              DAICc = 4.45,  weight = 0.090
#   Best Q1 + conn:       DAICc = 4.85,  weight = 0.074
#   Interaction strongly supported (weight = 0.836).
#
#   Connectivity × market gravity: b = -0.511, p = 0.004 **
#     Negative — connectivity moderates market gravity.
#     At low connectivity, market access associates
#     positively with biomass. At high connectivity,
#     market access depletes biomass — consistent with
#     connectivity facilitating commercial fishing access.
#   Market gravity: b = +0.407, p = 0.008 **
#     Positive at mean connectivity.
#   Connectivity:   b = +0.265, p = 0.048 *
#     Positive at mean market gravity.
#   Dispersion parameter = 20.5 — higher than baseline
#     (8.51), reflecting increased within-group variance
#     once the interaction structure is included.
#   Rugosity:       b = +0.240, p = 0.067 . — marginal
#   DHARMa diagnostics — best Q2 model (conn x press):
#   KS test:        p = 0.949 — excellent fit
#   Dispersion:     p = 0.822 — no significant deviation
#   Outlier test:   p = 1.000 — no outliers
#   Residuals vs predicted: no significant problems detected
#   Overall: excellent fit.

# ============================================================
#  Q3 — FORMAL PROTECTION
#
#  Scientific question:
#  Does MPA status explain additional variation in piscivore
#  biomass beyond the best Q2 model?
#
#  Q3 is a separate question — does not update best model.
#  MPA tested last — non-randomly placed with respect to
#  pressure and connectivity (see mpa_placement_checks.R).
#
#  Criterion: AICc weight for model support.
# ============================================================

p_q3_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc +
    connectivity_sc +
    log_market_gravity_sc:connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = pisc_model_data)

cat("\n--- Q3: MPA comparison ---\n")
print(make_aicc_df(list(
  "Best Q2"       = p_best_q2,
  "Best Q2 + MPA" = p_q3_mpa
)))

cat("\n--- Q3: MPA coefficients ---\n")
summary(p_q3_mpa)

cat("\n--- Q3: MPA 95% CIs ---\n")
print(confint(p_q3_mpa))

# ── Medium MPA fold difference ────────────────────────────────
b_med_mpa_p  <- fixef(p_q3_mpa)$cond["mpa_statusmedium"]
ci_mpa_p     <- confint(p_q3_mpa)["mpa_statusmedium",
                                  c("2.5 %", "97.5 %")]
cat(sprintf("\nMedium MPA fold difference: %.2fx\n",
            exp(b_med_mpa_p)))
cat(sprintf("95%% CI: %.2fx to %.2fx\n",
            exp(ci_mpa_p[1]), exp(ci_mpa_p[2])))

# ── Raw biomass by MPA — sanity check ────────────────────────
cat("\n--- Q3: Raw biomass by MPA status ---\n")
pisc_model_data %>%
  group_by(mpa_status) %>%
  summarise(
    n              = n(),
    mean_biomass   = round(mean(mean_biomass), 1),
    median_biomass = round(median(mean_biomass), 1),
    sd_biomass     = round(sd(mean_biomass), 1),
    mean_conn      = round(mean(connectivity_sc), 3),
    .groups = "drop"
  ) %>% print()

# ── DHARMa diagnostics — Q3 model ────────────────────────────
cat("\n--- DHARMa diagnostics: Q3 model ---\n")
p_q3_sim <- simulateResiduals(p_q3_mpa, n = 1000)
plot(p_q3_sim)
testOutliers(p_q3_sim)

# Results:
#   Best Q2 + MPA: AICc = 792.22, weight = 0.821
#   Best Q2:       DAICc =  3.05, weight = 0.179
#   MPA strongly supported beyond connectivity x pressure.
#
#   Medium MPA: b = +0.868, z = 2.96, p = 0.003 **
#     fold difference = 2.38x (95% CI: 1.34x to 4.23x)
#     Significant positive effect — piscivores benefit
#     directly from harvest exclusion.
#   Low MPA:    b = +0.197, p = 0.649 ns — not significant
#   Interaction: b = -0.578, z = -3.41, p = 0.001 ***
#     Strengthened and stable once MPA included
#     (was -0.511, p = 0.004).
#   Rugosity: b = +0.299, p = 0.015 * — now significant.
#
#   Contrasts with total biomass (MPA not supported,
#   DAICc = 3.01) — piscivores most sensitive to formal
#   protection of all functional groups.
#   Raw means: none = 579g, low = 754g, medium = 1284g
#   DHARMa diagnostics — Q3 model (conn x press + MPA):
#   KS test:        p = 0.869 — no significant deviation
#   Dispersion:     p = 0.788 — no significant deviation
#   Outlier test:   p = 1.000 — no outliers
#   Residuals vs predicted: quantile deviations detected
#     (red curve, lower quantile elevated at high fitted
#     values). Combined adjusted quantile test n.s. —
#     acceptable overall despite visual pattern.

# ── Best Q3 model ─────────────────────────────────────────────
# Q3 is a separate question — best Q2 model unchanged.
p_best_q3 <- p_best_q2

# ── Predicted vs observed ─────────────────────────────────────
pred_p <- predict(p_best_q3, type = "response")
cat(sprintf("\nPearson r (predicted vs observed): %.3f\n",
            cor(pred_p, pisc_model_data$mean_biomass)))

# ── Flagged sites ─────────────────────────────────────────────
# DHARMa diagnostics clean — 0 outliers, no influential
# sites identified. Gaussian log baseline flags sites 4,
# 15, 41 (zero/near-zero biomass, lower tail Q-Q) and
# sites 21, 42, 47 (leverage plot) but none exceed Cook's
# distance threshold and Tweedie DHARMa shows no problems.
# No sites warrant exclusion or special treatment.

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

pisc_model_data_coords <- pisc_model_data %>%
  left_join(site_coords, by = "site")

coords_mat_p <- cbind(pisc_model_data_coords$lon,
                      pisc_model_data_coords$lat)
listw5_p <- nb2listw(knn2nb(knearneigh(coords_mat_p, k = 5)),
                     style = "W")

# Warning: knearneigh — identical points found, kd_tree not
# available. 4 sub-graphs expected — discontinuous sampling
# design across four countries. See total_biomass.R.

cat("\n--- Spatial autocorrelation: piscivore best model ---\n")
print(moran.test(residuals(p_best_q3, type = "pearson"),
                 listw5_p))

# Moran's I = +0.052, p = 0.162 — no significant spatial
# autocorrelation in residuals.
# Consistent with corallivores (I = -0.021) and excavators
# (I = -0.015). Contrasts with total biomass (I = 0.140)
# and grazer-detritivores (I = 0.210).
# Connectivity × pressure interaction and MPA adequately
# capture spatial variation in piscivore biomass.
# No spatial error modelling required.


# ============================================================
#  SENSITIVITY ANALYSES
#
#  (a) Alternative pressure metrics
#      Only metrics within DAICc < 2 of best Q1 model
#      are evaluated. Market gravity best supported —
#      baseline is the next closest (DAICc = 0.96).
#      Settlement pop. (DAICc = 2.96) and settlement
#      gravity (DAICc = 3.19) outside threshold —
#      not evaluated.
#      Sensitivity (a) therefore evaluates whether the
#      Q2-Q3 conclusions hold when the baseline (no
#      pressure) is used as the Q1 reference instead.
#
#  (b) Transect-level replication
#      ~47% zeros at transect level — ZI Tweedie tested.
# ============================================================

# ── (a) Baseline as Q1 reference ─────────────────────────────
# Tests whether Q2-Q3 conclusions are robust when
# pressure is excluded — i.e. whether connectivity
# interaction and MPA effects hold without market gravity.
cat("\n--- Sensitivity (a): baseline as Q1 reference ---\n")

p_sens_conn_int <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data)

cat("\nQ2 — model comparison (from baseline):\n")
print(make_aicc_df(list(
  "Baseline"                = p_baseline,
  "Baseline + conn"         = p_sens_conn_int
)))

cat("\nQ2 — connectivity coefficients:\n")
print(summary(p_sens_conn_int)$coefficients$cond)
print(confint(p_sens_conn_int))

p_sens_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc + mpa_status,
  family = tweedie(link = "log"),
  data   = pisc_model_data)

cat("\nQ3 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + conn"       = p_sens_conn_int,
  "Baseline + conn + MPA" = p_sens_mpa
)))

cat("\nQ3 — MPA coefficients:\n")
print(summary(p_sens_mpa)$coefficients$cond)
print(confint(p_sens_mpa))

# Sensitivity (a) results:
# Q2:
#   Baseline + conn: AICc = 800.17, weight = 0.563
#   Baseline:        DAICc = 0.51,  weight = 0.437
#   Genuine uncertainty — connectivity marginally supported
#   without market gravity (b = +0.249, p = 0.077 .)
#   95% CI [-0.027, 0.524] — overlaps zero.
#   Direction consistent with primary Q2.
#
# Q3:
#   Baseline + conn + MPA: weight = 0.492, DAICc = 0.07
#   Baseline + conn:       weight = 0.508
#   Complete uncertainty — MPA not clearly supported
#   without market gravity in model.
#   Medium MPA: b = +0.727, p = 0.025 *
#     95% CI [0.089, 1.366] — direction consistent with
#     primary Q3 but weaker support overall.
#   Low MPA: b = +0.149, p = 0.757 ns
#
# Conclusion: primary Q2 and Q3 findings are partly
#   dependent on market gravity being in the model.
#   The connectivity × market gravity interaction is
#   the key result — without pressure, connectivity
#   main effect is marginal (p = 0.077) and MPA support
#   is uncertain (DAICc = 0.07). Market gravity selection
#   in Q1 is therefore consequential for downstream
#   conclusions.


# ── (b) Transect-level replication ───────────────────────────
# ~47% zeros at transect level — ZI Tweedie tested.

p_trans_tw_base <- glmmTMB(
  transect_pisc_biomass ~ rugosity_sc + log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = pisc_transect_data)

p_trans_tw_zi_base <- glmmTMB(
  transect_pisc_biomass ~ rugosity_sc + log_chla_sc +
    (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = pisc_transect_data)

p_trans_res    <- simulateResiduals(p_trans_tw_base,    n = 500)
p_trans_res_zi <- simulateResiduals(p_trans_tw_zi_base, n = 500)

plot(p_trans_res);    testZeroInflation(p_trans_res)
plot(p_trans_res_zi); testZeroInflation(p_trans_res_zi)

cat("\n--- Sensitivity (b): transect family selection ---\n")
print(make_aicc_df(list(
  "Tweedie"    = p_trans_tw_base,
  "ZI Tweedie" = p_trans_tw_zi_base
)))

# Results:
#   Standard Tweedie: AICc = 2557.28, weight = 0.743
#   ZI Tweedie:       DAICc = 2.12,   weight = 0.257
#   Standard Tweedie selected — consistent with site-level.

# ── Transect Q1-Q3 sequence ───────────────────────────────────
p_trans_null <- glmmTMB(
  transect_pisc_biomass ~ 1 + (1 | site),
  family = tweedie(link = "log"),
  data   = pisc_transect_data)

p_trans_baseline <- glmmTMB(
  transect_pisc_biomass ~ rugosity_sc + log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = pisc_transect_data)

p_trans_pressure <- glmmTMB(
  transect_pisc_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = pisc_transect_data)

p_trans_conn <- glmmTMB(
  transect_pisc_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc +
    connectivity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = pisc_transect_data)

p_trans_conn_int <- glmmTMB(
  transect_pisc_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc +
    connectivity_sc +
    log_market_gravity_sc:connectivity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = pisc_transect_data)

p_trans_mpa_int <- glmmTMB(
  transect_pisc_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc +
    connectivity_sc +
    log_market_gravity_sc:connectivity_sc +
    mpa_status + (1 | site),
  family = tweedie(link = "log"),
  data   = pisc_transect_data)

cat("\n--- Sensitivity (b): transect model comparison ---\n")
print(make_aicc_df(list(
  "Null"               = p_trans_null,
  "Baseline"           = p_trans_baseline,
  "Baseline + press"   = p_trans_pressure,
  "Baseline + conn"    = p_trans_conn,
  "Conn x press"       = p_trans_conn_int,
  "Conn x press + MPA" = p_trans_mpa_int
)))

cat("\n--- Sensitivity (b): interaction model coefficients ---\n")
summary(p_trans_conn_int)

# ── Site random intercept SD ──────────────────────────────────
vc_p      <- VarCorr(p_trans_baseline)
site_sd_p <- sqrt(as.numeric(vc_p$cond$site))
cat(sprintf("\nSite random intercept SD = %.3f\n", site_sd_p))

# Results:
#   Conn x press + MPA: AICc = 2546.70, weight = 0.893
#   Conn x press:       DAICc = 4.92,   weight = 0.076
#   All other models:   DAICc > 9 — not competitive
#
#   Primary Q2 and Q3 findings fully replicated at
#   transect level — stronger support than site level
#   (site weight = 0.821).
#
#   Interaction: b = -0.529, z = -2.68, p = 0.007 **
#     Direction and magnitude consistent with site-level
#     (b = -0.511, p = 0.004).
#   Market gravity: b = +0.429, z = 2.52, p = 0.012 * — stable
#   Connectivity:   b = +0.292, z = 1.99, p = 0.047 * — stable
#   Site random intercept SD = 0.893
#   
#   DHARMa diagnostics — transect baseline:
#   KS test:        p = 0.483 — no significant deviation
#   Dispersion:     p = 0.944 — no significant deviation
#   Outlier test:   p = 0.253 — no significant outliers
#   Residuals vs predicted: no significant problems detected
#   Gaussian lmer structure confirmed appropriate.

# DHARMa diagnostics — transect conn x press + MPA:
#   KS test:        p = 0.341 — no significant deviation
#   Dispersion:     p = 0.900 — no significant deviation
#   Outlier test:   p = 0.075 — no significant outliers
#   Residuals vs predicted: quantile deviations detected
#     (red curve, lower quantile below expected). Combined
#     adjusted quantile test significant — mild
#     heteroscedasticity at low predicted values consistent
#     with high within-site variance (site SD = 0.893) and
#     patchy predator distribution. KS and dispersion pass
#     — acceptable overall.


# ============================================================
#  RESULTS SUMMARY
# ============================================================

cat("\n--- Piscivore results summary ---\n")
tribble(
  ~Question,   ~Result,           ~Key_finding,
  "Q1",        "Market gravity",  "weight = 0.488, DAICc = 0.96 vs baseline — marginal, ecologically motivated",
  "Q2 int",    "Supported",       "weight = 0.836, DAICc = 4.45; conn x mkt gravity b = -0.511, p = 0.004 **",
  "Q3 MPA",    "Supported",       "weight = 0.821, DAICc = 3.05; medium b = +0.868, p = 0.003; 2.38x (CI: 1.34x–4.23x)",
  "Moran",     "No SAC",          "I = +0.052, p = 0.162",
  "Pearson r", "0.490",           "moderate fit",
  "Sens (a)",  "Partly dependent","conn marginal without MG (p = 0.077); MPA uncertain (DAICc = 0.07); interaction key result",
  "Sens (b)",  "Consistent",      "conn x press + MPA weight = 0.893; interaction b = -0.529 p = 0.007; site SD = 0.893"
) %>% print()


# ── End of script ─────────────────────────────────────────────

