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
#       One influential site removed prior to analysis
#       (mnemb; Kenya-Tanzania north ecoregion) —
#       unusually high piscivore biomass (1584g) relative
#       to structural complexity (rugosity = -2.80 SD) with
#       no formal protection. chumb retained — also flat
#       reef but medium MPA and ecologically interpretable
#       within Q3 framework.
#       Market gravity best supported in Q1 (weight = 0.490)
#       — carried forward. Interaction (market gravity ×
#       connectivity) strongly supported in Q2
#       (weight = 0.838). MPA strongly supported at Q3
#       (weight = 0.905). Medium MPA: 2.60x (p = 0.001).
#
#  Sensitivity analyses:
#       (a) Baseline as Q1 reference — evaluates whether
#           Q2-Q3 conclusions hold without market gravity.
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

# Piscivore transects: 243
# Sites: 54

# ── Transect-level dataset ────────────────────────────────────
pisc_transect_data <- pisc_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(pisc_transect_data$transect_pisc_biomass == 0),
    "/", nrow(pisc_transect_data),
    sprintf("(%.3f)\n",
            mean(pisc_transect_data$transect_pisc_biomass == 0)))

# Transect zeros: 110 / 243 (0.453)

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

# ── Remove influential site ───────────────────────────────────
# mnemb (Kenya-Tanzania north ecoregion) identified as an
# influential observation prior to analysis. The site exhibited
# unusually high piscivore biomass (1584g) relative to its
# structural complexity (rugosity = -2.80 SD, the lowest value
# in the dataset) with no formal protection, suggesting this
# flat reef habitat represents a structurally distinct reef
# type where the predictors in this model operate differently.
# chumb (same ecoregion, rugosity = -2.80 SD, biomass = 3199g)
# was retained — as a medium MPA site it is ecologically
# interpretable within the Q3 framework and its influence
# on MPA estimates is of direct scientific interest.

pisc_model_data <- pisc_model_data %>%
  filter(site != "mnemb")

cat("\nPiscivore model data after removing influential site:",
    nrow(pisc_model_data), "sites\n")

# 53 sites retained.

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

# NA check passed.
# Site zeros: 6 / 53 (0.113)
# Response: min = 0, median = 465g, mean = 809g, max = 4121g
# MPA: none = 29, low = 7, medium = 17


# ============================================================
#  MODEL FAMILY SELECTION
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
#     to extreme lower tail values. Heteroscedasticity
#     driven by zero/near-zero sites.
#
#   Tweedie (log link): SELECTED
#   DHARMa diagnostics (n = 1000):
#       KS test:        p = 0.987 — excellent fit
#       Dispersion:     p = 0.606, ratio = 1.131 — acceptable
#       Zero inflation: p = 1.000, ratio = 0.957 — not needed
#       Outlier test:   p = 1.000 — no outliers
#
#   Proceed: glmmTMB(family = tweedie(link = "log")) on
#   raw mean_biomass throughout all piscivore analyses.


# ============================================================
#  RANDOM EFFECT STRUCTURE
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

# Random effect structure:
#   No RE:            AICc = 783.59, weight = 0.782
#   (1 | ecoregion):  DAICc = 2.55,  weight = 0.218
#   Ecoregion RE not supported — consistent with all
#   previous functional groups. All models without RE.


# ============================================================
#  Q1 — HUMAN PRESSURE
#
#  Approach: AICc comparison of baseline vs baseline + each
#  pressure metric.
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
#   Baseline + mkt: AICc = 782.65, weight = 0.490 (BEST)
#   Baseline:       DAICc = 0.94,  weight = 0.306
#   Settlement pop: DAICc = 2.99,  weight = 0.110
#   Settl. gravity: DAICc = 3.28,  weight = 0.095
#
#   Market gravity best supported (weight = 0.490).
#   Ecological baseline competitive (DAICc = 0.94) —
#   genuine model selection uncertainty.
#   Market gravity carried forward as ecologically
#   motivated — piscivores are commercially targeted
#   large-bodied species.
#   Settlement pop. (DAICc = 2.99) and settlement gravity
#   (DAICc = 3.28) outside threshold — not evaluated
#   in sensitivity (a).

# ── Best Q1 model ─────────────────────────────────────────────
p_best_q1 <- p_q1_mktgrav

summary(p_q1_mktgrav)
print(confint(p_q1_mktgrav))


# ============================================================
#  Q2 — LARVAL CONNECTIVITY
#
#  Approach: three-model AICc comparison.
#  (1) Best Q1 — market gravity only
#  (2) Best Q1 + connectivity main effect
#  (3) Best Q1 + connectivity + market gravity × connectivity
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

cat("\n--- Q2: Connectivity ---\n")
print(make_aicc_df(list(
  "Best Q1"                      = p_best_q1,
  "Best Q1 + conn"               = p_q2_conn,
  "Best Q1 + conn + interaction" = p_q2_conn_int
)))

# Results:
#   Best Q1 + conn + int: AICc = 778.02, weight = 0.838
#   Best Q1:              DAICc = 4.63,  weight = 0.083
#   Best Q1 + conn:       DAICc = 4.73,  weight = 0.079
#
#   Interaction strongly supported (weight = 0.838,
#   DAICc = 4.63 vs market gravity alone).
#
#   Interaction: b = -0.507, CI [-0.859, -0.156], p = 0.005 **
#   Market gravity: b = +0.402, CI [+0.100, +0.704], p = 0.009 **
#   Connectivity:   b = +0.283, CI [+0.017, +0.548], p = 0.037 *
#   Rugosity:       b = +0.290, CI [+0.014, +0.565], p = 0.039 *
#   Chl-a:          b = +0.038, CI [-0.260, +0.336], p = 0.802
#
#   Interpretation: at low connectivity, greater market
#   access associates positively with piscivore biomass.
#   At high connectivity, greater market access associates
#   negatively — consistent with connectivity facilitating
#   commercial fishing access to otherwise well-replenished
#   stocks. Main effects of market gravity and connectivity
#   not interpretable independently of the interaction.
#
#   DHARMa diagnostics — best Q2 model:
#   Outlier test: p = 1.000 — no outliers

summary(p_q2_conn)
print(confint(p_q2_conn))

# ── Best Q2 model ─────────────────────────────────────────────
p_best_q2 <- p_q2_conn_int



# ── Best model summary ────────────────────────────────────────
cat("\n--- Best model: full summary ---\n")
summary(p_best_q2)


cat("\n--- Best model: 95% CIs ---\n")
print(confint(p_best_q2))

# ── DHARMa diagnostics ───────────────────────────────────────
cat("\n--- DHARMa diagnostics: best Q2 model ---\n")
p_best_q2_sim <- simulateResiduals(p_best_q2, n = 1000)
plot(p_best_q2_sim)
testOutliers(p_best_q2_sim)


# ============================================================
#  Q3 — FORMAL PROTECTION
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

# ── Raw biomass by MPA ────────────────────────────────────────
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

# ── DHARMa diagnostics ───────────────────────────────────────
cat("\n--- DHARMa diagnostics: Q3 model ---\n")
p_q3_sim <- simulateResiduals(p_q3_mpa, n = 1000)
plot(p_q3_sim)
testOutliers(p_q3_sim)

# Results:
#   Best Q2 + MPA: AICc = 773.51, weight = 0.905
#   Best Q2:       DAICc = 4.51,  weight = 0.095
#   MPA strongly supported (weight = 0.905).
#
#   Medium MPA: b = +0.957, CI [+0.381, +1.533], p = 0.001 **
#     fold difference = 2.60x (95% CI: 1.46x to 4.63x)
#   Low MPA:    b = +0.262, p = 0.543 ns
#   Interaction: b = -0.581, CI [-0.909, -0.253], p = 0.001 ***
#     Strengthened and stable with MPA included
#     (was -0.507, p = 0.005 in best Q2).
#   Market gravity: b = +0.448, p = 0.002 **
#   Rugosity:       b = +0.383, p = 0.004 **
#   Connectivity:   b = +0.139, p = 0.323 — attenuated with MPA
#
#   Raw means: none = 545g, low = 754g, medium = 1284g
#   DHARMa: outlier test p = 1.000 — no outliers

# ── Best Q3 model ─────────────────────────────────────────────
# Q3 is a separate question — best Q2 model unchanged.
p_best_q3 <- p_best_q2

# ── Predicted vs observed ─────────────────────────────────────
pred_p <- predict(p_best_q3, type = "response")
cat(sprintf("\nPearson r (predicted vs observed): %.3f\n",
            cor(pred_p, pisc_model_data$mean_biomass)))

# Pearson r = 0.492 — moderate fit.


# ============================================================
#  SPATIAL AUTOCORRELATION CHECK
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

# Spatial autocorrelation:
#   Moran's I = +0.022, p = 0.283 — no significant spatial
#   autocorrelation in residuals. No spatial error
#   modelling required.


# ============================================================
#  SENSITIVITY ANALYSES
#
#  (a) Baseline as Q1 reference
#      Market gravity wins Q1 but baseline competitive
#      (DAICc = 0.94). Sensitivity (a) evaluates whether
#      Q2-Q3 conclusions hold when baseline is carried
#      forward instead of market gravity.
#
#  (b) Transect-level replication (Tweedie GLMM)
# ============================================================

# ── (a) Baseline as Q1 reference ─────────────────────────────
cat("\n--- Sensitivity (a): baseline as Q1 reference ---\n")

p_sens_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data)

cat("\nQ2 — model comparison (from baseline):\n")
print(make_aicc_df(list(
  "Baseline"        = p_baseline,
  "Baseline + conn" = p_sens_conn
)))

cat("\nConnectivity coefficients:\n")
print(summary(p_sens_conn)$coefficients$cond)
print(confint(p_sens_conn))

p_sens_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = pisc_model_data)

cat("\nQ3 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + conn"       = p_sens_conn,
  "Baseline + conn + MPA" = p_sens_mpa
)))

cat("\nMPA coefficients:\n")
print(summary(p_sens_mpa)$coefficients$cond)
print(confint(p_sens_mpa))

# Sensitivity (a) results (baseline as Q1 reference):
#
#   Q2:
#   Baseline + conn: AICc = 782.76, weight = 0.602
#   Baseline:        DAICc = 0.82,  weight = 0.398
#   Connectivity only marginally supported without
#   market gravity (b = +0.265, CI [-0.014, +0.543],
#   p = 0.062 .) — weaker than primary Q2.
#   Direction consistent with primary analysis.
#
#   Q3:
#   Baseline + conn + MPA: AICc = 781.69, weight = 0.631
#   Baseline + conn:       DAICc = 1.07,  weight = 0.370
#   MPA competitive (DAICc = 1.07) without market gravity.
#   Medium MPA: b = +0.825, CI [+0.179, +1.471], p = 0.012 *
#     Direction consistent with primary Q3 (b = +0.957)
#     but weaker model support overall.
#   Low MPA: b = +0.214, p = 0.657 ns
#   Connectivity: b = +0.141, p = 0.367 — not significant
#
# Conclusion: primary Q2 and Q3 findings are partly
#   dependent on market gravity being in the model.
#   The interaction and strong MPA support are contingent
#   on market gravity selection in Q1. Direction of all
#   effects consistent across primary and sensitivity
#   analyses. Market gravity selection is consequential
#   for downstream model structure.


# ── (b) Transect-level replication ───────────────────────────
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

# Family selection (transect level):
#   Standard Tweedie: AICc = 2557.28, weight = 0.743
#   ZI Tweedie:       DAICc = 2.12,   weight = 0.257
#   Standard Tweedie selected — consistent with site-level.

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

# Sensitivity (b) results (transect level, n = 243, 54 sites):
#   Note: transect data includes all 54 sites — mnemb not
#   excluded at transect level since site-level exclusion
#   does not propagate to transect data.
#
#   Conn x press + MPA: AICc = 2546.70, weight = 0.893
#   Conn x press:       DAICc = 4.92,   weight = 0.076
#   Baseline + conn:    DAICc = 9.07,   weight = 0.010
#   All other models:   DAICc > 9 — not competitive
#
#   Primary Q2 and Q3 findings fully replicated at
#   transect level — stronger support than site level.
#   Interaction: b = -0.529, z = -2.68, p = 0.007 **
#     Direction and magnitude consistent with site-level
#     (b = -0.507, p = 0.005).
#   Market gravity: b = +0.429, z = 2.52, p = 0.012 *
#   Connectivity:   b = +0.292, z = 1.99, p = 0.047 *
#   Site random intercept SD = 0.893


# ============================================================
#  RESULTS SUMMARY
# ============================================================

cat("\n--- Piscivore results summary ---\n")
tribble(
  ~Question,   ~Result,            ~Key_finding,
  "Excl.",     "mnemb removed",    "rugosity = -2.80 SD; biomass = 1584g; unprotected — influential, distinct reef type. chumb retained (medium MPA, ecologically interpretable)",
  "Q1",        "Market gravity",   "w = 0.490, DAICc = 0.94 vs baseline; carried forward — ecologically motivated",
  "Q2 int",    "Supported",        "w = 0.838, DAICc = 4.63; interaction b = -0.507, CI [-0.859,-0.156], p = 0.005 **",
  "Q3 MPA",    "Strongly supp.",   "w = 0.905, DAICc = 4.51; medium b = +0.957, p = 0.001; 2.60x (CI: 1.46x-4.63x)",
  "Moran",     "No SAC",           "I = +0.022, p = 0.283",
  "Pearson r", "0.492",            "moderate fit",
  "Sens (a)",  "Partly dependent", "conn p = 0.062 without MG; MPA competitive (DAICc = 1.07, p = 0.012) — direction consistent",
  "Sens (b)",  "Consistent",       "conn x press + MPA w = 0.893; interaction b = -0.529, p = 0.007; site SD = 0.893"
) %>% print()


# ── End of script ─────────────────────────────────────────────