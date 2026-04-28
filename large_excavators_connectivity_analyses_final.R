# ============================================================
#  DRIVERS OF LARGE EXCAVATOR BIOMASS
#  Chapter 1 — Functional Group Analysis: Large Excavators
#
#  Large excavators are large-bodied parrotfish (Scaridae)
#  that excavate coral substrate. Heavily targeted by
#  spearfishing.
#
#  Analytical framework mirrors browser, corallivore, and
#  grazer-detritivore analyses.
#  Structure: Q1 (pressure) → Q2 (connectivity) → Q3 (MPA)
#
#  Key differences from other groups:
#    ~13% zeros at site level — Tweedie required.
#    ~58% zeros at transect level — ZI Tweedie tested.
#
#  Sensitivity analyses:
#    (a) Alternative pressure metrics
#    (b) Transect-level replication (Tweedie GLMM)
# ============================================================

source(here::here("data_preparation.R"))

# ============================================================
#  DATA AGGREGATION
# ============================================================

excavator_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_excavator_biomass = sum(
      ifelse(trophic_group == "large_excavators", tot_wt_g, 0),
      na.rm = TRUE
    ),
    country = first(country),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(
    site    = as.factor(site),
    country = as.factor(country)
  )

cat("Excavator transects:", nrow(excavator_transects), "\n")
cat("Sites:",               n_distinct(excavator_transects$site), "\n")

# ── Site-level dataset ────────────────────────────────────────
excavator_model_data <- excavator_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_excavator_biomass,
                                  na.rm = TRUE),
    n_transects            = n(),
    rugosity_sc            = first(rugosity_sc),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_chla_sc            = first(log_chla_sc),
    connectivity_sc        = first(connectivity_sc),
    mpa_status             = first(mpa_status),
    ecoregion              = first(ecoregion),
    log_settlement_pop_sc  = first(log_settlement_pop_sc),
    log_market_gravity_sc  = first(log_market_gravity_sc),
    .groups = "drop"
  ) %>%
  mutate(
    site       = as.factor(site),
    country    = as.factor(country),
    ecoregion  = as.factor(ecoregion),
    mpa_status = factor(mpa_status,
                        levels  = c("none", "low", "medium"),
                        ordered = FALSE)
  )

cat("\nExcavator model data:", nrow(excavator_model_data), "sites\n")

# ── Data checks ───────────────────────────────────────────────
excavator_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                log_chla_sc, connectivity_sc, mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

cat("\nZeros in mean_biomass:",
    sum(excavator_model_data$mean_biomass == 0), "\n")
cat("Site-level zero proportion:",
    round(mean(excavator_model_data$mean_biomass == 0), 3), "\n")
cat("\nResponse summary:\n")
print(summary(excavator_model_data$mean_biomass))

cat("\nMPA status counts:\n")
print(table(excavator_model_data$mpa_status))

# ── Transect-level dataset ────────────────────────────────────
excav_transect_data <- excavator_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(excav_transect_data$transect_excavator_biomass == 0),
    "/", nrow(excav_transect_data),
    "(", round(mean(excav_transect_data$transect_excavator_biomass == 0),
               3), ")\n")


# ============================================================
#  MODEL FAMILY SELECTION
#  Run on baseline model.
#  Zeros present — Tweedie required.
# ============================================================

excav_lm_baseline <- lm(
  log(mean_biomass + 0.01) ~ rugosity_sc +
    log_chla_sc,
  data = excavator_model_data
)

par(mfrow = c(2, 2))
plot(excav_lm_baseline, main = "Gaussian log — baseline")
par(mfrow = c(1, 1))

excav_tw_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data
)

excav_tw_res <- simulateResiduals(excav_tw_baseline, n = 1000)
plot(excav_tw_res)
testZeroInflation(excav_tw_res)
testDispersion(excav_tw_res)

# ── Family selection decision ─────────────────────────────────
#
# Gaussian log: REJECTED
#   Zero sites drive residuals to extreme values —
#   identical problem to browsers. No offset resolves
#   bimodal log distribution.
#
# Tweedie (log link): SELECTED
#   Handles zeros natively. DHARMa diagnostics: TBC.
#
#  Proceed: glmmTMB(family = tweedie(link = "log")) on
#  raw mean_biomass throughout.

# ============================================================
#  RANDOM EFFECT STRUCTURE
#  Tested on baseline model — same rationale as all
#  previous analyses. Tweedie family throughout.
# ============================================================

excav_re_null <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data
)

excav_re_ecoregion <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | ecoregion),
  family = tweedie(link = "log"),
  data   = excavator_model_data
)

cat("\n--- Excavator RE structure: ecoregion ---\n")
print(make_aicc_df(list(
  "No RE"           = excav_re_null,
  "(1 | ecoregion)" = excav_re_ecoregion
)))

# Random effect structure: ecoregion
# Tested on baseline model to avoid circularity.
# No RE: AICc = 686.05, weight = 0.781 (BEST)
# (1 | ecoregion): DAICc = 2.54,  weight = 0.220
# Ecoregion RE not supported — identical result to
# corallivores (DAICc = 2.54) and consistent with
# total biomass (DAICc = 2.25), browsers (TBC), and
# grazer-detritivores (DAICc = 2.93).

# ── Variance inflation factors ────────────────────────────────
# Note: no pressure stage — pressure not supported in Q1.
# Sequence: baseline → baseline + settgrav (Q1 check) →
#           connectivity → MPA (baseline + MPA only,
#           connectivity not in best Q3 model).

cat("\n--- VIF: baseline ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data))

cat("\n--- VIF: baseline + settlement gravity (Q1 check) ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data))

cat("\n--- VIF: baseline + connectivity ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + connectivity_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data))

cat("\n--- VIF: baseline + MPA ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + mpa_status,
  family = tweedie(link = "log"),
  data   = excavator_model_data))

# ============================================================
#  Q1 — HUMAN PRESSURE
#
#  Does human pressure explain variation in excavator
#  biomass beyond local ecological context?
#  Which metric best captures SSF exploitation intensity?
#
#  A priori prediction: settlement gravity outperforms
#  alternatives. If baseline best supported, no pressure
#  carried forward — consistent with corallivores and
#  grazer-detritivores.
# ============================================================

e_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data
)

e_q1_settgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data
)

e_q1_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data
)

e_q1_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_pop_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data
)

cat("\n--- Q1 Step 1: Excavator metric comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                      = e_baseline,
  "Baseline + settlement gravity" = e_q1_settgrav,
  "Baseline + market gravity"     = e_q1_mktgrav,
  "Baseline + settlement pop."    = e_q1_settpop
)))

cat("\n--- Q1 Step 2: Baseline coefficients ---\n")
summary(e_baseline)

# ── Q1: Rugosity range, fold difference and confint ──────────
cat("\n--- Excavator Q1: rugosity range and fold difference ---\n")
rug_range_e <- range(excavator_model_data$rugosity_sc, na.rm = TRUE)
cat("Range:", rug_range_e, "\n")
rug_span_e  <- diff(rug_range_e)
cat(sprintf("Span: %.3f SD units\n", rug_span_e))
b_rug_e     <- 0.506
fold_rug_e  <- exp(abs(b_rug_e * rug_span_e))
cat(sprintf("Fold difference (low vs high rugosity): %.2fx\n", fold_rug_e))

cat("\n--- Excavator Q1: baseline confint ---\n")
print(confint(e_baseline))

cat("\n--- Q1: Pressure metric direction checks ---\n")
cat("Settlement gravity:\n")
print(summary(e_q1_settgrav)$coefficients$cond[
  "log_settlement_grav_sc", ])
cat("\nMarket gravity:\n")
print(summary(e_q1_mktgrav)$coefficients$cond[
  "log_market_gravity_sc", ])
cat("\nSettlement population:\n")
print(summary(e_q1_settpop)$coefficients$cond[
  "log_settlement_pop_sc", ])

# Q1 results:
#   Baseline:           AICc = 686.05, weight = 0.494 (BEST)
#   Settlement pop.:    DAICc = 1.85,  weight = 0.196
#   Settlement gravity: DAICc = 2.13,  weight = 0.171
#   Market gravity:     DAICc = 2.53,  weight = 0.139
#
#   Baseline best supported — no pressure metric
#   outperforms ecological context for excavators.
#   Complete metric uncertainty — all within DAICc < 2.6.
#
# Baseline coefficients:
#   Rugosity: b = +0.506, p = 0.009 ** SIGNIFICANT
#     Strongest rugosity signal of all functional groups
#     — habitat structural complexity is the primary
#     driver of large excavator biomass. Excavators
#     depend on complex reef structure for shelter,
#     foraging substrate, and spawning aggregations.
#   Chla:     b = -0.300, p = 0.157 ns
#     Negative but not significant — retained as
#     baseline control.
#   Dispersion = 8.51 — high within-group variance
#     consistent with patchily distributed large-bodied
#     species.
#
# Pressure metric direction checks:
#   Settlement gravity: b = -0.223, p = 0.517 ns
#   Market gravity:     b = -0.021, p = 0.945 ns
#   Settlement pop.:    b = +0.269, p = 0.400 ns
#
#   Inconsistent directions — settlement gravity and
#   market gravity negative, settlement population
#   positive. All non-significant. No coherent pressure
#   signal for excavators. Fourth consecutive functional
#   group (corallivores, grazers, excavators) where
#   baseline outperforms all pressure metrics —
#   confirms pressure only detectable for total biomass
#   and browsers where SSF targeting is most direct.
#
# Best Q1 model: e_baseline (rugosity + chla)
# No pressure term carried forward.
e_best_q1 <- e_baseline

# ============================================================
#  Q2 — LARVAL CONNECTIVITY
#
#  Does connectivity explain additional variation beyond
#  the ecological baseline?
#
#  Main effect only — connectivity x pressure interaction
#  not testable because pressure is not a supported
#  predictor for excavators (Q1: baseline best supported,
#  weight = 0.494). Interaction requires pressure as a
#  main effect first.
#
#  Tweedie family throughout.
# ============================================================

e_q2_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data
)

cat("\n--- Q2: Excavator connectivity main effect ---\n")
print(make_aicc_df(list(
  "Best Q1"        = e_best_q1,
  "Best Q1 + conn" = e_q2_conn
)))

cat("\n--- Q2: Connectivity coefficients ---\n")
summary(e_q2_conn)

# ── Q2: Connectivity confint ─────────────────────────────────
cat("\n--- Excavator Q2: confint ---\n")
print(confint(e_q2_conn))

# Q2 results:
#   Best Q1:        AICc = 686.05, weight = 0.758 (BEST)
#   Best Q1 + conn: DAICc = 2.28,  weight = 0.242
#   Connectivity not supported (DAICc = 2.28).
#
#   Connectivity: b = -0.129, p = 0.614 ns
#     Not significant, near zero. No independent
#     effect on excavator biomass once habitat and
#     productivity controlled.
#   Rugosity:     b = +0.527, p = 0.008 ** — stable
#   Chla:         b = -0.320, p = 0.127 ns — stable
#
#   Contrasts with corallivores (b = -0.221, p = 0.015,
#   supported) and grazers (b = -0.196, p = 0.058,
#   marginal) where negative connectivity was detected.
#   Despite similar direction (negative), effect size
#   is small and clearly not supported by AICc.
#   High within-group variance (dispersion = 8.5)
#   likely reduces power to detect weak signals.

# ── Best Q2 model ─────────────────────────────────────────────
# Connectivity not supported — baseline retained.
e_best_q2 <- e_best_q1  # rugosity + chla only
# ============================================================
#  Q3 — FORMAL PROTECTION
#
#  Does MPA status explain additional variation beyond
#  the ecological baseline?
#
#  MPA x pressure not testable — pressure not supported
#  in Q1. MPA x connectivity not tested — MPA not
#  clearly supported as main effect (DAICc = 1.22).
#
#  Tweedie family throughout.
# ============================================================

e_q3_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = excavator_model_data
)

cat("\n--- Q3: Excavator MPA main effect ---\n")
print(make_aicc_df(list(
  "Best Q2"       = e_best_q2,
  "Best Q2 + MPA" = e_q3_mpa
)))

cat("\n--- Q3: MPA coefficients ---\n")
summary(e_q3_mpa)

# ── Q3: MPA placement check ──────────────────────────────────
cat("\n--- Excavator Q3: pressure and connectivity by MPA ---\n")
excavator_model_data %>%
  group_by(mpa_status) %>%
  summarise(
    n                    = n(),
    mean_biomass         = round(mean(mean_biomass), 1),
    mean_settlement_grav = round(mean(log_settlement_grav_sc), 3),
    mean_connectivity    = round(mean(connectivity_sc), 3),
    .groups = "drop"
  ) %>%
  print()

# Q3 results:
#   Best Q2:       AICc = 686.05, weight = 0.648 (BEST)
#   Best Q2 + MPA: DAICc = 1.22,  weight = 0.353
#   MPA not clearly supported — genuine model selection
#   uncertainty but best Q2 preferred.
#
#   Low MPA:    b = -1.760, p = 0.029 * — significant
#     but NEGATIVE and artefactual. Low MPA sites
#     confined to narrow high-connectivity range
#     (z = 0.66 to max, n = 7) with mean biomass
#     70.9g vs 621g for unprotected sites. Negative
#     coefficient reflects site characteristics, not
#     a genuine negative protection effect.
#   Medium MPA: b = -0.437, p = 0.427 ns
#     Not significant. Mean biomass (640g) essentially
#     identical to unprotected sites (621g).
#   Rugosity: b = +0.516, p = 0.014 * — stable.
#
#   MPA site distribution:
#     none:   n = 30, mean = 621g
#     low:    n = 7,  mean = 71g  — very low biomass,
#             high-connectivity sites (z = 0.66 to max)
#     medium: n = 17, mean = 640g — similar to none
#
#   MPA not interpreted as a driver for excavators.
#   Low MPA negative coefficient is same artefact
#   pattern as corallivores — low MPA sites have
#   atypical characteristics that depress biomass
#   independently of protection status.
#   Baseline model retained as best throughout Q3.
#
#   MPA x connectivity interaction not tested —
#   MPA not supported as main effect.

# ── Best Q3 model ─────────────────────────────────────────────
# Baseline retained — MPA not clearly supported.
e_best_q3 <- e_best_q2  # rugosity + chla only

cat("\n--- Excavator: predicted vs observed ---\n")
pred_e <- predict(e_best_q3, type = "response")
obs_e  <- excavator_model_data$mean_biomass
cat(sprintf("Pearson r: %.3f\n", cor(pred_e, obs_e)))

cat("\n--- Q3: DHARMa diagnostics ---\n")
excav_sim <- simulateResiduals(e_best_q3, n = 1000)
plot(excav_sim)
testOutliers(excav_sim)

# ============================================================
#  SPATIAL AUTOCORRELATION
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

cat("\n--- Spatial autocorrelation: excavator best model ---\n")
print(moran.test(residuals(e_best_q3, type = "pearson"),
                 listw5_e))

# Spatial autocorrelation: excavator best model
# (rugosity + chla)
# Moran's I = -0.015, p = 0.479 — no significant
# spatial autocorrelation in residuals.
#
# Consistent with browsers (I = -0.078, p = 0.811)
# and corallivores (I = -0.021, p = 0.509).
# Contrasts with total biomass (I = 0.140, p = 0.015)
# and grazer-detritivores (I = 0.210, p = 0.001).
# Rugosity adequately captures the spatial variation
# in excavator biomass without leaving a residual
# geographic signal. No spatial error modelling
# required.

# ============================================================
#  SENSITIVITY ANALYSIS
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
e_sens_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data
)

e_sens_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_pop_sc,
  family = tweedie(link = "log"),
  data   = excavator_model_data
)

cat("\n--- Sensitivity (a): excavator alternative metrics ---\n")
cat("Market gravity:\n")
print(summary(e_sens_mktgrav)$coefficients$cond)
cat("\nSettlement population:\n")
print(summary(e_sens_settpop)$coefficients$cond)

# Sensitivity (a) results:
#   Market gravity:      b = -0.021, p = 0.945 ns
#   Settlement pop.:     b = +0.269, p = 0.400 ns
#
#   Consistent with Q1 direction checks — inconsistent
#   directions (market gravity negative, settlement
#   population positive) and both non-significant.
#   Rugosity significant and stable across all metrics:
#     Settlement gravity: b = +0.506, p = 0.009 **
#     Market gravity:     b = +0.499, p = 0.022 *
#     Settlement pop.:    b = +0.585, p = 0.007 **
#   Rugosity is the sole robust predictor regardless
#   of pressure metric. Q1 null conclusion robust.

# ── (b) Transect-level replication ───────────────────────────
# ~58% zeros at transect level — ZI Tweedie tested.

e_trans_tw_base <- glmmTMB(
  transect_excavator_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = excav_transect_data
)

e_trans_tw_zi_base <- glmmTMB(
  transect_excavator_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = excav_transect_data
)

e_trans_res    <- simulateResiduals(e_trans_tw_base,    n = 500)
e_trans_res_zi <- simulateResiduals(e_trans_tw_zi_base, n = 500)

plot(e_trans_res);    testZeroInflation(e_trans_res)
plot(e_trans_res_zi); testZeroInflation(e_trans_res_zi)

cat("\n--- Sensitivity (b): transect family selection ---\n")
print(make_aicc_df(list(
  "Tweedie"    = e_trans_tw_base,
  "ZI Tweedie" = e_trans_tw_zi_base
)))

# Transect family selection:
#   Standard Tweedie: AICc = 1855.91, weight = 0.743
#   ZI Tweedie:       DAICc = 2.12,   weight = 0.257
#   Standard Tweedie selected — consistent with site-
#   level family selection. ZI not needed despite 58%
#   transect-level zeros.

e_trans_null <- glmmTMB(
  transect_excavator_biomass ~ 1 + (1 | site),
  family = tweedie(link = "log"),
  data   = excav_transect_data
)

e_trans_baseline <- glmmTMB(
  transect_excavator_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = excav_transect_data
)

# Pressure — for completeness, update if Q1 supports
e_trans_pressure <- glmmTMB(
  transect_excavator_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = excav_transect_data
)

# Connectivity — update predictors once Q1 known
e_trans_conn <- glmmTMB(
  transect_excavator_biomass ~ rugosity_sc +
    log_chla_sc +
    # log_settlement_grav_sc +  # uncomment if Q1 supports
    connectivity_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = excav_transect_data
)

# MPA — update predictors once Q2 known
e_trans_mpa <- glmmTMB(
  transect_excavator_biomass ~ rugosity_sc +
    log_chla_sc +
    # log_settlement_grav_sc +  # uncomment if Q1 supports
    # connectivity_sc +          # uncomment if Q2 supports
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = excav_transect_data
)

cat("\n--- Sensitivity (b): excavator transect comparison ---\n")
print(make_aicc_df(list(
  "Null"             = e_trans_null,
  "Baseline"         = e_trans_baseline,
  "Baseline + press" = e_trans_pressure,
  "Baseline + conn"  = e_trans_conn,
  "Best + MPA"       = e_trans_mpa
)))

cat("\n--- Sensitivity (b): baseline coefficients ---\n")
summary(e_trans_baseline)

vc_e <- VarCorr(e_trans_baseline)
site_sd_e <- sqrt(as.numeric(vc_e$cond$site))
cat(sprintf("\nSite random intercept SD = %.3f\n", site_sd_e))

# Sensitivity (b) results:
# AICc comparison (n = 243 transects, 54 sites):
#   Best + MPA:        AICc = 1855.63, weight = 0.308
#   Baseline:          DAICc = 0.28,   weight = 0.268
#   Null:              DAICc = 1.16,   weight = 0.173
#   Baseline + conn:   DAICc = 1.36,   weight = 0.156
#   Baseline + press:  DAICc = 2.36,   weight = 0.095
#
#   Complete model selection uncertainty at transect
#   level — all five models within DAICc < 2.5.
#   Null model third-ranked (weight = 0.173) confirms
#   excavator biomass is the most weakly structured
#   functional group at both analytical scales.
#   Pressure not supported (DAICc = 2.36) — consistent
#   with Q1. MPA marginally top-ranked but within
#   genuine uncertainty (DAICc = 0.28 vs baseline)
#   — consistent with Q3 ambiguity at site level.
#
# Baseline coefficients (REML):
#   Rugosity:   b = +0.605, t = 1.953 . (marginal)
#     Direction consistent with site-level (b = +0.506,
#     p = 0.009) — positive rugosity signal preserved
#     but significance lost at transect level due to
#     extremely high within-site variance.
#   Chla:       b = -0.405, t = -1.234 ns — stable
#   Site variance:    3.396 (SD = 1.843)
#   Dispersion:       26.1 — extremely high, highest
#     of all functional groups. Reflects patchy
#     distribution of large-bodied excavators within
#     sites — individual transects either catch large
#     parrotfish or not at all.
#   ICC = 3.396 / (3.396 + residual) — very high
#     between-site clustering despite weak fixed
#     effects, confirming (1|site) essential.
#
# Overall: rugosity direction consistent across scales
# but significance not replicated at transect level
# due to extremely high within-site variance (SD =
# 1.843). Site-level result (b = +0.506, p = 0.009)
# is the primary inference. Complete model uncertainty
# at transect level mirrors site-level pattern where
# null model was competitive. Excavators are the most
# weakly and patchily structured functional group.


# ============================================================
#  MARGINAL EFFECT PLOTS
#  Update once Q3 complete.
# ============================================================

best_model_e <- e_best_q3  # UPDATE after Q3 confirmed

# ── Rugosity effect ───────────────────────────────────────────
rug_grid_e <- data.frame(
  rugosity_sc     = seq(
    min(excavator_model_data$rugosity_sc),
    max(excavator_model_data$rugosity_sc),
    length.out = 200),
  log_chla_sc     = 0
  # add connectivity_sc = 0 if Q2 supported
  # add mpa_status if Q3 supported
)

rug_pred_e     <- predict(best_model_e,
                          newdata = rug_grid_e,
                          se.fit  = TRUE,
                          type    = "response",
                          re.form = NA)
rug_grid_e$fit <- rug_pred_e$fit
rug_grid_e$lwr <- rug_pred_e$fit - 1.96 * rug_pred_e$se.fit
rug_grid_e$upr <- rug_pred_e$fit + 1.96 * rug_pred_e$se.fit

p_e_rugosity <- ggplot(rug_grid_e,
                       aes(x = rugosity_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "#2c7bb6", alpha = 0.15) +
  geom_line(colour = "#2c7bb6", linewidth = 1.1) +
  geom_point(data = excavator_model_data,
             aes(x = rugosity_sc, y = mean_biomass),
             colour = "grey40", size = 1.5,
             alpha  = 0.5, inherit.aes = FALSE) +
  labs(x = "Rugosity (standardised)",
       y = "Excavator biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

print(p_e_rugosity)

# ============================================================
#  MARGINAL EFFECT PLOTS
#  Best model: rugosity + chla (e_best_q3 = e_baseline)
#  Rugosity is the sole significant predictor.
#  No connectivity or MPA plots — neither supported.
#  Predictions on response scale (raw biomass).
#  Tweedie — re.form = NA for population-level prediction.
# ============================================================

best_model_e <- e_best_q3  # rugosity + chla only

# ── Rugosity effect ───────────────────────────────────────────
rug_grid_e <- data.frame(
  rugosity_sc = seq(
    min(excavator_model_data$rugosity_sc),
    max(excavator_model_data$rugosity_sc),
    length.out = 200),
  log_chla_sc = 0
)

rug_pred_e     <- predict(best_model_e,
                          newdata = rug_grid_e,
                          se.fit  = TRUE,
                          type    = "response",
                          re.form = NA)
rug_grid_e$fit <- rug_pred_e$fit
rug_grid_e$lwr <- rug_pred_e$fit - 1.96 * rug_pred_e$se.fit
rug_grid_e$upr <- rug_pred_e$fit + 1.96 * rug_pred_e$se.fit

p_e_rugosity <- ggplot(rug_grid_e,
                       aes(x = rugosity_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "#2c7bb6", alpha = 0.15) +
  geom_line(colour = "#2c7bb6", linewidth = 1.1) +
  geom_point(data = excavator_model_data,
             aes(x = rugosity_sc, y = mean_biomass),
             colour = "grey40", size = 1.5,
             alpha  = 0.5, inherit.aes = FALSE) +
  labs(x = "Rugosity (standardised)",
       y = "Excavator biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

print(p_e_rugosity)

# jpeg("excavator_marginal_effects.jpg",
#      width = 11, height = 11, units = "cm", res = 300)
# print(p_e_rugosity)
# dev.off()


# ============================================================
#  RESULTS SUMMARY
# ============================================================

cat("\n--- Excavator results summary ---\n")
tribble(
  ~Question,   ~Result,           ~Key_finding,
  "Q1",        "Baseline best",   "weight = 0.494, all metrics ns, inconsistent directions",
  "Q2 conn",   "Not supported",   "DAICc = 2.28, weight = 0.242, b = -0.129, p = 0.614",
  "Q3 MPA",    "Not supported",   "DAICc = 1.22, weight = 0.353 — low MPA artefactual",
  "Spatial",   "Clean",           "Moran's I = -0.015, p = 0.479 — no autocorrelation",
  "Sens (a)",  "Consistent",      "market b = -0.021 ns, pop b = +0.269 ns — rugosity sole driver",
  "Sens (b)",  "Consistent",      "complete model uncertainty at transect level, rugosity direction stable"
) %>% print()


# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()