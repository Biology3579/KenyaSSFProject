# ============================================================
#  DRIVERS OF CORALLIVORE BIOMASS
#  Chapter 1 — Functional Group Analysis: Corallivores
#
#  Analytical framework mirrors total biomass and browsers
#  (Q1–Q3 + sensitivity).
#
#  Key differences from browsers:
#    Zero proportion: 0% at site level — Gaussian log
#    viable. Family selection determines approach.
#    Pressure metric: complete uncertainty in Q1
#    (all metrics within ΔAICc < 1) — settlement gravity
#    selected for consistency with other groups.
#    Ecological context: corallivores are coral-dependent
#    specialists not heavily targeted by SSF — expected
#    to show weaker pressure signal than browsers or
#    total biomass.
#
#  Sensitivity analysis:
#    (a) Alternative pressure metrics
#    (b) Transect-level replication
# ============================================================

# ── SOURCE SHARED DATA PREPARATION ───────────────────────────
# Loads: fish_data, scaled_predictors, final_predictors,
# total_transects, transect_model_data, total_model_data,
# make_aicc_df(), plot_effect(), and all packages.
source(here::here("data_preparation.R"))

# ============================================================
#  DATA AGGREGATION
# ============================================================

coralliv_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_coralliv_biomass = sum(
      ifelse(trophic_group == "corallivores", tot_wt_g, 0),
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

cat("Corallivore transects:", nrow(coralliv_transects), "\n")
cat("Sites:",                 n_distinct(coralliv_transects$site), "\n")
cat("Countries:",             n_distinct(coralliv_transects$country), "\n")

# ── Site-level dataset ────────────────────────────────────────
coralliv_model_data <- coralliv_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass              = mean(transect_coralliv_biomass,
                                     na.rm = TRUE),
    log_mean_biomass          = log(mean(transect_coralliv_biomass,
                                         na.rm = TRUE)),
    n_transects               = n(),
    rugosity_sc               = first(rugosity_sc),
    log_settlement_grav_sc    = first(log_settlement_grav_sc),
    log_chla_sc               = first(log_chla_sc),
    connectivity_sc           = first(connectivity_sc),
    mpa_status                = first(mpa_status),
    log_settlement_pop_sc     = first(log_settlement_pop_sc),
    log_market_gravity_sc     = first(log_market_gravity_sc),
    .groups = "drop"
  ) %>%
  mutate(
    site       = as.factor(site),
    country    = as.factor(country),
    mpa_status = factor(mpa_status,
                        levels  = c("none", "low", "medium"),
                        ordered = FALSE)
  )

cat("\nCorallivore model data:", nrow(coralliv_model_data), "sites\n")

# ── Data checks ───────────────────────────────────────────────
coralliv_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                log_chla_sc, connectivity_sc, mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

cat("\nZeros in mean_biomass:",
    sum(coralliv_model_data$mean_biomass == 0), "\n")
cat("-Inf in log_mean_biomass:",
    sum(is.infinite(coralliv_model_data$log_mean_biomass)), "\n")
cat("Site-level zero proportion:",
    round(mean(coralliv_model_data$mean_biomass == 0), 3), "\n")
cat("\nResponse summary:\n")
print(summary(coralliv_model_data$log_mean_biomass))

# ── MPA classification check ──────────────────────────────────
cat("\nMPA status counts:\n")
print(table(coralliv_model_data$mpa_status))

# ── Transect-level dataset ────────────────────────────────────
coralliv_transect_data <- coralliv_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(coralliv_transect_data$transect_coralliv_biomass == 0),
    "/", nrow(coralliv_transect_data),
    "(", round(mean(coralliv_transect_data$transect_coralliv_biomass == 0),
               3), ")\n")


# ============================================================
#  MODEL FAMILY SELECTION
#  Zero proportion determines viable families.
#  If site-level zeros = 0 and log distribution continuous:
#    Gaussian log (lm) is appropriate and preferred —
#    consistent with total biomass, varpart available.
#  If zeros present: Tweedie as per browsers.
# ============================================================

# ── Distribution plot ─────────────────────────────────────────
par(mfrow = c(1, 2))
hist(coralliv_model_data$mean_biomass,
     breaks = 30, main = "Raw", xlab = "Mean biomass (g)")
hist(coralliv_model_data$log_mean_biomass,
     breaks = 25, main = "Log-transformed",
     xlab = "log(Mean biomass)")
par(mfrow = c(1, 1))

# ── Gaussian log — global model ───────────────────────────────
coralliv_lm_global <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  data = coralliv_model_data
)

par(mfrow = c(2, 2))
plot(coralliv_lm_global, main = "Gaussian log")
par(mfrow = c(1, 1))

# ── Tweedie — global model ────────────────────────────────────
coralliv_tw_global <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

coralliv_tw_res <- simulateResiduals(coralliv_tw_global, n = 1000)
plot(coralliv_tw_res)
testZeroInflation(coralliv_tw_res)
testDispersion(coralliv_tw_res)

# ── Family selection decision ─────────────────────────────────
#
# Tweedie (log link): SELECTED
#   Zero-free response confirmed (0 site-level zeros) —
#   Tweedie not strictly required but selected for
#   consistency across all functional group analyses.
#   Using a common family enables direct comparison of
#   model structures and coefficients across groups
#   without justifying family differences.
#
#   DHARMa diagnostics (n = 1000):
#     KS test:        p = 0.976 — excellent fit
#     Dispersion:     p = 0.762 — acceptable
#     Zero inflation: ratio = 0, p = 1.000 — no zeros
#     Outlier test:   p = 1.000 — no outliers
#   No significant problems detected.
#
#   Gaussian log diagnostics adequate but not selected —
#   analytical consistency across functional groups
#   takes precedence over parsimony for a zero-free
#   response.
#
#  Proceed: glmmTMB(family = tweedie(link = "log")) on
#  raw mean_biomass throughout all corallivore analyses.
#  Consistent with browser analyses.


# ============================================================
#  RANDOM EFFECT STRUCTURE
#  Tweedie family — glmmTMB with and without (1|country).
#  Global predictor set — most demanding test.
#  Both models comparable via AICc directly (same family).
# ============================================================

coralliv_re_null <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

coralliv_re_country <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | country),
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

cat("\n--- Corallivore RE structure: country-level ---\n")
print(make_aicc_df(list(
  "No RE"         = coralliv_re_null,
  "(1 | country)" = coralliv_re_country
)))

# RE result: No RE: AICc = 557.55, weight = 0.820
# (1|country): ΔAICc = 3.03, weight = 0.180
# Country clustering not supported — consistent with
# total biomass (ΔAICc = 2.86) and browsers (ΔAICc = 3.03).
# All corallivore models fitted without random effects.

# ============================================================
#  Q1 — PRESSURE METRIC SELECTION
# ============================================================

# ── Without connectivity control ─────────────────────────────
c_q1_settgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

c_q1_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

c_q1_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_pop_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

cat("\n--- Q1: Corallivore pressure metric (without connectivity) ---\n")
print(make_aicc_df(list(
  "Settlement gravity"   = c_q1_settgrav,
  "Market gravity"       = c_q1_mktgrav,
  "Settlement pop. 25km" = c_q1_settpop
)))

# ── With connectivity control ─────────────────────────────────
c_q1_settgrav_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

c_q1_mktgrav_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

c_q1_settpop_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_pop_sc +
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

cat("\n--- Q1: Corallivore pressure metric (with connectivity) ---\n")
print(make_aicc_df(list(
  "Settlement gravity"   = c_q1_settgrav_conn,
  "Market gravity"       = c_q1_mktgrav_conn,
  "Settlement pop. 25km" = c_q1_settpop_conn
)))

cat("\n--- Q1: Corallivore coefficient summary ---\n")
summary(c_q1_settgrav)

# ── Q1 results ────────────────────────────────────────────────
# Without connectivity control:
#   Settlement gravity:    AICc = 562.09, weight = 0.394
#   Settlement pop. 25km:  ΔAICc = 0.41,  weight = 0.321
#   Market gravity:        ΔAICc = 0.65,  weight = 0.285
#   Complete metric uncertainty — all three within
#   ΔAICc < 1. No metric meaningfully preferred.
#
# With connectivity control:
#   Settlement gravity:    AICc = 558.50, weight = 0.467
#   Market gravity:        ΔAICc = 0.98,  weight = 0.286
#   Settlement pop. 25km:  ΔAICc = 1.28,  weight = 0.247
#   Settlement gravity marginally preferred but all
#   three still within ΔAICc < 2 — genuine uncertainty
#   persists even with spatial control.
#
# Complete metric uncertainty confirmed across both
# comparisons. Contrasts with total biomass (settlement
# gravity weight = 0.878) and browsers (weight = 0.581
# with connectivity). The absence of metric discrimination
# is ecologically informative — corallivores do not
# respond detectably to SSF pressure regardless of how
# it is measured, consistent with being a non-target
# functional group.
#
# Baseline coefficients (settlement gravity):
#   Rugosity:           β = +0.116, p = 0.197 ns
#   Settlement gravity: β = +0.088, p = 0.409 ns
#   Chla:               β = +0.239, p = 0.036 *
#
# Settlement gravity positive and non-significant —
# wrong direction for a fishing pressure effect.
# Confirms corallivore biomass is decoupled from
# exploitation intensity. Chla the only significant
# predictor — positive direction consistent with
# corallivores depending on coral prey which is
# supported by productive environments.
# Rugosity non-significant — corallivores less
# habitat-dependent than total biomass or browsers.
# Dispersion = 1.18 — much lower than browsers (13.3),
# reflecting less variance across sites.
#
# Primary metric: settlement gravity (log_settlement_grav_sc)
#   Selected for consistency with other functional groups.
#   No ecological justification for preferring any metric.
#   All three retained for sensitivity analysis.

# ============================================================
#  Q2 — ARE CONNECTIVITY AND MPA BASELINE DRIVERS?
#
#  Tweedie throughout — McFadden pseudo-R² used in place
#  of adj. R². Consistent with browser analyses.
# ============================================================

c_null <- glmmTMB(
  mean_biomass ~ 1,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

c_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

c_baseline_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

c_baseline_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

c_global_additive <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

c_q2_models <- list(
  "Null"                  = c_null,
  "Baseline"              = c_baseline,
  "Baseline + conn"       = c_baseline_conn,
  "Baseline + MPA"        = c_baseline_mpa,
  "Baseline + conn + MPA" = c_global_additive
)

cat("\n--- Q2: Corallivore model comparison (AICc ranked) ---\n")
print(make_aicc_df(c_q2_models))

# ── McFadden pseudo-R² relative to baseline ───────────────────
cat("\n--- Q2: Corallivore pseudo-R² relative to baseline ---\n")

null_ll_c     <- as.numeric(logLik(c_null))
baseline_r2_c <- 1 - (as.numeric(logLik(c_baseline)) / null_ll_c)

c_q2_models %>%
  imap_dfr(~ tibble(
    Model  = .y,
    McF_R2 = round(1 - (as.numeric(logLik(.x)) / null_ll_c), 3)
  )) %>%
  mutate(
    Delta_R2 = round(McF_R2 - baseline_r2_c, 3),
    Delta_R2 = ifelse(Model %in% c("Null", "Baseline"),
                      NA, Delta_R2)
  ) %>%
  print()

# ── Q2 results ────────────────────────────────────────────────
#
# AICc comparison:
#   Baseline + conn + MPA: AICc = 557.55, weight = 0.500 (BEST)
#   Baseline + conn:       ΔAICc = 0.95,  weight = 0.312
#   Null:                  ΔAICc = 2.69,  weight = 0.130
#   Baseline:              ΔAICc = 4.54,  weight = 0.052
#   Baseline + MPA:        ΔAICc = 8.85,  weight = 0.006
#
# McFadden pseudo-R² relative to baseline (R² = 0.010):
#   Baseline + conn:       ΔR² = +0.011
#   Baseline + MPA:        ΔR² = +0.002
#   Baseline + conn + MPA: ΔR² = +0.023
#
# KEY FINDINGS:
#
# Connectivity supported — Baseline + conn competitive
# (ΔAICc = 0.95, weight = 0.312) and combined weight
# of models including connectivity = 0.812. Connectivity
# adds more explanatory power than MPA (ΔR² = +0.011
# vs +0.002).
#
# MPA not independently supported — Baseline + MPA
# performs poorly (ΔAICc = 8.85, weight = 0.006),
# worse than the null model. MPA only contributes
# when connectivity is already in the model.
#
# Null model competitive (ΔAICc = 2.69, weight = 0.130)
# — there is genuine uncertainty about whether any
# predictor explains corallivore biomass reliably.
# Baseline performs worse than null (ΔAICc = 4.54) —
# rugosity, pressure and chla together explain less
# than nothing after penalty for parameters.
#
# This is the starkest contrast with total biomass
# and browsers — for corallivores the local baseline
# predictors have essentially no explanatory power.
# Corallivore biomass is weakly structured by
# connectivity and possibly MPA in combination, but
# the overall signal is weak (max McF R² = 0.033).
#
# Best-supported model: Baseline + conn + MPA
# (weight = 0.500). Used as reference for Q3.
# Note: Baseline + conn is nearly equally supported
# (ΔAICc = 0.95) — connectivity is the primary driver,
# MPA contribution is marginal.

cat("\n--- Q2: Baseline coefficients ---\n")
summary(c_baseline)

cat("\n--- Q2: Baseline + connectivity coefficients ---\n")
summary(c_baseline_conn)

cat("\n--- Q2: Baseline + conn + MPA coefficients ---\n")
summary(c_global_additive)

# ── Q2: Coefficient summaries ─────────────────────────────────
#
# Baseline:
#   Rugosity:           β = +0.116, p = 0.197 ns
#   Settlement gravity: β = +0.088, p = 0.409 ns
#   Chla:               β = +0.239, p = 0.036 *
#   Only chla significant — local baseline predictors
#   explain minimal variance for corallivores.
#
# Baseline + connectivity:
#   Rugosity:           β = +0.121, p = 0.159 ns
#   Settlement gravity: β = +0.115, p = 0.256 ns
#   Chla:               β = +0.182, p = 0.091 . (marginal)
#   Connectivity:       β = -0.231, p = 0.010 *
#   Connectivity significant and negative — well-connected
#   sites have lower corallivore biomass. Chla weakens
#   slightly when connectivity included. Settlement
#   gravity remains non-significant.
#
# Baseline + conn + MPA (best supported, weight = 0.500):
#   Rugosity:           β = +0.162, p = 0.059 . (marginal)
#   Settlement gravity: β = +0.201, p = 0.088 . (marginal)
#   Chla:               β = +0.193, p = 0.087 . (marginal)
#   Connectivity:       β = -0.359, p < 0.001 ***
#   MPA low:            β = +0.589, p = 0.032 *
#   MPA medium:         β = +0.453, p = 0.040 *
#
# Connectivity strengthens markedly when MPA included
# (β = -0.231 → -0.359) — MPA and connectivity share
# variance, and once MPA placement is accounted for
# the negative connectivity signal intensifies.
#
# Both MPA levels positive and significant in global
# model — contrasts with total biomass (MPA not
# supported) and browsers (only medium MPA supported).
# Low and medium protection both associated with higher
# corallivore biomass, suggesting MPAs benefit
# corallivores through mechanisms other than direct
# fishing pressure reduction — likely reduced destructive
# fishing practices (anchoring, blast fishing) that
# damage coral substrate and corallivore prey.
#
# Settlement gravity marginal (p = 0.088) in global
# model — positive direction throughout, confirming
# corallivores are not negatively affected by SSF
# pressure. Slight positive association may reflect
# productive coastal areas supporting both human
# settlement and coral prey.
#
# Coefficient stability across Q2 models:
#   connectivity:  — → -0.231* → -0.359***
#   mpa_low:       — →    —    → +0.589*
#   mpa_medium:    — →    —    → +0.453*
#   chla:       +0.239* → +0.182. → +0.193.
#
# Best-supported model: Baseline + conn + MPA
# (weight = 0.500). Used as reference for Q3.

# ============================================================
#  SPATIAL AUTOCORRELATION CHECK
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

cat("\n--- Spatial autocorrelation: corallivore best Q2 model ---\n")
print(moran.test(residuals(c_global_additive,
                           type = "pearson"), listw5_c))

# ── Spatial autocorrelation result ────────────────────────────
# Moran's I = -0.051, p = 0.669 — no significant spatial
# autocorrelation in corallivore best Q2 model residuals.
#
# Negative value indicates slight spatial dispersion —
# neighbouring sites tend to have more dissimilar
# residuals than expected by chance. Consistent with
# corallivore biomass varying discretely rather than
# showing continuous spatial gradients.
#
# Pattern consistent across all three analyses:
#   Total biomass: I = +0.140, p = 0.015 — weak positive
#   Browsers:      I = -0.074, p = 0.802 — no signal
#   Corallivores:  I = -0.051, p = 0.669 — no signal
#
# The positive spatial signal in total biomass is absent
# for both functional groups — suggesting it reflects
# broad-scale community-level geographic patterning
# rather than functional group-specific processes.
#
# Warnings: identical points and 4 sub-graphs — expected,
# as per total biomass and browsers.
#
# No spatial error modelling required.

## ── H1: MPA effectiveness depends on fishing intensity ────────
# connectivity_sc retained as main effect — present in reference
c_int_mpa_press <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status * log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

# ── H2: Connectivity buffers exploitation effects ─────────────
# mpa_status retained as main effect — present in reference
c_int_conn_press <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    mpa_status +
    connectivity_sc * log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

# ── H3: MPA effectiveness depends on larval supply ────────────
# mpa_status * connectivity_sc expands to include both
# main effects — log_settlement_grav_sc retained explicitly
c_int_mpa_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc +
    mpa_status * connectivity_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

c_q3_models <- list(
  "Reference (best Q2)"    = c_q3_reference,
  "H1: MPA × pressure"     = c_int_mpa_press,
  "H2: Conn × pressure"    = c_int_conn_press,
  "H3: MPA × connectivity" = c_int_mpa_conn
)

cat("\n--- Q3: Corallivore interaction comparison ---\n")
print(make_aicc_df(c_q3_models))

cat("\n--- Q3: H1 MPA × pressure coefficients ---\n")
summary(c_int_mpa_press)

cat("\n--- Q3: H2 Connectivity × pressure coefficients ---\n")
summary(c_int_conn_press)

cat("\n--- Q3: H3 MPA × connectivity coefficients ---\n")
summary(c_int_mpa_conn)

# ── Q3 results ────────────────────────────────────────────────
#
# AICc comparison (updated — all interactions include
# full reference main effects):
#   Reference (best Q2):     AICc = 557.55, weight = 0.668 (BEST)
#   H2: Conn × pressure:     ΔAICc = 2.85,  weight = 0.161
#   H1: MPA × pressure:      ΔAICc = 3.54,  weight = 0.114
#   H3: MPA × connectivity:  ΔAICc = 4.92,  weight = 0.057
#
# Reference model decisively preferred — no interaction
# improves on the additive Baseline + conn + MPA structure.
# Weight increases from 0.804 to 0.668 after correcting
# H1 and H2 to include full reference main effects —
# the previously reported result was slightly inflated
# by unfair model comparison. Conclusion unchanged.
#
# H1 (MPA × pressure): ΔAICc = 3.54 — not supported.
#   Interaction terms non-significant (low × pressure
#   β = +0.036, p = 0.913; medium × pressure β = -0.654,
#   p = 0.102). Connectivity strengthens when MPA ×
#   pressure included (β = -0.358***) — connectivity
#   signal is robust and independent of MPA structure.
#   Settlement gravity becomes significant as main effect
#   (β = +0.237, p = 0.048) once MPA interaction absorbs
#   variance — not a reliable signal given unsupported
#   interaction.
#
# H2 (Conn × pressure): ΔAICc = 2.85 — not supported.
#   Interaction term near-zero (β = -0.051, p = 0.673).
#   All main effects stable — connectivity β = -0.372***,
#   MPA low β = +0.612*, MPA medium β = +0.454*.
#   Additive structure confirmed — connectivity effect
#   does not depend on pressure level.
#
# H3 (MPA × connectivity): ΔAICc = 4.92 — not supported.
#   MPA low × connectivity non-significant (β = +0.414,
#   p = 0.807). MPA medium × connectivity non-significant
#   (β = -0.211, p = 0.294). Connectivity main effect
#   remains significant (β = -0.270, p = 0.044).
#
# Overall Q3 conclusion:
# No interaction supported for corallivores. The additive
# Baseline + conn + MPA structure from Q2 is the final
# best model — connectivity (negative) and MPA (positive
# for both levels) operate independently without
# modifying each other or the pressure-biomass
# relationship. Contrasts with browsers where MPA ×
# pressure was strongly supported (weight = 0.895).

# ============================================================
#  SENSITIVITY ANALYSIS
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
c_sens_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

c_sens_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_pop_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

cat("\n--- Sensitivity (a): corallivore alternative metrics ---\n")
cat("Market gravity:\n")
print(summary(c_sens_mktgrav)$coefficients$cond)
cat("\nSettlement population 25km:\n")
print(summary(c_sens_settpop)$coefficients$cond)

# ── Sensitivity (a) results ───────────────────────────────────
#
# Alternative pressure metrics — baseline structure retained,
# only pressure metric substituted.
#
# Market gravity (β = -0.019, p = 0.851):
#   Not significant, near-zero, wrong direction relative
#   to expectation. Chla marginal (β = +0.184, p = 0.071).
#   Rugosity non-significant. Pattern identical to
#   settlement gravity baseline.
#
# Settlement population 25km (β = -0.050, p = 0.601):
#   Not significant, near-zero. Chla marginal
#   (β = +0.182, p = 0.060). Rugosity non-significant.
#
# All three pressure metrics non-significant with
# coefficients near zero — complete metric-independence
# confirmed. Primary Q1 and Q2 conclusions robust:
# corallivore biomass is not structured by SSF pressure
# regardless of how it is measured.
#
# Chla consistently marginal across all three metrics
# (p = 0.036–0.091) — positive direction stable.
# This is the most robust predictor in the corallivore
# baseline, though significance varies with model
# structure.
#
# Rugosity non-significant across all metrics
# (β = 0.094–0.116, p > 0.16) — habitat complexity
# does not independently predict corallivore biomass,
# consistent with corallivores depending on coral
# prey availability rather than structural complexity.

# ── (b) Transect-level replication ───────────────────────────
# 17.7% zeros at transect level — standard Tweedie sufficient.
# ZI Tweedie not tested — zero proportion well below the
# threshold where zero inflation is a plausible concern
# (cf. browsers at 43.2% where ZI was also not needed).
# Standard Tweedie handles this proportion natively.

c_trans_tw <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data
)

c_trans_res <- simulateResiduals(c_trans_tw, n = 1000)
plot(c_trans_res)
testZeroInflation(c_trans_res)
testDispersion(c_trans_res)

# ── Transect Q2 sequence ──────────────────────────────────────
# Mirrors site-level Q2 — tests whether model ordering
# replicates at transect level.
# Key check: does negative connectivity signal persist
# when within-site variation is retained?

c_trans_null <- glmmTMB(
  transect_coralliv_biomass ~ 1 + (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data
)

c_trans_baseline <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data
)

c_trans_conn <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data
)

c_trans_mpa <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data
)

c_trans_global <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data
)

cat("\n--- Sensitivity (b): corallivore transect comparison ---\n")
print(make_aicc_df(list(
  "Null"                  = c_trans_null,
  "Baseline"              = c_trans_baseline,
  "Baseline + conn"       = c_trans_conn,
  "Baseline + MPA"        = c_trans_mpa,
  "Baseline + conn + MPA" = c_trans_global
)))

cat("\n--- Sensitivity (b) Coefficients  ---\n")
summary(c_trans_baseline)
summary(c_trans_global)

# ── Sensitivity (b) results ───────────────────────────────────
#
# AICc comparison (transect level):
#   Baseline + conn + MPA: AICc = 2319.21, weight = 0.779 (BEST)
#   Baseline + conn:       ΔAICc = 2.94,   weight = 0.179
#   Null:                  ΔAICc = 7.19,   weight = 0.021
#   Baseline:              ΔAICc = 7.57,   weight = 0.018
#   Baseline + MPA:        ΔAICc = 11.27,  weight = 0.003
#
# Model ordering consistent with site-level Q2 —
# Baseline + conn + MPA best-supported at both scales.
# Result clearer at transect level (weight = 0.779
# vs 0.500 at site level) — connectivity and MPA
# signals strengthen when within-site variation retained.
#
# Key findings replicated:
# (1) Connectivity supported — combined weight of models
#     including connectivity = 0.958. Signal robust to
#     aggregation level.
# (2) MPA contributes alongside connectivity — Baseline
#     + conn + MPA clearly preferred over Baseline + conn
#     (ΔAICc = 2.94).
# (3) Baseline no better than null (ΔAICc = 0.38) —
#     local predictors without connectivity explain
#     nothing for corallivores at either scale.
# (4) Baseline + MPA worst additive model (ΔAICc = 11.27)
#     — MPA without connectivity not supported, confirming
#     MPA contribution is conditional on connectivity.
#
# Baseline coefficients (n = 243 transects, 54 sites):
#   Rugosity:           β = +0.130, p = 0.161 ns
#   Settlement gravity: β = +0.081, p = 0.475 ns
#   Chla:               β = +0.234, p = 0.037 *
#   Site variance: 0.296 (SD = 0.544)
#   Identical pattern to site-level — chla only
#   significant predictor, local predictors uninformative.
#
# Best model coefficients (Baseline + conn + MPA):
#   Rugosity:           β = +0.153, p = 0.061 .
#   Settlement gravity: β = +0.170, p = 0.148 ns
#   Chla:               β = +0.164, p = 0.121 ns
#   Connectivity:       β = -0.383, p < 0.001 ***
#   MPA low:            β = +0.682, p = 0.013 *
#   MPA medium:         β = +0.427, p = 0.050 .
#   Site variance: 0.184 (SD = 0.429)
#
# Connectivity strengthens at transect level
# (β = -0.359 → -0.383) — genuine between-site
# effect, not aggregation artefact. Site variance
# reduces from 0.296 to 0.184 when connectivity
# and MPA included — both absorb meaningful
# between-site variance.
#
# Coefficient stability across scales:
#   connectivity: -0.359*** (site) → -0.383*** (transect)
#   mpa_low:      +0.589*   (site) → +0.682*   (transect)
#   mpa_medium:   +0.453*   (site) → +0.427.   (transect)
#
# Primary Q2 conclusions fully replicated —
# corallivore biomass structured by negative
# connectivity and positive MPA effects, independent
# of local ecological conditions at both scales.

# ============================================================
#  MARGINAL EFFECT PLOTS
#  Significant predictors in best Q2 model:
#  connectivity (negative) and MPA (positive).
#  Tweedie — predictions on response scale (raw biomass).
# ============================================================

# ── Connectivity effect ───────────────────────────────────────
conn_grid_c <- data.frame(
  connectivity_sc        = seq(
    min(coralliv_model_data$connectivity_sc),
    max(coralliv_model_data$connectivity_sc),
    length.out = 200),
  rugosity_sc            = 0,
  log_settlement_grav_sc = 0,
  log_chla_sc            = 0,
  mpa_status             = factor("none",
                                  levels = c("none", "low",
                                             "medium"))
)

conn_pred_c <- predict(c_global_additive,
                       newdata = conn_grid_c,
                       se.fit  = TRUE,
                       type    = "response",
                       re.form = NA)

conn_grid_c$fit <- conn_pred_c$fit
conn_grid_c$lwr <- conn_pred_c$fit - 1.96 * conn_pred_c$se.fit
conn_grid_c$upr <- conn_pred_c$fit + 1.96 * conn_pred_c$se.fit

p_c_conn <- ggplot(conn_grid_c,
                   aes(x = connectivity_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "#2c7bb6", alpha = 0.15) +
  geom_line(colour = "#2c7bb6", linewidth = 1.1) +
  geom_point(data = coralliv_model_data,
             aes(x = connectivity_sc, y = mean_biomass),
             colour = "grey40", size = 1.5,
             alpha = 0.5, inherit.aes = FALSE) +
  labs(x = "Connectivity (standardised)",
       y = "Corallivore biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ── MPA marginal means ────────────────────────────────────────
mpa_grid_c <- data.frame(
  mpa_status             = factor(c("none", "low", "medium"),
                                  levels = c("none", "low",
                                             "medium")),
  rugosity_sc            = 0,
  log_settlement_grav_sc = 0,
  log_chla_sc            = 0,
  connectivity_sc        = 0
)

mpa_pred_c <- predict(c_global_additive,
                      newdata = mpa_grid_c,
                      se.fit  = TRUE,
                      type    = "response",
                      re.form = NA)

mpa_grid_c$fit <- mpa_pred_c$fit
mpa_grid_c$lwr <- mpa_pred_c$fit - 1.96 * mpa_pred_c$se.fit
mpa_grid_c$upr <- mpa_pred_c$fit + 1.96 * mpa_pred_c$se.fit

p_c_mpa <- ggplot(mpa_grid_c,
                  aes(x = mpa_status, y = fit)) +
  geom_hline(yintercept = mpa_grid_c$fit[1],
             linetype   = "dashed",
             colour     = "grey70",
             linewidth  = 0.4) +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  colour    = "#0570b0",
                  linewidth = 0.7,
                  size      = 0.6) +
  scale_x_discrete(labels = c("none"   = "No MPA",
                              "low"    = "Low",
                              "medium" = "Medium")) +
  labs(x = "MPA status",
       y = "Corallivore biomass (g)") +
  theme_bw(base_size = 12) +
  theme(
    axis.title         = element_text(face = "bold"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

gridExtra::grid.arrange(p_c_conn, p_c_mpa, ncol = 2)


# ============================================================
#  RESULTS SUMMARY
# ============================================================

coralliv_results <- tribble(
  ~Question,  ~Best_model,           ~Key_finding,
  "Q1",       "Sett. gravity",       "complete metric uncertainty — all ΔAICc < 1",
  "Q2 conn",  "Baseline + conn",     "β = -0.231, p = 0.010; negative — well-connected sites lower biomass",
  "Q2 MPA",   "Baseline+conn+MPA",   "weight = 0.500; low β = +0.589*, medium β = +0.453*",
  "Q3 H1",    "Reference",           "ΔAICc = 12.41 — not supported",
  "Q3 H2",    "Reference",           "ΔAICc = 3.71 — not supported",
  "Q3 H3",    "Reference",           "ΔAICc = 4.92 — not supported"
)

cat("\n--- Corallivore results summary ---\n")
print(coralliv_results)


# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()

