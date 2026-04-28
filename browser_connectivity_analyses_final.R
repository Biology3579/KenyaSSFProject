# ============================================================
#  DRIVERS OF BROWSER BIOMASS
#  Chapter 1 — Functional Group Analysis: Browsers
#
#  Analytical framework mirrors total biomass (Q1–Q3 + 
#  sensitivity), with group-specific pressure metric
#  selection and family selection.
#
#  Scientific questions:
#
#  Q1 — Pressure metric validity
#       Does settlement gravity outperform market gravity
#       and settlement population for browser biomass
#       specifically? Functional groups may respond
#       differently to pressure metrics depending on
#       their vulnerability and targeting by SSFs.
#
#  Q2 — Are connectivity and MPA baseline drivers?
#       Does connectivity explain browser biomass
#       independently of pressure and ecological context?
#       Does MPA status explain additional variance
#       beyond the baseline?
#
#  Q3 — Do management and connectivity modify pressure
#       effects on browser biomass?
#       Three a priori interaction hypotheses tested
#       against best Q2 model.
#
#  Baseline model (fixed a priori):
#       biomass ~ rugosity + pressure_metric + chla
#       Identical structure to total biomass baseline.
#       Pressure metric selected in Q1.
#
#  Key difference from total biomass:
#       Browser biomass has ~11% zeros at site level.
#       Tweedie distribution selected over Gaussian log
#       (see family selection). All models use
#       glmmTMB(family = tweedie(link = "log")) on
#       raw mean_biomass throughout.
#
#  Sensitivity analysis:
#       (a) Alternative pressure metrics
#       (b) Transect-level replication (Tweedie GLMM)
# ============================================================

# ── SOURCE SHARED DATA PREPARATION ───────────────────────────
# Loads: fish_data, scaled_predictors, final_predictors,
# total_transects, transect_model_data, total_model_data,
# make_aicc_df(), plot_effect(), and all packages.
source(here::here("data_preparation.R"))


# ============================================================
#  DATA AGGREGATION
# ============================================================

browser_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_browser_biomass = sum(
      ifelse(trophic_group == "browsers", tot_wt_g, 0),
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

cat("Browser transects:", nrow(browser_transects), "\n")
cat("Sites:",             n_distinct(browser_transects$site), "\n")
cat("Countries:",         n_distinct(browser_transects$country), "\n")

# ── Site-level dataset ────────────────────────────────────────
browser_model_data <- browser_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_browser_biomass,
                                  na.rm = TRUE),
    n_transects            = n(),
    rugosity_sc            = first(rugosity_sc),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_chla_sc            = first(log_chla_sc),
    connectivity_sc        = first(connectivity_sc),
    mpa_status             = first(mpa_status),
    log_settlement_pop_sc  = first(log_settlement_pop_sc),
    log_market_gravity_sc  = first(log_market_gravity_sc),
    .groups = "drop"
  ) %>%
  mutate(
    site       = as.factor(site),
    country    = as.factor(country),
    mpa_status = factor(mpa_status,
                        levels  = c("none", "low", "medium"),
                        ordered = FALSE)
  )

cat("\nBrowser model data:", nrow(browser_model_data), "sites\n")

# ── Data checks ───────────────────────────────────────────────
browser_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                log_chla_sc, connectivity_sc, mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

cat("\nZeros in mean_biomass:",
    sum(browser_model_data$mean_biomass == 0), "\n")
cat("Site-level zero proportion:",
    round(mean(browser_model_data$mean_biomass == 0), 3), "\n")
cat("\nResponse summary:\n")
print(summary(browser_model_data$mean_biomass))

# ── MPA classification check ──────────────────────────────────
cat("\nMPA status counts:\n")
print(table(browser_model_data$mpa_status))

# ── Transect-level dataset ────────────────────────────────────
browser_transect_data <- browser_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(browser_transect_data$transect_browser_biomass == 0),
    "/", nrow(browser_transect_data),
    "(", round(mean(browser_transect_data$transect_browser_biomass == 0),
               3), ")\n")


# ============================================================
#  MODEL FAMILY SELECTION
#  Run once on global model structure.
#  Browser biomass has ~11% zeros at site level —
#  Gaussian log requires an offset constant for zeros.
#  Log-transformed distribution is bimodal: zero sites
#  cluster at log(offset) separated from the main
#  distribution by a gap reflecting genuine biological
#  absence. No offset constant resolves this.
#  Tweedie handles zeros natively without offset.
#
#  Two families tested:
#    Gaussian log: lm() on log(mean_biomass + offset)
#    Tweedie:      glmmTMB() on raw mean_biomass, log link
#
#  DHARMa used for Tweedie diagnostics.
#  Standard plot() for Gaussian log diagnostics.
# ============================================================

# ── Gaussian log — global model ───────────────────────────────
browser_lm_global <- lm(
  log(mean_biomass + 0.01) ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  data = browser_model_data
)

par(mfrow = c(2, 2))
plot(browser_lm_global, main = "Gaussian log")
par(mfrow = c(1, 1))

# ── Tweedie — global model ────────────────────────────────────
browser_tw_global <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

browser_tw_res <- simulateResiduals(browser_tw_global, n = 1000)
plot(browser_tw_res)
testZeroInflation(browser_tw_res)
testDispersion(browser_tw_res)

# ── Family selection decision ─────────────────────────────────
#
# Gaussian log: REJECTED
#   Residuals vs Fitted: strong downward curve at low
#     fitted values — zero sites (sites 27, 39, 42)
#     pulling residuals to -10, creating severe non-
#     linearity that no offset constant can resolve.
#   Q-Q: severe lower tail deviation — sites 27, 39
#     fall far from theoretical line.
#   Scale-Location: strong downward trend —
#     heteroscedasticity throughout fitted range.
#   Residuals vs Leverage: sites 27, 42, 32 approach
#     or exceed Cook's distance threshold — zero sites
#     unduly influential on fitted coefficients.
#
# Tweedie (log link): SELECTED
#   Handles zeros natively without offset.
#   DHARMa diagnostics (n = 1000 simulations):
#     KS test:        p = 0.762 — good fit
#     Dispersion:     p = 0.222, ratio = 1.615 — acceptable
#     Zero inflation: p = 0.902, ratio = 0.872 — not needed
#     Outlier test:   p = 0.102 — no significant outliers
#   Residuals vs predicted: one flagged outlier (red *)
#     at high predicted value — minor, does not affect
#     overall fit assessment. No significant problems
#     detected by DHARMa.
#
#  Proceed: glmmTMB(family = tweedie(link = "log")) on
#  raw mean_biomass throughout all browser analyses.

# ============================================================
#  RANDOM EFFECT STRUCTURE
#  Country RE tested as per total biomass procedure.
#  Global predictor set — most demanding test.
#  Tweedie family used throughout for consistency.
# ============================================================

browser_re_null <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

browser_re_country <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | country),
  family = tweedie(link = "log"),
  data   = browser_model_data
)

cat("\n--- Browser RE structure: country-level ---\n")
print(make_aicc_df(list(
  "No RE"         = browser_re_null,
  "(1 | country)" = browser_re_country
)))

# RE result: No RE: AICc = 829.31, weight = 0.820
# (1|country): ΔAICc = 3.03, weight = 0.180
# Country clustering not supported — consistent with
# total biomass (ΔAICc = 2.86). All browser models
# fitted without random effects throughout.


# ============================================================
#  Q1 — PRESSURE METRIC SELECTION
#
#  Browser-specific Q1 — pressure metric may differ from
#  total biomass. Browsers are a heavily targeted functional
#  group; the metric that best captures exploitation
#  intensity for the total community may not be optimal
#  for a selectively harvested group.
#
#  Three metrics compared: settlement gravity, market
#  gravity, and settlement population (25km radius).
#  Settlement population included here despite lacking
#  an access cost component — its performance relative
#  to gravity metrics is itself informative about
#  whether direct exploitation or broader human
#  footprint drives browser biomass.
#
#  Same structure as total biomass Q1:
#  Three models identical except pressure metric.
#  Baseline controls (rugosity + chla) held constant.
#  Robustness confirmed with and without connectivity.
#  Tweedie family used throughout.
# ============================================================

# ── Without connectivity control ─────────────────────────────
b_q1_settgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

b_q1_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

b_q1_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_pop_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

cat("\n--- Q1: Browser pressure metric (without connectivity) ---\n")
print(make_aicc_df(list(
  "Settlement gravity"  = b_q1_settgrav,
  "Market gravity"      = b_q1_mktgrav,
  "Settlement pop. 25km" = b_q1_settpop
)))

# ── With connectivity control ─────────────────────────────────
# Confirms ranking robust to inclusion of spatial structure
b_q1_settgrav_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

b_q1_mktgrav_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

b_q1_settpop_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_pop_sc +
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

cat("\n--- Q1: Browser pressure metric (with connectivity) ---\n")
print(make_aicc_df(list(
  "Settlement gravity"   = b_q1_settgrav_conn,
  "Market gravity"       = b_q1_mktgrav_conn,
  "Settlement pop. 25km" = b_q1_settpop_conn
)))

# ── Q1 results ────────────────────────────────────────────────
#
# Without connectivity control:
#   Settlement gravity:    AICc = 837.50, weight = 0.379
#   Market gravity:        ΔAICc = 0.04,  weight = 0.371
#   Settlement pop. 25km:  ΔAICc = 0.83,  weight = 0.250
#   All three metrics within ΔAICc < 1 — complete model
#   selection uncertainty without spatial control.
#
# With connectivity control:
#   Settlement gravity:    AICc = 832.03, weight = 0.581
#   Settlement pop. 25km:  ΔAICc = 1.66,  weight = 0.253
#   Market gravity:        ΔAICc = 2.51,  weight = 0.166
#   Settlement gravity more clearly preferred once spatial
#   structure is controlled.
#
# Key finding from 25km radius correction:
#   Settlement population at 25km (weight = 0.250/0.253)
#   no longer dominates — all three metrics are competitive
#   without connectivity, and settlement gravity leads with
#   connectivity included. This confirms that the earlier
#   50km result was a radius artefact — the inflated radius
#   was capturing broad coastal human footprint rather than
#   direct SSF pressure. With the ecologically appropriate
#   25km radius, settlement population performs similarly
#   to the gravity metrics, consistent with Cinner et al.
#   (2016) showing travel time components add explanatory
#   power beyond raw population size.
#
# Primary metric: settlement gravity (log_settlement_grav_sc)
#   Justified by: (1) best or competitive performance
#   across both comparisons, (2) clearest support once
#   spatial structure controlled (weight = 0.581),
#   (3) consistency with total biomass Q1 result,
#   (4) theoretically preferred over raw population
#   counts due to explicit travel cost component
#   (Cinner et al. 2016).
#
# Market gravity and settlement population (25km) retained
# for sensitivity analysis.

cat("\n--- Q1: Browser coefficient summary — selected metric ---\n")
summary(b_q1_settgrav)  

# ── Q1: Baseline model coefficients (settlement gravity) ──────
#
# glmmTMB(mean_biomass ~ rugosity_sc + log_settlement_grav_sc
#         + log_chla_sc, family = tweedie)
# n = 54 sites, dispersion = 13.3
#
# Rugosity:           β = +0.440, p = 0.025 *
#   Only significant predictor in the baseline —
#   habitat complexity drives browser biomass
#   independently of pressure or productivity.
#
# Settlement gravity: β = -0.221, p = 0.333 ns
#   Negative but not significant. Weaker pressure
#   signal than total biomass (β = -0.251, p = 0.012)
#   — browser biomass is not strongly structured by
#   the raw exploitation pressure gradient at baseline,
#   consistent with management context being the
#   primary driver of this heavily targeted group.
#
# Chla:               β = -0.277, p = 0.188 ns
#   Non-significant. Retained as baseline control.
#   Consistent with total biomass result (p = 0.208).

# ============================================================
#  Q2 — ARE CONNECTIVITY AND MPA BASELINE DRIVERS?
#
#  Identical structure to total biomass Q2.
#  Baseline: rugosity + [selected pressure metric] + chla
#  Tweedie family throughout.
#  McFadden pseudo-R² used in place of adj. R² —
#  Tweedie precludes standard R² calculation.
# ============================================================

# ── Null ──────────────────────────────────────────────────────
b_null <- glmmTMB(
  mean_biomass ~ 1,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

# ── Baseline (fixed a priori) ─────────────────────────────────
# Settlement gravity as the selected human pressure metric

b_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +  
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

# ── Baseline + connectivity ───────────────────────────────────
# Tests Q2: does connectivity add as main effect?
# Extends Warmuth et al. (2024) — fishing pressure
# included as covariate, biomass not abundance.
b_baseline_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +  
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

# ── Baseline + MPA ────────────────────────────────────────────
# Tests Q2: does formal protection add beyond baseline?
b_baseline_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

# ── Baseline + connectivity + MPA (global additive) ───────────
b_global_additive <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

# ── Q2 comparison ─────────────────────────────────────────────
b_q2_models <- list(
  "Null"                  = b_null,
  "Baseline"              = b_baseline,
  "Baseline + conn"       = b_baseline_conn,
  "Baseline + MPA"        = b_baseline_mpa,
  "Baseline + conn + MPA" = b_global_additive
)

cat("\n--- Q2: Browser model comparison (AICc ranked) ---\n")
print(make_aicc_df(b_q2_models))

# ── McFadden pseudo-R² relative to baseline ───────────────────
cat("\n--- Q2: Browser pseudo-R² relative to baseline ---\n")

null_ll     <- as.numeric(logLik(b_null))
baseline_r2 <- 1 - (as.numeric(logLik(b_baseline)) / null_ll)

b_q2_models %>%
  imap_dfr(~ tibble(
    Model  = .y,
    McF_R2 = round(1 - (as.numeric(logLik(.x)) / null_ll), 3)
  )) %>%
  mutate(
    Delta_R2 = round(McF_R2 - baseline_r2, 3),
    Delta_R2 = ifelse(Model %in% c("Null", "Baseline"),
                      NA, Delta_R2)
  ) %>%
  print()

# ── Q2 results ────────────────────────────────────────────────
#
# AICc comparison:
#   Baseline + MPA:        AICc = 827.63, weight = 0.644 (BEST)
#   Baseline + conn + MPA: ΔAICc = 1.68,  weight = 0.278
#   Baseline + conn:       ΔAICc = 4.40,  weight = 0.072
#   Baseline:              ΔAICc = 9.87,  weight = 0.005
#   Null:                  ΔAICc = 10.94, weight = 0.003
#
# McFadden pseudo-R² relative to baseline (R² = 0.010):
#   Baseline + conn:       ΔR² = +0.010
#   Baseline + MPA:        ΔR² = +0.018
#   Baseline + conn + MPA: ΔR² = +0.020
#
# MPA status is strongly supported for browsers —
# Baseline + MPA is the best-supported model (weight =
# 0.644) and Baseline + conn + MPA is competitive
# (ΔAICc = 1.68, weight = 0.278). Combined weight of
# models including MPA = 0.921 — MPA support is
# unambiguous.
#
# This contrasts sharply with total biomass where MPA
# was not supported (ΔAICc = 3.01, weight = 0.134).
# Browsers are demonstrably more sensitive to formal
# protection than the total fish community — consistent
# with their status as a heavily targeted functional
# group whose local abundance depends strongly on
# management regime.
#
# Connectivity: competitive when combined with MPA
# (Baseline + conn + MPA ΔAICc = 1.68) but not
# independently supported (Baseline + conn ΔAICc = 4.40).
# Connectivity does not add explanatory power beyond
# the baseline alone — it only appears competitive
# when MPA is already in the model. This suggests
# connectivity and MPA may share variance or that
# connectivity's effect on browsers operates through
# or alongside protection rather than independently.
#
# Baseline alone barely improves on null (ΔAICc = 1.07)
# — confirming the Q1 finding that settlement gravity
# is not a strong independent predictor of browser
# biomass. Browser biomass is primarily structured by
# governance and habitat complexity, not by the raw
# exploitation pressure gradient.
#
# McFadden R² values are low overall (max 0.030) —
# Tweedie pseudo-R² is not directly comparable to
# adj. R² from total biomass OLS. AICc is the primary
# comparison metric.
#
# Best-supported model: Baseline + MPA (weight = 0.644)
# Used as reference for Q3.

cat("\n--- Q2: Baseline + MPA coefficients ---\n")
summary(b_baseline_mpa)

cat("\n--- Q2: Baseline + conn + MPA coefficients ---\n")
summary(b_global_additive)

# ── Q2: Baseline + conn + MPA coefficients ────────────────────
#
# rugosity_sc:            β = +0.463, p = 0.004 **
# log_settlement_grav_sc: β = -0.012, p = 0.957 ns
# log_chla_sc:            β = -0.218, p = 0.274 ns
# connectivity_sc:        β = +0.211, p = 0.267 ns
# mpa_statuslow:          β = +0.003, p = 0.996 ns
# mpa_statusmedium:       β = +1.134, p = 0.005 **
#
# Connectivity does not reach significance when MPA is
# included (β = +0.211, p = 0.267) — confirming that
# connectivity's apparent competitiveness in Q2 AICc
# reflects shared variance with MPA rather than an
# independent effect. The medium MPA signal remains
# strong and stable (β = +1.134 vs +1.028 in baseline
# + MPA only) — not diminished by including connectivity.
# Settlement gravity collapses entirely (β = -0.012,
# p = 0.957) — consistent with browser biomass being
# governed by protection and habitat, not the raw
# exploitation gradient.
#
# Coefficient stability across Q2 models:
#   rugosity_sc:      +0.440* → +0.436* → +0.516** → +0.463**
#   settlement_grav:  -0.221  → -0.217  → -0.179   → -0.012
#   mpa_medium:           —       —      +1.028**  → +1.134**
#   connectivity:         —    +0.195   →    —      → +0.211
#
# Rugosity and medium MPA are the two stable, significant
# drivers of browser biomass. All other predictors
# including settlement gravity and connectivity are
# not independently supported.
#
# Best-supported model: Baseline + MPA (weight = 0.644).
# Used as reference for Q3.


# ============================================================
#  SPATIAL AUTOCORRELATION CHECK
#  Residuals from best Q2 model tested for spatial structure.
#  Moran's I on Pearson residuals from glmmTMB Tweedie.
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

cat("\n--- Spatial autocorrelation: browser best Q2 model ---\n")
print(moran.test(residuals(b_baseline_mpa,
                           type = "pearson"), listw5_b))

# ── Spatial autocorrelation result ────────────────────────────
# Moran's I = -0.074, p = 0.802 — no significant spatial
# autocorrelation in browser model residuals.
#
# This contrasts with total biomass (I = 0.140, p = 0.015)
# where weak but significant autocorrelation was detected.
# The absence of spatial structure in browser residuals
# suggests that the predictors in the best Q2 model
# (rugosity + settlement gravity + chla + MPA) adequately
# capture the spatial variation in browser biomass without
# leaving a residual geographic signal.
#
# The negative Moran's I indicates slight spatial
# dispersion — neighbouring sites tend to have more
# dissimilar residuals than expected by chance. This
# is consistent with browser biomass being structured
# by localised management differences (MPA status)
# that vary discretely between adjacent sites rather
# than varying smoothly across space.
#
# Warnings: identical points found — reflects replicate
# survey sites at the same location. 4 sub-graphs
# reflect the four geographically isolated country
# clusters, as per total biomass. Both are expected
# and do not affect interpretation.
#
# No spatial error modelling required for browser
# analysis.

# ============================================================
#  Q3 — DO MANAGEMENT AND CONNECTIVITY MODIFY
#       BROWSER PRESSURE EFFECTS?
#
#  Same three a priori hypotheses as total biomass.
#  Reference model: Baseline + MPA (best Q2 model,
#  weight = 0.644).
#
#  Interpretation note: MPA is supported as a main effect
#  in Q2 (weight = 0.644) — interactions involving MPA
#  are therefore confirmatory if supported. Connectivity
#  not supported as independent main effect (ΔAICc = 4.40
#  without MPA; p = 0.267 within global additive model)
#  — interactions involving connectivity are exploratory.
#
#  Tweedie family throughout.
#  All three run regardless of Q2 main-effect results.
# ============================================================

# Baseline + MPA confirmed as best Q2 model (weight = 0.644)
b_q3_reference <- b_baseline_mpa

# ── H1: MPA effectiveness depends on fishing intensity ────────
# Protected sites only detectable where external pressure
# is low enough for recovery to occur.
b_int_mpa_press <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    mpa_status * log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

# ── H2: Connectivity buffers exploitation effects ─────────────
# Well-connected sites sustain higher browser biomass
# under pressure through larval replenishment.
b_int_conn_press <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    mpa_status +
    connectivity_sc * log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

# ── H3: MPA effectiveness depends on larval supply ────────────
# Protected sites recover faster where connectivity
# is high enough to subsidise recruitment.
b_int_mpa_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc +
    mpa_status * connectivity_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

# ── Q3 comparison ─────────────────────────────────────────────
b_q3_models <- list(
  "Reference (best Q2)"    = b_q3_reference,
  "H1: MPA × pressure"     = b_int_mpa_press,
  "H2: Conn × pressure"    = b_int_conn_press,
  "H3: MPA × connectivity" = b_int_mpa_conn
)

cat("\n--- Q3: Browser interaction comparison ---\n")
print(make_aicc_df(b_q3_models))

cat("\n--- Q3: H1 MPA × pressure coefficients ---\n")
summary(b_int_mpa_press)

# ── Q3 results ────────────────────────────────────────────────
#
# AICc comparison:
#   H1: MPA × pressure:       AICc = 822.86, weight = 0.895 (BEST)
#   Reference (baseline+MPA): ΔAICc = 4.77,  weight = 0.083
#   H3: MPA × connectivity:   ΔAICc = 8.35,  weight = 0.014
#   H2: Conn × pressure:      ΔAICc = 9.33,  weight = 0.008
#
#
# H1 (MPA × pressure) overwhelmingly supported —
# weight = 0.895, ΔAICc = 4.77 vs additive reference.
# MPA was supported as main effect in Q2 — interaction
# is therefore confirmatory, not exploratory.
# H3 and H2 not supported.
#
# H1 coefficients:
#   Rugosity:                    β = +0.534, p < 0.001 ***
#   Chla:                        β = -0.355, p = 0.064 .
#   Settlement gravity (no MPA): β = -0.282, p = 0.220 ns
#   MPA low:                     β = -0.144, p = 0.820 ns
#   MPA medium:                  β = +1.954, p < 0.001 ***
#   MPA low × pressure:          β = +0.870, p = 0.202 ns
#   MPA medium × pressure:       β = +1.838, p = 0.001 **
#
# Medium MPA × pressure significant and in hypothesised
# direction. Effective pressure slope at medium MPA
# sites = -0.282 + 1.838 = +1.556 — browser biomass
# increases with pressure at protected sites, opposite
# to unprotected sites where no significant gradient
# exists. At mean pressure, medium MPA sites have
# e^1.954 = 7.1x higher browser biomass than
# unprotected sites.
#
# Interpretation: at unprotected sites browsers appear
# depleted across the full pressure gradient — no
# detectable variation remains. At medium MPA sites,
# higher surrounding pressure is associated with higher
# biomass, likely reflecting that MPAs in productive
# high-pressure areas provide the strongest refuge
# effect. Unlike the equivalent H1 pattern for total
# biomass (identified as a data structure artefact),
# this result is ecologically interpretable and
# confirmatory given Q2 MPA support.
#
# Overall Q3 conclusion:
# MPA effectiveness depends on external fishing pressure
# for browsers — medium protection reverses the pressure-
# biomass relationship. This is the primary management
# finding for this functional group and contrasts
# directly with total biomass where no interaction
# was supported.

# ── H1 exploratory checks ─────────────────────────────────────
# Verify interaction is genuine before interpreting as
# ecologically meaningful. Mirrors total biomass H1 checks.

# Check 1 — VIFs on interaction model
cat("\n--- H1 browser: VIF check ---\n")
check_collinearity(b_int_mpa_press)

#   All VIFs < 5 — no multicollinearity concern.
#   Interaction term inflation is modest and expected
#   given that interaction terms share variance with
#   their constituent main effects. Does not affect
#   reliability of coefficient estimates.

# Check 2 — medium MPA site distribution on pressure gradient
cat("\n--- H1 browser: medium MPA sites on pressure gradient ---\n")
browser_model_data %>%
  filter(mpa_status == "medium") %>%
  dplyr::select(site, log_settlement_grav_sc,
                mean_biomass) %>%
  arrange(log_settlement_grav_sc) %>%
  print()

# Check 3 — biomass comparison within observed pressure range
cat("\n--- H1 browser: biomass by MPA within observed range ---\n")
browser_model_data %>%
  filter(log_settlement_grav_sc < 0.20) %>%
  group_by(mpa_status) %>%
  summarise(
    n            = n(),
    mean_biomass = round(mean(mean_biomass), 1),
    sd_biomass   = round(sd(mean_biomass),   1),
    min_pressure = round(min(log_settlement_grav_sc), 3),
    max_pressure = round(max(log_settlement_grav_sc), 3),
    .groups = "drop"
  ) %>%
  print()

# Check 4 — raw data plot
ggplot(browser_model_data,
       aes(x = log_settlement_grav_sc,
           y = mean_biomass,
           colour = mpa_status,
           shape  = mpa_status)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_colour_manual(
    values = c("none"   = "#636363",
               "low"    = "#74a9cf",
               "medium" = "#0570b0"),
    labels = c("No MPA", "Low protection", "Medium protection")
  ) +
  scale_shape_manual(
    values = c("none" = 16, "low" = 17, "medium" = 15),
    labels = c("No MPA", "Low protection", "Medium protection")
  ) +
  labs(x      = "log(Settlement gravity) (standardised)",
       y      = "Browser biomass (g)",
       colour = "MPA status",
       shape  = "MPA status") +
  theme_bw(base_size = 12) +
  theme(axis.title        = element_text(face = "bold"),
        legend.position   = c(0.82, 0.85),
        legend.background = element_rect(fill      = "white",
                                         colour    = "grey80",
                                         linewidth = 0.3),
        panel.grid.minor  = element_blank())

# ── Biomass by MPA status within observed pressure range ──────
browser_model_data %>%
  filter(log_settlement_grav_sc < 0.20) %>%
  mutate(MPA = factor(mpa_status,
                      levels = c("none", "low", "medium"),
                      labels = c("No MPA", "Low", "Medium"))) %>%
  ggplot(aes(x = MPA, y = mean_biomass, fill = MPA)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5,
               alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.1, size = 1.8,
              alpha = 0.6, colour = "grey30") +
  scale_fill_manual(values = c("No MPA" = "#bdbdbd",
                               "Low"    = "#74a9cf",
                               "Medium" = "#0570b0")) +
  scale_y_continuous(labels = scales::comma) +
  labs(x       = "MPA status",
       y       = "Browser biomass (g)",
       caption = paste("Sites within observed medium MPA",
                       "pressure range (z < 0.20) only.",
                       "\nn = 18 (none), 3 (low), 17 (medium)")) +
  theme_bw(base_size = 12) +
  theme(
    axis.title         = element_text(face = "bold"),
    legend.position    = "none",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.caption       = element_text(colour = "grey50", size = 8)
  )

# ── Influence diagnostics — DHARMa ────────────────────────────
browser_h1_sim <- simulateResiduals(b_int_mpa_press, n = 1000)
plot(browser_h1_sim)
testOutliers(browser_h1_sim)

# ── Leave-one-out sensitivity: medium MPA × pressure ──────────
loo_results <- map_dfr(seq_len(nrow(browser_model_data)),
                       ~ {
                         fit <- glmmTMB(
                           mean_biomass ~ rugosity_sc + log_chla_sc +
                             mpa_status * log_settlement_grav_sc,
                           family = tweedie(link = "log"),
                           data   = browser_model_data[-.x, ]
                         )
                         coefs <- fixef(fit)$cond
                         tibble(
                           dropped_site     = browser_model_data$site[.x],
                           mpa_status       = browser_model_data$mpa_status[.x],
                           coef_med_x_press = coefs["mpa_statusmedium:log_settlement_grav_sc"]
                         )
                       }
)

cat("\n--- LOO: medium MPA × pressure coefficient stability ---\n")
loo_results %>%
  arrange(coef_med_x_press) %>%
  print(n = 20)

cat("\nFull model coefficient:",
    round(fixef(b_int_mpa_press)$cond[
      "mpa_statusmedium:log_settlement_grav_sc"], 3), "\n")
cat("LOO range:",
    round(min(loo_results$coef_med_x_press), 3), "to",
    round(max(loo_results$coef_med_x_press), 3), "\n")
cat("LOO all positive:",
    all(loo_results$coef_med_x_press > 0), "\n")

# ── H1 exploratory checks ─────────────────────────────────────
#
# DHARMa (n = 1000): KS p = 0.688, dispersion p = 0.202,
#   outlier p = 1.000 — good fit throughout.
#
# VIF: all < 5 — no multicollinearity concern.
#
# Medium MPA distribution: all 17 sites between
#   z = -0.95 and z = 0.17 — positive slope at high
#   pressure is geometric extrapolation. However the
#   biomass elevation within the observed range is
#   genuine: medium MPA mean = 2461g vs no MPA = 878g
#   (2.8x) at equivalent pressure — confirmed by raw
#   data plot and boxplot.
#
# LOO: coefficient positive across all 54 iterations
#   (range 0.644–2.455). Dindini (biomass = 13,386g)
#   has notable influence — dropping it reduces the
#   coefficient from 1.838 to 0.644. Direction robust;
#   magnitude sensitive to this one site.
#
# Conclusion: H1 interaction is ecologically real.
#   Unlike total biomass H1 (no biomass advantage at
#   medium MPAs), browsers show genuine protection
#   benefit. Slope magnitude should be interpreted
#   with caution given dindini's influence.

# ============================================================
#  SENSITIVITY ANALYSIS
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
# Mirrors baseline structure — only pressure metric swapped.
# Purpose: confirm Q1 and Q2 conclusions are not
# metric-dependent for browsers.
# Particularly important here given the higher metric
# uncertainty in Q1 (market gravity vs settlement gravity
# weights more similar than for total biomass).

b_sens_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

b_sens_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_pop_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

cat("\n--- Sensitivity (a): browser alternative metrics ---\n")
cat("Market gravity:\n")
print(summary(b_sens_mktgrav)$coefficients$cond)
cat("\nSettlement population:\n")
print(summary(b_sens_settpop)$coefficients$cond)

# ── Sensitivity (a) results ───────────────────────────────────
#
# Alternative pressure metrics — baseline structure retained,
# only pressure metric substituted. Purpose: confirm Q1
# and Q2 conclusions are not metric-dependent for browsers.
#
# Market gravity (β = +0.191, p = 0.347):
#   Not significant and in the wrong direction — positive
#   coefficient suggests higher market access associated
#   with higher browser biomass, contrary to expectation.
#   Rugosity remains significant (β = +0.540, p = 0.008).
#   Market gravity is a poor proxy for browser exploitation
#   intensity — confirms Q1 decision to use settlement
#   gravity.
#
# Settlement population 25km (β = -0.055, p = 0.761):
#   Not significant and near zero. Much weaker than the
#   old 50km version (β = -0.441, p = 0.009) — confirming
#   that the earlier significant result was a radius
#   artefact driven by inclusion of distant inland
#   populations beyond realistic SSF fishing range.
#   With the ecologically appropriate 25km radius,
#   settlement population adds no explanatory power
#   beyond rugosity and chla. Rugosity remains
#   significant (β = +0.460, p = 0.021).
#
# Both alternative metrics non-significant — settlement
# gravity is the best-supported pressure proxy for
# browsers. Rugosity is the only stable significant
# predictor across all three metrics (β = 0.440–0.540,
# p < 0.05) — habitat complexity signal robust to
# pressure metric choice.
#
# The 25km settlement population result is particularly
# informative: it directly validates the radius
# correction. The 50km metric was capturing broad
# coastal human footprint rather than direct SSF
# pressure — once restricted to the realistic fishing
# range, the population signal disappears entirely.
# This strengthens the case for gravity metrics as
# the appropriate SSF pressure proxy in this system.

# ── (b) Transect-level replication ───────────────────────────
# 43% zeros at transect level — Tweedie required.
# ZI Tweedie tested against standard Tweedie.
# (1 | site) random intercept throughout.

# ── Transect family selection ─────────────────────────────────
b_trans_tw <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_settlement_grav_sc +   
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

b_trans_tw_zi <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_settlement_grav_sc +   
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = browser_transect_data
)

b_trans_res    <- simulateResiduals(b_trans_tw,    n = 1000)
b_trans_res_zi <- simulateResiduals(b_trans_tw_zi, n = 1000)

plot(b_trans_res);    testZeroInflation(b_trans_res)
plot(b_trans_res_zi); testZeroInflation(b_trans_res_zi)

cat("\n--- Sensitivity (b): transect family selection ---\n")
print(make_aicc_df(list(
  "Tweedie"    = b_trans_tw,
  "ZI Tweedie" = b_trans_tw_zi
)))

# ── Transect family selection ─────────────────────────────────
#
# ZI Tweedie: AICc = 2654.34, weight = 0.643
# Tweedie:    ΔAICc = 1.18,   weight = 0.357
#
# ZI Tweedie marginally preferred by AICc but within
# ΔAICc < 2 — genuine uncertainty between the two.
#
# Zero inflation tests (DHARMa, n = 1000):
#   Tweedie:    ratio = 0.965, p = 0.684 — not significant
#   ZI Tweedie: ratio = 0.970, p = 0.738 — not significant
#   Neither model shows evidence of zero inflation —
#   standard Tweedie handles 43% transect-level zeros
#   adequately without an explicit ZI component.
#
# Note: step failure warning on standard Tweedie —
#   DHARMa simulation triggered a convergence warning
#   during GAM smoothing of residuals. This is a
#   DHARMa diagnostic issue, not a model fitting
#   failure — the model itself converged normally.
#   Residual plots should be inspected visually for
#   any systematic patterns.
#
# Proceed: standard Tweedie + (1|site) throughout
# transect-level browser analyses.

# ── Transect Q2 sequence ─────────────────────────────────────
b_trans_null <- glmmTMB(
  transect_browser_biomass ~ 1 + (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

b_trans_baseline <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_settlement_grav_sc +  
    log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

b_trans_conn <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_settlement_grav_sc + 
    log_chla_sc +
    connectivity_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

b_trans_mpa <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_settlement_grav_sc +   
    log_chla_sc +
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

b_trans_global <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_settlement_grav_sc +  
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

cat("\n--- Sensitivity (b): browser transect comparison ---\n")
print(make_aicc_df(list(
  "Null"                  = b_trans_null,
  "Baseline"              = b_trans_baseline,
  "Baseline + conn"       = b_trans_conn,
  "Baseline + MPA"        = b_trans_mpa,
  "Baseline + conn + MPA" = b_trans_global
)))

cat("\n--- Sensitivity (b): baseline coefficients ---\n")
summary(b_trans_baseline)

# ── Sensitivity (b) results ───────────────────────────────────
#
# AICc comparison (transect level):
#   Baseline + MPA:        AICc = 2653.73, weight = 0.441 (BEST)
#   Baseline + conn:       ΔAICc = 1.61,   weight = 0.197
#   Baseline + conn + MPA: ΔAICc = 1.79,   weight = 0.181
#   Baseline:              ΔAICc = 1.87,   weight = 0.174
#   Null:                  ΔAICc = 7.97,   weight = 0.008
#
# Baseline + MPA top-ranked at transect level —
# consistent with site-level Q2 result (weight = 0.644).
# However model selection uncertainty is higher at
# transect level — four models within ΔAICc < 2,
# compared to clearer support at site level. This is
# expected: site-level predictors (MPA status, rugosity,
# settlement gravity) vary between sites, not within
# them, so their signal is diluted by within-site
# transect variation.
#
# Key consistency: Baseline + MPA best-supported at
# both analytical scales — qualitative conclusion
# robust to aggregation level.
#
# Baseline coefficients (REML, n = 243 transects,
# 54 sites):
#   Rugosity:           β = +0.526, p = 0.005 **
#   Settlement gravity: β = -0.223, p = 0.326 ns
#   Chla:               β = -0.505, p = 0.026 *
#
# Random effects:
#   Site-level variance: 1.095 (SD = 1.046)
#   High ICC confirms substantial between-site
#   variation — (1|site) random intercept essential.
#
# Notable: chla reaches significance at transect level
# (β = -0.505, p = 0.026) but not site level
# (β = -0.277, p = 0.188). Greater power at transect
# level (n = 243 vs 54) detects a weaker signal.
# Direction consistent across both levels.
# Effect size substantially larger at transect level —
# may reflect within-site productivity variation
# being detectable at finer resolution.
#
# Rugosity and MPA support consistent across both
# analytical scales — primary Q2 conclusions robust
# to aggregation.


# ============================================================
#  MARGINAL EFFECT PLOTS
#  Generated for significant predictors in best Q2 model
#  (Baseline + MPA, weight = 0.644).
#  Tweedie — predictions on response scale (raw biomass).
#  All non-focal predictors held at 0 (their mean).
# ============================================================

# ── Rugosity effect ───────────────────────────────────────────
rug_grid_b <- data.frame(
  rugosity_sc            = seq(
    min(browser_model_data$rugosity_sc),
    max(browser_model_data$rugosity_sc),
    length.out = 200),
  log_settlement_grav_sc = 0,
  log_chla_sc            = 0,
  mpa_status             = factor("none",
                                  levels = c("none", "low",
                                             "medium"))
)

rug_pred_b <- predict(b_baseline_mpa,
                      newdata = rug_grid_b,
                      se.fit  = TRUE,
                      type    = "response",
                      re.form = NA)

rug_grid_b$fit <- rug_pred_b$fit
rug_grid_b$lwr <- rug_pred_b$fit - 1.96 * rug_pred_b$se.fit
rug_grid_b$upr <- rug_pred_b$fit + 1.96 * rug_pred_b$se.fit

p_b_rugosity <- ggplot(rug_grid_b,
                       aes(x = rugosity_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "#2c7bb6", alpha = 0.15) +
  geom_line(colour = "#2c7bb6", linewidth = 1.1) +
  geom_point(data = browser_model_data,
             aes(x = rugosity_sc, y = mean_biomass),
             colour = "grey40", size = 1.5,
             alpha = 0.5, inherit.aes = FALSE) +
  labs(x = "Rugosity (standardised)",
       y = "Browser biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ── MPA marginal means ────────────────────────────────────────
# From best Q2 model (Baseline + MPA).
# MPA strongly supported in Q2 (weight = 0.644).
# Shows genuine biomass difference — medium MPA sites
# have 2.8x higher browser biomass than unprotected
# sites at equivalent pressure values.

mpa_grid_b <- data.frame(
  mpa_status             = factor(c("none", "low", "medium"),
                                  levels = c("none", "low",
                                             "medium")),
  rugosity_sc            = 0,
  log_settlement_grav_sc = 0,
  log_chla_sc            = 0
)

mpa_pred_b <- predict(b_baseline_mpa,
                      newdata = mpa_grid_b,
                      se.fit  = TRUE,
                      type    = "response",
                      re.form = NA)

mpa_grid_b$fit <- mpa_pred_b$fit
mpa_grid_b$lwr <- mpa_pred_b$fit - 1.96 * mpa_pred_b$se.fit
mpa_grid_b$upr <- mpa_pred_b$fit + 1.96 * mpa_pred_b$se.fit

p_b_mpa <- ggplot(mpa_grid_b,
                  aes(x = mpa_status, y = fit)) +
  geom_hline(yintercept = mpa_grid_b$fit[1],
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
       y = "Browser biomass (g)") +
  theme_bw(base_size = 12) +
  theme(
    axis.title         = element_text(face = "bold"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

gridExtra::grid.arrange(p_b_rugosity, p_b_mpa, ncol = 2)


# ============================================================
#  RESULTS SUMMARY
# ============================================================

browser_results <- tribble(
  ~Question,  ~Best_model,          ~Key_finding,
  "Q1",       "Sett. gravity",      "weight = 0.581 with connectivity; all metrics competitive without",
  "Q2 conn",  "Baseline",           "β = +0.195, p = 0.163, ΔAICc = 4.40 — not supported",
  "Q2 MPA",   "Baseline + MPA",     "weight = 0.644; medium MPA β = +1.028, p = 0.001",
  "Q3 H1",    "MPA × pressure",     "weight = 0.900, ΔAICc = 4.77 — confirmed genuine",
  "Q3 H2",    "Reference",          "ΔAICc = 11.54 — not supported",
  "Q3 H3",    "Reference",          "ΔAICc = 8.35 — not supported"
)

cat("\n--- Browser results summary ---\n")
print(browser_results)


# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()