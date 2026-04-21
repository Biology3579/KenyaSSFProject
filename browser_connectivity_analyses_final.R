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

# Result: [update after running]
# Expected: country RE not supported — consistent with
# total biomass result (ΔAICc = 2.86).
# All browser models fitted without RE throughout.


# ============================================================
#  Q1 — PRESSURE METRIC SELECTION
#
#  Browser-specific Q1 — pressure metric may differ from
#  total biomass. Browsers are a heavily targeted functional
#  group; the metric that best captures exploitation
#  intensity for the total community may not be optimal
#  for a selectively harvested group.
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
  "Settlement gravity" = b_q1_settgrav,
  "Market gravity"     = b_q1_mktgrav,
  "Settlement pop."    = b_q1_settpop
)))

# ── With connectivity control ─────────────────────────────────
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
  "Settlement gravity" = b_q1_settgrav_conn,
  "Market gravity"     = b_q1_mktgrav_conn,
  "Settlement pop."    = b_q1_settpop_conn
)))

# ── Q1 results ────────────────────────────────────────────────
# [Update after running]
#
# Prior result (old analysis, univariate):
#   Market gravity:     AICc = 297.48, weight = 0.552 — best
#   Settlement gravity: ΔAICc = 1.47,  weight = 0.266
#   Settlement pop.:    ΔAICc = 2.22,  weight = 0.182
#
# Market gravity previously preferred over settlement gravity
# for browsers — contrasts with total biomass where settlement
# gravity was clearly preferred (weight = 0.878).
#
# This may reflect browsers being a higher-value commercial
# target (parrotfish), making market access a stronger
# predictor of exploitation intensity for this group than
# for the total community. However, model selection
# uncertainty was high (weights 0.552 vs 0.266) — check
# whether result is consistent across both model versions
# (with and without connectivity) before confirming metric.
#
# Primary metric: [update after running]
# Used throughout Q2, Q3, and sensitivity analyses.

cat("\n--- Q1: Browser coefficient summary — selected metric ---\n")
summary(b_q1_settgrav)  # update to selected metric


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
# Update pressure metric after Q1 decision
# Placeholder: settlement gravity — replace if Q1 selects
# different metric
b_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +   # ← update after Q1
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
    log_settlement_grav_sc +   # ← update after Q1
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

# ── Baseline + MPA ────────────────────────────────────────────
# Tests Q2: does formal protection add beyond baseline?
b_baseline_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +   # ← update after Q1
    log_chla_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

# ── Baseline + connectivity + MPA (global additive) ───────────
b_global_additive <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +   # ← update after Q1
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

# ── Coefficient summaries ─────────────────────────────────────
cat("\n--- Q2: Baseline coefficients ---\n")
summary(b_baseline)

cat("\n--- Q2: Baseline + connectivity coefficients ---\n")
summary(b_baseline_conn)

cat("\n--- Q2: Baseline + MPA coefficients ---\n")
summary(b_baseline_mpa)

# ── Q2 results ────────────────────────────────────────────────
# [Update after running]
#
# Key contrast expected from prior analysis:
#   MPA was strongly supported for browsers (old Local + MPA
#   weight = 0.693) — very different from total biomass
#   where MPA was not supported (ΔAICc = 3.01).
#   Connectivity also showed moderate support for browsers
#   (old ΔAICc = 3.32) vs no support for total biomass.
#
# Note: McFadden R² not directly comparable to adj. R²
# from total biomass. Use AICc as primary comparison metric
# across functional groups. McFadden R² used only for
# within-group ΔR² increments.
#
# Best-supported model: [update after running]
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
print(moran.test(residuals(b_baseline,
                           type = "pearson"), listw5_b))

# Note: use type = "pearson" for glmmTMB Tweedie residuals.
# Moran's I interpretation as per total biomass — weak
# signal expected given country RE not supported and
# predictors absorb between-cluster structure.


# ============================================================
#  Q3 — DO MANAGEMENT AND CONNECTIVITY MODIFY
#       BROWSER PRESSURE EFFECTS?
#
#  Same three a priori hypotheses as total biomass.
#  Reference model: best Q2 model.
#  Tweedie family throughout.
#  All three run regardless of Q2 main-effect results.
#  Interpretation conditional on Q2.
# ============================================================

# Reference: best Q2 model — update after Q2
b_q3_reference <- b_baseline   # replace after Q2

# ── H1: MPA effectiveness depends on fishing intensity ────────
b_int_mpa_press <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    mpa_status * log_settlement_grav_sc,  # ← update metric
  family = tweedie(link = "log"),
  data   = browser_model_data
)

# ── H2: Connectivity buffers exploitation effects ─────────────
b_int_conn_press <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc * log_settlement_grav_sc,  # ← update metric
  family = tweedie(link = "log"),
  data   = browser_model_data
)

# ── H3: MPA effectiveness depends on larval supply ────────────
b_int_mpa_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc +        # ← update metric
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

cat("\n--- Q3: H2 Connectivity × pressure coefficients ---\n")
summary(b_int_conn_press)

cat("\n--- Q3: H3 MPA × connectivity coefficients ---\n")
summary(b_int_mpa_conn)

# ── Q3 results ────────────────────────────────────────────────
# [Update after running]
#
# Prior analysis showed no interaction supported for browsers
# when using the additive Local + MPA as reference:
#   MPA × connectivity: ΔAICc = 2.98
#   MPA × pressure:     ΔAICc = 3.47
#   Conn × pressure:    ΔAICc = 5.92
#
# This may change with the updated baseline and metric.
# Interpret interactions conditional on Q2 support —
# if MPA and/or connectivity supported in Q2, interactions
# are confirmatory; if not, exploratory only.

# ── Overall Q3 conclusion ─────────────────────────────────────
# [Update after running]


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

# Key checks:
# (1) Rugosity direction and significance stable?
# (2) Pressure direction consistent across metrics?
# (3) Do Q2 conclusions about MPA and connectivity
#     hold regardless of metric?

# ── (b) Transect-level replication ───────────────────────────
# 43% zeros at transect level — Tweedie required.
# ZI Tweedie tested against standard Tweedie.
# (1 | site) random intercept throughout.

# ── Transect family selection ─────────────────────────────────
b_trans_tw <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_settlement_grav_sc +   # ← update after Q1
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

b_trans_tw_zi <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_settlement_grav_sc +   # ← update after Q1
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

# Prior result: standard Tweedie selected over ZI Tweedie
# (ZI test not significant for either; ZI Tweedie AICc
# marginally lower but within ΔAICc < 2).
# Proceed with standard Tweedie + (1|site) if confirmed.

# ── Transect Q2 sequence ─────────────────────────────────────
b_trans_null <- glmmTMB(
  transect_browser_biomass ~ 1 + (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

b_trans_baseline <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_settlement_grav_sc +   # ← update after Q1
    log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

b_trans_conn <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_settlement_grav_sc +   # ← update after Q1
    log_chla_sc +
    connectivity_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

b_trans_mpa <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_settlement_grav_sc +   # ← update after Q1
    log_chla_sc +
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

b_trans_global <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_settlement_grav_sc +   # ← update after Q1
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
# [Update after running]
# Key check: does the Q2 model ranking replicate at
# transect level? Particularly whether MPA support
# persists — if it does, the site-level finding is
# not an artefact of aggregation.


# ============================================================
#  MARGINAL EFFECT PLOTS
#  Generated for significant predictors in best Q2 model.
#  Tweedie — predictions on response scale (raw biomass).
# ============================================================

# ── Rugosity effect ───────────────────────────────────────────
rug_grid_b <- data.frame(
  rugosity_sc            = seq(
    min(browser_model_data$rugosity_sc),
    max(browser_model_data$rugosity_sc),
    length.out = 200),
  log_settlement_grav_sc = 0,   # ← update metric name
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
# Generated from best Q2 model including MPA.
# Update model object after Q2 result confirmed.

mpa_grid_b <- data.frame(
  mpa_status             = factor(c("none", "low", "medium"),
                                  levels = c("none", "low",
                                             "medium"),
                                  labels = c("None", "Low",
                                             "Medium")),
  rugosity_sc            = 0,
  log_settlement_grav_sc = 0,   # ← update metric name
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
#  Quick reference — update after running all models.
# ============================================================

browser_results <- tribble(
  ~Question,  ~Best_model,  ~Key_finding,
  "Q1",       "[update]",   "[update after running]",
  "Q2 conn",  "[update]",   "[update after running]",
  "Q2 MPA",   "[update]",   "[update after running]",
  "Q3 H1",    "[update]",   "[update after running]",
  "Q3 H2",    "[update]",   "[update after running]",
  "Q3 H3",    "[update]",   "[update after running]"
)

cat("\n--- Browser results summary ---\n")
print(browser_results)


# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()