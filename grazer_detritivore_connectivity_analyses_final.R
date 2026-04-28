# ============================================================
#  DRIVERS OF GRAZER-DETRITIVORE BIOMASS
#  Chapter 1 — Functional Group Analysis: Grazers/Detritivores
#
#  Analytical framework mirrors browser and corallivore
#  analyses. Structure: Q1 (pressure) → Q2 (connectivity)
#  → Q3 (MPA)
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
#       Connectivity x pressure interaction only testable
#       if pressure is supported in Q1.
#
#  Q3 — Formal protection
#       Does MPA status explain additional variation beyond
#       the best Q2 model?
#       MPA x pressure not testable if pressure not in
#       model. MPA x connectivity explored if MPA
#       supported as main effect.
#
#  Baseline model (fixed a priori, never tested):
#       log(biomass) ~ rugosity_sc + log_chla_sc
#
#  Key difference from browsers/corallivores:
#       Check zero proportion at site level (below).
#       If zero-free → lm() on log_mean_biomass.
#       If zeros present → Tweedie glmmTMB.
#       Family confirmed in family selection section.
#
#  Sensitivity analyses:
#       (a) Alternative pressure metrics
#       (b) Transect-level replication
# ============================================================

source(here::here("data_preparation.R"))


# ============================================================
#  DATA AGGREGATION
# ============================================================

grazer_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_grazer_biomass = sum(
      ifelse(trophic_group %in% c("grazers", "detritivores",
                                  "grazer-detritivores"),
             tot_wt_g, 0),
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

cat("Grazer transects:", nrow(grazer_transects), "\n")
cat("Sites:",            n_distinct(grazer_transects$site), "\n")

# ── Site-level dataset ────────────────────────────────────────
grazer_model_data <- grazer_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_grazer_biomass,
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
                        ordered = FALSE),
    log_mean_biomass = log(mean_biomass)
  )

cat("\nGrazer model data:", nrow(grazer_model_data), "sites\n")

# ── Data checks ───────────────────────────────────────────────
grazer_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                log_chla_sc, connectivity_sc, mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

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

# ── Transect-level dataset ────────────────────────────────────
grazer_transect_data <- grazer_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(grazer_transect_data$transect_grazer_biomass == 0),
    "/", nrow(grazer_transect_data),
    "(", round(mean(grazer_transect_data$transect_grazer_biomass == 0),
               3), ")\n")


# ============================================================
#  MODEL FAMILY SELECTION
#  Run on baseline model.
#  If zero-free: lm() on log_mean_biomass (Gaussian log).
#  If zeros present: Tweedie glmmTMB.
#  Check zero proportion above before proceeding.
# ============================================================

# ── Gaussian log — baseline ───────────────────────────────────
grazer_lm_baseline <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_chla_sc,
  data = grazer_model_data
)

par(mfrow = c(2, 2))
plot(grazer_lm_baseline, main = "Gaussian log — baseline")
par(mfrow = c(1, 1))

# ── Family selection decision ─────────────────────────────────
#
# Gaussian log: SELECTED — lm() on log_mean_biomass.
#   Zero-free response (0 site-level zeros) — log
#   transformation valid without offset.
#   Diagnostics (baseline model): TBC — update after
#   reviewing plots above.
#
# Tweedie: REJECTED — convergence failure on baseline
#   model (non-positive-definite Hessian). Consistent
#   with total biomass where Tweedie was not selected.
#   Zero-free response makes Gaussian log the correct
#   and sufficient choice.
#
#  Proceed: lm() on log_mean_biomass throughout all
#  grazer-detritivore site-level analyses.

# ============================================================
#  RANDOM EFFECT STRUCTURE
#  Tested on baseline model — same rationale as all
#  previous analyses.
# ============================================================

g_re_null <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_chla_sc,
  data = grazer_model_data
)

g_re_ecoregion <- lmer(
  log_mean_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | ecoregion),
  data = grazer_model_data,
  REML = FALSE
)

cat("\n--- Grazer RE structure: ecoregion ---\n")
print(make_aicc_df(list(
  "No RE"           = g_re_null,
  "(1 | ecoregion)" = g_re_ecoregion
)))

# Random effect structure: ecoregion
# Tested on baseline model to avoid circularity.
# (1 | ecoregion): AICc = 121.36, weight = 0.813 (BEST)
# No RE:           DAICc = 2.93,  weight = 0.187
#
# Despite stronger AICc support, RE not pursued:
# (1) Only 4 ecoregions — insufficient group-level
#     replication for reliable variance component
#     estimation (Gelman & Hill 2007).
# (2) Severely uneven group sizes:
#     Kenya-Tanzania north:       2 sites
#     Comoros:                    8 sites
#     Madagascar:                 9 sites
#     Tanzania south-Mozambique: 35 sites
#     Kenya-Tanzania north (n = 2) cannot support a
#     meaningful group-level intercept estimate.
#
# Between-ecoregion variation acknowledged as a
# limitation — more prominent for this functional
# group than total biomass. May reflect differential
# grazer-detritivore communities across ecoregions
# beyond what rugosity and chla capture.
# All models fitted as lm() throughout.

# ============================================================
#  VIF CHECKS - can be added to total bioamss only or data prep. 
#   VIF identical to total biomass analysis — same scaled
#   predictor set used across all functional group analyses.
#   Baseline VIF = 1.000, max VIF = 1.49 at baseline +
#   pressure + connectivity. No multicollinearity concern.
#   See 00_data_preparation.R for full VIF output.
# ============================================================

cat("\n--- VIF: baseline ---\n")
vif(lm(log_mean_biomass ~ rugosity_sc +
         log_chla_sc,
       data = grazer_model_data))

cat("\n--- VIF: baseline + pressure ---\n")
vif(lm(log_mean_biomass ~ rugosity_sc +
         log_chla_sc +
         log_settlement_grav_sc,
       data = grazer_model_data))

cat("\n--- VIF: baseline + pressure + connectivity ---\n")
vif(lm(log_mean_biomass ~ rugosity_sc +
         log_chla_sc +
         log_settlement_grav_sc +
         connectivity_sc,
       data = grazer_model_data))

cat("\n--- VIF: baseline + connectivity + MPA ---\n")
vif(lm(log_mean_biomass ~ rugosity_sc +
         log_chla_sc +
         connectivity_sc +
         mpa_status,
       data = grazer_model_data))

# VIF checks:
#   Baseline: all VIF = 1.000 — orthogonal
#   Baseline + pressure: max VIF = 1.44 — acceptable
#   Baseline + pressure + connectivity: max VIF = 1.49
#   All VIFs < 2.0 — no multicollinearity concern.
#   Consistent with total biomass VIF results.

# ============================================================
#  Q1 — HUMAN PRESSURE
#
#  Does human pressure explain variation in grazer-
#  detritivore biomass beyond local ecological context,
#  and which spatial metric best captures SSF exploitation?
#
#  Two steps:
#  Step 1 — AICc comparison: does any metric outperform
#            ecological baseline?
#  Step 2 — Coefficient check: direction and significance
#            of selected metric (or direction checks only
#            if baseline best supported).
# ============================================================

g_baseline <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_chla_sc,
  data = grazer_model_data
)

# ── Q1 Step 1: Metric comparison ─────────────────────────────
g_q1_settgrav <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc,
  data = grazer_model_data
)

g_q1_mktgrav <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_market_gravity_sc,
  data = grazer_model_data
)

g_q1_settpop <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_pop_sc,
  data = grazer_model_data
)

cat("\n--- Q1 Step 1: Grazer metric comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                      = g_baseline,
  "Baseline + settlement gravity" = g_q1_settgrav,
  "Baseline + market gravity"     = g_q1_mktgrav,
  "Baseline + settlement pop."    = g_q1_settpop
)))

# ── Q1 Step 2: Coefficient check ─────────────────────────────
cat("\n--- Q1 Step 2: Baseline coefficients ---\n")
summary(g_baseline)

# ── Q1: Baseline confint ─────────────────────────────────────
cat("\n--- Grazer Q1: baseline confint ---\n")
print(confint(g_baseline))


cat("\n--- Q1: Pressure metric direction checks ---\n")
cat("Settlement gravity:\n")
print(summary(g_q1_settgrav)$coefficients["log_settlement_grav_sc", ])
cat("\nMarket gravity:\n")
print(summary(g_q1_mktgrav)$coefficients["log_market_gravity_sc", ])
cat("\nSettlement population:\n")
print(summary(g_q1_settpop)$coefficients["log_settlement_pop_sc", ])

# Q1 results:
#   Baseline:           AICc = 124.30, weight = 0.497 (BEST)
#   Market gravity:     DAICc = 1.81,  weight = 0.201
#   Settlement gravity: DAICc = 2.34,  weight = 0.154
#   Settlement pop.:    DAICc = 2.43,  weight = 0.147
#
#   Baseline best supported — no pressure metric
#   outperforms ecological context for grazer-detritivores.
#   All three metrics positive and non-significant —
#   consistent absence of a pressure signal and all
#   in the wrong direction relative to exploitation
#   hypothesis.
#
# Baseline coefficients:
#   Rugosity: b = +0.245, p = 0.017 * SIGNIFICANT
#     Habitat structural complexity is a significant
#     positive driver — consistent with grazers
#     depending on complex reef structure for foraging
#     substrate and shelter.
#   Chla:     b = -0.069, p = 0.506 ns
#     Non-significant, retained as baseline control.
#
# Pressure metric direction checks:
#   Settlement gravity: b = +0.037, p = 0.766 ns
#   Market gravity:     b = +0.082, p = 0.449 ns
#   Settlement pop.:    b = +0.005, p = 0.960 ns
#
#   All three positive — opposite to expected direction
#   if pressure depletes grazer-detritivores. Effect
#   sizes near zero throughout. This is biologically
#   interpretable — grazer-detritivores are less
#   selectively targeted by SSFs than browsers or
#   piscivores, and may even benefit indirectly from
#   fishing through reduced predation pressure or
#   competitive release. Absence of pressure signal
#   robust across all three metrics.
#
# Best Q1 model: g_baseline (rugosity + chla)
# No pressure term carried forward.
g_best_q1 <- g_baseline


# ============================================================
#  Q2 — LARVAL CONNECTIVITY
#
#  Does connectivity explain additional variation beyond
#  the ecological baseline?
#
#  Main effect only — connectivity x pressure interaction
#  not testable because pressure is not a supported
#  predictor for grazer-detritivores (Q1: baseline best
#  supported, weight = 0.497). Interaction requires
#  pressure as a main effect first.
#
#  lm() throughout — consistent with family selection.
# ============================================================

g_q2_conn <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc,
  data = grazer_model_data
)

cat("\n--- Q2: Grazer connectivity main effect ---\n")
print(make_aicc_df(list(
  "Best Q1"        = g_best_q1,
  "Best Q1 + conn" = g_q2_conn
)))

cat("\n--- Q2: Connectivity coefficients ---\n")
summary(g_q2_conn)

# ── Q2: Connectivity range and fold difference ───────────────
cat("\n--- Grazer Q2: connectivity range and fold difference ---\n")
conn_range_g <- range(grazer_model_data$connectivity_sc, na.rm = TRUE)
cat("Range:", conn_range_g, "\n")
conn_span_g  <- diff(conn_range_g)
cat(sprintf("Span: %.3f SD units\n", conn_span_g))
b_conn_g     <- -0.196
fold_conn_g  <- exp(abs(b_conn_g * conn_span_g))
cat(sprintf("Fold difference (low vs high connectivity): %.2fx\n", fold_conn_g))

# ── Q2: Connectivity confint ──────────────────────────────────
cat("\n--- Grazer Q2: confint ---\n")
print(confint(g_q2_conn))


# Q2 results:
#   Best Q1 + conn: AICc = 122.82, weight = 0.677 (BEST)
#   Best Q1:        DAICc = 1.48,  weight = 0.323
#   Connectivity marginally supported (weight = 0.677,
#   DAICc = 1.48 — below conventional DAICc > 2 threshold
#   but weight clearly favours connectivity model).
#
#   Connectivity: b = -0.196, p = 0.058 .
#     Marginal negative effect — well-connected sites
#     tend to have lower grazer-detritivore biomass.
#     Same direction as corallivores (b = -0.221,
#     p = 0.015) — contrasts with browsers (positive).
#     Below conventional significance threshold (p = 0.05)
#     but consistent in direction with corallivores.
#   Rugosity:     b = +0.241, p = 0.016 * — stable
#   Chla:         b = -0.125, p = 0.235 ns — retained
#
#   Interpretation: negative connectivity signal for
#   grazer-detritivores mirrors corallivores and contrasts
#   with browsers. Well-connected sites may experience
#   higher predator or competitor recruitment, or
#   connectivity may reflect site accessibility and
#   associated exploitation pressure not captured by
#   settlement gravity. Effect is marginal — report
#   with caution. DAICc = 1.48 means genuine model
#   selection uncertainty between connectivity and
#   baseline models.
#
#   Given marginal support (weight = 0.677, p = 0.058),
#   connectivity carried forward into Q3 as best
#   supported model but result reported as suggestive.

# ── Best Q2 model ─────────────────────────────────────────────
g_best_q2 <- g_q2_conn  # rugosity + chla + connectivity

# ============================================================
#  Q3 — FORMAL PROTECTION
#
#  Does MPA status explain additional variation beyond
#  the connectivity baseline?
#
#  MPA x pressure not testable — pressure not supported
#  in Q1. MPA x connectivity not tested — MPA not
#  supported as main effect.
#
#  lm() throughout.
# ============================================================

g_q3_mpa <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  data = grazer_model_data
)

cat("\n--- Q3: Grazer MPA main effect ---\n")
print(make_aicc_df(list(
  "Best Q2"       = g_best_q2,
  "Best Q2 + MPA" = g_q3_mpa
)))

cat("\n--- Q3: MPA coefficients ---\n")
summary(g_q3_mpa)

# ── Q3: MPA placement check ──────────────────────────────────
cat("\n--- Grazer Q3: pressure and connectivity by MPA ---\n")
grazer_model_data %>%
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
#   Best Q2:       AICc = 122.82, weight = 0.925 (BEST)
#   Best Q2 + MPA: DAICc = 5.03,  weight = 0.075
#   MPA clearly not supported — connectivity model
#   strongly preferred.
#
#   Low MPA:    b = +0.062, p = 0.856 ns
#   Medium MPA: b = -0.053, p = 0.822 ns
#   Both near zero and non-significant. Raw means
#   confirm no protection effect: none = 3818g,
#   low = 2187g, medium = 2750g — unprotected sites
#   have highest mean biomass, consistent with MPA
#   placement in lower-biomass areas or low MPA sites'
#   high-connectivity characteristics (z = 0.66 to max).
#
#   Connectivity: b = -0.199, p = 0.096 . — marginally
#     weakens slightly once MPA included but direction
#     consistent.
#   Rugosity: b = +0.241, p = 0.019 * — stable.
#
#   MPA x connectivity interaction not tested —
#   MPA not supported as main effect (DAICc = 5.03).
#
#   Consistent with total biomass (DAICc = 3.01) and
#   corallivores (DAICc = 0.39). Grazer-detritivore
#   biomass is not structured by formal protection
#   designation in this system.
#
# Best Q3 model: rugosity + chla + connectivity
g_best_q3 <- g_best_q2

cat("\n--- Grazer: predicted vs observed ---\n")
pred_g <- predict(g_best_q3)
obs_g  <- grazer_model_data$log_mean_biomass
cat(sprintf("Pearson r: %.3f\n", cor(pred_g, obs_g)))

# ============================================================
#  SPATIAL AUTOCORRELATION
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

cat("\n--- Spatial autocorrelation: grazer best model ---\n")
print(moran.test(residuals(g_best_q3), listw5_g))

# Spatial autocorrelation: grazer best model
# (rugosity + chla + connectivity)
# Moran's I = 0.210, p = 0.001 — significant spatial
# autocorrelation in residuals.
#
# Stronger than total biomass (I = 0.140, p = 0.015)
# and contrasts with browsers (I = -0.078, p = 0.811)
# and corallivores (I = -0.021, p = 0.509).
#
# Spatial error modelling not pursued — same reasons
# as total biomass:
# (1) Only 4 ecoregions with severely uneven group
#     sizes (n = 2, 8, 9, 35) — RE not estimable.
# (2) Discontinuous sampling design — k-NN weights
#     bridge geographically isolated country clusters,
#     not reflecting true within-region covariance
#     structure (Kissling & Carl 2008).
#
# Residual spatial structure more pronounced for
# grazer-detritivores than other groups — consistent
# with the stronger ecoregion RE support at baseline
# (DAICc = 2.93, weight = 0.813). Between-ecoregion
# variation in grazer-detritivore communities is not
# fully captured by rugosity, chla, and connectivity.
# Acknowledged as a limitation — may inflate type I
# error rates for connectivity coefficient (b = -0.196,
# p = 0.058). Marginal connectivity result should be
# interpreted with additional caution given residual
# spatial structure.


# ============================================================
#  SENSITIVITY ANALYSIS
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
# Confirms Q1 null result is not metric-dependent.
# All three metrics expected non-significant and positive
# — consistent with Q1 direction checks.

g_sens_mktgrav <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_market_gravity_sc,
  data = grazer_model_data
)

g_sens_settpop <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_pop_sc,
  data = grazer_model_data
)

cat("\n--- Sensitivity (a): grazer alternative metrics ---\n")
cat("Market gravity:\n")
print(summary(g_sens_mktgrav)$coefficients)
cat("\nSettlement population:\n")
print(summary(g_sens_settpop)$coefficients)

# Sensitivity (a) results:
#   Market gravity:      b = +0.082, p = 0.449 ns
#   Settlement pop.:     b = +0.005, p = 0.960 ns
#
#   Both positive and non-significant — consistent with
#   Q1 direction checks (settlement gravity b = +0.037,
#   market gravity b = +0.082, settlement pop. b = +0.005).
#   All three metrics positive across Q1 and sensitivity —
#   absence of pressure signal is not metric-dependent
#   and direction is consistently opposite to exploitation
#   hypothesis. Rugosity significant and stable across
#   all three metrics (b = +0.245-0.262, p < 0.05).
#   Q1 null conclusion robust.


# ── (b) Transect-level replication ───────────────────────────
# 4 zero-biomass transects (1.6%) produce -Inf on log
# scale — filtered before transect-level analysis.
# Consistent with site-level lm() on log_mean_biomass.
# All 54 sites retained.

grazer_transect_data_nz <- grazer_transect_data %>%
  filter(transect_grazer_biomass > 0) %>%
  mutate(log_transect_biomass = log(transect_grazer_biomass))

cat("Transects retained:", nrow(grazer_transect_data_nz),
    "of", nrow(grazer_transect_data), "\n")

g_trans_null <- lmer(
  log_transect_biomass ~ 1 + (1 | site),
  data = grazer_transect_data_nz,   # _nz
  REML = TRUE
)

g_trans_baseline <- lmer(
  log_transect_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | site),
  data = grazer_transect_data_nz,   # _nz
  REML = TRUE
)

g_trans_pressure <- lmer(
  log_transect_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc +
    (1 | site),
  data = grazer_transect_data_nz,   # _nz
  REML = TRUE
)

g_trans_conn <- lmer(
  log_transect_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc +
    (1 | site),
  data = grazer_transect_data_nz,   # _nz
  REML = TRUE
)

g_trans_mpa <- lmer(
  log_transect_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | site),
  data = grazer_transect_data_nz,   # _nz
  REML = TRUE
)

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

cat("\n--- Sensitivity (b): grazer transect comparison ---\n")
print(make_aicc_df(list(
  "Null"             = g_trans_null_ml,
  "Baseline"         = g_trans_baseline_ml,
  "Baseline + press" = g_trans_pressure_ml,
  "Baseline + conn"  = g_trans_conn_ml,
  "Best + MPA"       = g_trans_mpa_ml
)))

# Sensitivity (b) results:
#   Baseline + conn: AICc = 625.50, weight = 0.453 (BEST)
#   Baseline:        DAICc = 0.71,  weight = 0.317
#   Baseline + press: DAICc = 2.39, weight = 0.137
#   Best + MPA:      DAICc = 4.14,  weight = 0.057
#   Null:            DAICc = 5.04,  weight = 0.036
#
#   Connectivity model top-ranked at transect level —
#   consistent with site-level Q2 (weight = 0.677).
#   Model selection uncertainty higher at transect level
#   (connectivity weight = 0.453 vs 0.677 at site level)
#   but direction of support consistent across scales.
#   Pressure not supported (DAICc = 2.39) — consistent
#   with Q1 null result.
#   MPA not supported (DAICc = 4.14) — consistent with
#   Q3 (DAICc = 5.03).
#
#   Qualitative conclusions robust to aggregation level.
#   Connectivity is the primary driver at both scales,
#   though marginal at both (site: DAICc = 1.48,
#   weight = 0.677; transect: DAICc = 0.71, weight = 0.453).

cat("\n--- Sensitivity (b): connectivity model coefficients ---\n")
summary(g_trans_conn)

# ── ICC calculation ───────────────────────────────────────────
vc_g <- as.data.frame(VarCorr(g_trans_conn))
site_var_g     <- vc_g$vcov[vc_g$grp == "site"]
residual_var_g <- vc_g$vcov[vc_g$grp == "Residual"]
icc_g          <- site_var_g / (site_var_g + residual_var_g)

cat(sprintf("\nICC = %.3f — %.1f%% of variance attributable
to between-site differences\n", icc_g, icc_g * 100))

# Sensitivity (b) — connectivity model coefficients
# (REML, n = 239 transects after removing 4 zeros,
# 54 sites):
#   Connectivity: b = -0.146, t = -1.637
#     Consistent direction with site-level (b = -0.196)
#     — negative signal robust to aggregation level.
#     Weaker at transect level, consistent with marginal
#     site-level result (p = 0.058).
#   Rugosity:     b = +0.235, t = 2.729 — stable
#   Chla:         b = -0.147, t = -1.568 ns — retained
#
# Random effects:
#   Site variance:    0.256 (SD = 0.506)
#   Residual variance: 0.613 (SD = 0.783)
#   ICC = 0.306 — 30.6% of variance attributable to
#   between-site differences beyond fixed predictors.
#   Confirms (1 | site) random intercept is justified.
#   Higher ICC than total biomass (0.205) — consistent
#   with stronger ecoregion RE support for grazers
#   at baseline (DAICc = 2.93).
#
# Overall sensitivity (b) conclusion:
#   Connectivity negative at both site and transect
#   levels — direction consistent, magnitude stable.
#   Marginal support at both scales confirms suggestive
#   rather than conclusive interpretation. Pressure
#   and MPA not supported at either level — consistent
#   with Q1 and Q3 conclusions. Site-level aggregation
#   does not alter qualitative inference.

# ============================================================
#  MARGINAL EFFECT PLOTS
#  Best model: rugosity + chla + connectivity (g_best_q3)
#  MPA not supported — no MPA plot.
#  Two plots: rugosity and connectivity.
#  Response on log scale — lm() on log_mean_biomass.
# ============================================================

best_model_g <- g_best_q3  # rugosity + chla + connectivity

# ── Rugosity effect ───────────────────────────────────────────
rug_grid_g <- data.frame(
  rugosity_sc     = seq(
    min(grazer_model_data$rugosity_sc),
    max(grazer_model_data$rugosity_sc),
    length.out = 200),
  log_chla_sc     = 0,
  connectivity_sc = 0
)

rug_pred_g     <- predict(best_model_g,
                          newdata = rug_grid_g,
                          se.fit  = TRUE)
rug_grid_g$fit <- rug_pred_g$fit
rug_grid_g$lwr <- rug_pred_g$fit - 1.96 * rug_pred_g$se.fit
rug_grid_g$upr <- rug_pred_g$fit + 1.96 * rug_pred_g$se.fit

p_g_rugosity <- ggplot(rug_grid_g,
                       aes(x = rugosity_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "#2c7bb6", alpha = 0.15) +
  geom_line(colour = "#2c7bb6", linewidth = 1.1) +
  geom_point(data = grazer_model_data,
             aes(x = rugosity_sc, y = log_mean_biomass),
             colour = "grey40", size = 1.5,
             alpha  = 0.5, inherit.aes = FALSE) +
  labs(x = "Rugosity (standardised)",
       y = "log(Grazer-detritivore biomass)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())


# ── Connectivity effect ───────────────────────────────────────
# Connectivity marginally supported in Q2 (weight = 0.677,
# DAICc = 1.48, b = -0.196, p = 0.058) — negative direction
# consistent with corallivores (b = -0.221, p = 0.015).
# Reported as suggestive given DAICc < 2 and significant
# residual spatial autocorrelation (Moran's I = 0.210,
# p = 0.001). Direction consistent at transect level
# (b = -0.146, t = -1.637).

conn_grid_g <- data.frame(
  connectivity_sc = seq(
    min(grazer_model_data$connectivity_sc),
    max(grazer_model_data$connectivity_sc),
    length.out = 200),
  rugosity_sc = 0,
  log_chla_sc = 0
)

conn_pred_g      <- predict(best_model_g,
                            newdata = conn_grid_g,
                            se.fit  = TRUE)
conn_grid_g$fit  <- conn_pred_g$fit
conn_grid_g$lwr  <- conn_pred_g$fit - 1.96 * conn_pred_g$se.fit
conn_grid_g$upr  <- conn_pred_g$fit + 1.96 * conn_pred_g$se.fit

p_g_conn <- ggplot(conn_grid_g,
                   aes(x = connectivity_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "#2c7bb6", alpha = 0.15) +
  geom_line(colour = "#2c7bb6", linewidth = 1.1) +
  geom_point(data = grazer_model_data,
             aes(x = connectivity_sc, y = log_mean_biomass),
             colour = "grey40", size = 1.5,
             alpha  = 0.5, inherit.aes = FALSE) +
  labs(x = "Connectivity (standardised)",
       y = "log(Grazer-detritivore biomass)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())


# ── Arrange plots ─────────────────────────────────────────────
gridExtra::grid.arrange(p_g_rugosity, p_g_conn, ncol = 2)

# jpeg("grazer_marginal_effects.jpg",
#      width = 22, height = 11, units = "cm", res = 300)
# gridExtra::grid.arrange(p_g_rugosity, p_g_conn, ncol = 2)
# dev.off()


# ============================================================
#  RESULTS SUMMARY
# ============================================================

cat("\n--- Grazer-detritivore results summary ---\n")
tribble(
  ~Question,   ~Result,            ~Key_finding,
  "Q1",        "Baseline best",    "weight = 0.497, all metrics positive ns — no pressure signal",
  "Q2 conn",   "Suggestive",       "weight = 0.677, DAICc = 1.48, b = -0.196, p = 0.058",
  "Q3 MPA",    "Not supported",    "DAICc = 5.03, weight = 0.075, connectivity model retained",
  "Spatial",   "Significant",      "Moran's I = 0.210, p = 0.001 — strongest of all groups",
  "Sens (a)",  "Consistent",       "market b = +0.082 ns, pop b = +0.005 ns — all positive, Q1 null robust",
  "Sens (b)",  "Consistent",       "connectivity top-ranked at transect level (weight = 0.453), ICC = 0.306"
) %>% print()


# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()