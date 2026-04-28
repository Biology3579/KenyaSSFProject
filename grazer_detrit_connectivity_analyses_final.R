# ============================================================
#  DRIVERS OF GRAZER/DETRITIVORE BIOMASS
#  Chapter 1 — Functional Group Analysis: Grazers
#
#  Trophic group definition:
#    Grazers, detritivores, and grazer-detritivores treated
#    as a single functional guild. Roles overlap substantially
#    and combining reduces classification uncertainty.
#    Consistent with standard practice in Indo-Pacific
#    reef fish ecology.
#
#  Analytical framework mirrors total biomass and browsers
#  (Q1–Q3 + sensitivity).
#
#  Key differences from other functional groups:
#    Family: Gaussian log (lm) — Tweedie showed systematic
#    convergence failure across all model structures,
#    likely due to extremely high biomass variance in
#    this abundant guild. Gaussian log appropriate given
#    zero-free response and consistent with total biomass.
#    Zero proportion: 0% at site level — no offset needed.
#    Transect zeros: 2.1% — Tweedie at transect level.
#    Ecological context: heavily targeted guild —
#    rugosity and human pressure expected as primary
#    drivers, consistent with total biomass.
#
#  Sensitivity analysis:
#    (a) Alternative pressure metrics
#    (b) Transect-level replication (Tweedie GLMM)
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
cat("Sites:",           n_distinct(grazer_transects$site), "\n")
cat("Countries:",       n_distinct(grazer_transects$country), "\n")

# ── Site-level dataset ────────────────────────────────────────
grazer_model_data <- grazer_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass              = mean(transect_grazer_biomass,
                                     na.rm = TRUE),
    log_mean_biomass          = log(mean(transect_grazer_biomass,
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

cat("\nGrazer model data:", nrow(grazer_model_data), "sites\n")

# ── Data checks ───────────────────────────────────────────────
grazer_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                log_chla_sc, connectivity_sc, mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

cat("\nZeros in mean_biomass:",
    sum(grazer_model_data$mean_biomass == 0), "\n")
cat("-Inf in log_mean_biomass:",
    sum(is.infinite(grazer_model_data$log_mean_biomass)), "\n")
cat("Site-level zero proportion:",
    round(mean(grazer_model_data$mean_biomass == 0), 3), "\n")
cat("\nResponse summary:\n")
print(summary(grazer_model_data$log_mean_biomass))

# ── MPA classification check ──────────────────────────────────
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
#  Zero-free response — Gaussian log appropriate.
#  Tweedie attempted but showed systematic convergence
#  failure across all model structures — not used.
# ============================================================

grazer_lm_global <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  data = grazer_model_data
)

par(mfrow = c(2, 2))
plot(grazer_lm_global, main = "Gaussian log")
par(mfrow = c(1, 1))

# ── Family selection decision ─────────────────────────────────
#
# Gaussian log: SELECTED — lm() on log_mean_biomass.
#   Zero-free response (0 site-level zeros) — log
#   transformation valid without offset.
#
#   Diagnostics (global model):
#     Residuals vs Fitted: broadly flat with minor curve —
#       acceptable for n = 54. Sites 6, 48 elevated,
#       site 21 low but none extreme.
#     Q-Q: good through main body — site 6 deviates at
#       upper tail (standardised residual ~4), site 21
#       at lower tail (~-2.1). Both within acceptable
#       range for a zero-free continuous response.
#     Scale-Location: broadly flat — homoscedasticity
#       broadly met. Sites 6, 21, 48 show elevated
#       sqrt(standardised residuals) but no systematic
#       trend.
#     Residuals vs Leverage: site 6 has highest
#       standardised residual (~3.8) but leverage ~0.08
#       — well within Cook's distance threshold.
#       Site 21 approaches Cook's threshold at leverage
#       ~0.31 — monitor in Q2 and Q3 models.
#       No sites exceed 0.5 Cook's threshold.
#
#   Overall: Gaussian log diagnostics adequate.
#   Consistent with total biomass family selection.
#
# Tweedie: REJECTED — systematic convergence failure
#   across all model structures including baseline.
#   Non-positive-definite Hessian and false convergence
#   (8). Gaussian log is the appropriate substitute
#   for this zero-free high-variance guild.
#
#  Proceed: lm() on log_mean_biomass throughout all
#  grazer site-level analyses.
#  adj. R² used for variance explained.

# ============================================================
#  RANDOM EFFECT STRUCTURE
#  lm() vs lmer() with (1|country).
#  Baseline predictor set used — global model convergence
#  issues with Tweedie preclude global set test.
#  lmer() fitted with REML = FALSE for AICc comparison.
# ============================================================

grazer_re_null <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  data = grazer_model_data
)

grazer_re_country <- lmer(
  log_mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | country),
  data = grazer_model_data,
  REML = FALSE
)

cat("\n--- Grazer RE structure: country-level ---\n")
print(make_aicc_df(list(
  "No RE"         = grazer_re_null,
  "(1 | country)" = grazer_re_country
)))

# RE result: No RE: AICc = 130.50, weight = 0.809
# (1|country): ΔAICc = 2.89, weight = 0.191
# Country clustering not supported — consistent with
# all other analyses
# Remarkably consistent ΔAICc = 2.86–3.03 across all
# five analyses — confirms between-country differences
# are systematically absorbed by fixed predictors
# regardless of functional group or model family.
# All grazer models fitted as lm() without random
# effects throughout.
# ============================================================
#  Q1 — PRESSURE METRIC SELECTION
#
#  lm() on log_mean_biomass — consistent with family
#  selection. Three metrics compared.
#  Prior result showed complete metric uncertainty
#  (all ΔAICc < 0.30).
# ============================================================

# ── Without connectivity control ─────────────────────────────
g_q1_settgrav <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc,
  data = grazer_model_data
)

g_q1_mktgrav <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc,
  data = grazer_model_data
)

g_q1_settpop <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_settlement_pop_sc +
    log_chla_sc,
  data = grazer_model_data
)

cat("\n--- Q1: Grazer pressure metric (without connectivity) ---\n")
print(make_aicc_df(list(
  "Settlement gravity"   = g_q1_settgrav,
  "Market gravity"       = g_q1_mktgrav,
  "Settlement pop. 25km" = g_q1_settpop
)))

# ── With connectivity control ─────────────────────────────────
g_q1_settgrav_conn <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc,
  data = grazer_model_data
)

g_q1_mktgrav_conn <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc +
    connectivity_sc,
  data = grazer_model_data
)

g_q1_settpop_conn <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_settlement_pop_sc +
    log_chla_sc +
    connectivity_sc,
  data = grazer_model_data
)

cat("\n--- Q1: Grazer pressure metric (with connectivity) ---\n")
print(make_aicc_df(list(
  "Settlement gravity"   = g_q1_settgrav_conn,
  "Market gravity"       = g_q1_mktgrav_conn,
  "Settlement pop. 25km" = g_q1_settpop_conn
)))

cat("\n--- Q1: Grazer coefficient summary ---\n")
summary(g_q1_settgrav)

# ── Q1 results ────────────────────────────────────────────────
#
# Without connectivity control:
#   Market gravity:        AICc = 126.11, weight = 0.400
#   Settlement gravity:    ΔAICc = 0.53,  weight = 0.307
#   Settlement pop. 25km:  ΔAICc = 0.62,  weight = 0.293
#   All three within ΔAICc < 1 — uncertainty without
#   spatial control.
#
# With connectivity control:
#   Market gravity:        AICc = 123.65, weight = 0.520
#   Settlement gravity:    ΔAICc = 1.48,  weight = 0.248
#   Settlement pop. 25km:  ΔAICc = 1.61,  weight = 0.232
#   Market gravity has twice the weight of settlement
#   gravity — clearest metric differentiation of any
#   non-browser functional group.
#
# Market gravity selected as primary metric:
#   (1) Best or competitive across both comparisons
#   (2) Weight doubles settlement gravity once spatial
#       structure controlled (0.520 vs 0.248)
#   (3) Ecologically interpretable — grazers/detritivores
#       include commercially targeted species (surgeonfish,
#       rabbitfish) sold in urban markets. Market access
#       may capture exploitation pressure more accurately
#       for this partially commercial guild than
#       residential proximity.
#   (4) Only functional group where a metric other than
#       settlement gravity is clearly preferred — reflects
#       genuine group-specific differences in how SSF
#       pressure operates.
#
# Contrasts with total biomass (settlement gravity
# weight = 0.878), browsers (0.581), corallivores,
# and excavators (all complete uncertainty) — grazers
# show the most ecologically interpretable metric
# differentiation of any functional group.
#
# Settlement gravity and settlement population retained
# for sensitivity analysis.

cat("\n--- Q1: Grazer coefficient summary — market gravity ---\n")
summary(g_q1_mktgrav)

# ── Q1: Baseline model coefficients (market gravity) ──────────
#
# lm(log_mean_biomass ~ rugosity_sc + log_market_gravity_sc
#    + log_chla_sc)
# n = 54 sites, adj. R² = 0.072, F(3,50) = 2.364, p = 0.082
#
# Rugosity:         β = +0.262, p = 0.013 *
# Market gravity:   β = +0.082, p = 0.449 ns
# Chla:             β = -0.043, p = 0.697 ns
#
# Rugosity the only significant predictor — consistent
# with all other non-browser functional groups.
# Market gravity positive and non-significant — wrong
# direction for a fishing pressure effect, consistent
# with grazers not being structured by direct SSF
# pressure at baseline level despite market gravity
# being the preferred metric in Q1.
# Overall baseline explains minimal variance (adj.
# R² = 0.072) — consistent with grazers being weakly
# structured by local predictors.
# Slightly better fit than settlement gravity baseline
# (adj. R² = 0.063) — consistent with Q1 AICc result.

# ============================================================
#  Q2 — ARE CONNECTIVITY AND MPA BASELINE DRIVERS?
#
#  lm() on log_mean_biomass — consistent with family.
#  adj. R² used for variance explained.
#  Market gravity as primary pressure metric.
# ============================================================

g_null <- lm(log_mean_biomass ~ 1,
             data = grazer_model_data)

g_baseline <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc,
  data = grazer_model_data
)

g_baseline_conn <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc +
    connectivity_sc,
  data = grazer_model_data
)

g_baseline_mpa <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc +
    mpa_status,
  data = grazer_model_data
)

g_global_additive <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  data = grazer_model_data
)

g_q2_models <- list(
  "Null"                  = g_null,
  "Baseline"              = g_baseline,
  "Baseline + conn"       = g_baseline_conn,
  "Baseline + MPA"        = g_baseline_mpa,
  "Baseline + conn + MPA" = g_global_additive
)

cat("\n--- Q2: Grazer model comparison (AICc ranked) ---\n")
print(make_aicc_df(g_q2_models))

cat("\n--- Q2: Variance explained relative to baseline ---\n")
baseline_r2_g <- summary(g_baseline)$adj.r.squared

g_q2_models %>%
  imap_dfr(~ tibble(
    Model  = .y,
    Adj_R2 = round(summary(.x)$adj.r.squared, 3)
  )) %>%
  mutate(
    Delta_R2 = round(Adj_R2 - baseline_r2_g, 3),
    Delta_R2 = ifelse(Model %in% c("Null", "Baseline"),
                      NA, Delta_R2)
  ) %>%
  print()

cat("\n--- Q2: Baseline + connectivity coefficients ---\n")
summary(g_baseline_conn)


# ── Q2 results ────────────────────────────────────────────────
#
# AICc comparison:
#   Baseline + conn:       AICc = 123.65, weight = 0.598 (BEST)
#   Baseline:              ΔAICc = 2.45,  weight = 0.176
#   Null:                  ΔAICc = 2.60,  weight = 0.163
#   Baseline + conn + MPA: ΔAICc = 5.35,  weight = 0.041
#   Baseline + MPA:        ΔAICc = 6.60,  weight = 0.022
#
# Variance explained relative to baseline (adj. R² = 0.072):
#   Baseline + conn:       ΔR² = +0.064
#   Baseline + conn + MPA: ΔR² = +0.029
#   Baseline + MPA:        ΔR² = -0.021
#
# Connectivity clearly supported — Baseline + conn best-
# supported model (weight = 0.598). Combined weight of
# models including connectivity = 0.639. ΔR² = +0.064 —
# largest single predictor contribution for grazers.
#
# MPA not supported — Baseline + MPA worse than baseline
# and null (ΔR² = -0.021). Contrasts with browsers
# (MPA strongly supported) and corallivores (MPA
# significant in global model).
#
# Baseline + conn coefficients (best model):
#   Rugosity:         β = +0.268, p = 0.009 **
#   Market gravity:   β = +0.133, p = 0.217 ns
#   Chla:             β = -0.091, p = 0.402 ns
#   Connectivity:     β = -0.224, p = 0.034 *
#   adj. R² = 0.136, F(4,49) = 3.092, p = 0.024
#
# Connectivity significant and negative — well-connected
# sites have lower grazer biomass. Consistent with
# corallivores (β = -0.231*) — negative connectivity
# signal emerging as a pattern across functional groups
# not directly targeted by SSF.
#
# Rugosity strengthens when connectivity included
# (β = +0.262 → +0.268, p = 0.013 → 0.009) — habitat
# complexity signal more precisely estimated once
# spatial structure is controlled.
#
# Market gravity remains non-significant (p = 0.217)
# — consistent with Q1 baseline. Pressure does not
# independently structure grazer biomass once
# connectivity is accounted for.
#
# Best-supported model: Baseline + conn (weight = 0.598).
# Used as reference for Q3.

# ============================================================
#  SPATIAL AUTOCORRELATION CHECK
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

cat("\n--- Spatial autocorrelation: grazer best Q2 model ---\n")
print(moran.test(residuals(g_baseline_conn), listw5_g))

# ── Spatial autocorrelation result ────────────────────────────
# Moran's I = +0.187, p = 0.002 — significant spatial
# autocorrelation in grazer best Q2 model residuals.
#
# This is the strongest spatial signal of any functional
# group and the only one where p < 0.05 besides total
# biomass (I = +0.140, p = 0.015).
#
# Pattern across all analyses:
#   Total biomass: I = +0.140, p = 0.015 — significant
#   Browsers:      I = -0.074, p = 0.802 — no signal
#   Corallivores:  I = -0.051, p = 0.669 — no signal
#   Excavators:    I = -0.026, p = 0.537 — no signal
#   Grazers:       I = +0.187, p = 0.002 — significant
#
# Positive Moran's I indicates spatial clustering —
# neighbouring sites have more similar residuals than
# expected by chance. The signal is stronger than total
# biomass, suggesting grazer biomass has spatial
# structure not captured by rugosity, market gravity,
# chla, and connectivity.
#
# Spatial error modelling not pursued for the same
# reasons as total biomass:
# (1) Country-level RE tested and not supported
#     (ΔAICc = 2.89) — between-country clustering
#     absorbed by fixed predictors.
# (2) Discontinuous sampling design — k-NN weights
#     bridge across four isolated country clusters,
#     not reflecting true spatial covariance structure
#     (Kissling & Carl 2008, Dormann et al. 2007).
#
# Residual spatial autocorrelation acknowledged as a
# limitation for grazer analysis specifically. May
# inflate type I error rates for connectivity and
# rugosity coefficients more than for other groups.
# Qualitative conclusions remain valid but p-values
# should be interpreted with additional caution.

# ============================================================
#  Q3 — DO MANAGEMENT AND CONNECTIVITY MODIFY
#       GRAZER PRESSURE EFFECTS?
#
#  Reference model: Baseline + conn
#  lm() throughout. Market gravity as pressure metric.
# ============================================================

# Reference: update after Q2
g_q3_reference <- g_baseline_conn  

# ── H1: MPA effectiveness depends on fishing intensity ────────
# Connectivity retained as main effect — present in reference
g_int_mpa_press <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status * log_market_gravity_sc,
  data = grazer_model_data
)

# ── H2: Connectivity buffers exploitation effects ─────────────
# Connectivity × pressure interaction — tests whether the
# relationship between connectivity and biomass depends
# on pressure level, beyond the additive main effect
g_int_conn_press <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc * log_market_gravity_sc,
  data = grazer_model_data
)

# ── H3: MPA effectiveness depends on larval supply ────────────
# Connectivity retained as main effect in interaction term
g_int_mpa_conn <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_market_gravity_sc +
    mpa_status * connectivity_sc,
  data = grazer_model_data
)

g_q3_models <- list(
  "Reference (best Q2)"    = g_q3_reference,
  "H1: MPA × pressure"     = g_int_mpa_press,
  "H2: Conn × pressure"    = g_int_conn_press,
  "H3: MPA × connectivity" = g_int_mpa_conn
)

cat("\n--- Q3: Grazer interaction comparison ---\n")
print(make_aicc_df(g_q3_models))

cat("\n--- Q3: H2 Connectivity × pressure coefficients ---\n")
summary(g_int_conn_press)

# ── H2 exploratory checks ────────────────────────────────────
# Check 1 — VIF
cat("\n--- H2 grazer: VIF check ---\n")
check_collinearity(g_int_conn_press)

# Check 2 — raw data plot
ggplot(grazer_model_data,
       aes(x = log_market_gravity_sc,
           y = log_mean_biomass,
           colour = connectivity_sc)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_colour_viridis_c(name = "Connectivity (sc)") +
  labs(x = "log(Market gravity) (standardised)",
       y = "log(Grazer biomass)") +
  theme_bw(base_size = 12)

# Check 3 — interaction plot at low/medium/high connectivity
conn_quantiles <- quantile(grazer_model_data$connectivity_sc,
                           c(0.10, 0.50, 0.90))

pred_grid <- expand.grid(
  log_market_gravity_sc = seq(
    min(grazer_model_data$log_market_gravity_sc),
    max(grazer_model_data$log_market_gravity_sc),
    length.out = 100),
  connectivity_sc = conn_quantiles,
  rugosity_sc     = 0,
  log_chla_sc     = 0
)

pred_grid$fit <- predict(g_int_conn_press, newdata = pred_grid)
pred_grid$conn_label <- factor(
  round(pred_grid$connectivity_sc, 2),
  labels = c("Low connectivity (10th)",
             "Medium connectivity (50th)",
             "High connectivity (90th)"))

ggplot(pred_grid,
       aes(x = log_market_gravity_sc,
           y = fit,
           colour = conn_label,
           group  = conn_label)) +
  geom_line(linewidth = 1.1) +
  labs(x      = "log(Market gravity) (standardised)",
       y      = "Predicted log(Grazer biomass)",
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top")

# Check 4 — diagnostic plots
par(mfrow = c(2, 2))
plot(g_int_conn_press)
par(mfrow = c(1, 1))

# ── H2 exploratory checks ─────────────────────────────────────
#
# VIF: all < 2 — no collinearity concern.
#
# Diagnostics: broadly adequate. Residuals vs Fitted
#   flat — no systematic pattern. Q-Q good through
#   main body. Site 6 high standardised residual (~4)
#   but low leverage (~0.08) — not influential. Site 21
#   moderate leverage (~0.34) but within Cook's threshold.
#   Homoscedasticity broadly met.
#
# Interaction plot: clear fan pattern across full
#   observed pressure range — not extrapolation.
#   Low connectivity (10th): strong positive slope
#   High connectivity (90th): negative/flat slope
#   Lines cross near z = -0.8 within observed data.
#
# Raw data: pattern visible in observed data —
#   yellow/green points (high connectivity) flat or
#   negative across pressure gradient. Purple points
#   (low connectivity) show positive trend. Not a
#   model artefact.
#
# Biological interpretation: connectivity modifies
#   the market gravity-biomass relationship for grazers.
#   At isolated (low connectivity) sites, proximity to
#   markets associates positively with grazer biomass —
#   possibly reflecting productive coastal areas that
#   support both human settlement and grazer habitat.
#   At well-connected sites, market access associates
#   negatively — connectivity may amplify exploitation
#   pressure by linking sites into broader fishing
#   networks, making well-connected high-market sites
#   more vulnerable to depletion.
#
# Conclusion: H2 grazer interaction is ecologically
#   real and interpretable. The fan pattern spans the
#   full observed data range, VIF is clean, and the
#   biological mechanism is coherent. This is the only
#   functional group where a connectivity × pressure
#   interaction is supported — contrasts with browsers
#   (MPA × pressure) and corallivores/excavators
#   (no interaction).

# ── Q3 results ────────────────────────────────────────────────
#
# AICc comparison:
#   H2: Conn × pressure:     AICc = 121.25, weight = 0.764 (BEST)
#   Reference (best Q2):     ΔAICc = 2.41,  weight = 0.229
#   H3: MPA × connectivity:  ΔAICc = 10.53, weight = 0.004
#   H1: MPA × pressure:      ΔAICc = 11.13, weight = 0.003
#
# H2 (Conn × pressure) clearly supported (weight = 0.764,
# ΔAICc = 2.41 vs reference). H1 and H3 not supported.
#
# H2 coefficients:
#   Rugosity:                    β = +0.348, p = 0.001 **
#   Chla:                        β = -0.056, p = 0.596 ns
#   Connectivity:                β = -0.221, p = 0.031 *
#   Market gravity:              β = +0.232, p = 0.044 *
#   Conn × market gravity:       β = -0.289, p = 0.035 *
#   adj. R² = 0.197, F(5,48) = 3.602, p = 0.008
#
# The interaction is negative — the positive market
# gravity effect weakens and reverses as connectivity
# increases. At low-connectivity sites, higher market
# access associates positively with grazer biomass.
# At high-connectivity sites, this relationship
# reverses — well-connected sites near markets show
# lower grazer biomass, consistent with connectivity
# amplifying exploitation pressure by linking sites
# into broader regional fishing networks.
#
# Rugosity strengthens in H2 (β = +0.348**) compared
# to reference (β = +0.268*) — habitat complexity
# signal more precisely estimated once the interaction
# structure captures pressure × connectivity variance.
#
# Exploratory checks confirm genuine signal:
# VIF all < 2, fan pattern spans full observed range,
# pattern visible in raw data. See H2 checks above.
#
# Overall Q3 conclusion:
# Connectivity modifies the market gravity-biomass
# relationship for grazers — the only functional group
# where a connectivity × pressure interaction is
# supported. This contrasts with browsers (MPA ×
# pressure) and corallivores/excavators (no interaction).
# The result suggests grazer depletion risk is highest
# at well-connected sites with high market access —
# a management-relevant finding for reef fish conservation
# in the WIO.

# ============================================================
#  SENSITIVITY ANALYSIS
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
# Primary metric: market gravity.
# Sensitivity tests settlement gravity and settlement
# population to confirm Q2 conclusions not metric-dependent.

g_sens_settgrav <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc,
  data = grazer_model_data
)

g_sens_settpop <- lm(
  log_mean_biomass ~ rugosity_sc +
    log_settlement_pop_sc +
    log_chla_sc,
  data = grazer_model_data
)

cat("\n--- Sensitivity (a): grazer alternative metrics ---\n")
cat("Settlement gravity:\n")
print(summary(g_sens_settgrav)$coefficients)
cat("\nSettlement population 25km:\n")
print(summary(g_sens_settpop)$coefficients)

# ── Sensitivity (a) results ───────────────────────────────────
#
# Alternative pressure metrics — baseline structure retained,
# only pressure metric substituted.
#
# Settlement gravity (β = +0.037, p = 0.766):
#   Not significant, near-zero. Rugosity remains
#   significant (β = +0.249, p = 0.017). Chla
#   non-significant. Pattern identical to market
#   gravity baseline — rugosity sole driver.
#
# Settlement population 25km (β = +0.005, p = 0.960):
#   Not significant, essentially zero. Rugosity
#   remains significant (β = +0.246, p = 0.019).
#   Chla non-significant.
#
# Primary Q2 conclusions robust across all metrics:
# (1) Rugosity significant and positive throughout
#     (β = 0.246–0.262, p < 0.02 across all three)
# (2) Pressure metrics non-significant regardless
#     of metric — market gravity β = +0.082,
#     settlement gravity β = +0.037, settlement
#     pop. β = +0.005. All near-zero and positive —
#     no fishing pressure signal for grazers at
#     baseline level.
# (3) The market gravity preference from Q1 is
#     confirmed as metric-independent at baseline —
#     the market gravity signal only emerges when
#     connectivity is included (Q2) and interacted
#     (Q3), not in isolation.
#
# Connectivity signal from Q2 not directly testable
# here — sensitivity (a) uses baseline structure
# without connectivity. Key robustness confirmed:
# rugosity is stable, pressure non-significant
# across all three metrics.

# ── Transect family selection ─────────────────────────────────
# lmer() on log-transformed transect biomass — zero
# transects removed (n = 5, 2.1%) to avoid log(0)
# without offset. Removal negligible given n = 243.
# Gaussian log consistent with site-level family.
# (1|site) accounts for within-site correlation.

grazer_transect_data_nozero <- grazer_transect_data %>%
  filter(transect_grazer_biomass > 0) %>%
  mutate(log_transect_biomass = log(transect_grazer_biomass))

cat("Transects after zero removal:",
    nrow(grazer_transect_data_nozero), "\n")
cat("Sites retained:",
    n_distinct(grazer_transect_data_nozero$site), "\n")

g_trans_null <- lmer(
  log_transect_biomass ~ 1 + (1 | site),
  data = grazer_transect_data_nozero,
  REML = FALSE
)

g_trans_baseline <- lmer(
  log_transect_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc +
    (1 | site),
  data = grazer_transect_data_nozero,
  REML = FALSE
)

g_trans_conn <- lmer(
  log_transect_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc +
    connectivity_sc +
    (1 | site),
  data = grazer_transect_data_nozero,
  REML = FALSE
)

g_trans_mpa <- lmer(
  log_transect_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc +
    mpa_status +
    (1 | site),
  data = grazer_transect_data_nozero,
  REML = FALSE
)

g_trans_global <- lmer(
  log_transect_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | site),
  data = grazer_transect_data_nozero,
  REML = FALSE
)

# ── Diagnostics ───────────────────────────────────────────────
par(mfrow = c(1, 2))
plot(fitted(g_trans_global), residuals(g_trans_global),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2, col = "grey60")
lines(lowess(fitted(g_trans_global),
             residuals(g_trans_global)), col = "red")
qqnorm(residuals(g_trans_global),
       main = "Q-Q Residuals")
qqline(residuals(g_trans_global), col = "red")
par(mfrow = c(1, 1))

# ── Transect diagnostics ──────────────────────────────────────
# lmer() on log(transect_grazer_biomass), zeros removed
# (n = 5), n = 238 transects retained.
#
# Residuals vs Fitted: broadly flat — minor downward
#   trend at high fitted values, acceptable for n = 238.
#   No systematic pattern suggesting model misspecification.
# Q-Q: follows theoretical line well through main body.
#   Minor deviations at both tails — acceptable.
#   Upper tail (2 points ~2.5) and lower tail (~-2.5)
#   within expected range for this sample size.
#
# lmer() confirmed as adequate family for grazer
# transect-level sensitivity analysis.
# Proceed with AICc comparison and coefficient summaries.

cat("\n--- Sensitivity (b): grazer transect comparison (lmer) ---\n")
print(make_aicc_df(list(
  "Null"                  = g_trans_null,
  "Baseline"              = g_trans_baseline,
  "Baseline + conn"       = g_trans_conn,
  "Baseline + MPA"        = g_trans_mpa,
  "Baseline + conn + MPA" = g_trans_global
)))

cat("\n--- Sensitivity (b): baseline coefficients ---\n")
summary(g_trans_baseline)

cat("\n--- Sensitivity (b): best model coefficients ---\n")
summary(g_trans_conn)


# ── Sensitivity (b) results ───────────────────────────────────
#
# AICc comparison (transect level, lmer, n = 238):
#   Baseline + conn:       AICc = 627.24, weight = 0.500 (BEST)
#   Baseline:              ΔAICc = 1.02,  weight = 0.300
#   Null:                  ΔAICc = 3.29,  weight = 0.096
#   Baseline + conn + MPA: ΔAICc = 4.21,  weight = 0.061
#   Baseline + MPA:        ΔAICc = 4.87,  weight = 0.044
#
# Model ordering fully consistent with site-level Q2:
#   Site level:      Baseline + conn best (weight = 0.598)
#   Transect level:  Baseline + conn best (weight = 0.500)
#
# MPA not supported at either scale. Primary Q2
# conclusion robust to aggregation level.
#
# Best model coefficients (Baseline + conn):
#   Rugosity:         β = +0.248, t = 2.908 — positive
#   Market gravity:   β = +0.057, t = 0.614 — ns
#   Chla:             β = -0.132, t = -1.417 — ns
#   Connectivity:     β = -0.159, t = -1.802 — marginal
#   Site variance:    0.226 (SD = 0.475)
#   Residual:         0.612 (SD = 0.783)
#
# Connectivity direction consistent with site level
# (β = -0.224* at site, β = -0.159 marginal at transect).
# Signal weakens slightly at transect level — within-site
# variance in grazer biomass reduces precision — but
# direction fully preserved. Site variance reduces from
# 0.248 (baseline) to 0.226 (baseline + conn) —
# connectivity absorbs modest between-site variance.
#
# Rugosity consistent across scales (β = +0.241 baseline,
# +0.248 with connectivity at transect vs +0.268** at
# site level) — direction and magnitude stable.
#
# Primary Q2 and Q3 conclusions robust:
# connectivity signal replicates directionally at
# transect level despite slight attenuation.


# ============================================================
#  MARGINAL EFFECT PLOTS
#  Significant predictors in best Q2 model:
#  rugosity (positive) and connectivity (negative).
#  lm() — predictions on log scale, back-transformed
#  for display.
# ============================================================

# ── Rugosity effect ───────────────────────────────────────────
rug_grid_g <- data.frame(
  rugosity_sc           = seq(
    min(grazer_model_data$rugosity_sc),
    max(grazer_model_data$rugosity_sc),
    length.out = 200),
  log_market_gravity_sc = 0,
  log_chla_sc           = 0,
  connectivity_sc       = 0
)

rug_pred_g        <- predict(g_baseline_conn,
                             newdata  = rug_grid_g,
                             interval = "confidence")
rug_grid_g$fit    <- rug_pred_g[, "fit"]
rug_grid_g$lwr    <- rug_pred_g[, "lwr"]
rug_grid_g$upr    <- rug_pred_g[, "upr"]

p_g_rugosity <- ggplot(rug_grid_g,
                       aes(x = rugosity_sc, y = exp(fit))) +
  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)),
              fill = "#2c7bb6", alpha = 0.15) +
  geom_line(colour = "#2c7bb6", linewidth = 1.1) +
  geom_point(data = grazer_model_data,
             aes(x = rugosity_sc, y = mean_biomass),
             colour = "grey40", size = 1.5,
             alpha = 0.5, inherit.aes = FALSE) +
  labs(x = "Rugosity (standardised)",
       y = "Grazer biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ── Connectivity effect ───────────────────────────────────────
conn_grid_g <- data.frame(
  connectivity_sc       = seq(
    min(grazer_model_data$connectivity_sc),
    max(grazer_model_data$connectivity_sc),
    length.out = 200),
  rugosity_sc           = 0,
  log_market_gravity_sc = 0,
  log_chla_sc           = 0
)

conn_pred_g        <- predict(g_baseline_conn,
                              newdata  = conn_grid_g,
                              interval = "confidence")
conn_grid_g$fit    <- conn_pred_g[, "fit"]
conn_grid_g$lwr    <- conn_pred_g[, "lwr"]
conn_grid_g$upr    <- conn_pred_g[, "upr"]

p_g_conn <- ggplot(conn_grid_g,
                   aes(x = connectivity_sc, y = exp(fit))) +
  geom_ribbon(aes(ymin = exp(lwr), ymax = exp(upr)),
              fill = "#2c7bb6", alpha = 0.15) +
  geom_line(colour = "#2c7bb6", linewidth = 1.1) +
  geom_point(data = grazer_model_data,
             aes(x = connectivity_sc, y = mean_biomass),
             colour = "grey40", size = 1.5,
             alpha = 0.5, inherit.aes = FALSE) +
  labs(x = "Connectivity (standardised)",
       y = "Grazer biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

gridExtra::grid.arrange(p_g_rugosity, p_g_conn, ncol = 2)


# ============================================================
#  RESULTS SUMMARY
# ============================================================

grazer_results <- tribble(
  ~Question,  ~Best_model,          ~Key_finding,
  "Q1",       "Market gravity",     "weight = 0.520 with connectivity — only group preferring non-settlement metric",
  "Q2 conn",  "Baseline + conn",    "β = -0.224, p = 0.034; negative — consistent with corallivores",
  "Q2 MPA",   "Baseline",           "ΔAICc = 6.60 — not supported; ΔR² = -0.021",
  "Q3 H1",    "Reference",          "ΔAICc = 11.13 — not supported",
  "Q3 H2",    "H2 (conn×pressure)", "weight = 0.764; β = -0.289, p = 0.035 — genuine interaction confirmed",
  "Q3 H3",    "Reference",          "ΔAICc = 10.53 — not supported"
)

cat("\n--- Grazer results summary ---\n")
print(grazer_results)


# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()