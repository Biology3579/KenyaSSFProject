# ============================================================
#  TOTAL BIOMASS ANALYSIS
#  Created by: Candela Ferrer Diez
#  Date: 24/01/2026
# ============================================================
#
#  Scientific questions:
#  Q1 — Pressure metric validity
#       Does settlement gravity outperform market gravity and
#       settlement population as a proxy for small-scale
#       fisheries pressure on total reef fish biomass?
#       Tested via AICc comparison — four models identical
#       in structure, differing only in pressure metric.
#
#  Q2 — Larval connectivity
#       Does larval connectivity explain additional variation
#       in total biomass beyond the human pressure baseline,
#       and does it modify the relationship between human
#       pressure and total biomass?
#       Tested as main effect and connectivity × pressure
#       interaction (a priori hypothesis).
#
#  Q3 — Formal protection
#       Does MPA status explain additional variation in
#       total biomass beyond the best-supported pressure
#       and connectivity model?
#       Separate question — does not update the best model.
#
#  Rationale for sequence:
#       Pressure first as primary driver of interest,
#       connectivity second to extend Warmuth et al. (2024),
#       MPA last as governance response downstream of both.
#
#  Baseline model (fixed a priori, never tested):
#       log(biomass) ~ rugosity_sc + log_chla_sc
#
#  Model family:
#       Gaussian lm on log-transformed biomass throughout
#       (see family selection section below).
#       No zeros in total biomass at site level.
#
#  Sensitivity analysis:
#       Transect-level replication — confirms site-level
#       findings not an artefact of aggregation.
# ============================================================

# ── Source predictors ───────────────────────────
# Loads raw data, and prepared predictors.
# Also loads key functions

source(here::here("predictor_preparation.R"))

# ============================================================
#  RESPONSE VARIABLE — DATA AGGREGATION AND CHECKS
# ============================================================

# ── Site-level biomass summary ────────────────────────────────
site_data <- total_transects %>%
  group_by(site) %>%
  summarise(mean_biomass = mean(transect_total_biomass,
                                na.rm = TRUE),
            n_transects  = n(),
            .groups = "drop") %>%
  mutate(site = as.factor(site))

cat("Sites:", nrow(site_data), "\n")
cat("Proportion of zeros:",
    round(mean(site_data$mean_biomass == 0), 3), "\n")

# ── Raw distribution ──────────────────────────────────────────
( p_raw <- ggplot(site_data, aes(x = mean_biomass)) +
    geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
    labs(x = "Mean total biomass per site (g)", y = "Frequency",
         title = "Raw site-level biomass") +
    theme_bw() )
# Strong right skew — log transformation required

# ── Normality checks ──────────────────────────────────────────
cat("\n--- Shapiro-Wilk: raw biomass ---\n")
shapiro.test(site_data$mean_biomass)
# Expected: W << 1, p < 0.001 — strongly non-normal

# ── Box-Cox ───────────────────────────────────────────────────
MASS::boxcox(
  lm(mean_biomass ~ 1, data = site_data),
  lambda = seq(-2, 2, 0.1)
)
# Lambda CI includes 0 — log transformation appropriate

# ── Log transformation ────────────────────────────────────────
site_data <- site_data %>%
  mutate(log_mean_biomass = log(mean_biomass))

( p_log <- ggplot(site_data, aes(x = log_mean_biomass)) +
    geom_histogram(bins = 25, fill = "#1a9641", colour = "white") +
    labs(x = "log(Mean biomass)", y = "Frequency",
         title = "Log-transformed site-level biomass") +
    theme_bw() )

qqnorm(site_data$log_mean_biomass,
       main = "Q-Q plot: log(mean biomass per site)")
qqline(site_data$log_mean_biomass, col = "red")

cat("\n--- Shapiro-Wilk: log(biomass) ---\n")
shapiro.test(site_data$log_mean_biomass)
# No strong evidence of departure from normality —
# log transformation sufficient

# ── Biomass by site ───────────────────────────────────────────
ggplot(site_data,
       aes(x = reorder(site, mean_biomass), y = mean_biomass)) +
  geom_col(fill = "#2c7bb6", alpha = 0.7) +
  coord_flip() +
  labs(x = NULL, y = "Mean total biomass (g)",
       title = "Mean biomass by site (raw)") +
  theme_bw(base_size = 9)


# ============================================================
#  ANALYSIS DATASETS
# ============================================================

# ── Transect-level dataset ────────────────────────────────────
transect_model_data <- total_transects %>%
  left_join(final_predictors, by = "site") %>%
  mutate(log_transect_biomass = log(transect_total_biomass))

cat("\nTransect data:", nrow(transect_model_data), "rows,",
    n_distinct(transect_model_data$site), "sites\n")
cat("Biomass zeros:",
    sum(transect_model_data$transect_total_biomass == 0), "\n")

# ── Site-level dataset ────────────────────────────────────────
total_model_data <- transect_model_data %>%
  group_by(site) %>%
  summarise(
    log_mean_biomass       = log(mean(transect_total_biomass,
                                      na.rm = TRUE)),
    mean_biomass           = mean(transect_total_biomass,
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
    ecoregion = as.factor(ecoregion)
  )

cat("\nSite data:", nrow(total_model_data), "sites,",
    n_distinct(total_model_data$ecoregion), "ecoregions\n")

# ── Data checks ───────────────────────────────────────────────
total_model_data %>%
  dplyr::select(site, ecoregion, rugosity_sc, log_chla_sc,
                log_settlement_grav_sc, connectivity_sc,
                mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  { if (nrow(.) > 0) { warning("NAs in predictors:"); print(.) }
    else cat("NA check passed.\n") }

cat("\nZeros in mean_biomass:",
    sum(total_model_data$mean_biomass == 0), "\n")
cat("-Inf in log_mean_biomass:",
    sum(is.infinite(total_model_data$log_mean_biomass)), "\n")
cat("\nMPA status counts:\n")
print(table(total_model_data$mpa_status))
cat("\nResponse variable summary:\n")
print(summary(total_model_data$log_mean_biomass))


# ============================================================
#  MODEL FAMILY SELECTION
#
#  Assessed on baseline model — if family holds for baseline
#  it holds for all models built on it.
# ============================================================

lm_gaussian_raw <- lm(mean_biomass ~ rugosity_sc + log_chla_sc,
                      data = total_model_data)
lm_gaussian_log <- lm(log_mean_biomass ~ rugosity_sc + log_chla_sc,
                      data = total_model_data)
glm_gamma       <- glm(mean_biomass ~ rugosity_sc + log_chla_sc,
                       family = Gamma(link = "log"),
                       data   = total_model_data)

par(mfrow = c(2, 2))
plot(lm_gaussian_raw, main = "Gaussian raw")
plot(lm_gaussian_log, main = "Gaussian log")
plot(glm_gamma,       main = "Gamma log-link")
par(mfrow = c(1, 1))

# Family selection:
#   Gaussian (raw):   rejected — severe heteroscedasticity,
#                     U-shaped residuals vs fitted, heavy
#                     upper Q-Q tail.
#   Gamma (log link): rejected — systematic Q-Q deviation
#                     at upper tail despite flat residuals.
#   Gaussian (log):   SELECTED — broadly flat residuals,
#                     Q-Q closely follows theoretical line,
#                     minor upper tail deviation only.
#                     No sites exceed Cook's distance threshold.
#
#   Proceed: lm() on log_mean_biomass throughout.


# ============================================================
#  RANDOM EFFECT STRUCTURE
#
#  Tested on baseline model — adding unvalidated predictors
#  would presuppose Q1-Q3 outcomes and absorb between-
#  ecoregion variance, biasing the test.
#  ML (REML = FALSE) for AICc comparison.
# ============================================================

re_null <- lm(log_mean_biomass ~ rugosity_sc + log_chla_sc,
              data = total_model_data)

re_ecoregion <- lmer(log_mean_biomass ~ rugosity_sc + log_chla_sc +
                       (1 | ecoregion),
                     data = total_model_data,
                     REML = FALSE)

cat("\n--- Random effect structure: ecoregion ---\n")
print(make_aicc_df(list(
  "No RE"           = re_null,
  "(1 | ecoregion)" = re_ecoregion
)))

# Ecoregion RE not supported (DAICc = 2.25, weight = 0.245 vs 0.75 for the no RE).

# Between-ecoregion variation acknowledged as a limitation.
# All models fitted as lm() throughout.


# ============================================================
#  SST BASELINE CHECK (SENSITIVITY)
#
#  SST evaluated as additional baseline covariate.
#  Excluded if: (1) strong collinearity with connectivity
#  (confirmed in 00_data_preparation.R, r = 0.70) AND
#  (2) negligible improvement in model fit (DAICc < 2).
#  Both conditions must hold to justify exclusion.
# ============================================================

# Join SST to model data
sst_model <- sst_sites %>%
  mutate(site = as.character(site)) %>%
  dplyr::select(site, sst_sc)

total_model_data <- total_model_data %>%
  left_join(sst_model, by = "site")

m_baseline     <- lm(log_mean_biomass ~ rugosity_sc + log_chla_sc,
                     data = total_model_data)
m_baseline_sst <- lm(log_mean_biomass ~ rugosity_sc + log_chla_sc +
                       sst_sc,
                     data = total_model_data)

cat("\n--- SST baseline check ---\n")
print(make_aicc_df(list(
  "Baseline"       = m_baseline,
  "Baseline + SST" = m_baseline_sst
)))

# DAICc = 1.80 — SST not supported.
# Combined with r = 0.70 with connectivity
# (predictor_preparation.R), SST excluded from all subsequent analyses.


# ============================================================
#  Q1 — HUMAN PRESSURE
#
#  Scientific question:
#  Does human pressure explain variation in reef fish biomass
#  beyond local ecological context, and which spatial metric
#  best captures SSF exploitation intensity in the WIO?
#
#  A priori prediction:
#  Settlement gravity will outperform market gravity because
#  residential proximity better captures the spatial footprint
#  of subsistence-oriented exploitation than commercial market
#  access in this SSF-dominated system (Cinner et al. 2016,
#  Samoilys et al. 2019).
#
#  Approach: AICc comparison of baseline vs baseline + each
#  pressure metric. Best metric = highest AICc weight AND
#  outperforms baseline.
# ============================================================

q1_settgrav <- lm(log_mean_biomass ~ rugosity_sc + log_chla_sc +
                    log_settlement_grav_sc,
                  data = total_model_data)

q1_mktgrav  <- lm(log_mean_biomass ~ rugosity_sc + log_chla_sc +
                    log_market_gravity_sc,
                  data = total_model_data)

q1_settpop  <- lm(log_mean_biomass ~ rugosity_sc + log_chla_sc +
                    log_settlement_pop_sc,
                  data = total_model_data)

cat("\n--- Q1: Pressure metric comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                      = m_baseline,
  "Baseline + settlement gravity" = q1_settgrav,
  "Baseline + market gravity"     = q1_mktgrav,
  "Baseline + settlement pop."    = q1_settpop
)))

# Results:
#   Settlement gravity: AICc = 101.36, weight = 0.826
#   Baseline:           DAICc = 4.39,  weight = 0.092
#   Market gravity:     DAICc = 5.86,  weight = 0.044
#   Settlement pop.:    DAICc = 6.17,  weight = 0.038
#
#   Settlement gravity decisively selected (weight = 0.826).
#   Outperforms ecological baseline alone (DAICc = 4.39),
#   confirming human pressure explains meaningful variance
#   beyond habitat and productivity context.
#   Market gravity and settlement population not competitive
#   (combined weight = 0.082).

# ── Best Q1 model ─────────────────────────────────────────────
m_best_q1 <- q1_settgrav


# ============================================================
#  Q2 — LARVAL CONNECTIVITY
#
#  Scientific question:
#  Does larval connectivity explain additional variation in
#  reef fish biomass beyond the human pressure baseline, and
#  does it modify the relationship between fishing pressure
#  and biomass?
#
#  Directly extends Warmuth et al. (2024) — tests whether
#  the connectivity signal persists once pressure is controlled.
#
#  Approach: three-model AICc comparison.
#  (1) Best Q1 — pressure only
#  (2) Best Q1 + connectivity main effect
#  (3) Best Q1 + connectivity + connectivity × pressure
#
#  The interaction (connectivity × pressure) is evaluated
#  alongside the main effect because connectivity may modify
#  the pressure-biomass relationship rather than act
#  independently
#
#  Best model identified here is carried forward to Q3.
# ============================================================

m_q2_conn <- lm(log_mean_biomass ~ rugosity_sc + log_chla_sc +
                  log_settlement_grav_sc +
                  connectivity_sc,
                data = total_model_data)

m_q2_conn_int <- lm(log_mean_biomass ~ rugosity_sc + log_chla_sc +
                      log_settlement_grav_sc +
                      connectivity_sc +
                      log_settlement_grav_sc:connectivity_sc,
                    data = total_model_data)

cat("\n--- Q2: Connectivity main effect and interaction ---\n")
print(make_aicc_df(list(
  "Best Q1"                      = m_best_q1,
  "Best Q1 + conn"               = m_q2_conn,
  "Best Q1 + conn + interaction" = m_q2_conn_int
)))

summary(m_q2_conn_int)
confint(m_q2_conn_int)

# Results:
#   Best Q1:              AICc = 101.36, weight = 0.529
#   Best Q1 + conn + int: DAICc =  1.06, weight = 0.312
#   Best Q1 + conn:       DAICc =  2.41, weight = 0.158
#
#   Connectivity as a main effect not supported (DAICc = 2.41,
#   weight = 0.158).
#   Interaction model carries moderate weight (0.312,
#   DAICc = 1.06), but the interaction term is marginal
#   (b = -0.225, p = 0.060, 95% CI [-0.460, 0.010]) with
#   CIs overlapping zero. No directional prediction was made
#   for this interaction — given marginal statistical support
#   and wide uncertainty, it is not interpreted as ecologically
#   meaningful. Best Q1 model retained as m_best_q2.

# ── Best Q2 model ─────────────────────────────────────────────
m_best_q2 <- m_best_q1

# ── Best model summary — reported in results ─────────────────
# Full coefficient reporting for the best model identified
# across Q1-Q2. This is the primary reported model.
# Q3 is a separate question and does not update this model.

# ── Model diagnostics — best Q2 model ────────────────────────
par(mfrow = c(2, 2))
plot(m_best_q2, main = "Best Q2 model diagnostics")
par(mfrow = c(1, 1))

cat("\n--- Shapiro-Wilk: best Q2 model residuals ---\n")
shapiro.test(residuals(m_best_q2))

# Diagnostics — best Q2 model (rugosity + chla + settlement gravity):
#   Residuals vs Fitted: broadly flat with minor downward
#     curve — acceptable, no systematic non-linearity.
#     Sites 6, 25, 27 flagged as potential outliers.
#   Q-Q: excellent — points follow theoretical line closely
#     across full range. Site 27 minor lower tail deviation,
#     site 6 minor upper tail. No cause for concern.
#   Scale-Location: slight upward trend — minor heteroscedasticity
#     at high fitted values. Sites 6, 27, 35 flagged.
#     Acceptable for n = 54.
#   Residuals vs Leverage: no sites exceed Cook's distance
#     threshold. Sites 6, 35, 48 have moderate leverage
#     but standardised residuals within acceptable range.
#   Shapiro-Wilk: W = 0.984, p = 0.705 — no evidence of
#     departure from normality.
#   Overall: Gaussian lm assumptions adequately met.


cat("\n--- Best model: full summary ---\n")
summary(m_best_q2)

cat("\n--- Best model: 95% CIs ---\n")
print(confint(m_best_q2))

cat("\n--- Best model: R² ---\n")
cat(sprintf("R²     = %.3f\n", summary(m_best_q2)$r.squared))
cat(sprintf("Adj R² = %.3f\n", summary(m_best_q2)$adj.r.squared))

# ── Effect sizes: fold differences across observed range ──────
b_sg  <- coef(m_best_q2)["log_settlement_grav_sc"]
b_rug <- coef(m_best_q2)["rugosity_sc"]

sg_span  <- diff(range(total_model_data$log_settlement_grav_sc,
                       na.rm = TRUE))
rug_span <- diff(range(total_model_data$rugosity_sc,
                       na.rm = TRUE))

cat(sprintf("\nSettlement gravity — span: %.3f SD units\n", sg_span))
cat(sprintf("Fold difference (low vs high pressure): %.2fx\n",
            exp(abs(b_sg * sg_span))))

cat(sprintf("\nRugosity — span: %.3f SD units\n", rug_span))
cat(sprintf("Fold difference (low vs high rugosity): %.2fx\n",
            exp(abs(b_rug * rug_span))))

# Results:
#   Rugosity:           b = +0.216, t(50) =  2.707, p = 0.009 **
#     95% CI [0.056, 0.376]
#     Fold difference: 3.04x across observed range
#     (span = 5.151 SD units)
#   Settlement gravity: b = -0.251, t(50) = -2.596, p = 0.012 *
#     95% CI [-0.445, -0.057]
#     Fold difference: 2.73x across observed range
#     (span = 4.001 SD units)
#   Chla:               b = -0.125, t(50) = -1.275, p = 0.208 ns
#     95% CI [-0.322, 0.072]
#   Model: F(3,50) = 5.253, p = 0.003,
#     R² = 0.240, adj. R² = 0.194

# ============================================================
#  Q3 — FORMAL PROTECTION
#
#  Scientific question:
#  Does MPA status explain additional variation in reef fish
#  biomass beyond the best-supported model from Q1-Q2?
#
#  Q3 is a separate question — it does not update the best
#  model. The best Q2 model is used as the reference, and
#  MPA is evaluated as an additional candidate predictor.
#
#  MPA tested last because MPAs are non-randomly placed —
#  medium-protection sites face lower fishing pressure and
#  higher connectivity than low and unprotected sites
#  (see 00_mpa_placement.R). Apparent protection effects
#  are confounded by pressure and connectivity unless
#  both are controlled first.
#
#  Criterion: AICc weight for model support.
#  p-values for coefficient direction if supported.
# ============================================================

m_q3_mpa <- lm(log_mean_biomass ~ rugosity_sc + log_chla_sc +
                 log_settlement_grav_sc +
                 mpa_status,
               data = total_model_data)

cat("\n--- Q3: MPA model comparison ---\n")
print(make_aicc_df(list(
  "Best Q2"       = m_best_q2,
  "Best Q2 + MPA" = m_q3_mpa
)))

# ── Model diagnostics — Q3 model ─────────────────────────────
par(mfrow = c(2, 2))
plot(m_q3_mpa, main = "Q3 MPA model diagnostics")
par(mfrow = c(1, 1))

cat("\n--- Shapiro-Wilk: Q3 model residuals ---\n")
shapiro.test(residuals(m_q3_mpa))

# Diagnostics — Q3 model (rugosity + chla + settlement gravity + MPA):
#   Residuals vs Fitted: broadly flat — similar to Q2 model.
#     Sites 6, 25, 27 flagged consistently.
#   Q-Q: excellent — closely follows theoretical line.
#     Site 27 minor lower tail, site 6 minor upper tail —
#     consistent with Q2 model, not introduced by MPA addition.
#   Scale-Location: slight upward trend consistent with Q2 —
#     MPA addition does not worsen heteroscedasticity.
#   Residuals vs Leverage: site 35 leverage increases slightly
#     with MPA addition — MPA category may be influential
#     for this site. Sites 6, 19, 35 flagged. No sites
#     exceed Cook's distance threshold.
#   Shapiro-Wilk: W = 0.989, p = 0.893 — excellent,
#     slight improvement over Q2 model.
#   Overall: Gaussian assumptions well met. Residual structure
#     consistent across Q2 and Q3 models — MPA addition does
#     not introduce new violations.

cat("\n--- Q3: MPA coefficients ---\n")
summary(m_q3_mpa)

cat("\n--- Q3: MPA model 95% CIs ---\n")
print(confint(m_q3_mpa))

# Results:
#   Best Q2:       AICc = 101.36, weight = 0.819
#   Best Q2 + MPA: DAICc =  3.01, weight = 0.182
#   MPA status not supported beyond pressure baseline.
#
#   Low MPA:    b = -0.342, t(48) = -1.390, p = 0.171 ns
#     95% CI [-0.838, 0.153]
#   Medium MPA: b = -0.084, t(48) = -0.430, p = 0.669 ns
#     95% CI [-0.478, 0.310]
#   Neither protection level significantly predicts biomass
#   once habitat and pressure are controlled.
#   Negative coefficients reflect preferential MPA placement
#   in lower-pressure, higher-connectivity areas (see
#   00_mpa_placement.R) — not a genuine negative effect
#   of protection.
#
#   Rugosity and pressure coefficients stable throughout:
#     Rugosity:           b = +0.211, p = 0.011 *
#     Settlement gravity: b = -0.239, p = 0.036 *

# ── Identify flagged sites ────────────────────────────────────
# Row indices from diagnostic plots: 6, 25, 27, 35
flagged_rows <- c(6, 25, 27, 35)

total_model_data %>%
  slice(flagged_rows) %>%
  dplyr::select(site, mpa_status, ecoregion,
                mean_biomass, log_mean_biomass,
                rugosity_sc, log_settlement_grav_sc,
                connectivity_sc) %>%
  print()

# Flagged sites (row indices 6, 25, 27, 35):
#   Row 6  — ankao_s (none, Madagascar):
#     mean_biomass = 76,467g — highest biomass in the dataset.
#     Unprotected remote site in Madagascar with very high
#     rugosity likely driving extreme biomass. Genuine outlier
#     in the upper tail — real ecological signal, not an error.
#
#   Row 25 — metundo_ne (medium, Tanzania S-Mozambique):
#     mean_biomass = 7,011g — medium MPA site, moderate biomass.
#     Flagged for leverage in Q3 model — MPA category
#     influential given site characteristics.
#
#   Row 27 — metundo_nw (medium, Tanzania S-Mozambique):
#     mean_biomass = 5,239g — sister site to metundo_ne,
#     similar profile. Pair likely influential together
#     in Q3 model given shared MPA status and proximity.
#
#   Row 35 — mnemb (none, Kenya-Tanzania north):
#     mean_biomass = 26,207g — high biomass unprotected site.
#     Only 2 sites in Kenya-Tanzania north ecoregion —
#     this site has high leverage partly due to ecoregion
#     isolation. Consistent with acknowledged between-ecoregion
#     limitation.
#
# None of these sites warrant exclusion — all represent
# genuine ecological variation. ankao_s warrants a note
# as the highest-biomass site; its influence on the
# rugosity coefficient should be monitored but removal
# would not be scientifically justified.

# ============================================================
#  SPATIAL AUTOCORRELATION CHECK
#
#  Residuals from best model tested for spatial structure.
#  Reported as a diagnostic — not part of model selection.
# ============================================================

site_coords <- location_data %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(lon = first(longitude),
            lat = first(latitude),
            .groups = "drop")

total_model_data_coords <- total_model_data %>%
  left_join(site_coords, by = "site")

coords_mat <- cbind(total_model_data_coords$lon,
                    total_model_data_coords$lat)
listw5 <- nb2listw(knn2nb(knearneigh(coords_mat, k = 5)),
                   style = "W")

# Warning: knearneigh — identical points found, kd_tree not available.
# neighbour object has 4 sub-graphs.
# This is expected given the discontinuous sampling design —
# sites from four countries form geographically isolated
# clusters. k-NN weights bridge across these clusters,
# which is a known limitation already acknowledged in the
# spatial autocorrelation comment above. Warning does not
# invalidate the Moran's I result but reinforces why spatial
# error modelling was not pursued.

cat("\n--- Spatial autocorrelation: best model residuals ---\n")
print(moran.test(residuals(m_best_q2), listw5))

# Moran's I = 0.140, p = 0.015 — weak but significant.
#
# Spatial error modelling not pursued:
# (1) Ecoregion RE not supported (DAICc = 2.25) and only
#     4 groups with severely uneven sizes (n = 2, 8, 9, 35).
# (2) Discontinuous sampling design — k-NN weights bridge
#     across geographically isolated country clusters,
#     producing a weights matrix that does not reflect true
#     within-region spatial covariance structure
#     (Kissling & Carl 2008; Dormann et al. 2007).
#
# Weak residual autocorrelation acknowledged as a limitation.
# May slightly inflate type I error for pressure and rugosity.
# Substantive conclusions unlikely to change given effect
# sizes and consistency across model specifications.


# ============================================================
#  SENSITIVITY ANALYSIS — TRANSECT-LEVEL REPLICATION
#
#  Purpose: confirm site-level Q1-Q3 conclusions are not an
#  artefact of spatial aggregation. Replicates Q1-Q3 model
#  sequence at transect level using lmer with a site-level
#  random intercept.
#
#  Models fitted with REML = TRUE for coefficient inference;
#  refitted with REML = FALSE for AICc comparison.
# ============================================================

# ── Fit transect models (REML = TRUE) ────────────────────────
sens_t_null <- lmer(log_transect_biomass ~ 1 + (1 | site),
                    data = transect_model_data, REML = TRUE)

sens_t_baseline <- lmer(log_transect_biomass ~ rugosity_sc +
                          log_chla_sc + (1 | site),
                        data = transect_model_data, REML = TRUE)

sens_t_pressure <- lmer(log_transect_biomass ~ rugosity_sc +
                          log_chla_sc +
                          log_settlement_grav_sc + (1 | site),
                        data = transect_model_data, REML = TRUE)

sens_t_conn <- lmer(log_transect_biomass ~ rugosity_sc +
                      log_chla_sc +
                      log_settlement_grav_sc +
                      connectivity_sc + (1 | site),
                    data = transect_model_data, REML = TRUE)

sens_t_mpa <- lmer(log_transect_biomass ~ rugosity_sc +
                     log_chla_sc +
                     log_settlement_grav_sc +
                     mpa_status + (1 | site),
                   data = transect_model_data, REML = TRUE)

# ── Diagnostics on baseline model ────────────────────────────
par(mfrow = c(1, 2))
plot(fitted(sens_t_baseline), residuals(sens_t_baseline),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Transect baseline: Residuals vs Fitted")
abline(h = 0, lty = 2, col = "grey60")
lines(lowess(fitted(sens_t_baseline),
             residuals(sens_t_baseline)), col = "red")
qqnorm(residuals(sens_t_baseline),
       main = "Transect baseline: Q-Q Residuals")
qqline(residuals(sens_t_baseline), col = "red")
par(mfrow = c(1, 1))

# ── AICc comparison — ML refits ──────────────────────────────
sens_t_null_ml     <- update(sens_t_null,     REML = FALSE)
sens_t_baseline_ml <- update(sens_t_baseline, REML = FALSE)
sens_t_pressure_ml <- update(sens_t_pressure, REML = FALSE)
sens_t_conn_ml     <- update(sens_t_conn,     REML = FALSE)
sens_t_mpa_ml      <- update(sens_t_mpa,      REML = FALSE)

cat("\n--- Sensitivity: transect-level model comparison ---\n")
print(make_aicc_df(list(
  "Null"                = sens_t_null_ml,
  "Baseline"            = sens_t_baseline_ml,
  "Baseline + pressure" = sens_t_pressure_ml,
  "Best + conn"         = sens_t_conn_ml,
  "Best + MPA"          = sens_t_mpa_ml
)))

# Results:
#   Baseline + pressure: AICc = 634.93, weight = 0.626
#   Best + conn:         DAICc = 2.10,  weight = 0.219
#   Best + MPA:          DAICc = 3.08,  weight = 0.134
#   Baseline:            DAICc = 6.84,  weight = 0.021
#   Null:                DAICc = 13.05, weight = 0.001
#
#   Model ordering identical to site-level Q1-Q3.
#   Connectivity and MPA not supported — consistent
#   with site-level results.

# ── Coefficient summary — REML pressure model ────────────────
cat("\n--- Sensitivity: pressure model coefficients (REML) ---\n")
summary(sens_t_pressure)
print(confint(sens_t_pressure))

# Results:
#   Rugosity:           b = +0.240, t = 3.100 **
#     Consistent with site-level (b = +0.216)
#   Settlement gravity: b = -0.285, t = -3.006 **
#     Consistent with site-level (b = -0.251)
#   Chla:               b = -0.205, t = -2.187 *
#     Reaches significance at transect level (n = 243)
#     but not site level (n = 54) — power difference only,
#     direction consistent.

# ── ICC ───────────────────────────────────────────────────────
vc           <- as.data.frame(VarCorr(sens_t_pressure))
site_var     <- vc$vcov[vc$grp == "site"]
residual_var <- vc$vcov[vc$grp == "Residual"]
icc          <- site_var / (site_var + residual_var)

cat(sprintf("\nICC = %.3f — %.1f%% of variance attributable to
between-site differences beyond fixed predictors\n",
            icc, icc * 100))
# ICC = 0.205 — (1 | site) random intercept justified.

# ── Diagnostics on pressure model ────────────────────────────
par(mfrow = c(1, 2))
plot(fitted(sens_t_pressure), residuals(sens_t_pressure),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Transect pressure model: Residuals vs Fitted")
abline(h = 0, lty = 2, col = "grey60")
lines(lowess(fitted(sens_t_pressure),
             residuals(sens_t_pressure)), col = "red")
qqnorm(residuals(sens_t_pressure),
       main = "Transect pressure model: Q-Q Residuals")
qqline(residuals(sens_t_pressure), col = "red")
par(mfrow = c(1, 1))


# ============================================================
#  FIGURE 1 — Total biomass marginal effects
#
#  Best model: rugosity + chla + settlement gravity
#  All non-focal predictors held at mean (z = 0).
#  Observed data overlaid on fitted lines.
#  Shared y-axis across panels.
# ============================================================

col_pressure <- "#C0392B"
col_rugosity <- "#27AE60"

# ── Prediction grids ──────────────────────────────────────────
sg_grid <- data.frame(
  log_settlement_grav_sc = seq(
    min(total_model_data$log_settlement_grav_sc),
    max(total_model_data$log_settlement_grav_sc),
    length.out = 200),
  rugosity_sc = 0,
  log_chla_sc = 0)

sg_pred     <- predict(m_best_q2, newdata = sg_grid, se.fit = TRUE)
sg_grid$fit <- sg_pred$fit
sg_grid$lwr <- sg_pred$fit - 1.96 * sg_pred$se.fit
sg_grid$upr <- sg_pred$fit + 1.96 * sg_pred$se.fit

rug_grid <- data.frame(
  rugosity_sc = seq(
    min(total_model_data$rugosity_sc),
    max(total_model_data$rugosity_sc),
    length.out = 200),
  log_settlement_grav_sc = 0,
  log_chla_sc = 0)

rug_pred      <- predict(m_best_q2, newdata = rug_grid,
                         se.fit = TRUE)
rug_grid$fit  <- rug_pred$fit
rug_grid$lwr  <- rug_pred$fit - 1.96 * rug_pred$se.fit
rug_grid$upr  <- rug_pred$fit + 1.96 * rug_pred$se.fit

# ── Shared y-axis limits ──────────────────────────────────────
y_min <- min(c(sg_grid$lwr, rug_grid$lwr,
               total_model_data$log_mean_biomass), na.rm = TRUE)
y_max <- max(c(sg_grid$upr, rug_grid$upr,
               total_model_data$log_mean_biomass), na.rm = TRUE)
y_pad <- (y_max - y_min) * 0.05
y_lim <- c(y_min - y_pad, y_max + y_pad)

# ── Shared theme ──────────────────────────────────────────────
theme_fig1 <- theme_bw(base_size = 12) +
  theme(axis.title.x     = element_text(face = "bold"),
        axis.title.y     = element_blank(),
        panel.grid.minor = element_blank(),
        plot.tag         = element_text(face = "bold", size = 13),
        plot.margin      = margin(5, 10, 5, 2))

# ── Panel (a): Settlement gravity ─────────────────────────────
p1a <- ggplot(sg_grid,
              aes(x = log_settlement_grav_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = col_pressure, alpha = 0.15) +
  geom_line(colour = col_pressure, linewidth = 1.1) +
  geom_point(data = total_model_data,
             aes(x = log_settlement_grav_sc,
                 y = log_mean_biomass),
             colour = "grey40", size = 1.5, alpha = 0.5,
             inherit.aes = FALSE) +
  coord_cartesian(ylim = y_lim) +
  labs(x   = "Fishing pressure\n(settlement gravity, standardised)",
       tag = "(a)") +
  theme_fig1

# ── Panel (b): Rugosity ───────────────────────────────────────
p1b <- ggplot(rug_grid, aes(x = rugosity_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = col_rugosity, alpha = 0.15) +
  geom_line(colour = col_rugosity, linewidth = 1.1) +
  geom_point(data = total_model_data,
             aes(x = rugosity_sc, y = log_mean_biomass),
             colour = "grey40", size = 1.5, alpha = 0.5,
             inherit.aes = FALSE) +
  coord_cartesian(ylim = y_lim) +
  labs(x   = "Habitat complexity\n(rugosity, standardised)",
       tag = "(b)") +
  theme_fig1

# ── Shared y-axis label ───────────────────────────────────────
y_label <- textGrob(
  label = "Total reef\nfish biomass\n(log g)",
  rot = 0, vjust = 0.5, hjust = 0.5,
  gp = gpar(fontsize = 11, fontface = "bold"))

fig_title <- textGrob(
  label = paste("Figure 1. Marginal effects of fishing pressure",
                "and habitat complexity on total reef fish biomass"),
  x = 0, hjust = 0,
  gp = gpar(fontsize = 10, fontface = "plain"))

grid.arrange(
  fig_title,
  arrangeGrob(y_label, p1a, p1b,
              ncol = 3, widths = c(0.08, 0.46, 0.46)),
  nrow = 2, heights = c(0.06, 0.94))

# jpeg("figure1_total_biomass.jpg",
#      width = 22, height = 11, units = "cm", res = 300)
# grid.arrange(
#   fig_title,
#   arrangeGrob(y_label, p1a, p1b,
#               ncol = 3, widths = c(0.08, 0.46, 0.46)),
#   nrow = 2, heights = c(0.06, 0.94))
# dev.off()


# ============================================================
#  RESULTS SUMMARY
#  Quick reference for writing — verify all values match
#  reported results before writing up.
# ============================================================

results_summary <- tribble(
  ~Question,      ~Result,           ~Key_finding,
  "Q1 metric",    "Sett. gravity",   "weight = 0.826, DAICc = 4.39 vs baseline",
  "Q1 effect",    "Significant",     "b = -0.251, p = 0.012, 95% CI [-0.445, -0.057], fold diff = 2.73x",
  "Q1 rugosity",  "Significant",     "b = +0.216, p = 0.009, 95% CI [0.056, 0.376], fold diff = 3.04x",
  "Q1 chla",      "ns",              "b = -0.125, p = 0.208, 95% CI [-0.322, 0.072]",
  "Q1 fit",       "R² = 0.240",      "adj. R² = 0.194, F(3,50) = 5.253, p = 0.003",
  "Q2 main",      "Not supported",   "DAICc = 2.41, weight = 0.158",
  "Q2 int",       "Not supported",   "DAICc = 1.06, weight = 0.312, b = -0.225, p = 0.060, CI [-0.460, 0.010] — marginal, CIs overlap zero",
  "Q3 MPA",       "Not supported",   "DAICc = 3.01, weight = 0.182, low: b = -0.342 ns, medium: b = -0.084 ns",
  "Sensitivity",  "Consistent",      "Model ordering identical; rugosity b = +0.240, pressure b = -0.285, ICC = 0.205"
)

cat("\n--- Results summary ---\n")
print(results_summary)

# ── End of script ─────────────────────────────────────────────