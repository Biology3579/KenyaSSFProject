# ============================================================
#  DRIVERS OF CORALLIVORE BIOMASS
#  Chapter 1 — Functional Group Analysis: Corallivores
#
#  Analytical framework mirrors browser analysis.
#  Structure: Q1 (pressure) → Q2 (connectivity) → Q3 (MPA)
#
#  Scientific questions:
#
#  Q1 — Human pressure
#       Does human pressure explain variation in corallivore
#       biomass beyond local ecological context, and which
#       spatial metric best captures SSF exploitation
#       intensity for this functional group?
#
#  Q2 — Larval connectivity
#       Does larval connectivity explain additional variation
#       in corallivore biomass beyond the best Q1 model?
#       Connectivity x pressure interaction only testable
#       if pressure is supported in Q1.
#
#  Q3 — Formal protection
#       Does MPA status explain additional variation in
#       corallivore biomass beyond the best Q2 model?
#       MPA x pressure not testable if pressure not in
#       model. MPA x connectivity explored if MPA
#       supported as main effect.
#
#  Baseline model (fixed a priori, never tested):
#       biomass ~ rugosity_sc + log_chla_sc
#
#  Key difference from total biomass:
#       Corallivore biomass may have zeros at site level
#       (check below). Tweedie used for consistency
#       across all functional group analyses regardless
#       of zero proportion.
#
#  Sensitivity analyses:
#       (a) Alternative pressure metrics
#       (b) Transect-level replication (Tweedie GLMM)
# ============================================================

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

# ── Site-level dataset ────────────────────────────────────────
coralliv_model_data <- coralliv_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_coralliv_biomass,
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

cat("\nCorallivore model data:", nrow(coralliv_model_data), "sites\n")

# ── Data checks ───────────────────────────────────────────────
coralliv_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                log_chla_sc, connectivity_sc, mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

cat("\nZeros in mean_biomass:",
    sum(coralliv_model_data$mean_biomass == 0), "\n")
cat("Site-level zero proportion:",
    round(mean(coralliv_model_data$mean_biomass == 0), 3), "\n")
cat("\nResponse summary:\n")
print(summary(coralliv_model_data$mean_biomass))

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
#  Run on baseline model — consistent with all functional
#  group analyses. Tweedie used for consistency across
#  groups regardless of zero proportion.
# ============================================================

coralliv_lm_baseline <- lm(
  log(mean_biomass + 0.01) ~ rugosity_sc +
    log_chla_sc,
  data = coralliv_model_data
)

par(mfrow = c(2, 2))
plot(coralliv_lm_baseline, main = "Gaussian log — baseline")
par(mfrow = c(1, 1))

coralliv_tw_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

coralliv_tw_res <- simulateResiduals(coralliv_tw_baseline,
                                     n = 1000)
plot(coralliv_tw_res)
testZeroInflation(coralliv_tw_res)
testDispersion(coralliv_tw_res)

# ── Family selection decision ─────────────────────────────────
#
# Gaussian log:
#   Residuals vs Fitted: broadly flat, minor curve at
#     low fitted values. Sites 19 and 13 pull residuals
#     to -2 — lower tail deviation visible.
#   Q-Q: site 13 falls below theoretical line at lower
#     tail — moderate deviation but not severe.
#   Scale-Location: slight downward trend, sites 19 and
#     13 elevated. Acceptable but not ideal.
#   Residuals vs Leverage: site 13 approaches Cook's
#     distance threshold (leverage ~0.12). Site 9
#     influential but within bounds.
#   Adequate but lower tail sites warrant caution.
#
# Tweedie (log link): SELECTED
#   DHARMa diagnostics (n = 1000):
#     KS test:        p = 0.795 — excellent fit
#     Dispersion:     p = 0.804, ratio = 0.892 — acceptable
#     Zero inflation: ratio = 0, p = 1.000 — no zeros
#     Outlier test:   p = 1.000 — no outliers
#   QQ plot follows theoretical line closely throughout.
#   Residuals vs predicted: slight negative slope in
#     smoothing line but well within confidence band.
#     No significant problems detected by DHARMa.
#   Selected for fit quality and consistency across
#   all functional group analyses.
#
#  Proceed: glmmTMB(family = tweedie(link = "log")) on
#  raw mean_biomass throughout all corallivore analyses.

# ============================================================
#  RANDOM EFFECT STRUCTURE
#  Tested on baseline model — same rationale as total
#  biomass and browser analyses.
# ============================================================

coralliv_re_null <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

coralliv_re_ecoregion <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | ecoregion),
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

cat("\n--- Corallivore RE structure: ecoregion ---\n")
print(make_aicc_df(list(
  "No RE"           = coralliv_re_null,
  "(1 | ecoregion)" = coralliv_re_ecoregion
)))

# Random effect structure: ecoregion
# Tested on baseline model to avoid circularity.
# No RE: AICc = 560.23, weight = 0.781
# (1 | ecoregion): DAICc = 2.54, weight = 0.220
# Ecoregion RE not supported (DAICc = 2.54, weight = 0.220).
# Not pursued — only 4 ecoregions with severely uneven
# group sizes (n = 2, 8, 9, 35; Gelman & Hill 2007).
# Consistent with total biomass (DAICc = 2.25) and
# browser (DAICc = 2.13) results.
# All corallivore models fitted without RE throughout.

# ── Variance inflation factors ────────────────────────────────
# Confirms pairwise correlations do not translate into
# meaningful variance inflation in multivariate models.
# Checked at each stage predictors are added.
# Note: no pressure stage — pressure not supported in Q1.
# Sequence: baseline → baseline + settgrav (Q1 check) →
#           connectivity → connectivity + MPA.

cat("\n--- VIF: baseline ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data))

cat("\n--- VIF: baseline + settlement gravity (Q1 check) ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data))

cat("\n--- VIF: baseline + connectivity ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + connectivity_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data))

cat("\n--- VIF: baseline + connectivity + MPA ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + connectivity_sc + mpa_status,
  family = tweedie(link = "log"),
  data   = coralliv_model_data))

# ============================================================
#  Q1 — HUMAN PRESSURE
#
#  Scientific question:
#  Does human pressure explain variation in corallivore
#  biomass beyond local ecological context, and which
#  spatial metric best captures SSF exploitation intensity?
#
#  A priori prediction: settlement gravity outperforms
#  alternatives — consistent with total biomass Q1.
#  If baseline best supported (as for browsers), no
#  pressure term carried forward.
#
#  Two steps:
#  Step 1 — AICc comparison: does any metric outperform
#            ecological baseline?
#  Step 2 — Coefficient check: direction and significance
#            of selected metric (or direction checks only
#            if baseline best supported).
# ============================================================

c_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

# ── Q1 Step 1: Metric comparison ─────────────────────────────
c_q1_settgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

c_q1_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

c_q1_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_pop_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

cat("\n--- Q1 Step 1: Corallivore metric comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                      = c_baseline,
  "Baseline + settlement gravity" = c_q1_settgrav,
  "Baseline + market gravity"     = c_q1_mktgrav,
  "Baseline + settlement pop."    = c_q1_settpop
)))

# ── Q1 Step 2: Coefficient check ─────────────────────────────
cat("\n--- Q1 Step 2: Baseline coefficients ---\n")
summary(c_baseline)

# ── Q1: Rugosity and chla ranges ─────────────────────────────
cat("\n--- Corallivore Q1: baseline confint ---\n")
print(confint(c_baseline))

cat("\n--- Q1: Pressure metric direction checks ---\n")
cat("Settlement gravity:\n")
print(summary(c_q1_settgrav)$coefficients$cond[
  "log_settlement_grav_sc", ])
cat("\nMarket gravity:\n")
print(summary(c_q1_mktgrav)$coefficients$cond[
  "log_market_gravity_sc", ])
cat("\nSettlement population:\n")
print(summary(c_q1_settpop)$coefficients$cond[
  "log_settlement_pop_sc", ])

# Q1 results:
#   Baseline:           AICc = 560.23, weight = 0.499 (BEST)
#   Settlement gravity: DAICc = 1.86,  weight = 0.197
#   Settlement pop.:    DAICc = 2.26,  weight = 0.161
#   Market gravity:     DAICc = 2.50,  weight = 0.143
#
#   Baseline best supported — no pressure metric
#   outperforms ecological context for corallivores.
#   High model selection uncertainty across all metrics
#   (all within DAICc < 2.5) — no best metric identifiable.
#
# Baseline coefficients:
#   Rugosity: b = +0.105, p = 0.244 ns
#     Positive but not significant — habitat structural
#     complexity is a weak driver of corallivore biomass.
#     Corallivores depend more on live coral availability
#     than reef architectural complexity per se.
#   Chla:     b = +0.190, p = 0.049 *
#     Significant and positive — productive waters
#     support higher corallivore biomass, consistent
#     with corallivores' dependence on live coral which
#     is associated with clear, nutrient-rich conditions.
#     Contrasts with total biomass (chla negative once
#     pressure controlled) and browsers (chla negative,
#     ns) — ecologically coherent group-specific pattern.
#
# Pressure metric direction checks:
#   Settlement gravity: b = +0.088, p = 0.409 ns
#   Market gravity:     b = -0.019, p = 0.851 ns
#   Settlement population: b = -0.050, p = 0.601 ns
#
#   Directions inconsistent across metrics — settlement
#   gravity positive, market gravity and settlement
#   population weakly negative. All non-significant.
#   No coherent pressure signal for corallivores —
#   consistent with this group being relatively less
#   targeted by SSFs than other functional groups.
#   Absence of pressure effect not metric-dependent.
#
# Best Q1 model: c_baseline (rugosity + chla)
# No pressure term carried forward.
c_best_q1 <- c_baseline
# ============================================================
#  Q2 — LARVAL CONNECTIVITY
#
#  Does connectivity explain additional variation in
#  corallivore biomass beyond the ecological baseline?
#
#  Main effect only — connectivity x pressure interaction
#  not testable because pressure is not a supported
#  predictor for corallivores (Q1: baseline best supported,
#  weight = 0.499). Interaction requires pressure as a
#  main effect first.
#
#  Tweedie family throughout.
# ============================================================

c_q2_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

cat("\n--- Q2: Corallivore connectivity main effect ---\n")
print(make_aicc_df(list(
  "Best Q1"        = c_best_q1,
  "Best Q1 + conn" = c_q2_conn
)))

cat("\n--- Q2: Connectivity coefficients ---\n")
summary(c_q2_conn)

# ── Q2: Connectivity range and fold difference ───────────────
cat("\n--- Corallivore Q2: connectivity range and fold difference ---\n")
conn_range_c <- range(coralliv_model_data$connectivity_sc, na.rm = TRUE)
cat("Range:", conn_range_c, "\n")
conn_span_c  <- diff(conn_range_c)
cat(sprintf("Span: %.3f SD units\n", conn_span_c))
b_conn_c     <- -0.221
fold_conn_c  <- exp(abs(b_conn_c * conn_span_c))
cat(sprintf("Fold difference (low vs high connectivity): %.2fx\n", fold_conn_c))

cat("\n--- Corallivore Q2: confint ---\n")
print(confint(c_q2_conn))

# Q2 results:
#   Best Q1 + conn: AICc = 557.13, weight = 0.825 (BEST)
#   Best Q1:        DAICc = 3.10,  weight = 0.175
#   Connectivity strongly supported (weight = 0.825).
#
#   Connectivity: b = -0.221, p = 0.015 *
#     Significant NEGATIVE effect — well-connected sites
#     have LOWER corallivore biomass.
#     Opposite direction to browsers (b = +0.421).
#   Rugosity:     b = +0.108, p = 0.214 ns — stable
#   Chla:         b = +0.121, p = 0.196 ns — weakens
#     slightly from baseline (was +0.190, p = 0.049)
#     once connectivity included.
#
#   Negative connectivity signal is biologically
#   interpretable — corallivores may be more abundant
#   at isolated reefs where reduced larval connectivity
#   means lower predator and competitor recruitment,
#   or where isolated reefs experience lower fishing
#   pressure due to reduced accessibility. The negative
#   direction contrasts sharply with browsers (positive)
#   and total biomass (near zero) — connectivity effects
#   are group-specific and reflect different ecological
#   dependencies across functional groups.
#
#   Note: connectivity correlates positively with market
#   gravity (r = +0.30) — well-connected sites tend to
#   be more accessible. For corallivores this may mean
#   the connectivity signal partly reflects accessibility
#   rather than pure larval replenishment. Interpret
#   with caution.
#
# Best Q2 model: rugosity + chla + connectivity
c_best_q2 <- c_q2_conn

# ============================================================
#  Q3 — FORMAL PROTECTION
#
#  Does MPA status explain additional variation in
#  corallivore biomass beyond the ecological baseline
#  and connectivity?
#
#  MPA tested despite absence of pressure signal in Q1 —
#  formal harvest exclusion operates independently of
#  the broader exploitation pressure gradient.
#
#  MPA x pressure not testable — pressure not supported
#  in Q1. MPA x connectivity explored if MPA supported
#  as main effect — do protected sites show stronger
#  connectivity effects?
#
#  Tweedie family throughout.
# ============================================================

c_q3_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

cat("\n--- Q3: Corallivore MPA main effect ---\n")
print(make_aicc_df(list(
  "Best Q2"       = c_best_q2,
  "Best Q2 + MPA" = c_q3_mpa
)))

cat("\n--- Q3: MPA coefficients ---\n")
summary(c_q3_mpa)

# ── Q3: MPA placement check ──────────────────────────────────
cat("\n--- Corallivore Q3: pressure and connectivity by MPA ---\n")
coralliv_model_data %>%
  group_by(mpa_status) %>%
  summarise(
    n                    = n(),
    mean_biomass         = round(mean(mean_biomass), 1),
    mean_settlement_grav = round(mean(log_settlement_grav_sc), 3),
    mean_connectivity    = round(mean(connectivity_sc), 3),
    .groups = "drop"
  ) %>%
  print()

# ── MPA data structure check ──────────────────────────────────
cat("\n--- Q3: MPA site distribution ---\n")
coralliv_model_data %>%
  group_by(mpa_status) %>%
  summarise(
    n                 = n(),
    mean_biomass      = round(mean(mean_biomass), 1),
    mean_connectivity = round(mean(connectivity_sc), 3),
    min_connectivity  = round(min(connectivity_sc),  3),
    max_connectivity  = round(max(connectivity_sc),  3),
    .groups = "drop"
  ) %>%
  print()

# Q3 results:
#   Best Q2:       AICc = 557.13, weight = 0.549
#   Best Q2 + MPA: DAICc = 0.39,  weight = 0.451
#   Genuine model selection uncertainty — MPA does not
#   substantially improve on connectivity baseline.
#
#   Low MPA:    b = +0.622, p = 0.027 * — significant
#     but artefactual. Low MPA sites confined to narrow
#     high-connectivity range (z = 0.66 to max, n = 7)
#     with characteristics that elevate biomass
#     independently of protection status.
#   Medium MPA: b = +0.248, p = 0.192 ns
#     Not significant. Raw mean (68.4g) essentially
#     identical to unprotected sites (67.6g).
#   Connectivity: b = -0.337, p = 0.001 ***
#     Strengthens and remains highly significant once
#     MPA included — negative connectivity signal robust
#     and not confounded by MPA placement.
#
#   Raw biomass differences negligible across MPA
#   categories (none = 67.6g, low = 78.9g,
#   medium = 68.4g) — no biologically meaningful
#   protection signal for corallivores.
#
#   Conclusion: MPA not supported as a driver of
#   corallivore biomass (DAICc = 0.39, weight = 0.451).
#   Low MPA coefficient reflects site characteristics
#   not genuine protection benefit. Connectivity model
#   retained as best supported throughout Q3.

# ── Q3 MPA × connectivity interaction ────────────────────────
c_q3_mpa_conn_int <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    mpa_status * connectivity_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

cat("\n--- Q3: MPA x connectivity interaction ---\n")
print(make_aicc_df(list(
  "Best Q2"              = c_best_q2,
  "Best Q2 + MPA"        = c_q3_mpa,
  "Best Q2 + MPA x conn" = c_q3_mpa_conn_int
)))

cat("\n--- Q3: MPA x connectivity coefficients ---\n")
summary(c_q3_mpa_conn_int)

# Q3 MPA x connectivity interaction:
#   Best Q2:              AICc = 557.13, weight = 0.525
#   Best Q2 + MPA:        DAICc = 0.39,  weight = 0.432
#   Best Q2 + MPA x conn: DAICc = 4.97,  weight = 0.044
#
#   Interaction not supported (DAICc = 4.97, weight = 0.044)
#   — adding the interaction substantially worsens model
#   fit relative to both simpler models. MPA effectiveness
#   does not vary with connectivity for corallivores.
#   Confirms connectivity model as best supported
#   throughout Q1-Q3.
#
# Final model sequence:
#   Q1: baseline best (weight = 0.499) — no pressure
#   Q2: connectivity supported (weight = 0.825,
#       DAICc = 3.10, b = -0.221, p = 0.015)
#   Q3: MPA not supported (DAICc = 0.39) — connectivity
#       model retained throughout
#
# Best model: rugosity + chla + connectivity
c_best_q3 <- c_best_q2

cat("\n--- Corallivore: predicted vs observed ---\n")
pred_c <- predict(c_best_q3, type = "response")
obs_c  <- coralliv_model_data$mean_biomass
cat(sprintf("Pearson r: %.3f\n", cor(pred_c, obs_c)))

cat("\n--- Q3: DHARMa diagnostics ---\n")
coralliv_sim <- simulateResiduals(c_best_q3, n = 1000)
plot(coralliv_sim)
testOutliers(coralliv_sim)

# ============================================================
#  SPATIAL AUTOCORRELATION
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

cat("\n--- Spatial autocorrelation: corallivore best model ---\n")
print(moran.test(residuals(c_best_q3, type = "pearson"),
                 listw5_c))

# Spatial autocorrelation: corallivore best model
# (rugosity + chla + connectivity)
# Moran's I = -0.021, p = 0.509 — no significant
# spatial autocorrelation in residuals.
#
# Consistent with browser result (I = -0.078, p = 0.811)
# and contrasts with total biomass (I = 0.140, p = 0.015).
# Connectivity predictor appears to adequately capture
# the spatial structure in corallivore biomass without
# leaving a residual geographic signal.
# No spatial error modelling required.


# ============================================================
#  SENSITIVITY ANALYSIS
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
c_sens_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

c_sens_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_pop_sc,
  family = tweedie(link = "log"),
  data   = coralliv_model_data
)

cat("\n--- Sensitivity (a): corallivore alternative metrics ---\n")
cat("Market gravity:\n")
print(summary(c_sens_mktgrav)$coefficients$cond)
cat("\nSettlement population:\n")
print(summary(c_sens_settpop)$coefficients$cond)

# Sensitivity (a): alternative pressure metrics
# Purpose: confirm Q1 null result is not metric-dependent.
#
# Market gravity:      b = -0.019, p = 0.851 ns
# Settlement pop.:     b = -0.050, p = 0.601 ns
#
# Both non-significant and near zero — consistent with
# Q1 result (settlement gravity b = +0.088, p = 0.409).
# Note: settlement gravity positive, market gravity and
# settlement population negative — inconsistent directions
# confirm genuine absence of a coherent pressure signal
# for corallivores rather than a detection issue.
# Chla remains positive across all three metric models
# (b = 0.181-0.190) — productivity signal stable
# regardless of pressure metric choice.
# Rugosity non-significant throughout (b = 0.094-0.108).
#
# Q1 null conclusion robust to metric choice.


# ── (b) Transect-level replication ───────────────────────────
# Mirrors site-level Q1-Q3 sequence at transect level.
# Standard Tweedie + (1|site) throughout — consistent
# with site-level family selection.

c_trans_null <- glmmTMB(
  transect_coralliv_biomass ~ 1 + (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data
)

c_trans_baseline <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data
)

# Pressure — for completeness, not expected to be supported
c_trans_pressure <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data
)

# Mirrors Q2 — connectivity against baseline
c_trans_conn <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data
)

# Mirrors Q3 — MPA against baseline + connectivity
c_trans_mpa <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data
)

cat("\n--- Sensitivity (b): corallivore transect comparison ---\n")
print(make_aicc_df(list(
  "Null"             = c_trans_null,
  "Baseline"         = c_trans_baseline,
  "Baseline + press" = c_trans_pressure,
  "Baseline + conn"  = c_trans_conn,
  "Best + MPA"       = c_trans_mpa
)))

cat("\n--- Sensitivity (b): connectivity coefficients ---\n")
summary(c_trans_conn)

cat("\n--- Sensitivity (b): connectivity + MPA coefficients ---\n")
summary(c_trans_mpa)

vc_c <- VarCorr(c_trans_conn)
site_sd_c <- sqrt(as.numeric(vc_c$cond$site))
cat(sprintf("\nSite random intercept SD = %.3f\n", site_sd_c))

# ============================================================
#  MARGINAL EFFECT PLOTS
#  Best model: rugosity + chla + connectivity (c_best_q3)
#  MPA not supported (DAICc = 0.39) — MPA plot omitted.
#  Rugosity and connectivity plotted only.
#  All non-focal predictors held at 0 (their mean).
#  Predictions on response scale (raw biomass).
# ============================================================

best_model_c <- c_best_q3  # rugosity + chla + connectivity

# ── Rugosity effect ───────────────────────────────────────────
rug_grid_c <- data.frame(
  rugosity_sc     = seq(
    min(coralliv_model_data$rugosity_sc),
    max(coralliv_model_data$rugosity_sc),
    length.out = 200),
  log_chla_sc     = 0,
  connectivity_sc = 0
)

rug_pred_c     <- predict(best_model_c,
                          newdata = rug_grid_c,
                          se.fit  = TRUE,
                          type    = "response",
                          re.form = NA)
rug_grid_c$fit <- rug_pred_c$fit
rug_grid_c$lwr <- rug_pred_c$fit - 1.96 * rug_pred_c$se.fit
rug_grid_c$upr <- rug_pred_c$fit + 1.96 * rug_pred_c$se.fit

p_c_rugosity <- ggplot(rug_grid_c,
                       aes(x = rugosity_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "#2c7bb6", alpha = 0.15) +
  geom_line(colour = "#2c7bb6", linewidth = 1.1) +
  geom_point(data = coralliv_model_data,
             aes(x = rugosity_sc, y = mean_biomass),
             colour = "grey40", size = 1.5,
             alpha  = 0.5, inherit.aes = FALSE) +
  labs(x = "Rugosity (standardised)",
       y = "Corallivore biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())


# ── Connectivity effect ───────────────────────────────────────
# Connectivity supported in Q2 (weight = 0.825, b = -0.221)
# Negative — well-connected sites have lower corallivore
# biomass. Contrasts with browsers (positive direction).

conn_grid_c <- data.frame(
  connectivity_sc = seq(
    min(coralliv_model_data$connectivity_sc),
    max(coralliv_model_data$connectivity_sc),
    length.out = 200),
  rugosity_sc = 0,
  log_chla_sc = 0
)

conn_pred_c      <- predict(best_model_c,
                            newdata = conn_grid_c,
                            se.fit  = TRUE,
                            type    = "response",
                            re.form = NA)
conn_grid_c$fit  <- conn_pred_c$fit
conn_grid_c$lwr  <- conn_pred_c$fit - 1.96 * conn_pred_c$se.fit
conn_grid_c$upr  <- conn_pred_c$fit + 1.96 * conn_pred_c$se.fit

p_c_conn <- ggplot(conn_grid_c,
                   aes(x = connectivity_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "#2c7bb6", alpha = 0.15) +
  geom_line(colour = "#2c7bb6", linewidth = 1.1) +
  geom_point(data = coralliv_model_data,
             aes(x = connectivity_sc, y = mean_biomass),
             colour = "grey40", size = 1.5,
             alpha  = 0.5, inherit.aes = FALSE) +
  labs(x = "Connectivity (standardised)",
       y = "Corallivore biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ── Arrange plots ─────────────────────────────────────────────
gridExtra::grid.arrange(p_c_rugosity, p_c_conn, ncol = 2)

# jpeg("corallivore_marginal_effects.jpg",
#      width = 22, height = 11, units = "cm", res = 300)
# gridExtra::grid.arrange(p_c_rugosity, p_c_conn, ncol = 2)
# dev.off()


# ============================================================
#  RESULTS SUMMARY
# ============================================================

cat("\n--- Corallivore results summary ---\n")
tribble(
  ~Question,   ~Result,          ~Key_finding,
  "Q1",        "Baseline best",  "weight = 0.499, no pressure metric outperforms baseline; chla b = +0.190, p = 0.049",
  "Q2 conn",   "Supported",      "weight = 0.825, DAICc = 3.10, b = -0.221, p = 0.015 — negative, contrasts with browsers",
  "Q3 MPA",    "Not supported",  "DAICc = 0.39, connectivity model retained; low MPA artefactual",
  "Sens (a)",  "Consistent",     "all metrics ns, inconsistent directions — Q1 null robust",
  "Sens (b)",  "Consistent",     "connectivity negative at both levels (b = -0.231, p = 0.009); MPA ambiguous"
) %>% print()


# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()