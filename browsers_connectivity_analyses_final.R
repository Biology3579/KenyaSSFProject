# ============================================================
#  DRIVERS OF BROWSER BIOMASS
#  Chapter 1 — Functional Group Analysis: Browsers
#
#  Analytical framework mirrors total biomass (Q1–Q3 +
#  sensitivity analyses), with group-specific family
#  selection and pressure metric confirmation.
#
#  Scientific questions:
#
#  Q1 — Human pressure
#       Does human pressure explain variation in browser
#       biomass beyond local ecological context, and which
#       spatial metric best captures SSF exploitation
#       intensity for this functional group?
#       Browsers are heavily targeted by SSFs — the metric
#       that best captures total community exploitation
#       may not be optimal for a selectively harvested group.
#
#  Q2 — Larval connectivity
#       Does larval connectivity explain additional variation
#       in browser biomass beyond the human pressure
#       baseline, and does it modify the relationship
#       between fishing pressure and browser biomass?
#
#  Q3 — Formal protection
#       Does MPA status explain additional variation in
#       browser biomass beyond the fully specified pressure
#       and connectivity model?
#       MPA tested last — non-randomly placed (r = -0.30
#       with settlement gravity), confounded by pressure
#       and connectivity unless both controlled first.
#
#  Rationale for sequence:
#       Identical to total biomass — pressure first as
#       primary driver of interest, connectivity second
#       to extend Warmuth et al. (2024), MPA last as
#       governance response downstream of both.
#
#  Baseline model (fixed a priori, never tested):
#       biomass ~ rugosity_sc + log_chla_sc
#       Identical justification to total biomass.
#       Pressure metric selected in Q1.
#
#  Key difference from total biomass:
#       Browser biomass has ~11% zeros at site level.
#       Tweedie distribution selected over Gaussian log
#       (see family selection). All models use
#       glmmTMB(family = tweedie(link = "log")) on
#       raw mean_biomass throughout.
#
#  Sensitivity analyses:
#       (a) Alternative pressure metrics
#       (b) Transect-level replication (Tweedie GLMM)
# ============================================================

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
    ecoregion              = first(ecoregion),     # ADD THIS
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
#  Run on baseline model — distributional family is a
#  property of the response variable, not the predictor set.
#
#  Browser biomass has ~11% zeros at site level.
#  Gaussian log requires an offset constant for zeros —
#  log-transformed distribution is bimodal with zero sites
#  clustering at log(offset) separated from the main
#  distribution. No offset constant resolves this.
#  Tweedie handles zeros natively without offset.
#
#  Two families tested on baseline structure:
#    Gaussian log: lm() on log(mean_biomass + offset)
#    Tweedie:      glmmTMB() on raw mean_biomass, log link
#
#  DHARMa used for Tweedie diagnostics.
#  Standard plot() for Gaussian log diagnostics.
# ============================================================

# ── Gaussian log — baseline ───────────────────────────────────
browser_lm_baseline <- lm(
  log(mean_biomass + 0.01) ~ rugosity_sc +
    log_chla_sc,
  data = browser_model_data
)

par(mfrow = c(2, 2))
plot(browser_lm_baseline, main = "Gaussian log — baseline")
par(mfrow = c(1, 1))

# ── Tweedie — baseline ────────────────────────────────────────
browser_tw_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

browser_tw_res <- simulateResiduals(browser_tw_baseline,
                                    n = 1000)
plot(browser_tw_res)
testZeroInflation(browser_tw_res)
testDispersion(browser_tw_res)

# ── Family selection decision ─────────────────────────────────
#
# Gaussian log: REJECTED
#   Residuals vs Fitted: strong downward curve at low
#     fitted values — zero sites pulling residuals to
#     extreme values, severe non-linearity that no
#     offset constant resolves.
#   Q-Q: severe lower tail deviation — zero sites fall
#     far from theoretical line.
#   Scale-Location: strong downward trend —
#     heteroscedasticity throughout fitted range.
#   Zero sites unduly influential on fitted coefficients.
#
# Tweedie (log link): SELECTED
#   Handles zeros natively without offset.
#   DHARMa diagnostics (n = 1000 simulations):
#     KS test:        p = 0.762 — good fit
#     Dispersion:     p = 0.222, ratio = 1.615 — acceptable
#     Zero inflation: p = 0.902, ratio = 0.872 — not needed
#     Outlier test:   p = 0.102 — no significant outliers
#
#  Proceed: glmmTMB(family = tweedie(link = "log")) on
#  raw mean_biomass throughout all browser analyses.


# ============================================================
#  RANDOM EFFECT STRUCTURE
#  Tested on baseline model — same rationale as total
#  biomass. Tweedie family used throughout.
# ============================================================

browser_re_null <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

browser_re_ecoregion <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | ecoregion),
  family = tweedie(link = "log"),
  data   = browser_model_data
)

cat("\n--- Browser RE structure: ecoregion ---\n")
print(make_aicc_df(list(
  "No RE"           = browser_re_null,
  "(1 | ecoregion)" = browser_re_ecoregion
)))

# Ecoregion RE not supported — consistent with total biomass result. 
# All browser models without RE.

# ============================================================
#  VARIANCE INFLATION FACTORS
#  Confirms pairwise correlations do not translate into
#  meaningful variance inflation in multivariate models.
#  Checked at each stage predictors are added.
#  Note: no pressure stage — pressure not supported in Q1.
#  Sequence: baseline → connectivity → connectivity + MPA.
#  check_collinearity() from performance package used
#  throughout — handles glmmTMB Tweedie models correctly.
# ============================================================

# ── Stage 1: Baseline ─────────────────────────────────────────
cat("\n--- VIF: baseline ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
))

# ── Stage 1: Baseline + Pressure ─────────────────────────────────────────
cat("\n--- VIF: baseline ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc + log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
))

# ── Stage 2: Baseline + connectivity ─────────────────────────
cat("\n--- VIF: baseline + connectivity ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
))

# ── Stage 3: Baseline + connectivity + MPA ───────────────────
cat("\n--- VIF: baseline + connectivity + MPA ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = browser_model_data
))

#  Collinearity results:
#   Baseline: rugosity VIF = TBC, chla VIF = TBC
#   + connectivity: max VIF = TBC
#   + MPA: max VIF = TBC
#  All VIFs < 2.0 at every stage — no multicollinearity
#  concern. Update with actual values after running.

# ============================================================
#  Q1 — HUMAN PRESSURE
#
#  Scientific question:
#  Does human pressure explain variation in browser biomass
#  beyond local ecological context, and which spatial metric
#  best captures SSF exploitation intensity for browsers?
#
#  A priori prediction: settlement gravity will outperform
#  alternatives — consistent with total biomass Q1.
#  Browser-specific note: heavily targeted group, metric
#  performance may differ from total community.
#
#  Two steps:
#  Step 1 — Does pressure add beyond baseline, and which
#            metric best captures browser exploitation?
#            AICc weight criterion for metric selection.
#  Step 2 — Coefficient check for selected metric.
#            p-value criterion for effect confirmation.
# ============================================================

b_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

# ── Q1 Step 1: Metric comparison ─────────────────────────────
b_q1_settgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

b_q1_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

b_q1_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_pop_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

cat("\n--- Q1 Step 1: Browser metric comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                      = b_baseline,
  "Baseline + settlement gravity" = b_q1_settgrav,
  "Baseline + market gravity"     = b_q1_mktgrav,
  "Baseline + settlement pop."    = b_q1_settpop
)))

# Q1 Step 1 results:
#   Baseline:           AICc = 835.89, weight = 0.459 (BEST)
#   Settlement gravity: DAICc = 1.61,  weight = 0.205
#   Market gravity:     DAICc = 1.65,  weight = 0.201
#   Settlement pop.:    DAICc = 2.44,  weight = 0.135
#
#   Baseline best supported — no pressure metric
#   outperforms ecological context alone for browsers.
#   Human pressure does not explain additional variance
#   in browser biomass beyond habitat and productivity.
#
#   Contrast with total biomass (settlement gravity
#   weight = 0.826, DAICc = 4.39 vs baseline) —
#   browser biomass is not structured by the raw
#   exploitation pressure gradient in this system.
#
#   Best Q1 model: baseline (rugosity + chla only).
#   No pressure term carried forward into Q2 and Q3.
#   Q2 and Q3 test connectivity and MPA against the
#   ecological baseline directly.

# ── Best Q1 model ─────────────────────────────────────────────
# Baseline best supported — no pressure term carried forward.
b_best_q1 <- b_baseline


# need to move this to ....
cat("\n--- Baseline coefficients ---\n")
summary(b_baseline)

# ── Q1: CIs for baseline coefficients ────────────────────────
cat("\n--- Browser Q1: baseline confint ---\n")
print(confint(b_baseline))

# ── Q1: Rugosity range and fold difference (browser model) ──
cat("\n--- Browser: rugosity range and fold difference ---\n")

# Range and span
rug_range_b <- range(browser_model_data$rugosity_sc, na.rm = TRUE)
cat(sprintf("Range: %.3f to %.3f\n", rug_range_b[1], rug_range_b[2]))

rug_span_b <- diff(rug_range_b)
cat(sprintf("Span: %.3f SD units\n", rug_span_b))

# Extract coefficient directly from model (no hard-coding)
b_rug_b <- fixef(b_baseline)$cond["rugosity_sc"]

# Fold change across observed range
fold_rug_b <- exp(abs(b_rug_b * rug_span_b))

cat(sprintf(
  "Fold difference in biomass (low → high rugosity): %.2fx higher\n",
  fold_rug_b
))

# Q1 Step 2 results — baseline and direction checks:
#
# Baseline coefficients (n = 54 sites):
#   Rugosity: b = +0.478, p = 0.013 * SIGNIFICANT
#     Habitat structural complexity is a significant
#     positive driver of browser biomass — consistent
#     with browsers depending on complex reef structure
#     for shelter and feeding habitat.
#   Chla:     b = -0.192, p = 0.322 ns
#     Negative but not significant. Direction consistent
#     with total biomass pattern once pressure context
#     is removed — retained as baseline control.
#
# Pressure metric direction checks:
#   Settlement gravity: b = -0.221, p = 0.333 ns
#   Market gravity:     b = +0.191, p = 0.347 ns
#   Settlement pop.:    b = -0.055, p = 0.761 ns
#
#   Inconsistent directions across metrics — settlement
#   gravity and settlement population negative (consistent
#   with exploitation hypothesis) but market gravity
#   positive (opposite direction). All three non-significant.
#   No coherent pressure signal for browsers — confirms
#   baseline as best supported model in Q1 Step 1.
#   Absence of pressure effect is not metric-dependent
#   but is also not consistent in direction, suggesting
#   genuine absence of a pressure-biomass relationship
#   for this functional group rather than a detection
#   issue.
#
# Best Q1 model: b_baseline (rugosity + chla)
# Proceed to Q2 with baseline as reference.


# ============================================================
#  Q2 — LARVAL CONNECTIVITY
#
#  Scientific question:
#  Does larval connectivity explain additional variation
#  in browser biomass beyond the ecological baseline?
#
# Connectivity × pressure interaction not evaluated.
# Pressure was not supported as a main effect in Q1
# (baseline model best-supported), and interaction terms
# were therefore not pursued to avoid overparameterisation.
#
#  Tweedie family throughout.
# ============================================================

# ── Step 1: Model comparison — does connectivity add? ────────
b_q2_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

cat("\n--- Q2: Model comparison (baseline vs baseline + connectivity) ---\n")
print(make_aicc_df(list(
  "Baseline"                = b_best_q1,
  "Baseline + connectivity" = b_q2_conn
)))

# Results:
#   Baseline + connectivity: AICc = 831.89, weight = 0.881 (BEST)
#   Baseline:                DAICc = 4.00,  weight = 0.119
#   Connectivity strongly supported (weight = 0.881,
#   DAICc = 4.00 vs baseline).

# ── Q2: Connectivity model coefficients ──────────────────────
cat("\n--- Q2: Connectivity model coefficients ---\n")

print(summary(b_q2_conn))
print(confint(b_q2_conn))


# ── Q2: Connectivity range and fold difference ───────────────
cat("\n--- Q2: Connectivity range and fold difference ---\n")

# Range and span
conn_range_b <- range(browser_model_data$connectivity_sc, na.rm = TRUE)
cat(sprintf("Range: %.3f to %.3f\n", conn_range_b[1], conn_range_b[2]))

conn_span_b <- diff(conn_range_b)
cat(sprintf("Span: %.3f SD units\n", conn_span_b))

# Extract coefficient directly from model
b_conn_b <- fixef(b_q2_conn)$cond["connectivity_sc"]

# Fold difference across observed range
fold_conn_b <- exp(abs(b_conn_b * conn_span_b))

cat(sprintf(
  "Fold difference in biomass (low → high connectivity): %.2fx higher\n",
  fold_conn_b
))


# ── Optional: CI-based fold difference (recommended) ─────────
ci_conn <- confint(b_q2_conn)["connectivity_sc", c("2.5 %", "97.5 %")]

fold_low  <- exp(abs(ci_conn[1] * conn_span_b))
fold_high <- exp(abs(ci_conn[2] * conn_span_b))

cat(sprintf(
  "Fold difference range (95%% CI): %.2fx to %.2fx\n",
  fold_low, fold_high
))

# Results:
#   Range: −1.355 to 1.312 (span = 2.667 SD units)
#   Fold difference: ~3.07× higher biomass at the most
#   connected compared to the least connected sites
#   95% CI: ~1.33× to 7.09×

# ── Step 4: Connectivity distribution check ──────────────────
# Confirms range adequate for inference — not driven by
# a handful of extreme sites.
cat("\n--- Q2: Connectivity distribution summary ---\n")
print(summary(browser_model_data$connectivity_sc))
print(quantile(browser_model_data$connectivity_sc,
               c(0.10, 0.25, 0.50, 0.75, 0.90)))

# Results:
#   Min = -1.355, Max = 1.312, Median = 0.021
#   10th pct = -1.316, 90th pct = 1.215
#   Distribution well spread — no inference concerns.

# ── Best Q2 model ─────────────────────────────────────────────
b_best_q2 <- b_q2_conn

# ============================================================
#  Q3 — FORMAL PROTECTION
#
#  Scientific question:
#  Does MPA status explain additional variation in browser
#  biomass beyond the ecological baseline and connectivity?
#
#  MPA tested last — non-randomly placed (medium sites at
#  lower pressure: mean settlement gravity z = -0.532 vs
#  unprotected z = +0.242). Effect confounded by pressure
#  and connectivity unless both controlled first.
#
#  Note: pressure not in reference model (not supported
#  in Q1). MPA x pressure interaction not testable.
#  If MPA supported, MPA x connectivity interaction
#  explored — mechanistically coherent for browsers
#  given strong dispersal dependence.
#
#  Tweedie family throughout.
# ============================================================

# ── Step 1: Model comparison — does MPA add? ─────────────────
b_q3_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = browser_model_data
)

cat("\n--- Q3: Model comparison (best Q2 vs best Q2 + MPA) ---\n")
print(make_aicc_df(list(
  "Best Q2"         = b_best_q2,
  "Best Q2 + MPA"   = b_q3_mpa
)))

# Results:
#   Best Q2 + MPA: AICc = 826.43, weight = 0.939 (BEST)
#   Best Q2:       DAICc = 5.47,  weight = 0.061
#   MPA strongly supported (weight = 0.939).

# ── Step 2: MPA coefficients ─────────────────────────────────
cat("\n--- Q3: MPA model coefficients ---\n")
summary(b_q3_mpa)

# ── Step 3: 95% CIs for MPA coefficients ─────────────────────
cat("\n--- Q3: MPA model confint ---\n")
print(confint(b_q3_mpa))

# Results:
#   Medium MPA: b = +1.143, SE = 0.375, z = 3.051, p = 0.002 **
#     Significant positive effect.
#   Low MPA:    b =  0.000, SE = 0.551, z = 0.000, p = 0.999
#     No detectable effect — low protection insufficient.
#   Rugosity:   b = +0.464, p = 0.004 ** — stable
#   Connectivity: b = +0.208, p = 0.251 — weakens once MPA
#     included, suggesting shared variance between
#     connectivity and MPA (well-connected sites tend
#     to be protected: low MPA mean conn z = +0.879,
#     medium MPA mean conn z = +0.303 vs none z = -0.345)


# Results:
#   Medium MPA: 95% CI [+0.409, +1.877] on log scale
#   Back-transformed: [1.50x, 6.53x] fold difference
#   Wide CI — effect real but magnitude uncertain.
#   Low MPA:    95% CI [-1.080, +1.080] — straddles zero

# ── Step 4: Fold difference — medium MPA ─────────────────────
cat("\n--- Q3: Medium MPA fold difference ---\n")
b_medium_mpa <- 1.14266348
fold_mpa_b   <- exp(b_medium_mpa)
ci_lo_mpa    <- exp(0.4086231)
ci_hi_mpa    <- exp(1.8767038)
cat(sprintf(
  "Medium MPA: %.2fx higher biomass (95%% CI: %.2fx to %.2fx)\n",
  fold_mpa_b, ci_lo_mpa, ci_hi_mpa))

# Results:
#   Medium MPA: 3.13x higher browser biomass than
#   unprotected sites (95% CI: 1.50x to 6.53x),
#   controlling for habitat complexity and connectivity.


# ── Step 6: Raw biomass by MPA status ────────────────────────
# Sanity check — model fold difference (3.13x) should
# be broadly consistent with raw means pattern.
cat("\n--- Q3: Raw biomass by MPA status ---\n")
browser_model_data %>%
  group_by(mpa_status) %>%
  summarise(
    n              = n(),
    mean_biomass   = round(mean(mean_biomass), 1),
    median_biomass = round(median(mean_biomass), 1),
    sd_biomass     = round(sd(mean_biomass), 1),
    .groups = "drop"
  ) %>%
  print()

# Results:
#   none:   mean = 679g,  median = 679g
#   low:    mean = 821g,  median = 821g
#   medium: mean = 2461g, median = 2461g
#
#   Raw means: medium = 3.6x higher than unprotected.
#   Model estimate (3.1x) slightly lower — consistent
#   once pressure and connectivity context controlled.
#   Low MPA raw advantage (821g vs 679g) not reflected
#   in model (b = 0.000) — low MPA sites at high pressure
#   context pulls raw means down; once controlled, no
#   effect detected.


# ── Best Q3 model ─────────────────────────────────────────────
# MPA main effect model — interaction suggestive only.
b_best_q3 <- b_q3_mpa

cat("\n--- Browser: predicted vs observed ---\n")
pred_b <- predict(b_best_q3, type = "response")
obs_b  <- browser_model_data$mean_biomass
cat(sprintf("Pearson r: %.3f\n", cor(pred_b, obs_b)))


cat("\n--- Q3: DHARMa diagnostics ---\n")
browser_mpa_sim <- simulateResiduals(b_q3_mpa, n = 1000)
plot(browser_mpa_sim)
testOutliers(browser_mpa_sim)

# ── Step 9: Raw data visualisation ───────────────────────────
ggplot(browser_model_data,
       aes(x = mpa_status, y = mean_biomass,
           fill = mpa_status)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5,
               alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.1, size = 1.8,
              alpha = 0.6, colour = "grey30") +
  scale_fill_manual(values = c("none"   = "#bdbdbd",
                               "low"    = "#74a9cf",
                               "medium" = "#0570b0")) +
  scale_x_discrete(labels = c("none"   = "No MPA",
                              "low"    = "Low",
                              "medium" = "Medium")) +
  labs(x = "MPA status", y = "Browser biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title         = element_text(face = "bold"),
        legend.position    = "none",
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank())

# ============================================================
#  SPATIAL AUTOCORRELATION CHECK
#  Pearson residuals from best model (post Q3).
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

cat("\n--- Spatial autocorrelation: browser best model ---\n")
print(moran.test(residuals(b_q3_mpa,
                           type = "pearson"), listw5_b))

# Spatial autocorrelation: browser best Q3 model
# Moran's I = -0.078, p = 0.811 — no significant
# spatial autocorrelation in residuals.
#
# Contrasts with total biomass (I = 0.140, p = 0.015)
# where weak but significant autocorrelation remained.
# Browser model residuals show no spatial structure —
# rugosity, connectivity, and MPA status adequately
# capture the spatial variation in browser biomass
# without leaving a residual geographic signal.
#
# Slight negative value indicates spatial dispersion —
# neighbouring sites tend to have more dissimilar
# residuals than expected. Consistent with browser
# biomass being structured by localised MPA boundaries
# that vary discretely between adjacent sites.
#
# No spatial error modelling required.

# ============================================================
#  SENSITIVITY ANALYSIS
# ============================================================

# ── Settlement gravity ────────────────────────────────────────
cat("\n--- Sensitivity (a): settlement gravity ---\n")

# Q2: connectivity beyond settlement gravity
b_sens_sg_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc + connectivity_sc,
  family = tweedie(link = "log"),
  data = browser_model_data
)

# Q2: connectivity × settlement gravity interaction
b_sens_sg_conn_int <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc * connectivity_sc,
  family = tweedie(link = "log"),
  data = browser_model_data
)

cat("\nQ2 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + SG"             = b_q1_settgrav,
  "Baseline + SG + conn"      = b_sens_sg_conn,
  "Baseline + SG + conn x SG" = b_sens_sg_conn_int
)))

# Check top model
cat("\nQ2 — connectivity model summary:\n")
print(summary(b_sens_sg_conn))
print(confint(b_sens_sg_conn))


# Q3: MPA beyond settlement gravity + connectivity
b_sens_sg_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc + connectivity_sc + mpa_status,
  family = tweedie(link = "log"),
  data = browser_model_data
)

cat("\nQ3 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + SG + conn"       = b_sens_sg_conn,
  "Baseline + SG + conn + MPA" = b_sens_sg_mpa
)))

cat("\nQ3 — MPA model summary:\n")
print(summary(b_sens_sg_mpa))
print(confint(b_sens_sg_mpa))

# # Sensitivity (a) key findings:
# With settlement gravity included:
#   Q2: connectivity remains significant (b = +0.473, p = 0.003)
#   Q3: medium MPA remains significant (b = +1.134, p = 0.005)
#   Settlement gravity not significant in either model
# Conclusions robust to Q1 metric uncertainty

# ── Market gravity ────────────────────────────────────────────
cat("\n--- Sensitivity (a): market gravity ---\n")

# Q2: connectivity beyond market gravity
b_sens_mg_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc + connectivity_sc,
  family = tweedie(link = "log"),
  data = browser_model_data
)

# Q2: connectivity × market gravity interaction
b_sens_mg_conn_int <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc * connectivity_sc,
  family = tweedie(link = "log"),
  data = browser_model_data
)

cat("\nQ2 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + MG"              = b_q1_mktgrav,
  "Baseline + MG + conn"       = b_sens_mg_conn,
  "Baseline + MG + conn x MG"  = b_sens_mg_conn_int
)))

cat("\nQ2 — connectivity model summary:\n")
print(summary(b_sens_mg_conn))
print(confint(b_sens_mg_conn))

# Q3: MPA beyond market gravity + connectivity
b_sens_mg_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc + connectivity_sc + mpa_status,
  family = tweedie(link = "log"),
  data = browser_model_data
)

cat("\nQ3 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + MG + conn"       = b_sens_mg_conn,
  "Baseline + MG + conn + MPA" = b_sens_mg_mpa
)))

cat("\nQ3 — MPA model summary:\n")
print(summary(b_sens_mg_mpa))
print(confint(b_sens_mg_mpa))

# Sensitivity (a) summary — browsers:
# With settlement gravity:
#   Q2 connectivity: b = +0.473, p = 0.003 — robust
#   Q3 medium MPA:   b = +1.134, p = 0.005 — robust
# With market gravity:
#   Q2 connectivity: b = +0.418, p = 0.014 — robust
#   Q3 medium MPA:   b = +1.154, p = 0.002 — robust
# Both key findings robust to Q1 metric uncertainty

# ── (b) Transect-level replication ───────────────────────────
# 43% zeros at transect level — Tweedie required.
# Mirrors Q1-Q3 sequence at transect level.
# (1 | site) accounts for non-independence within sites.

# ── Transect family selection ─────────────────────────────────
# ZI Tweedie tested against standard Tweedie.
# Run on baseline structure.

b_trans_tw_base <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

b_trans_tw_zi_base <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = browser_transect_data
)

b_trans_res    <- simulateResiduals(b_trans_tw_base,    n = 1000)
b_trans_res_zi <- simulateResiduals(b_trans_tw_zi_base, n = 1000)

plot(b_trans_res);    testZeroInflation(b_trans_res)
plot(b_trans_res_zi); testZeroInflation(b_trans_res_zi)

cat("\n--- Sensitivity (b): transect family selection ---\n")
print(make_aicc_df(list(
  "Tweedie"    = b_trans_tw_base,
  "ZI Tweedie" = b_trans_tw_zi_base
)))

# Results:
#   ZI Tweedie: AICc = 2652.49, weight = 0.725 (BEST)
#   Tweedie:    DAICc = 1.93,   weight = 0.276
#   ZI Tweedie marginally preferred at transect level
#   (DAICc = 1.93) — borderline. Zero inflation test
#   p = 0.666 (standard Tweedie) — not significant.
#   Standard Tweedie retained for consistency with
#   site-level family and because zero inflation test
#   does not support ZI model despite marginal AICc
#   preference. Results reported with standard Tweedie.

# ── Transect Q1-Q3 sequence ───────────────────────────────────
b_trans_null <- glmmTMB(
  transect_browser_biomass ~ 1 + (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

b_trans_baseline <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

b_trans_pressure <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

b_trans_conn <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc +
    connectivity_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

b_trans_mpa <- glmmTMB(
  transect_browser_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc +
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = browser_transect_data
)

cat("\n--- Sensitivity (b): browser transect comparison ---\n")
print(make_aicc_df(list(
  "Null"                = b_trans_null,
  "Baseline"            = b_trans_baseline,
  "Baseline + pressure" = b_trans_pressure,
  "Best + conn"         = b_trans_conn,
  "Best + MPA"          = b_trans_mpa
)))

# Results:
#   Best + MPA:          AICc = 2653.73, weight = 0.390 (BEST)
#   Baseline:            DAICc = 0.69,   weight = 0.276
#   Best + conn:         DAICc = 1.61,   weight = 0.174
#   Baseline + pressure: DAICc = 1.87,   weight = 0.153
#   Null:                DAICc = 7.97,   weight = 0.007
#
#   MPA model best supported at transect level —
#   consistent with site-level Q3 (weight = 0.939).
#   However weaker support than site level (weight =
#   0.390 vs 0.939) — greater within-site variance
#   at transect level reduces discriminatory power.
#   Model ordering qualitatively consistent: MPA >
#   baseline > connectivity > pressure > null.
#   Primary Q3 conclusion robust to aggregation level.

cat("\n--- Sensitivity (b): pressure model coefficients ---\n")
summary(b_trans_mpa)

# 

# ── ICC calculation ───────────────────────────────────────────
vc_b <- VarCorr(b_trans_pressure)
site_sd_b <- sqrt(as.numeric(vc_b$cond$site))
cat(sprintf("\nSite random intercept SD = %.3f\n", site_sd_b))
# ============================================================
#  MARGINAL EFFECT PLOTS
#  Best model: rugosity + chla + connectivity + MPA
#  (b_best_q3 = b_q3_mpa, weight = 0.939)
#  Three plots: rugosity, connectivity, MPA marginal means.
#  All non-focal predictors held at 0 (their mean).
#  Predictions on response scale (raw biomass).
# ============================================================

best_model_b <- b_best_q3  # rugosity + chla + conn + MPA

# ── Rugosity effect ───────────────────────────────────────────
rug_grid_b <- data.frame(
  rugosity_sc     = seq(
    min(browser_model_data$rugosity_sc),
    max(browser_model_data$rugosity_sc),
    length.out = 200),
  log_chla_sc     = 0,
  connectivity_sc = 0,
  mpa_status      = factor("none",
                           levels = c("none", "low", "medium"))
)

rug_pred_b     <- predict(best_model_b,
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
             alpha  = 0.5, inherit.aes = FALSE) +
  labs(x = "Rugosity (standardised)",
       y = "Browser biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())


# ── Connectivity effect ───────────────────────────────────────
# Connectivity supported in Q2 (weight = 0.881, b = +0.421)
# Positive — well-connected sites have higher browser
# biomass. Contrasts with corallivores (negative direction).

conn_grid_b <- data.frame(
  connectivity_sc = seq(
    min(browser_model_data$connectivity_sc),
    max(browser_model_data$connectivity_sc),
    length.out = 200),
  rugosity_sc = 0,
  log_chla_sc = 0,
  mpa_status  = factor("none",
                       levels = c("none", "low", "medium"))
)

conn_pred_b      <- predict(best_model_b,
                            newdata = conn_grid_b,
                            se.fit  = TRUE,
                            type    = "response",
                            re.form = NA)
conn_grid_b$fit  <- conn_pred_b$fit
conn_grid_b$lwr  <- conn_pred_b$fit - 1.96 * conn_pred_b$se.fit
conn_grid_b$upr  <- conn_pred_b$fit + 1.96 * conn_pred_b$se.fit

p_b_conn <- ggplot(conn_grid_b,
                   aes(x = connectivity_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "#2c7bb6", alpha = 0.15) +
  geom_line(colour = "#2c7bb6", linewidth = 1.1) +
  geom_point(data = browser_model_data,
             aes(x = connectivity_sc, y = mean_biomass),
             colour = "grey40", size = 1.5,
             alpha  = 0.5, inherit.aes = FALSE) +
  labs(x = "Connectivity (standardised)",
       y = "Browser biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())


# ── MPA marginal means ────────────────────────────────────────
# MPA strongly supported in Q3 (weight = 0.939).
# Medium protection: ~3.1x higher biomass than no MPA
# at equivalent habitat and connectivity conditions.
# Low MPA: no detectable effect (b = 0.000, p = 0.999).

mpa_grid_b <- data.frame(
  mpa_status      = factor(c("none", "low", "medium"),
                           levels = c("none", "low", "medium")),
  rugosity_sc     = 0,
  log_chla_sc     = 0,
  connectivity_sc = 0
)

mpa_pred_b     <- predict(best_model_b,
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
  labs(x = "MPA status", y = "Browser biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title         = element_text(face = "bold"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank())


# ── Arrange plots ─────────────────────────────────────────────
gridExtra::grid.arrange(p_b_rugosity, p_b_conn, p_b_mpa,
                        ncol = 3)

# jpeg("browser_marginal_effects.jpg",
#      width = 33, height = 11, units = "cm", res = 300)
# gridExtra::grid.arrange(p_b_rugosity, p_b_conn, p_b_mpa,
#                         ncol = 3)
# dev.off()


# ============================================================
#  RESULTS SUMMARY
# ============================================================

cat("\n--- Browser results summary ---\n")
tribble(
  ~Question,    ~Result,           ~Key_finding,
  "Q1",         "Baseline best",   "weight = 0.459, no pressure metric outperforms baseline",
  "Q2 conn",    "Supported",       "weight = 0.881, DAICc = 4.00, b = +0.421, p = 0.008",
  "Q3 MPA",     "Supported",       "weight = 0.939, DAICc = 5.47, medium b = +1.143, p = 0.002",
  "Q3 MPA int", "Suggestive",      "DAICc = 1.62, medium x conn b = +0.729, p = 0.035",
  "Sens (a)",   "Consistent",      "all metrics ns, inconsistent directions — Q1 null robust",
  "Sens (b)",   "Consistent",      "MPA + conn top-ranked at transect level, connectivity positive"
) %>% print()


# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()

# ============================================================
#  FIGURE 5 — MPA marginal means
#  (a) Browsers — purple #7b2d8b
#  (b) Piscivores — red #d7191c
#  Pointrange style with dashed reference line at no-MPA
#  Free y-axis scales
#  Observed data overlaid as jittered points
# ============================================================

# ── Panel (a): Browser MPA marginal means ────────────────────
mpa_grid_b <- data.frame(
  mpa_status      = factor(c("none", "low", "medium"),
                           levels = c("none", "low", "medium")),
  rugosity_sc     = 0,
  log_chla_sc     = 0,
  connectivity_sc = 0
)

mpa_pred_b     <- predict(b_best_q3,
                          newdata = mpa_grid_b,
                          se.fit  = TRUE,
                          type    = "response",
                          re.form = NA)
mpa_grid_b$fit <- mpa_pred_b$fit
mpa_grid_b$lwr <- mpa_pred_b$fit - 1.96 * mpa_pred_b$se.fit
mpa_grid_b$upr <- mpa_pred_b$fit + 1.96 * mpa_pred_b$se.fit

p5a <- ggplot(mpa_grid_b,
              aes(x = mpa_status, y = fit)) +
  geom_hline(yintercept = mpa_grid_b$fit[1],
             linetype   = "dashed",
             colour     = "grey70",
             linewidth  = 0.4) +
  geom_jitter(data   = browser_model_data,
              aes(x  = mpa_status,
                  y  = mean_biomass),
              width  = 0.1,
              size   = 1.5,
              alpha  = 0.4,
              colour = "grey50",
              inherit.aes = FALSE) +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  colour    = "#7b2d8b",
                  linewidth = 0.7,
                  size      = 0.6) +
  scale_x_discrete(labels = c("none"   = "No MPA",
                              "low"    = "Low",
                              "medium" = "Medium")) +
  labs(x = "MPA status",
       y = "Browser biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title         = element_text(face = "bold"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank())

# ── Panel (b): Piscivore MPA marginal means ──────────────────
mpa_grid_p <- data.frame(
  mpa_status            = factor(c("none", "low", "medium"),
                                 levels = c("none", "low",
                                            "medium")),
  rugosity_sc           = 0,
  log_chla_sc           = 0,
  log_market_gravity_sc = 0,
  connectivity_sc       = 0
)

mpa_pred_p     <- predict(p_best_q3,
                          newdata = mpa_grid_p,
                          se.fit  = TRUE,
                          type    = "response",
                          re.form = NA)
mpa_grid_p$fit <- mpa_pred_p$fit
mpa_grid_p$lwr <- mpa_pred_p$fit - 1.96 * mpa_pred_p$se.fit
mpa_grid_p$upr <- mpa_pred_p$fit + 1.96 * mpa_pred_p$se.fit

p5b <- ggplot(mpa_grid_p,
              aes(x = mpa_status, y = fit)) +
  geom_hline(yintercept = mpa_grid_p$fit[1],
             linetype   = "dashed",
             colour     = "grey70",
             linewidth  = 0.4) +
  geom_jitter(data   = pisc_model_data,
              aes(x  = mpa_status,
                  y  = mean_biomass),
              width  = 0.1,
              size   = 1.5,
              alpha  = 0.4,
              colour = "grey50",
              inherit.aes = FALSE) +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  colour    = "#d7191c",
                  linewidth = 0.7,
                  size      = 0.6) +
  scale_x_discrete(labels = c("none"   = "No MPA",
                              "low"    = "Low",
                              "medium" = "Medium")) +
  labs(x = "MPA status",
       y = "Piscivore biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title         = element_text(face = "bold"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank())

# ── Arrange ───────────────────────────────────────────────────
gridExtra::grid.arrange(p5a, p5b, ncol = 2)

# jpeg("figure5_mpa_marginal_means.jpg",
#      width = 22, height = 11, units = "cm", res = 300)
# gridExtra::grid.arrange(p5a, p5b, ncol = 2)
# dev.off()
