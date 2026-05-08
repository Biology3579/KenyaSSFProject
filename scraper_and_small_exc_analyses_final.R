# ============================================================
#  DRIVERS OF SCRAPER/SMALL EXCAVATOR BIOMASS
#  Functional Group Analysis: Scrapers/Small Excavators
# ============================================================
#
#  Scientific questions:
#
#  Q1 — Human pressure
#       Does human pressure explain variation in scraper
#       biomass beyond local ecological context, and which
#       spatial metric best captures SSF exploitation
#       intensity for this functional group?
#       Scrapers are subsistence-targeted — settlement
#       gravity expected to outperform market gravity.
#
#  Q2 — Larval connectivity
#       Does larval connectivity explain additional variation
#       beyond the best Q1 model?
#       Connectivity × pressure interaction only evaluated
#       if connectivity main effect supported in Q2.
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
#       ~2% zeros at site level — Tweedie selected despite
#       near-zero proportion due to extreme influential
#       observation (itsan site) causing Gaussian log
#       instability.
#       ~8.5% zeros at transect level — standard Tweedie.
#
#  Sensitivity analyses:
#       (a) Alternative pressure metrics — only metrics
#           within DAICc < 2 of best Q1 model evaluated
#       (b) Transect-level replication (Tweedie GLMM)
# ============================================================

source(here::here("predictor_preparation.R"))


# ============================================================
#  RESPONSE VARIABLE — DATA AGGREGATION AND CHECKS
# ============================================================

# ── Transect-level aggregation ────────────────────────────────
scraper_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_scraper_biomass = sum(
      ifelse(trophic_group %in% c("scrapers", "small_excavators"),
             tot_wt_g, 0),
      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(site = as.factor(site))

cat("Scraper transects:", nrow(scraper_transects), "\n")
cat("Sites:",             n_distinct(scraper_transects$site), "\n")

# ── Transect-level dataset ────────────────────────────────────
scraper_transect_data <- scraper_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(scraper_transect_data$transect_scraper_biomass == 0),
    "/", nrow(scraper_transect_data),
    sprintf("(%.3f)\n",
            mean(scraper_transect_data$transect_scraper_biomass == 0)))

# ── Site-level dataset ────────────────────────────────────────
scraper_model_data <- scraper_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    mean_biomass           = mean(transect_scraper_biomass,
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

cat("\nScraper model data:", nrow(scraper_model_data), "sites\n")

# ── Data checks ───────────────────────────────────────────────
scraper_model_data %>%
  dplyr::select(site, rugosity_sc, log_chla_sc,
                log_settlement_grav_sc, connectivity_sc,
                mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  { if (nrow(.) > 0) { warning("NAs in predictors:"); print(.) }
    else cat("NA check passed.\n") }

cat("\nZeros in mean_biomass:",
    sum(scraper_model_data$mean_biomass == 0), "\n")
cat("Site-level zero proportion:",
    round(mean(scraper_model_data$mean_biomass == 0), 3), "\n")
cat("\nResponse summary:\n")
print(summary(scraper_model_data$mean_biomass))
cat("\nMPA status counts:\n")
print(table(scraper_model_data$mpa_status))


# ============================================================
#  MODEL FAMILY SELECTION
#
#  ~2% zeros at site level but Tweedie selected over
#  Gaussian log — extreme influential observation (itsan,
#  Comoros) causes Gaussian log instability.
# ============================================================

scr_lm_baseline <- lm(
  log(mean_biomass + 0.01) ~ rugosity_sc + log_chla_sc,
  data = scraper_model_data)

par(mfrow = c(2, 2))
plot(scr_lm_baseline, main = "Gaussian log — baseline")
par(mfrow = c(1, 1))

scr_tw_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data)

scr_tw_res <- simulateResiduals(scr_tw_baseline, n = 1000)
plot(scr_tw_res)
testZeroInflation(scr_tw_res)
testDispersion(scr_tw_res)

# Family selection:
#   Gaussian log: REJECTED
#     Site 12 (itsan, Comoros) — standardised residual
#     = -4.5, exceeds Cook's distance threshold. Severe
#     lower tail Q-Q deviation. Scale-location shows
#     elevated sqrt(standardised residual) at low fitted
#     values. No offset constant resolves this — itsan
#     is a genuine biological observation (extreme low
#     biomass, ~130x below median) that Gaussian log
#     cannot accommodate without coefficient distortion.
#
#   Tweedie (log link): SELECTED
#     Handles extreme low values natively.
#     DHARMa diagnostics (n = 1000):
#       KS test:        p = 0.589 — good fit
#       Dispersion:     p = 0.240, ratio = 1.359 — acceptable
#       Zero inflation: p = 1.000 — not needed
#       Outlier test:   p = 1.000 — no outliers
#     Minor lower quantile deviation in residuals vs
#     predicted — formal tests all non-significant.
#
#   Proceed: glmmTMB(family = tweedie(link = "log")) on
#   raw mean_biomass throughout all scraper analyses.


# ============================================================
#  RANDOM EFFECT STRUCTURE
#
#  Tested on baseline model — same rationale as all
#  previous analyses. Tweedie family throughout.
# ============================================================

scr_re_null <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data)

scr_re_ecoregion <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + (1 | ecoregion),
  family = tweedie(link = "log"),
  data   = scraper_model_data)

cat("\n--- Random effect structure: ecoregion ---\n")
print(make_aicc_df(list(
  "No RE"           = scr_re_null,
  "(1 | ecoregion)" = scr_re_ecoregion
)))

# Ecoregion RE not supported (DAICc = 2.10, weight = 0.259)
# — consistent with all previous functional groups.
# All scraper models fitted without RE throughout.


# ============================================================
#  Q1 — HUMAN PRESSURE
#
#  Scientific question:
#  Does human pressure explain variation in scraper biomass
#  beyond local ecological context, and which spatial metric
#  best captures SSF exploitation intensity?
#
#  Approach: AICc comparison of baseline vs baseline + each
#  pressure metric. Best metric = highest AICc weight AND
#  outperforms baseline. Coefficient directions checked
#  across all metrics regardless of support.
# ============================================================

s_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data)

s_q1_settgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data)

s_q1_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data)

s_q1_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_pop_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data)

cat("\n--- Q1: Pressure metric comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                      = s_baseline,
  "Baseline + settlement gravity" = s_q1_settgrav,
  "Baseline + market gravity"     = s_q1_mktgrav,
  "Baseline + settlement pop."    = s_q1_settpop
)))

# Results:
#   Settlement gravity: AICc = 1007.32, weight = 0.562 (BEST)
#   Market gravity:     DAICc = [update],  weight = [update]
#   Baseline:           DAICc = 2.85,  weight = [update]
#   Settlement pop.:    DAICc = [update],  weight = [update]
#
#   Settlement gravity best supported — ecologically
#   motivated as scrapers are subsistence-targeted.
#   Market gravity within DAICc < 2 threshold —
#   evaluated in sensitivity (a).
#   [Update all values after running]

# ── Best Q1 model ─────────────────────────────────────────────
# Settlement gravity best supported — carried forward.
s_best_q1 <- s_q1_settgrav

# ── Best model summary — reported in results ─────────────────
cat("\n--- Best model: full summary ---\n")
summary(s_best_q1)

cat("\n--- Best model: 95% CIs ---\n")
print(confint(s_best_q1))

# ── Effect sizes: fold differences across observed range ──────
b_sg_s   <- fixef(s_best_q1)$cond["log_settlement_grav_sc"]
b_chla_s <- fixef(s_best_q1)$cond["log_chla_sc"]
b_rug_s  <- fixef(s_best_q1)$cond["rugosity_sc"]

sg_span_s   <- diff(range(scraper_model_data$log_settlement_grav_sc,
                          na.rm = TRUE))
chla_span_s <- diff(range(scraper_model_data$log_chla_sc,
                          na.rm = TRUE))
rug_span_s  <- diff(range(scraper_model_data$rugosity_sc,
                          na.rm = TRUE))

cat(sprintf("\nSettlement gravity span: %.3f SD units\n", sg_span_s))
cat(sprintf("Fold difference (low vs high pressure): %.2fx\n",
            exp(abs(b_sg_s * sg_span_s))))

cat(sprintf("\nChla span: %.3f SD units\n", chla_span_s))
cat(sprintf("Fold difference (low vs high chla): %.2fx\n",
            exp(abs(b_chla_s * chla_span_s))))

cat(sprintf("\nRugosity span: %.3f SD units\n", rug_span_s))
cat(sprintf("Fold difference (low vs high rugosity): %.2fx\n",
            exp(abs(b_rug_s * rug_span_s))))

# ── CI-based fold differences ─────────────────────────────────
ci_sg_s <- confint(s_best_q1)["log_settlement_grav_sc",
                              c("2.5 %", "97.5 %")]
cat(sprintf("Settlement gravity fold diff 95%% CI: %.2fx to %.2fx\n",
            exp(abs(ci_sg_s[1] * sg_span_s)),
            exp(abs(ci_sg_s[2] * sg_span_s))))

# ── DHARMa diagnostics — best Q1/Q2 model ────────────────────
cat("\n--- DHARMa diagnostics: best Q1 model ---\n")
s_best_q1_sim <- simulateResiduals(s_best_q1, n = 1000)
plot(s_best_q1_sim)
testOutliers(s_best_q1_sim)

# Results:
#   Settlement gravity: b = -0.261, z = [update], p = 0.017 *
#     95% CI: [update]
#     Fold difference: [update]x across observed range
#     (span = [update] SD units)
#     95% CI on fold difference: [update]x to [update]x
#   Chla: b = +0.204, p = 0.037 *
#     Fold difference: [update]x
#   Rugosity: b = +0.143, p = 0.126 ns
#   DHARMa: [update from plots]


# ============================================================
#  Q2 — LARVAL CONNECTIVITY
#
#  Scientific question:
#  Does larval connectivity explain additional variation
#  in scraper biomass beyond the pressure baseline?
#
#  Pressure supported in Q1 — connectivity × pressure
#  interaction evaluated in step 2 if main effect supported.
#
#  Criterion: AICc weight for model support.
# ============================================================

s_q2_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data)

cat("\n--- Q2: Connectivity comparison ---\n")
print(make_aicc_df(list(
  "Best Q1"                = s_best_q1,
  "Best Q1 + connectivity" = s_q2_conn
)))

# ── Best Q2 model ─────────────────────────────────────────────
# Connectivity not supported — best Q1 retained.
s_best_q2 <- s_best_q1

# Results:
#   Best Q1:                AICc = 1007.32, weight = 0.655
#   Best Q1 + connectivity: DAICc = 1.28,   weight = 0.345
#   Connectivity not supported.
#
#   Connectivity: b = -0.105, p = 0.238 ns
#     Negative direction — consistent with corallivores
#     (b = -0.221, p = 0.015) and grazer-detritivores
#     (b = -0.196, p = 0.058). Not significant.
#   Settlement gravity: b = -0.252, p = 0.020 * — stable
#   Chla: b = +0.171, p = 0.087 . — weakens slightly
#   Connectivity × pressure interaction not tested —
#   main effect not supported.


# ============================================================
#  Q3 — FORMAL PROTECTION
#
#  Scientific question:
#  Does MPA status explain additional variation in scraper
#  biomass beyond the best Q2 model?
#
#  Q3 is a separate question — does not update best model.
#  MPA tested last — non-randomly placed with respect to
#  pressure and connectivity (see mpa_placement_checks.R).
#
#  MPA × pressure interaction tested if MPA supported —
#  both components in model.
#
#  Criterion: AICc weight for model support.
# ============================================================

s_q3_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = scraper_model_data)

s_q3_mpa_press_int <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    mpa_status * log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data)

cat("\n--- Q3: MPA comparison ---\n")
print(make_aicc_df(list(
  "Best Q2"                  = s_best_q2,
  "Best Q2 + MPA"            = s_q3_mpa,
  "Best Q2 + MPA x pressure" = s_q3_mpa_press_int
)))

cat("\n--- Q3: MPA coefficients ---\n")
summary(s_q3_mpa)

cat("\n--- Q3: MPA 95% CIs ---\n")
print(confint(s_q3_mpa))

# ── Raw biomass by MPA — sanity check ────────────────────────
cat("\n--- Q3: Raw biomass by MPA status ---\n")
scraper_model_data %>%
  group_by(mpa_status) %>%
  summarise(
    n              = n(),
    mean_biomass   = round(mean(mean_biomass), 1),
    median_biomass = round(median(mean_biomass), 1),
    sd_biomass     = round(sd(mean_biomass), 1),
    mean_conn      = round(mean(connectivity_sc), 3),
    .groups = "drop"
  ) %>% print()

# ── DHARMa diagnostics — Q3 model ────────────────────────────
cat("\n--- DHARMa diagnostics: Q3 model ---\n")
s_q3_sim <- simulateResiduals(s_q3_mpa, n = 1000)
plot(s_q3_sim)
testOutliers(s_q3_sim)

# Results:
#   Best Q2:               AICc = 1007.32, weight = [update]
#   Best Q2 + MPA:         DAICc =  0.89,  weight = [update]
#   Best Q2 + MPA x press: DAICc =  5.58,  weight = [update]
#   MPA not clearly supported — genuine uncertainty
#   (DAICc < 2) but best Q2 preferred.
#   MPA × pressure interaction not supported (DAICc = 5.58).
#
#   Medium MPA: b = -0.422, p = 0.034 * — significant but
#     NEGATIVE and artefactual. Raw means confirm:
#     none = 5267g, low = 4192g, medium = 3975g —
#     unprotected sites have highest mean biomass.
#     MPA placement bias: medium MPA sites at moderately
#     high connectivity (mean z = +0.303), likely in
#     lower-productivity or higher-pressure areas.
#   Low MPA: b = +0.034, p = 0.897 ns
#   Settlement gravity: b = -0.363, p = 0.002 ** —
#     strengthens when MPA included.
#
#   Conclusion: MPA not supported. Best Q2 model retained.
#   DHARMa: [update from plots]

# ── Best Q3 model ─────────────────────────────────────────────
# Q3 is a separate question — best Q2 model unchanged.
s_best_q3 <- s_best_q2

# ── Predicted vs observed ─────────────────────────────────────
pred_s <- predict(s_best_q3, type = "response")
cat(sprintf("\nPearson r (predicted vs observed): %.3f\n",
            cor(pred_s, scraper_model_data$mean_biomass)))

# ── Flagged sites ─────────────────────────────────────────────
# Update row indices from diagnostic plots after running
flagged_rows_s <- c(12)  # itsan (Comoros) — extreme low biomass
# confirmed in family selection; Tweedie handles adequately

scraper_model_data %>%
  slice(flagged_rows_s) %>%
  dplyr::select(site, mpa_status, ecoregion,
                mean_biomass, rugosity_sc,
                log_settlement_grav_sc, connectivity_sc) %>%
  print()


# ============================================================
#  SPATIAL AUTOCORRELATION CHECK
#  Pearson residuals from best Q3 model.
# ============================================================

site_coords <- location_data %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(lon = first(longitude),
            lat = first(latitude),
            .groups = "drop")

scraper_model_data_coords <- scraper_model_data %>%
  left_join(site_coords, by = "site")

coords_mat_scr <- cbind(scraper_model_data_coords$lon,
                        scraper_model_data_coords$lat)
listw5_scr <- nb2listw(knn2nb(knearneigh(coords_mat_scr, k = 5)),
                       style = "W")

# Warning: knearneigh — identical points found, kd_tree not
# available. 4 sub-graphs expected — discontinuous sampling
# design across four countries. See total_biomass.R.

cat("\n--- Spatial autocorrelation: scraper best model ---\n")
print(moran.test(residuals(s_best_q3, type = "pearson"),
                 listw5_scr))

# Moran's I = +0.033, p = 0.230 — no significant spatial
# autocorrelation in residuals.
# Positive but non-significant — similar to piscivores
# (I = +0.052, p = 0.162). Settlement gravity and chla
# adequately capture spatial variation in scraper biomass.
# No spatial error modelling required.


# ============================================================
#  SENSITIVITY ANALYSES
#
#  (a) Alternative pressure metrics
#      Only metrics within DAICc < 2 of best Q1 model
#      are evaluated. Settlement gravity best supported.
#      Market gravity within threshold — evaluated.
#      Settlement pop. outside threshold — not evaluated.
#
#  (b) Transect-level replication (Tweedie GLMM)
# ============================================================

# ── (a) Market gravity ────────────────────────────────────────
cat("\n--- Sensitivity (a): market gravity ---\n")

s_sens_mg_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc + connectivity_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data)

cat("\nQ2 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + MG"        = s_q1_mktgrav,
  "Baseline + MG + conn" = s_sens_mg_conn
)))

cat("\nQ2 — connectivity coefficients:\n")
print(summary(s_sens_mg_conn)$coefficients$cond)
print(confint(s_sens_mg_conn))

s_sens_mg_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc +
    log_market_gravity_sc + connectivity_sc + mpa_status,
  family = tweedie(link = "log"),
  data   = scraper_model_data)

cat("\nQ3 — model comparison:\n")
print(make_aicc_df(list(
  "Baseline + MG + conn"       = s_sens_mg_conn,
  "Baseline + MG + conn + MPA" = s_sens_mg_mpa
)))

cat("\nQ3 — MPA coefficients:\n")
print(summary(s_sens_mg_mpa)$coefficients$cond)
print(confint(s_sens_mg_mpa))

# Results: [update after running]
# Expected: market gravity negative and significant
# (b = -0.198, p = 0.049 from Q1 direction check);
# connectivity and MPA conclusions consistent.


# ── (b) Transect-level replication ───────────────────────────
# ~8.5% zeros at transect level — ZI Tweedie tested.

scr_trans_tw_base <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc + log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data)

scr_trans_tw_zi_base <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc + log_chla_sc +
    (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = scraper_transect_data)

scr_trans_res    <- simulateResiduals(scr_trans_tw_base,    n = 500)
scr_trans_res_zi <- simulateResiduals(scr_trans_tw_zi_base, n = 500)

plot(scr_trans_res);    testZeroInflation(scr_trans_res)
plot(scr_trans_res_zi); testZeroInflation(scr_trans_res_zi)

cat("\n--- Sensitivity (b): transect family selection ---\n")
print(make_aicc_df(list(
  "Tweedie"    = scr_trans_tw_base,
  "ZI Tweedie" = scr_trans_tw_zi_base
)))

# Results:
#   Standard Tweedie: AICc = 4335.13, weight = 0.862
#   ZI Tweedie:       DAICc = 3.66,   weight = 0.138
#   Standard Tweedie clearly selected.

# ── Transect Q1-Q3 sequence ───────────────────────────────────
scr_trans_null <- glmmTMB(
  transect_scraper_biomass ~ 1 + (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data)

scr_trans_baseline <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc + log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data)

scr_trans_pressure <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data)

scr_trans_conn <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc +
    connectivity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data)

scr_trans_mpa <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc + log_chla_sc +
    log_settlement_grav_sc +
    mpa_status + (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data)

cat("\n--- Sensitivity (b): transect model comparison ---\n")
print(make_aicc_df(list(
  "Null"                = scr_trans_null,
  "Baseline"            = scr_trans_baseline,
  "Baseline + pressure" = scr_trans_pressure,
  "Baseline + conn"     = scr_trans_conn,
  "Best + MPA"          = scr_trans_mpa
)))

cat("\n--- Sensitivity (b): pressure model coefficients ---\n")
summary(scr_trans_pressure)

# ── Site random intercept SD ──────────────────────────────────
vc_s      <- VarCorr(scr_trans_pressure)
site_sd_s <- sqrt(as.numeric(vc_s$cond$site))
cat(sprintf("\nSite random intercept SD = %.3f\n", site_sd_s))

# Results:
#   Baseline + press: AICc = 4330.61, weight = 0.452
#   Baseline + conn:  DAICc = 0.84,   weight = 0.296
#   Best + MPA:       DAICc = 1.58,   weight = 0.205
#   Baseline:         DAICc = 4.53,   weight = 0.047
#   Null:             DAICc = 15.61,  weight = 0.000
#
#   Pressure model best supported — consistent with
#   site-level Q1. Null decisively rejected (DAICc =
#   15.61) — strongest transect-level replication
#   across all functional groups.
#
#   Settlement gravity: b = -0.283, p = 0.009 **
#     Consistent with site-level (b = -0.261, p = 0.017).
#   Rugosity: b = +0.178, p = 0.042 * — significant at
#     transect level (marginal at site level, p = 0.126).
#   Chla: b = +0.195, p = 0.054 . — consistent direction.
#   Site random intercept SD = [update]


# ============================================================
#  RESULTS SUMMARY
# ============================================================

cat("\n--- Scraper results summary ---\n")
tribble(
  ~Question,   ~Result,          ~Key_finding,
  "Q1",        "Sett. gravity",  "weight = 0.562, DAICc = 2.85 vs baseline; b = -0.261, p = 0.017; fold diff = [update]x",
  "Q1 chla",   "Significant",    "b = +0.204, p = 0.037; fold diff = [update]x — most robust baseline predictor",
  "Q2 conn",   "Not supported",  "DAICc = 1.28, weight = 0.345; b = -0.105, p = 0.238 ns",
  "Q3 MPA",    "Not supported",  "DAICc = 0.89, weight = [update] — medium MPA negative, artefactual; best Q2 retained",
  "Q3 int",    "Not supported",  "MPA x pressure DAICc = 5.58",
  "Moran",     "No SAC",         "I = +0.033, p = 0.230",
  "Pearson r", "[update]",       "[update after running]",
  "Sens (a)",  "[update]",       "MG within threshold (DAICc < 2); [update after running]",
  "Sens (b)",  "Consistent",     "pressure best at transect level (weight = 0.452); null DAICc = 15.61; site SD = [update]"
) %>% print()


# ── End of script ─────────────────────────────────────────────

