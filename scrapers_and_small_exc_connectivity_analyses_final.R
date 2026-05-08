# ============================================================
#  DRIVERS OF SCRAPER/SMALL EXCAVATOR BIOMASS
#  Chapter 1 — Functional Group Analysis: Scrapers
#
#  Scrapers and small excavators combined as a single
#  functional guild — both process reef substrate through
#  grazing and minor excavation.
#
#  Analytical framework mirrors all other functional groups.
#  Structure: Q1 (pressure) → Q2 (connectivity) → Q3 (MPA)
#
#  Prior analysis found settlement gravity decisively
#  preferred in Q1 (weight = 0.954) — strongest metric
#  selection result across all functional groups.
#  Confirm with current data.
#
#  Key differences from other groups:
#    ~2% zeros at site level — Tweedie selected despite
#    near-zero proportion due to extreme influential
#    observation (itsan site) causing Gaussian log
#    instability.
#    ~8.5% zeros at transect level — standard Tweedie.
#
#  Sensitivity analyses:
#    (a) Alternative pressure metrics
#    (b) Transect-level replication (Tweedie GLMM)
# ============================================================

source(here::here("data_preparation.R"))

# ============================================================
#  DATA AGGREGATION
# ============================================================

scraper_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_scraper_biomass = sum(
      ifelse(trophic_group %in% c("scrapers", "small_excavators"),
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

cat("Scraper transects:", nrow(scraper_transects), "\n")
cat("Sites:",             n_distinct(scraper_transects$site), "\n")

# ── Site-level dataset ────────────────────────────────────────
scraper_model_data <- scraper_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_scraper_biomass,
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

cat("\nScraper model data:", nrow(scraper_model_data), "sites\n")

# ── Data checks ───────────────────────────────────────────────
scraper_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                log_chla_sc, connectivity_sc, mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

cat("\nZeros in mean_biomass:",
    sum(scraper_model_data$mean_biomass == 0), "\n")
cat("Site-level zero proportion:",
    round(mean(scraper_model_data$mean_biomass == 0), 3), "\n")
cat("\nResponse summary:\n")
print(summary(scraper_model_data$mean_biomass))

cat("\nMPA status counts:\n")
print(table(scraper_model_data$mpa_status))

# ── Transect-level dataset ────────────────────────────────────
scraper_transect_data <- scraper_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(scraper_transect_data$transect_scraper_biomass == 0),
    "/", nrow(scraper_transect_data),
    "(", round(mean(scraper_transect_data$transect_scraper_biomass == 0),
               3), ")\n")


# ============================================================
#  MODEL FAMILY SELECTION
#  Run on baseline model.
#  Despite near-zero site-level zero proportion (~2%),
#  Tweedie expected over Gaussian log due to extreme
#  influential observation (itsan site).
# ============================================================

scr_lm_baseline <- lm(
  log(mean_biomass + 0.01) ~ rugosity_sc +
    log_chla_sc,
  data = scraper_model_data
)

par(mfrow = c(2, 2))
plot(scr_lm_baseline, main = "Gaussian log — baseline")
par(mfrow = c(1, 1))

scr_tw_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

scr_tw_res <- simulateResiduals(scr_tw_baseline, n = 1000)
plot(scr_tw_res)
testZeroInflation(scr_tw_res)
testDispersion(scr_tw_res)

# ── Family selection decision ─────────────────────────────────
#
# Gaussian log: REJECTED
#   Site 12 (itsan, Comoros) — standardised residual
#   = -4.5, exceeds Cook's distance = 1 threshold.
#   Both outlier and highly influential observation.
#   Severe lower tail Q-Q deviation. Scale-location
#   shows elevated sqrt(standardised residual) at low
#   fitted values. No offset constant resolves this —
#   itsan is a genuine biological observation (extreme
#   low biomass, ~130x below median) that Gaussian log
#   cannot accommodate without coefficient distortion.
#
# Tweedie (log link): SELECTED
#   Handles extreme low values natively through
#   variance structure. DHARMa diagnostics (n = 1000):
#     KS test:        p = 0.589 — good fit
#     Dispersion:     p = 0.240, ratio = 1.359 — acceptable
#     Zero inflation: p = 1.000 — not needed
#     Outlier test:   p = 1.000 — no outliers
#   Minor lower quantile deviation in residuals vs
#   predicted — formal tests all non-significant.
#   Substantially better than Gaussian log for this
#   dataset.
#
#  Proceed: glmmTMB(family = tweedie(link = "log")) on
#  raw mean_biomass throughout all scraper analyses.

# ============================================================
#  RANDOM EFFECT STRUCTURE
#  Tested on baseline model. Tweedie family throughout.
# ============================================================

scr_re_null <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

scr_re_ecoregion <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | ecoregion),
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

cat("\n--- Scraper RE structure: ecoregion ---\n")
print(make_aicc_df(list(
  "No RE"           = scr_re_null,
  "(1 | ecoregion)" = scr_re_ecoregion
)))

# Random effect structure: ecoregion
# Tested on baseline model to avoid circularity.
# No RE:           AICc = 1010.17, weight = 0.741 (BEST)
# (1 | ecoregion): DAICc = 2.10,   weight = 0.259
# Ecoregion RE not supported — consistent with all
# previous analyses. Note: prior analysis found country
# RE more marginal (DAICc = 0.91) — ecoregion RE here
# shows clearer non-support (DAICc = 2.10), consistent
# with total biomass (DAICc = 2.25) and corallivores
# (DAICc = 2.54).
# All scraper models fitted without RE throughout.

# ── Variance inflation factors ────────────────────────────────
# Sequence: baseline → baseline + settlement gravity →
#           baseline + settlement gravity + connectivity →
#           baseline + settlement gravity + MPA

cat("\n--- VIF: baseline ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data))

cat("\n--- VIF: baseline + settlement gravity ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data))

cat("\n--- VIF: baseline + settlement gravity + connectivity ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + log_settlement_grav_sc + connectivity_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data))

cat("\n--- VIF: baseline + settlement gravity + MPA ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + log_settlement_grav_sc + mpa_status,
  family = tweedie(link = "log"),
  data   = scraper_model_data))

# ============================================================
#  Q1 — HUMAN PRESSURE
#
#  Does human pressure explain variation in scraper
#  biomass beyond local ecological context?
#  Which metric best captures SSF exploitation intensity?
#
#  Prior analysis found settlement gravity decisively
#  preferred — ecologically motivated as scrapers are
#  subsistence-targeted. Confirm with current data.
# ============================================================

s_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

s_q1_settgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

s_q1_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

s_q1_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_pop_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

cat("\n--- Q1 Step 1: Scraper metric comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                      = s_baseline,
  "Baseline + settlement gravity" = s_q1_settgrav,
  "Baseline + market gravity"     = s_q1_mktgrav,
  "Baseline + settlement pop."    = s_q1_settpop
)))

cat("\n--- Q1 Step 2: Baseline coefficients ---\n")
summary(s_baseline)

# ── Q1: Baseline confint ─────────────────────────────────────
cat("\n--- Scraper Q1: baseline confint ---\n")
print(confint(s_baseline))

cat("\n--- Q1: Pressure metric direction checks ---\n")
cat("Settlement gravity:\n")
summary(s_q1_settgrav)

# ── Q1: Settlement gravity range and fold difference ─────────
cat("\n--- Scraper Q1: settlement gravity range and fold difference ---\n")
sg_range_s <- range(scraper_model_data$log_settlement_grav_sc, na.rm = TRUE)
cat("Range:", sg_range_s, "\n")
sg_span_s  <- diff(sg_range_s)
cat(sprintf("Span: %.3f SD units\n", sg_span_s))
b_sg_s     <- -0.261
fold_sg_s  <- exp(abs(b_sg_s * sg_span_s))
cat(sprintf("Fold difference (low vs high pressure): %.2fx\n", fold_sg_s))

# ── Q1: Settlement gravity confint ────────────────────────────
cat("\n--- Scraper Q1: settlement gravity model confint ---\n")
print(confint(s_q1_settgrav))

# ── Q1: Chla range and fold difference ───────────────────────
cat("\n--- Scraper Q1: chla range and fold difference ---\n")
chla_range_s <- range(scraper_model_data$log_chla_sc, na.rm = TRUE)
chla_span_s  <- diff(chla_range_s)
cat(sprintf("Chla span: %.3f SD units\n", chla_span_s))
b_chla_s     <- 0.204
fold_chla_s  <- exp(abs(b_chla_s * chla_span_s))
cat(sprintf("Fold difference (low vs high chla): %.2fx\n", fold_chla_s))

cat("\nMarket gravity:\n")
print(summary(s_q1_mktgrav)$coefficients$cond[
  "log_market_gravity_sc", ])
cat("\nSettlement population:\n")
print(summary(s_q1_settpop)$coefficients$cond[
  "log_settlement_pop_sc", ])

# Q1 Step 2 results:
#
# Baseline coefficients:
#   Rugosity: b = +0.163, p = 0.088 . marginal
#   Chla:     b = +0.320, p < 0.001 *** SIGNIFICANT
#     Chla strongly significant even before pressure
#     included — the most robust baseline predictor
#     across all functional groups for scrapers.
#     Positive direction — higher productivity supports
#     scraper biomass, consistent with algal-grazing
#     guild benefiting from productive environments.
#   Dispersion = 3.28 — moderate, lower than piscivores
#     and excavators, reflecting less patchy distribution.
#
# Settlement gravity coefficients (best Q1 model):
#   Settlement gravity: b = -0.261, p = 0.017 * SIGNIFICANT
#     Negative — higher proximity-weighted SSF pressure
#     reduces scraper biomass. Expected direction,
#     consistent with scrapers being directly targeted
#     by subsistence fishing. The only functional group
#     besides total biomass where pressure is significant
#     at this stage.
#   Chla:     b = +0.204, p = 0.037 * — weakens slightly
#     but remains significant. Direction consistent.
#   Rugosity: b = +0.143, p = 0.126 ns — positive but
#     not significant once pressure included.
#   Dispersion = 2.88 — improves from baseline (3.28).
#
# Pressure metric direction checks:
#   Market gravity:     b = -0.198, p = 0.049 * significant
#     Consistent direction and significant — pressure
#     signal detectable with market gravity too, though
#     weaker than settlement gravity.
#   Settlement pop.:    b = -0.076, p = 0.430 ns
#     Negative but not significant — population size
#     alone insufficient proxy for SSF pressure.
#
#   All three metrics negative — direction consistent
#   throughout. Settlement gravity gives strongest
#   signal (p = 0.017), consistent with Q1 metric
#   selection result.

# ── Best Q1 model ─────────────────────────────────────────────
s_best_q1 <- s_q1_settgrav  # (rugosity + chla + settlement gravity)

# Given that market grvaity falls within 2 DAICc, we include both models in the paper
# Must recalculate Weights to incldue these 
print(make_aicc_df(list(
  "Settlement gravity" = s_q1_settgrav,
  "Market gravity"     = s_q1_mktgrav
)))


# ============================================================
#  Q2 — LARVAL CONNECTIVITY
#
#  Does connectivity explain additional variation beyond
#  the best Q1 model?
#
#  Pressure supported in Q1 (settlement gravity, weight =
#  0.562) — connectivity x pressure interaction testable
#  in Step 2 if connectivity main effect supported.
#
#  Tweedie family throughout.
# ============================================================

s_q2_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

cat("\n--- Q2: Scraper connectivity main effect ---\n")
print(make_aicc_df(list(
  "Best Q1"        = s_best_q1,
  "Best Q1 + conn" = s_q2_conn
)))

cat("\n--- Q2: Connectivity coefficients ---\n")
summary(s_q2_conn)

# Q2 results:
#   Best Q1:        AICc = 1007.32, weight = 0.655 (BEST)
#   Best Q1 + conn: DAICc = 1.28,   weight = 0.345
#   Connectivity not supported (DAICc = 1.28) — best
#   Q1 model retained.
#
#   Connectivity: b = -0.105, p = 0.238 ns
#     Negative direction — consistent with corallivores
#     (b = -0.221, p = 0.015) and grazer-detritivores
#     (b = -0.196, p = 0.058). Not significant for
#     scrapers.
#   Settlement gravity: b = -0.252, p = 0.020 * — stable
#   Chla:               b = +0.171, p = 0.087 . — weakens
#     slightly but direction consistent.
#   Rugosity:           b = +0.143, p = 0.123 ns — stable
#
#   Connectivity × pressure interaction not tested —
#   main effect not supported and prior analysis found
#   MPA × connectivity confound (VIF > 50) suggesting
#   data structure limitations for interaction models.

# ── Best Q2 model ─────────────────────────────────────────────
# Connectivity not supported — best Q1 retained.
s_best_q2 <- s_best_q1  # rugosity + chla + settlement gravity

# ============================================================
#  Q3 — FORMAL PROTECTION
#
#  Does MPA status explain additional variation beyond
#  the best Q2 model (rugosity + chla + settlement gravity)?
#
#  MPA x pressure interaction testable — both components
#  in model. Tested if MPA supported as main effect.
#
#  MPA x connectivity NOT tested — connectivity not
#  supported as main effect in Q2 (DAICc = 1.28).
#  Interaction requires both components justified
#  as main effects first.
#
#  Tweedie family throughout.
# ============================================================

s_q3_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

cat("\n--- Q3: Scraper MPA main effect ---\n")
print(make_aicc_df(list(
  "Best Q2"       = s_best_q2,
  "Best Q2 + MPA" = s_q3_mpa
)))

cat("\n--- Q3: MPA coefficients ---\n")
summary(s_q3_mpa)

# ── Q3: MPA placement check ──────────────────────────────────
cat("\n--- Scraper Q3: pressure and connectivity by MPA ---\n")
scraper_model_data %>%
  group_by(mpa_status) %>%
  summarise(
    n                    = n(),
    mean_biomass         = round(mean(mean_biomass), 1),
    mean_settlement_grav = round(mean(log_settlement_grav_sc), 3),
    mean_connectivity    = round(mean(connectivity_sc), 3),
    .groups = "drop"
  ) %>%
  print()

# Q3 MPA main effect results:
#   Best Q2:       AICc = 1007.32, weight = 0.610 (BEST)
#   Best Q2 + MPA: DAICc = 0.89,   weight = 0.390
#   MPA not clearly supported — genuine model selection
#   uncertainty but best Q2 preferred.
#
#   Medium MPA: b = -0.422, p = 0.034 * — significant
#     but NEGATIVE. Counterintuitive — medium protection
#     associated with lower scraper biomass than
#     unprotected sites. Check site distribution below.
#   Low MPA:    b = +0.034, p = 0.897 ns — not significant.
#   Settlement gravity: b = -0.363, p = 0.002 **
#     Strengthens when MPA included — MPA absorbs
#     additional governance variance, clarifying the
#     pressure signal.
#   Chla:       b = +0.147, p = 0.127 ns — weakens
#     slightly but direction consistent.

# ── MPA data structure check ──────────────────────────────────
cat("\n--- Q3: MPA site distribution ---\n")
scraper_model_data %>%
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

# ── VIF check on MPA main effect ─────────────────────────────
cat("\n--- Q3: VIF check on MPA model ---\n")
check_collinearity(s_q3_mpa)

# VIF clean (all < 2) — no collinearity concern in
# MPA main effect model. Medium MPA negative coefficient
# reflects site characteristics, not collinearity.

# ── Q3 MPA × pressure interaction ────────────────────────────
# Both pressure and MPA are in the model — interaction
# testable. Does MPA effectiveness depend on fishing
# intensity? Ecologically coherent — protection may be
# more effective at high-pressure sites where scraper
# biomass is most depleted.

s_q3_mpa_press_int <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    mpa_status * log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

cat("\n--- Q3: MPA x pressure interaction ---\n")
print(make_aicc_df(list(
  "Best Q2"                  = s_best_q2,
  "Best Q2 + MPA"            = s_q3_mpa,
  "Best Q2 + MPA x pressure" = s_q3_mpa_press_int
)))

cat("\n--- Q3: MPA x pressure coefficients ---\n")
summary(s_q3_mpa_press_int)

# Q3 results:
#
# MPA main effect:
#   Best Q2:       AICc = 1007.32, weight = 0.588 (BEST)
#   Best Q2 + MPA: DAICc = 0.89,   weight = 0.376
#   MPA not clearly supported — best Q2 preferred but
#   genuine model selection uncertainty (DAICc < 2).
#
#   Medium MPA: b = -0.422, p = 0.034 * — significant
#     but negative. Counterintuitive direction.
#   Low MPA:    b = +0.034, p = 0.897 ns
#   VIF: all < 2 — no collinearity concern.
#
#   MPA site distribution:
#     none:   n = 30, mean = 5267g, conn = -0.345
#     low:    n = 7,  mean = 4192g, conn = +0.879
#     medium: n = 17, mean = 3975g, conn = +0.303
#
#   Raw means confirm negative direction — unprotected
#   sites have highest mean scraper biomass (5267g vs
#   3975g for medium). This likely reflects MPA placement
#   bias: medium MPA sites are at moderately high
#   connectivity (mean z = +0.303) — same pattern seen
#   across other functional groups where MPA placement
#   confounds the protection signal.
#   Settlement gravity strengthens (b = -0.363, p = 0.002)
#   when MPA included — MPA placement in lower-pressure
#   areas absorbs variance that clarifies the pressure
#   signal once controlled.
#
# MPA x pressure interaction:
#   Best Q2 + MPA x press: DAICc = 5.58, weight = 0.036
#   Clearly not supported — interaction terms both ns
#   (low: p = 0.600, medium: p = 0.313). MPA
#   effectiveness does not depend on fishing intensity
#   for scrapers.
#
# Overall Q3 conclusion:
#   MPA not supported for scrapers. Best Q2 model
#   (rugosity + chla + settlement gravity) retained.
#   Negative medium MPA coefficient reflects placement
#   bias — not a genuine negative protection effect.
#   Baseline pressure and productivity drivers are
#   the primary findings.
#
# Best Q3 model: s_best_q2 — unchanged throughout Q3.
s_best_q3 <- s_best_q2  # rugosity + chla + settlement gravity

cat("\n--- Scraper: predicted vs observed ---\n")
pred_s <- predict(s_best_q3, type = "response")
obs_s  <- scraper_model_data$mean_biomass
cat(sprintf("Pearson r: %.3f\n", cor(pred_s, obs_s)))

cat("\n--- Best model: DHARMa diagnostics ---\n")
scr_best_sim <- simulateResiduals(s_best_q3, n = 1000)
plot(scr_best_sim)
testOutliers(scr_best_sim)

# ── Q3 exploratory checks ────────────────────────────────────
cat("\n--- Q3: Biomass by MPA status ---\n")
scraper_model_data %>%
  group_by(mpa_status) %>%
  summarise(
    n            = n(),
    mean_biomass = round(mean(mean_biomass), 1),
    sd_biomass   = round(sd(mean_biomass),   1),
    mean_conn    = round(mean(connectivity_sc), 3),
    .groups = "drop"
  ) %>%
  print()

scr_mpa_sim <- simulateResiduals(s_q3_mpa, n = 1000)
plot(scr_mpa_sim)
testOutliers(scr_mpa_sim)

ggplot(scraper_model_data,
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
  labs(x = "MPA status", y = "Scraper biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title         = element_text(face = "bold"),
        legend.position    = "none",
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank())

# ============================================================
#  SPATIAL AUTOCORRELATION
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
listw5_scr <- nb2listw(knn2nb(knearneigh(coords_mat_scr,
                                         k = 5)),
                       style = "W")

cat("\n--- Spatial autocorrelation: scraper best model ---\n")
print(moran.test(residuals(s_best_q3, type = "pearson"),
                 listw5_scr))

# Spatial autocorrelation: scraper best model
# (rugosity + chla + settlement gravity)
# Moran's I = +0.033, p = 0.230 — no significant
# spatial autocorrelation in residuals.
#
# Positive but non-significant — identical magnitude
# to piscivores (I = +0.034, p = 0.630). No spatial
# error modelling required.
# Settlement gravity and chla adequately capture the
# spatial variation in scraper biomass without leaving
# a residual geographic signal.

# ============================================================
#  SENSITIVITY ANALYSIS
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
s_sens_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

s_sens_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_pop_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

cat("\n--- Sensitivity (a): scraper alternative metrics ---\n")
cat("Market gravity:\n")
print(summary(s_sens_mktgrav)$coefficients$cond)
cat("\nSettlement population:\n")
print(summary(s_sens_settpop)$coefficients$cond)

# Sensitivity (a) results:
#   Market gravity:      b = -0.198, p = 0.049 * significant
#   Settlement pop.:     b = -0.076, p = 0.430 ns
#
#   Market gravity significant and negative — pressure
#   signal detectable with both gravity metrics.
#   Settlement population not significant — population
#   size alone insufficient proxy for SSF pressure
#   on scrapers.
#
#   Chla positive and significant across all metrics:
#     Settlement gravity: b = -0.261*, chla b = +0.204*
#     Market gravity:     b = -0.198*, chla b = +0.283**
#     Settlement pop.:    b = -0.076 ns, chla b = +0.313***
#   Chla is the most robust predictor — significant
#   regardless of pressure metric. Pressure negative
#   in direction throughout — both gravity metrics
#   significant, settlement population not.
#   Consistent with Q1 metric selection: gravity-based
#   framework outperforms population size alone.
#   Settlement gravity retained as primary metric —
#   strongest signal and most ecologically appropriate
#   for subsistence SSF on scrapers.

# ── (b) Transect-level replication ───────────────────────────

# ── Family selection ──────────────────────────────────────────
scr_trans_tw_base <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data
)

scr_trans_tw_zi_base <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = scraper_transect_data
)

scr_trans_res    <- simulateResiduals(scr_trans_tw_base,    n = 500)
scr_trans_res_zi <- simulateResiduals(scr_trans_tw_zi_base, n = 500)

plot(scr_trans_res);    testZeroInflation(scr_trans_res)
plot(scr_trans_res_zi); testZeroInflation(scr_trans_res_zi)

cat("\n--- Sensitivity (b): transect family selection ---\n")
print(make_aicc_df(list(
  "Tweedie"    = scr_trans_tw_base,
  "ZI Tweedie" = scr_trans_tw_zi_base
)))

# Tweedie:    AICc = 4335.13, weight = 0.862 — SELECTED
# ZI Tweedie: DAICc = 3.66,   weight = 0.138
# Standard Tweedie clearly preferred.

# ── Transect models mirroring Q1-Q3 sequence ─────────────────
scr_trans_null <- glmmTMB(
  transect_scraper_biomass ~ 1 + (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data
)

scr_trans_baseline <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data
)

# Mirrors Q1 — pressure baseline
scr_trans_pressure <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data
)

# Mirrors Q2 — connectivity added to pressure baseline
scr_trans_conn <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc +
    connectivity_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data
)

# Mirrors Q3 — MPA added to pressure baseline
scr_trans_mpa <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc +
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data
)

cat("\n--- Sensitivity (b): scraper transect comparison ---\n")
print(make_aicc_df(list(
  "Null"             = scr_trans_null,
  "Baseline"         = scr_trans_baseline,
  "Baseline + press" = scr_trans_pressure,
  "Baseline + conn"  = scr_trans_conn,
  "Best + MPA"       = scr_trans_mpa
)))

cat("\n--- Sensitivity (b): pressure model coefficients ---\n")
summary(scr_trans_pressure)

vc_s <- VarCorr(scr_trans_pressure)
site_sd_s <- sqrt(as.numeric(vc_s$cond$site))
cat(sprintf("\nSite random intercept SD = %.3f\n", site_sd_s))

# Sensitivity (b) results:
#
# Transect family selection:
#   Standard Tweedie: AICc = 4335.13, weight = 0.862
#   ZI Tweedie:       DAICc = 3.66,   weight = 0.138
#   Standard Tweedie clearly selected.
#
# AICc comparison (n = 243 transects, 54 sites):
#   Baseline + press: AICc = 4330.61, weight = 0.452
#   Baseline + conn:  DAICc = 0.84,   weight = 0.296
#   Best + MPA:       DAICc = 1.58,   weight = 0.205
#   Baseline:         DAICc = 4.53,   weight = 0.047
#   Null:             DAICc = 15.61,  weight = 0.000
#
#   Pressure model best supported — consistent with
#   site-level Q1 result. Null decisively rejected
#   (DAICc = 15.61) — stronger rejection than any
#   other functional group at transect level.
#   Model ordering mirrors site level exactly:
#   pressure > connectivity > MPA > baseline > null.
#   This is the strongest transect-level replication
#   across all functional groups.
#
# Pressure model coefficients (n = 243 transects):
#   Settlement gravity: b = -0.283, p = 0.009 **
#     Consistent with site-level (b = -0.261, p = 0.017).
#     Direction and significance stable across scales.
#   Rugosity:           b = +0.178, p = 0.042 *
#     Reaches significance at transect level — greater
#     statistical power (n = 243 vs 54). Direction
#     consistent with site level (b = +0.143, p = 0.126).
#   Chla:               b = +0.195, p = 0.054 .
#     Marginal at transect level — consistent positive
#     direction (site level: b = +0.204, p = 0.037).
#   Site variance: 0.224 (SD = 0.474) — between-site
#     clustering present, (1|site) justified.
#   Dispersion: 65.1 — high within-site variance
#     consistent with scraper aggregation behaviour.
#
# Overall: primary findings fully replicated at
#   transect level. Settlement gravity negative and
#   significant at both scales. Chla positive and
#   consistent. Null decisively rejected — scrapers
#   are the most strongly structured functional group
#   after total biomass.

# ============================================================
#  MARGINAL EFFECT PLOTS
#  Best model: rugosity + chla + settlement gravity
#  (s_best_q3 = s_q1_settgrav)
#  Two plots: settlement gravity (negative) and chla
#  (positive) — the two significant predictors.
#  No connectivity or MPA plots — neither supported.
#  Predictions on response scale (raw biomass).
# ============================================================

best_model_s <- s_best_q3  # rugosity + chla + settlement gravity

# ── Settlement gravity effect ─────────────────────────────────
# Primary result — negative effect of SSF pressure.
# Chla and rugosity held at mean (z = 0).

grav_grid_scr <- data.frame(
  log_settlement_grav_sc = seq(
    min(scraper_model_data$log_settlement_grav_sc),
    max(scraper_model_data$log_settlement_grav_sc),
    length.out = 200),
  rugosity_sc = 0,
  log_chla_sc = 0
)

grav_pred_scr     <- predict(best_model_s,
                             newdata = grav_grid_scr,
                             se.fit  = TRUE,
                             type    = "response",
                             re.form = NA)
grav_grid_scr$fit <- grav_pred_scr$fit
grav_grid_scr$lwr <- grav_pred_scr$fit - 1.96 * grav_pred_scr$se.fit
grav_grid_scr$upr <- grav_pred_scr$fit + 1.96 * grav_pred_scr$se.fit

p_scr_grav <- ggplot(grav_grid_scr,
                     aes(x = log_settlement_grav_sc,
                         y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "#2c7bb6", alpha = 0.15) +
  geom_line(colour = "#2c7bb6", linewidth = 1.1) +
  geom_point(data = scraper_model_data,
             aes(x = log_settlement_grav_sc,
                 y = mean_biomass),
             colour = "grey40", size = 1.5,
             alpha  = 0.5, inherit.aes = FALSE) +
  labs(x = "log(Settlement gravity) (standardised)",
       y = "Scraper biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())


# ── Chla effect ───────────────────────────────────────────────
# Positive effect — higher productivity supports
# scraper biomass. Most robust predictor across all
# metrics and scales.
# Settlement gravity and rugosity held at mean (z = 0).

chla_grid_scr <- data.frame(
  log_chla_sc            = seq(
    min(scraper_model_data$log_chla_sc),
    max(scraper_model_data$log_chla_sc),
    length.out = 200),
  rugosity_sc            = 0,
  log_settlement_grav_sc = 0
)

chla_pred_scr     <- predict(best_model_s,
                             newdata = chla_grid_scr,
                             se.fit  = TRUE,
                             type    = "response",
                             re.form = NA)
chla_grid_scr$fit <- chla_pred_scr$fit
chla_grid_scr$lwr <- chla_pred_scr$fit - 1.96 * chla_pred_scr$se.fit
chla_grid_scr$upr <- chla_pred_scr$fit + 1.96 * chla_pred_scr$se.fit

p_scr_chla <- ggplot(chla_grid_scr,
                     aes(x = log_chla_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "#1a9641", alpha = 0.15) +
  geom_line(colour = "#1a9641", linewidth = 1.1) +
  geom_point(data = scraper_model_data,
             aes(x = log_chla_sc, y = mean_biomass),
             colour = "grey40", size = 1.5,
             alpha  = 0.5, inherit.aes = FALSE) +
  labs(x = "log(Chlorophyll-a) (standardised)",
       y = "Scraper biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())


# ── Arrange plots ─────────────────────────────────────────────
gridExtra::grid.arrange(p_scr_grav, p_scr_chla, ncol = 2)

# jpeg("scraper_marginal_effects.jpg",
#      width = 22, height = 11, units = "cm", res = 300)
# gridExtra::grid.arrange(p_scr_grav, p_scr_chla, ncol = 2)
# dev.off()


# ============================================================
#  RESULTS SUMMARY
# ============================================================

cat("\n--- Scraper results summary ---\n")
tribble(
  ~Question,   ~Result,           ~Key_finding,
  "Q1",        "Sett. gravity",   "weight = 0.562, DAICc = 2.85 vs baseline, b = -0.261, p = 0.017",
  "Q2 conn",   "Not supported",   "DAICc = 1.28, weight = 0.345, b = -0.105, p = 0.238",
  "Q3 MPA",    "Not supported",   "DAICc = 0.89, weight = 0.390 — medium MPA negative, artefact",
  "Q3 int",    "Not supported",   "MPA x pressure DAICc = 5.58 — not supported",
  "Spatial",   "Clean",           "Moran's I = +0.033, p = 0.230",
  "Sens (a)",  "Consistent",      "market gravity b = -0.198*, pop b = -0.076 ns — gravity metrics both significant",
  "Sens (b)",  "Consistent",      "pressure best at transect level (weight = 0.452), null DAICc = 15.61"
) %>% print()


# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()

# ============================================================
#  FIGURE 2 — Scraper marginal effects
#  (a) Settlement gravity (negative) — red #d73027
#  (b) Chlorophyll-a (positive) — green #1a9641
#  Predictions on response scale (Tweedie glmmTMB)
#  Built separately then arranged
# ============================================================

# ── Panel (a): Settlement gravity ────────────────────────────
grav_grid_scr <- data.frame(
  log_settlement_grav_sc = seq(
    min(scraper_model_data$log_settlement_grav_sc),
    max(scraper_model_data$log_settlement_grav_sc),
    length.out = 200),
  rugosity_sc = 0,
  log_chla_sc = 0
)

grav_pred_scr     <- predict(s_best_q3,
                             newdata = grav_grid_scr,
                             se.fit  = TRUE,
                             type    = "response",
                             re.form = NA)
grav_grid_scr$fit <- grav_pred_scr$fit
grav_grid_scr$lwr <- grav_pred_scr$fit - 1.96 * grav_pred_scr$se.fit
grav_grid_scr$upr <- grav_pred_scr$fit + 1.96 * grav_pred_scr$se.fit

p2a <- ggplot(grav_grid_scr,
              aes(x = log_settlement_grav_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "#d73027", alpha = 0.15) +
  geom_line(colour = "#d73027", linewidth = 1.1) +
  geom_point(data = scraper_model_data,
             aes(x = log_settlement_grav_sc, y = mean_biomass),
             colour = "grey40", size = 1.5, alpha = 0.5,
             inherit.aes = FALSE) +
  labs(x = "Settlement gravity (standardised)",
       y = "Scraper biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ── Panel (b): Chlorophyll-a ──────────────────────────────────
chla_grid_scr <- data.frame(
  log_chla_sc            = seq(
    min(scraper_model_data$log_chla_sc),
    max(scraper_model_data$log_chla_sc),
    length.out = 200),
  rugosity_sc            = 0,
  log_settlement_grav_sc = 0
)

chla_pred_scr     <- predict(s_best_q3,
                             newdata = chla_grid_scr,
                             se.fit  = TRUE,
                             type    = "response",
                             re.form = NA)
chla_grid_scr$fit <- chla_pred_scr$fit
chla_grid_scr$lwr <- chla_pred_scr$fit - 1.96 * chla_pred_scr$se.fit
chla_grid_scr$upr <- chla_pred_scr$fit + 1.96 * chla_pred_scr$se.fit

p2b <- ggplot(chla_grid_scr,
              aes(x = log_chla_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "#1a9641", alpha = 0.15) +
  geom_line(colour = "#1a9641", linewidth = 1.1) +
  geom_point(data = scraper_model_data,
             aes(x = log_chla_sc, y = mean_biomass),
             colour = "grey40", size = 1.5, alpha = 0.5,
             inherit.aes = FALSE) +
  labs(x = "Chlorophyll-a (standardised)",
       y = "Scraper biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ── Arrange ───────────────────────────────────────────────────
gridExtra::grid.arrange(p2a, p2b, ncol = 2)

# jpeg("figure2_scrapers.jpg",
#      width = 22, height = 11, units = "cm", res = 300)
# gridExtra::grid.arrange(p2a, p2b, ncol = 2)
# dev.off()