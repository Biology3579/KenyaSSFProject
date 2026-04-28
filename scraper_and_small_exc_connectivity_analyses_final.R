# ============================================================
#  DRIVERS OF SCRAPER/SMALL EXCAVATOR BIOMASS
#  Chapter 1 — Functional Group Analysis: Scrapers
#
#  Scrapers and small excavators combined as a single
#  functional guild — both process reef substrate through
#  grazing and minor excavation. Distinguished from large
#  excavators by body size and substrate removal depth.
#
#  Analytical framework mirrors total biomass and other
#  functional groups (Q1–Q3 + sensitivity).
#
#  Key differences from other groups:
#    Zero proportion: ~2% at site level (1 site) —
#    Tweedie selected despite near-zero proportion due
#    to extreme influential observation (itsan, 130×
#    below median) causing Gaussian log instability.
#    Transect zeros: ~8.5% — standard Tweedie.
#    Pressure metric: settlement gravity decisively
#    preferred (weight = 0.954) — strongest metric
#    selection of any functional group. Contrasts
#    directly with piscivores (market gravity) reflecting
#    different fishing pressure mechanisms.
#    Ecological context: subsistence/artisanal targeted
#    guild — settlement gravity captures proximity-
#    weighted SSF pressure appropriately.
#
#  Sensitivity analysis:
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
cat("Countries:",         n_distinct(scraper_transects$country), "\n")

# ── Site-level dataset ────────────────────────────────────────
scraper_model_data <- scraper_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass              = mean(transect_scraper_biomass,
                                     na.rm = TRUE),
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

# ── MPA classification check ──────────────────────────────────
cat("\nMPA status counts:\n")
print(table(scraper_model_data$mpa_status))

# ── Distribution plot ─────────────────────────────────────────
par(mfrow = c(1, 2))
hist(scraper_model_data$mean_biomass,
     breaks = 30, main = "Raw",
     xlab = "Mean scraper biomass (g)")
hist(log(scraper_model_data$mean_biomass + 0.01),
     breaks = 25, main = "Log-transformed",
     xlab = "log(Mean biomass + 0.01)")
abline(v = log(0.01), col = "red", lty = 2)
par(mfrow = c(1, 1))

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
# ============================================================

# ── Gaussian log ──────────────────────────────────────────────
scr_lm_global <- lm(
  log(mean_biomass + 0.01) ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  data = scraper_model_data
)

par(mfrow = c(2, 2))
plot(scr_lm_global, main = "Gaussian log")
par(mfrow = c(1, 1))

# ── Tweedie ───────────────────────────────────────────────────
scr_tw_global <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

scr_tw_res <- simulateResiduals(scr_tw_global, n = 500)
plot(scr_tw_res)
testZeroInflation(scr_tw_res)
testDispersion(scr_tw_res)

# ── Family selection decision ─────────────────────────────────
# Tweedie (log link): SELECTED
#   Gaussian log rejected — site itsan (Comoros,
#   mean biomass = 30g, ~130× below median) exceeds
#   Cook's distance threshold (standardised residual
#   ~-4.5, leverage ~0.20). Itsan is a genuine
#   biological observation — not an error — but
#   Gaussian log cannot accommodate its extreme
#   value without substantial coefficient distortion.
#
#   Gaussian log diagnostics:
#     Site 12 (itsan): standardised residual ~-4.5,
#       exceeds Cook's distance 0.5 threshold —
#       both outlier and influential observation.
#     Site 11: high leverage (~0.20), positive
#       residual — within threshold.
#     Site 47: highest leverage (~0.28), small
#       residual — not influential.
#     Q-Q: severe lower tail deviation at site 12.
#     Scale-Location: heteroscedasticity at low end.
#
#   Tweedie DHARMa diagnostics (n = 500):
#     KS test:        [update from plot]
#     Dispersion:     p = 0.160, ratio = 1.457 — acceptable
#     Zero inflation: p = 1.000 — not needed
#   Tweedie handles itsan's extreme value natively
#   through its variance structure without distortion.
#
#  Proceed: glmmTMB(family = tweedie(link = "log")) on
#  raw mean_biomass throughout all scraper analyses.

# ============================================================
#  RANDOM EFFECT STRUCTURE
#  Tweedie family throughout.
#  Global predictor set — most demanding test.
# ============================================================

scr_re_null <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

scr_re_country <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | country),
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

cat("\n--- Scraper RE structure: country-level ---\n")
print(make_aicc_df(list(
  "No RE"         = scr_re_null,
  "(1 | country)" = scr_re_country
)))

# RE result: No RE: AICc = 1010.74, weight = 0.612
# (1|country): ΔAICc = 0.91, weight = 0.389
#
# Country clustering marginally not supported —
# ΔAICc = 0.91 is the weakest RE penalty across
# all functional groups:
#   Total biomass:  ΔAICc = 2.86
#   Browsers:       ΔAICc = 3.03
#   Corallivores:   ΔAICc = 3.03
#   Excavators:     ΔAICc = 3.03
#   Grazers:        ΔAICc = 2.89
#   Piscivores:     ΔAICc = 3.04
#   Scrapers:       ΔAICc = 0.91
#
# Unlike all other groups where country RE was clearly
# not supported (ΔAICc > 2), scrapers show genuine
# model selection uncertainty (ΔAICc < 2, weight =
# 0.389). Some between-country clustering in scraper
# biomass may exist — Comoros (itsan) drives much of
# the between-country variance and the Tweedie may
# not fully absorb this signal.
#
# Decision: no RE retained for consistency with all
# other functional group analyses and because ΔAICc
# does not reach the threshold for clear RE support.
# Acknowledged as a limitation — between-country
# variance in scraper biomass may be partially
# confounded with fixed predictor estimates.
# Sensitivity: if key Q2/Q3 results change when
# (1|country) is included, note this explicitly.

# ============================================================
#  Q1 — PRESSURE METRIC SELECTION
#
#  Settlement gravity decisively preferred in prior
#  analysis (weight = 0.954, ΔAICc = 6.79 over market
#  gravity) — strongest metric selection result across
#  all functional groups. Confirm with updated data.
#  Ecological motivation: scrapers/small excavators
#  targeted by subsistence SSF — settlement gravity
#  (proximity-weighted population pressure) is the
#  appropriate proxy. Contrasts with piscivores
#  (market gravity) and grazers (market gravity).
# ============================================================

# ── Without connectivity control ─────────────────────────────
scr_q1_settgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

scr_q1_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

scr_q1_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_pop_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

cat("\n--- Q1: Scraper pressure metric (without connectivity) ---\n")
print(make_aicc_df(list(
  "Settlement gravity"   = scr_q1_settgrav,
  "Market gravity"       = scr_q1_mktgrav,
  "Settlement pop. 25km" = scr_q1_settpop
)))

# ── With connectivity control ─────────────────────────────────
scr_q1_settgrav_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

scr_q1_mktgrav_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

scr_q1_settpop_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_pop_sc +
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

cat("\n--- Q1: Scraper pressure metric (with connectivity) ---\n")
print(make_aicc_df(list(
  "Settlement gravity"   = scr_q1_settgrav_conn,
  "Market gravity"       = scr_q1_mktgrav_conn,
  "Settlement pop. 25km" = scr_q1_settpop_conn
)))

cat("\n--- Q1: Scraper coefficient summary ---\n")
summary(scr_q1_settgrav)

# ── Q1: Baseline model coefficients (settlement gravity) ──────
#
# glmmTMB(mean_biomass ~ rugosity_sc + log_settlement_grav_sc
#         + log_chla_sc, family = tweedie)
# n = 54 sites, dispersion = 2.88
#
# Rugosity:           β = +0.143, p = 0.126 ns
# Settlement gravity: β = -0.261, p = 0.017 *
# Chla:               β = +0.204, p = 0.037 *
#
# Two significant predictors — the clearest baseline
# signal of any functional group alongside total biomass.
# Contrasts with corallivores, grazers, excavators, and
# piscivores where baseline predictors were largely
# non-significant.
#
# Settlement gravity negative and significant —
# higher proximity-weighted SSF pressure reduces
# scraper biomass. Expected direction, consistent
# with scrapers being directly targeted by subsistence
# fishing. The only functional group besides total
# biomass where the pressure effect is significant
# at baseline level.
#
# Chla positive and significant — higher productivity
# supports scraper biomass. Consistent with scrapers
# exploiting algal resources that are enhanced by
# productive environments. Also consistent with the
# original analysis where chla was a significant
# predictor.
#
# Rugosity non-significant (p = 0.126) — habitat
# complexity not independently important once pressure
# and productivity are included. Positive direction
# maintained throughout.
#
# Dispersion = 2.88 — moderate, lower than piscivores
# (18.2) and excavators (8.5), reflecting less extreme
# variance in scraper biomass across sites.

# ============================================================
#  Q2 — ARE CONNECTIVITY AND MPA BASELINE DRIVERS?
#
#  Tweedie throughout — McFadden pseudo-R².
#  Identical structure to total biomass Q2.
# ============================================================

scr_null <- glmmTMB(
  mean_biomass ~ 1,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

scr_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

scr_baseline_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

scr_baseline_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

scr_global_additive <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

scr_q2_models <- list(
  "Null"                  = scr_null,
  "Baseline"              = scr_baseline,
  "Baseline + conn"       = scr_baseline_conn,
  "Baseline + MPA"        = scr_baseline_mpa,
  "Baseline + conn + MPA" = scr_global_additive
)

cat("\n--- Q2: Scraper model comparison (AICc ranked) ---\n")
print(make_aicc_df(scr_q2_models))

# ── McFadden pseudo-R² relative to baseline ───────────────────
cat("\n--- Q2: Scraper pseudo-R² relative to baseline ---\n")

null_ll_scr     <- as.numeric(logLik(scr_null))
baseline_r2_scr <- 1 - (as.numeric(logLik(scr_baseline)) /
                          null_ll_scr)

scr_q2_models %>%
  imap_dfr(~ tibble(
    Model  = .y,
    McF_R2 = round(1 - (as.numeric(logLik(.x)) / null_ll_scr),
                   3)
  )) %>%
  mutate(
    Delta_R2 = round(McF_R2 - baseline_r2_scr, 3),
    Delta_R2 = ifelse(Model %in% c("Null", "Baseline"),
                      NA, Delta_R2)
  ) %>%
  print()

cat("\n--- Q2: Baseline coefficients ---\n")
summary(scr_baseline)

cat("\n--- Q2: Baseline + connectivity coefficients ---\n")
summary(scr_baseline_conn)

cat("\n--- Q2: Baseline + MPA coefficients ---\n")
summary(scr_baseline_mpa)

# ── Q2 results ────────────────────────────────────────────────
#
# AICc comparison:
#   Baseline:              AICc = 1007.32, weight = 0.426 (BEST)
#   Baseline + MPA:        ΔAICc = 0.89,   weight = 0.272
#   Baseline + conn:       ΔAICc = 1.28,   weight = 0.224
#   Baseline + conn + MPA: ΔAICc = 3.42,   weight = 0.077
#   Null:                  ΔAICc = 12.11,  weight = 0.001
#
# McFadden pseudo-R² relative to baseline (R² = 0.019):
#   Baseline + conn:       ΔR² = +0.002
#   Baseline + MPA:        ΔR² = +0.005
#   Baseline + conn + MPA: ΔR² = +0.005
#
# Baseline best-supported (weight = 0.426) but Baseline
# + MPA and Baseline + conn both competitive (combined
# weight of top three = 0.922). Null strongly rejected
# (ΔAICc = 12.11) — scrapers are the only non-browser
# functional group where the null is decisively rejected.
#
# Baseline + connectivity coefficients:
#   Settlement gravity: β = -0.252, p = 0.020 *
#   Chla:               β = +0.171, p = 0.087 . (weakens)
#   Connectivity:       β = -0.105, p = 0.238 ns
#   Connectivity negative but non-significant — adds
#   minimal variance (ΔR² = +0.002). Settlement gravity
#   signal preserved. Chla weakens when connectivity
#   included — modest collinearity between the two.
#
# Baseline + MPA coefficients:
#   Settlement gravity: β = -0.363, p = 0.002 **
#     Strengthens when MPA included — MPA absorbs
#     governance variance, clarifying the pressure signal.
#   Chla:               β = +0.147, p = 0.127 ns (weakens)
#   MPA low:            β = +0.034, p = 0.897 ns
#   MPA medium:         β = -0.422, p = 0.034 *
#     NEGATIVE medium MPA — counterintuitive. Medium
#     protection associated with lower scraper biomass
#     than unprotected sites. Requires exploratory
#     checks in Q3. Possible explanations: (1) medium
#     MPAs in this dataset coincide with low-habitat-
#     quality sites; (2) scraper biomass responds to
#     total benthic community structure rather than
#     protection level directly; (3) data structure
#     artefact — check medium MPA biomass distribution.
#
# Coefficient stability across Q2 models:
#   Settlement gravity: -0.261* → -0.252* → -0.363**
#   Direction robust throughout — strengthens when
#   MPA absorbs additional variance.
#
# Baseline selected as Q3 reference (weight = 0.426,
# most parsimonious). Negative medium MPA noted —
# investigate in Q3 interaction models.

# ============================================================
#  SPATIAL AUTOCORRELATION CHECK
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

cat("\n--- Spatial autocorrelation: scraper best Q2 model ---\n")
print(moran.test(residuals(scr_baseline,
                           type = "pearson"), listw5_scr))

# ── Spatial autocorrelation result ────────────────────────────
# Moran's I = +0.033, p = 0.230 — no significant spatial
# autocorrelation in scraper baseline model residuals.
#
# Pattern across all analyses:
#   Total biomass: I = +0.140, p = 0.015 — significant
#   Browsers:      I = -0.074, p = 0.802 — no signal
#   Corallivores:  I = -0.051, p = 0.669 — no signal
#   Excavators:    I = -0.026, p = 0.537 — no signal
#   Grazers:       I = +0.187, p = 0.002 — significant
#   Piscivores:    I = +0.034, p = 0.235 — no signal
#   Scrapers:      I = +0.033, p = 0.230 — no signal
#
# Positive but non-significant — similar magnitude to
# piscivores (I = +0.034). No spatial error modelling
# required.
#
# Warnings: identical points and 4 sub-graphs —
# expected, as per all other analyses.

# ============================================================
#  Q3 — DO MANAGEMENT AND CONNECTIVITY MODIFY
#       SCRAPER PRESSURE EFFECTS?
#
#  Reference model: Baseline (best Q2 model, weight = 0.426).
#  Baseline contains: rugosity + settlement_gravity + chla.
#
#  Note: Baseline + MPA competitive (ΔAICc = 0.89) and
#  medium MPA negative in Q2 (β = -0.422, p = 0.034) —
#  H1 (MPA × pressure) and H3 (MPA × connectivity) are
#  of particular interest to understand whether the
#  negative MPA medium effect is conditional on pressure
#  or connectivity.
#
#  Tweedie throughout.
# ============================================================

scr_q3_reference <- scr_baseline

# ── H1: MPA effectiveness depends on fishing intensity ────────
scr_int_mpa_press <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    mpa_status * log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

# ── H2: Connectivity buffers exploitation effects ─────────────
scr_int_conn_press <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc * log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

# ── H3: MPA effectiveness depends on larval supply ────────────
scr_int_mpa_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc +
    mpa_status * connectivity_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

scr_q3_models <- list(
  "Reference (best Q2)"    = scr_q3_reference,
  "H1: MPA × pressure"     = scr_int_mpa_press,
  "H2: Conn × pressure"    = scr_int_conn_press,
  "H3: MPA × connectivity" = scr_int_mpa_conn
)

cat("\n--- Q3: Scraper interaction comparison ---\n")
print(make_aicc_df(scr_q3_models))

cat("\n--- Q3: H3 MPA × connectivity coefficients ---\n")
summary(scr_int_mpa_conn)

# ── H3 exploratory checks ────────────────────────────────────
# Check 1 — VIF
cat("\n--- H3 scraper: VIF check ---\n")
check_collinearity(scr_int_mpa_conn)

# Check 2 — MPA group biomass vs connectivity
cat("\n--- H3 scraper: MPA biomass by group ---\n")
scraper_model_data %>%
  group_by(mpa_status) %>%
  summarise(
    n              = n(),
    mean_biomass   = round(mean(mean_biomass), 1),
    mean_conn      = round(mean(connectivity_sc), 3),
    min_conn       = round(min(connectivity_sc), 3),
    max_conn       = round(max(connectivity_sc), 3),
    .groups = "drop"
  ) %>%
  print()

# Check 3 — interaction plot
conn_seq <- seq(min(scraper_model_data$connectivity_sc),
                max(scraper_model_data$connectivity_sc),
                length.out = 200)

pred_grid_scr <- expand.grid(
  connectivity_sc        = conn_seq,
  mpa_status             = factor(c("none", "low", "medium"),
                                  levels = c("none", "low",
                                             "medium")),
  rugosity_sc            = 0,
  log_chla_sc            = 0,
  log_settlement_grav_sc = 0
)

pred_grid_scr$fit <- predict(scr_int_mpa_conn,
                             newdata = pred_grid_scr,
                             type    = "response")

ggplot(pred_grid_scr,
       aes(x      = connectivity_sc,
           y      = fit,
           colour = mpa_status,
           group  = mpa_status)) +
  geom_line(linewidth = 1.1) +
  geom_rug(data = scraper_model_data,
           aes(x = connectivity_sc),
           inherit.aes = FALSE,
           alpha = 0.4, sides = "b") +
  scale_colour_manual(
    values = c("none"   = "#636363",
               "low"    = "#74a9cf",
               "medium" = "#0570b0"),
    labels = c("No MPA", "Low protection",
               "Medium protection")
  ) +
  labs(x      = "Connectivity (standardised)",
       y      = "Scraper biomass (g)",
       colour = "MPA status") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        legend.position  = "top",
        panel.grid.minor = element_blank())

# Check 4 — DHARMa
scr_h3_sim <- simulateResiduals(scr_int_mpa_conn, n = 500)
plot(scr_h3_sim)
testOutliers(scr_h3_sim)

# ── H3 exploratory checks ─────────────────────────────────────
#
# VIF: CRITICAL CONCERN
#   mpa_status:              VIF = 52.80 — severe
#   mpa_status:connectivity: VIF = 72.36 — severe
#   Near-perfect collinearity between MPA levels and
#   connectivity. Interaction coefficient estimates
#   are unreliable.
#
# Data structure: LOW MPA × CONNECTIVITY CONFOUND
#   No MPA:     mean connectivity = -0.345 (n = 30)
#   Low MPA:    mean connectivity = +0.879 (n = 7)
#               range: 0.662–1.00 — ALL at high connectivity
#   Medium MPA: mean connectivity = +0.303 (n = 17)
#
#   All 7 low MPA sites cluster at high connectivity —
#   no low MPA sites at low connectivity to independently
#   estimate the MPA low main effect. The model cannot
#   separate the effect of low MPA status from the
#   effect of high connectivity. Large opposing
#   coefficients (MPA low: β = +4.718; MPA low ×
#   connectivity: β = -5.396) mathematically cancel
#   but are not independently identified.
#
# DHARMa (n = 500):
#   KS test:        p = 0.596 — good fit
#   Dispersion:     p = 0.100 — acceptable
#   Outlier test:   p = 0.194 — no significant outliers
#   Residuals vs predicted: systematic downward trend
#     across all three quantile lines — model overpredicts
#     at high fitted values, consistent with collinearity
#     distorting predictions.
#
# Conclusion: H3 scraper interaction is a DATA
# STRUCTURE ARTEFACT driven by confounding between
# low MPA status and high connectivity. VIF > 50
# for MPA terms confirms severity. Despite decisive
# AICc support (weight = 0.910, ΔAICc = 5.11),
# interaction coefficients are not interpretable.
# H3 not supported as a genuine ecological signal.
#
# Differs from other artefacts identified:
#   Total biomass H1: extrapolation beyond data range
#   Excavator H1:     single-site influence (kifinge)
#   Scraper H3:       group-level confounding between
#                     categorical predictor and modifier
#
# Overall Q3 conclusion:
# No interaction supported for scrapers. Baseline
# remains the most defensible model. Settlement
# gravity (negative) and chla (positive) are the
# primary drivers. Neither management nor connectivity
# modifies the pressure-biomass relationship in a
# manner reliably estimable from the current data.
# The low MPA × connectivity confound is a fundamental
# data limitation — resolving it would require either
# more low MPA sites at varying connectivity levels
# or a different analytical approach.


# ============================================================
#  SENSITIVITY ANALYSIS
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
# Settlement gravity decisively preferred (weight = 0.954)
# — alternative metrics are genuinely poor proxies for
# subsistence SSF pressure on scrapers. Sensitivity
# analysis documents metric-dependence but primary
# robustness test is transect-level replication (b).

scr_sens_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

scr_sens_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_settlement_pop_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

cat("\n--- Sensitivity (a): scraper alternative metrics ---\n")
cat("Market gravity:\n")
print(summary(scr_sens_mktgrav)$coefficients$cond)
cat("\nSettlement population 25km:\n")
print(summary(scr_sens_settpop)$coefficients$cond)

# ── Sensitivity (a) results ───────────────────────────────────
#
# Alternative pressure metrics — baseline structure retained,
# only pressure metric substituted.
#
# Market gravity (β = -0.198, p = 0.049 *):
#   Significant and negative — pressure effect preserved
#   with market gravity despite it being a weaker proxy
#   for subsistence SSF pressure on scrapers.
#   Chla strengthens (β = +0.283, p = 0.001**) —
#   more precisely estimated when market gravity
#   substituted for settlement gravity.
#   Rugosity non-significant (p = 0.141).
#
# Settlement population 25km (β = -0.076, p = 0.430):
#   Not significant — settlement population fails to
#   capture SSF pressure on scrapers. Chla remains
#   strongly significant (β = +0.313, p < 0.001**).
#   Consistent with Q1 result where settlement population
#   was strongly rejected (ΔAICc > 4).
#
# Key patterns across all three metrics:
# (1) Chla positive and significant across all metrics
#     (β = +0.171–0.313, p = 0.037–0.001) — the most
#     robust predictor for scrapers, metric-independent.
# (2) Pressure metrics negative in direction throughout
#     but significance varies: settlement gravity p = 0.017*,
#     market gravity p = 0.049*, settlement pop. p = 0.430.
#     The two gravity metrics both significant — settlement
#     population not. Consistent with gravity-based
#     framework outperforming population size alone.
# (3) Rugosity non-significant across all metrics —
#     habitat complexity not independently important
#     once pressure and productivity included.
#
# Primary Q2 conclusion robust: negative pressure
# effect (significant with both gravity metrics) and
# positive chla effect (consistent throughout).
# Settlement gravity retained as primary metric —
# strongest signal and most ecologically appropriate
# for subsistence SSF on scrapers.

# ── (b) Transect-level replication ───────────────────────────
# ~8.5% zeros at transect level — standard Tweedie.
# ZI Tweedie tested for completeness.
# Prior result: standard Tweedie strongly preferred
# (weight = 0.928, ΔAICc = 5.12).

scr_trans_tw <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data
)

scr_trans_tw_zi <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = scraper_transect_data
)

scr_trans_res    <- simulateResiduals(scr_trans_tw,    n = 500)
scr_trans_res_zi <- simulateResiduals(scr_trans_tw_zi, n = 500)

plot(scr_trans_res);    testZeroInflation(scr_trans_res)
plot(scr_trans_res_zi); testZeroInflation(scr_trans_res_zi)

cat("\n--- Sensitivity (b): transect family selection ---\n")
print(make_aicc_df(list(
  "Tweedie"    = scr_trans_tw,
  "ZI Tweedie" = scr_trans_tw_zi
)))

# ── Transect family selection ─────────────────────────────────
#
# Tweedie:    AICc = 4334.06, weight = 0.941 (SELECTED)
# ZI Tweedie: ΔAICc = 5.53,   weight = 0.059
#
# Zero inflation tests (DHARMa, n = 500):
#   Tweedie:    ratio = 0.967, p = 0.956 — not significant
#   ZI Tweedie: ratio = 0.990, p = 1.000 — not significant
#
# Standard Tweedie strongly preferred — ΔAICc = 5.53,
# the largest family selection margin across all
# transect-level analyses:
#   Browsers:    ΔAICc = [update]
#   Excavators:  ΔAICc = 2.19
#   Piscivores:  ΔAICc = 2.19
#   Scrapers:    ΔAICc = 5.53 — strongest
#
# ZI clearly not warranted — zero inflation non-
# significant for both models and ZI Tweedie strongly
# penalised by AICc. ~8.5% transect zeros handled
# natively by standard Tweedie compound Poisson
# structure.
#
# Proceed: standard Tweedie + (1|site) throughout
# scraper transect-level analyses.

# ── Transect Q2 sequence + Q3 best model ─────────────────────
scr_trans_null <- glmmTMB(
  transect_scraper_biomass ~ 1 + (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data
)

scr_trans_baseline <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data
)

scr_trans_conn <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data
)

scr_trans_mpa <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data
)

scr_trans_global <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    connectivity_sc +
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data
)

# Q3 best model replicated at transect level
scr_trans_mpa_conn <- glmmTMB(
  transect_scraper_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    mpa_status * connectivity_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_transect_data
)

cat("\n--- Sensitivity (b): scraper transect comparison ---\n")
print(make_aicc_df(list(
  "Null"                  = scr_trans_null,
  "Baseline"              = scr_trans_baseline,
  "Baseline + conn"       = scr_trans_conn,
  "Baseline + MPA"        = scr_trans_mpa,
  "Baseline + conn + MPA" = scr_trans_global,
  "MPA × connectivity"    = scr_trans_mpa_conn
)))

cat("\n--- Sensitivity (b): best model coefficients ---\n")
summary(scr_trans_mpa_conn)

# ── Sensitivity (b) results ───────────────────────────────────
#
# AICc comparison (transect level, n = 243):
#   MPA × connectivity:    AICc = 4325.44, weight = 0.853 (BEST)
#   Baseline:              ΔAICc = 5.17,   weight = 0.064
#   Baseline + conn:       ΔAICc = 6.01,   weight = 0.042
#   Baseline + MPA:        ΔAICc = 6.75,   weight = 0.029
#   Baseline + conn + MPA: ΔAICc = 8.62,   weight = 0.011
#   Null:                  ΔAICc = 20.78,  weight = 0.000
#
# MPA × connectivity dominant at transect level
# (weight = 0.853) — even stronger than site level
# (weight = 0.910 vs 0.853 — similar). Null strongly
# rejected (ΔAICc = 20.78). Model ordering replicates
# site level exactly.
#
# Best model coefficients (transect level):
#   Rugosity:                    β = +0.200, p = 0.011 *
#   Settlement gravity:          β = -0.342, p = 0.002 **
#   Chla:                        β = +0.199, p = 0.044 *
#   MPA low:                     β = +4.825, p < 0.001 ***
#   MPA medium:                  β = -0.326, p = 0.097 .
#   Connectivity:                β = -0.119, p = 0.311 ns
#   MPA low × connectivity:      β = -5.439, p < 0.001 ***
#   MPA medium × connectivity:   β = +0.255, p = 0.158 ns
#   Site variance: 0.129 (SD = 0.360)
#   Dispersion: 65.4 — high within-site variance
#
# Coefficient stability across scales:
#   Settlement gravity: -0.340** (site) → -0.342** (transect)
#   Chla:               +0.201*  (site) → +0.199*  (transect)
#   MPA low:            +4.718** (site) → +4.825** (transect)
#   MPA low × conn:     -5.396** (site) → -5.439** (transect)
#   Near-identical across scales — within-site variation
#   does not alter the pattern.
#
# HOWEVER: the collinearity problem identified at site
# level (VIF > 50 for MPA terms) applies equally here.
# All 7 low MPA sites still cluster at high connectivity
# regardless of analytical scale. The transect-level
# replication confirms the pattern is not a site-level
# aggregation artefact — but it cannot resolve the
# fundamental data structure confound. The interaction
# coefficients remain unreliable for the same reason.
#
# The transect result strengthens confidence that
# something real is happening with low MPA sites at
# high connectivity — but the collinearity prevents
# reliable estimation of the separate MPA and
# connectivity effects. This is a data limitation
# that cannot be resolved analytically.
#
# Baseline predictors fully replicated at transect
# level (settlement gravity negative**, chla positive*)
# — these are the robust findings for scrapers.

# ============================================================
#  MARGINAL EFFECT PLOTS
#  No interaction supported — Q3 artefact confirmed.
#  Plot baseline significant predictors:
#  settlement gravity (negative) and chla (positive).
#  Tweedie — predictions on response scale (raw biomass).
# ============================================================

# ── Settlement gravity effect ─────────────────────────────────
grav_grid_scr <- data.frame(
  log_settlement_grav_sc = seq(
    min(scraper_model_data$log_settlement_grav_sc),
    max(scraper_model_data$log_settlement_grav_sc),
    length.out = 200),
  rugosity_sc = 0,
  log_chla_sc = 0
)

grav_pred_scr <- predict(scr_baseline,
                         newdata = grav_grid_scr,
                         se.fit  = TRUE,
                         type    = "response")

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
             alpha = 0.5, inherit.aes = FALSE) +
  labs(x = "log(Settlement gravity) (standardised)",
       y = "Scraper biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# ── Chla effect ───────────────────────────────────────────────
chla_grid_scr <- data.frame(
  log_chla_sc            = seq(
    min(scraper_model_data$log_chla_sc),
    max(scraper_model_data$log_chla_sc),
    length.out = 200),
  rugosity_sc            = 0,
  log_settlement_grav_sc = 0
)

chla_pred_scr <- predict(scr_baseline,
                         newdata = chla_grid_scr,
                         se.fit  = TRUE,
                         type    = "response")

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
             alpha = 0.5, inherit.aes = FALSE) +
  labs(x = "log(Chla) (standardised)",
       y = "Scraper biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())

gridExtra::grid.arrange(p_scr_grav, p_scr_chla, ncol = 2)


# ============================================================
#  RESULTS SUMMARY
# ============================================================

scraper_results <- tribble(
  ~Question,  ~Best_model,      ~Key_finding,
  "Q1",       "Sett. gravity",  "weight = 0.650 — clear preference; contrasts with piscivores",
  "Q2 conn",  "Baseline",       "ΔAICc = 1.28, competitive but ΔR² = +0.002 only",
  "Q2 MPA",   "Baseline",       "ΔAICc = 0.89, competitive; medium MPA β = -0.422* negative",
  "Q3 H1",    "Reference",      "ΔAICc = 10.69 — not supported",
  "Q3 H2",    "Reference",      "ΔAICc = 8.26 — not supported",
  "Q3 H3",    "H3 (artefact)",  "weight = 0.910 but VIF > 50 — low MPA/connectivity confound"
)

cat("\n--- Scraper results summary ---\n")
print(scraper_results)


# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()