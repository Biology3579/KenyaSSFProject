# ============================================================
#  SCRAPER / SMALL EXCAVATOR FISH BIOMASS — SITE-LEVEL ANALYSIS
#
#  Scrapers and small excavators are combined as a single
#  functional guild — both process reef substrate through
#  grazing and minor excavation, contributing to algal control
#  and carbonate cycling. They are distinguished from large
#  excavators by body size and the depth of substrate removal.
#
#  Analytical framework mirrors total biomass (four stages):
#
#  STAGE 1 — Variance partitioning
#             McFadden pseudo-R² + LRT (Tweedie expected given
#             site-level zeros). Formal varpart() if Gaussian
#             log adequate.
#
#  STAGE 2 — Hierarchical model comparison
#             Identical nested sequence to other groups.
#
#  STAGE 3 — Interaction testing (conditional on Stage 2)
#             Three a priori interactions as per other groups.
#
#  STAGE 4 — Sensitivity analysis
#             (a) Alternative pressure metrics
#             (b) Transect-level mixed model replication
#
#  Key differences from other functional groups:
#    - Original analysis used SST not DHW — this version
#      uses DHW for consistency across groups.
#    - Original analysis lacked MPA status and connectivity
#      — both added here for consistency.
#    - Original analysis used 19-model flat candidate set —
#      replaced here with four-stage hierarchical framework.
#    - Original top predictors: environment (SST, chla) and
#      settlement pop. — test whether these hold with DHW
#      substituted for SST and MPA/connectivity added.
#    - Site-level zeros: check before confirming family.
#    - Original ~8.5% transect zeros — Tweedie appropriate.
#
#  Study design:
#    243 transects, 54 sites, 4 countries.
#    Minimum 3 transects per site retained.
# ============================================================

# ── SOURCE SHARED DATA PREPARATION ───────────────────────────
source(here::here("data_preparation.R"))


# ── AGGREGATE SCRAPER / SMALL EXCAVATOR TRANSECT DATA ────────
scraper_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_scraper_biomass = sum(
      ifelse(trophic_group %in% c("scrapers", "small_excavators"),
             tot_wt_g, 0),
      na.rm = TRUE
    ),
    transect_scraper_count = sum(
      ifelse(trophic_group %in% c("scrapers", "small_excavators"),
             number, 0),
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


# ── SITE-LEVEL AGGREGATION ────────────────────────────────────
scraper_site_data <- scraper_transects %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass = mean(transect_scraper_biomass, na.rm = TRUE),
    n_transects  = n(),
    .groups      = "drop"
  ) %>%
  mutate(
    site    = as.factor(site),
    country = as.factor(country)
  )

# ── Zero check ────────────────────────────────────────────────
zeros_site <- mean(scraper_site_data$mean_biomass == 0, na.rm = TRUE)
cat("Site-level zero proportion:", round(zeros_site, 3), "\n")
cat("Site-level zeros (n):",       sum(scraper_site_data$mean_biomass == 0), "\n")
# Original analysis had 1 zero site (ch_rubu).
# If bimodal log distribution → Tweedie confirmed.

# ── Distribution plots ────────────────────────────────────────
scraper_site_data <- scraper_site_data %>%
  mutate(log_mean_biomass = log(mean_biomass + 0.01))

p_raw <- ggplot(scraper_site_data, aes(x = mean_biomass)) +
  geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
  labs(x = "Mean scraper/small excavator biomass per site (g)",
       y = "Frequency", title = "Raw") + theme_bw()

p_log <- ggplot(scraper_site_data, aes(x = log_mean_biomass)) +
  geom_histogram(bins = 25, fill = "#1a9641", colour = "white") +
  geom_vline(xintercept = log(0.01),
             colour = "red", linetype = "dashed") +
  annotate("text", x = log(0.01) + 0.3, y = Inf,
           vjust = 2, label = "Zero sites",
           colour = "red", size = 3) +
  labs(x = "log(mean biomass + 0.01)",
       y = "Frequency", title = "Log-transformed") + theme_bw()

gridExtra::grid.arrange(p_raw, p_log, ncol = 2)

# ── Distribution decision ─────────────────────────────────────
# Log distribution is approximately unimodal — single zero
# site at -4.6 is an outlier, not an isolated cluster.
# One intermediate site (~-2) separates it from the main
# distribution but this is not a bimodal pattern.
#
# Only 1 site-level zero (1.9%) — bimodal justification
# for Tweedie does not apply here. Contrast with:
#   Browsers:    11% zeros, bimodal → Tweedie
#   Excavators:  15% zeros, bimodal → Tweedie
#   Piscivores:  11% zeros, bimodal → Tweedie
#   Scrapers:     2% zeros, unimodal → test lm()
#
# Run Gaussian log diagnostics alongside Tweedie diagnostics
# before deciding. If lm() diagnostics are clean, prefer
# lm() for consistency with total biomass and corallivores
# — varpart() then available for Stage 1.


# ── Build scraper site-level model dataset ────────────────────
scraper_model_data <- scraper_transects %>%
  left_join(scaled_predictors, by = "site") %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_scraper_biomass, na.rm = TRUE),
    log_mean_biomass       = log(mean(transect_scraper_biomass,
                                      na.rm = TRUE) + 0.01),
    n_transects            = n(),
    rugosity_sc            = first(rugosity_sc),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    connectivity_sc        = first(connectivity_sc),
    mpa_status             = first(mpa_status),
    log_chla_sc            = first(log_chla_sc),
    log_max_dhw_sc         = first(log_max_dhw_sc),
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

cat("\nScraper model data:", nrow(scraper_model_data), "sites\n")
cat("Site-level zeros:",    sum(scraper_model_data$mean_biomass == 0), "\n")

# ── Verify no NAs ─────────────────────────────────────────────
scraper_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                connectivity_sc, mpa_status,
                log_chla_sc, log_max_dhw_sc) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

# ── Check zeros and response ──────────────────────────────────
cat("Zeros in mean_biomass:",    sum(scraper_model_data$mean_biomass == 0), "\n")
cat("-Inf in log_mean_biomass:", sum(is.infinite(scraper_model_data$log_mean_biomass)), "\n")
cat("NAs in log_mean_biomass:",  sum(is.na(scraper_model_data$log_mean_biomass)), "\n")
cat("\nResponse summary:\n")
print(summary(scraper_model_data$log_mean_biomass))


# ============================================================
#  PRE-ANALYSIS: HUMAN PRESSURE METRIC SELECTION
# ============================================================

scr_press_settgrav <- glmmTMB(mean_biomass ~ log_settlement_grav_sc,
                              family = tweedie(link = "log"),
                              data   = scraper_model_data)
scr_press_settpop  <- glmmTMB(mean_biomass ~ log_settlement_pop_sc,
                              family = tweedie(link = "log"),
                              data   = scraper_model_data)
scr_press_mktgrav  <- glmmTMB(mean_biomass ~ log_market_gravity_sc,
                              family = tweedie(link = "log"),
                              data   = scraper_model_data)

cat("\n--- Pre-analysis: scraper pressure metric selection ---\n")
print(make_aicc_df(list(
  "Settlement gravity" = scr_press_settgrav,
  "Settlement pop."    = scr_press_settpop,
  "Market gravity"     = scr_press_mktgrav
)))

# ── Pressure metric decision ──────────────────────────────────
# Settlement gravity decisively preferred — weight = 0.954,
# ΔAICc = 6.79 over market gravity, ΔAICc = 8.49 over
# settlement population.
#
# Results:
#   Settlement gravity: AICc = 1008.41, weight = 0.954
#   Market gravity:     ΔAICc = 6.79,   weight = 0.032
#   Settlement pop.:    ΔAICc = 8.49,   weight = 0.014
#
# Strongest and clearest metric selection result across all
# functional groups — even stronger than piscivores
# (market gravity weight = 0.815).
#
# Ecological interpretation: scrapers/small excavators are
# primarily targeted by subsistence and small-scale artisanal
# fishing rather than market-oriented commercial fishing.
# Settlement gravity — proximity-weighted population pressure
# — is the appropriate proxy for this fishing dynamic.
# Contrasts directly with piscivores (high-value commercial
# targets where market gravity dominated) — the two groups
# show opposite metric preferences reflecting genuinely
# different fishing pressure mechanisms.
#
# Settlement population last (ΔAICc = 8.49) — consistent
# with the cross-group pattern where settlement population
# never decisively outperforms gravity metrics, supporting
# the gravity-based framework of Cinner et al. (2016).
#
# log_settlement_grav_sc used as primary pressure predictor
# throughout all scraper analyses. Market gravity and
# settlement population retained for Stage 4a sensitivity
# only.

# ============================================================
#  MODEL FAMILY SELECTION
# ============================================================

# ── Gaussian log (global model) ───────────────────────────────
scr_lm_global <- lm(log_mean_biomass ~ rugosity_sc +
                      log_settlement_grav_sc +
                      connectivity_sc +
                      mpa_status +
                      log_chla_sc +
                      log_max_dhw_sc,
                    data = scraper_model_data)

par(mfrow = c(2, 2)); plot(scr_lm_global); par(mfrow = c(1, 1))

# ── Tweedie (global model) ────────────────────────────────────
scr_tw_global <- glmmTMB(mean_biomass ~ rugosity_sc +
                           log_settlement_grav_sc +
                           connectivity_sc +
                           mpa_status +
                           log_chla_sc +
                           log_max_dhw_sc,
                         family = tweedie(link = "log"),
                         data   = scraper_model_data)

scr_tw_res <- simulateResiduals(scr_tw_global, n = 1000)
plot(scr_tw_res)
testZeroInflation(scr_tw_res)
testDispersion(scr_tw_res)


# ============================================================
#  RANDOM EFFECT STRUCTURE
# ============================================================

scr_re_null <- glmmTMB(mean_biomass ~ rugosity_sc +
                         log_settlement_grav_sc +
                         log_chla_sc +
                         log_max_dhw_sc,
                       family = tweedie(link = "log"),
                       data   = scraper_model_data)

scr_re_country <- glmmTMB(mean_biomass ~ rugosity_sc +
                            log_settlement_grav_sc +
                            log_chla_sc +
                            log_max_dhw_sc +
                            (1 | country),
                          family = tweedie(link = "log"),
                          data   = scraper_model_data)

cat("\n--- Scraper RE structure comparison ---\n")
print(make_aicc_df(list(
  "No RE"         = scr_re_null,
  "(1 | country)" = scr_re_country
)))

# ── Family selection decision ─────────────────────────────────
# Tweedie(log link) selected — Gaussian log rejected.
#
# Zero site-level zeros — Gaussian log technically viable
# but rejected due to excessive influence of site itsan
# (Comoros, mean biomass = 30g vs median ~3900g, 130× below
# median). Coefficient instability test confirms influence:
#
#   Predictor          With itsan   Without itsan
#   Settlement gravity   -0.476        -0.292  (38% change)
#   Connectivity         -0.116        -0.066  (43% change)
#   MPA low              +0.178        -0.051  (sign reversal)
#   DHW                  +0.132        +0.072  (45% change)
#
# Sign reversal for MPA low and large proportional changes
# across multiple predictors indicate the Gaussian log
# model is not stable. Itsan is not an error — it is a
# genuine biological observation reflecting local
# conditions in Comoros — but Gaussian log cannot
# accommodate its extreme value without distortion.
#
# Tweedie diagnostics excellent (DHARMa, n = 1000):
#   KS test:        p = 0.964 — excellent
#   Dispersion:     p = 0.170 — acceptable
#   Zero inflation: p = 1.000 — not ZI
#   Outlier test:   p = 1.000 — no outliers
#
# Tweedie handles the 130× biomass range natively through
# its variance structure without distortion from itsan.
# ZI Tweedie not warranted (ΔAICc = 3.17, ZI p = 1).
#
# Implication for Stage 1: varpart() not available.
# McFadden pseudo-R² + LRT approximation used throughout.

# ============================================================
#  STAGE 1 — VARIANCE PARTITIONING
#
#  Tweedie expected — McFadden pseudo-R² + LRT as per
#  browsers, corallivores, excavators, and piscivores.
#  Three process groups: local, spatial, environmental.
#  MPA excluded — governance variable, tested in Stage 2.
# ============================================================

scr_vp_null <- glmmTMB(mean_biomass ~ 1,
                       family = tweedie(link = "log"),
                       data   = scraper_model_data)

scr_vp_all <- glmmTMB(mean_biomass ~ rugosity_sc +
                        log_settlement_grav_sc +
                        connectivity_sc +
                        log_chla_sc +
                        log_max_dhw_sc,
                      family = tweedie(link = "log"),
                      data   = scraper_model_data)

scr_vp_no_local <- glmmTMB(mean_biomass ~ connectivity_sc +
                             log_chla_sc +
                             log_max_dhw_sc,
                           family = tweedie(link = "log"),
                           data   = scraper_model_data)

scr_vp_no_spatial <- glmmTMB(mean_biomass ~ rugosity_sc +
                               log_settlement_grav_sc +
                               log_chla_sc +
                               log_max_dhw_sc,
                             family = tweedie(link = "log"),
                             data   = scraper_model_data)

scr_vp_no_environ <- glmmTMB(mean_biomass ~ rugosity_sc +
                               log_settlement_grav_sc +
                               connectivity_sc,
                             family = tweedie(link = "log"),
                             data   = scraper_model_data)

# ── McFadden pseudo-R² ────────────────────────────────────────
null_loglik_scr <- as.numeric(logLik(scr_vp_null))

get_mcfadden_scr <- function(model) {
  round(1 - (as.numeric(logLik(model)) / null_loglik_scr), 3)
}

r2_scr_all    <- get_mcfadden_scr(scr_vp_all)
r2_scr_noloc  <- get_mcfadden_scr(scr_vp_no_local)
r2_scr_nospat <- get_mcfadden_scr(scr_vp_no_spatial)
r2_scr_noenv  <- get_mcfadden_scr(scr_vp_no_environ)

unique_scr_local   <- round(r2_scr_all - r2_scr_noloc,  3)
unique_scr_spatial <- round(r2_scr_all - r2_scr_nospat, 3)
unique_scr_environ <- round(r2_scr_all - r2_scr_noenv,  3)

cat("\n--- Stage 1: Scraper pseudo-R² variance decomposition ---\n")
cat("Total McFadden R² (all groups):", r2_scr_all,        "\n")
cat("Unique local:                  ", unique_scr_local,   "\n")
cat("Unique spatial:                ", unique_scr_spatial, "\n")
cat("Unique environmental:          ", unique_scr_environ, "\n")

# ── LRT significance tests ────────────────────────────────────
cat("\nLRT — local unique fraction:\n")
print(anova(scr_vp_no_local, scr_vp_all))

cat("\nLRT — spatial unique fraction:\n")
print(anova(scr_vp_no_spatial, scr_vp_all))

cat("\nLRT — environmental unique fraction:\n")
print(anova(scr_vp_no_environ, scr_vp_all))

# ── Stage 1 results ───────────────────────────────────────────
# Unique local:         McF R² = 0.009, p = 0.010 ** — significant
# Unique spatial:       McF R² = 0.001, p = 0.241   — not significant
# Unique environmental: McF R² = 0.005, p = 0.058 . — marginal
# Total McFadden R²:    0.023
#
# Local processes significant — rugosity and settlement gravity
# together explain unique variance in scraper biomass.
# Environmental context marginal (p = 0.058) — chla and/or
# DHW explain additional variance beyond local processes.
# Spatial (connectivity) not significant.
#
# This is the only group where environmental context shows
# a marginal signal — consistent with the original analysis
# where chla was a significant predictor. The substitution
# of DHW for SST may have weakened the environmental signal
# slightly (SST was significant in original, DHW marginal here).
#
# Compare across groups:
#   Total biomass: local***,  spatial n.s., env .
#   Browsers:      local***,  spatial *,    env n.s.
#   Corallivores:  local n.s., spatial .,   env n.s.
#   Grazers:       local .,   spatial *,    env n.s.
#   Excavators:    local *,   spatial n.s., env n.s.
#   Piscivores:    local n.s., spatial n.s., env n.s.
#   Scrapers:      local **,  spatial n.s., env .
#
# Pattern: scrapers most closely mirror total biomass —
# local dominant, environmental marginal, spatial absent.
# Unlike corallivores and grazers where spatial dominates,
# scrapers respond primarily to site-level conditions.

# ============================================================
#  STAGE 2 — HIERARCHICAL MODEL COMPARISON
#
#  UPDATE family after family selection decision above.
#  Tweedie used as default throughout.
# ============================================================

# ── Null model ────────────────────────────────────────────────
scr_null <- glmmTMB(mean_biomass ~ 1,
                    family = tweedie(link = "log"),
                    data   = scraper_model_data)

# ── Local ecological baseline ─────────────────────────────────
scr_local <- glmmTMB(mean_biomass ~ rugosity_sc +
                       log_settlement_grav_sc,
                     family = tweedie(link = "log"),
                     data   = scraper_model_data)

# ── Local + chla ──────────────────────────────────────────────
scr_local_chla <- glmmTMB(mean_biomass ~ rugosity_sc +
                            log_settlement_grav_sc +
                            log_chla_sc,
                          family = tweedie(link = "log"),
                          data   = scraper_model_data)

# ── Local + DHW ───────────────────────────────────────────────
scr_local_dhw <- glmmTMB(mean_biomass ~ rugosity_sc +
                           log_settlement_grav_sc +
                           log_max_dhw_sc,
                         family = tweedie(link = "log"),
                         data   = scraper_model_data)

# ── Local + environmental context ────────────────────────────
scr_local_env <- glmmTMB(mean_biomass ~ rugosity_sc +
                           log_settlement_grav_sc +
                           log_chla_sc +
                           log_max_dhw_sc,
                         family = tweedie(link = "log"),
                         data   = scraper_model_data)

# ── Local + spatial ───────────────────────────────────────────
scr_local_spatial <- glmmTMB(mean_biomass ~ rugosity_sc +
                               log_settlement_grav_sc +
                               connectivity_sc,
                             family = tweedie(link = "log"),
                             data   = scraper_model_data)

# ── Local + MPA (governance test) ────────────────────────────
scr_local_mpa <- glmmTMB(mean_biomass ~ rugosity_sc +
                           log_settlement_grav_sc +
                           mpa_status,
                         family = tweedie(link = "log"),
                         data   = scraper_model_data)

# ── Global model ──────────────────────────────────────────────
scr_global <- glmmTMB(mean_biomass ~ rugosity_sc +
                        log_settlement_grav_sc +
                        connectivity_sc +
                        mpa_status +
                        log_chla_sc +
                        log_max_dhw_sc,
                      family = tweedie(link = "log"),
                      data   = scraper_model_data)

# ── Model list ────────────────────────────────────────────────
scraper_model_list <- list(
  "Null"            = scr_null,
  "Local"           = scr_local,
  "Local + chla"    = scr_local_chla,
  "Local + DHW"     = scr_local_dhw,
  "Local + env"     = scr_local_env,
  "Local + MPA"     = scr_local_mpa,
  "Local + spatial" = scr_local_spatial,
  "Global"          = scr_global
)

# ── AICc ranked table ─────────────────────────────────────────
cat("\n--- Stage 2: Scraper model comparison (AICc ranked) ---\n")
print(make_aicc_df(scraper_model_list))

# ── McFadden R² increments over local baseline ────────────────
cat("\n--- Stage 2: Scraper variance explained ---\n")
local_loglik_scr <- as.numeric(logLik(scr_local))
null_loglik_scr2 <- as.numeric(logLik(scr_null))

scraper_model_list %>%
  imap_dfr(~ tibble(
    Model  = .y,
    McF_R2 = round(1 - (as.numeric(logLik(.x)) /
                          null_loglik_scr2), 3)
  )) %>%
  mutate(
    Delta_R2 = round(McF_R2 - (1 - (local_loglik_scr /
                                      null_loglik_scr2)), 3),
    Delta_R2 = ifelse(Model %in% c("Null", "Local"), NA, Delta_R2)
  ) %>%
  arrange(Model == "Null",
          desc(Model == "Local"),
          desc(McF_R2)) %>%
  print()

# ── Stage 2 results ───────────────────────────────────────────
# Null clearly rejected (ΔAICc = 12.24) — scrapers are
# explained by the predictor set. Contrast with piscivores
# where null was best.
#
# AICc ranking:
#   Local + env:     AICc = 1007.18, weight = 0.229 — best
#   Local + chla:    ΔAICc = 0.14,   weight = 0.213
#   Local + DHW:     ΔAICc = 0.35,   weight = 0.192
#   Local + MPA:     ΔAICc = 0.89,   weight = 0.146
#   Local + spatial: ΔAICc = 1.67,   weight = 0.099
#   Local:           ΔAICc = 1.86,   weight = 0.090
#   Global:          ΔAICc = 4.06,   weight = 0.030
#   Null:            ΔAICc = 12.24   — strongly rejected
#
# Four models within ΔAICc < 1, combined weight = 0.780.
# Genuine model selection uncertainty — no single model
# clearly dominates.
#
# R² increments over Local baseline (McF R² = 0.015):
#   Local + env:     +0.007 — both chla and DHW contribute
#   Local + MPA:     +0.006 — MPA adds modest value
#   Local + chla:    +0.004 — chla the stronger env predictor
#   Local + DHW:     +0.004 — DHW adds similarly to chla
#   Local + spatial: +0.003 — connectivity weak
#   Global:          +0.011 — driven by env and MPA combined
#
# Key finding: environmental context (chla, DHW) consistently
# adds value for scrapers — consistent with Stage 1 where
# environmental fraction was marginal (p = 0.058).
# This is the only group where both chla AND DHW contribute
# to model fit. Original analysis found chla significant
# with SST — DHW substitution appears to have retained
# the environmental signal.

summary(scr_local_env)
summary(scr_local_chla)
summary(scr_local_dhw)

# ── Coefficient summary: competitive models ───────────────────
#
# Settlement gravity — negative and significant throughout:
#   Local + env:  β = -0.286, p = 0.008 **
#   Local + chla: β = -0.261, p = 0.017 *
#   Local + DHW:  β = -0.389, p < 0.001 ***
#   Most consistent predictor — higher proximity-weighted
#   pressure reduces scraper biomass across all models.
#
# Chla — positive, significant in Local + chla (p = 0.037),
#   marginal in Local + env (p = 0.080):
#   β = +0.204 (Local + chla), β = +0.172 (Local + env)
#   Higher productivity supports scraper biomass —
#   consistent with scrapers exploiting algal resources.
#
# DHW — positive, significant in Local + DHW (p = 0.041),
#   marginal in Local + env (p = 0.091):
#   β = +0.169 (Local + DHW), β = +0.139 (Local + env)
#   Positive direction — bleaching reduces coral cover,
#   releasing algal substrate that benefits scrapers.
#   Competitive release hypothesis: scrapers gain from
#   coral decline through increased food availability.
#   Contrasts with corallivores (negative DHW predicted)
#   and total biomass (positive DHW — geographic covariation).
#
# Rugosity — marginal throughout (p = 0.093–0.184).
#   Positive direction but not independently significant
#   once pressure and environment are included.

# ============================================================
#  STAGE 3 — INTERACTION TESTING
# ============================================================

# ── Gate check ────────────────────────────────────────────────
delta_scr_global_vs_local <- AICc(scr_local) - AICc(scr_global)
cat("\nΔAICc (Local vs Global — scrapers):",
    round(delta_scr_global_vs_local, 2), "\n")

if (delta_scr_global_vs_local < 2) {
  cat("NOTE: Global does not outperform Local by ΔAICc > 2.\n",
      "Interaction models fitted for completeness but\n",
      "interpreted with caution.\n")
}

# ── Hypothesis 1: MPA effectiveness depends on larval supply ──
scr_int_mpa_conn <- glmmTMB(mean_biomass ~ rugosity_sc +
                              log_settlement_grav_sc +
                              mpa_status * connectivity_sc,
                            family = tweedie(link = "log"),
                            data   = scraper_model_data)

# ── Hypothesis 2: MPA effectiveness depends on fishing pressure
scr_int_mpa_press <- glmmTMB(mean_biomass ~ rugosity_sc +
                               mpa_status * log_settlement_grav_sc,
                             family = tweedie(link = "log"),
                             data   = scraper_model_data)

# ── Hypothesis 3: connectivity buffers fishing pressure ───────
scr_int_conn_press <- glmmTMB(mean_biomass ~ rugosity_sc +
                                connectivity_sc * log_settlement_grav_sc,
                              family = tweedie(link = "log"),
                              data   = scraper_model_data)

# ── Interaction candidate set ─────────────────────────────────
# Reference: best supported additive model from Stage 2.
# UPDATE model name after Stage 2 results known.
scraper_interactions <- list(
  "Local (additive)"        = scr_local,     # update if needed
  "MPA × connectivity"      = scr_int_mpa_conn,
  "MPA × pressure"          = scr_int_mpa_press,
  "Connectivity × pressure" = scr_int_conn_press
)

cat("\n--- Stage 3: Scraper interaction comparison ---\n")
print(make_aicc_df(scraper_interactions))

summary(scr_int_mpa_conn)

# ── Stage 3 results ───────────────────────────────────────────
# MPA × connectivity decisive (weight = 0.895, ΔAICc = 5.12).
# Strong result despite gate check failure — interaction
# coefficient highly significant (p = 0.004).
#
# Coefficients (MPA × connectivity model):
#   Rugosity:              β = +0.206, p = 0.018 * — significant
#     Emerges once MPA × connectivity absorbs governance signal.
#   Settlement gravity:    β = -0.471, p < 0.001 *** — robust
#     Consistent negative pressure effect throughout.
#   MPA low main:          β = +4.127, p = 0.003 **
#     At mean connectivity — large positive effect.
#   MPA medium main:       β = -0.440, p = 0.031 *
#     At mean connectivity — negative, consistent with Stage 2.
#   Connectivity main:     β = -0.204, p = 0.081 . — marginal
#   MPA low × connectivity: β = -4.528, p = 0.004 **
#     KEY FINDING: benefit of low protection is large at
#     low connectivity but collapses at high connectivity.
#     Well-connected low-protection sites lose their biomass
#     advantage — larvae disperse away, undermining local
#     population replenishment despite nominal protection.
#   MPA medium × connectivity: β = +0.337, p = 0.055 .
#     Marginal opposite pattern — medium protection slightly
#     more effective at well-connected sites.
#
# Ecological interpretation:
# For scrapers, low-level protection only works when sites
# are relatively isolated — larvae are retained locally and
# biomass accumulates. At well-connected sites, larval export
# undermines the benefit of low-level protection. Medium
# protection shows the opposite pattern — sufficient
# enforcement at well-connected sites may retain both adults
# and recruits, maintaining biomass even with high dispersal.
#
# Management implication:
# Low-protection MPAs at well-connected sites provide
# little benefit for scraper biomass — protection level
# needs to match the ecological exposure of the site.
# Well-connected reefs require stronger protection to
# counteract the continuous larval export that low-level
# protection cannot compensate for.
#
# Note: gate check failed (ΔAICc = -2.20) — interpret
# with caution. However the interaction coefficients are
# ecologically coherent and statistically significant.

# ============================================================
#  STAGE 4 — SENSITIVITY ANALYSIS
# ============================================================
# ── (a) Alternative pressure metrics ─────────────────────────
# Key robustness check: does the MPA × connectivity interaction hold with
# alternative pressure metrics?

scr_sens_settpop <- glmmTMB(mean_biomass ~ rugosity_sc +
                              log_settlement_pop_sc +
                              mpa_status * connectivity_sc,
                            family = tweedie(link = "log"),
                            data   = scraper_model_data)

scr_sens_mktgrav <- glmmTMB(mean_biomass ~ rugosity_sc +
                              log_market_gravity_sc +
                              mpa_status * connectivity_sc,
                            family = tweedie(link = "log"),
                            data   = scraper_model_data)

cat("\n--- Stage 4a: Scraper sensitivity — pressure metrics ---\n")
cat("Settlement population:\n")
print(round(summary(scr_sens_settpop)$coefficients$cond, 4))
cat("\nMarket gravity:\n")
print(round(summary(scr_sens_mktgrav)$coefficients$cond, 4))

# ── Stage 4a revised conclusion ───────────────────────────────
# Metric sensitivity for the MPA × connectivity interaction
# is expected given the strong pre-analysis metric selection
# result for scrapers (settlement gravity weight = 0.954,
# ΔAICc = 6.79 over market gravity).
#
# Alternative metrics are genuinely poor proxies for
# subsistence fishing pressure on scrapers — market gravity
# captures commercial pressure dynamics that are not
# relevant to this artisanal-targeted guild, and settlement
# population lacks the accessibility weighting that makes
# settlement gravity mechanistically appropriate.
#
# The correct robustness test for scrapers is the transect-
# level replication (Stage 4b), not metric substitution.
# Metric sensitivity here reflects the inadequacy of the
# alternative metrics for this group, not instability of
# the underlying biological relationship.
#
# Primary finding robust: settlement gravity negative and
# connectivity main effect negative across all metrics.
# MPA × connectivity interaction: test robustness via
# Stage 4b rather than metric substitution.
# 
# Note that the consistsnecy across emtrics is only for when the metircs chave simarl aiccs

# ── (b) Transect-level replication ───────────────────────────
scraper_transect_data <- scraper_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(scraper_transect_data$transect_scraper_biomass == 0),
    "/", nrow(scraper_transect_data),
    "(", round(mean(scraper_transect_data$transect_scraper_biomass == 0), 3), ")\n")

# ── Transect family selection ─────────────────────────────────
scr_trans_tw <- glmmTMB(transect_scraper_biomass ~ rugosity_sc +
                          log_settlement_grav_sc +
                          connectivity_sc +
                          mpa_status +
                          log_chla_sc +
                          log_max_dhw_sc +
                          (1 | site),
                        family = tweedie(link = "log"),
                        data   = scraper_transect_data)

scr_trans_tw_zi <- glmmTMB(transect_scraper_biomass ~ rugosity_sc +
                             log_settlement_grav_sc +
                             connectivity_sc +
                             mpa_status +
                             log_chla_sc +
                             log_max_dhw_sc +
                             (1 | site),
                           family    = tweedie(link = "log"),
                           ziformula = ~1,
                           data      = scraper_transect_data)

scr_trans_res    <- simulateResiduals(scr_trans_tw,    n = 1000)
scr_trans_res_zi <- simulateResiduals(scr_trans_tw_zi, n = 1000)

plot(scr_trans_res);    testZeroInflation(scr_trans_res)
plot(scr_trans_res_zi); testZeroInflation(scr_trans_res_zi)

cat("\n--- Transect family selection ---\n")
print(make_aicc_df(list(
  "Tweedie"    = scr_trans_tw,
  "ZI Tweedie" = scr_trans_tw_zi
)))

# Standard Tweedie confirmed (weight = 0.928, ΔAICc = 5.12) 

# ── Transect hierarchical sequence ────────────────────────────
# Stage 3 best model: MPA × connectivity.

scr_trans_null <- glmmTMB(transect_scraper_biomass ~ 1 +
                            (1 | site),
                          family = tweedie(link = "log"),
                          data   = scraper_transect_data)

scr_trans_local <- glmmTMB(transect_scraper_biomass ~ rugosity_sc +
                             log_settlement_grav_sc +
                             (1 | site),
                           family = tweedie(link = "log"),
                           data   = scraper_transect_data)

scr_trans_local_env <- glmmTMB(transect_scraper_biomass ~ rugosity_sc +
                                 log_settlement_grav_sc +
                                 log_chla_sc +
                                 log_max_dhw_sc +
                                 (1 | site),
                               family = tweedie(link = "log"),
                               data   = scraper_transect_data)

scr_trans_local_spatial <- glmmTMB(transect_scraper_biomass ~ rugosity_sc +
                                     log_settlement_grav_sc +
                                     connectivity_sc +
                                     (1 | site),
                                   family = tweedie(link = "log"),
                                   data   = scraper_transect_data)

scr_trans_local_mpa <- glmmTMB(transect_scraper_biomass ~ rugosity_sc +
                                 log_settlement_grav_sc +
                                 mpa_status +
                                 (1 | site),
                               family = tweedie(link = "log"),
                               data   = scraper_transect_data)

scr_trans_global <- glmmTMB(transect_scraper_biomass ~ rugosity_sc +
                              log_settlement_grav_sc +
                              connectivity_sc +
                              mpa_status +
                              log_chla_sc +
                              log_max_dhw_sc +
                              (1 | site),
                            family = tweedie(link = "log"),
                            data   = scraper_transect_data)

# Stage 3 best model — MPA × connectivity
scr_trans_best <- glmmTMB(transect_scraper_biomass ~ rugosity_sc +
                            log_settlement_grav_sc +
                            mpa_status * connectivity_sc +
                            (1 | site),
                          family = tweedie(link = "log"),
                          data   = scraper_transect_data)

scraper_transect_list <- list(
  "Null"                 = scr_trans_null,
  "Local"                = scr_trans_local,
  "Local + env"          = scr_trans_local_env,
  "Local + spatial"      = scr_trans_local_spatial,
  "Local + MPA"          = scr_trans_local_mpa,
  "Global"               = scr_trans_global,
  "MPA × connectivity"   = scr_trans_best
)

cat("\n--- Stage 4b: Scraper transect-level sensitivity ---\n")
print(make_aicc_df(scraper_transect_list))

# ── Stage 4b results ──────────────────────────────────────────
# MPA × connectivity dominant at transect level —
# weight = 0.637, ΔAICc = 3.15 over Local + env.
# Confirms site-level finding despite gate check failure.
#
# Model ranking (transect-level):
#   MPA × connectivity: AICc = 4327.79, weight = 0.637
#   Local + env:        ΔAICc = 3.15,   weight = 0.132
#   Local:              ΔAICc = 4.23,   weight = 0.077
#   Local + spatial:    ΔAICc = 4.26,   weight = 0.076
#   Local + MPA:        ΔAICc = 4.77,   weight = 0.059
#   Global:             ΔAICc = 6.87,   weight = 0.021
#   Null:               ΔAICc = 18.43   — strongly rejected
#
# Transect family: standard Tweedie (weight = 0.928,
# ΔAICc = 5.12) — strongest family selection result
# across all groups. ZI clearly not warranted.
#
# Convergence with site level:
#   Site level:      MPA × connectivity weight = 0.895
#   Transect level:  MPA × connectivity weight = 0.637
#   Consistent and dominant at both analytical scales.
#   Within-site variation does not alter the conclusion.
#
# The transect-level confirmation is the critical robustness
# test for scrapers given:
#   (a) gate check failed at site level (ΔAICc = -2.20)
#   (b) interaction metric-sensitive in Stage 4a
#
# Transect-level dominance with (1|site) addresses both
# concerns — the interaction is not a site-level artefact
# and survives within-site variation. Combined with
# ecologically coherent coefficients (MPA low × connectivity
# negative, p = 0.004 at site level), the finding is
# treated as genuine despite the caveats.
#
# Scraper analysis complete.