# ============================================================
#  PISCIVORE FISH BIOMASS — SITE-LEVEL ANALYSIS
#
#  Piscivores are fish-eating predators occupying high trophic
#  positions on coral reefs. They are among the most heavily
#  targeted functional groups by artisanal and recreational
#  fishing, making them sensitive indicators of exploitation.
#
#  Analytical framework mirrors total biomass (four stages):
#
#  STAGE 1 — Variance partitioning
#             McFadden pseudo-R² + LRT (Tweedie expected).
#             Formal varpart() if Gaussian log adequate.
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
#    - Site-level zeros expected (~17% in original analysis).
#    - ~47% transect-level zeros — Tweedie required at
#      transect level regardless of site-level family.
#    - Original analysis found rugosity + market gravity
#      as top predictors — expect rugosity to hold but
#      test whether MPA and connectivity add value.
#    - Piscivores are heavily targeted — MPA × pressure
#      interaction ecologically motivated as per browsers.
#
#  Study design:
#    243 transects, 54 sites, 4 countries.
#    Minimum 3 transects per site retained.
# ============================================================

# ── SOURCE SHARED DATA PREPARATION ───────────────────────────
# Loads: fish_data, scaled_predictors, final_predictors,
# total_transects, transect_model_data, total_model_data,
# make_aicc_df(), plot_effect(), and all packages.
source(here::here("data_preparation.R"))


# ── AGGREGATE PISCIVORE TRANSECT DATA ────────────────────────
pisc_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_pisc_biomass = sum(
      ifelse(trophic_group == "piscivores", tot_wt_g, 0),
      na.rm = TRUE
    ),
    transect_pisc_count = sum(
      ifelse(trophic_group == "piscivores", number, 0),
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

cat("Piscivore transects:", nrow(pisc_transects), "\n")
cat("Sites:",               n_distinct(pisc_transects$site), "\n")
cat("Countries:",           n_distinct(pisc_transects$country), "\n")


# ── SITE-LEVEL AGGREGATION ────────────────────────────────────
pisc_site_data <- pisc_transects %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass = mean(transect_pisc_biomass, na.rm = TRUE),
    n_transects  = n(),
    .groups      = "drop"
  ) %>%
  mutate(
    site    = as.factor(site),
    country = as.factor(country)
  )

# ── Zero check ────────────────────────────────────────────────
zeros_site <- mean(pisc_site_data$mean_biomass == 0, na.rm = TRUE)
cat("Site-level zero proportion:", round(zeros_site, 3), "\n")
cat("Site-level zeros (n):",       sum(pisc_site_data$mean_biomass == 0), "\n")
# Zeros expected — piscivores absent from some sites.
# If bimodal log distribution → Tweedie confirmed.

# ── Distribution plots ────────────────────────────────────────
pisc_site_data <- pisc_site_data %>%
  mutate(log_mean_biomass = log(mean_biomass + 0.01))

p_raw <- ggplot(pisc_site_data, aes(x = mean_biomass)) +
  geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
  labs(x = "Mean piscivore biomass per site (g)",
       y = "Frequency", title = "Raw") + theme_bw()

p_log <- ggplot(pisc_site_data, aes(x = log_mean_biomass)) +
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
# Bimodal — 6 zero sites form isolated cluster at -4.6,
# separated from main distribution by clear gap.
# Tweedie selected — consistent with family selection
# diagnostics (ZI p = 1, ΔAICc = 3.17). varpart() not available.

# ── Build piscivore site-level model dataset ──────────────────
pisc_model_data <- pisc_transects %>%
  left_join(scaled_predictors, by = "site") %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_pisc_biomass, na.rm = TRUE),
    log_mean_biomass       = log(mean(transect_pisc_biomass,
                                      na.rm = TRUE) + 0.01),
    n_transects            = n(),
    rugosity_sc            = first(rugosity_sc),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_settlement_pop_sc = first(log_settlement_pop_sc),
    connectivity_sc        = first(connectivity_sc),
    mpa_status             = first(mpa_status),
    log_chla_sc            = first(log_chla_sc),
    log_max_dhw_sc         = first(log_max_dhw_sc),
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

cat("\nPiscivore model data:", nrow(pisc_model_data), "sites\n")
cat("Site-level zeros:",      sum(pisc_model_data$mean_biomass == 0), "\n")

# ── Verify no NAs ─────────────────────────────────────────────
pisc_model_data %>%
  dplyr::select(site, rugosity_sc, log_market_gravity_sc ,
                connectivity_sc, mpa_status,
                log_chla_sc, log_max_dhw_sc) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

# ── Check zeros and response ──────────────────────────────────
cat("Zeros in mean_biomass:",    sum(pisc_model_data$mean_biomass == 0), "\n")
cat("-Inf in log_mean_biomass:", sum(is.infinite(pisc_model_data$log_mean_biomass)), "\n")
cat("NAs in log_mean_biomass:",  sum(is.na(pisc_model_data$log_mean_biomass)), "\n")
cat("\nResponse summary:\n")
print(summary(pisc_model_data$log_mean_biomass))


# ============================================================
#  PRE-ANALYSIS: HUMAN PRESSURE METRIC SELECTION
# ============================================================

pisc_press_settgrav <- glmmTMB(mean_biomass ~ log_settlement_grav_sc,
                               family = tweedie(link = "log"),
                               data   = pisc_model_data)
pisc_press_settpop  <- glmmTMB(mean_biomass ~ log_settlement_pop_sc,
                               family = tweedie(link = "log"),
                               data   = pisc_model_data)
pisc_press_mktgrav  <- glmmTMB(mean_biomass ~ log_market_gravity_sc,
                               family = tweedie(link = "log"),
                               data   = pisc_model_data)

cat("\n--- Pre-analysis: piscivore pressure metric selection ---\n")
print(make_aicc_df(list(
  "Settlement gravity" = pisc_press_settgrav,
  "Settlement pop."    = pisc_press_settpop,
  "Market gravity"     = pisc_press_mktgrav
)))


# ── Pressure metric decision ──────────────────────────────────
# Market gravity clearly preferred — weight = 0.815, ΔAICc = 2.96.
# First functional group where one metric decisively dominates
# (all other groups had ΔAICc < 1 across metrics).
#
# Results:
#   Market gravity:     AICc = 795.95, weight = 0.815
#   Settlement gravity: ΔAICc = 2.96,  weight = 0.185
#
# Ecological interpretation: piscivores are large-bodied,
# high-value commercial species specifically targeted by
# market-oriented fishing. Market gravity — which captures
# both market size and travel time — is a more appropriate
# proxy for commercial fishing pressure than settlement
# gravity, which reflects subsistence-scale proximity.
# This is consistent with the original analysis.
#
# log_market_gravity_sc used as primary pressure predictor
# throughout all piscivore analyses. Settlement gravity
# retained for Stage 4a sensitivity check only.

# ============================================================
#  MODEL FAMILY SELECTION
# ============================================================

# ── Tweedie (global model) ────────────────────────────────────
pisc_tw_global <- glmmTMB(mean_biomass ~ rugosity_sc +
                            log_market_gravity_sc  +
                            connectivity_sc +
                            mpa_status +
                            log_chla_sc +
                            log_max_dhw_sc,
                          family = tweedie(link = "log"),
                          data   = pisc_model_data)

pisc_tw_res <- simulateResiduals(pisc_tw_global, n = 1000)
plot(pisc_tw_res)
testZeroInflation(pisc_tw_res)
testDispersion(pisc_tw_res)

# ── ZI Tweedie ────────────────────────────────────────────────
pisc_tw_zi_global <- glmmTMB(mean_biomass ~ rugosity_sc +
                               log_market_gravity_sc  +
                               connectivity_sc +
                               mpa_status +
                               log_chla_sc +
                               log_max_dhw_sc,
                             family    = tweedie(link = "log"),
                             ziformula = ~1,
                             data      = pisc_model_data)

pisc_tw_zi_res <- simulateResiduals(pisc_tw_zi_global, n = 1000)
plot(pisc_tw_zi_res)
testZeroInflation(pisc_tw_zi_res)

cat("\n--- Family selection: Tweedie vs ZI Tweedie ---\n")
print(make_aicc_df(list(
  "Tweedie"    = pisc_tw_global,
  "ZI Tweedie" = pisc_tw_zi_global
)))

# ── Family selection decision ─────────────────────────────────
# Standard Tweedie(log link) selected.
#
# DHARMa diagnostics (global model, n = 1000):
#   Dispersion:     p = 0.706, ratio = 1.049 — excellent
#   Zero inflation: p = 1.000, ratio = 0.996 — not ZI
#
# ZI Tweedie not supported:
#   AICc improves by only 3.17 — below the >2 threshold
#   ZI test non-significant (p = 1) for both models
#   Neither condition for ZI adoption met
#
# Bimodal log distribution confirmed from plots — 6 zero
# sites form an isolated cluster, consistent with Tweedie
# being the appropriate family. Tweedie handles these zeros
# natively without a log + constant approximation.

# # ============================================================
# #  RANDOM EFFECT STRUCTURE
# # ============================================================
# 
# pisc_re_null <- glmmTMB(mean_biomass ~ rugosity_sc +
#                           log_market_gravity_sc  +
#                           log_chla_sc +
#                           log_max_dhw_sc,
#                         family = tweedie(link = "log"),
#                         data   = pisc_model_data)
# 
# pisc_re_country <- glmmTMB(mean_biomass ~ rugosity_sc +
#                              log_market_gravity_sc  +
#                              log_chla_sc +
#                              log_max_dhw_sc +
#                              (1 | country),
#                            family = tweedie(link = "log"),
#                            data   = pisc_model_data)
# 
# cat("\n--- Piscivore RE structure comparison ---\n")
# print(make_aicc_df(list(
#   "No RE"         = pisc_re_null,
#   "(1 | country)" = pisc_re_country
# )))

# Not really a signficiant difference but either way no RE kept

# ============================================================
#  STAGE 1 — VARIANCE PARTITIONING
#
#  Tweedie expected — McFadden pseudo-R² + LRT as per
#  browsers, corallivores, and excavators.
#  Three process groups: local, spatial, environmental.
#  MPA excluded — governance variable, tested in Stage 2.
# ============================================================

p_vp_null <- glmmTMB(mean_biomass ~ 1,
                     family = tweedie(link = "log"),
                     data   = pisc_model_data)

p_vp_all <- glmmTMB(mean_biomass ~ rugosity_sc +
                      log_market_gravity_sc  +
                      connectivity_sc +
                      log_chla_sc +
                      log_max_dhw_sc,
                    family = tweedie(link = "log"),
                    data   = pisc_model_data)

p_vp_no_local <- glmmTMB(mean_biomass ~ connectivity_sc +
                           log_chla_sc +
                           log_max_dhw_sc,
                         family = tweedie(link = "log"),
                         data   = pisc_model_data)

p_vp_no_spatial <- glmmTMB(mean_biomass ~ rugosity_sc +
                             log_market_gravity_sc  +
                             log_chla_sc +
                             log_max_dhw_sc,
                           family = tweedie(link = "log"),
                           data   = pisc_model_data)

p_vp_no_environ <- glmmTMB(mean_biomass ~ rugosity_sc +
                             log_market_gravity_sc  +
                             connectivity_sc,
                           family = tweedie(link = "log"),
                           data   = pisc_model_data)

# ── McFadden pseudo-R² ────────────────────────────────────────
null_loglik_p <- as.numeric(logLik(p_vp_null))

get_mcfadden_p <- function(model) {
  round(1 - (as.numeric(logLik(model)) / null_loglik_p), 3)
}

r2_p_all    <- get_mcfadden_p(p_vp_all)
r2_p_noloc  <- get_mcfadden_p(p_vp_no_local)
r2_p_nospat <- get_mcfadden_p(p_vp_no_spatial)
r2_p_noenv  <- get_mcfadden_p(p_vp_no_environ)

unique_p_local   <- round(r2_p_all - r2_p_noloc,  3)
unique_p_spatial <- round(r2_p_all - r2_p_nospat, 3)
unique_p_environ <- round(r2_p_all - r2_p_noenv,  3)

cat("\n--- Stage 1: Piscivore pseudo-R² variance decomposition ---\n")
cat("Total McFadden R² (all groups):", r2_p_all,        "\n")
cat("Unique local:                  ", unique_p_local,   "\n")
cat("Unique spatial:                ", unique_p_spatial, "\n")
cat("Unique environmental:          ", unique_p_environ, "\n")

# ── LRT significance tests ────────────────────────────────────
cat("\nLRT — local unique fraction:\n")
print(anova(p_vp_no_local, p_vp_all))

cat("\nLRT — spatial unique fraction:\n")
print(anova(p_vp_no_spatial, p_vp_all))

cat("\nLRT — environmental unique fraction:\n")
print(anova(p_vp_no_environ, p_vp_all))

# ── Stage 1 results ───────────────────────────────────────────
# Unique local:         McF R² = 0.003, p = 0.316 — not significant
# Unique spatial:       McF R² = 0.003, p = 0.133 — not significant
# Unique environmental: McF R² = 0.000, p = 0.927 — not significant
# Total McFadden R²:    0.008
#
# Complete null result — no process group explains significant
# unique variance in piscivore biomass.
#
# This is the starkest Stage 1 result across all functional
# groups. Even the local fraction (rugosity + market gravity),
# which was significant for total biomass (p = 0.001),
# browsers (p = 0.001), and excavators (p = 0.042), and
# marginal for grazers (p = 0.073), is completely
# non-significant for piscivores (p = 0.316).
#
# Compare across groups:
#   Total biomass: local***,  spatial n.s., env .
#   Browsers:      local***,  spatial *,    env n.s.
#   Corallivores:  local n.s., spatial .,   env n.s.
#   Grazers:       local .,   spatial *,    env n.s.
#   Excavators:    local *,   spatial n.s., env n.s.
#   Piscivores:    local n.s., spatial n.s., env n.s.
#
# Piscivore biomass is not structured by any of the measured
# process groups at this spatial scale. This is consistent
# with piscivores being determined by historical exploitation,
# body size dynamics, and large home ranges rather than
# contemporary site-level environmental gradients.
#
# Note: spatial fraction (p = 0.133) is the closest to
# significance — connectivity may play a role but signal
# is too weak to detect at n = 54.

# ============================================================
#  STAGE 2 — HIERARCHICAL MODEL COMPARISON
#
#  UPDATE family after family selection decision above.
#  Tweedie used as default throughout.
# ============================================================

# ── Null model ────────────────────────────────────────────────
p_null <- glmmTMB(mean_biomass ~ 1,
                  family = tweedie(link = "log"),
                  data   = pisc_model_data)

# ── Local ecological baseline ─────────────────────────────────
p_local <- glmmTMB(mean_biomass ~ rugosity_sc +
                     log_market_gravity_sc ,
                   family = tweedie(link = "log"),
                   data   = pisc_model_data)

# ── Local + chla ──────────────────────────────────────────────
p_local_chla <- glmmTMB(mean_biomass ~ rugosity_sc +
                          log_market_gravity_sc  +
                          log_chla_sc,
                        family = tweedie(link = "log"),
                        data   = pisc_model_data)

# ── Local + DHW ───────────────────────────────────────────────
p_local_dhw <- glmmTMB(mean_biomass ~ rugosity_sc +
                         log_market_gravity_sc  +
                         log_max_dhw_sc,
                       family = tweedie(link = "log"),
                       data   = pisc_model_data)

# ── Local + environmental context ────────────────────────────
p_local_env <- glmmTMB(mean_biomass ~ rugosity_sc +
                         log_market_gravity_sc  +
                         log_chla_sc +
                         log_max_dhw_sc,
                       family = tweedie(link = "log"),
                       data   = pisc_model_data)

# ── Local + spatial ───────────────────────────────────────────
p_local_spatial <- glmmTMB(mean_biomass ~ rugosity_sc +
                             log_market_gravity_sc  +
                             connectivity_sc,
                           family = tweedie(link = "log"),
                           data   = pisc_model_data)

# ── Local + MPA (governance test) ────────────────────────────
p_local_mpa <- glmmTMB(mean_biomass ~ rugosity_sc +
                         log_market_gravity_sc  +
                         mpa_status,
                       family = tweedie(link = "log"),
                       data   = pisc_model_data)

# ── Global model ──────────────────────────────────────────────
p_global <- glmmTMB(mean_biomass ~ rugosity_sc +
                      log_market_gravity_sc  +
                      connectivity_sc +
                      mpa_status +
                      log_chla_sc +
                      log_max_dhw_sc,
                    family = tweedie(link = "log"),
                    data   = pisc_model_data)

# ── Model list ────────────────────────────────────────────────
pisc_model_list <- list(
  "Null"            = p_null,
  "Local"           = p_local,
  "Local + chla"    = p_local_chla,
  "Local + DHW"     = p_local_dhw,
  "Local + env"     = p_local_env,
  "Local + MPA"     = p_local_mpa,
  "Local + spatial" = p_local_spatial,
  "Global"          = p_global
)

# ── AICc ranked table ─────────────────────────────────────────
cat("\n--- Stage 2: Piscivore model comparison (AICc ranked) ---\n")
print(make_aicc_df(pisc_model_list))

# ── McFadden R² increments over local baseline ────────────────
cat("\n--- Stage 2: Piscivore variance explained ---\n")
local_loglik_p <- as.numeric(logLik(p_local))
null_loglik_p2 <- as.numeric(logLik(p_null))

pisc_model_list %>%
  imap_dfr(~ tibble(
    Model  = .y,
    McF_R2 = round(1 - (as.numeric(logLik(.x)) /
                          null_loglik_p2), 3)
  )) %>%
  mutate(
    Delta_R2 = round(McF_R2 - (1 - (local_loglik_p /
                                      null_loglik_p2)), 3),
    Delta_R2 = ifelse(Model %in% c("Null", "Local"), NA, Delta_R2)
  ) %>%
  arrange(Model == "Null",
          desc(Model == "Local"),
          desc(McF_R2)) %>%
  print()

# ── Stage 2 results ───────────────────────────────────────────
# Null best (weight = 0.273) — Local + MPA essentially tied
# (ΔAICc = 0.10, weight = 0.260). Four models within ΔAICc < 2
# with combined weight = 0.771 — high model selection uncertainty.
#
# AICc ranking:
#   Null:            AICc = 796.61, weight = 0.273
#   Local + MPA:     ΔAICc = 0.10,  weight = 0.260
#   Local:           ΔAICc = 0.92,  weight = 0.173
#   Local + spatial: ΔAICc = 1.00,  weight = 0.166
#   All others:      ΔAICc > 3.11  — not supported
#
# R² increments over Local baseline (McF R² = 0.005):
#   Local + MPA:     +0.007 — marginal
#   Local + spatial: +0.003 — negligible
#   Chla, DHW, env:  +0.000 — nothing
#   Global:          +0.009 — driven by MPA and spatial
#
# Model sequence is appropriate — the eight models cover the
# full range of ecologically motivated predictors. The null
# result is not an artefact of testing the wrong models.
# Neither habitat, pressure, connectivity, MPA, nor
# environmental context explains piscivore biomass reliably.
#
# Local + MPA marginally competitive — MPA adds a small
# positive increment (+0.007) but this is within model
# selection noise given the null is best overall.
# Examine summary(p_local_mpa) to check MPA coefficient
# direction before proceeding to Stage 3.

summary(p_local_mpa)
summary(p_local)

# ── Coefficient summary: Local and Local + MPA ────────────────
# Local model:
#   Rugosity:       β = +0.118, p = 0.354 — not significant
#   Market gravity: β = +0.298, p = 0.047 * — significant,
#     POSITIVE direction — counterintuitive. Higher market
#     gravity associated with higher piscivore biomass.
#     Likely reflects spatial confounding: high-gravity sites
#     near urban centres may coincide with historically
#     productive or better-monitored reefs rather than
#     indicating a true positive fishing effect.
#     Interpret with caution throughout.
#
# Local + MPA model:
#   Rugosity:       β = +0.141, p = 0.245 — not significant
#   Market gravity: β = +0.306, p = 0.034 * — significant,
#     positive direction, consistent with Local model
#   MPA low:        β = +0.081, p = 0.852 — not significant
#   MPA medium:     β = +0.707, p = 0.015 * — significant
#     Medium protection sites have substantially higher
#     piscivore biomass than unprotected sites. Most
#     ecologically coherent signal in Stage 2 — consistent
#     with piscivores being heavily targeted by fishing and
#     benefiting from protection.
#
# NOTE: positive market gravity coefficient is inconsistent
# with the expected negative fishing pressure effect.
# This was flagged in the original analysis as likely
# spatial confounding. Connectivity × pressure interaction
# (Stage 3, weight = 0.894) may help resolve this —
# the market gravity effect may be conditional on
# connectivity rather than a direct pressure effect.

# ============================================================
#  STAGE 3 — INTERACTION TESTING
# ============================================================

# ── Gate check ────────────────────────────────────────────────
delta_p_global_vs_local <- AICc(p_local) - AICc(p_global)
cat("\nΔAICc (Local vs Global — piscivores):",
    round(delta_p_global_vs_local, 2), "\n")

if (delta_p_global_vs_local < 2) {
  cat("NOTE: Global does not outperform Local by ΔAICc > 2.\n",
      "Interaction models fitted for completeness but\n",
      "interpreted with caution.\n")
}

# ── Hypothesis 1: MPA effectiveness depends on larval supply ──
p_int_mpa_conn <- glmmTMB(mean_biomass ~ rugosity_sc +
                            log_market_gravity_sc  +
                            mpa_status * connectivity_sc,
                          family = tweedie(link = "log"),
                          data   = pisc_model_data)

# ── Hypothesis 2: MPA effectiveness depends on fishing pressure
p_int_mpa_press <- glmmTMB(mean_biomass ~ rugosity_sc +
                             mpa_status * log_market_gravity_sc ,
                           family = tweedie(link = "log"),
                           data   = pisc_model_data)

# ── Hypothesis 3: connectivity buffers fishing pressure ───────
p_int_conn_press <- glmmTMB(mean_biomass ~ rugosity_sc +
                              connectivity_sc * log_market_gravity_sc ,
                            family = tweedie(link = "log"),
                            data   = pisc_model_data)

# ── Interaction candidate set ─────────────────────────────────
# Reference: best supported additive model from Stage 2.
# UPDATE model name after Stage 2 results known.
pisc_interactions <- list(
  "Local (additive)"        = p_local,      # update if needed
  "MPA × connectivity"      = p_int_mpa_conn,
  "MPA × pressure"          = p_int_mpa_press,
  "Connectivity × pressure" = p_int_conn_press
)

cat("\n--- Stage 3: Piscivore interaction comparison ---\n")
print(make_aicc_df(pisc_interactions))

summary(p_int_conn_press)

# ── Stage 3 results ───────────────────────────────────────────
# Connectivity × pressure decisive (weight = 0.894, ΔAICc = 5.02
# over Local additive). Strongest interaction result after
# browsers (weight = 0.965).
#
# Coefficients (connectivity × pressure model):
#   Rugosity:              β = +0.239, p = 0.066 . — marginal
#   Connectivity:          β = +0.264, p = 0.047 * — significant
#   Market gravity:        β = +0.405, p = 0.007 ** — significant
#   Conn × market gravity: β = -0.509, p = 0.004 ** — significant
#
# Interaction interpretation:
# The effect of market gravity on piscivore biomass depends
# strongly on connectivity. At low-connectivity sites, market
# gravity is positively associated with biomass — likely
# reflecting spatial confounding (accessible sites near
# markets may also be historically productive). At high-
# connectivity sites, this relationship reverses — well-
# connected sites with high market gravity have lower
# piscivore biomass, consistent with the expected negative
# fishing pressure effect.
#
# This resolves the counterintuitive positive market gravity
# direction observed in Stage 2. The positive main effect
# was masking a conditional relationship — market gravity
# only appears positive when connectivity is not accounted
# for. At high connectivity, fishing pressure operates as
# expected: higher market access reduces piscivore biomass.
#
# Ecological interpretation:
# Well-connected reefs are more accessible to fishing
# vessels — connectivity correlates with open-water
# exposure and network centrality, which also determines
# how easily fishers can reach a site. At high-connectivity
# sites, market gravity therefore translates more directly
# into fishing mortality, reducing piscivore biomass.
# At isolated sites, market gravity is a weaker proxy for
# actual exploitation because distance and accessibility
# constrain fishing effort regardless of market size.
#
# Management implication:
# Well-connected reefs with high market gravity are the
# most vulnerable sites for piscivore depletion —
# both ecologically exposed (high larval export, open
# network position) and economically accessible. These
# sites should be priority targets for protection.
# This complements the browser finding (MPA × pressure)
# and strengthens the case for connectivity-informed
# MPA placement.
#
# MPA × pressure:     ΔAICc = 6.96 — not supported
# MPA × connectivity: ΔAICc = 10.19 — not supported

# ============================================================
#  STAGE 4 — SENSITIVITY ANALYSIS
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
# Mirrors connectivity × pressure — best Stage 3 model.
# Only pressure metric swapped. Key robustness check:
# does the interaction hold with settlement gravity?

p_sens_settgrav <- glmmTMB(mean_biomass ~ rugosity_sc +
                             connectivity_sc * log_settlement_grav_sc,
                           family = tweedie(link = "log"),
                           data   = pisc_model_data)

cat("\n--- Stage 4a: Piscivore sensitivity — pressure metrics ---\n")
cat("Settlement gravity:\n")
print(round(summary(p_sens_settgrav)$coefficients$cond, 4))

# ── Stage 4a results ──────────────────────────────────────────
# Connectivity × pressure interaction metric-sensitive:
#
#   Market gravity (primary):
#     Connectivity × pressure: β = -0.509, p = 0.004 **
#     Market gravity main:     β = +0.405, p = 0.007 **
#     Connectivity main:       β = +0.264, p = 0.047 *
#
#   Settlement gravity (sensitivity):
#     Connectivity × pressure: β = -0.385, p = 0.082 .
#     Settlement gravity main: β = +0.108, p = 0.570 n.s.
#     Connectivity main:       β = +0.214, p = 0.131 n.s.
#
# Negative interaction direction consistent across both
# metrics — the pattern holds qualitatively. However
# significance is lost with settlement gravity, and the
# settlement gravity main effect is non-significant
# throughout (p = 0.570).
#
# This metric-sensitivity is expected given market gravity
# was clearly the better-performing metric for piscivores
# (weight = 0.815, ΔAICc = 2.96). Settlement gravity is a
# weaker proxy for commercial fishing pressure on piscivores
# — it captures subsistence-scale proximity rather than
# market-oriented exploitation.
#
# Conclusion: connectivity × pressure interaction is
# supported with the ecologically appropriate metric
# (market gravity) but not fully robust to metric
# substitution. Interpret as a suggestive finding
# requiring replication with direct fishing effort data.
# The direction is consistent — this is not a spurious
# result, but it is metric-dependent.


# ── (b) Transect-level replication ───────────────────────────
# ~47% zeros at transect level — Tweedie required.
# ZI Tweedie tested given moderate zero proportion.

pisc_transect_data <- pisc_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(pisc_transect_data$transect_pisc_biomass == 0),
    "/", nrow(pisc_transect_data),
    "(", round(mean(pisc_transect_data$transect_pisc_biomass == 0), 3), ")\n")

# ── Transect family selection ─────────────────────────────────
p_trans_tw <- glmmTMB(transect_pisc_biomass ~ rugosity_sc +
                        log_market_gravity_sc  +
                        connectivity_sc +
                        mpa_status +
                        log_chla_sc +
                        log_max_dhw_sc +
                        (1 | site),
                      family = tweedie(link = "log"),
                      data   = pisc_transect_data)

p_trans_tw_zi <- glmmTMB(transect_pisc_biomass ~ rugosity_sc +
                           log_market_gravity_sc  +
                           connectivity_sc +
                           mpa_status +
                           log_chla_sc +
                           log_max_dhw_sc +
                           (1 | site),
                         family    = tweedie(link = "log"),
                         ziformula = ~1,
                         data      = pisc_transect_data)

p_trans_res    <- simulateResiduals(p_trans_tw,    n = 1000)
p_trans_res_zi <- simulateResiduals(p_trans_tw_zi, n = 1000)

plot(p_trans_res);    testZeroInflation(p_trans_res)
plot(p_trans_res_zi); testZeroInflation(p_trans_res_zi)

cat("\n--- Transect family selection ---\n")
print(make_aicc_df(list(
  "Tweedie"    = p_trans_tw,
  "ZI Tweedie" = p_trans_tw_zi
)))

# ── Transect family selection decision ────────────────────────
# Standard Tweedie selected (weight = 0.752, ΔAICc = 2.21).
# ZI test non-significant for both models (p = 0.742, p = 0.804).
# Despite 45.3% transect zeros, standard Tweedie handles
# these natively — no evidence of excess zeros beyond what
# the Tweedie distribution expects.
# Consistent with site-level family selection.

# ── Transect hierarchical sequence ────────────────────────────
# Stage 3 best model used — UPDATE after Stage 3.

p_trans_null <- glmmTMB(transect_pisc_biomass ~ 1 +
                          (1 | site),
                        family = tweedie(link = "log"),
                        data   = pisc_transect_data)

p_trans_local <- glmmTMB(transect_pisc_biomass ~ rugosity_sc +
                           log_market_gravity_sc  +
                           (1 | site),
                         family = tweedie(link = "log"),
                         data   = pisc_transect_data)

p_trans_local_env <- glmmTMB(transect_pisc_biomass ~ rugosity_sc +
                               log_market_gravity_sc  +
                               log_chla_sc +
                               log_max_dhw_sc +
                               (1 | site),
                             family = tweedie(link = "log"),
                             data   = pisc_transect_data)

p_trans_local_spatial <- glmmTMB(transect_pisc_biomass ~ rugosity_sc +
                                   log_market_gravity_sc  +
                                   connectivity_sc +
                                   (1 | site),
                                 family = tweedie(link = "log"),
                                 data   = pisc_transect_data)

p_trans_local_mpa <- glmmTMB(transect_pisc_biomass ~ rugosity_sc +
                               log_market_gravity_sc  +
                               mpa_status +
                               (1 | site),
                             family = tweedie(link = "log"),
                             data   = pisc_transect_data)

p_trans_global <- glmmTMB(transect_pisc_biomass ~ rugosity_sc +
                            log_market_gravity_sc  +
                            connectivity_sc +
                            mpa_status +
                            log_chla_sc +
                            log_max_dhw_sc +
                            (1 | site),
                          family = tweedie(link = "log"),
                          data   = pisc_transect_data)

# Stage 3 best model at transect level — connectivity × pressure
p_trans_best <- glmmTMB(transect_pisc_biomass ~ rugosity_sc +
                          connectivity_sc * log_market_gravity_sc +
                          (1 | site),
                        family = tweedie(link = "log"),
                        data   = pisc_transect_data)

pisc_transect_list <- list(
  "Null"                       = p_trans_null,
  "Local"                      = p_trans_local,
  "Local + env"                = p_trans_local_env,
  "Local + spatial"            = p_trans_local_spatial,
  "Local + MPA"                = p_trans_local_mpa,
  "Global"                     = p_trans_global,
  "Connectivity × pressure"    = p_trans_best
)

cat("\n--- Stage 4b: Piscivore transect-level sensitivity ---\n")
print(make_aicc_df(pisc_transect_list))

# ── Stage 4b results ──────────────────────────────────────────
# Connectivity × pressure dominant at transect level —
# weight = 0.742, ΔAICc = 4.20 over Local + MPA.
# Stronger and clearer than site-level result.
#
# Model ranking (transect-level):
#   Connectivity × pressure: AICc = 2549.73, weight = 0.742
#   Local + MPA:             ΔAICc = 4.20,   weight = 0.091
#   Local + spatial:         ΔAICc = 4.80,   weight = 0.067
#   Local:                   ΔAICc = 5.80,   weight = 0.041
#   Null:                    ΔAICc = 6.05,   weight = 0.036
#   Global:                  ΔAICc = 8.08,   weight = 0.013
#   Local + env:             ΔAICc = 8.71,   weight = 0.010
#
# Transect family selection:
#   Standard Tweedie selected (weight = 0.752, ΔAICc = 2.21).
#   ZI test non-significant (p = 0.742) — no excess zeros
#   despite 45.3% transect zeros.
#
# Convergence with site level:
#   Site level:      Connectivity × pressure weight = 0.894
#   Transect level:  Connectivity × pressure weight = 0.742
#   Consistent and decisive at both analytical scales.
#   Within-site variation does not alter the conclusion.
#
# The connectivity × pressure interaction is the most
# robust finding in the piscivore analysis — significant
# at site level (p = 0.004), dominant at transect level
# (weight = 0.742), and directionally consistent with
# settlement gravity (p = 0.082).
#
# Piscivore analysis complete.