# ============================================================
#  LARGE EXCAVATOR FISH BIOMASS — SITE-LEVEL ANALYSIS
#
#  Large excavators are large-bodied parrotfish (Scaridae)
#  that excavate coral substrate, contributing significantly
#  to bioerosion and carbonate cycling. They are among the
#  most important functional groups for reef structural
#  processes and are heavily targeted by spearfishing.
#
#  Analytical framework mirrors total biomass (four stages):
#
#  STAGE 1 — Variance partitioning
#             McFadden pseudo-R² + LRT if Tweedie selected.
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
#    - Original analysis had ~15% site-level zeros (10/64
#      sites) — check zero proportion in current dataset.
#    - Rugosity was the dominant predictor in original
#      analysis — expect this to hold.
#    - ~57.8% transect-level zeros — Tweedie required at
#      transect level regardless of site-level family.
#
#  Study design:
#    243 transects, 54 sites, 4 countries.
#    Minimum 3 transects per site retained.
# ============================================================

# ── SOURCE SHARED DATA PREPARATION ───────────────────────────
source(here::here("data_preparation.R"))


# ── AGGREGATE LARGE EXCAVATOR TRANSECT DATA ──────────────────
excavator_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_excavator_biomass = sum(
      ifelse(trophic_group == "large_excavators", tot_wt_g, 0),
      na.rm = TRUE
    ),
    transect_excavator_count = sum(
      ifelse(trophic_group == "large_excavators", number, 0),
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

cat("Excavator transects:", nrow(excavator_transects), "\n")
cat("Sites:",               n_distinct(excavator_transects$site), "\n")
cat("Countries:",           n_distinct(excavator_transects$country), "\n")


# ── SITE-LEVEL AGGREGATION ────────────────────────────────────
excavator_site_data <- excavator_transects %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass = mean(transect_excavator_biomass, na.rm = TRUE),
    n_transects  = n(),
    .groups      = "drop"
  ) %>%
  mutate(
    site    = as.factor(site),
    country = as.factor(country)
  )

# ── Zero check ────────────────────────────────────────────────
zeros_site <- mean(excavator_site_data$mean_biomass == 0, na.rm = TRUE)
cat("Site-level zero proportion:", round(zeros_site, 3), "\n")
cat("Site-level zeros (n):",       sum(excavator_site_data$mean_biomass == 0), "\n")
# Large excavators are absent from some sites — zeros expected.
# If > 5%: Tweedie required. If bimodal distribution: confirm Tweedie.

# ── Distribution plots ────────────────────────────────────────
excavator_site_data <- excavator_site_data %>%
  mutate(log_mean_biomass = log(mean_biomass + 0.01))

p_raw <- ggplot(excavator_site_data, aes(x = mean_biomass)) +
  geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
  labs(x = "Mean excavator biomass per site (g)",
       y = "Frequency", title = "Raw") + theme_bw()

p_log <- ggplot(excavator_site_data, aes(x = log_mean_biomass)) +
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
# There are a lot of zeros, which means a tweedie model will likely be needed

# ── Build full excavator site-level model dataset ─────────────
excavator_model_data <- excavator_transects %>%
  left_join(scaled_predictors, by = "site") %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_excavator_biomass, na.rm = TRUE),
    log_mean_biomass       = log(mean(transect_excavator_biomass,
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

cat("\nExcavator model data:", nrow(excavator_model_data), "sites\n")
cat("Site-level zeros:",      sum(excavator_model_data$mean_biomass == 0), "\n")

# ── Verify no NAs ─────────────────────────────────────────────
excavator_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                connectivity_sc, mpa_status,
                log_chla_sc, log_max_dhw_sc) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

# ── Check zeros and response ──────────────────────────────────
cat("Zeros in mean_biomass:",    sum(excavator_model_data$mean_biomass == 0), "\n")
cat("-Inf in log_mean_biomass:", sum(is.infinite(excavator_model_data$log_mean_biomass)), "\n")
cat("NAs in log_mean_biomass:",  sum(is.na(excavator_model_data$log_mean_biomass)), "\n")
cat("\nResponse summary:\n")
print(summary(excavator_model_data$log_mean_biomass))


# ============================================================
#  PRE-ANALYSIS: HUMAN PRESSURE METRIC SELECTION
# ============================================================

excav_press_settgrav <- glmmTMB(mean_biomass ~ log_settlement_grav_sc,
                                family = tweedie(link = "log"),
                                data   = excavator_model_data)
excav_press_settpop  <- glmmTMB(mean_biomass ~ log_settlement_pop_sc,
                                family = tweedie(link = "log"),
                                data   = excavator_model_data)
excav_press_mktgrav  <- glmmTMB(mean_biomass ~ log_market_gravity_sc,
                                family = tweedie(link = "log"),
                                data   = excavator_model_data)

cat("\n--- Pre-analysis: excavator pressure metric selection ---\n")
print(make_aicc_df(list(
  "Settlement gravity" = excav_press_settgrav,
  "Settlement pop."    = excav_press_settpop,
  "Market gravity"     = excav_press_mktgrav
)))

# ── Pressure metric decision ──────────────────────────────────
# Complete metric uncertainty — all three within ΔAICc < 0.35.
#
# Results:
#   Market gravity:     AICc = 690.62, weight = 0.366 (best)
#   Settlement pop.:    ΔAICc = 0.23,  weight = 0.326
#   Settlement gravity: ΔAICc = 0.34,  weight = 0.309
#
# Consistent with corallivores and grazers — metric selection
# indeterminate for functional groups where fishing pressure
# is not the primary driver. Settlement gravity selected for
# consistency with total biomass and browsers.
# All three retained for sensitivity analysis 

# ============================================================
#  MODEL FAMILY SELECTION
#
#  Large excavators have site-level zeros — Gaussian log
#  with offset requires checking for bimodality as per browsers.
#  Tweedie likely required. Run both diagnostics.
# ============================================================

# ── Gaussian log (global model) ───────────────────────────────
excav_lm_global <- lm(log_mean_biomass ~ rugosity_sc +
                        log_settlement_grav_sc +
                        connectivity_sc +
                        mpa_status +
                        log_chla_sc +
                        log_max_dhw_sc,
                      data = excavator_model_data)

par(mfrow = c(2, 2)); plot(excav_lm_global); par(mfrow = c(1, 1))

# ── Tweedie (global model) ────────────────────────────────────
excav_tw_global <- glmmTMB(mean_biomass ~ rugosity_sc +
                             log_settlement_grav_sc +
                             connectivity_sc +
                             mpa_status +
                             log_chla_sc +
                             log_max_dhw_sc,
                           family = tweedie(link = "log"),
                           data   = excavator_model_data)

excav_tw_res <- simulateResiduals(excav_tw_global, n = 1000)
plot(excav_tw_res)
testZeroInflation(excav_tw_res)
testDispersion(excav_tw_res)

# ── ZI Tweedie (test only if zero proportion warrants) ────────
excav_tw_zi_global <- glmmTMB(mean_biomass ~ rugosity_sc +
                                log_settlement_grav_sc +
                                connectivity_sc +
                                mpa_status +
                                log_chla_sc +
                                log_max_dhw_sc,
                              family    = tweedie(link = "log"),
                              ziformula = ~1,
                              data      = excavator_model_data)

excav_tw_zi_res <- simulateResiduals(excav_tw_zi_global, n = 1000)
plot(excav_tw_zi_res)
testZeroInflation(excav_tw_zi_res)

# ── Family selection decision ─────────────────────────────────
# Tweedie (log link) selected.
#
# Gaussian log rejected — bimodal log distribution confirmed
# (7 zero sites at log(0.01) = -4.6, separated from main
# distribution). Zero sites exert undue leverage.
#
# Tweedie diagnostics (DHARMa, n = 1000):
#   KS test:        p = 0.985 — excellent fit
#   Dispersion:     p = 0.662, ratio = 1.019 — near perfect
#   Zero inflation: p = 0.984, ratio = 0.905 — not significant
#
# ZI Tweedie rejected on convergence grounds:
#   "false convergence (8)" — model did not converge
#   reliably. Standard Tweedie sufficient given non-
#   significant zero inflation test.
#
# Proceed: glmmTMB(family = tweedie(link = "log")) on
# raw mean_biomass throughout excavator site-level analyses.
# varpart() not available — Stage 1 uses pseudo-R²
# approximation.

# # ============================================================
# #  RANDOM EFFECT STRUCTURE
# # ============================================================
# 
# excav_re_null <- glmmTMB(mean_biomass ~ rugosity_sc +
#                            log_settlement_grav_sc +
#                            log_chla_sc +
#                            log_max_dhw_sc,
#                          family = tweedie(link = "log"),
#                          data   = excavator_model_data)
# 
# excav_re_country <- glmmTMB(mean_biomass ~ rugosity_sc +
#                               log_settlement_grav_sc +
#                               log_chla_sc +
#                               log_max_dhw_sc +
#                               (1 | country),
#                             family = tweedie(link = "log"),
#                             data   = excavator_model_data)
# 
# cat("\n--- Excavator RE structure comparison ---\n")
# print(make_aicc_df(list(
#   "No RE"         = excav_re_null,
#   "(1 | country)" = excav_re_country
# )))
# 
# # No RE supported, consistent with other groups.

# ============================================================
#  STAGE 1 — VARIANCE PARTITIONING
#
#  Tweedie selected — vegan::varpart() not available.
#  McFadden pseudo-R² + LRT approximation as per browsers
#  and corallivores.
#  Three process groups: local, spatial, environmental.
#  MPA excluded — governance variable, tested in Stage 2.
# ============================================================

e_vp_null <- glmmTMB(mean_biomass ~ 1,
                     family = tweedie(link = "log"),
                     data   = excavator_model_data)

e_vp_all <- glmmTMB(mean_biomass ~ rugosity_sc +
                      log_settlement_grav_sc +
                      connectivity_sc +
                      log_chla_sc +
                      log_max_dhw_sc,
                    family = tweedie(link = "log"),
                    data   = excavator_model_data)

e_vp_no_local <- glmmTMB(mean_biomass ~ connectivity_sc +
                           log_chla_sc +
                           log_max_dhw_sc,
                         family = tweedie(link = "log"),
                         data   = excavator_model_data)

e_vp_no_spatial <- glmmTMB(mean_biomass ~ rugosity_sc +
                             log_settlement_grav_sc +
                             log_chla_sc +
                             log_max_dhw_sc,
                           family = tweedie(link = "log"),
                           data   = excavator_model_data)

e_vp_no_environ <- glmmTMB(mean_biomass ~ rugosity_sc +
                             log_settlement_grav_sc +
                             connectivity_sc,
                           family = tweedie(link = "log"),
                           data   = excavator_model_data)

# ── McFadden pseudo-R² ────────────────────────────────────────
null_loglik_e <- as.numeric(logLik(e_vp_null))

get_mcfadden_e <- function(model) {
  round(1 - (as.numeric(logLik(model)) / null_loglik_e), 3)
}

r2_e_all    <- get_mcfadden_e(e_vp_all)
r2_e_noloc  <- get_mcfadden_e(e_vp_no_local)
r2_e_nospat <- get_mcfadden_e(e_vp_no_spatial)
r2_e_noenv  <- get_mcfadden_e(e_vp_no_environ)

unique_e_local   <- round(r2_e_all - r2_e_noloc,  3)
unique_e_spatial <- round(r2_e_all - r2_e_nospat, 3)
unique_e_environ <- round(r2_e_all - r2_e_noenv,  3)

cat("\n--- Stage 1: Excavator pseudo-R² variance decomposition ---\n")
cat("Total McFadden R² (all groups):", r2_e_all,        "\n")
cat("Unique local:                  ", unique_e_local,   "\n")
cat("Unique spatial:                ", unique_e_spatial, "\n")
cat("Unique environmental:          ", unique_e_environ, "\n")

# ── LRT significance tests ────────────────────────────────────
cat("\nLRT — local unique fraction:\n")
print(anova(e_vp_no_local, e_vp_all))

cat("\nLRT — spatial unique fraction:\n")
print(anova(e_vp_no_spatial, e_vp_all))

cat("\nLRT — environmental unique fraction:\n")
print(anova(e_vp_no_environ, e_vp_all))

# ── Results ───────────────────────────────────────────────────
# ── Stage 1 results ───────────────────────────────────────────
# Unique local:         McF R² = 0.009, p = 0.042 * — significant
# Unique spatial:       McF R² = 0.000, p = 0.705  — not significant
# Unique environmental: McF R² = 0.004, p = 0.299  — not significant
# Total McFadden R²:    0.012
#
# Local processes are the only significant contributor —
# consistent with total biomass and browsers.
# Rugosity and fishing pressure explain unique variance
# independently of spatial and environmental context.
#
# Spatial and environmental both non-significant — excavators
# are not structured by connectivity or thermal/productivity
# context at this scale.
#
# Compare across groups:
#   Total biomass: local***,  spatial n.s., env .
#   Browsers:      local***,  spatial *,    env n.s.
#   Corallivores:  local n.s., spatial .,   env n.s.
#   Grazers:       local .,   spatial *,    env n.s.
#   Excavators:    local *,   spatial n.s., env n.s.
#
# Excavators most closely mirror total biomass — local
# dominant, spatial absent. Consistent with rugosity being
# the primary driver in the original analysis.
#
# NOTE: McFadden R² values are low throughout (total = 0.012)
# — consistent with other Tweedie groups. LRT significance
# is the primary inferential tool, not R² magnitude.

# ============================================================
#  STAGE 2 — HIERARCHICAL MODEL COMPARISON
#
#  UPDATE family after family selection decision above.
#  Using Tweedie as default — update to lm() if Gaussian
#  log selected.
# ============================================================

# ── Null model ────────────────────────────────────────────────
e_null <- glmmTMB(mean_biomass ~ 1,
                  family = tweedie(link = "log"),
                  data   = excavator_model_data)

# ── Local ecological baseline ─────────────────────────────────
e_local <- glmmTMB(mean_biomass ~ rugosity_sc +
                     log_settlement_grav_sc,
                   family = tweedie(link = "log"),
                   data   = excavator_model_data)

# ── Local + chla ──────────────────────────────────────────────
e_local_chla <- glmmTMB(mean_biomass ~ rugosity_sc +
                          log_settlement_grav_sc +
                          log_chla_sc,
                        family = tweedie(link = "log"),
                        data   = excavator_model_data)

# ── Local + DHW ───────────────────────────────────────────────
e_local_dhw <- glmmTMB(mean_biomass ~ rugosity_sc +
                         log_settlement_grav_sc +
                         log_max_dhw_sc,
                       family = tweedie(link = "log"),
                       data   = excavator_model_data)

# ── Local + environmental context ────────────────────────────
e_local_env <- glmmTMB(mean_biomass ~ rugosity_sc +
                         log_settlement_grav_sc +
                         log_chla_sc +
                         log_max_dhw_sc,
                       family = tweedie(link = "log"),
                       data   = excavator_model_data)

# ── Local + spatial ───────────────────────────────────────────
e_local_spatial <- glmmTMB(mean_biomass ~ rugosity_sc +
                             log_settlement_grav_sc +
                             connectivity_sc,
                           family = tweedie(link = "log"),
                           data   = excavator_model_data)

# ── Local + MPA (governance test) ────────────────────────────
e_local_mpa <- glmmTMB(mean_biomass ~ rugosity_sc +
                         log_settlement_grav_sc +
                         mpa_status,
                       family = tweedie(link = "log"),
                       data   = excavator_model_data)

# ── Global model ──────────────────────────────────────────────
e_global <- glmmTMB(mean_biomass ~ rugosity_sc +
                      log_settlement_grav_sc +
                      connectivity_sc +
                      mpa_status +
                      log_chla_sc +
                      log_max_dhw_sc,
                    family = tweedie(link = "log"),
                    data   = excavator_model_data)

# ── Model list ────────────────────────────────────────────────
excav_model_list <- list(
  "Null"            = e_null,
  "Local"           = e_local,
  "Local + chla"    = e_local_chla,
  "Local + DHW"     = e_local_dhw,
  "Local + env"     = e_local_env,
  "Local + MPA"     = e_local_mpa,
  "Local + spatial" = e_local_spatial,
  "Global"          = e_global
)

# ── AICc ranked table ─────────────────────────────────────────
cat("\n--- Stage 2: Excavator model comparison (AICc ranked) ---\n")
print(make_aicc_df(excav_model_list))

# ── McFadden R² increments over local baseline ────────────────
cat("\n--- Stage 2: Excavator variance explained ---\n")
local_loglik_e <- as.numeric(logLik(e_local))
null_loglik_e2 <- as.numeric(logLik(e_null))

excav_model_list %>%
  imap_dfr(~ tibble(
    Model  = .y,
    McF_R2 = round(1 - (as.numeric(logLik(.x)) /
                          null_loglik_e2), 3)
  )) %>%
  mutate(
    Delta_R2 = round(McF_R2 - (1 - (local_loglik_e /
                                      null_loglik_e2)), 3),
    Delta_R2 = ifelse(Model %in% c("Null", "Local"), NA, Delta_R2)
  ) %>%
  arrange(Model == "Null",
          desc(Model == "Local"),
          desc(McF_R2)) %>%
  print()

# ── Stage 2 results ───────────────────────────────────────────
# Local best (weight = 0.260) but Null competitive
# (ΔAICc = 0.88, weight = 0.168) — local baseline barely
# outperforms intercept-only model.
#
# AICc ranking:
#   Local:        AICc = 687.87, weight = 0.260 — best
#   Local + chla: ΔAICc = 0.30,  weight = 0.224 — competitive
#   Null:         ΔAICc = 0.88,  weight = 0.168 — competitive
#   Local + MPA:  ΔAICc = 1.41,  weight = 0.129 — marginal
#   All others:   ΔAICc > 2.45  — not supported
#
# R² increments:
#   Local + chla: +0.004 — minimal
#   Local + MPA:  +0.006 — minimal
#   Spatial, DHW: +0.000 — nothing
#   Global:       +0.009 — driven by chla and MPA
#
# Overall signal is weak — excavator biomass is largely
# unexplained by the current predictor set.
# Local processes (rugosity + pressure) explain the only
# supported unique variance fraction (Stage 1, p = 0.042)
# but the effect is small and the null model competitive.
#
# Consistent with original analysis where habitat (rugosity)
# was the dominant predictor but explained modest variance.
# The addition of MPA and connectivity does not improve
# inference substantially.

summary(e_local)
summary(e_local_chla)

# ── Coefficients: competitive models ─────────────────────────
#
#                      Local        Local + chla
# Rugosity:            +0.492**     +0.481*
# Settlement gravity:  +0.056 n.s.  -0.223 n.s.
# Chla:                —            -0.383 n.s.
#
# Rugosity is the only significant predictor — positive and
# stable across both competitive models. Structurally complex
# reefs support higher excavator biomass, consistent with
# large-bodied parrotfish being strongly habitat-dependent.
#
# Settlement gravity non-significant throughout — fishing
# pressure has no detectable independent effect on excavator
# biomass. Consistent with complete pressure metric uncertainty
# at pre-analysis selection and the weak pressure signal
# across all functional groups except browsers.
#
# Chla non-significant (p = 0.115) — negative direction
# possibly reflects inverse relationship between productivity
# and reef structural complexity. Not a reliable signal.
#
# Primary conclusion: rugosity is the sole robust predictor
# of large excavator biomass. Consistent with original
# analysis (rugosity best predictor, weight = 0.33).

# ============================================================
#  STAGE 3 — INTERACTION TESTING
# ============================================================

# ── Gate check ────────────────────────────────────────────────
delta_e_global_vs_local <- AICc(e_local) - AICc(e_global)
cat("\nΔAICc (Local vs Global — excavators):",
    round(delta_e_global_vs_local, 2), "\n")

if (delta_e_global_vs_local < 2) {
  cat("NOTE: Global does not outperform Local by ΔAICc > 2.\n",
      "Interaction models fitted for completeness but\n",
      "interpreted with caution.\n")
}

# ── Hypothesis 1: MPA effectiveness depends on larval supply ──
e_int_mpa_conn <- glmmTMB(mean_biomass ~ rugosity_sc +
                            log_settlement_grav_sc +
                            mpa_status * connectivity_sc,
                          family = tweedie(link = "log"),
                          data   = excavator_model_data)

# ── Hypothesis 2: MPA effectiveness depends on fishing pressure
e_int_mpa_press <- glmmTMB(mean_biomass ~ rugosity_sc +
                             mpa_status * log_settlement_grav_sc,
                           family = tweedie(link = "log"),
                           data   = excavator_model_data)

# ── Hypothesis 3: connectivity buffers fishing pressure ───────
e_int_conn_press <- glmmTMB(mean_biomass ~ rugosity_sc +
                              connectivity_sc * log_settlement_grav_sc,
                            family = tweedie(link = "log"),
                            data   = excavator_model_data)

# ── Interaction candidate set ─────────────────────────────────
# Reference: best supported additive model from Stage 2.
# UPDATE model name after Stage 2 results known.
excav_interactions <- list(
  "Local (additive)"        = e_local,      # update if needed
  "MPA × connectivity"      = e_int_mpa_conn,
  "MPA × pressure"          = e_int_mpa_press,
  "Connectivity × pressure" = e_int_conn_press
)

cat("\n--- Stage 3: Excavator interaction comparison ---\n")
print(make_aicc_df(excav_interactions))

# ── Stage 3 results ───────────────────────────────────────────
# Local additive dominates (weight = 0.848).
# No interaction improves on the additive structure.
#
# MPA × pressure:          ΔAICc = 4.92 — not supported
# Connectivity × pressure: ΔAICc = 5.13 — not supported
# MPA × connectivity:      ΔAICc = 8.17 — not supported
#
# Gate check failed (Global worse than Local by 7.84 units)
# — result consistent with expectation. MPA and connectivity
# show no additive support so interactions are untestable.
#
# Conclusion: no interaction supported for excavators.
# Rugosity is the sole predictor of large excavator biomass.
# Neither spatial management nor connectivity modifies the
# habitat-biomass relationship at this scale.
#
# Contrast with browsers (MPA × pressure decisive, weight =
# 0.965) — excavators do not show the same management-
# dependent pressure response despite being spearfishing
# targets. Likely reflects different vulnerability profiles:
# browsers are more accessible to artisanal fishing than
# large-bodied parrotfish.

# ============================================================
#  STAGE 4 — SENSITIVITY ANALYSIS
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
# Mirrors Local model — best supported from Stage 3.
# Only pressure metric swapped. Rugosity retained as the
# key predictor whose robustness we are testing.

e_sens_settpop <- glmmTMB(mean_biomass ~ rugosity_sc +
                            log_settlement_pop_sc,
                          family = tweedie(link = "log"),
                          data   = excavator_model_data)

e_sens_mktgrav <- glmmTMB(mean_biomass ~ rugosity_sc +
                            log_market_gravity_sc,
                          family = tweedie(link = "log"),
                          data   = excavator_model_data)

cat("\n--- Stage 4a: Excavator sensitivity — pressure metrics ---\n")
cat("Settlement population:\n")
print(round(summary(e_sens_settpop)$coefficients$cond, 4))
cat("\nMarket gravity:\n")
print(round(summary(e_sens_mktgrav)$coefficients$cond, 4))

# ── Sensitivity alternative pressure metrics: ────────────────────────────────
# Rugosity robust across all three metrics:
#   Settlement gravity:  β = +0.492, p = 0.009
#   Settlement pop.:     β = +0.590, p = 0.015
#   Market gravity:      β = +0.513, p = 0.016
# Significant and positive throughout — rugosity is the
# only robust predictor of excavator biomass regardless
# of pressure metric. Pressure metrics themselves
# non-significant across all three (p = 0.472–0.788).

# ── (b) Transect-level replication ───────────────────────────
# ~57.8% zeros at transect level — Tweedie required.
# ZI Tweedie tested given high zero proportion.

excav_transect_data <- excavator_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(excav_transect_data$transect_excavator_biomass == 0),
    "/", nrow(excav_transect_data),
    "(", round(mean(excav_transect_data$transect_excavator_biomass == 0), 3), ")\n")

# ── Transect family selection ─────────────────────────────────
e_trans_tw <- glmmTMB(transect_excavator_biomass ~ rugosity_sc +
                        log_settlement_grav_sc +
                        connectivity_sc +
                        mpa_status +
                        log_chla_sc +
                        log_max_dhw_sc +
                        (1 | site),
                      family = tweedie(link = "log"),
                      data   = excav_transect_data)

e_trans_tw_zi <- glmmTMB(transect_excavator_biomass ~ rugosity_sc +
                           log_settlement_grav_sc +
                           connectivity_sc +
                           mpa_status +
                           log_chla_sc +
                           log_max_dhw_sc +
                           (1 | site),
                         family    = tweedie(link = "log"),
                         ziformula = ~1,
                         data      = excav_transect_data)

e_trans_res    <- simulateResiduals(e_trans_tw,    n = 1000)
e_trans_res_zi <- simulateResiduals(e_trans_tw_zi, n = 1000)

plot(e_trans_res);    testZeroInflation(e_trans_res)
plot(e_trans_res_zi); testZeroInflation(e_trans_res_zi)

cat("\n--- Transect family selection ---\n")
print(make_aicc_df(list(
  "Tweedie"    = e_trans_tw,
  "ZI Tweedie" = e_trans_tw_zi
)))

# Stage 4b — Transect family selection:
# Standard Tweedie selected (weight = 0.752, ΔAICc = 2.21).
# ZI test non-significant for both models (p > 0.68).
# Despite 56.4% transect zeros, standard Tweedie sufficient.


# ── Transect hierarchical sequence ────────────────────────────
# Stage 3 additive Local model is best — no interaction
# supported. Transect best model mirrors site-level Local.

e_trans_null <- glmmTMB(transect_excavator_biomass ~ 1 +
                          (1 | site),
                        family = tweedie(link = "log"),
                        data   = excav_transect_data)

e_trans_local <- glmmTMB(transect_excavator_biomass ~ rugosity_sc +
                           log_settlement_grav_sc +
                           (1 | site),
                         family = tweedie(link = "log"),
                         data   = excav_transect_data)

e_trans_local_env <- glmmTMB(transect_excavator_biomass ~ rugosity_sc +
                               log_settlement_grav_sc +
                               log_chla_sc +
                               log_max_dhw_sc +
                               (1 | site),
                             family = tweedie(link = "log"),
                             data   = excav_transect_data)

e_trans_local_spatial <- glmmTMB(transect_excavator_biomass ~ rugosity_sc +
                                   log_settlement_grav_sc +
                                   connectivity_sc +
                                   (1 | site),
                                 family = tweedie(link = "log"),
                                 data   = excav_transect_data)

e_trans_local_mpa <- glmmTMB(transect_excavator_biomass ~ rugosity_sc +
                               log_settlement_grav_sc +
                               mpa_status +
                               (1 | site),
                             family = tweedie(link = "log"),
                             data   = excav_transect_data)

e_trans_global <- glmmTMB(transect_excavator_biomass ~ rugosity_sc +
                            log_settlement_grav_sc +
                            connectivity_sc +
                            mpa_status +
                            log_chla_sc +
                            log_max_dhw_sc +
                            (1 | site),
                          family = tweedie(link = "log"),
                          data   = excav_transect_data)

excav_transect_list <- list(
  "Null"            = e_trans_null,
  "Local"           = e_trans_local,
  "Local + env"     = e_trans_local_env,
  "Local + spatial" = e_trans_local_spatial,
  "Local + MPA"     = e_trans_local_mpa,
  "Global"          = e_trans_global
)

cat("\n--- Stage 4b: Excavator transect-level sensitivity ---\n")
print(make_aicc_df(excav_transect_list))

# ── Interpretation ────────────────────────────────────────────
# Stage 4b — Transect hierarchical comparison:
# Null best (weight = 0.268) — Local competitive (ΔAICc = 0.38).
# Combined weight of Null + Local = 0.489.
# No model clearly supported over the null at transect level.
#
# This is a weaker result than site level where Local was
# best (weight = 0.260, Null weight = 0.168). At transect
# level the signal is essentially absent — within-site
# variation in excavator biomass swamps the rugosity signal.
#
# Primary conclusion: rugosity signal is robust in
# coefficient direction and significance across pressure
# metrics (Stage 4a) but weakens at transect level where
# within-site variation is large. Site-level result is
# the primary inference — transect level confirms direction
# but not significance.