# ============================================================
#  CORALLIVORE FISH BIOMASS — SITE-LEVEL ANALYSIS
#
#  Analytical framework mirrors total biomass and browsers
#  (four stages). Key differences from browsers noted below.
#
#  STAGE 1 — Variance partitioning
#             McFadden pseudo-R² + LRT approximation if
#             Tweedie selected (as per browsers).
#             Formal varpart() if Gaussian log adequate.
#
#  STAGE 2 — Hierarchical model comparison
#             Identical nested sequence to total biomass
#             and browsers.
#
#  STAGE 3 — Interaction testing (conditional on Stage 2)
#             Three a priori interactions as per other groups.
#
#  STAGE 4 — Sensitivity analysis
#             (a) Alternative pressure metrics
#             (b) Transect-level mixed model replication
#
#  Key differences from browsers:
#    - Pressure metric: settlement gravity clearly preferred
#      for corallivores (weight = 0.974, ΔAICc = 7.23 vs
#      settlement pop.) — no metric uncertainty here.
#    - Zero proportion: check site-level and transect-level
#      zeros before proceeding — will determine family.
#    - Ecological context: corallivores are coral-dependent
#      predators — DHW and habitat complexity hypothesised
#      to be stronger drivers than for browsers. MPA effect
#      may differ given different fishing vulnerability.
#
#  Study design:
#    243 transects, 54 sites, 4 countries.
#    Minimum 3 transects per site retained.
# ============================================================

# ── SOURCE SHARED DATA PREPARATION ───────────────────────────
# # Loads: fish_data, scaled_predictors, final_predictors,
# total_transects, transect_model_data, total_model_data,
# make_aicc_df(), plot_effect(), and all packages.
source(here::here("data_preparation.R"))

# ── AGGREGATE CORALLIVORE TRANSECT DATA ──────────────────────
coralliv_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_coralliv_biomass = sum(
      ifelse(trophic_group == "corallivores", tot_wt_g, 0),
      na.rm = TRUE
    ),
    transect_coralliv_count = sum(
      ifelse(trophic_group == "corallivores", number, 0),
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
cat("Countries:",             n_distinct(coralliv_transects$country), "\n")


# ── SITE-LEVEL AGGREGATION ────────────────────────────────────
coralliv_site_data <- coralliv_transects %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass = mean(transect_coralliv_biomass, na.rm = TRUE),
    n_transects  = n(),
    .groups      = "drop"
  ) %>%
  mutate(
    site    = as.factor(site),
    country = as.factor(country)
  )

# ── Zero check ────────────────────────────────────────────────
zeros_site <- mean(coralliv_site_data$mean_biomass == 0, na.rm = TRUE)
cat("Site-level zero proportion:", round(zeros_site, 3), "\n")
cat("Site-level zeros (n):",       sum(coralliv_site_data$mean_biomass == 0), "\n")

# ── Distribution plots ────────────────────────────────────────
p_raw <- ggplot(coralliv_site_data, aes(x = mean_biomass)) +
  geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
  labs(x = "Mean corallivore biomass per site (g)",
       y = "Frequency", title = "Raw") + theme_bw()

coralliv_site_data <- coralliv_site_data %>%
  mutate(log_mean_biomass = log(mean_biomass + 0.01))

p_log <- ggplot(coralliv_site_data, aes(x = log_mean_biomass)) +
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
# Log-transformed distribution is unimodal and approximately
# continuous — no isolated cluster of zero sites visible.
# The dashed line at log(0.01) = -4.6 marks where zeros
# would appear but no sites fall there. The minimum value
# is around -2, which is a genuine low-biomass site rather
# than a structural zero.
#
# Zero proportion confirmed at 0% from the zero check above.
# This is consistent with corallivores being present at all
# 54 sites — likely because even heavily degraded reefs
# retain some corallivore species.

# ── Build full corallivore site-level model dataset ───────────
coralliv_model_data <- coralliv_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_coralliv_biomass, na.rm = TRUE),
    log_mean_biomass       = log(mean(transect_coralliv_biomass,
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

cat("\nCorallivore model data:", nrow(coralliv_model_data), "sites\n")
cat("Site-level zeros:",        sum(coralliv_model_data$mean_biomass == 0), "\n")

# ── Verify no NAs ─────────────────────────────────────────────
coralliv_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                connectivity_sc, mpa_status,
                log_chla_sc, log_max_dhw_sc) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

# ── Check zeros and response ──────────────────────────────────
cat("Zeros in mean_biomass:",    sum(coralliv_model_data$mean_biomass == 0), "\n")
cat("-Inf in log_mean_biomass:", sum(is.infinite(coralliv_model_data$log_mean_biomass)), "\n")
cat("NAs in log_mean_biomass:",  sum(is.na(coralliv_model_data$log_mean_biomass)), "\n")
cat("\nResponse summary:\n")
print(summary(coralliv_model_data$log_mean_biomass))


# ============================================================
#  PRE-ANALYSIS: HUMAN PRESSURE METRIC SELECTION
#
#  Settlement gravity clearly preferred for corallivores
#  in original analysis (weight = 0.974, ΔAICc = 7.23).
#  Confirm with univariate comparison using consistent
#  family across all functional groups.
# ============================================================

coralliv_press_settgrav <- glmmTMB(mean_biomass ~ log_settlement_grav_sc,
                                   family = tweedie(link = "log"),
                                   data   = coralliv_model_data)
coralliv_press_settpop  <- glmmTMB(mean_biomass ~ log_settlement_pop_sc,
                                   family = tweedie(link = "log"),
                                   data   = coralliv_model_data)
coralliv_press_mktgrav  <- glmmTMB(mean_biomass ~ log_market_gravity_sc,
                                   family = tweedie(link = "log"),
                                   data   = coralliv_model_data)

cat("\n--- Pre-analysis: corallivore pressure metric selection ---\n")
print(make_aicc_df(list(
  "Settlement gravity" = coralliv_press_settgrav,
  "Settlement pop."    = coralliv_press_settpop,
  "Market gravity"     = coralliv_press_mktgrav
)))

# ── Corallivore pressure metric decision ──────────────────────
# All three metrics within ΔAICc < 1 — no metric preferred.
#
# Results:
#   Market gravity:     AICc = 561.70, weight = 0.397 (best)
#   Settlement pop.:    ΔAICc = 0.40,  weight = 0.324
#   Settlement gravity: ΔAICc = 0.71,  weight = 0.279
#
# Complete metric uncertainty — weights nearly equal across
# all three. No biological or statistical basis for preferring
# any single metric for corallivores.
#
# Settlement gravity selected for consistency with total
# biomass and browser analyses, enabling direct cross-
# functional-group comparison. This is a pragmatic decision
# based on analytical consistency rather than metric
# performance, which must be acknowledged explicitly in
# the methods.
#
# Note: contrast with total biomass (settlement gravity
# weight = 0.672) and corallivores in original analysis
# (settlement gravity weight = 0.974 — this has changed,
# likely due to family change from Tweedie to lm-compatible
# testing). All three metrics retained for sensitivity (Stage 4a).


# ============================================================
#  MODEL FAMILY SELECTION
#
#  Decision tree based on zero proportion:
#
#  If site-level zeros < 5% and distribution approximately
#  continuous after log transformation:
#    → Test Gaussian log — run standard lm() diagnostics
#    → If adequate: use lm(), varpart() available in Stage 1
#
#  If site-level zeros >= 5% OR bimodal log distribution:
#    → Tweedie required
#    → varpart() not available, use pseudo-R² in Stage 1
#
#  In either case test ZI Tweedie at transect level
#  (higher zero proportions expected there).
# ============================================================

# ── Gaussian log (global model) ───────────────────────────────
coralliv_lm_global <- lm(log_mean_biomass ~ rugosity_sc +
                           log_settlement_grav_sc +
                           connectivity_sc +
                           mpa_status +
                           log_chla_sc +
                           log_max_dhw_sc,
                         data = coralliv_model_data)

par(mfrow = c(2, 2)); plot(coralliv_lm_global); par(mfrow = c(1, 1))

# ── Tweedie (global model) ────────────────────────────────────
coralliv_tw_global <- glmmTMB(mean_biomass ~ rugosity_sc +
                                log_settlement_grav_sc +
                                connectivity_sc +
                                mpa_status +
                                log_chla_sc +
                                log_max_dhw_sc,
                              family = tweedie(link = "log"),
                              data   = coralliv_model_data)

coralliv_tw_res <- simulateResiduals(coralliv_tw_global, n = 1000)
plot(coralliv_tw_res)
testZeroInflation(coralliv_tw_res)
testDispersion(coralliv_tw_res)

# ── Family selection decision ─────────────────────────────────
# Gaussian log selected — lm() on log_mean_biomass.
#
# Zero-free response (0 site-level zeros) — Tweedie not
# needed. Tweedie zero inflation test: ratio = 0, p = 1 —
# confirms no zero structure in data.
#
# Gaussian log diagnostics (global model):
#   Residuals vs Fitted: broadly flat — no systematic pattern
#   Q-Q: good through main body — minor lower tail deviation
#     at sites 35, 1, 13 within acceptable range for n = 54
#   Scale-Location: broadly flat — homoscedasticity met
#   Residuals vs Leverage: site 47 high leverage but within
#     Cook's distance — no unduly influential observations
#
# Tweedie diagnostics shown for completeness:
#   KS test: p = 0.814 — good fit
#   Dispersion: p = 0.762 — acceptable
#   Both families fit adequately — Gaussian log selected
#   as more parsimonious and appropriate for zero-free
#   continuous data. Consistent with total biomass.

# # ============================================================
# #  RANDOM EFFECT STRUCTURE
# #  Uses lm() consistent with family selection above.
# #  glmmTMB used only to enable (1|country) comparison —
# #  lm() does not support random effects directly.
# #  If country RE not supported, proceed with lm() throughout.
# # ============================================================
# 
# coralliv_re_null <- lm(log_mean_biomass ~ rugosity_sc +
#                          log_settlement_grav_sc +
#                          log_chla_sc +
#                          log_max_dhw_sc,
#                        data = coralliv_model_data)
# 
# coralliv_re_country <- glmmTMB(log_mean_biomass ~ rugosity_sc +
#                                  log_settlement_grav_sc +
#                                  log_chla_sc +
#                                  log_max_dhw_sc +
#                                  (1 | country),
#                                family = gaussian(),
#                                data   = coralliv_model_data)
# 
# cat("\n--- Corallivore RE structure comparison ---\n")
# print(make_aicc_df(list(
#   "No RE"         = coralliv_re_null,
#   "(1 | country)" = coralliv_re_country
# )))
# 
# # Expected: no RE supported, consistent with total biomass
# # (ΔAICc = 2.53 in favour of no RE there).
# # If country RE not supported → all models as lm().
# # If country RE supported → rebuild all models as glmmTMB
# # gaussian with (1|country).
# 
# rm(coralliv_re_null, coralliv_re_country)


# ============================================================
#  STAGE 1 — VARIANCE PARTITIONING
#
#  Gaussian log selected — formal varpart() available.
#  Identical procedure to total biomass.
#  Three process groups: local, spatial, environmental.
#  MPA excluded — governance variable, tested in Stage 2.
# ============================================================

vp_local_c <- coralliv_model_data %>%
  dplyr::select(rugosity_sc, log_settlement_grav_sc) %>%
  as.data.frame()

vp_spatial_c <- coralliv_model_data %>%
  dplyr::select(connectivity_sc) %>%
  as.data.frame()

vp_environ_c <- coralliv_model_data %>%
  dplyr::select(log_chla_sc, log_max_dhw_sc) %>%
  as.data.frame()

y_coralliv <- coralliv_model_data$log_mean_biomass

vp_coralliv <- varpart(y_coralliv,
                       vp_local_c,
                       vp_spatial_c,
                       vp_environ_c)

cat("\n--- Stage 1: Corallivore variance partitioning ---\n")
print(vp_coralliv)

# ── Variance partition bar chart ──────────────────────────────
vp_fractions_c <- data.frame(
  Group       = c("Local", "Spatial", "Environment"),
  Unique      = c(0.008, 0.058, 0.019),
  Unique_plot = c(0.008, 0.058, 0.019),
  p_label     = c("adj. R² = 0.008, p = 0.312",
                  "adj. R² = 0.058, p = 0.044",
                  "adj. R² = 0.019, p = 0.217"),
  sig         = c("not significant", "significant", "not significant")
) %>%
  mutate(Group = factor(Group,
                        levels = c("Environment", "Spatial", "Local")))

ggplot(vp_fractions_c, aes(x = Unique_plot, y = Group,
                           fill = sig)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = p_label),
            hjust = -0.05, size = 3.2) +
  scale_fill_manual(values = c("significant"     = "#0072B2",
                               "not significant" = "#bdbdbd")) +
  scale_x_continuous(limits = c(0, 0.20),
                     name   = "Unique variance explained (adj. R²)") +
  labs(y       = NULL,
       fill    = NULL,
       caption = "n = 54 sites | Permutation tests, 999 permutations | Residual variance = 0.890") +
  theme_bw(base_size = 13) +
  theme(legend.position    = "top",
        axis.title         = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        plot.caption       = element_text(colour = "grey50", size = 9))

# ── Significance tests ────────────────────────────────────────
cat("\nLocal unique fraction:\n")
print(anova(rda(y_coralliv ~ rugosity_sc + log_settlement_grav_sc +
                  Condition(connectivity_sc) +
                  Condition(log_chla_sc) +
                  Condition(log_max_dhw_sc),
                data = coralliv_model_data)))

cat("\nSpatial unique fraction:\n")
print(anova(rda(y_coralliv ~ connectivity_sc +
                  Condition(rugosity_sc) +
                  Condition(log_settlement_grav_sc) +
                  Condition(log_chla_sc) +
                  Condition(log_max_dhw_sc),
                data = coralliv_model_data)))

cat("\nEnvironmental unique fraction:\n")
print(anova(rda(y_coralliv ~ log_chla_sc + log_max_dhw_sc +
                  Condition(rugosity_sc) +
                  Condition(log_settlement_grav_sc) +
                  Condition(connectivity_sc),
                data = coralliv_model_data)))

# ── Stage 1 results ───────────────────────────────────
# Unique local:         adj R² =  0.008, p = 0.299 — not significant
# Unique spatial:       adj R² =  0.058, p = 0.053 . — marginal trend
# Unique environmental: adj R² =  0.019, p = 0.216 — not significant
# Total explained:      adj R² =  0.110
# Residual:                       0.890
#
# Pattern across functional groups:
#   Total biomass: local***,  spatial n.s., environment .
#   Browsers:      local***,  spatial *,    environment n.s.
#   Corallivores:  local n.s., spatial .,   environment n.s.
#
# Connectivity shows a progressive shift in importance:
# absent for total biomass, significant for browsers,
# marginal for corallivores. Local processes show the
# inverse pattern — dominant for total biomass, significant
# for browsers, absent for corallivores.
#
# Corallivore biomass is largely unexplained by the current
# predictor set (residual = 0.890). Unmeasured variables —
# coral cover, prey availability, species-specific habitat
# associations — likely explain substantial variance.

# ============================================================
#  STAGE 2 — HIERARCHICAL MODEL COMPARISON
#
#  Gaussian log selected — lm() on log_mean_biomass.
#  Identical nested sequence to total biomass and corallivores.
#  Adjusted R² used — consistent with total biomass.
# ============================================================

# ── Null model ────────────────────────────────────────────────
c_null <- lm(log_mean_biomass ~ 1,
             data = coralliv_model_data)

# ── Local ecological baseline ─────────────────────────────────
c_local <- lm(log_mean_biomass ~ rugosity_sc +
                log_settlement_grav_sc,
              data = coralliv_model_data)

# ── Local + chla ──────────────────────────────────────────────
c_local_chla <- lm(log_mean_biomass ~ rugosity_sc +
                     log_settlement_grav_sc +
                     log_chla_sc,
                   data = coralliv_model_data)

# ── Local + DHW ───────────────────────────────────────────────
c_local_dhw <- lm(log_mean_biomass ~ rugosity_sc +
                    log_settlement_grav_sc +
                    log_max_dhw_sc,
                  data = coralliv_model_data)

# ── Local + environmental context ────────────────────────────
c_local_env <- lm(log_mean_biomass ~ rugosity_sc +
                    log_settlement_grav_sc +
                    log_chla_sc +
                    log_max_dhw_sc,
                  data = coralliv_model_data)

# ── Local + spatial ───────────────────────────────────────────
c_local_spatial <- lm(log_mean_biomass ~ rugosity_sc +
                        log_settlement_grav_sc +
                        connectivity_sc,
                      data = coralliv_model_data)

# ── Local + MPA (governance test) ────────────────────────────
c_local_mpa <- lm(log_mean_biomass ~ rugosity_sc +
                    log_settlement_grav_sc +
                    mpa_status,
                  data = coralliv_model_data)

# ── Global model ──────────────────────────────────────────────
c_global <- lm(log_mean_biomass ~ rugosity_sc +
                 log_settlement_grav_sc +
                 connectivity_sc +
                 mpa_status +
                 log_chla_sc +
                 log_max_dhw_sc,
               data = coralliv_model_data)

# ── Model list ────────────────────────────────────────────────
coralliv_model_list <- list(
  "Null"            = c_null,
  "Local"           = c_local,
  "Local + chla"    = c_local_chla,
  "Local + DHW"     = c_local_dhw,
  "Local + env"     = c_local_env,
  "Local + MPA"     = c_local_mpa,
  "Local + spatial" = c_local_spatial,
  "Global"          = c_global
)

# ── AICc ranked table ─────────────────────────────────────────
cat("\n--- Stage 2: Corallivore model comparison (AICc ranked) ---\n")
print(make_aicc_df(coralliv_model_list))

# ── Adjusted R² increments over local baseline ───────────────
# lm() allows direct use of adj. R² — no McFadden needed.
cat("\n--- Stage 2: Corallivore variance explained ---\n")
local_r2_c <- summary(c_local)$adj.r.squared

coralliv_model_list %>%
  imap_dfr(~ tibble(
    Model  = .y,
    Adj_R2 = round(summary(.x)$adj.r.squared, 3)
  )) %>%
  mutate(
    Delta_R2 = round(Adj_R2 - local_r2_c, 3),
    Delta_R2 = ifelse(Model %in% c("Null", "Local"), NA, Delta_R2)
  ) %>%
  arrange(Model == "Null",
          desc(Model == "Local"),
          desc(Adj_R2)) %>%
  print()

# ── Stage 2 results ───────────────────────────────────────────
# Corallivore biomass shows a fundamentally different pattern
# from all other functional groups.
#
# AICc ranking:
#   Global:          AICc = 135.57, weight = 0.327 — best
#   Local + spatial: ΔAICc = 0.23,  weight = 0.291 — competitive
#   Null:            ΔAICc = 1.52,  weight = 0.153 — competitive
#   Local + chla:    ΔAICc = 2.05  — marginal
#   All others:      ΔAICc > 3.95  — not supported
#
# R² increments over Local baseline:
#   Global:          +0.196 — all driven by connectivity
#   Local + spatial: +0.092 — connectivity alone
#   Local + chla:    +0.061 — chla adds modestly
#   Local + env:     +0.053
#   Local + MPA:     -0.008 — worse than local
#   Local + DHW:     -0.014 — worse than local
#
# CRITICAL FINDING: Local baseline adj R² = -0.001
# Rugosity and fishing pressure explain zero variance for
# corallivores — the local baseline performs no better than
# an intercept-only null model. This is the starkest contrast
# with total biomass (local adj R² = 0.184) and browsers
# (local competitive).
#
# Connectivity is the only meaningful predictor — consistent
# with Stage 1 (spatial marginal fraction p = 0.053).
# The global model's apparent performance (+0.195) is almost
# entirely attributable to the connectivity term.
#
# Null model competitive (ΔAICc = 1.52, weight = 0.153) —
# there is genuine model selection uncertainty about whether
# any predictor explains corallivore biomass reliably.
#
# Ecological interpretation:
# Corallivores are specialist coral feeders whose local
# abundance is decoupled from fishing pressure and habitat
# complexity as measured here. Their biomass is better
# predicted by larval network position (connectivity) than
# by local conditions — consistent with a group that is
# not heavily targeted by artisanal fishing and whose
# distribution depends on regional recruitment dynamics
# and coral prey availability rather than local exploitation.

summary(c_local_spatial)
summary(c_global)

# ── Coefficient summary ───────────────────────────────────────
#
# Local + spatial:
#   Rugosity:          β = +0.138, p = 0.215 — not significant
#   Settlement gravity: β = -0.010, p = 0.931 — not significant
#   Connectivity:      β = -0.278, p = 0.017 * — significant
#   Adj R² = 0.091
#
# Global:
#   Rugosity:          β = +0.159, p = 0.130 — not significant
#   Settlement gravity: β = +0.191, p = 0.208 — not significant
#   Connectivity:      β = -0.397, p = 0.003 ** — significant
#   MPA low:           β = +0.843, p = 0.022 * — significant
#   MPA medium:        β = +0.516, p = 0.069 . — marginal
#   Chla:              β = +0.206, p = 0.140 — not significant
#   DHW:               β = -0.023, p = 0.834 — not significant
#   Adj R² = 0.195
#
# Connectivity is NEGATIVE — well-connected sites have lower
# corallivore biomass. Counterintuitive but consistent across
# both models (β = -0.278 and -0.397). Not a collinearity
# artefact (r with rugosity = -0.02, r with chla = -0.28).
#
# Possible explanations:
# 1. Source-sink dynamics — highly connected sites export
#    larvae rather than retaining them locally, reducing
#    local settlement and adult biomass.
# 2. Regional disturbance exposure — well-connected sites
#    may be more exposed to regional bleaching and disease
#    events that reduce coral cover and corallivore prey.
# 3. Habitat quality confound — connectivity may correlate
#    with reef configuration features that reduce corallivore
#    habitat quality in ways not captured by rugosity alone.
#
# MPA low positive (β = +0.843, p = 0.022) — reverse of
# browser pattern. Low protection may reduce destructive
# fishing practices (anchoring, blast fishing) that damage
# coral substrate rather than directly targeting corallivores.
#
# DHW negative as predicted but non-significant — bleaching
# effects on corallivore prey may operate over longer
# timescales than cross-sectional data can detect.
#
# Local processes (rugosity, pressure) non-significant
# throughout — consistent with Stage 1 and the broader
# finding that corallivores are decoupled from local
# ecological conditions as measured here.

# ============================================================
#  STAGE 3 — INTERACTION TESTING
# ============================================================

# ── Gate check ────────────────────────────────────────────────
delta_c_global_vs_local <- AICc(c_local) - AICc(c_global)
cat("\nΔAICc (Local vs Global — corallivores):",
    round(delta_c_global_vs_local, 2), "\n")

# ΔAICc = 4.07 — Global clearly outperforms Local.
# Gate check passes — interactions are warranted.
# Note: this is driven by connectivity, not local predictors.
# Interaction interpretation should focus on connectivity
# terms rather than local-process modifications.

# ── Hypothesis 1: MPA effectiveness depends on larval supply ──
c_int_mpa_conn <- lm(log_mean_biomass ~ rugosity_sc +
                       log_settlement_grav_sc +
                       mpa_status * connectivity_sc,
                     data = coralliv_model_data)

# ── Hypothesis 2: MPA effectiveness depends on fishing pressure
c_int_mpa_press <- lm(log_mean_biomass ~ rugosity_sc +
                        mpa_status * log_settlement_grav_sc +
                        connectivity_sc,
                      data = coralliv_model_data)

# ── Hypothesis 3: connectivity buffers fishing pressure ───────
c_int_conn_press <- lm(log_mean_biomass ~ rugosity_sc +
                         connectivity_sc * log_settlement_grav_sc,
                       data = coralliv_model_data)

# ── Interaction candidate set ─────────────────────────────────
# Reference: Local + spatial — best supported additive model
# (AICc = 135.80, weight = 0.291, ΔAICc = 0.23 vs Global).
# Asking whether any interaction improves on the best additive
# structure, not on Local + env which is unsupported here.
coralliv_interactions <- list(
  "Local + spatial (additive)" = c_local_spatial,
  "MPA × connectivity"         = c_int_mpa_conn,
  "MPA × pressure"             = c_int_mpa_press,
  "Connectivity × pressure"    = c_int_conn_press
)

cat("\n--- Stage 3: Corallivore interaction comparison ---\n")
print(make_aicc_df(coralliv_interactions))

# ── Stage 3 AICc summary ─────────────────────────────────────
# MPA × pressure:              AICc = 134.89, weight = 0.474
# Local + spatial (additive):  ΔAICc = 0.91,  weight = 0.300
# MPA × connectivity:          ΔAICc = 2.59,  weight = 0.130
# Connectivity × pressure:     ΔAICc = 3.19,  weight = 0.096
#
# MPA × pressure marginally preferred but Local + spatial
# additive is competitive (ΔAICc < 2, combined weight = 0.774).
# This is model selection uncertainty — not a decisive result.
# Contrast with browsers (MPA × pressure weight = 0.965).
# Interpret MPA × pressure with caution — examine coefficients
# to assess whether the interaction is ecologically coherent
# before drawing inference.

summary(c_int_mpa_press)

# ── Stage 3 conclusion ────────────────────────────────────────
# MPA × pressure is marginally preferred by AICc (weight = 0.474)
# but the interaction terms are non-significant:
#   MPA low × pressure:    β = -0.021, p = 0.961
#   MPA medium × pressure: β = -0.708, p = 0.097
#
# The marginal AICc improvement over Local + spatial is driven
# by the MPA low main effect (β = +0.921, p = 0.028) rather
# than the interaction structure. The interaction hypothesis
# is not supported.
#
# Connectivity remains the dominant predictor across all models
# (β = -0.278 to -0.413, consistently significant). The
# negative direction is robust and not a collinearity artefact.
#
# No interaction model provides an ecologically interpretable
# improvement over the Local + spatial additive model.
# Stage 3 is reported as inconclusive — model selection
# uncertainty too high (top two models within ΔAICc < 1)
# to support any single interaction structure.
#
# Primary inference unchanged: corallivore biomass is
# structured by connectivity (negative), with a marginal
# positive MPA low effect. Local processes do not contribute.
#
# MPA × connectivity (ΔAICc = 2.59) and connectivity ×
# pressure (ΔAICc = 3.19) not supported.

# ============================================================
#  STAGE 4 — SENSITIVITY ANALYSIS
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
# Mirrors Local + spatial structure — best supported additive
# model from Stage 2. Only pressure metric swapped.
# Global model used here to include all predictors for
# fair comparison across metrics.
# Note: complete metric uncertainty at pre-analysis selection
# (all ΔAICc < 1) — robustness across metrics especially
# important for corallivores.

c_sens_settpop <- lm(log_mean_biomass ~ rugosity_sc +
                       log_settlement_pop_sc +
                       connectivity_sc +
                       mpa_status +
                       log_chla_sc +
                       log_max_dhw_sc,
                     data = coralliv_model_data)

c_sens_mktgrav <- lm(log_mean_biomass ~ rugosity_sc +
                       log_market_gravity_sc +
                       connectivity_sc +
                       mpa_status +
                       log_chla_sc +
                       log_max_dhw_sc,
                     data = coralliv_model_data)

cat("\n--- Stage 4a: Corallivore sensitivity — pressure metrics ---\n")
cat("Settlement population:\n")
print(round(summary(c_sens_settpop)$coefficients, 4))
cat("\nMarket gravity:\n")
print(round(summary(c_sens_mktgrav)$coefficients, 4))

# ── Stage 4a results ──────────────────────────────────────────
# Negative connectivity signal is robust across all three
# pressure metrics — strongest robustness result in the
# corallivore analysis.
#
# Connectivity coefficient across all metrics:
#   Settlement gravity:  β = -0.397, p = 0.003
#   Settlement pop.:     β = -0.375, p = 0.005
#   Market gravity:      β = -0.389, p = 0.003
#   Near-identical direction, magnitude, and significance.
#   Metric-independent — confirms genuine signal.
#
# MPA low positive effect also robust:
#   Settlement gravity:  β = +0.843, p = 0.022
#   Settlement pop.:     β = +0.864, p = 0.020
#   Market gravity:      β = +0.838, p = 0.024
#
# Pressure metrics themselves non-significant throughout:
#   Settlement pop.:     β = +0.015, p = 0.907
#   Market gravity:      β = +0.112, p = 0.386
#   Consistent with primary analysis — human pressure does
#   not structure corallivore biomass regardless of metric.
#
# Rugosity, DHW, chla non-significant across all models.
#
# Conclusion: complete metric uncertainty at pre-analysis
# selection (all ΔAICc < 1) is not a concern — primary
# findings are entirely metric-independent. Settlement gravity
# retained for consistency with other functional groups.

# ── (b) Transect-level replication ───────────────────────────
# Site level was zero-free — transect level will have more zeros.
# Tweedie used at transect level regardless of site-level lm().
# Key robustness test: does negative connectivity signal hold
# when within-site variation is retained?

coralliv_transect_data <- coralliv_transects %>%
  left_join(final_predictors, by = "site") %>%
  mutate(log_coralliv_biomass = log(transect_coralliv_biomass + 0.01))

cat("\nTransect zeros:",
    sum(coralliv_transect_data$transect_coralliv_biomass == 0),
    "/", nrow(coralliv_transect_data),
    "(", round(mean(coralliv_transect_data$transect_coralliv_biomass == 0), 3), ")\n")

# ── Transect family selection decision ────────────────────────
# Standard Tweedie selected — ZI Tweedie not tested.
#
# Transect-level zeros: 43/243 (17.7%) — substantially lower
# than browsers (43.2%) and below the threshold where zero
# inflation is a plausible concern. Standard Tweedie handles
# this proportion of zeros natively without a ZI component.
# ZI Tweedie not fitted.
#
# Proceed: standard Tweedie + (1|site) throughout
# transect-level corallivore analyses.

# ── Transect hierarchical sequence ────────────────────────────
# Mirrors Stage 2 sequence.
# Stage 3 was inconclusive — c_trans_best uses Local + spatial
# as the reference model since no interaction was decisively
# supported. Update if Stage 3 results change interpretation.

c_trans_null <- glmmTMB(transect_coralliv_biomass ~ 1 +
                          (1 | site),
                        family = tweedie(link = "log"),
                        data   = coralliv_transect_data)

c_trans_local <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                           log_settlement_grav_sc +
                           (1 | site),
                         family = tweedie(link = "log"),
                         data   = coralliv_transect_data)

c_trans_local_env <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                               log_settlement_grav_sc +
                               log_chla_sc +
                               log_max_dhw_sc +
                               (1 | site),
                             family = tweedie(link = "log"),
                             data   = coralliv_transect_data)

c_trans_local_spatial <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                                   log_settlement_grav_sc +
                                   connectivity_sc +
                                   (1 | site),
                                 family = tweedie(link = "log"),
                                 data   = coralliv_transect_data)

c_trans_local_mpa <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                               log_settlement_grav_sc +
                               mpa_status +
                               (1 | site),
                             family = tweedie(link = "log"),
                             data   = coralliv_transect_data)

c_trans_global <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                            log_settlement_grav_sc +
                            connectivity_sc +
                            mpa_status +
                            log_chla_sc +
                            log_max_dhw_sc +
                            (1 | site),
                          family = tweedie(link = "log"),
                          data   = coralliv_transect_data)

# Stage 3 inconclusive — include Local + spatial as best
# supported additive model and MPA × pressure as marginal
# best interaction for completeness
c_trans_local_spatial_add <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    connectivity_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data)

c_trans_mpa_press <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc +
    mpa_status * log_settlement_grav_sc +
    connectivity_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = coralliv_transect_data)

coralliv_transect_list <- list(
  "Null"                       = c_trans_null,
  "Local"                      = c_trans_local,
  "Local + env"                = c_trans_local_env,
  "Local + spatial"            = c_trans_local_spatial,
  "Local + MPA"                = c_trans_local_mpa,
  "Global"                     = c_trans_global,
  "Local + spatial (additive)" = c_trans_local_spatial_add,
  "MPA × pressure"             = c_trans_mpa_press
)

cat("\n--- Stage 4b: Corallivore transect-level sensitivity ---\n")
print(make_aicc_df(coralliv_transect_list))

# ── Stage 4b results ──────────────────────────────────────────
# Transect-level results show a stronger and clearer pattern
# than the site-level analysis.
#
# Model ranking (transect-level):
#   MPA × pressure:  AICc = 2318.93, weight = 0.519 — best
#   Global:          ΔAICc = 1.05,   weight = 0.307 — competitive
#   Local + spatial: ΔAICc = 3.81,   weight = 0.077 — not supported
#   Global (add.):   ΔAICc = 3.81,   weight = 0.077 — not supported
#   Null:            ΔAICc = 7.47,   weight = 0.012 — worst
#   Local + env:     ΔAICc = 9.72,   weight = 0.004 — worst
#   Local:           ΔAICc = 9.94,   weight = 0.004 — worst
#   Local + MPA:     ΔAICc = 13.15,  weight = 0.001 — worst
#
# Key findings:
# 1. MPA × pressure is the best supported model at transect
#    level (weight = 0.519) — this is a stronger result than
#    at site level where the interaction was only marginal.
#    The interaction strengthens rather than weakens when
#    within-site variation is retained.
#
# 2. Global model competitive (ΔAICc = 1.05, weight = 0.307)
#    — connectivity and environmental terms contribute
#    alongside MPA × pressure at transect resolution.
#
# 3. Local + spatial drops out of the competitive set
#    (ΔAICc = 3.81) — the negative connectivity signal
#    from the site-level analysis does not dominate at
#    transect level once the MPA × pressure interaction
#    is included.
#
# 4. Local + MPA performs poorly (ΔAICc = 13.15) —
#    the additive MPA effect alone is insufficient at
#    transect level; the interaction with pressure is
#    needed to capture the pattern.
#
# Conclusion: the MPA × pressure interaction is the most
# consistent finding across analytical scales for
# corallivores, strengthening from marginal at site level
# to clearly supported at transect level. This supports
# treating it as a genuine biological signal rather than
# an artefact of site-level aggregation.

# --------------------------------------------------------------------------
# ============================================================
#  PART 3 — TRANSECT-LEVEL COUNTS (COMPLEMENTARY ANALYSIS)
# ============================================================

cat("Transects:", nrow(transect_model_data), "\n")
cat("Zero counts:", sum(transect_model_data$transect_coralliv_count == 0), "\n")
cat("Proportion zeros:", round(mean(transect_model_data$transect_coralliv_count == 0), 3), "\n")

summary(transect_model_data$transect_coralliv_count)

ggplot(transect_model_data, aes(x = transect_coralliv_count)) +
  geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
  labs(x = "Total corallivore count per transect", y = "Frequency") +
  theme_bw()

transect_model_data %>%
  group_by(site) %>%
  summarise(mean_count = mean(transect_coralliv_count),
            var_count = var(transect_coralliv_count),
            .groups = "drop") %>%
  ggplot(aes(x = mean_count, y = var_count)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Site mean count", y = "Site variance",
       title = "Mean-variance (red = Poisson expectation)") +
  theme_bw()

# ── Family selection ──────────────────────────────────────────

m_count_poisson <- glmmTMB(
  transect_coralliv_count ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    log_max_dhw_sc +
    (1 | site),
  family = poisson(link = "log"),
  data = transect_model_data
)

m_count_nb2 <- glmmTMB(
  transect_coralliv_count ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    log_max_dhw_sc +
    (1 | site),
  family = nbinom2(link = "log"),
  data = transect_model_data
)

m_count_nb1 <- glmmTMB(
  transect_coralliv_count ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    log_max_dhw_sc +
    (1 | site),
  family = nbinom1(link = "log"),
  data = transect_model_data
)

res_poisson <- simulateResiduals(m_count_poisson, n = 1000)
res_nb2 <- simulateResiduals(m_count_nb2, n = 1000)
res_nb1 <- simulateResiduals(m_count_nb1, n = 1000)

jpeg("dharma_coralliv_count_poisson.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_poisson, main = "DHARMa — Poisson"); dev.off()

jpeg("dharma_coralliv_count_nb2.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_nb2, main = "DHARMa — NB2"); dev.off()

jpeg("dharma_coralliv_count_nb1.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_nb1, main = "DHARMa — NB1"); dev.off()

plot(res_poisson); testDispersion(res_poisson); testZeroInflation(res_poisson)
plot(res_nb2); testDispersion(res_nb2); testZeroInflation(res_nb2)
plot(res_nb1); testDispersion(res_nb1); testZeroInflation(res_nb1)

cat("\n--- Family selection: corallivore count models ---\n")
print(make_aicc_df(list(
  "Poisson" = m_count_poisson,
  "NB2" = m_count_nb2,
  "NB1" = m_count_nb1
)))
# Update decision after running.

# ── Random effect structure ───────────────────────────────────

count_family <- nbinom2(link = "log") # update after family selection

re_c_null <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                       log_settlement_grav_sc +
                       log_chla_sc +
                       log_max_dhw_sc,
                     family = count_family, data = transect_model_data)

re_c_site <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                       log_settlement_grav_sc +
                       log_chla_sc +
                       log_max_dhw_sc +
                       (1 | site),
                     family = count_family, data = transect_model_data)

cat("\n--- RE structure: corallivore count models ---\n")
print(make_aicc_df(list(
  "No RE" = re_c_null,
  "(1 | site)" = re_c_site
)))
# (1 | site) expected to be supported — retain for consistency
# with all other count analyses regardless of delta.

# ── Candidate models ──────────────────────────────────────────

c_m1_hab <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                      (1 | site),
                    family = count_family, data = transect_model_data)

c_m2_hab_press <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                            log_settlement_grav_sc +
                            (1 | site),
                          family = count_family, data = transect_model_data)

c_m3_hab_press_mpa <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                                log_settlement_grav_sc +
                                mpa_status +
                                (1 | site),
                              family = count_family, data = transect_model_data)

c_m4_conn <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                       log_settlement_grav_sc +
                       mpa_status +
                       connectivity_sc +
                       (1 | site),
                     family = count_family, data = transect_model_data)

c_m5_chla <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                       log_settlement_grav_sc +
                       mpa_status +
                       connectivity_sc +
                       log_chla_sc +
                       (1 | site),
                     family = count_family, data = transect_model_data)

c_m6_dhw <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                      log_settlement_grav_sc +
                      mpa_status +
                      connectivity_sc +
                      log_max_dhw_sc +
                      (1 | site),
                    family = count_family, data = transect_model_data)

c_m7_mpa_conn <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                           log_settlement_grav_sc +
                           mpa_status * connectivity_sc +
                           (1 | site),
                         family = count_family, data = transect_model_data)

c_m8_mpa_press <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                            mpa_status * log_settlement_grav_sc +
                            connectivity_sc +
                            (1 | site),
                          family = count_family, data = transect_model_data)

c_m9_conn_press <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                             mpa_status +
                             connectivity_sc * log_settlement_grav_sc +
                             (1 | site),
                           family = count_family, data = transect_model_data)

c_sens_settpop <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                            log_settlement_pop_sc +
                            mpa_status +
                            connectivity_sc +
                            (1 | site),
                          family = count_family, data = transect_model_data)

c_sens_mktgrav <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                            log_market_gravity_sc +
                            mpa_status +
                            connectivity_sc +
                            (1 | site),
                          family = count_family, data = transect_model_data)

model_list_counts <- list(
  "Habitat" = c_m1_hab,
  "Habitat + pressure" = c_m2_hab_press,
  "Habitat + pressure + MPA" = c_m3_hab_press_mpa,
  "Above + connectivity" = c_m4_conn,
  "Above + chla" = c_m5_chla,
  "Above + DHW" = c_m6_dhw,
  "MPA x connectivity" = c_m7_mpa_conn,
  "MPA x pressure" = c_m8_mpa_press,
  "Connectivity x pressure" = c_m9_conn_press,
  "Settlement pop. (sensitivity)" = c_sens_settpop,
  "Market gravity (sensitivity)" = c_sens_mktgrav
)

cat("\n--- AICc: Corallivore count models ---\n")
print(make_aicc_df(model_list_counts))

# ── Results — update after running ───────────────────────────
# Key comparison: do count results converge with biomass?
# Divergence indicates whether predictors operate through
# fish abundance, body size, or both.

# ============================================================
#  SYNTHESIS: CORALLIVORE BIOMASS vs COUNT CONCLUSIONS
# ============================================================
# Update after running all three parts.
# Key questions:
#   1. Is the null model dominant as in the old analysis?
#   2. Do MPA interactions appear for corallivores?
#   3. Does DHW affect corallivore counts?
#   4. How do corallivores differ from browsers?

cat("\n=== CORALLIVORE — Site-level biomass, best model ===\n")
# summary(best model — update after running)

cat("\n=== CORALLIVORE — Transect-level biomass, best model ===\n")
# summary(best model — update after running)

cat("\n=== CORALLIVORE — Count models, top models ===\n")
# summary(best model — update after running)