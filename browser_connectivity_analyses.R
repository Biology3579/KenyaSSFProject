# ============================================================
#  BROWSER FISH BIOMASS — SITE-LEVEL ANALYSIS
#
#  Analytical framework mirrors total biomass (four stages):
#
#  STAGE 1 — Variance partitioning
#             Quantifies unique and shared variance attributable
#             to local ecological, spatial, and environmental
#             process groups. MPA status excluded (governance
#             variable — tested separately in Stage 2).
#
#  STAGE 2 — Hierarchical model comparison
#             Nested model sequence adds process groups
#             progressively. Tests incremental explanatory value
#             of each group beyond the local baseline.
#
#  STAGE 3 — Interaction testing (conditional on Stage 2)
#             Three a priori interactions test whether spatial
#             management and connectivity modify local-scale
#             relationships.
#
#  STAGE 4 — Sensitivity analysis
#             (a) Alternative pressure metrics
#             (b) Transect-level mixed model replication
#
#  Key difference from total biomass:
#    Browser biomass has ~11% zeros at site level and ~43%
#    at transect level. Family selection therefore tests
#    Tweedie in addition to Gaussian log. If Tweedie is
#    selected, varpart is not directly available (vegan::varpart
#    requires lm-compatible models). In that case Stage 1 is
#    approximated using pseudo-R² from glmmTMB models.
#    See family selection decision below.
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

# ── AGGREGATE BROWSER TRANSECT DATA ──────────────────────
browser_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_browser_biomass = sum(
      ifelse(trophic_group == "browsers", tot_wt_g, 0),
      na.rm = TRUE
    ),
    transect_browser_count = sum(
      ifelse(trophic_group == "browsers", number, 0),
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
cat("Sites:",                 n_distinct(browser_transects$site), "\n")
cat("Countries:",             n_distinct(browser_transects$country), "\n")

# ── SITE-LEVEL AGGREGATION ────────────────────────────────────
browser_site_data <- browser_transects %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass = mean(transect_browser_biomass, na.rm = TRUE),
    n_transects  = n(),
    .groups      = "drop"
  ) %>%
  mutate(
    site    = as.factor(site),
    country = as.factor(country)
  )

# ── Zero check ────────────────────────────────────────────────
zeros_site <- mean(browser_site_data$mean_biomass == 0, na.rm = TRUE)
cat("Site-level zero proportion:", round(zeros_site, 3), "\n")
# If > 0: log(x + offset) required for Gaussian
# If substantial (> 0.15): Tweedie warranted

# ── Distribution plots ────────────────────────────────────────
p_raw <- ggplot(browser_site_data, aes(x = mean_biomass)) +
  geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
  labs(x = "Mean browser biomass per site (g)", y = "Frequency",
       title = "Raw") + theme_bw()

# Apply log(x + 0.01) transformation — offset chosen as 1%
# of minimum non-zero value to minimise distortion
browser_site_data <- browser_site_data %>%
  mutate(log_mean_biomass = log(mean_biomass + 0.01))

p_log <- ggplot(browser_site_data, aes(x = log_mean_biomass)) +
  geom_histogram(bins = 25, fill = "#1a9641", colour = "white") +
  labs(x = "log(mean biomass + 0.01)", y = "Frequency",
       title = "Log-transformed") + theme_bw()

gridExtra::grid.arrange(p_raw, p_log, ncol = 2)

# ── Distribution decision ─────────────────────────────────────
# Log-transformed distribution is bimodal — 6 zero sites
# (11%) form an isolated cluster at log(0.01) = -4.6,
# separated from the main distribution by a gap of ~10 log
# units. This pattern is invariant to the choice of offset
# constant: no log transformation can bridge the gap between
# genuine absence and presence because the two states are
# biologically distinct, not simply different points on a
# continuous scale. A tweedie model will liekly be needed.

# ── Box-Cox on non-zero values ────────────────────────────────
browser_nonzero <- browser_site_data %>% filter(mean_biomass > 0)
MASS::boxcox(lm(mean_biomass ~ 1, data = browser_nonzero),
             lambda = seq(-2, 2, 0.1))


# ── Build full browser site-level model dataset ───────────────
# Uses scaled_predictors from total biomass preparation —
# predictors are identical across functional groups.
# MPA status: unordered factor, reference = "none"

browser_model_data <- browser_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_browser_biomass, na.rm = TRUE),
    log_mean_biomass       = log(mean(transect_browser_biomass,
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
    site    = as.factor(site),
    country = as.factor(country),
    mpa_status = factor(mpa_status,
                        levels  = c("none", "low", "medium"),
                        ordered = FALSE)
  )

cat("\nBrowser model data:", nrow(browser_model_data), "sites\n")
cat("Site-level zeros:",   sum(browser_model_data$mean_biomass == 0), "\n")

# ── Verify no NAs ─────────────────────────────────────────────
browser_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                connectivity_sc, mpa_status,
                log_chla_sc, log_max_dhw_sc) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

# ── Check zeros and response ──────────────────────────────────
cat("Zeros in mean_biomass:",       sum(browser_model_data$mean_biomass == 0), "\n")
cat("-Inf in log_mean_biomass:",    sum(is.infinite(browser_model_data$log_mean_biomass)), "\n")
cat("NAs in log_mean_biomass:",     sum(is.na(browser_model_data$log_mean_biomass)), "\n")
cat("\nResponse summary:\n")
print(summary(browser_model_data$log_mean_biomass))



# ============================================================
#  PRE-ANALYSIS: HUMAN PRESSURE METRIC SELECTION
#  Mirrors total biomass procedure.
#  Tests both Gaussian log and Tweedie to match family
#  selection — use Gaussian log here for comparability with
#  total biomass univariate selection.
# ============================================================

browser_press_settgrav <- lm(log_mean_biomass ~ log_settlement_grav_sc,
                             data = browser_model_data)
browser_press_settpop  <- lm(log_mean_biomass ~ log_settlement_pop_sc,
                             data = browser_model_data)
browser_press_mktgrav  <- lm(log_mean_biomass ~ log_market_gravity_sc,
                             data = browser_model_data)

cat("\n--- Pre-analysis: browser pressure metric selection ---\n")
print(make_aicc_df(list(
  "Settlement gravity" = browser_press_settgrav,
  "Settlement pop."    = browser_press_settpop,
  "Market gravity"     = browser_press_mktgrav
)))


# ── Browser pressure metric decision ──────────────────────────
# Univariate AICc comparison of three candidate metrics.
# Identical procedure to total biomass pre-analysis selection.
#
# Results:
#   Market gravity:     AICc = 297.48, weight = 0.552 — best
#   Settlement gravity: ΔAICc = 1.47,  weight = 0.266
#   Settlement pop.:    ΔAICc = 2.22,  weight = 0.182
#
# Gravity metrics are clearly preferred — both within
# ΔAICc < 2, indicating genuine uncertainty in which
# pressure proxy best captures fishing impacts on browsers.
# Market gravity selected as primary metric: it has the
# lowest AICc and greatest weight (more than twice that of the other two).

# This differs from total biomass where settlement gravity was clearly preferred
# (weight = 0.672, ΔAICc > 2 for alternatives) — the
# weaker metric discrimination for browsers likely reflects
# the more complex, management-mediated relationships
# governing this heavily targeted functional group.
#
# All three metrics retained for sensitivity analysis (Stage 4a).

# ============================================================
#  MODEL FAMILY SELECTION
#
#  Browser biomass differs from total biomass:
#    Site-level zeros: ~11% — log(x + 0.01) offsets these
#    Tweedie handles zeros natively without offset
#
#  Three families tested on global model structure:
#    Gaussian log:  lm() on log(mean_biomass + 0.01)
#    Tweedie:       glmmTMB() on mean_biomass, log link
#
#  DHARMa residual simulation used for Tweedie diagnostics.
#  Standard plot() used for Gaussian log diagnostics.
#
#  Decision rule:
#    If Gaussian log diagnostics adequate → use lm() as per
#    total biomass for consistency and varpart compatibility.
#    If Gaussian log fails → use Tweedie; note that varpart
#    is not directly available and Stage 1 will use an
#    approximation (see below).
# ============================================================

# ── Gaussian log (global model) ───────────────────────────────
browser_lm_global <- lm(log_mean_biomass ~ rugosity_sc +
                          log_market_gravity_sc +
                          connectivity_sc +
                          mpa_status +
                          log_chla_sc +
                          log_max_dhw_sc,
                        data = browser_model_data)

par(mfrow = c(2, 2)); plot(browser_lm_global); par(mfrow = c(1, 1))

# ── Tweedie (global model) ────────────────────────────────────
browser_tw_global <- glmmTMB(mean_biomass ~ rugosity_sc +
                               log_market_gravity_sc +
                               connectivity_sc +
                               mpa_status +
                               log_chla_sc +
                               log_max_dhw_sc,
                             family = tweedie(link = "log"),
                             data   = browser_model_data)

browser_tw_res <- simulateResiduals(browser_tw_global, n = 1000)
plot(browser_tw_res)
testZeroInflation(browser_tw_res)
testDispersion(browser_tw_res)

# ============================================================
#  MODEL FAMILY SELECTION
#
#  Browser biomass has 6 zero sites (11%) at the site level.
#  Log-transformed distribution is bimodal — zero sites form
#  an isolated cluster at log(0.01) = -4.6, separated from
#  the main distribution by ~10 log units. No offset constant
#  resolves this as the gap reflects genuine biological
#  absence, not low biomass.
#
#  Two families tested on the global model structure:
#
#  Gaussian log: REJECTED
#    Residuals vs Fitted: strong downward curve driven by
#      zero sites — clear non-linearity.
#    Q-Q: severe lower tail deviation — zero sites fall
#      off the theoretical line entirely.
#    Scale-Location: strong downward trend —
#      heteroscedasticity throughout.
#    Residuals vs Leverage: sites 5, 42, 32 outside Cook's
#      distance — zero sites unduly influential.
#
#  Tweedie (log link): SELECTED
#    Handles zeros natively without offset.
#    DHARMa diagnostics (n = 1000 simulations):
#      KS test:         p = 0.785 — good fit
#      Dispersion:      p = 0.218, ratio = 1.605 — acceptable
#      Zero inflation:  p = 0.984 — ZI component not needed
#      Outlier test:    p = 1.000 — no outliers
#      Residuals vs predicted: no systematic pattern
#
#  Proceed: glmmTMB(family = tweedie(link = "log")) on
#  raw mean_biomass throughout browser site-level analyses.
#
#  Implication for Stage 1:
#  vegan::varpart() requires lm-compatible models and cannot
#  be used with glmmTMB Tweedie. Stage 1 uses a sequential
#  pseudo-R² approximation instead (see below).
# ============================================================

# ============================================================
#  RANDOM EFFECT STRUCTURE ---- make sure to update this or decide on this...
#  Country RE tested as per total biomass procedure.
#  Uses Tweedie on raw mean_biomass — consistent with
#  family selection decision above.
# ============================================================

browser_re_null <- glmmTMB(mean_biomass ~ rugosity_sc +
                             log_market_gravity_sc +
                             log_chla_sc +
                             log_max_dhw_sc,
                           family = tweedie(link = "log"),
                           data   = browser_model_data)

browser_re_country <- glmmTMB(mean_biomass ~ rugosity_sc +
                                log_market_gravity_sc +
                                log_chla_sc +
                                log_max_dhw_sc +
                                (1 | country),
                              family = tweedie(link = "log"),
                              data   = browser_model_data)

cat("\n--- Browser RE structure comparison ---\n")
print(make_aicc_df(list(
  "No RE"         = browser_re_null,
  "(1 | country)" = browser_re_country
)))

# ── Random effect structure ───────────────────────────────────
# Country RE tested as per total biomass procedure.
#
# Results:
#   (1 | country): AICc = 838.77, weight = 0.510
#   No RE:         ΔAICc = 0.08,  weight = 0.490
#
# Weights are split almost exactly 50/50 — the two models
# are statistically indistinguishable. Country RE adds no
# meaningful explanatory value, consistent with the total
# biomass result. Proceed without country RE throughout
# all browser site-level analyses.


# ============================================================
#  STAGE 1 — VARIANCE PARTITIONING
#
#  Tweedie selected for browser biomass — vegan::varpart()
#  not available as it requires lm-compatible models.
#
#  Sequential pseudo-R² decomposition used instead:
#  Fit models adding each process group in turn, compute
#  Nagelkerke R² at each step using the performance package.
#  Unique fractions approximated as the increment in R²
#  when each group is added while conditioning on others.
#
#  This is not equivalent to formal variance partitioning —
#  shared fractions cannot be cleanly decomposed and negative
#  unique fractions are not estimable. Results are reported
#  as approximate contributions and interpreted alongside
#  Stage 2 hierarchical model comparison rather than as
#  standalone variance decomposition.
#
#  MPA excluded — governance variable, tested in Stage 2.
#  Three process groups as per total biomass:
#    Local:       rugosity + settlement gravity
#    Spatial:     connectivity
#    Environmental: chla + DHW
# ============================================================

# ── Null model ────────────────────────────────────────────────
vp_b_null <- glmmTMB(mean_biomass ~ 1,
                     family = tweedie(link = "log"),
                     data   = browser_model_data)

# ── Local only ────────────────────────────────────────────────
vp_b_local <- glmmTMB(mean_biomass ~ rugosity_sc +
                        log_market_gravity_sc,
                      family = tweedie(link = "log"),
                      data   = browser_model_data)

# ── Spatial only ──────────────────────────────────────────────
vp_b_spatial <- glmmTMB(mean_biomass ~ connectivity_sc,
                        family = tweedie(link = "log"),
                        data   = browser_model_data)

# ── Environmental only ────────────────────────────────────────
vp_b_environ <- glmmTMB(mean_biomass ~ log_chla_sc +
                          log_max_dhw_sc,
                        family = tweedie(link = "log"),
                        data   = browser_model_data)

# ── All groups combined (no MPA) ──────────────────────────────
vp_b_all <- glmmTMB(mean_biomass ~ rugosity_sc +
                      log_market_gravity_sc +
                      connectivity_sc +
                      log_chla_sc +
                      log_max_dhw_sc,
                    family = tweedie(link = "log"),
                    data   = browser_model_data)

# ── Unique fractions: condition each group on the others ──────
# Local | spatial + environmental
vp_b_local_unique <- glmmTMB(mean_biomass ~ rugosity_sc +
                               log_market_gravity_sc +
                               connectivity_sc +
                               log_chla_sc +
                               log_max_dhw_sc,
                             family = tweedie(link = "log"),
                             data   = browser_model_data)

vp_b_no_local <- glmmTMB(mean_biomass ~ connectivity_sc +
                           log_chla_sc +
                           log_max_dhw_sc,
                         family = tweedie(link = "log"),
                         data   = browser_model_data)

vp_b_no_spatial <- glmmTMB(mean_biomass ~ rugosity_sc +
                             log_market_gravity_sc +
                             log_chla_sc +
                             log_max_dhw_sc,
                           family = tweedie(link = "log"),
                           data   = browser_model_data)

vp_b_no_environ <- glmmTMB(mean_biomass ~ rugosity_sc +
                             log_market_gravity_sc +
                             connectivity_sc,
                           family = tweedie(link = "log"),
                           data   = browser_model_data)

# ── Extract McFadden pseudo-R² ────────────────────────────────
# McFadden R² = 1 - (logLik(model) / logLik(null))
# Comparable across Tweedie models with same response and data.

get_mcfadden <- function(model, null_model) {
  round(1 - (as.numeric(logLik(model)) / 
               as.numeric(logLik(null_model))), 3)
}

r2_all     <- get_mcfadden(vp_b_all,        vp_b_null)
r2_noloc   <- get_mcfadden(vp_b_no_local,   vp_b_null)
r2_nospat  <- get_mcfadden(vp_b_no_spatial, vp_b_null)
r2_noenv   <- get_mcfadden(vp_b_no_environ, vp_b_null)

# Unique fractions = R²(all) - R²(all without that group)
unique_local   <- round(r2_all - r2_noloc,  3)
unique_spatial <- round(r2_all - r2_nospat, 3)
unique_environ <- round(r2_all - r2_noenv,  3)

cat("\n--- Stage 1: Browser pseudo-R² variance decomposition ---\n")
cat("Total McFadden R² (all groups):", r2_all,        "\n")
cat("Unique local:                  ", unique_local,   "\n")
cat("Unique spatial:                ", unique_spatial, "\n")
cat("Unique environmental:          ", unique_environ, "\n")
cat("\nNOTE: McFadden R² used in place of adjusted R² —",
    "\nTweedie family precludes vegan::varpart().",
    "\nUnique fractions approximate only — shared variance",
    "\nnot decomposable with this method.\n")

cat("\nLRT — local unique fraction:\n")
print(anova(vp_b_no_local, vp_b_all))

cat("\nLRT — spatial unique fraction:\n")
print(anova(vp_b_no_spatial, vp_b_all))

cat("\nLRT — environmental unique fraction:\n")
print(anova(vp_b_no_environ, vp_b_all))

# ── Stage 1 results ───────────────────────────────────────────
# McFadden R² fractions (approximate — but not comparable to
# adjusted R² from total biomass varpart):
#   Total R²:             0.018
#   Unique local:         0.010
#   Unique spatial:       0.006
#   Unique environmental: 0.002
#
# LRT significance of unique fractions:
#   Local:         χ²(2) = 7.55, p = 0.023 * — significant
#   Spatial:       χ²(1) = 4.93, p = 0.026 * — significant
#   Environmental: χ²(2) = 0.91, p = 0.635   — not significant
#
# KEY CONTRAST WITH TOTAL BIOMASS:
#   Total biomass: local significant, spatial p = 0.904
#   Browser biomass: BOTH local and spatial significant
#
# Connectivity explains a significant unique fraction of
# browser biomass independently of local and environmental
# processes. This suggests browsers as a functional group
# are more sensitive to larval network position than the
# total fish community — consistent with browsers being
# heavily targeted by fishing and dependent on external
# recruitment to maintain populations at exploited sites.
#
# Environmental context not significant for browsers —
# consistent with total biomass result (p = 0.077 there,
# p = 0.635 here).
#
# NOTE: McFadden R² magnitudes are not directly comparable
# to adjusted R² from total biomass OLS varpart. LRT
# significance tests are the primary inferential tool here.
# Magnitude comparisons between functional groups should
# be made cautiously and only qualitatively.

# ============================================================
#  STAGE 2 — HIERARCHICAL MODEL COMPARISON
#
#  Identical nested sequence to total biomass.
#  All models fitted using glmmTMB(tweedie) on raw
#  mean_biomass — consistent with family selection.
#
#  AICc comparison uses make_aicc_df() as per total biomass.
#  R² increments use McFadden pseudo-R² over null model
#  (not adjusted R² — Tweedie precludes this).
#  Delta_R2 computed relative to Local baseline as per
#  total biomass procedure.
# ============================================================

# ── Null model ────────────────────────────────────────────────
b_null <- glmmTMB(mean_biomass ~ 1,
                  family = tweedie(link = "log"),
                  data   = browser_model_data)

# ── Local ecological baseline ─────────────────────────────────
b_local <- glmmTMB(mean_biomass ~ rugosity_sc +
                     log_market_gravity_sc,
                   family = tweedie(link = "log"),
                   data   = browser_model_data)

# ── Local + chla ──────────────────────────────────────────────
b_local_chla <- glmmTMB(mean_biomass ~ rugosity_sc +
                          log_market_gravity_sc +
                          log_chla_sc,
                        family = tweedie(link = "log"),
                        data   = browser_model_data)

# ── Local + DHW ───────────────────────────────────────────────
b_local_dhw <- glmmTMB(mean_biomass ~ rugosity_sc +
                         log_market_gravity_sc +
                         log_max_dhw_sc,
                       family = tweedie(link = "log"),
                       data   = browser_model_data)

# ── Local + environmental context ────────────────────────────
b_local_env <- glmmTMB(mean_biomass ~ rugosity_sc +
                         log_market_gravity_sc +
                         log_chla_sc +
                         log_max_dhw_sc,
                       family = tweedie(link = "log"),
                       data   = browser_model_data)

# ── Local + spatial ───────────────────────────────────────────
b_local_spatial <- glmmTMB(mean_biomass ~ rugosity_sc +
                             log_market_gravity_sc +
                             connectivity_sc,
                           family = tweedie(link = "log"),
                           data   = browser_model_data)

# ── Local + MPA (governance test) ────────────────────────────
b_local_mpa <- glmmTMB(mean_biomass ~ rugosity_sc +
                         log_market_gravity_sc +
                         mpa_status,
                       family = tweedie(link = "log"),
                       data   = browser_model_data)

# ── Global model ──────────────────────────────────────────────
b_global <- glmmTMB(mean_biomass ~ rugosity_sc +
                      log_market_gravity_sc +
                      connectivity_sc +
                      mpa_status +
                      log_chla_sc +
                      log_max_dhw_sc,
                    family = tweedie(link = "log"),
                    data   = browser_model_data)

# ── Model list ────────────────────────────────────────────────
browser_model_list <- list(
  "Null"            = b_null,
  "Local"           = b_local,
  "Local + chla"    = b_local_chla,
  "Local + DHW"     = b_local_dhw,
  "Local + env"     = b_local_env,
  "Local + MPA"     = b_local_mpa,
  "Local + spatial" = b_local_spatial,
  "Global"          = b_global
)

# ── AICc ranked table ─────────────────────────────────────────
cat("\n--- Stage 2: Browser model comparison (AICc ranked) ---\n")
print(make_aicc_df(browser_model_list))

# ── McFadden R² increments over local baseline ────────────────
cat("\n--- Stage 2: Browser variance explained ---\n")

local_loglik_b <- as.numeric(logLik(b_local))
null_loglik_b  <- as.numeric(logLik(b_null))

browser_model_list %>%
  imap_dfr(~ tibble(
    Model      = .y,
    McF_R2     = round(1 - (as.numeric(logLik(.x)) /
                              null_loglik_b), 3)
  )) %>%
  mutate(
    Delta_R2 = round(McF_R2 - (1 - (local_loglik_b /
                                      null_loglik_b)), 3),
    Delta_R2 = ifelse(Model %in% c("Null", "Local"), NA, Delta_R2)
  ) %>%
  arrange(Model == "Null",
          desc(Model == "Local"),
          desc(McF_R2)) %>%
  print()

# ── Stage 2 results ───────────────────────────────────────────
# Browser biomass shows a fundamentally different structuring
# pattern from total biomass.
#
# AICc ranking:
#   Local + MPA:     AICc = 828.87, weight = 0.693 — best
#   Local + spatial: ΔAICc = 3.32,  weight = 0.132 — moderate support
#   Global:          ΔAICc = 3.46,  weight = 0.123 — moderate support
#   Local:           ΔAICc = 6.84,  weight = 0.023 — not supported
#   Local + DHW:     ΔAICc = 8.50,  weight = 0.010 — not supported
#   Local + chla:    ΔAICc = 8.68,  weight = 0.009 — not supported
#   Null:            ΔAICc = 9.71,  weight = 0.005 — worst
#   Local + env:     ΔAICc = 9.99,  weight = 0.005 — worst
#
# McFadden R² increments over Local baseline:
#   MPA:         +0.015 — largest single contribution
#   Spatial:     +0.007 — second largest
#   Environment: +0.003 — negligible
#   Global:      +0.021 — driven by MPA and spatial
#
# KEY CONTRASTS WITH TOTAL BIOMASS:
#
# 1. MPA status is the strongest predictor for browsers
#    (best model weight = 0.693) — not supported for total
#    biomass (ΔAICc = 3.52).
#
# 2. Connectivity adds moderate support for browsers
#    (ΔAICc = 3.32) — made models worse for total biomass
#    (ΔR² = -0.011). Consistent with Stage 1 LRT
#    (spatial p = 0.026).
#
# 3. Local baseline barely improves on null for browsers
#    (ΔAICc = 0.87 vs null) — for total biomass Local was
#    strongly supported (ΔAICc = 8.49 vs null).
#
# 4. Environmental context adds nothing for browsers
#    (ΔR² = +0.003) — consistent with Stage 1 (p = 0.635).
#
# Interpretation: browser biomass is structured primarily by
# governance and spatial processes rather than local ecological
# conditions. Protection status is the dominant driver —
# consistent with browsers being a heavily targeted, mobile
# functional group whose local abundance depends on management
# regime and regional replenishment rather than habitat
# carrying capacity alone.

summary(b_local_mpa)

# ── Coefficients: best supported model (Local + MPA) ─────────
#
# Rugosity:          β = +0.513, p = 0.003 **
#   Positive and significant — structurally complex reefs
#   support higher browser biomass. Effect larger than total
#   biomass (β = +0.209), suggesting browsers are more
#   habitat-dependent than the total community.
#
# Market gravity:    β = +0.178, p = 0.322 — not significant
#   No detectable independent pressure effect once MPA
#   status is included. Browsers may be depleted across
#   the full pressure gradient at unprotected sites,
#   leaving no residual pressure signal.
#
# MPA low:           β = +0.136, p = 0.793 — not significant
#   Low protection indistinguishable from none — likely
#   reflects insufficient enforcement and small sample
#   (n = 7 low-protection sites).
#
# MPA medium:        β = +1.176, p < 0.001 ***
#   Highly significant. Medium-protection sites have
#   approximately 3.2x higher browser biomass than
#   unprotected sites (e^1.176 = 3.24) at mean rugosity
#   and pressure. The none→medium contrast is the dominant
#   signal in the browser analysis.
#
# Key interpretation:
#   Browser biomass is primarily structured by protection
#   regime and habitat complexity. The significant medium
#   MPA effect — absent for total biomass — indicates that
#   browsers are particularly sensitive to management
#   intervention, consistent with their status as a heavily
#   targeted functional group. The non-significant pressure
#   main effect suggests MPAs mediate the fishing pressure
#   signal rather than pressure operating independently.

# ============================================================
#  STAGE 3 — INTERACTION TESTING
#
#  Reference model: Local + MPA (best supported in Stage 2,
#  weight = 0.601). Interactions tested against this baseline.
#
#  Gate check: Global outperforms Local by ΔAICc = 7.70 —
#  well above the ΔAICc > 2 threshold. Spatial and governance
#  terms are clearly supported in additive models, making
#  interaction testing both warranted and meaningful.
#
#  This contrasts with total biomass where the gate check
#  failed (Global worse than Local) — for browsers, Stage 3
#  is genuinely interpretable.
# ============================================================

# ── Gate check ────────────────────────────────────────────────
delta_b_global_vs_local <- AICc(b_local) - AICc(b_global)
cat("\nΔAICc (Local vs Global — browsers):",
    round(delta_b_global_vs_local, 2), "\n")

# ── Gate check result ─────────────────────────────────────────
# ΔAICc (Local vs Global) = 6.84
# Global outperforms Local by 6.84 AICc units — well above
# the ΔAICc > 2 threshold. Governance terms add meaningful
# explanatory value beyond the local baseline.
# Stage 3 interactions are warranted and interpretable.
#
# Contrast with total biomass (ΔAICc = -3.88) where global
# was worse than local — interactions were untestable there.
# For browsers, the gate check passes clearly.

# ── Hypothesis 1: MPA effectiveness depends on larval supply ──
# Well-connected protected sites should recover faster and
# maintain higher browser biomass than isolated protected sites
b_int_mpa_conn <- glmmTMB(mean_biomass ~ rugosity_sc +
                            log_market_gravity_sc +
                            mpa_status * connectivity_sc,
                          family = tweedie(link = "log"),
                          data   = browser_model_data)

# ── Hypothesis 2: MPA effectiveness depends on fishing pressure
# MPAs should only be detectable where external pressure is
# low enough that protection translates into biomass difference
b_int_mpa_press <- glmmTMB(mean_biomass ~ rugosity_sc +
                             mpa_status * log_market_gravity_sc,
                           family = tweedie(link = "log"),
                           data   = browser_model_data)

# ── Hypothesis 3: connectivity buffers fishing pressure ───────
# Well-connected sites should be more resilient to exploitation
# through sustained larval replenishment offsetting mortality
b_int_conn_press <- glmmTMB(mean_biomass ~ rugosity_sc +
                              log_market_gravity_sc +
                              connectivity_sc * log_market_gravity_sc,
                            family = tweedie(link = "log"),
                            data   = browser_model_data)

# ── Interaction candidate set ─────────────────────────────────
# Reference: Local + MPA (best supported additive model)
browser_interactions <- list(
  "Local + MPA (additive)"  = b_local_mpa,
  "MPA × connectivity"      = b_int_mpa_conn,
  "MPA × pressure"          = b_int_mpa_press,
  "Connectivity × pressure" = b_int_conn_press
)

cat("\n--- Stage 3: Browser interaction comparison ---\n")
print(make_aicc_df(browser_interactions))

# ── AICc summary ─────────────────────────────────────────────
# Local + MPA (additive):  AICc = 828.87, weight = 0.688 — best
# MPA × connectivity:      ΔAICc = 2.98,  weight = 0.155 — not supported
# MPA × pressure:          ΔAICc = 3.47,  weight = 0.121 — not supported
# Connectivity × pressure: ΔAICc = 5.92,  weight = 0.036 — not supported
#
# No interaction is supported — the additive Local + MPA
# model is the best supported model (weight = 0.688) and
# all interaction models perform worse. ΔAICc > 2 for all
# interactions relative to the additive baseline.
# Stage 3 interactions do not improve on the additive model
# and are not carried forward.
#
# Proceed with Local + MPA (b_local_mpa) as the final
# best-supported model for browser site-level biomass.

# ── Marginal effect plots — Local + MPA (best model) ─────────

# ── Rugosity effect (at mean pressure, reference MPA = none) ──
rug_grid <- data.frame(
  rugosity_sc           = seq(min(browser_model_data$rugosity_sc),
                              max(browser_model_data$rugosity_sc),
                              length.out = 100),
  log_market_gravity_sc = 0,
  mpa_status            = factor("none", levels = c("none", "low", "medium"))
)

rug_pred <- predict(b_local_mpa,
                    newdata = rug_grid,
                    se.fit  = TRUE,
                    type    = "response",
                    re.form = NA)

rug_grid$fit <- rug_pred$fit
rug_grid$lwr <- rug_pred$fit - 1.96 * rug_pred$se.fit
rug_grid$upr <- rug_pred$fit + 1.96 * rug_pred$se.fit

ggplot(rug_grid, aes(x = rugosity_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "#2c7bb6", alpha = 0.2) +
  geom_line(colour = "#2c7bb6", linewidth = 1.1) +
  labs(x = "Rugosity (standardised)", y = "Browser biomass (g)") +
  theme_bw(base_size = 13) +
  theme(axis.title = element_text(face = "bold"))

# ── MPA status effect (at mean rugosity and pressure) ─────────
mpa_grid <- data.frame(
  mpa_status            = factor(c("none", "low", "medium"),
                                 levels = c("none", "low", "medium")),
  rugosity_sc           = 0,
  log_market_gravity_sc = 0
)

mpa_pred <- predict(b_local_mpa,
                    newdata = mpa_grid,
                    se.fit  = TRUE,
                    type    = "response",
                    re.form = NA)

mpa_grid$fit <- mpa_pred$fit
mpa_grid$lwr <- mpa_pred$fit - 1.96 * mpa_pred$se.fit
mpa_grid$upr <- mpa_pred$fit + 1.96 * mpa_pred$se.fit

ggplot(mpa_grid, aes(x = mpa_status, y = fit)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  colour = "#0072B2", size = 0.8) +
  labs(x = "MPA status", y = "Browser biomass (g)") +
  theme_bw(base_size = 13) +
  theme(axis.title = element_text(face = "bold"))

# ── Interpretation ────────────────────────────────────────────
# No interaction is supported — all three interaction models
# perform worse than the additive Local + MPA baseline
# (ΔAICc > 2 for all). The additive model is retained as
# the final best-supported model for browser site-level
# biomass.
#
# MPA × connectivity (ΔAICc = 2.98): the most biologically
#   motivated interaction — MPAs more effective when
#   well-connected — is not supported. MPA and connectivity
#   appear to operate as independent additive processes
#   rather than synergistically.
#
# MPA × pressure (ΔAICc = 3.47): no evidence that MPA
#   effectiveness depends on fishing intensity at the
#   site level when market gravity is used as the pressure
#   metric. Note this interaction was strongly supported
#   in an earlier version of this analysis using settlement
#   gravity — the metric-dependence of this result suggests
#   caution in interpreting any pressure × MPA signal.
#
# Connectivity × pressure (ΔAICc = 5.92): no support for
#   connectivity buffering fishing pressure through larval
#   replenishment at the site level.
#
# Final model: Local + MPA (b_local_mpa)
#   Rugosity and MPA medium status are the two significant
#   drivers of browser biomass. Their effects are additive
#   and do not depend on spatial context or pressure regime.

# ============================================================
#  STAGE 4 — SENSITIVITY ANALYSIS
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
# Best supported model from Stage 3 is MPA × pressure.
# Sensitivity models mirror this structure — only pressure
# metric swapped. MPA × pressure interaction retained to
# test whether the interaction finding is robust to
# pressure proxy choice.

b_sens_settgrav <- glmmTMB(mean_biomass ~ rugosity_sc +
                             mpa_status * log_settlement_grav_sc,
                           family = tweedie(link = "log"),
                           data   = browser_model_data)

b_sens_settpop <- glmmTMB(mean_biomass ~ rugosity_sc +
                            mpa_status * log_settlement_pop_sc,
                          family = tweedie(link = "log"),
                          data   = browser_model_data)

cat("\n--- Stage 4a: Browser sensitivity — pressure metrics ---\n")
cat("Settlement population:\n")
print(round(summary(b_sens_settgrav)$coefficients$cond, 4))
cat("\nMarket gravity:\n")
print(round(summary(b_sens_settpop)$coefficients$cond, 4))

# ── Stage 4a results ──────────────────────────────────────────
# Sensitivity analysis tests MPA × pressure interaction
# structure with two alternative pressure metrics.
# Primary metric (market gravity) uses additive Local + MPA
# as best model — sensitivity models test whether interactions
# emerge under alternative metrics.
#
# Rugosity: significant and positive across all three metrics
#   Market gravity (primary): β = +0.513, p < 0.001
#   Settlement gravity:       β = +0.581, p < 0.001
#   Settlement pop.:          β = +0.396, p = 0.029
#   Most robust predictor — unaffected by metric choice.
#
# MPA medium main effect:
#   Market gravity (primary): β = +1.176, p < 0.001 ✓
#   Settlement gravity:       β = +1.942, p < 0.001 ✓
#   Settlement pop.:          β = +0.838, p = 0.073 — marginal
#   Direction consistent but weakens with settlement pop.
#
# MPA medium × pressure interaction:
#   Market gravity (primary): not tested — additive model best
#   Settlement gravity:       β = +1.958, p < 0.001 ✓
#   Settlement pop.:          β = -0.695, p = 0.145 — not significant,
#                             direction reverses
#
# The MPA × pressure interaction emerges only with settlement
# gravity — it is absent with market gravity (primary metric,
# no interaction supported in Stage 3) and with settlement
# pop. (non-significant, direction reverses). This metric-
# dependence confirms the decision to use the additive model
# as the primary result.
#
# Conservative conclusion: the positive effect of medium
# protection on browser biomass is robust across all three
# metrics. The MPA × pressure interaction is not a general
# finding — it is specific to settlement gravity and should
# not be reported as a primary result. The additive MPA
# medium effect is the reliable biological signal.

# ── (b) Transect-level replication ───────────────────────────
# 43% zeros at transect level — Tweedie required.
# ZI Tweedie tested: if zero inflation test significant
# AND AICc improves by > 2, ZI Tweedie adopted.
# Transect sequence mirrors Stage 2 hierarchical structure
# plus Stage 3 best model (MPA × pressure).

browser_transect_data <- browser_transects %>%
  left_join(final_predictors, by = "site") %>%
  mutate(log_browser_biomass = log(transect_browser_biomass + 0.01))

cat("\nTransect zeros:",
    sum(browser_transect_data$transect_browser_biomass == 0),
    "/", nrow(browser_transect_data),
    "(", round(mean(browser_transect_data$transect_browser_biomass == 0), 3), ")\n")

# ── Transect family selection ─────────────────────────────────
b_trans_tw <- glmmTMB(transect_browser_biomass ~ rugosity_sc +
                        log_market_gravity_sc +
                        connectivity_sc +
                        mpa_status +
                        log_chla_sc +
                        log_max_dhw_sc +
                        (1 | site),
                      family = tweedie(link = "log"),
                      data   = browser_transect_data)

b_trans_tw_zi <- glmmTMB(transect_browser_biomass ~ rugosity_sc +
                           log_market_gravity_sc +
                           connectivity_sc +
                           mpa_status +
                           log_chla_sc +
                           log_max_dhw_sc +
                           (1 | site),
                         family    = tweedie(link = "log"),
                         ziformula = ~1,
                         data      = browser_transect_data)

b_trans_res    <- simulateResiduals(b_trans_tw,    n = 1000)
b_trans_res_zi <- simulateResiduals(b_trans_tw_zi, n = 1000)

plot(b_trans_res);    testZeroInflation(b_trans_res)
plot(b_trans_res_zi); testZeroInflation(b_trans_res_zi)

cat("\n--- Transect family selection ---\n")
print(make_aicc_df(list(
  "Tweedie"    = b_trans_tw,
  "ZI Tweedie" = b_trans_tw_zi
)))

# ── Transect family selection decision ────────────────────────
# Standard Tweedie selected over ZI Tweedie.
#
# Zero inflation tests:
#   Tweedie:    ratio = 0.964, p = 0.680 — not significant
#   ZI Tweedie: ratio = 0.971, p = 0.734 — not significant
#   Neither model shows evidence of zero inflation —
#   standard Tweedie handles 43% zeros adequately.
#
# AICc:
#   ZI Tweedie: AICc = 2656.11, weight = 0.633
#   Tweedie:    ΔAICc = 1.09,   weight = 0.367
#   ZI Tweedie marginally preferred by AICc but within
#   ΔAICc < 2 — genuine uncertainty.
#
# Proceed: standard Tweedie + (1|site) throughout
# transect-level browser analyses.

# ── Transect hierarchical sequence ───────────────────────────
# Mirrors Stage 2 sequence plus Stage 3 interaction model.
# Key test: does the additive Local + MPA result from the
# site-level analysis replicate at transect level?

b_trans_null <- glmmTMB(transect_browser_biomass ~ 1 +
                          (1 | site),
                        family = tweedie(link = "log"),
                        data   = browser_transect_data)

b_trans_local <- glmmTMB(transect_browser_biomass ~ rugosity_sc +
                           log_market_gravity_sc +
                           (1 | site),
                         family = tweedie(link = "log"),
                         data   = browser_transect_data)

b_trans_local_env <- glmmTMB(transect_browser_biomass ~ rugosity_sc +
                               log_market_gravity_sc +
                               log_chla_sc +
                               log_max_dhw_sc +
                               (1 | site),
                             family = tweedie(link = "log"),
                             data   = browser_transect_data)

b_trans_local_spatial <- glmmTMB(transect_browser_biomass ~ rugosity_sc +
                                   log_market_gravity_sc +
                                   connectivity_sc +
                                   (1 | site),
                                 family = tweedie(link = "log"),
                                 data   = browser_transect_data)

b_trans_local_mpa <- glmmTMB(transect_browser_biomass ~ rugosity_sc +
                               log_market_gravity_sc +
                               mpa_status +
                               (1 | site),
                             family = tweedie(link = "log"),
                             data   = browser_transect_data)

b_trans_global <- glmmTMB(transect_browser_biomass ~ rugosity_sc +
                            log_market_gravity_sc +
                            connectivity_sc +
                            mpa_status +
                            log_chla_sc +
                            log_max_dhw_sc +
                            (1 | site),
                          family = tweedie(link = "log"),
                          data   = browser_transect_data)

# ── Stage 3 best model at transect level ─────────────────────
b_trans_mpa_press <- glmmTMB(transect_browser_biomass ~ rugosity_sc +
                               mpa_status * log_market_gravity_sc +
                               (1 | site),
                             family = tweedie(link = "log"),
                             data   = browser_transect_data)

browser_transect_list <- list(
  "Null"            = b_trans_null,
  "Local"           = b_trans_local,
  "Local + env"     = b_trans_local_env,
  "Local + spatial" = b_trans_local_spatial,
  "Local + MPA"     = b_trans_local_mpa,
  "Global"          = b_trans_global,
  "MPA × pressure"  = b_trans_mpa_press
)

cat("\n--- Stage 4b: Browser transect-level sensitivity ---\n")
print(make_aicc_df(browser_transect_list))

# ── Stage 4b results ──────────────────────────────────────────
# Results are consistent with the site-level finding —
# Local + MPA is again the top-ranked model — but model
# selection uncertainty is high at transect level.
#
# Model ranking (transect-level):
#   Local + MPA:     AICc = 2655.98, weight = 0.238 — best
#   Local + spatial: ΔAICc = 0.15,   weight = 0.221 — equivalent
#   Local:           ΔAICc = 0.33,   weight = 0.202 — equivalent
#   Local + env:     ΔAICc = 1.05,   weight = 0.140 — equivalent
#   Global:          ΔAICc = 1.21,   weight = 0.130 — equivalent
#   MPA × pressure:  ΔAICc = 2.86,   weight = 0.057 — not supported
#   Null:            ΔAICc = 5.72,   weight = 0.014 — worst
#
# Five models fall within ΔAICc < 2 — no single model is
# clearly preferred at transect level. This is a weaker
# result than the site-level analysis (Local + MPA weight
# = 0.693) and reflects the additional within-site noise
# retained at transect resolution diluting between-site
# predictor signals.
#
# Key points:
# 1. Local + MPA top-ranked at both levels — consistent
#    qualitative conclusion across analytical scales.
# 2. MPA × pressure not supported at transect level
#    (ΔAICc = 2.86) — further confirms the interaction
#    is not a robust finding.
# 3. High model uncertainty at transect level is expected:
#    site-level predictors (rugosity, MPA status, market
#    gravity) vary between sites, not within them, so
#    their signal is diluted when within-site transect
#    variation is retained.
#
# Conclusion: the site-level result (Local + MPA best
# supported) is qualitatively replicated at transect level.
# The sensitivity analysis supports the primary finding
# despite higher uncertainty at finer resolution.


# ------------------------------------------------------------------------------

# ============================================================
#  PART 3 — TRANSECT-LEVEL COUNTS (COMPLEMENTARY ANALYSIS)
# ============================================================

cat("Transects:", nrow(transect_model_data), "\n")
cat("Zero counts:", sum(transect_model_data$transect_browser_count == 0), "\n")
cat("Proportion zeros:", round(mean(transect_model_data$transect_browser_count == 0), 3), "\n")

summary(transect_model_data$transect_browser_count)

ggplot(transect_model_data, aes(x = transect_browser_count)) +
  geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
  labs(x = "Total browser count per transect", y = "Frequency") +
  theme_bw()

transect_model_data %>%
  group_by(site) %>%
  summarise(mean_count = mean(transect_browser_count),
            var_count = var(transect_browser_count),
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
  transect_browser_count ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc +
    log_max_dhw_sc +
    (1 | site),
  family = poisson(link = "log"),
  data = transect_model_data
)

m_count_nb2 <- glmmTMB(
  transect_browser_count ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc +
    log_max_dhw_sc +
    (1 | site),
  family = nbinom2(link = "log"),
  data = transect_model_data
)

m_count_nb1 <- glmmTMB(
  transect_browser_count ~ rugosity_sc +
    log_market_gravity_sc +
    log_chla_sc +
    log_max_dhw_sc +
    (1 | site),
  family = nbinom1(link = "log"),
  data = transect_model_data
)

res_poisson <- simulateResiduals(m_count_poisson, n = 1000)
res_nb2 <- simulateResiduals(m_count_nb2, n = 1000)
res_nb1 <- simulateResiduals(m_count_nb1, n = 1000)

jpeg("dharma_browser_count_poisson.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_poisson, main = "DHARMa — Poisson"); dev.off()

jpeg("dharma_browser_count_nb2.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_nb2, main = "DHARMa — NB2"); dev.off()

jpeg("dharma_browser_count_nb1.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_nb1, main = "DHARMa — NB1"); dev.off()

plot(res_poisson); testDispersion(res_poisson); testZeroInflation(res_poisson)
plot(res_nb2); testDispersion(res_nb2); testZeroInflation(res_nb2)
plot(res_nb1); testDispersion(res_nb1); testZeroInflation(res_nb1)

cat("\n--- Family selection: browser count models ---\n")
print(make_aicc_df(list(
  "Poisson" = m_count_poisson,
  "NB2" = m_count_nb2,
  "NB1" = m_count_nb1
)))

# ── Random effect structure ───────────────────────────────────

count_family <- nbinom1(link = "log") # update after family selection

re_c_null <- glmmTMB(transect_browser_count ~ rugosity_sc +
                       log_market_gravity_sc +
                       log_chla_sc +
                       log_max_dhw_sc,
                     family = count_family, data = transect_model_data)

re_c_site <- glmmTMB(transect_browser_count ~ rugosity_sc +
                       log_market_gravity_sc +
                       log_chla_sc +
                       log_max_dhw_sc +
                       (1 | site),
                     family = count_family, data = transect_model_data)

cat("\n--- RE structure: browser count models ---\n")
print(make_aicc_df(list(
  "No RE" = re_c_null,
  "(1 | site)" = re_c_site
)))
# (1 | site) retained for consistency, despite the no RE delta being slightly higher 
#  Although this does say something about abundance...

# ── Candidate models ──────────────────────────────────────────

c_m1_hab <- glmmTMB(transect_browser_count ~ rugosity_sc +
                      (1 | site),
                    family = count_family, data = transect_model_data)

c_m2_hab_press <- glmmTMB(transect_browser_count ~ rugosity_sc +
                            log_market_gravity_sc +
                            (1 | site),
                          family = count_family, data = transect_model_data)

c_m3_hab_press_mpa <- glmmTMB(transect_browser_count ~ rugosity_sc +
                                log_market_gravity_sc +
                                mpa_status +
                                (1 | site),
                              family = count_family, data = transect_model_data)

c_m4_conn <- glmmTMB(transect_browser_count ~ rugosity_sc +
                       log_market_gravity_sc +
                       mpa_status +
                       connectivity_sc +
                       (1 | site),
                     family = count_family, data = transect_model_data)

c_m5_chla <- glmmTMB(transect_browser_count ~ rugosity_sc +
                       log_market_gravity_sc +
                       mpa_status +
                       connectivity_sc +
                       log_chla_sc +
                       (1 | site),
                     family = count_family, data = transect_model_data)

c_m6_dhw <- glmmTMB(transect_browser_count ~ rugosity_sc +
                      log_market_gravity_sc +
                      mpa_status +
                      connectivity_sc +
                      log_max_dhw_sc +
                      (1 | site),
                    family = count_family, data = transect_model_data)

c_m7_mpa_conn <- glmmTMB(transect_browser_count ~ rugosity_sc +
                           log_market_gravity_sc +
                           mpa_status * connectivity_sc +
                           (1 | site),
                         family = count_family, data = transect_model_data)

c_m8_mpa_press <- glmmTMB(transect_browser_count ~ rugosity_sc +
                            mpa_status * log_market_gravity_sc +
                            connectivity_sc +
                            (1 | site),
                          family = count_family, data = transect_model_data)

c_m9_conn_press <- glmmTMB(transect_browser_count ~ rugosity_sc +
                             mpa_status +
                             connectivity_sc * log_market_gravity_sc +
                             (1 | site),
                           family = count_family, data = transect_model_data)

c_sens_settpop <- glmmTMB(transect_browser_count ~ rugosity_sc +
                            log_settlement_pop_sc +
                            mpa_status +
                            connectivity_sc +
                            (1 | site),
                          family = count_family, data = transect_model_data)

c_sens_mktgrav <- glmmTMB(transect_browser_count ~ rugosity_sc +
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

cat("\n--- AICc: Browser count models ---\n")
print(make_aicc_df(model_list_counts))

# ── Browser count model results ───────────────────────────────
#
#                          Model   AICc Delta Weight
#                 MPA x pressure 897.30  0.00  0.40 — BEST
#           Above + connectivity 899.03  1.73  0.17 — equivalent
#                    Above + DHW 899.68  2.38  0.12
#        Connectivity x pressure 899.92  2.62  0.11
#                   Above + chla 900.59  3.29  0.08
#             MPA x connectivity 900.97  3.67  0.06
#       Habitat + pressure + MPA 901.76  4.46  0.04
#             Habitat + pressure 904.15  6.85  0.01
#   Market gravity (sensitivity) 907.24  9.94  0.00
#  Settlement pop. (sensitivity) 912.94 15.64  0.00
#                        Habitat 916.19 18.89  0.00
#
# PRIMARY EVIDENCE THRESHOLD (delta < 2):
#   1. MPA x pressure        (delta 0.00, weight 0.40) — best
#   2. Above + connectivity  (delta 1.73, weight 0.17) — equivalent
#
# NOTABLE: these results are identical to the previous analysis
# using settlement pop. as primary metric — count model ladder
# is completely robust to pressure metric choice. This confirms
# the count signal is genuine and not sensitive to how human
# pressure is measured.
#
# SENSITIVITY METRICS:
#   Market gravity (delta 9.94) and settlement pop. (delta 15.64)
#   not supported — settlement gravity confirmed as appropriate
#   metric for browser counts as for biomass.
#
summary(c_m8_mpa_press)
summary(c_m4_conn)

# ============================================================
#  SYNTHESIS: BROWSER BIOMASS vs COUNT CONCLUSIONS
# ============================================================
#
# Site-level biomass best model:      MPA x pressure (weight 0.74)
# Transect-level biomass best model:  MPA x pressure (weight 0.58)
# Count best models (delta < 2):      MPA x pressure (weight 0.40)
#                                     Above + connectivity (weight 0.17)
#
# ── Rugosity ─────────────────────────────────────────────────
# Strong positive effect on biomass at both levels
# (site β = 0.57***, transect β = 0.62***) but not significant
# for counts (β = 0.10-0.11, p > 0.18). Complex reefs support
# larger-bodied browsers, not more browsers. Contrast with
# total biomass where rugosity was also the dominant driver of
# counts — browsers are less numerically responsive to habitat
# complexity than the total community.
#
# ── MPA effects ──────────────────────────────────────────────
# MPA x pressure dominates browser biomass at all three levels —
# the strongest and most consistent management signal in the
# dataset. Medium protection drives the biomass interaction
# (β ≈ 1.83-1.86**, both levels) while low protection drives
# the count interaction (β = +0.63*, p = 0.017). This suggests:
#   - Low MPAs recover browser numbers
#   - Medium MPAs recover browser body size
# Potentially reflecting sequential stages of community
# recovery under increasing protection intensity.
#
# ── Pressure direction reversal ──────────────────────────────
# Biomass: negative pressure slope at unprotected sites
#          (β = -0.21, n.s.) — no detectable effect, browsers
#          likely already depleted across full gradient
# Counts:  positive pressure slope at unprotected sites
#          (β = +0.29, p = 0.002) — more browsers where
#          pressure is higher
#
# This reversal is biologically meaningful. High-pressure
# unprotected sites have more but smaller-bodied browsers —
# consistent with size-selective fishing removing large
# individuals while small individuals persist or increase
# through release from predation or competitive compensation.
# MPAs reverse this by recovering body size at high-pressure
# sites — medium protection slope for counts (+0.29 + 0.39
# = +0.68, n.s.) vs biomass (+1.61, p = 0.001).
#
# ── Connectivity ─────────────────────────────────────────────
# Significant for counts in both top models:
#   MPA x pressure model: β = +0.26, p = 0.011
#   Above + connectivity: β = +0.23, p = 0.023
# Absent from biomass models at both levels (delta > 5).
# Larval supply supports browser abundance via recruitment
# but does not translate into biomass recovery — connectivity
# delivers small recruits rather than recovering large-bodied
# fish.
#
# ── Pressure positive for counts — interpretation ────────────
# Settlement gravity positively associated with browser counts
# in both top count models (β = +0.29-0.38, p < 0.01).
# More browsers numerically where fishing pressure is higher —
# but these are smaller individuals (biomass negative at
# unprotected sites). Fishing pressure restructures the size
# distribution of browsers without necessarily reducing numbers.
#
# ── Collapsed RE in MPA x pressure count model ───────────────
# Site RE variance ≈ 0 in c_m8_mpa_press — the MPA x pressure
# interaction absorbs all site-level variation. Connectivity
# inference more reliable from c_m4_conn where RE is stable
# (variance = 0.011). Both models agree on connectivity
# direction and magnitude.
#
# ── DHW ──────────────────────────────────────────────────────
# Marginal for counts (delta 2.38) — some evidence thermal
# stress reduces browser numbers, consistent with total fish
# count result. Not supported for biomass at either level.
#
# ── Contrast with total biomass ──────────────────────────────
# Total biomass: habitat + pressure best model, MPA marginal
# Browser biomass: MPA x pressure overwhelmingly dominant,
#                  habitat not independently supported
# Browsers are the functional group most sensitive to the
# protection x pressure interaction — consistent with being
# heavily targeted by fishing and highly responsive to
# exclusion of fishing effort.
#
# ── Pressure metric note ─────────────────────────────────────
# Settlement gravity used as primary pressure metric throughout.
# Settlement pop. tested as primary metric but produced diffuse
# results in full candidate set (6 models within delta < 2,
# MPA x pressure drops to delta 5.29). Settlement gravity
# produces a more interpretable result and is consistent with
# total biomass analyses. Settlement pop. retained as sensitivity
# only. Count model results identical under both metrics —
# confirms count signal is robust to pressure metric choice.
# ============================================================

cat("\n=== BROWSER — Site-level biomass, best model ===\n")
summary(m8_mpa_press)

cat("\n=== BROWSER — Transect-level biomass, best model ===\n")
summary(transect_m8_mpa_press)

cat("\n=== BROWSER — Count models, top two ===\n")
summary(c_m8_mpa_press)
summary(c_m4_conn)

# maybe update models to also test 