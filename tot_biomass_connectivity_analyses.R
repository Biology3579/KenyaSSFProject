# ============================================================
#  DRIVERS OF REEF FISH BIOMASS
#  Chapter 1 — Site-level Analysis
#
#  Analytical framework (four sequential stages):
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
#             Reports both ΔAICc and ΔR² at each step.
#
#  STAGE 3 — Interaction testing (conditional on Stage 2)
#             Three a priori interactions test whether spatial
#             management and connectivity modify local-scale
#             relationships. Run only if Global outperforms
#             Local baseline by ΔAICc > 2.
#
#  STAGE 4 — Sensitivity analysis
#             (a) Alternative pressure metrics — best Stage 2
#                 model refitted with settlement pop. and market
#                 gravity substituted for settlement gravity.
#             (b) Transect-level mixed model — confirms
#                 site-level findings not an artefact of
#                 aggregation.
#
#  Study design:
#    Transects nested within stations, stations within sites,
#    sites within locations, locations within countries.
#
#  Response:
#    log(mean total fish biomass per site)  [site-level LM]
#    log(transect total biomass)            [sensitivity, GLMM]
#
#  Primary pressure metric: settlement gravity (selected via
#    pre-analysis univariate AICc comparison — see Section 2.3)
# ============================================================

options(scipen = 999)

# ── PACKAGES ─────────────────────────────────────────────────
library(tidyverse)
library(sf)
library(glmmTMB)
library(DHARMa)
library(MuMIn)
library(AICcmodavg)
library(vegan)       # varpart()
library(ggcorrplot)
library(corrplot)
library(gridExtra)
library(MASS)
library(here)
library(spdep) # spatial autocorrelation

# ── FUNCTIONS ─────────────────────────────────────────────────

# AICc comparison table from a named list of models
make_aicc_df <- function(model_list) {
  aicc_v  <- sapply(model_list, AICc)
  delta_v <- aicc_v - min(aicc_v)
  wt_v    <- exp(-0.5 * delta_v) / sum(exp(-0.5 * delta_v))
  data.frame(
    Model  = names(model_list),
    AICc   = round(aicc_v,  2),
    Delta  = round(delta_v, 2),
    Weight = round(wt_v,    4),
    row.names = NULL
  ) %>% arrange(AICc)
}


# Marginal effect plot for a single focal predictor
# All other scaled predictors held at 0 (their mean)
plot_effect <- function(model, data, focal_var,
                        x_label,
                        y_label = "Fitted log(biomass)",
                        colour  = "#2c7bb6",
                        n = 200) {
  scaled_vars <- names(data)[endsWith(names(data), "_sc")]
  grid <- as.data.frame(
    matrix(0, nrow = n, ncol = length(scaled_vars),
           dimnames = list(NULL, scaled_vars))
  )
  grid[[focal_var]] <- seq(
    min(data[[focal_var]], na.rm = TRUE),
    max(data[[focal_var]], na.rm = TRUE),
    length.out = n
  )
  
  is_lm <- inherits(model, "lm") && !inherits(model, "glmmTMB")
  pred  <- if (is_lm) {
    predict(model, newdata = grid, se.fit = TRUE)
  } else {
    predict(model, newdata = grid, type = "response",
            se.fit = TRUE, re.form = NA)
  }
  
  grid$fit <- pred$fit
  grid$lwr <- pred$fit - 1.96 * pred$se.fit
  grid$upr <- pred$fit + 1.96 * pred$se.fit
  
  ggplot(grid, aes(x = .data[[focal_var]])) +
    geom_ribbon(aes(ymin = lwr, ymax = upr),
                fill = colour, alpha = 0.2) +
    geom_line(aes(y = fit), colour = colour, linewidth = 1.1) +
    labs(x = x_label, y = y_label) +
    theme_bw(base_size = 13) +
    theme(axis.title = element_text(face = "bold"))
}


# ============================================================
#  DATA LOADING
# ============================================================

fish_data      <- readr::read_rds(here("processed_data", "clean_fish_connectivity.rds"))
gravity_data   <- readr::read_rds(here("city_data", "locations_with_grav_combined.rds"))
chla_data      <- read.csv(here("processed_data", "locations_with_chla_2009.csv"))
rugosity_data  <- readr::read_rds(here("processed_data", "clean_dive_details_connectivity.rds"))
location_data  <- readr::read_rds(here("processed_data", "clean_location_connectivity.rds"))
dhw_data       <- readr::read_rds(here("processed_data", "locations_with_dhw_2009.rds"))


# ============================================================
#  DATA AGGREGATION
#  Minimum 3 transects per site retained.
# ============================================================

total_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_total_biomass = sum(tot_wt_g, na.rm = TRUE),
    transect_total_count   = sum(number,   na.rm = TRUE),
    country  = first(country),
    .groups  = "drop"
  ) %>%
  group_by(site) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(
    site    = as.factor(site),
    country = as.factor(country)
  )

cat("Transects:", nrow(total_transects), "\n")
cat("Sites:",     n_distinct(total_transects$site), "\n")
cat("Countries:", n_distinct(total_transects$country), "\n")


# ============================================================
#  PREDICTOR PREPARATION
# ============================================================

# ── Human gravity metrics ─────────────────────────────────────
gravity_sites <- gravity_data %>%
  st_drop_geometry() %>%
  group_by(site) %>%
  summarise(
    market_gravity  = mean(market_grav,          na.rm = TRUE),
    settlement_pop  = mean(settlement_tot_pop,   na.rm = TRUE),
    settlement_grav = mean(nearest_pop75_grav,   na.rm = TRUE),
    .groups = "drop"
  )

# ── Chlorophyll-a ─────────────────────────────────────────────
chla_sites <- chla_data %>%
  group_by(site) %>%
  summarise(mean_annual_chla = mean(chla_annual_mean, na.rm = TRUE),
            .groups = "drop")

# ── Degree heating weeks ──────────────────────────────────────
dhw_sites <- dhw_data %>%
  filter(!is.na(max_dhw)) %>%
  group_by(site) %>%
  summarise(mean_annual_dhw = mean(max_dhw, na.rm = TRUE),
            .groups = "drop")

# ── Rugosity ──────────────────────────────────────────────────
rugosity_sites <- rugosity_data %>%
  group_by(site) %>%
  summarise(rugosity = mean(rugosity, na.rm = TRUE),
            .groups = "drop")

# ── MPA status and connectivity ───────────────────────────────
location_sites <- location_data %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(
    mpa_status   = first(mpa_status),
    connectivity = mean(prop_connectivity, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
#  TRANSFORMATIONS, SCALING, AND PREDICTOR CHECKS
# ============================================================

#Join all predictors
raw_predictors <- location_sites %>%
  left_join(chla_sites,     by = "site") %>%
  left_join(dhw_sites,      by = "site") %>%
  left_join(rugosity_sites, by = "site") %>%
  left_join(gravity_sites,  by = "site")

# Visualise
#...

# ── Apply transformations ─────────────────────────────────────
# Rugosity:         no transformation (approximately normal)
# Gravity metrics:  log (right-skewed)
# Chla:             log (right-skewed)
# DHW:              log(x + 1) (right-skewed with zeros)
# MPA status:       ordered factor (governance gradient)
# Connectivity:     no transformation (...)

transformed_predictors <- raw_predictors %>%
  transmute(
    site                = site,
    rugosity            = rugosity,
    log_settlement_grav = log(settlement_grav),
    log_settlement_pop  = log(settlement_pop),
    log_market_gravity  = log(market_gravity),
    log_chla            = log(mean_annual_chla),
    log_max_dhw         = log(mean_annual_dhw + 1),
    connectivity        = connectivity,
    mpa_status          = factor(mpa_status,
                                 levels  = c("none", "low", "medium"),
                                 ordered = FALSE)   # ← change here
  )

# ── Standardise continuous predictors ────────────────────────
# z-score scaling (mean = 0, SD = 1) enables direct comparison
# of effect sizes across predictors with different units.

scaled_predictors <- transformed_predictors %>%
  transmute(
    site                   = site,
    rugosity_sc            = as.numeric(scale(rugosity)),
    log_settlement_grav_sc = as.numeric(scale(log_settlement_grav)),
    log_settlement_pop_sc  = as.numeric(scale(log_settlement_pop)),
    log_market_gravity_sc  = as.numeric(scale(log_market_gravity)),
    connectivity_sc        = as.numeric(scale(connectivity)),
    log_chla_sc            = as.numeric(scale(log_chla)),
    log_max_dhw_sc         = as.numeric(scale(log_max_dhw)),
    mpa_status             = mpa_status   # ordered factor — not scaled
  )

# ── Predictor correlation matrix ──────────────────────────────
# Check for blocking collinearity before modelling.
# Rule of thumb: |r| > 0.70 warrants caution; > 0.80 is problematic.

corr_matrix <- scaled_predictors %>%
  mutate(mpa_numeric = as.numeric(mpa_status)) %>%
  dplyr::select(ends_with("_sc"), mpa_numeric) %>%
  rename(
    "Market gravity"     = log_market_gravity_sc,
    "Settlement gravity" = log_settlement_grav_sc,
    "Settlement pop."    = log_settlement_pop_sc,
    "Chlorophyll-a"      = log_chla_sc,
    "Max DHW"            = log_max_dhw_sc,
    "Rugosity"           = rugosity_sc,
    "Connectivity"       = connectivity_sc,
    "MPA status"         = mpa_numeric
  ) %>%
  cor(use = "complete.obs")

# ── Correlation matrix with direction and magnitude ───────────
# Red = positive correlation, blue = negative correlation
# Opacity scales with magnitude — stronger correlations more opaque

corrplot(corr_matrix,           # use signed matrix, not abs()
         method      = "square",
         type        = "lower",
         tl.col      = "black",
         tl.srt      = 0,
         tl.offset   = 0.5,
         addCoef.col = "black",
         number.cex  = 0.8,
         col         = colorRampPalette(
           c("#2c7bb6",  # strong negative — deep blue
             "#abd9e9",  # weak negative — light blue
             "white",    # zero — white
             "#f4a7a3",  # weak positive — light orange
             "#d7191c")  # strong positive — deep red
         )(200),
         is.corr     = TRUE,    # treat as correlation — enables opacity scaling
         mar         = c(0, 0, 4, 2))

# ── Collinearity summary ──────────────────────────────────────
#
# GRAVITY METRICS (settlement gravity / pop / market): r = 0.53–0.54
#   Moderate positive correlations — all three proxy the same
#   underlying construct (local human pressure). Use one per
#   model only. Primary metric selected via pre-analysis
#   AICc comparison (see below).
#
# CHLA vs pressure metrics: r = -0.57 (settlement gravity),
#   -0.23 (settlement pop), -0.33 (market gravity)
#   Negative — productive sites tend to be less fished.
#   Reflects geographic covariation between offshore
#   productivity and remoteness from human settlements.
#   Monitor chla coefficient in models where settlement
#   gravity is included simultaneously.
#
# MPA vs settlement gravity: r = -0.36
# MPA vs settlement pop:     r = -0.43
#   Negative — protected sites tend to have lower nearby
#   human pressure, consistent with MPAs being preferentially
#   placed in more remote locations. MPA and settlement
#   gravity partially capture the same geographic gradient.
#
# DHW vs market gravity: r = 0.42
#   Moderate positive — warmer, higher-DHW sites tend to be
#   near larger markets. Acceptable, monitor in models.
#
# CONNECTIVITY: all |r| < 0.30 across all predictors.
#   Genuinely orthogonal to other predictors in this dataset.
#
# RUGOSITY: all |r| < 0.20 across all predictors.
#   Clean and independent of all other predictors.
#
# No blocking collinearity in primary predictor set.
# All pairwise |r| < 0.60 — threshold for concern.
# Gravity metrics handled by single-metric-per-model rule.


# ============================================================
#  PRE-ANALYSIS: HUMAN PRESSURE METRIC SELECTION
#
#  Three gravity metrics all proxy local human pressure.
#  Select the best-performing metric via univariate AICc
#  before entering the main candidate set. This is a
#  pre-analysis selection step, not part of main inference.
# ============================================================

# ── Build site-level dataset for selection ───────────────────
pressure_selection_data <- total_transects %>%
  left_join(scaled_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    log_mean_biomass       = log(mean(transect_total_biomass, na.rm = TRUE)),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_settlement_pop_sc  = first(log_settlement_pop_sc),
    log_market_gravity_sc  = first(log_market_gravity_sc),
    .groups = "drop"
  )

# ── Univariate models ────────────────────────────────────────
press_settgrav <- lm(log_mean_biomass ~ log_settlement_grav_sc,
                     data = pressure_selection_data)
press_settpop  <- lm(log_mean_biomass ~ log_settlement_pop_sc,
                     data = pressure_selection_data)
press_mktgrav  <- lm(log_mean_biomass ~ log_market_gravity_sc,
                     data = pressure_selection_data)

cat("\n--- Pre-analysis: human pressure metric selection ---\n")
print(make_aicc_df(list(
  "Settlement gravity" = press_settgrav,
  "Settlement pop."    = press_settpop,
  "Market gravity"     = press_mktgrav
)))

# ── Decision ─────────────────────────────────────────────────
# Settlement gravity: AICc = 105.87, weight = 0.67 (SELECTED)
# Settlement pop.:    ΔAICc = 2.41,  weight = 0.20
# Market gravity:     ΔAICc = 3.34,  weight = 0.12
#
# Settlement gravity is the primary pressure metric throughout.
# Settlement pop. and market gravity retained for sensitivity
# analysis only (Stage 4).

# ============================================================
#  ANALYSIS DATASETS
# ============================================================

# ── Combine all predictors ────────────────────────────────────
final_predictors <- scaled_predictors %>%
  dplyr::select(
    site,
    rugosity_sc,                # LOCAL      — habitat complexity
    log_settlement_grav_sc,     # LOCAL      — human pressure (primary)
    connectivity_sc,            # SPATIAL    — larval network position
    mpa_status,                 # GOVERNANCE — protection gradient
    log_chla_sc,                # ENV        — primary productivity
    log_max_dhw_sc,             # ENV        — thermal stress history
    log_settlement_pop_sc,      # SENSITIVITY — alternative pressure
    log_market_gravity_sc       # SENSITIVITY — alternative pressure
  )

# ── Transect-level dataset (sensitivity check) ───────────────
transect_model_data <- total_transects %>%
  left_join(final_predictors, by = "site") %>%
  mutate(log_transect_biomass = log(transect_total_biomass))

cat("\nTransect model data:", nrow(transect_model_data), "rows,",
    n_distinct(transect_model_data$site), "sites\n")
cat("Biomass zeros:", sum(transect_model_data$transect_total_biomass == 0), "\n")

# ── Site-level dataset (primary analysis) ────────────────────
total_model_data <- transect_model_data %>%
  group_by(site, country) %>%
  summarise(
    log_mean_biomass       = log(mean(transect_total_biomass, na.rm = TRUE)),
    mean_biomass           = mean(transect_total_biomass, na.rm = TRUE),
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
    country = as.factor(country)
  )

cat("\nSite model data:", nrow(total_model_data), "sites,",
    n_distinct(total_model_data$country), "countries\n")

# ── Confirm pressure metric selection in multivariate context ─ maybe do this...
# Each model includes the full primary predictor set —
# only the pressure metric differs. This ensures the comparison
# reflects metric performance under the same conditions as the
# main candidate models, not a reduced subset.
# 
# press_multi_settgrav <- lm(log_mean_biomass ~ rugosity_sc +
#                              log_settlement_grav_sc +
#                              connectivity_sc +
#                              mpa_status +
#                              log_chla_sc +
#                              log_max_dhw_sc,
#                            data = total_model_data)
# 
# press_multi_settpop  <- lm(log_mean_biomass ~ rugosity_sc +
#                              log_settlement_pop_sc +
#                              connectivity_sc +
#                              mpa_status +
#                              log_chla_sc +
#                              log_max_dhw_sc,
#                            data = total_model_data)
# 
# press_multi_mktgrav  <- lm(log_mean_biomass ~ rugosity_sc +
#                              log_market_gravity_sc +
#                              connectivity_sc +
#                              mpa_status +
#                              log_chla_sc +
#                              log_max_dhw_sc,
#                            data = total_model_data)
# 
# cat("\n--- Pressure metric selection: multivariate confirmation ---\n")
# print(make_aicc_df(list(
#   "Settlement gravity" = press_multi_settgrav,
#   "Settlement pop."    = press_multi_settpop,
#   "Market gravity"     = press_multi_mktgrav
# )))


# Multivariate confirmation (full predictor set):
#   Settlement gravity: AICc = 104.53, weight = 0.795 (best)
#   Market gravity:     ΔAICc = 3.32,  weight = 0.152
#   Settlement pop.:    ΔAICc = 5.40,  weight = 0.053
#
# Settlement gravity selected as primary pressure metric.
# Selection is consistent and strengthens in multivariate
# context — weight increases from 0.672 to 0.795 when all
# predictors are included simultaneously.
# Market gravity and settlement population retained for
# sensitivity analysis only.

# ── Verify no NAs or zeros in primary predictors and response ─
total_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                connectivity_sc, mpa_status,
                log_chla_sc, log_max_dhw_sc) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

# Check for zeros in response variable
# log(0) = -Inf so zeros in mean_biomass would produce
# -Inf in log_mean_biomass — check both
cat("\nZeros in mean_biomass:", 
    sum(total_model_data$mean_biomass == 0), "\n")
cat("-Inf in log_mean_biomass:", 
    sum(is.infinite(total_model_data$log_mean_biomass)), "\n")
cat("NAs in log_mean_biomass:", 
    sum(is.na(total_model_data$log_mean_biomass)), "\n")

# Summary of response distribution
cat("\nResponse variable summary:\n")
print(summary(total_model_data$log_mean_biomass))


# ============================================================
#  MODEL FAMILY SELECTION (documented, not re-run each time)
#
#  Purpose: select the appropriate error distribution and
#  response transformation before fitting any candidate models.
#
#  Family selection uses the global model (all predictors)
#  to test the distribution under the most demanding conditions
#  the data will face. If the chosen distribution holds with
#  all predictors included, it holds for all reduced models.
#  MPA status included as it is part of the global model —
#  its categorical structure is present in the actual analyses
#  and should therefore be present in the diagnostic test.
#
#  Three candidate families tested on identical predictor set:
#
#  Gaussian (raw):   rejected — heteroscedasticity, non-normality,
#                    curved residual pattern.
#  Gamma (log link): rejected — Q-Q shows systematic deviation
#                    from half-normal line.
#  Gaussian (log):   SELECTED — flat residuals vs fitted,
#                    Q-Q closely follows theoretical line,
#                    homoscedasticity met. Minor upper tail
#                    deviation at 2 sites, within acceptable range.
#
#  AICc (raw scale only — Gaussian log not comparable because
#  response is transformed, changing the likelihood scale):
#    Gamma:         AICc = 1146.93, weight = 1.00
#    Gaussian raw:  AICc = 1191.14, ΔAICc = 44.21
#  Gamma outperforms raw Gaussian, confirming transformation
#  is required. Gaussian log selected on diagnostic quality
#  over Gamma despite Gamma's better raw-scale AICc — the
#  Q-Q deviation in Gamma indicates the error distribution
#  is not well supported by these data regardless of AICc.
#
#  Proceed: lm() on log_mean_biomass throughout.
# ============================================================

# ── Gaussian on raw mean biomass ─────────────────────────────
lm_gaussian_raw <- lm(mean_biomass ~ rugosity_sc +
                        log_settlement_grav_sc +
                        connectivity_sc +
                        mpa_status +
                        log_chla_sc +
                        log_max_dhw_sc,
                      data = total_model_data)

par(mfrow = c(2, 2)); plot(lm_gaussian_raw); par(mfrow = c(1, 1))


# ── Gaussian on log-transformed mean biomass ─────────────────
lm_gaussian_log <- lm(log_mean_biomass ~ rugosity_sc +
                        log_settlement_grav_sc +
                        connectivity_sc +
                        mpa_status +
                        log_chla_sc +
                        log_max_dhw_sc,
                      data = total_model_data)

par(mfrow = c(2, 2)); plot(lm_gaussian_log); par(mfrow = c(1, 1))

# ── Gamma (log link) on raw mean biomass ─────────────────────
glm_gamma <- glm(mean_biomass ~ rugosity_sc +
                   log_settlement_grav_sc +
                   connectivity_sc +
                   mpa_status +
                   log_chla_sc +
                   log_max_dhw_sc,
                 family = Gamma(link = "log"),
                 data   = total_model_data)

par(mfrow = c(2, 2)); plot(glm_gamma); par(mfrow = c(1, 1))

# ── AICc comparison (raw scale only) ─────────────────────────
cat("\n--- Family selection: AICc (raw scale) ---\n")
print(make_aicc_df(list(
  "Gaussian (raw)" = lm_gaussian_raw,
  "Gamma"          = glm_gamma
)))

# ── Clean up ──────────────────────────────────────────────────
# Family selection diagnostic summary:
#
# 1. Gaussian (raw):
#    Residuals vs Fitted: strong curved pattern,
#      heteroscedasticity clearly visible.
#    Q-Q: systematic deviation, heavy upper tail.
#    Scale-Location: strong upward trend — variance
#      increases with fitted values. REJECTED.
#
# 2. Gaussian (log): SELECTED
#    Residuals vs Fitted: flat, no systematic pattern.
#    Q-Q: follows theoretical line closely — minor
#      deviations at sites 27 (lower) and 48 (upper)
#      within acceptable range for n = 54.
#    Scale-Location: broadly flat.
#    No sites outside Cook's distance threshold.
#
# 3. Gamma (log link):
#    Residuals vs Fitted: flatter than raw Gaussian
#      but residual pattern remains.
#    Q-Q: systematic deviation in upper tail —
#      Gamma error distribution not well supported.
#    Scale-Location: slight upward trend remains.
#    Sites 6, 48, 35 approach Cook's distance.
#    REJECTED despite better raw-scale AICc (ΔAICc =
#    46.05 vs Gaussian raw) — diagnostic quality
#    does not support Gamma error structure.
# ============================================================
#  RANDOM EFFECT STRUCTURE (documented, not re-run each time) -- need to run for geomorpholgoy and lcoations
#
#  Country-level RE tested using glmmTMB (Gaussian):
#    No RE:          AICc weight = 0.78
#    (1 | country):  AICc weight = 0.22, ΔAICc = 2.53
#
#  Country-level clustering not supported once environmental
#  and human pressure predictors are included.
#  All site-level models fitted as lm() — no random effects.
# ============================================================


# ============================================================
#  STAGE 1 — VARIANCE PARTITIONING
#
#  Quantifies unique and shared variance attributable to three
#  a priori ecological process groups.
#
#  Groups defined before analysis based on mechanistic role:
#    Local:       rugosity + settlement gravity
#    Spatial:     connectivity
#    Environment: chla + DHW
#
#  MPA status excluded — governance variable, not an ecological
#  filter. Tested independently in Stage 2.
#
#  Uses vegan::varpart() with adjusted R² throughout.
# ============================================================

# ── Build matrix inputs for varpart() ────────────────────────
# varpart() requires predictor matrices, not a formula.
# MPA status is excluded here (see rationale above).

vp_local <- total_model_data %>%
  dplyr::select(rugosity_sc, log_settlement_grav_sc) %>%
  as.data.frame()

vp_spatial <- total_model_data %>%
  dplyr::select(connectivity_sc) %>%
  as.data.frame()

vp_environ <- total_model_data %>%
  dplyr::select(log_chla_sc, log_max_dhw_sc) %>%
  as.data.frame()

y_biomass <- total_model_data$log_mean_biomass

# ── Run variance partition ────────────────────────────────────
vp_result <- varpart(y_biomass,
                     vp_local,
                     vp_spatial,
                     vp_environ)

cat("\n--- Stage 1: Variance partitioning ---\n")
print(vp_result)

# ── Variance plot ──────────────────────────────────────────────
vp_fractions <- data.frame(
  Group      = c("Local", "Spatial", "Environment"),
  Unique     = c(0.242, -0.016, 0.054),
  Unique_plot = c(0.242, 0.001, 0.054),  # floor spatial at 0.001 for visibility
  p_label    = c("adj. R² = 0.242, p = 0.001", 
                 "adj. R² < 0, p = 0.904",
                 "adj. R² = 0.054, p = 0.077"),
  sig        = c("significant", "not significant", "marginal")
) %>%
  mutate(Group = factor(Group, 
                        levels = c("Environment", "Spatial", "Local")))

#chaneg colours
ggplot(vp_fractions, aes(x = Unique_plot, y = Group, 
                         fill = sig)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = p_label),
            hjust = -0.05, size = 3.2) +
  scale_fill_manual(values = c("significant"     = "#0072B2",  # blue
                               "marginal"        = "#009E73",  # green
                               "not significant" = "#bdbdbd")) + # grey 
  scale_x_continuous(limits = c(0, 0.42),
                     name   = "Unique variance explained (adj. R²)") +
  labs(y       = NULL,
       fill    = NULL,
       caption = "n = 54 sites | Permutation tests, 999 permutations | Residual variance = 0.773") +
  theme_bw(base_size = 13) +
  theme(legend.position    = "top",
        axis.title         = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        plot.caption       = element_text(colour = "grey50", size = 9))

# ── Significance tests for individual fractions ───────────────
# Tests the unique fraction of each group (conditioned on others)
# using permutation-based RDA.
cat("\n--- Significance of unique fractions ---\n")
cat("Local unique fraction:\n")
print(anova(rda(y_biomass ~ rugosity_sc + log_settlement_grav_sc +
                  Condition(connectivity_sc) +
                  Condition(log_chla_sc) +
                  Condition(log_max_dhw_sc),
                data = total_model_data)))

cat("Spatial unique fraction:\n")
print(anova(rda(y_biomass ~ connectivity_sc +
                  Condition(rugosity_sc) +
                  Condition(log_settlement_grav_sc) +
                  Condition(log_chla_sc) +
                  Condition(log_max_dhw_sc),
                data = total_model_data)))

cat("Environmental unique fraction:\n")
print(anova(rda(y_biomass ~ log_chla_sc + log_max_dhw_sc +
                  Condition(rugosity_sc) +
                  Condition(log_settlement_grav_sc) +
                  Condition(connectivity_sc),
                data = total_model_data)))

# ── Variance partition results and interpretation ─────────────
#
# Unique local fraction    [a] =  0.242  p = 0.001 ***
# Unique spatial fraction  [b] = -0.016  p = 0.904
# Unique environ fraction  [c] =  0.054  p = 0.077 .
# Total explained (adj R²)     =  0.227
# Residual                     =  0.773
#
# Local processes are the only statistically significant
# contributor to reef fish biomass variance. Rugosity and
# human pressure together explain 24.2% of variance uniquely —
# i.e. independently of connectivity and environmental context.
#
# Connectivity explains no unique variance (adj R² < 0,
# reported as ~0; F = 0.017, p = 0.904). This is not a
# borderline result — the F-statistic indicates a near-complete
# absence of independent spatial signal at this scale.
#
# Environmental context (chla + DHW) explains a modest
# independent fraction (5.4%) that does not reach conventional
# significance (p = 0.077). Reported as a trend only.
#
# All shared fractions [d–g] are near zero or negative,
# indicating that the three process groups are largely
# orthogonal in this dataset. The unique local fraction is
# therefore a reliable estimate, not a conservative lower
# bound inflated by collinearity.
#
# Negative adjusted R² values arise when the penalty for
# degrees of freedom exceeds the raw explained variance.
# They indicate the fraction is indistinguishable from zero
# and are reported as ~0 throughout.
#
# CAUTION: The absence of a spatial signal should be
# interpreted as scale-dependent, not as evidence that
# connectivity is ecologically irrelevant. At the spatial
# grain of this study, biomass is more strongly constrained
# by local carrying capacity and exploitation intensity than
# by larval network position.

# ============================================================
#  STAGE 2 — HIERARCHICAL MODEL COMPARISON
#
#  Nested sequence adds process groups progressively.
#  Tests incremental explanatory value of each group beyond
#  the local ecological baseline.
#
#  Sequence:
#    Null            → intercept only
#    Local           → habitat + pressure
#    Local + env     → adds environmental context
#    Local + spatial → adds larval connectivity
#    Local + MPA     → adds governance gradient (independent test)
#    Global          → all predictors
#
#  Both ΔAICc and ΔR² reported at each step.
# ============================================================

# ── Null model ────────────────────────────────────────────────
m_null <- lm(log_mean_biomass ~ 1,
             data = total_model_data)

# ── Local ecological baseline (Tier 1) ───────────────────────
# Habitat complexity + human pressure only.
# Represents the hypothesis that local processes are sufficient
# to explain biomass without spatial or environmental context.
m_local <- lm(log_mean_biomass ~ rugosity_sc +
                log_settlement_grav_sc,
              data = total_model_data)

# ── Local + environmental context ────────────────────────────
# Adds productivity and thermal stress to the local baseline.
# Tests whether background abiotic conditions explain additional
# variance beyond local habitat and pressure.
m_local_env <- lm(log_mean_biomass ~ rugosity_sc +
                    log_settlement_grav_sc +
                    log_chla_sc +
                    log_max_dhw_sc,
                  data = total_model_data)

# And separate
m_local_chla <- lm(log_mean_biomass ~ rugosity_sc +
                     log_settlement_grav_sc +
                     log_chla_sc,
                   data = total_model_data)

m_local_dhw  <- lm(log_mean_biomass ~ rugosity_sc +
                     log_settlement_grav_sc +
                     log_max_dhw_sc,
                   data = total_model_data)

# ── Local + spatial processes ────────────────────────────────
# Adds larval connectivity to the local baseline.
# Tests whether network position explains additional variance
# beyond local habitat and pressure.
m_local_spatial <- lm(log_mean_biomass ~ rugosity_sc +
                        log_settlement_grav_sc +
                        connectivity_sc,
                      data = total_model_data)

# ── Local + MPA status (governance test) ─────────────────────
# Adds the protection gradient to the local baseline.
# Key test: does management regime explain variance beyond
# local ecological processes alone?
# MPA treated as ordered factor: none < low < medium.
m_local_mpa <- lm(log_mean_biomass ~ rugosity_sc +
                    log_settlement_grav_sc +
                    mpa_status,
                  data = total_model_data)

# ── Global model (all processes) ─────────────────────────────
# Full additive model including all hypothesised drivers.
# Upper bound on explained variance for this predictor set.
m_global <- lm(log_mean_biomass ~ rugosity_sc +
                 log_settlement_grav_sc +
                 connectivity_sc +
                 mpa_status +
                 log_chla_sc +
                 log_max_dhw_sc,
               data = total_model_data)

# ── Hierarchical comparison table ────────────────────────────
# Listed in intended hierarchical order — ΔR² is the increment
# relative to the previous model in the sequence.
model_list <- list(
  "Null"            = m_null,
  "Local"           = m_local,
  "Local + chla"    = m_local_chla,
  "Local + DHW"     = m_local_dhw,
  "Local + env"     = m_local_env,
  "Local + MPA"     = m_local_mpa,
  "Local + spatial" = m_local_spatial,
  "Global"          = m_global
)

# AICc-ranked table — for model selection
cat("\n--- Stage 2: Model comparison (AICc ranked) ---\n")
print(make_aicc_df(model_list))

# Hierarchical table — for R² increments, in sequence order
cat("\n--- Stage 2: Variance explained (hierarchical sequence) ---\n")
local_r2 <- summary(m_local)$adj.r.squared

model_list %>%
  imap_dfr(~ tibble(
    Model  = .y,
    Adj_R2 = round(summary(.x)$adj.r.squared, 3)
  )) %>%
  mutate(Delta_R2 = round(Adj_R2 - local_r2, 3),
         Delta_R2 = ifelse(Model %in% c("Null", "Local"), NA, Delta_R2)) %>%
  arrange(Model == "Null",        # Null always last
          desc(Model == "Local"), # Local always first
          desc(Adj_R2)) %>%       # everything else by R²
  print()

# ── Key contrasts ─────────────────────────────────────────────
# Local vs Null:           Local baseline explains 18.4% of variance
#                          (ΔAICc = 8.49) — local processes clearly
#                          supported over intercept-only model.
#
# Local + env vs Local:    Adding environmental context (chla + DHW)
#                          improves fit by ΔR² = 0.059 — best supported
#                          model (ΔAICc = 0.00, weight = 0.30).
#
# Local + DHW vs Local:    DHW alone adds ΔR² = 0.038 — thermal stress
#                          is the stronger environmental driver.
#
# Local + chla vs Local:   Chla alone adds ΔR² = 0.010 — productivity
#                          contributes modestly and independently.
#
# Local + MPA vs Local:    MPA status adds ΔR² = 0.008 (ΔAICc = 3.52)
#                          — governance gradient not supported beyond
#                          local ecological baseline.
#
# Local + spatial vs Local: Connectivity adds ΔR² = -0.011 (ΔAICc = 3.30)
#                          — spatial processes add no explanatory value.
#
# Global vs Local:         Full model adds ΔR² = 0.056 over Local but
#                          ΔAICc = 5.11 — gains driven by environmental
#                          terms only; MPA and connectivity contribute
#                          dead-weight parameters. Global is the worst
#                          supported model after Null.
#
# Best-supported model: Local + env (AICc = 99.43, weight = 0.30)
# Model selection uncertainty: Local + DHW competitive (ΔAICc = 0.02)
# Combined weight of Local + env and Local + DHW = 0.59
# ΔAICc > 2 threshold: only Local, Local + env, Local + DHW supported

# ── Coefficient summary: best-supported model (Local + env) ───
summary(m_local_env)   # already done
summary(m_local_dhw)   
summary(m_local)

# ── Results: best-supported model (Local + env) ───────────────
# lm(log_mean_biomass ~ rugosity_sc + log_settlement_grav_sc +
#    log_chla_sc + log_max_dhw_sc)
# n = 54 sites, 4 countries
#
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)
# (Intercept)             9.6539     0.0765  126.23  < 0.001 ***
# rugosity_sc             0.2095     0.0773    2.71   0.0092 **
# log_settlement_grav_sc -0.2794     0.0947   -2.95   0.0049 **
# log_chla_sc            -0.1476     0.0958   -1.54   0.1297
# log_max_dhw_sc          0.1611     0.0783    2.06   0.0451 *
#
# Residual standard error: 0.561 on 49 df
# Adjusted R² = 0.243
# F(4, 49) = 5.25, p = 0.001
#
# Coefficient stability across competitive models (ΔAICc < 2):
#
#                      Local      Local+DHW   Local+env
# rugosity_sc         +0.221**   +0.217**    +0.209**
# settlement_grav_sc  -0.184*    -0.198*     -0.279**
# log_max_dhw_sc          —      +0.147.     +0.161*
# log_chla_sc             —          —       -0.148 ns
#
# Rugosity and settlement gravity are stable in direction and
# magnitude across all three competitive models, confirming
# that the local signal is robust and does not depend on which
# environmental terms are included.
#
# Interpretation:
#
# Rugosity (β = +0.209, p = 0.009):
#   Positive effect — greater structural complexity associated
#   with higher biomass, consistent with rugosity setting local
#   carrying capacity through provision of refuge and foraging
#   habitat. Effect is significant and stable across all
#   competitive models.
#
# Settlement gravity (β = -0.279, p = 0.005):
#   Negative effect — higher human pressure associated with
#   lower biomass, consistent with exploitation reducing
#   standing stock. Effect strengthens slightly when
#   environmental terms are included, suggesting modest
#   suppression by collinearity with chla (r = 0.57).
#   Significant and stable across all competitive models.
#
# Chla (β = -0.148, p = 0.130):
#   Non-significant. Negative direction is counterintuitive
#   and likely reflects moderate collinearity with settlement
#   gravity (r = 0.57), which destabilises the chla coefficient
#   when both are included. Interpret with caution.
#
# DHW (β = +0.161, p = 0.045):
#   Positive and significant — counterintuitive given expected
#   negative effect of thermal stress on coral and fish biomass.
#   Not driven by zeros (15 zero-DHW sites have near-identical
#   mean biomass to non-zero sites: 9.62 vs 9.67) and not
#   confounded with fishing pressure (r = 0.10 with settlement
#   gravity). Most likely reflects geographic covariation —
#   higher DHW sites tend to be in warmer, more tropical
#   environments that support greater reef fish biomass
#   independently of bleaching history. Interpret as a
#   geographic context variable rather than a direct
#   thermal stress effect. Address in discussion.

# ============================================================
#  SPATIAL AUTOCORRELATION CHECK
#  Tests residuals from best-supported Stage 2 model.
#  Run once — result reported as a limitation in methods.
# ============================================================

# ── Site coordinates ──────────────────────────────────────────
site_coords <- location_data %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(lon = first(longitude),
            lat = first(latitude),
            .groups = "drop")

total_model_data_coords <- total_model_data %>%
  left_join(site_coords, by = "site")

# ── Spatial weights (k = 5 nearest neighbours) ────────────────
coords_mat <- cbind(total_model_data_coords$lon,
                    total_model_data_coords$lat)
listw5 <- nb2listw(knn2nb(knearneigh(coords_mat, k = 5)),
                   style = "W")

# ── Moran's I on best-supported model residuals ───────────────
cat("\n--- Moran's I: residuals from Local + env (best Stage 2) ---\n")
print(moran.test(residuals(m_local_env), listw5))

# ── Result ────────────────────────────────────────────────────
# I = 0.108, p = 0.043 — weak but statistically significant
# positive spatial autocorrelation in residuals.
#
# Warnings noted:
#   4 sub-graphs: reflects 4 geographically separate
#     country clusters. k-NN weights bridge clusters
#     correctly — acceptable for this dataset.
#
# Effect is modest and decreases from Local alone
# (I = 0.142, p = 0.015), confirming predictors absorb
# some spatial structure. 
# Acknowledged as a limitation — # n = 54 precludes reliable spatial error modelling.

# ============================================================
#  STAGE 3 — INTERACTION TESTING
#
#  Tests whether spatial management (MPA) and connectivity
#  modify local-scale relationships — the mechanistic
#  hypothesis underlying the Tier 2 structure.
#
#  Run only if Global outperforms Local by ΔAICc > 2.
#  If the Tier 2 variables add little, interactions are
#  untestable and this stage is reported as not supported.
#
#  Three a priori interactions, each representing a distinct
#  mechanistic hypothesis:
#    MPA × connectivity:  effectiveness depends on larval supply
#    MPA × pressure:      effectiveness depends on fishing intensity
#    Connectivity × pressure: connectivity buffers fishing impact
#
#  Compared against the additive Global model as reference.
#  Interaction preferred only if ΔAICc > 2 vs Global additive.
# ============================================================

# ── Check Stage 2 before proceeding ──────────────────────────
delta_global_vs_local <- AICc(m_local) - AICc(m_global)
cat("\nΔAICc (Local vs Global):", round(delta_global_vs_local, 2), "\n")

# ── Gate check result ─────────────────────────────────────────
# ΔAICc (Local vs Global) = -3.88
# Global is WORSE than Local by 3.88 AICc units.
# MPA and connectivity add no explanatory value as main effects.
# Interaction models are fitted for completeness and to satisfy
# the a priori analytical framework, but any significant
# interaction should be interpreted with considerable caution —
# it would be difficult to justify a context-dependent effect
# of variables whose main effects are not supported.
# Results reported in supplementary material if no interaction
# is supported; in main text only if ΔAICc > 2 vs Local + env.

# ── Hypothesis 1: MPA effectiveness depends on larval supply ──
# Protected sites with high connectivity should recover faster
# and maintain higher biomass than isolated protected sites
m_int_mpa_conn <- lm(log_mean_biomass ~ rugosity_sc +
                       log_settlement_grav_sc +
                       log_chla_sc +
                       log_max_dhw_sc +
                       mpa_status * connectivity_sc,
                     data = total_model_data)

# ── Hypothesis 2: MPA effectiveness depends on fishing intensity
# MPAs should only be detectable where external pressure is low
# enough that protection translates into a biomass difference
m_int_mpa_press <- lm(log_mean_biomass ~ rugosity_sc +
                        log_chla_sc +
                        log_max_dhw_sc +
                        mpa_status * log_settlement_grav_sc,
                      data = total_model_data)

# ── Hypothesis 3: connectivity buffers fishing pressure ───────
# Well-connected sites should be more resilient to exploitation
# through sustained larval replenishment offsetting mortality
m_int_conn_press <- lm(log_mean_biomass ~ rugosity_sc +
                         log_chla_sc +
                         log_max_dhw_sc +
                         connectivity_sc * log_settlement_grav_sc,
                       data = total_model_data)

# ── Interaction candidate set ─────────────────────────────────
model_list_interactions <- list(
  "Local + env (additive)"  = m_local_env,
  "MPA × connectivity"      = m_int_mpa_conn,
  "MPA × pressure"          = m_int_mpa_press,
  "Connectivity × pressure" = m_int_conn_press
)

cat("\n--- Stage 3: Interaction model comparison ---\n")
print(make_aicc_df(model_list_interactions))


# ── Interaction candidate set ─────────────────────────────────
# Global additive model is the reference — we ask whether any
# interaction improves on the additive structure.
model_list_interactions <- list(
  "Global (additive)"       = m_global,
  "MPA × connectivity"      = m_int_mpa_conn,
  "MPA × pressure"          = m_int_mpa_press,
  "Connectivity × pressure" = m_int_conn_press
)

cat("\n--- Stage 3: Interaction model comparison ---\n")
print(make_aicc_df(model_list_interactions))

# ── Interpretation ────────────────────────────────────────────
# An interaction model is preferred only if ΔAICc > 2 vs Global
# additive. If no interaction is supported, report this as
# evidence that Tier 2 effects are additive rather than
# context-dependent at this spatial scale.

summary(m_int_mpa_press)
summary(m_int_conn_press)

# ── Stage 3 conclusion ────────────────────────────────────────
# No interaction model provides a convincing, ecologically
# interpretable improvement over Local + env.
#
# MPA × pressure (ΔAICc = 1.29) is statistically competitive
# but the interaction pattern is not consistent with the a
# priori hypothesis. The main effect of settlement gravity
# collapses entirely within the interaction model, suggesting
# the signal reflects geographic covariation rather than a
# genuine management-by-pressure effect.
#
# Connectivity × pressure (ΔAICc = 1.81, p = 0.079) is
# marginal and not supported at conventional thresholds.
#
# MPA × connectivity (ΔAICc = 5.13) is not supported.
#
# Primary inference is unchanged: reef fish biomass is
# structured by local ecological processes. Spatial management
# and connectivity do not modify this relationship in a
# statistically or ecologically robust way at this scale.

# ============================================================
#  STAGE 4 — SENSITIVITY ANALYSIS
#
#  (a) Alternative pressure metrics
#      Best Stage 2 model refitted substituting settlement pop.
#      and market gravity for settlement gravity. All other
#      predictors retained to ensure comparability.
#      Results reported in Supplementary Table S2.
#
#  (b) Transect-level replication
#      Site-level hierarchical sequence repeated at transect
#      level using GLMM with (1 | site). Confirms site-level
#      findings are not an artefact of averaging.
#      Results reported in Supplementary Table S1.
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
# Mirrors the Global model structure — only pressure metric swapped.
# NOTE: earlier sensitivity models dropped env. variables,
# making them non-comparable. These are corrected here.

sens_settpop <- lm(log_mean_biomass ~ rugosity_sc +
                     log_settlement_pop_sc +      
                     connectivity_sc +
                     mpa_status +
                     log_chla_sc +
                     log_max_dhw_sc,
                   data = total_model_data)

sens_mktgrav <- lm(log_mean_biomass ~ rugosity_sc +
                     log_market_gravity_sc +      
                     connectivity_sc +
                     mpa_status +
                     log_chla_sc +
                     log_max_dhw_sc,
                   data = total_model_data)

cat("\n--- Stage 4a: Sensitivity — alternative pressure metrics ---\n")
cat("(Coefficients only — these are not compared via AICc\n",
    "because they test robustness, not process importance)\n\n")
cat("Settlement population:\n"); print(summary(sens_settpop)$coefficients)
cat("\nMarket gravity:\n");      print(summary(sens_mktgrav)$coefficients)

# ── Stage 4a results ──────────────────────────────────────────
# Robustness check: primary conclusions hold across alternative
# pressure metrics.
#
# Rugosity coefficient stability:
#   Primary (settlement gravity): β = +0.209, p = 0.009
#   Settlement population:        β = +0.210, p = 0.014
#   Market gravity:               β = +0.192, p = 0.024
#
# Rugosity is significant and stable in direction and magnitude
# across all three pressure metrics. The habitat complexity
# signal is robust to pressure proxy choice.
#
# Pressure metric performance:
#   Settlement population: β = -0.110, p = 0.258 — not significant
#   Market gravity:        β = -0.175, p = 0.081 — marginal
#
# Neither alternative reaches significance, confirming
# settlement gravity as the strongest pressure proxy.
#
# All other predictors (connectivity, MPA, chla) non-significant
# across both models — consistent with Stage 2.
#
# DHW positive in both models (β = +0.144 and +0.209) —
# consistent with primary model interpretation.
#
# Conclusion: sensitivity analysis supports primary inference.
# Rugosity is the most robust predictor regardless of pressure
# metric. Human pressure direction is consistent but weakens
# with alternative metrics, further justifying settlement
# gravity selection.

# ── (b) Transect-level sensitivity ───────────────────────────
# Replicates the Stage 2 hierarchical sequence at transect level.
# (1 | site) accounts for non-independence of transects within sites.

sens_t_null <- glmmTMB(log_transect_biomass ~ 1 +
                         (1 | site),
                       family = gaussian(),
                       data   = transect_model_data)

sens_t_local <- glmmTMB(log_transect_biomass ~ rugosity_sc +
                          log_settlement_grav_sc +
                          (1 | site),
                        family = gaussian(),
                        data   = transect_model_data)

sens_t_local_env <- glmmTMB(log_transect_biomass ~ rugosity_sc +
                              log_settlement_grav_sc +
                              log_chla_sc +
                              log_max_dhw_sc +
                              (1 | site),
                            family = gaussian(),
                            data   = transect_model_data)

sens_t_local_spatial <- glmmTMB(log_transect_biomass ~ rugosity_sc +
                                  log_settlement_grav_sc +
                                  connectivity_sc +
                                  (1 | site),
                                family = gaussian(),
                                data   = transect_model_data)

sens_t_local_mpa <- glmmTMB(log_transect_biomass ~ rugosity_sc +
                              log_settlement_grav_sc +
                              mpa_status +
                              (1 | site),
                            family = gaussian(),
                            data   = transect_model_data)

sens_t_global <- glmmTMB(log_transect_biomass ~ rugosity_sc +
                           log_settlement_grav_sc +
                           connectivity_sc +
                           mpa_status +
                           log_chla_sc +
                           log_max_dhw_sc +
                           (1 | site),
                         family = gaussian(),
                         data   = transect_model_data)

model_list_transect <- list(
  "Null"            = sens_t_null,
  "Local"           = sens_t_local,
  "Local + env"     = sens_t_local_env,
  "Local + spatial" = sens_t_local_spatial,
  "Local + MPA"     = sens_t_local_mpa,
  "Global"          = sens_t_global
)

cat("\n--- Stage 4b: Sensitivity — transect-level ---\n")
print(make_aicc_df(model_list_transect))

# ── Stage 4b results ──────────────────────────────────────────
# Transect-level mixed model (1 | site) replicates site-level
# findings — conclusions are not an artefact of aggregation.
#
# Model ranking (transect-level):
#   Local + env:     AICc = 634.92, weight = 0.663 — best supported
#   Local:           ΔAICc = 2.83,  weight = 0.161 — competitive
#   Local + spatial: ΔAICc = 4.63  — not supported
#   Global:          ΔAICc = 4.82  — not supported
#   Local + MPA:     ΔAICc = 5.16  — not supported
#   Null:            ΔAICc = 13.06 — worst supported
#
# Comparison with site-level ranking:
#   Site-level:      Local + env best (weight = 0.297)
#   Transect-level:  Local + env best (weight = 0.663)
#
# Model ordering is identical across both levels of analysis.
# Local + env is best supported at both levels. Local + spatial,
# Local + MPA, and Global are all unsupported at both levels.
# The stronger weight at transect level reflects greater
# statistical power from 243 observations vs 54 sites.
#
# Conclusion: site-level aggregation does not alter qualitative
# inference. Local ecological processes (habitat complexity and
# human pressure) with environmental context (chla + DHW)
# provide the best explanation of reef fish biomass at both
# the transect and site level. Spatial and governance variables
# add no explanatory value at either level of analysis.

# ============================================================
#  MARGINAL EFFECT PLOTS
#  Generated for primary predictors from best supported model:
#  Local + env (AICc = 99.43, weight = 0.297)
#  All other standardised predictors held at their mean (= 0).
# ============================================================

best_model <- m_local_env

# ── Continuous predictors ─────────────────────────────────────
p_rugosity <- plot_effect(best_model, total_model_data,
                          "rugosity_sc",
                          "Rugosity (standardised)")

p_pressure <- plot_effect(best_model, total_model_data,
                          "log_settlement_grav_sc",
                          "log(Settlement gravity) (standardised)")

p_chla <- plot_effect(best_model, total_model_data,
                      "log_chla_sc",
                      "log(Chlorophyll-a) (standardised)")

p_dhw <- plot_effect(best_model, total_model_data,
                     "log_max_dhw_sc",
                     "log(Max DHW + 1) (standardised)")

# Connectivity excluded from main marginal effect plots —
# not included in best supported model and adds no explanatory
# value (Stage 2: ΔAICc = 3.30, Stage 1: unique fraction ≈ 0)

# jpeg("marginal_effects_main.jpg",
#      width = 33, height = 11, units = "cm", res = 300)
# gridExtra::grid.arrange(p_rugosity, p_pressure,
#                         p_chla,    p_dhw,
#                         ncol = 4)
# dev.off()

# ── MPA marginal means ────────────────────────────────────────
# MPA status not included in best supported model (Local + env).
# Marginal means generated from m_local_mpa for completeness
# and to report the direction of the protection gradient.
# Interpret with caution — MPA not supported in Stage 2
# (ΔAICc = 3.52 vs Local).
# Unordered factor — coefficients are dummy contrasts:
#   low vs none, medium vs none.

mpa_grid <- data.frame(
  mpa_status             = factor(c("none", "low", "medium"),
                                  levels  = c("none", "low", "medium"),
                                  ordered = FALSE),
  rugosity_sc            = 0,
  log_settlement_grav_sc = 0,
  log_chla_sc            = 0,
  log_max_dhw_sc         = 0
)

mpa_pred <- predict(m_local_mpa, newdata = mpa_grid, se.fit = TRUE)
mpa_grid$fit <- mpa_pred$fit
mpa_grid$lwr <- mpa_pred$fit - 1.96 * mpa_pred$se.fit
mpa_grid$upr <- mpa_pred$fit + 1.96 * mpa_pred$se.fit

p_mpa <- ggplot(mpa_grid, aes(x = mpa_status, y = fit)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  colour = "#2c7bb6", size = 0.8) +
  labs(x = "MPA status",
       y = "Fitted log(biomass)",
       caption = "From m_local_mpa — MPA not supported in Stage 2") +
  theme_bw(base_size = 13) +
  theme(axis.title   = element_text(face = "bold"),
        plot.caption = element_text(colour = "grey50", size = 9))

# jpeg("marginal_effect_mpa.jpg",
#      width = 12, height = 12, units = "cm", res = 300)
# print(p_mpa)
# dev.off()

# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()