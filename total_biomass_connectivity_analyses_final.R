# Driver of reef fish biomass - Total Biomass Analyses ---
# Created by: Candela Ferrer Diez
# Date: 24/01/2026
# ──────────────────────────────────────────────────────────────────────────────
# 
# Introduction ----
# Analytical framework mirrors total biomass (Q1–Q3 + sensitivity analyses),
# with group-specific family selection and pressure metric confirmation.
#
# Scientific questions:
#  Q1 — Pressure metric validity
#       Does settlement gravity outperform market gravity and
#       settlement population as a proxy for small-scale
#       fisheries pressure on reef fish biomass?
#       Tested via multivariate AICc comparison — three models
#       identical in structure, differing only in pressure
#       metric.
#
# Q2 — Larval connectivity
#       Does larval connectivity explain additional variation
#       in total biomass beyond the human pressure
#       baseline, and does it modify the relationship
#       between human pressure and total biomass?
#
# Q3 — Formal protection
#       Does MPA status explain additional variation in
#       total biomass beyond the fully specified pressure
#       and connectivity model?
#
# Rationale for sequence:
#       Pressure first as primary driver of interest, connectivity second
#       to extend Warmuth et al. (2024), MPA last as
#       governance response downstream of both.
#
# Baseline model (fixed a priori, never tested):
#       biomass ~ rugosity_sc + log_chla_sc
#       Identical justification to total biomass.
#       Pressure metric selected in Q1.
#
# Key difference from total biomass:
#       Browser biomass has ~11% zeros at site level.
#       Tweedie distribution selected over Gaussian log
#       (see family selection). All models use
#       glmmTMB(family = tweedie(link = "log")) on
#       raw mean_biomass throughout.
#
# Sensitivity analyses:
#       (a) Alternative pressure metrics
#       (b) Transect-level replication GLMM — confirms site-level findings
#           not an artefact of aggregation
# ──────────────────────────────────────────────────────────────────────────────

# Set up ----
options(scipen = 999) # avoids scientific notation (for ease of interpretation)

## Packages ----- 
library(tidyverse)
library(sf)
library(MuMIn)
library(gridExtra)
library(here)
library(spdep)
library(ggplot2)
library(lme4)
library(car)

## Functions -----

# AICc comparison table
make_aicc_df <- function(model_list) {
  aicc_v <- sapply(model_list, AICc)
  delta_v <- aicc_v - min(aicc_v)
  wt_v <- exp(-0.5 * delta_v) / sum(exp(-0.5 * delta_v))
  data.frame( Model = names(model_list),
              AICc = round(aicc_v, 2),
              Delta = round(delta_v, 2),
              Weight = round(wt_v, 4),
              row.names = NULL) %>% 
  arrange(AICc)
}

# Marginal effect plots 
plot_effect <- function(model, data, focal_var, x_label, y_label = "Fitted log(biomass)",
                        colour  = "#2c7bb6", n = 200) {
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
  is_lm <- inherits(model, "lm")
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
                fill = colour, alpha = 0.15) +
    geom_line(aes(y = fit), colour = colour, linewidth = 1.1) +
    geom_point(data = data,
               aes(x = .data[[focal_var]],
                   y = log_mean_biomass),
               colour = "grey40", size = 1.5,
               alpha  = 0.5, inherit.aes = FALSE) +
    labs(x = x_label, y = y_label) +
    theme_bw(base_size = 13) +
    theme(axis.title      = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

## Load data ----
fish_data <- readr::read_rds(here("processed_data", "clean_fish_connectivity.rds"))
gravity_data <- readr::read_rds(here("city_data", "locations_with_grav_combined.rds"))
rugosity_data <- readr::read_rds(here("processed_data", "clean_dive_details_connectivity.rds"))
location_data <- readr::read_rds(here("processed_data", "clean_location_connectivity.rds"))
chla_data <- read.csv(here("processed_data", "locations_with_chla_2009.csv"))
# ──────────────────────────────────────────────────────────────────────────────

# Exploring fish data ----
# 
## Aggregate data at transect level ----
# Set a minimum of 3 transects per site.
total_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(transect_total_biomass = sum(tot_wt_g, na.rm = TRUE),
            transect_total_count= sum(number, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(site) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate( site= as.factor(site))

cat("Transects:", nrow(total_transects), "\n")
cat("Sites:", n_distinct(total_transects$site), "\n")

## Aggregate data at site level (mean ts_biomass per site) ----
site_data <- total_transects %>%
  group_by(site) %>%
  summarise(mean_biomass = mean(transect_total_biomass, na.rm = TRUE),
            n_transects  = n(),
            .groups = "drop") %>%
  mutate(site = as.factor(site))

## Basic summary ----
summary(site_data$mean_biomass)

zeros <- mean(site_data$mean_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")

## Raw distribution ----
( site_raw <- ggplot(site_data, aes(x = mean_biomass)) +
    geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
    labs(x = "Mean total biomass per site (g)", y = "Frequency",
         title = "Raw Site-Level Biomass") +
    theme_bw() )
# Great right-skewdness

## Normality checks ----
## Normality checks on log-transformed response ----
qqnorm(site_data$mean_biomass,
       main = "Q-Q plot: Mean biomass per site)")
qqline(site_data$mean_biomass, col = "red")
shapiro.test(site_data$mean_biomass)
# Raw biomass
cat("\n--- Shapiro-Wilk: raw biomass ---\n")
shapiro.test(site_data$mean_biomass)

#Raw biomass was strongly right-skewed (Shapiro-Wilk: W = 0.764, p < 0.001) and 
#was log-transformed prior to analysis."

## Box-Cox ----
## What power transformation does the data suggest?
MASS::boxcox(
  lm(mean_biomass ~ 1, data = site_data),
  lambda = seq(-2, 2, 0.1)
)
# Lambda has a wide range but includes 0, therefore log transformation appropriate
# Square-root also tried.

## Apply transformation and visualise ----
site_data <- site_data %>%
  mutate(log_mean_biomass  = log(mean_biomass))

( site_log <- ggplot(site_data, aes(x = log_mean_biomass)) +
    geom_histogram(bins = 25, fill = "#1a9641", colour = "white") +
    labs(x = "log(mean biomass)", y = "Frequency",
         title = "Log-transformed Site-Level Biomass") +
    theme_bw() )
# Right skew mostly corrected

## Normality checks on log-transformed response ----
qqnorm(site_data$log_mean_biomass,
       main = "Q-Q plot: log(mean biomass per site)")
qqline(site_data$log_mean_biomass, col = "red")
cat("\n--- Shapiro-Wilk: log(biomass) ---\n")
shapiro.test(site_data$log_mean_biomass)

# Confirms that there is no strong evidence that your log_mean_biomass deviates 
# from normality, therefore follows the assumptions of the lm. 

## Check variation by site ----
ggplot(site_data, aes(x = reorder(site, mean_biomass, median),
                      y = mean_biomass)) +
  geom_col(fill = "#2c7bb6", alpha = 0.7) +
  coord_flip() +
  labs(x = NULL, y = "Mean total biomass (g)",
       title = "Mean biomass by site (raw)") +
  theme_bw(base_size = 9)
# Great variability across sites
# ──────────────────────────────────────────────────────────────────────────────

# Preparing predictors ----
#
## Site-level preparatiom ----

###  Human gravity metrics ----
# Three metrics proxy small-scale fisheries pressure
# All retained here — primary metric selected in Q1
gravity_sites <- gravity_data %>%
  st_drop_geometry() %>%
  group_by(site) %>%
  summarise(market_gravity = mean(market_grav, na.rm = TRUE),
            settlement_pop  = mean(settlement_tot_pop, na.rm = TRUE),
            settlement_grav = mean(nearest_pop75_grav, na.rm = TRUE),
            .groups = "drop")

### Chlorophyll-a ----
# Baseline covariate: background productivity
chla_sites <- chla_data %>%
  group_by(site) %>%
  summarise(mean_annual_chla = mean(chla_annual_mean, na.rm = TRUE),
            .groups = "drop")

### Rugosity ----
# Baseline covariate: habitat structural complexity.
rugosity_sites <- rugosity_data %>%
  group_by(site) %>%
  summarise(rugosity = mean(rugosity, na.rm = TRUE),
            .groups = "drop")

### Connectivity and MPA status ----
# Connectivity: candidate driver tested in Q2 (main effect)
#   and Q3 (interaction with pressure and MPA).
# MPA: governance modifier — tested in Q3,
location_sites <- location_data %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(
    mpa_status = first(mpa_status),
    connectivity = mean(prop_connectivity, na.rm = TRUE),
    ecoregion = first(ecoregion),
    .groups = "drop"
  )

## Join predictors ----
raw_predictors <- location_sites %>%
  left_join(rugosity_sites, by = "site") %>%
  left_join(chla_sites, by = "site") %>%
  left_join(gravity_sites, by = "site")

## Check for transformations ----
predictor_labels <- c("Rugosity", "Chlorophyll-a", "Settlement gravity", "Market gravity",
                      "Settlement pop.", "Connectivity")

predictor_order_raw  <- c("rugosity", "mean_annual_chla", "settlement_grav", "market_gravity", 
                           "settlement_pop", "connectivity")
# Plot raw predictor histograms
( p_pred_raw <- raw_predictors %>%
    dplyr::select(all_of(predictor_order_raw)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, predictor_order_raw, predictor_labels)) %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
    facet_wrap(~ variable, scales = "free") +
    labs(title = "Raw predictors") + theme_bw() )

## Transformations ----
# Rugosity: no transformation (approximately normal)
# Gravity metrics: log (right-skewed)
# Chla: log (right-skewed)
# Connectivity: no transformation
# MPA status: unordered factor (none / low / medium)

transformed_predictors <- raw_predictors %>%
  transmute(
    site = site,
    ecoregion = ecoregion,
    rugosity = rugosity,
    log_chla = log(mean_annual_chla),
    log_settlement_grav = log(settlement_grav),
    log_settlement_pop = log(settlement_pop),
    log_market_gravity = log(market_gravity),
    connectivity = connectivity,
    mpa_status = factor(mpa_status, levels  = c("none", "low", "medium"),
                        ordered = FALSE))

predictor_order_tran <- c("rugosity", "log_settlement_grav", "log_market_gravity",
                          "log_settlement_pop", "log_chla", "connectivity")

# Plot transformed histograms
( p_pred_tran <- transformed_predictors %>%
    dplyr::select(all_of(predictor_order_tran)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, predictor_order_tran, predictor_labels)) %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
    facet_wrap(~ variable, scales = "free") +
    labs(title = "Transformed predictors") + theme_bw() )
# Much healthier distributions

## Standardise continuous predictors ----
# z-score scaling (mean = 0, SD = 1) enables direct comparison
# of effect sizes across predictors with different units
scaled_predictors <- transformed_predictors %>%
  transmute(
    site = site,
    ecoregion = ecoregion,
    rugosity_sc = as.numeric(scale(rugosity)),
    log_settlement_grav_sc = as.numeric(scale(log_settlement_grav)),
    log_settlement_pop_sc = as.numeric(scale(log_settlement_pop)),
    log_market_gravity_sc = as.numeric(scale(log_market_gravity)),
    log_chla_sc = as.numeric(scale(log_chla)),
    connectivity_sc = as.numeric(scale(connectivity)),
    mpa_status = mpa_status)
# ──────────────────────────────────────────────────────────────────────────────

# Predictor correlation checks ----
# Check for blocking collinearity before modelling.
# Note that r > 0.70 requires caution if modeeling these together

corr_matrix <- scaled_predictors %>%
  mutate(mpa_numeric = as.numeric(mpa_status)) %>%
  dplyr::select(ends_with("_sc"), mpa_numeric) %>%
  rename(
    "Market gravity" = log_market_gravity_sc,
    "Settlement gravity" = log_settlement_grav_sc,
    "Settlement pop." = log_settlement_pop_sc,
    "Chlorophyll-a" = log_chla_sc,
    "Rugosity" = rugosity_sc,
    "Connectivity" = connectivity_sc,
    "MPA status" = mpa_numeric
  ) %>%
  cor(use = "complete.obs")

vars <- colnames(corr_matrix)

corr_long <- as.data.frame(as.table(corr_matrix)) %>%
  rename(x = Var1, y = Var2, corr = Freq) %>%
  mutate(
    x = factor(x, levels = vars),
    y = factor(y, levels = rev(vars))
  ) %>%
  filter(as.numeric(x) < length(vars) - as.numeric(y) + 1)

ggplot(corr_long, aes(x = x, y = y, fill = corr)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", corr)),
            size = 3.2, color = "black") +
  scale_fill_gradient2(
    low = "#2c7bb6",
    mid = "white",
    high = "#d7191c",
    limits = c(-1, 1),
    name = "Correlation"
  ) +
  scale_x_discrete(limits = vars, drop = FALSE, position = "bottom") +
  scale_y_discrete(limits = rev(vars), drop = FALSE, position = "left") +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1,
                                     vjust = 1, size = 10),
    axis.text.y = element_text(hjust = 1, size = 9.5),
    legend.position = c(1.015, 0.48),
    legend.key.height = unit(1.75, "cm"),
    legend.key.width = unit(0.5,  "cm"),
    legend.title = element_text(size = 9.75, hjust = 1.5,
                                     vjust = 1.5, colour = "grey40"),
    legend.text = element_text(size = 9, vjust = 0.3,
                                     colour = "grey40"),
    plot.margin = margin(10, 10, 10, 10)
  )

# Collinearity notes
#
## Gravity metrics r = 0.66–0.54
#   Moderate positive, but this was expected as they all proxy the same construct
#   Only one used per model so no issues here
#
# Chlorophyll-a vs settlement gravity: r = -0.57
#   Strongest pairwise correlation in the dataset (that will be included in the same model)
#   Sites with higher chlorophyll-a (often a proxy for productivity) tend to have lower settlement gravity.
#   Note that this is typically for coastal sites.
#   Negative — productive sites tend to be more remote and less fished. 
#   Monitor chla coefficient in models including settlement gravity simultaneously.
#
# MPA vs settlement pop.: r = -0.41
#   Negative association with population-based metric than with settlement gravity 
#   Consistent again with MPAs placed away from densely populated areas.
#
# CONNECTIVITY vs market gravity: r = 0.30
#   Weak positive — more connected sites tend to be near
#   larger markets. Acceptable; monitor in models including
#   both terms.
#
# CONNECTIVITY vs chla: r = -0.30
#   Weak negative — more connected sites tend to be less
#   productive. No inference concern.
#
# RUGOSITY: all |r| ≤ 0.19 — orthogonal to all predictors.
# CONNECTIVITY: all |r| ≤ 0.34 — no blocking collinearity.
#
# No blocking collinearity. All pairwise |r| < 0.60.
# Gravity metrics handled by single-metric-per-model rule.


# ============================================================
#  ANALYSIS DATASETS
# ============================================================

final_predictors <- scaled_predictors %>%
  dplyr::select(
    site,
    ecoregion,                  # ADD THIS
    rugosity_sc,                # BASELINE    — habitat complexity
    log_settlement_grav_sc,     # Q1/Q2/Q3   — primary pressure metric
    log_chla_sc,                # BASELINE    — productivity
    connectivity_sc,            # Q2/Q3      — larval connectivity
    mpa_status,                 # Q2/Q3      — governance modifier
    log_settlement_pop_sc,      # SENSITIVITY
    log_market_gravity_sc       # SENSITIVITY
  )

# ── Transect-level dataset ────────────────────────────────────
transect_model_data <- total_transects %>%
  left_join(final_predictors, by = "site") %>%
  mutate(log_transect_biomass = log(transect_total_biomass))

cat("\nTransect data:", nrow(transect_model_data), "rows,",
    n_distinct(transect_model_data$site), "sites\n")
cat("Biomass zeros:",
    sum(transect_model_data$transect_total_biomass == 0), "\n")

# ── Site-level dataset ────────────────────────────────────────
total_model_data <- transect_model_data %>%
  group_by(site, ecoregion) %>%
  summarise(
    log_mean_biomass       = log(mean(transect_total_biomass,
                                      na.rm = TRUE)),
    mean_biomass           = mean(transect_total_biomass,
                                  na.rm = TRUE),
    n_transects            = n(),
    rugosity_sc            = first(rugosity_sc),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_chla_sc            = first(log_chla_sc),
    connectivity_sc        = first(connectivity_sc),
    mpa_status             = first(mpa_status),
    log_settlement_pop_sc  = first(log_settlement_pop_sc),
    log_market_gravity_sc  = first(log_market_gravity_sc),
    .groups = "drop"
  ) %>%
  mutate(
    site    = as.factor(site),
    ecoregion = as.factor(ecoregion)
  )

cat("\nSite data:", nrow(total_model_data), "sites,",
    n_distinct(total_model_data$ecoregion), "ecoregions\n")

# ── Data checks ───────────────────────────────────────────────
total_model_data %>%
  dplyr::select(site, ecoregion, rugosity_sc, log_settlement_grav_sc,
                log_chla_sc, connectivity_sc, mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

cat("\nZeros in mean_biomass:",
    sum(total_model_data$mean_biomass == 0), "\n")
cat("-Inf in log_mean_biomass:",
    sum(is.infinite(total_model_data$log_mean_biomass)), "\n")
cat("NAs in log_mean_biomass:",
    sum(is.na(total_model_data$log_mean_biomass)), "\n")
cat("\nResponse variable summary:\n")
print(summary(total_model_data$log_mean_biomass))

# ── MPA classification check ──────────────────────────────────
# Verify MPA counts after any raw data corrections.
cat("\nMPA status counts:\n")
print(table(total_model_data$mpa_status))


# ── Variance inflation factors ────────────────────────────────
# Confirms pairwise correlations do not translate into
# meaningful variance inflation in multivariate models.
# Checked at each stage predictors are added.

cat("\n--- VIF: baseline ---\n")
vif(lm(log_mean_biomass ~ rugosity_sc +
         log_chla_sc,
       data = total_model_data))

cat("\n--- VIF: baseline + pressure ---\n")
vif(lm(log_mean_biomass ~ rugosity_sc +
         log_chla_sc +
         log_settlement_grav_sc,
       data = total_model_data))

cat("\n--- VIF: baseline + pressure + connectivity ---\n")
vif(lm(log_mean_biomass ~ rugosity_sc +
         log_chla_sc +
         log_settlement_grav_sc +
         connectivity_sc,
       data = total_model_data))

cat("\n--- VIF: baseline + pressure + MPA ---\n")
vif(lm(log_mean_biomass ~ rugosity_sc +
         log_chla_sc +
         log_settlement_grav_sc +
         mpa_status,
       data = total_model_data))

#  Collinearity:
#  Pairwise correlations checked via correlation matrix.
#  VIF checked at three stages:
#    Baseline (rugosity + chla):
#      rugosity VIF = 1.000, chla VIF = 1.000
#      Baseline predictors completely orthogonal.
#    Baseline + pressure:
#      rugosity VIF = 1.011, chla VIF = 1.427,
#      settlement gravity VIF = 1.438
#    Baseline + pressure + connectivity:
#      rugosity VIF = 1.012, chla VIF = 1.486,
#      settlement gravity VIF = 1.446,
#      connectivity VIF = 1.090
#    Baseline + pressure + MPA:
#      rugosity GVIF^(1/(2*Df)) = 1.011 — no concern
#      chla GVIF^(1/(2*Df)) = 1.259 — no concern
#      settlement gravity GVIF^(1/(2*Df)) = 1.377 — no concern
#      mpa_status GVIF^(1/(2*Df)) = 1.073 — no concern
#      (GVIF reported for factor predictors — equivalent
#      to VIF for continuous predictors when Df = 1)
#  All VIFs/GVIFs < 2.0 at every stage — no multicollinearity
#  concern throughout Q1–Q3 model sequence.

# ============================================================
#  MODEL FAMILY SELECTION
#  Baseline sufficient to assess error distribution before human terms introduced. 
#  If family holds for baseline it holds for all models built on it.
# ============================================================

lm_gaussian_raw <- lm(mean_biomass ~ rugosity_sc +
                        log_chla_sc,
                      data = total_model_data)

par(mfrow = c(2, 2))
plot(lm_gaussian_raw, main = "Gaussian raw")

lm_gaussian_log <- lm(log_mean_biomass ~ rugosity_sc +
                        log_chla_sc,
                      data = total_model_data)

plot(lm_gaussian_log, main = "Gaussian log")

glm_gamma <- glm(mean_biomass ~ rugosity_sc +
                   log_chla_sc,
                 family = Gamma(link = "log"),
                 data   = total_model_data)

plot(glm_gamma, main = "Gamma log-link")

# Family selection --------------------------------------------
#
#  Gaussian (raw):   rejected — severe heteroscedasticity,
#                    U-shaped residuals vs fitted, heavy
#                    upper Q-Q tail. Sites 6, 48, 19
#                    influential.
#
#  Gamma (log link): rejected — systematic Q-Q deviation
#                    at upper tail despite flat residuals.
#                    Sites 6, 48, 19 pull from half-normal
#                    line.
#
#  Gaussian (log):   SELECTED — broadly flat residuals,
#                    Q-Q closely follows theoretical line
#                    across full range, minor upper tail
#                    deviation only. Scale-location flat.
#                    No sites exceed Cook's distance
#                    threshold.
#
#  Proceed: lm() on log_mean_biomass throughout.

# ============================================================
#  RANDOM EFFECT STRUCTURE
#  Tested on baseline model — adding unvalidated predictors
#  would presuppose Q1-Q3 outcomes and absorb between-
#  ecoregion variance, biasing the test.
#  ML (REML = FALSE) for AICc comparison.
# ============================================================

re_null <- lm(log_mean_biomass ~ rugosity_sc +
                         log_chla_sc,
                       data = total_model_data)

re_ecoregion <- lmer(log_mean_biomass ~ rugosity_sc +
                       log_chla_sc +
                       (1 | ecoregion),
                     data = total_model_data,
                     REML = FALSE)

cat("\n--- Random effect structure: ecoregion ---\n")
print(make_aicc_df(list(
  "No RE"           = re_null,
  "(1 | ecoregion)" = re_ecoregion
)))

# Random effect structure: ecoregion
# Tested on baseline model to avoid circularity —
# predictor inclusion not yet established at this stage.
# Ecoregion RE not supported (DAICc = 2.25, weight = 0.245).
# Not pursued for two additional reasons:
# (1) Only 4 ecoregions — insufficient group-level
#     replication for reliable variance component
#     estimation (Gelman & Hill 2007).
# (2) Severely uneven group sizes:
#     Kenya-Tanzania north:       2 sites
#     Comoros:                    8 sites
#     Madagascar:                 9 sites
#     Tanzania south-Mozambique: 35 sites
#     Kenya-Tanzania north (n = 2) cannot support a
#     meaningful group-level intercept estimate.
# Between-ecoregion variation acknowledged as limitation.
# All models fitted as lm() throughout.

# ============================================================
#  Q1 — HUMAN PRESSURE
#
#  Scientific question:
#  Does human pressure explain variation in reef fish biomass beyond local 
#  ecological context, and which spatial metric best captures SSF exploitation 
#  intensity in the WIO?
#
#  A priori prediction:
#  Settlement gravity will outperform market gravity because residential 
#  proximity better captures the spatial footprint of subsistence-oriented 
#  exploitation than commercial market access in this SSF-dominated system
#  (Cinner et al. 2016, Samoilys et al. 2019).
#
#  Baseline (fixed a priori):
#  log(biomass) ~ rugosity_sc + log_chla_sc
#  Justified by ...
#  Retained regardless of significance.
#
#  Two steps:
#  Step 1 — Does pressure add beyond baseline, and which metric best captures SSF exploitation?
#            AICc weight criterion for metric selection.
#  Step 2 — Coefficient check for selected metric.
#            p-value criterion for effect confirmation.
# ============================================================

m_baseline <- lm(log_mean_biomass ~ rugosity_sc + log_chla_sc,
                 data = total_model_data)
summary(m_baseline)

# ── Q1 Step 1: Does pressure add? Which metric is best? ──────
# Baseline vs baseline + each pressure metric.
# Single comparison answers both questions simultaneously:
# whether pressure matters AND which metric captures it best.
# Best metric: highest AICc weight AND outperforms baseline.
#
# Three candidate metrics:
# Settlement gravity — weights populations by proximity,
#   captures subsistence exploitation footprint (primary)
# Market gravity — weights population centres by accessibility, 
# validated globally by Cinner et al. 2016
# Settlement population — aggregate population size,
#   sensitivity check only, not in primary table

q1_settgrav <- lm(log_mean_biomass ~ rugosity_sc +
                    log_chla_sc +
                    log_settlement_grav_sc,
                  data = total_model_data)

q1_mktgrav  <- lm(log_mean_biomass ~ rugosity_sc +
                    log_chla_sc +
                    log_market_gravity_sc,
                  data = total_model_data)

q1_settpop  <- lm(log_mean_biomass ~ rugosity_sc +
                    log_chla_sc +
                    log_settlement_pop_sc,
                  data = total_model_data)

cat("\n--- Q1 Step 1: Metric comparison ---\n")
print(make_aicc_df(list(
    "Baseline"                      = m_baseline,
    "Baseline + settlement gravity" = q1_settgrav,
    "Baseline + settlement pop."    = q1_settpop,
    "Baseline + market gravity"     = q1_mktgrav
)))

# Results:
#   Settlement gravity: AICc = 101.36, weight = 0.826
#   Baseline:           DAICc = 4.39,  weight = 0.092
#   Market gravity:     DAICc = 5.86,  weight = 0.044
#   Settlement pop.:    DAICc = 6.17,  weight = 0.038
#
#   Settlement gravity decisively selected (weight = 0.826).
#   Outperforms ecological baseline alone (DAICc = 4.39),
#   confirming human pressure explains meaningful variance
#   beyond habitat and productivity context.
#   Market gravity and settlement population perform
#   similarly poorly and are not competitive alternatives
#   (combined weight = 0.082).
#   Residential proximity weighted by distance better
#   captures SSF exploitation intensity than either
#   market access or aggregate population size alone.

# ── Q1 Step 2: Coefficient check — selected metric ───────────
# Confirms direction and significance of selected metric.
# p-values are the criterion here — is the pressure effect
# real and in the predicted direction?
# Also confirms stability of baseline predictors across
# models — rugosity and chla should remain consistent.

# --- Model summary: coefficients, SE, t, p ---
cat("\n--- Q1 Step 2: Selected metric coefficients ---\n")
summary(q1_settgrav)

# --- Settlement gravity: back-transform effect size across observed range ---
cat("\n--- Settlement gravity: standardised range ---\n")
sg_range <- range(total_model_data$log_settlement_grav_sc, na.rm = TRUE)
print(sg_range)
sg_span <- diff(sg_range)
cat(sprintf("Span: %.3f SD units\n", sg_span))
b_sg <- -0.251 # coefficient estimate

# Calculate fold difference in biomass across observed settlement gravity range
fold_diff_sg <- exp(abs(b_sg * sg_span))
cat(sprintf("Fold difference in biomass (low vs high pressure): %.2fx\n", fold_diff_sg))

# --- Rugosity: back-transform effect size across observed range ---
cat("\n--- Rugosity: standardised range ---\n")
rug_range <- range(total_model_data$rugosity_sc, na.rm = TRUE)
rug_span  <- diff(rug_range)
cat(sprintf("Span: %.3f SD units\n", rug_span))
b_rug <- 0.216 # coefficient estimate

# Calculate fold difference in biomass across observed rugosity range
fold_diff_rug <- exp(abs(b_rug * rug_span))
cat(sprintf("Fold difference in biomass (low vs high rugosity): %.2fx\n", fold_diff_rug))

# --- 95% CIs: uncertainty around key coefficient estimates ---
cat("\n--- 95% CIs: best Q1 model ---\n")
print(confint(q1_settgrav))

# Results:
#   Settlement gravity: b = -0.251, t(50) = -2.596, p = 0.012, 95% CI [-0.445, -0.057] *
#     Fold difference in biomass (low vs high pressure): 2.74x
#     [span = 4.001 SD units; exp(0.251 × 4.001) = 2.74]
#   Rugosity:           b = +0.216, t(50) =  2.707, p = 0.009, 95% CI [ 0.056,  0.376] **
#     Fold difference in biomass (low vs high rugosity): 3.04x
#     [span = 5.151 SD units; exp(0.216 × 5.151) = 3.04]
#   Chla:               b = -0.125, t(50) = -1.275, p = 0.208, 95% CI [-0.322,  0.072] ns
#     Shifts negative as expected — ecological covariation
#     with pressure resolved once both included simultaneously
#     (r = -0.54, confirmed non-spatial in diagnostics)
#   Model: F(3,50) = 5.253, p = 0.003, R² = 0.240, adj. R² = 0.194

# ── Direction check for market gravity ───────────────────────
# Confirms pressure signal consistent in direction across
# metrics — failure to reach significance reflects metric
# choice, not absence of a pressure effect.
# Full sensitivity in Sensitivity (a).

cat("\n--- Q1: Market gravity direction check ---\n")
summary(q1_mktgrav)$coefficients["log_market_gravity_sc", ]

# Results:  b = -0.086, p = 0.346 ns
#   Consistent direction with settlement gravity
#   but effect size one third as large (-0.086 vs -0.251)
#   and not significant. Confirms settlement gravity
#   is the correct metric choice — the pressure signal
#   is detectable with residential proximity but not
#   with commercial market access in this SSF system.
#   Kept for sensititvity later. 

# ── Best Q1 model ─────────────────────────────────────────────
# rugosity + chla + settlement gravity
# Reference model for Q2 and Q3.
m_best_q1 <- q1_settgrav


# ============================================================
#  Q2 — LARVAL CONNECTIVITY
#
#  Scientific question:
#  Does larval connectivity explain additional variation in
#  reef fish biomass beyond the human pressure baseline, and
#  does it modify the relationship between fishing pressure
#  and biomass?
#
#  Directly extends Warmuth et al. (2024) — their models
#  included connectivity without fishing pressure as a
#  covariate. This analysis tests whether the connectivity
#  signal persists once pressure is controlled.
#
#  Two steps:
#  Step 1 — connectivity as main effect (additive)
#  Step 2 — connectivity × pressure interaction
#
#  Criterion for both: AICc weight for model support,
#  p-value for effect confirmation if supported.
# ============================================================

# ── Q2 Step 1: Does connectivity add beyond pressure? ─────────
# Tests whether sites that are more connected to the broader reef network have 
# higher biomass independently of local habitat and fishing pressure.

m_q2_conn <- lm(log_mean_biomass ~ rugosity_sc +
                  log_chla_sc +
                  log_settlement_grav_sc +
                  connectivity_sc,
                data = total_model_data)

cat("\n--- Q2 Step 1: Connectivity main effect ---\n")
print(make_aicc_df(list(
  "Best Q1"              = m_best_q1,
  "Best Q1 + conn"       = m_q2_conn
)))

cat("\n--- Q2 Step 1: Connectivity coefficients ---\n")
summary(m_q2_conn)

# Q2 Step 1 results:
#   Best Q1:         AICc = 101.36, weight = 0.770
#   Best Q1 + conn:  DAICc = 2.41,  weight = 0.230
#   Connectivity not supported as main effect.
#
#   Connectivity: b = +0.028, p = 0.739 — near zero,
#     not significant. No independent effect on biomass
#     once habitat and pressure controlled.
#   Rugosity:     b = +0.216, p = 0.010 — stable
#   Pressure:     b = -0.253, p = 0.013 — stable
#
#   Directly extends Warmuth et al. (2024): connectivity
#   signal on herbivore abundance does not translate to
#   total biomass once fishing pressure included as
#   covariate. Larval supply does not independently
#   subsidise total reef fish biomass at this spatial
#   scale in this SSF-dominated system.


# ── Q2 Step 2: Does connectivity modify pressure effects? ─────
# Tests whether the negative pressure-biomass relationship
# depends on connectivity — specifically whether well-connected
# sites sustain higher biomass under fishing pressure through
# larval replenishment.
#
# A priori hypothesis: connectivity × pressure interaction
# positive — high connectivity buffers exploitation effects.
# Interaction not supported if: DAICc > 2 vs best Q1 AND/OR
# direction inconsistent with hypothesis.
#
# Note: interaction model tested against best Q1 (not Q2
# main effect model) because if connectivity main effect
# not supported, reference should be the best supported
# model not a poorly supported intermediate.

m_q2_conn_int <- lm(log_mean_biomass ~ rugosity_sc +
                      log_chla_sc +
                      connectivity_sc *
                      log_settlement_grav_sc,
                    data = total_model_data)

cat("\n--- Q2 Step 2: Connectivity x pressure interaction ---\n")
print(make_aicc_df(list(
  "Best Q1"                  = m_best_q1,
  "Best Q1 + conn x pressure" = m_q2_conn_int
)))

cat("\n--- Q2 Step 2: Interaction coefficients ---\n")
summary(m_q2_conn_int)

# --- 95% CIs: uncertainty around key coefficient estimates ---
cat("\n--- 95% CIs: Connectivity x pressure interaction ---\n")
print(confint(m_q2_conn_int))



# Expected:
#   Interaction not supported (DAICc > 2 vs best Q1)
#   Interaction term: b = -0.225, p = 0.060 — marginal
#   and in WRONG direction — negative sign suggests
#   connectivity amplifies rather than buffers pressure,
#   contrary to hypothesis. Not interpretable as
#   buffering mechanism.

# ── Best Q2 model ─────────────────────────────────────────────
# Connectivity not supported as main effect or interaction.
# Best Q2 model = best Q1 model (unchanged).
# Reference model for Q3.
m_best_q2 <- m_best_q1


# ============================================================
#  Q3 — FORMAL PROTECTION
#
#  Scientific question:
#  Does MPA status explain additional variation in reef fish
#  biomass beyond the fully specified pressure and
#  connectivity model?
#
#  MPA tested last because MPAs are non-randomly placed —
#  preferentially in lower-pressure areas (r = -0.30 with
#  settlement gravity). Effect confounded by pressure,
#  habitat, and connectivity unless all three controlled.
#
#  Criterion: AICc weight for model support.
#  p-values for coefficient direction if supported.
#
#  Interpretation note: MPA variable captures nominal
#  protection designation only — does not incorporate
#  enforcement or compliance information. Null result
#  should be interpreted as nominal designation not
#  predicting biomass, not as protection having no effect
#  (Cinner et al. 2016).
# ============================================================

m_q3_mpa <- lm(log_mean_biomass ~ rugosity_sc +
                 log_chla_sc +
                 log_settlement_grav_sc +
                 mpa_status,
               data = total_model_data)

cat("\n--- Q3: MPA main effect ---\n")
print(make_aicc_df(list(
  "Best Q2"        = m_best_q2,
  "Best Q2 + MPA"  = m_q3_mpa
)))

cat("\n--- Q3: MPA coefficients ---\n")
summary(m_q3_mpa)

cat("\n--- Q3: MPA model confint ---\n")
print(confint(m_q3_mpa))

# Q3 results:
#   Best Q2:       AICc = 101.36, weight = 0.819
#   Best Q2 + MPA: DAICc = 3.01,  weight = 0.182
#   MPA status not supported beyond pressure baseline.
#
#   Low MPA:    b = -0.342, p = 0.171 ns
#   Medium MPA: b = -0.084, p = 0.669 ns
#   Neither protection level significantly predicts
#   biomass once habitat and pressure controlled.
#   Negative coefficients reflect preferential MPA
#   placement in lower-pressure areas (r = -0.30 with
#   settlement gravity) — not a genuine negative effect
#   of protection.
#
#   Rugosity: b = +0.211, p = 0.011 — stable throughout
#   Pressure: b = -0.239, p = 0.036 — stable throughout
#
#
#   Best model throughout: rugosity + chla +
#   settlement gravity (m_best_q1 = m_best_q2 = m_best_q3)

# ============================================================
#  SPATIAL AUTOCORRELATION CHECK
#  Residuals from best model tested for spatial structure.
#  Reported as diagnostic — not part of model selection.
# ============================================================

site_coords <- location_data %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(lon = first(longitude),
            lat = first(latitude),
            .groups = "drop")

total_model_data_coords <- total_model_data %>%
  left_join(site_coords, by = "site")

coords_mat <- cbind(total_model_data_coords$lon,
                    total_model_data_coords$lat)
listw5 <- nb2listw(knn2nb(knearneigh(coords_mat, k = 5)),
                   style = "W")

cat("\n--- Spatial autocorrelation: best model residuals ---\n")
print(moran.test(residuals(m_best_q2), listw5))

# Spatial autocorrelation: best model residuals
# Moran's I = 0.140, p = 0.015 — weak but significant.
#
# Spatial error modelling not pursued for two reasons:
# (1) Ecoregion RE not supported (DAICc = 2.25 at
#     baseline) and only 4 groups with severely uneven
#     sizes (n = 2, 8, 9, 35) — between-ecoregion
#     clustering not correctable with available data.
# (2) Discontinuous sampling design — k-NN spatial
#     weights bridge across geographically isolated
#     country clusters, producing a weights matrix
#     that does not reflect true within-region spatial
#     covariance structure (Kissling & Carl 2008,
#     Dormann et al. 2007).
#
# Weak residual autocorrelation acknowledged as a
# limitation. May slightly inflate type I error rates
# for pressure and rugosity coefficients. Substantive
# conclusions unlikely to change given effect sizes
# and consistency across model specifications.

# ============================================================
#  SENSITIVITY ANALYSIS
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
# Mirrors baseline structure — only pressure metric swapped.
# Purpose: robustness check, not model selection.
# AICc comparison not reported here — see Q1 for metric
# selection. 
# Coefficients confirm conclusions are not metric-dependent.

sens_mktgrav <- lm(log_mean_biomass ~ rugosity_sc +
                     log_market_gravity_sc +
                     log_chla_sc,
                   data = total_model_data)

sens_settpop  <- lm(log_mean_biomass ~ rugosity_sc +
                      log_settlement_pop_sc +
                      log_chla_sc,
                    data = total_model_data)

cat("\n--- Sensitivity (a): alternative pressure metrics ---\n")
cat("Market gravity:\n")
print(summary(sens_mktgrav)$coefficients)
cat("\nSettlement population:\n")
print(summary(sens_settpop)$coefficients)

# Sensitivity (a): alternative pressure metrics
# Purpose: confirm Q1 conclusions are not metric-dependent.
# Both alternatives substituted into baseline structure.
# AICc comparison not repeated — metric selected in Q1.
# Coefficients confirm direction consistent across metrics.
#
# Market gravity:       b = -0.086, p = 0.346 ns
# Settlement pop.:      b = -0.067, p = 0.439 ns
#
# Both negative — pressure signal consistent in direction.
# Neither significant — pressure only detectable with
# settlement gravity. Effect size approximately one third
# (market) and one quarter (population) of settlement
# gravity — confirms residential proximity weighted by
# distance captures SSF exploitation footprint more
# precisely than market access or aggregate population.
# Rugosity significant and stable across all three metrics
# (b = 0.220-0.225, p < 0.05) — habitat signal robust
# to pressure metric choice.

# ============================================================
#  SENSITIVITY ANALYSIS (b) — TRANSECT-LEVEL REPLICATION
#
#  Purpose: confirm that site-level Q1-Q3 conclusions are
#  not an artefact of spatial aggregation. Replicates the
#  Q1-Q3 model sequence at transect level using lmer with
#  a site-level random intercept.
#
#  (1 | site) accounts for non-independence of transects
#  within sites — conventional choice for Gaussian mixed
#  models with simple nested structure.
#
#  If model ordering and coefficient directions are
#  consistent with site-level results, site-level
#  aggregation does not alter qualitative inference.
#
#  Models fitted with REML = TRUE for coefficient
#  inference. Refitted with REML = FALSE for AICc
#  comparison — REML likelihoods not valid for comparing
#  models with different fixed effect structures.
# ============================================================

# ── Fit transect models (REML = TRUE) ────────────────────────

# Null — site clustering only, no predictors
sens_t_null <- lmer(log_transect_biomass ~ 1 +
                      (1 | site),
                    data = transect_model_data,
                    REML = TRUE)

# Baseline — rugosity + chla fixed a priori
# Mirrors site-level baseline
sens_t_baseline <- lmer(log_transect_biomass ~ rugosity_sc +
                          log_chla_sc +
                          (1 | site),
                        data = transect_model_data,
                        REML = TRUE)

# Baseline + pressure — mirrors Q1
# Does settlement gravity add beyond ecological baseline?
sens_t_pressure <- lmer(log_transect_biomass ~ rugosity_sc +
                          log_chla_sc +
                          log_settlement_grav_sc +
                          (1 | site),
                        data = transect_model_data,
                        REML = TRUE)

# Best + connectivity — mirrors Q2 main effect
# Does connectivity add beyond pressure baseline?
sens_t_conn <- lmer(log_transect_biomass ~ rugosity_sc +
                      log_chla_sc +
                      log_settlement_grav_sc +
                      connectivity_sc +
                      (1 | site),
                    data = transect_model_data,
                    REML = TRUE)

# Best + MPA — mirrors Q3
# Does MPA add beyond pressure baseline?
sens_t_mpa <- lmer(log_transect_biomass ~ rugosity_sc +
                     log_chla_sc +
                     log_settlement_grav_sc +
                     mpa_status +
                     (1 | site),
                   data = transect_model_data,
                   REML = TRUE)


# ── Diagnostics on baseline model ────────────────────────────
# Run on true baseline (rugosity + chla) — consistent with
# site-level approach. Confirms Gaussian lmer appropriate
# before introducing pressure terms.

par(mfrow = c(1, 2))
plot(fitted(sens_t_baseline), residuals(sens_t_baseline),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Transect baseline: Residuals vs Fitted")
abline(h = 0, lty = 2, col = "grey60")
lines(lowess(fitted(sens_t_baseline),
             residuals(sens_t_baseline)), col = "red")
qqnorm(residuals(sens_t_baseline),
       main = "Transect baseline: Q-Q Residuals")
qqline(residuals(sens_t_baseline), col = "red")
par(mfrow = c(1, 1))


# ── AICc comparison — ML refits ──────────────────────────────
sens_t_null_ml     <- update(sens_t_null,     REML = FALSE)
sens_t_baseline_ml <- update(sens_t_baseline, REML = FALSE)
sens_t_pressure_ml <- update(sens_t_pressure, REML = FALSE)
sens_t_conn_ml     <- update(sens_t_conn,     REML = FALSE)
sens_t_mpa_ml      <- update(sens_t_mpa,      REML = FALSE)

cat("\n--- Sensitivity (b): transect-level model comparison ---\n")
print(make_aicc_df(list(
  "Null"                = sens_t_null_ml,
  "Baseline"            = sens_t_baseline_ml,
  "Baseline + pressure" = sens_t_pressure_ml,
  "Best + conn"         = sens_t_conn_ml,
  "Best + MPA"          = sens_t_mpa_ml
)))

# Sensitivity (b) results:
#
# Diagnostics — transect baseline:
#   Residuals vs fitted: minor negative trend in lowess
#     line but acceptable for n = 243 transects.
#   Q-Q: excellent — points follow theoretical line
#     closely, minor lower tail deviations only.
#   Gaussian lmer structure confirmed appropriate.
#
# AICc comparison (ML):
#   Baseline + pressure: AICc = 634.93, weight = 0.626
#   Best + conn:         DAICc = 2.10,  weight = 0.219
#   Best + MPA:          DAICc = 3.08,  weight = 0.134
#   Baseline:            DAICc = 6.84,  weight = 0.021
#   Null:                DAICc = 13.05, weight = 0.001
#
#   Model ordering identical to site-level Q1-Q3.
#   Pressure model best supported at both levels of
#   analysis — site-level aggregation does not alter
#   qualitative inference.
#   Connectivity not supported (DAICc = 2.10) —
#     consistent with Q2.
#   MPA not supported (DAICc = 3.08) —
#     consistent with Q3.


# ── Coefficient summary — REML-fitted pressure model ─────────
# Primary interest: direction and magnitude of pressure
# and rugosity coefficients relative to site-level results.
# Chla expected to reach significance at transect level
# (n = 243 vs n = 54) due to greater statistical power —
# direction should be consistent.

cat("\n--- Sensitivity (b): pressure model coefficients (REML) ---\n")
summary(sens_t_pressure)

# ── Diagnostics on pressure model (supplementary) ────────────
# Confirms Gaussian lmer remains appropriate after addition
# of settlement gravity. Reported in supplementary materials.

par(mfrow = c(1, 2))
plot(fitted(sens_t_pressure), residuals(sens_t_pressure),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Transect pressure model: Residuals vs Fitted")
abline(h = 0, lty = 2, col = "grey60")
lines(lowess(fitted(sens_t_pressure),
             residuals(sens_t_pressure)), col = "red")
qqnorm(residuals(sens_t_pressure),
       main = "Transect pressure model: Q-Q Residuals")
qqline(residuals(sens_t_pressure), col = "red")
par(mfrow = c(1, 1))

# Sensitivity (b): pressure model coefficients (REML)
#
# Fixed effects:
#   Rugosity:           b = +0.240, t = 3.100 — significant, stable
#   Chla:               b = -0.205, t = -2.187 — significant at
#     transect level (n = 243) as predicted — greater power
#     resolves effect that was marginal at site level (p = 0.208)
#   Settlement gravity: b = -0.285, t = -3.006 — significant,
#     consistent direction with site-level result (b = -0.251)
#
# Random effects:
#   Site variance:     0.168
#   Residual variance: 0.651
#   ICC = 0.205 — 20.5% of residual variance attributable to
#     between-site differences — (1 | site) random intercept
#     justified (ICC well above 0.10 threshold)
#
# Qualitative conclusions unchanged from site-level analysis:
#   pressure negative and significant, rugosity positive and
#   significant, chla now reaches significance as predicted.

# ── ICC calculation ───────────────────────────────────────────
vc <- as.data.frame(VarCorr(sens_t_pressure))
site_var     <- vc$vcov[vc$grp == "site"]
residual_var <- vc$vcov[vc$grp == "Residual"]
icc          <- site_var / (site_var + residual_var)

cat(sprintf("\nICC = %.3f — %.1f%% of variance attributable
to between-site differences beyond fixed predictors\n",
            icc, icc * 100))

# Results
#
# Coefficient summary (REML, n = 243 transects, 54 sites):
#   Rugosity:           b = +0.240, t = 3.100 **
#     Consistent with site-level (b = +0.216) —
#     habitat signal stable across aggregation levels.
#   Settlement gravity: b = -0.285, t = -3.006 **
#     Consistent with site-level (b = -0.251) —
#     pressure signal robust to aggregation.
#   Chla:               b = -0.205, t = -2.187 *
#     Reaches significance at transect level (p < 0.05)
#     but not at site level (b = -0.125, p = 0.208).
#     Greater statistical power (n = 243 vs 54) allows
#     detection of a weaker signal. Direction consistent
#     across both levels — effect is real but site-level
#     sample is underpowered to detect it clearly.
#
# Random effects:
#   Site variance:    0.168 (SD = 0.410)
#   Residual variance: 0.651 (SD = 0.807)
#   ICC = 0.205 — 20.5% of variance attributable to
#   between-site differences beyond fixed predictors.
#   Confirms (1 | site) random intercept is justified.
#
# Overall conclusion:
#   Direction and magnitude of rugosity and pressure
#   coefficients consistent across site and transect
#   levels. Chla direction consistent, stronger signal
#   at transect level due to power difference only.
#   Site-level aggregation does not alter qualitative
#   inference for any Q1-Q3 conclusion.
# ============================================================
#  MARGINAL EFFECT PLOTS
#  Best model: rugosity + chla + settlement gravity
#  All predictors held at mean (z = 0) except focal var.
#  Observed data overlaid on fitted lines.
#  Chla non-significant at site level (p = 0.208) but
#  retained as baseline control — reaches significance
#  at transect level (b = -0.205, t = -2.187), direction
#  consistent across both levels.
# ============================================================

best_model <- m_best_q1  # rugosity + chla + settlement gravity

p_rugosity <- plot_effect(best_model, total_model_data,
                          "rugosity_sc",
                          "Rugosity (standardised)")

p_pressure <- plot_effect(best_model, total_model_data,
                          "log_settlement_grav_sc",
                          "log(Settlement gravity) (standardised)")

p_chla     <- plot_effect(best_model, total_model_data,
                          "log_chla_sc",
                          "log(Chlorophyll-a) (standardised)")

gridExtra::grid.arrange(p_rugosity, p_pressure, p_chla,
                        ncol = 3)

# jpeg("marginal_effects_best_model.jpg",
#      width = 33, height = 11, units = "cm", res = 300)
# gridExtra::grid.arrange(p_rugosity, p_pressure, p_chla,
#                         ncol = 3)
# dev.off()

# ============================================================
#  RESULTS SUMMARY
#  Quick reference for writing — not part of analysis.
#  Verify all values match reported results before writing.
# ============================================================

results_summary <- tribble(
  ~Question,      ~Result,          ~Key_finding,
  "Q1 metric",    "Sett. gravity",  "weight = 0.826, DAICc = 4.39 vs baseline",
  "Q1 effect",    "Significant",    "b = -0.251, p = 0.012",
  "Q1 rugosity",  "Significant",    "b = +0.216, p = 0.009 — stable throughout",
  "Q2 main",      "Not supported",  "DAICc = 2.41, weight = 0.230, b = +0.028, p = 0.739",
  "Q2 int",       "Not supported",  "DAICc = 1.06, weight = 0.312, wrong direction b = -0.225, p = 0.060",
  "Q3 MPA",       "Not supported",  "DAICc = 3.01, weight = 0.182, low ns, medium ns",
  "Sensitivity a","Consistent",     "market b = -0.086 ns, pop b = -0.067 ns — direction consistent",
  "Sensitivity b","Consistent",     "model ordering identical, rugosity and pressure stable"
)

cat("\n--- Results summary ---\n")
print(results_summary)


# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()

