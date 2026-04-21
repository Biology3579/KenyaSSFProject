# ============================================================
#  00_DATA_PREPARATION.R
#
#  Shared data preparation for all biomass analyses.
#  Source at the top of each analysis script.
#
#  Objects created:
#    fish_data            — raw transect fish data
#    location_data        — site-level location and MPA data
#    scaled_predictors    — all scaled predictors (site-level)
#    final_predictors     — primary + sensitivity predictors
#    total_transects      — aggregated transect data (all species)
#    transect_model_data  — transect-level dataset with predictors
#    total_model_data     — site-level dataset with predictors
#
#  Usage:
#    source(here::here("00_data_preparation.R"))
#
#  Functional group aggregation handled within each
#  functional group script — not here.
# ============================================================

options(scipen = 999)

# ── PACKAGES ──────────────────────────────────────────────────
library(tidyverse)
library(sf)
library(glmmTMB)     # functional group analyses
library(DHARMa)      # functional group diagnostics
library(MuMIn)       # AICc throughout
library(gridExtra)   # plot arrangement
library(here)        # file paths
library(spdep)       # spatial autocorrelation
library(lme4)        # RE structure tests + transect sensitivity
library(car)         # VIFs


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
# Observed data points overlaid for transparency
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
  # For glmmTMB models — set site to first level to allow
  # population-level prediction with re.form = NA
  if (inherits(model, "glmmTMB") && "site" %in% names(data)) {
    grid$site <- levels(data$site)[1]
  }
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
                   y = if ("log_mean_biomass" %in% names(data))
                     log_mean_biomass else mean_biomass),
               colour = "grey40", size = 1.5,
               alpha  = 0.5, inherit.aes = FALSE) +
    labs(x = x_label, y = y_label) +
    theme_bw(base_size = 13) +
    theme(axis.title       = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}


# ============================================================
#  DATA LOADING
# ============================================================

fish_data     <- readr::read_rds(here("processed_data",
                                      "clean_fish_connectivity.rds"))
gravity_data  <- readr::read_rds(here("city_data",
                                      "locations_with_grav_combined.rds"))
chla_data     <- read.csv(here("processed_data",
                               "locations_with_chla_2009.csv"))
rugosity_data <- readr::read_rds(here("processed_data",
                                      "clean_dive_details_connectivity.rds"))
location_data <- readr::read_rds(here("processed_data",
                                      "clean_location_connectivity.rds"))

cat("Data loaded successfully.\n")


# ============================================================
#  TRANSECT AGGREGATION (ALL SPECIES — TOTAL BIOMASS)
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

cat("Total transects:", nrow(total_transects), "\n")
cat("Sites:",           n_distinct(total_transects$site), "\n")
cat("Countries:",       n_distinct(total_transects$country), "\n")


# ============================================================
#  PREDICTOR PREPARATION
# ============================================================

# ── Human gravity metrics ─────────────────────────────────────
# Three metrics retained — primary selected per analysis in Q1
gravity_sites <- gravity_data %>%
  st_drop_geometry() %>%
  group_by(site) %>%
  summarise(
    market_gravity  = mean(market_grav,        na.rm = TRUE),
    settlement_pop  = mean(settlement_tot_pop, na.rm = TRUE),
    settlement_grav = mean(nearest_pop75_grav, na.rm = TRUE),
    .groups = "drop"
  )

# ── Chlorophyll-a ─────────────────────────────────────────────
# Baseline covariate: background productivity
# (Warmuth et al. 2024, Samoilys et al. 2025)
chla_sites <- chla_data %>%
  group_by(site) %>%
  summarise(mean_annual_chla = mean(chla_annual_mean,
                                    na.rm = TRUE),
            .groups = "drop")

# ── Rugosity ──────────────────────────────────────────────────
# Baseline covariate: habitat structural complexity
# (Samoilys et al. 2025, Darling et al. 2017)
rugosity_sites <- rugosity_data %>%
  group_by(site) %>%
  summarise(rugosity = mean(rugosity, na.rm = TRUE),
            .groups = "drop")

# ── MPA status and connectivity ───────────────────────────────
# MPA: governance modifier — Q2 and Q3 only
# Connectivity: candidate driver — Q2 and Q3
location_sites <- location_data %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(
    mpa_status   = first(mpa_status),
    connectivity = mean(prop_connectivity, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
#  TRANSFORMATIONS AND SCALING
# ============================================================

raw_predictors <- location_sites %>%
  left_join(chla_sites,     by = "site") %>%
  left_join(rugosity_sites, by = "site") %>%
  left_join(gravity_sites,  by = "site")

# ── Transformations ───────────────────────────────────────────
# Rugosity:     no transformation (approximately normal)
# Gravity:      log (right-skewed)
# Chla:         log (right-skewed)
# Connectivity: no transformation
# MPA status:   unordered factor, reference = "none"

transformed_predictors <- raw_predictors %>%
  transmute(
    site                = site,
    rugosity            = rugosity,
    log_settlement_grav = log(settlement_grav),
    log_settlement_pop  = log(settlement_pop),
    log_market_gravity  = log(market_gravity),
    log_chla            = log(mean_annual_chla),
    connectivity        = connectivity,
    mpa_status          = factor(mpa_status,
                                 levels  = c("none", "low",
                                             "medium"),
                                 ordered = FALSE)
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
    log_chla_sc            = as.numeric(scale(log_chla)),
    connectivity_sc        = as.numeric(scale(connectivity)),
    mpa_status             = mpa_status
  )

# ── Final predictor set ───────────────────────────────────────
# Used by all functional group analyses.
# Primary predictors + sensitivity alternatives.
# DHW excluded — study period (2009–2016) predates the
# 2016 El Niño bleaching event; insufficient thermal
# stress variation to justify inclusion.

final_predictors <- scaled_predictors %>%
  dplyr::select(
    site,
    rugosity_sc,                # BASELINE    — habitat complexity
    log_settlement_grav_sc,     # Q1/Q2/Q3   — primary pressure (total biomass)
    log_chla_sc,                # BASELINE    — productivity
    connectivity_sc,            # Q2/Q3      — larval connectivity
    mpa_status,                 # Q2/Q3      — governance modifier
    log_settlement_pop_sc,      # SENSITIVITY
    log_market_gravity_sc       # SENSITIVITY
  )

cat("Predictor preparation complete.\n")
cat("Sites with predictors:", nrow(final_predictors), "\n")


# ============================================================
#  ANALYSIS DATASETS (TOTAL BIOMASS)
# ============================================================

# ── Transect-level dataset ────────────────────────────────────
transect_model_data <- total_transects %>%
  left_join(final_predictors, by = "site") %>%
  mutate(log_transect_biomass = log(transect_total_biomass))

cat("\nTransect model data:", nrow(transect_model_data), "rows,",
    n_distinct(transect_model_data$site), "sites\n")
cat("Biomass zeros:",
    sum(transect_model_data$transect_total_biomass == 0), "\n")

# ── Site-level dataset ────────────────────────────────────────
total_model_data <- transect_model_data %>%
  group_by(site, country) %>%
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
    country = as.factor(country)
  )

cat("\nSite model data:", nrow(total_model_data), "sites,",
    n_distinct(total_model_data$country), "countries\n")

# ── NA check ──────────────────────────────────────────────────
na_check <- total_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                log_chla_sc, connectivity_sc, mpa_status) %>%
  filter(if_any(everything(), is.na))

if (nrow(na_check) > 0) {
  warning("NAs found in primary predictors:")
  print(na_check)
} else {
  cat("NA check passed — no missing values in primary predictors.\n")
}

cat("\n--- Data preparation complete ---\n")
cat("Objects available:\n")
cat("  fish_data, gravity_data, chla_data,\n")
cat("  rugosity_data, location_data\n")
cat("  scaled_predictors, final_predictors\n")
cat("  total_transects, transect_model_data,\n")
cat("  total_model_data\n")
cat("  make_aicc_df(), plot_effect()\n")