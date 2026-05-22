# ============================================================
#  Predictor variable preparation
#
#  Shared data preparation and predictor diagnostics for
#  all biomass analyses. Sourced at the top of each fish group script.
#   source(here::here("predictor_preparation.R"))
# 
# ============================================================

options(scipen = 999)

# ── Packages ──────────────────────────────────────────────────
library(tidyverse)
library(sf)
library(glmmTMB)     # functional group analyses
library(DHARMa)      # functional group diagnostics
library(MuMIn)       # AICc throughout
library(gridExtra)   # plot arrangement
library(grid)        # textGrob for figure labels
library(here)        # file paths
library(spdep)       # spatial autocorrelation
library(lme4)        # RE structure tests + transect sensitivity
library(car)         # VIFs


# ── Functions ─────────────────────────────────────────────────

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

# Marginal effect plot for a single focal predictor.
# All other scaled predictors held at 0 (their mean).
# Observed data points overlaid for transparency.
# Works for lm and glmmTMB models.
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
  # glmmTMB: set site to first level for population-level prediction
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
  
  # Observed y: use log_mean_biomass if available, else mean_biomass
  obs_y <- if ("log_mean_biomass" %in% names(data)) {
    data$log_mean_biomass
  } else {
    data$mean_biomass
  }
  
  ggplot(grid, aes(x = .data[[focal_var]])) +
    geom_ribbon(aes(ymin = lwr, ymax = upr),
                fill = colour, alpha = 0.15) +
    geom_line(aes(y = fit), colour = colour, linewidth = 1.1) +
    geom_point(data = data.frame(data, obs_y = obs_y),
               aes(x = .data[[focal_var]], y = obs_y),
               colour = "grey40", size = 1.5,
               alpha  = 0.5, inherit.aes = FALSE) +
    labs(x = x_label, y = y_label) +
    theme_bw(base_size = 13) +
    theme(axis.title       = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}


# ============================================================
#  1. DATA LOADING
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


# ============================================================
#  2. TRANSECT AGGREGATION (ALL SPECIES — TOTAL BIOMASS)
#
#  Minimum 3 transects per site retained.
#  Functional group scripts filter fish_data to relevant
#  taxa and repeat this aggregation step independently.
# ============================================================

total_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_total_biomass = sum(tot_wt_g, na.rm = TRUE),
    transect_total_count   = sum(number,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(site = as.factor(site))

cat("Total transects:", nrow(total_transects), "\n")
cat("Sites:",           n_distinct(total_transects$site), "\n")


# ============================================================
#  3. PREDICTOR PREPARATION
# ============================================================

# ── Human gravity metrics ─────────────────────────────────────
# Three metrics retained — primary metric selected per
# response variable in Q1 via AICc comparison.
# Settlement gravity: weights residential populations by
#   travel-cost distance — captures SSF exploitation footprint.
# Market gravity: weights urban centres by accessibility —
#   validated globally (Cinner et al. 2016).
# Settlement population: unweighted sum within 25 km radius —
#   distance-independent measure of human proximity.
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
# Baseline covariate: background productivity.
# Mean annual Aqua-MODIS composite for calendar year
# preceding each survey.
chla_sites <- chla_data %>%
  group_by(site) %>%
  summarise(mean_annual_chla = mean(chla_annual_mean, na.rm = TRUE),
            .groups = "drop")

# ── Rugosity ──────────────────────────────────────────────────
# Baseline covariate: habitat structural complexity.
# Visual scale 0–5 (Wilson et al. 2007).
rugosity_sites <- rugosity_data %>%
  group_by(site) %>%
  summarise(rugosity = mean(rugosity, na.rm = TRUE),
            .groups = "drop")

# ── MPA status and connectivity ───────────────────────────────
# MPA: three-level governance modifier (none / low / medium),
#   reflecting realised protection outcomes (Osuka et al. 2021).
# Connectivity: proportional oceanographic connectivity
#   (Warmuth et al. 2024) — fraction of WIO reef polygons
#   connected to each site within a 30-day dispersal window.
location_sites <- location_data %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(
    mpa_status   = first(mpa_status),
    connectivity = mean(prop_connectivity, na.rm = TRUE),
    ecoregion    = first(ecoregion),
    .groups = "drop"
  )

# ── Join all predictors ───────────────────────────────────────
raw_predictors <- location_sites %>%
  left_join(chla_sites,     by = "site") %>%
  left_join(rugosity_sites, by = "site") %>%
  left_join(gravity_sites,  by = "site")


# ============================================================
#  4. PREDICTOR DISTRIBUTION CHECKS AND TRANSFORMATIONS
# ============================================================

predictor_labels <- c("Rugosity", "Chlorophyll-a",
                      "Settlement gravity", "Market gravity",
                      "Settlement pop.", "Connectivity")

predictor_order_raw <- c("rugosity", "mean_annual_chla",
                         "settlement_grav", "market_gravity",
                         "settlement_pop", "connectivity")

# ── Raw distributions ─────────────────────────────────────────
( p_pred_raw <- raw_predictors %>%
    dplyr::select(all_of(predictor_order_raw)) %>%
    pivot_longer(everything(),
                 names_to  = "variable",
                 values_to = "value") %>%
    mutate(variable = factor(variable,
                             predictor_order_raw,
                             predictor_labels)) %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
    facet_wrap(~ variable, scales = "free") +
    labs(title = "Raw predictors") +
    theme_bw() )

# ── Transformations ───────────────────────────────────────────
# Rugosity:     no transformation (approximately normal)
# Gravity:      log (right-skewed)
# Chla:         log (right-skewed)
# Connectivity: no transformation (approximately uniform)
# MPA status:   unordered factor, reference level = "none"

transformed_predictors <- raw_predictors %>%
  transmute(
    site                = site,
    ecoregion           = ecoregion,
    rugosity            = rugosity,
    log_settlement_grav = log(settlement_grav),
    log_settlement_pop  = log(settlement_pop),
    log_market_gravity  = log(market_gravity),
    log_chla            = log(mean_annual_chla),
    connectivity        = connectivity,
    mpa_status          = factor(mpa_status,
                                 levels  = c("none", "low", "medium"),
                                 ordered = FALSE)
  )

predictor_order_tran <- c("rugosity", "log_chla",
                          "log_settlement_grav",
                          "log_market_gravity",
                          "log_settlement_pop",
                          "connectivity")

# ── Transformed distributions ─────────────────────────────────
( p_pred_tran <- transformed_predictors %>%
    dplyr::select(all_of(predictor_order_tran)) %>%
    pivot_longer(everything(),
                 names_to  = "variable",
                 values_to = "value") %>%
    mutate(variable = factor(variable,
                             predictor_order_tran,
                             predictor_labels)) %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
    facet_wrap(~ variable, scales = "free") +
    labs(title = "Transformed predictors") +
    theme_bw() )
# Gravity metrics and chla: right skew corrected by log transformation.
# Rugosity and connectivity: distributions acceptable without transformation.


# ============================================================
#  5. STANDARDISATION
#
#  z-score scaling (mean = 0, SD = 1) applied to all
#  continuous predictors. Enables direct comparison of
#  effect sizes across predictors with different units.
#  Scaling applied once here — consistent across all analyses.
# ============================================================

scaled_predictors <- transformed_predictors %>%
  transmute(
    site                   = site,
    ecoregion              = ecoregion,
    rugosity_sc            = as.numeric(scale(rugosity)),
    log_settlement_grav_sc = as.numeric(scale(log_settlement_grav)),
    log_settlement_pop_sc  = as.numeric(scale(log_settlement_pop)),
    log_market_gravity_sc  = as.numeric(scale(log_market_gravity)),
    log_chla_sc            = as.numeric(scale(log_chla)),
    connectivity_sc        = as.numeric(scale(connectivity)),
    mpa_status             = mpa_status
  )

# ── Final predictor set ───────────────────────────────────────
# Joined to response data in each analysis script via site.

final_predictors <- scaled_predictors %>%
  dplyr::select(
    site,
    ecoregion,
    rugosity_sc,                # BASELINE    — habitat complexity
    log_chla_sc,                # BASELINE    — productivity
    log_settlement_grav_sc,     # Q1/Q2/Q3   — primary pressure metric
    connectivity_sc,            # Q2/Q3      — larval connectivity
    mpa_status,                 # Q2/Q3      — governance modifier
    log_settlement_pop_sc,      # SENSITIVITY
    log_market_gravity_sc       # SENSITIVITY
  )

# ============================================================
#  6. COLLINEARITY CHECKS
#
#  Pairwise correlation matrix
#  Variance inflation factors (VIFs)
#
#  VIFs are a property of the predictor set, not the response
#  variable or model family. Computed once here on lm() —
#  applies equally to Tweedie GLMs since collinearity is
#  determined by the design matrix alone.
#
#  Gravity metrics are never modelled simultaneously
#  (single-metric-per-model rule), so VIFs are checked
#  separately for each at the baseline + pressure stage.
#  All subsequent stages use settlement gravity as the
#  selected metric (established in Q1 of total biomass and
#  confirmed across functional groups).
#
#  Log total biomass response is used for
#  VIF computation — the response does not affect VIF values.
# ============================================================

# ── Pairwise correlation matrix ──────────────────────────

corr_matrix <- scaled_predictors %>%
  mutate(mpa_numeric = as.numeric(mpa_status)) %>%
  dplyr::select(ends_with("_sc"), mpa_numeric) %>%
  rename(
    "Market gravity"     = log_market_gravity_sc,
    "Settlement gravity" = log_settlement_grav_sc,
    "Settlement pop."    = log_settlement_pop_sc,
    "Chlorophyll-a"      = log_chla_sc,
    "Rugosity"           = rugosity_sc,
    "Connectivity"       = connectivity_sc,
    "MPA status"         = mpa_numeric
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
    low    = "#2c7bb6",
    mid    = "white",
    high   = "#d7191c",
    limits = c(-1, 1),
    name   = "Correlation"
  ) +
  scale_x_discrete(limits = vars, drop = FALSE,
                   position = "bottom") +
  scale_y_discrete(limits = rev(vars), drop = FALSE,
                   position = "left") +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid           = element_blank(),
    axis.title           = element_blank(),
    axis.text.x          = element_text(angle = 45, hjust = 1,
                                        vjust = 1, size = 10),
    axis.text.y          = element_text(hjust = 1, size = 9.5),
    legend.position      = c(1.015, 0.48),
    legend.key.height    = unit(1.75, "cm"),
    legend.key.width     = unit(0.50, "cm"),
    legend.title         = element_text(size = 9.75,
                                        colour = "grey40"),
    legend.text          = element_text(size = 9,
                                        colour = "grey40"),
    plot.margin          = margin(10, 10, 10, 10)
  )

# Collinearity notes:
#
# Gravity metrics (never modelled together):
#   Settlement grav vs settlement pop:  r =  0.66
#   Settlement grav vs market gravity:  r =  0.54
#   Settlement pop  vs market gravity:  r =  0.65
#   Moderate positive — expected, all proxy human pressure.
#   Single-metric-per-model rule removes any concern.
#
# Chlorophyll-a vs settlement gravity: r = -0.57
#   Strongest pairwise correlation among predictors that
#   appear in the same model. Productive sites tend to be
#   more remote and less fished. Monitor chla coefficient
#   in models including settlement gravity simultaneously.
#
# All other pairwise |r| < 0.50 — no blocking collinearity.


# ── Variance inflation factors ───────────────────────────

vif_data <- total_transects %>%
  group_by(site) %>%
  summarise(mean_biomass = mean(transect_total_biomass,
                                na.rm = TRUE),
            .groups = "drop") %>%
  mutate(log_mean_biomass = log(mean_biomass),
         site = as.character(site)) %>%
  left_join(final_predictors, by = "site")

cat("\n--- VIF: baseline ---\n")
vif(lm(log_mean_biomass ~ rugosity_sc + log_chla_sc,
       data = vif_data))

cat("\n--- VIF: baseline + settlement gravity ---\n")
vif(lm(log_mean_biomass ~ rugosity_sc + log_chla_sc +
         log_settlement_grav_sc,
       data = vif_data))

cat("\n--- VIF: baseline + market gravity ---\n")
vif(lm(log_mean_biomass ~ rugosity_sc + log_chla_sc +
         log_market_gravity_sc,
       data = vif_data))

cat("\n--- VIF: baseline + settlement population ---\n")
vif(lm(log_mean_biomass ~ rugosity_sc + log_chla_sc +
         log_settlement_pop_sc,
       data = vif_data))

cat("\n--- VIF: baseline + settlement gravity + connectivity ---\n")
vif(lm(log_mean_biomass ~ rugosity_sc + log_chla_sc +
         log_settlement_grav_sc + connectivity_sc,
       data = vif_data))

cat("\n--- VIF: baseline + settlement gravity + MPA ---\n")
vif(lm(log_mean_biomass ~ rugosity_sc + log_chla_sc +
         log_settlement_grav_sc + mpa_status,
       data = vif_data))

# Results:
#  Baseline (rugosity + chla):
#    rugosity VIF = 1.000, chla VIF = 1.000
#    Baseline predictors completely orthogonal.
#
#  Baseline + settlement gravity:
#    rugosity VIF = 1.011, chla VIF = 1.427,
#    settlement gravity VIF = 1.438
#
#  Baseline + market gravity:
#    rugosity VIF = 1.050, chla VIF = 1.111,
#    market gravity VIF = 1.161
#
#  Baseline + settlement population:
#    rugosity VIF = 1.036, chla VIF = 1.031,
#    settlement population VIF = 1.067
#
#  Baseline + settlement gravity + connectivity:
#    rugosity VIF = 1.012, chla VIF = 1.486,
#    settlement gravity VIF = 1.446,
#    connectivity VIF = 1.090
#
#  Baseline + settlement gravity + MPA:
#    rugosity GVIF^(1/(2*Df)) = 1.011 — no concern
#    chla GVIF^(1/(2*Df)) = 1.259 — no concern
#    settlement gravity GVIF^(1/(2*Df)) = 1.377 — no concern
#    mpa_status GVIF^(1/(2*Df)) = 1.073 — no concern
#    (GVIF reported for factor predictors — equivalent
#    to VIF for continuous predictors when Df = 1)
#
#  All VIFs/GVIFs < 2.0 at every stage and across all
#  three pressure metrics — no multicollinearity concern
#  throughout Q1-Q3 model sequence for any response variable.


# ============================================================
#  7. SST SENSITIVITY CHECK
#
#  SST evaluated as a potential additional baseline covariate.
#  Excluded from primary analyses due to strong collinearity
#  with larval connectivity (r = 0.70) and negligible
#  improvement in model fit (DAICc = 1.82) — see total
#  biomass script for details.
#
#  Assessed here separately from the primary correlation
#  matrix to keep the predictor diagnostics clean.
#  SST is not carried forward into any analysis script.
# ============================================================

# ── Load SST ─────────────────────────────────────────────────
sst_data <- read.csv(here("processed_data",
                          "locations_with_sst_2009.csv"))

sst_sites <- sst_data %>%
  group_by(site) %>%
  summarise(mean_annual_sst = mean(sst_annual_mean, na.rm = TRUE),
            .groups = "drop")

# ── Distribution check ────────────────────────────────────────
hist(sst_sites$mean_annual_sst,
     breaks = 20, col = "#2c7bb6", border = "white",
     main = "SST: raw", xlab = "Mean annual SST (°C)")

cat("\n--- Shapiro-Wilk: raw SST ---\n")
shapiro.test(sst_sites$mean_annual_sst)

# Distribution is left-skewed/bimodal — likely reflects
# geographic structure across the WIO (lower cluster ~27.0–27.5°C
# corresponds to higher-latitude sites in Madagascar/Mozambique;
# upper cluster ~27.8–28.2°C to equatorial sites in
# Comoros/Tanzania north). Log transformation not appropriate
# as it would increase rather than reduce skew.
# Raw SST used for correlation check.

sst_sites <- sst_sites %>%
  mutate(sst_sc = as.numeric(scale(mean_annual_sst)))

# ── SST correlation matrix ────────────────────────────────────
sst_corr_data <- scaled_predictors %>%
  left_join(sst_sites %>% dplyr::select(site, sst_sc),
            by = "site") %>%
  mutate(mpa_numeric = as.numeric(mpa_status)) %>%
  dplyr::select(ends_with("_sc"), mpa_numeric) %>%
  rename(
    "Market gravity"     = log_market_gravity_sc,
    "Settlement gravity" = log_settlement_grav_sc,
    "Settlement pop."    = log_settlement_pop_sc,
    "Chlorophyll-a"      = log_chla_sc,
    "Rugosity"           = rugosity_sc,
    "Connectivity"       = connectivity_sc,
    "MPA status"         = mpa_numeric,
    "SST"                = sst_sc
  )

sst_corr_matrix <- cor(sst_corr_data, use = "complete.obs")
sst_vars        <- colnames(sst_corr_matrix)

sst_corr_long <- as.data.frame(as.table(sst_corr_matrix)) %>%
  rename(x = Var1, y = Var2, corr = Freq) %>%
  mutate(
    x = factor(x, levels = sst_vars),
    y = factor(y, levels = rev(sst_vars))
  ) %>%
  filter(as.numeric(x) < length(sst_vars) - as.numeric(y) + 1)

ggplot(sst_corr_long, aes(x = x, y = y, fill = corr)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", corr)),
            size = 3.2, color = "black") +
  scale_fill_gradient2(
    low    = "#2c7bb6",
    mid    = "white",
    high   = "#d7191c",
    limits = c(-1, 1),
    name   = "Correlation"
  ) +
  scale_x_discrete(limits = sst_vars, drop = FALSE,
                   position = "bottom") +
  scale_y_discrete(limits = rev(sst_vars), drop = FALSE,
                   position = "left") +
  coord_equal() +
  labs(title = "Sensitivity check: SST included") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid        = element_blank(),
    axis.title        = element_blank(),
    axis.text.x       = element_text(angle = 45, hjust = 1,
                                     vjust = 1, size = 10),
    axis.text.y       = element_text(hjust = 1, size = 9.5),
    legend.position   = c(1.015, 0.48),
    legend.key.height = unit(1.75, "cm"),
    legend.key.width  = unit(0.50, "cm"),
    legend.title      = element_text(size = 9.75,
                                     colour = "grey40"),
    legend.text       = element_text(size = 9,
                                     colour = "grey40"),
    plot.margin       = margin(10, 10, 10, 10)
  )

# SST collinearity notes:
#   SST vs connectivity:       r =  0.70
#     Strong collinearity — primary reason SST excluded.
#     SST and connectivity are likely capturing overlapping
#     spatial gradients across the WIO (equatorial vs
#     higher-latitude sites differ in both thermal regime
#     and larval connectivity patterns).
#   SST vs settlement gravity: r =  0.26 — no concern
#   SST vs market gravity:     r =  0.14 — no concern
#   SST vs settlement pop.:    r =  0.14 — no concern
#   SST vs chlorophyll-a:      r = -0.14 — no concern
#   SST vs rugosity:           r = -0.05 — no concern
#   SST vs MPA status:         r =  0.19 — no concern
#
# Conclusion: SST excluded from all analyses due to
# collinearity with connectivity (r = 0.70). Including
# both would confound interpretation of connectivity
# effects. Negligible improvement in model fit when added
# to baseline (DAICc = 1.82 — see total_biomass.R).
# Reported in methods as justification for exclusion.

# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()

# end of script
