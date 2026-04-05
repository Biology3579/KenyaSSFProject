# ============================================================
#  LARGE EXCAVATOR BIOMASS
#  Response:   Large excavator biomass (continuous, zero-inflated)
#  Predictors: SST, Chl-a, Human gravity, DHW, Connectivity, Rugosity
#  Random fx:  Transect nested in site  (1 | site/transect)
#  Families:   Tweedie, Zero-inflated Gamma/Lognormal (glmmTMB)
# ============================================================

options(scipen = 999) #turn off scientific notation

# ── PACKAGES ─────────────────────────────────────────────────
library(tidyverse)
library(sf)
library(glmmTMB)       # Tweedie, ZI models
library(DHARMa)        # Residual diagnostics
library(MuMIn)         # AICc
library(AICcmodavg)    # aictab()
library(ggcorrplot)    # Correlation heatmaps
library(corrplot)      # Alternative corrplot
library(MASS)          # boxcox
library(here)          # Reproducible paths

# ── PREPARE LARGE EXCAVATOR MODEL DATA ─────────────────────────────────────────────────
# Replace these paths with your own files.
# Expected: one row per transect, columns for biomass per group + all predictors.

fish_2009 <- readr::read_rds(here::here("processed_data", "clean_fish_2009.rds"))
location_2009 <- readr::read_rds(here::here("processed_data", "clean_location_2009.rds"))
gravity_2009 <- readr::read_rds(here::here("city_data", "locations_with_grav_combined.rds"))
old_grav_2009 <- read.csv(here::here("processed_data", "locations_with_gravity.csv"))
chla_2009 <- read.csv(here::here("processed_data", "locations_with_chla_2009.csv"))
sst_2009 <- read.csv(here::here("processed_data", "locations_with_sst_2009.csv"))
rugosity_2009 <- readr::read_rds(here::here("processed_data", "clean_dive_details_2009.rds"))

# ============================================================
#  LARGE EXCAVATOR BIOMASS PER TRANSECT
# ============================================================

excavator_transects <- fish_2009 %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_excavator_biomass = sum(
      ifelse(trophic_group %in% c("large_excavators"), tot_wt_g, 0),
      na.rm = TRUE
    ),
    country = first(country),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(site = as.factor(site))

cat("Transects after filtering:", nrow(excavator_transects), "\n")
cat("Sites after filtering:",     n_distinct(excavator_transects$site), "\n")


# ============================================================
#  LARGE EXCAVATOR BIOMASS EXPLORATION
# ============================================================

# ── Basic summary ─────────────────────────────────────────────
summary(excavator_transects$transect_excavator_biomass)

zeros <- mean(excavator_transects$transect_excavator_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")
# Check proportion — strong zero inflation may favour Tweedie or zero-inflated family

# ── Raw distribution ──────────────────────────────────────────
( excavator_raw <- ggplot(excavator_transects, aes(x = transect_excavator_biomass)) +
    geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
    labs(x = "Large excavator biomass per transect (g)", y = "Frequency",
         title = "Raw Large Excavator Biomass") +
    theme_bw() )

# ── Log transformation ────────────────────────────────────────
# Small constant (+0.01) handles zeros; check sensitivity to constant choice
excavator_transects <- excavator_transects %>%
  mutate(
    log_excavator_biomass  = log(transect_excavator_biomass + 0.01),
    sqrt_excavator_biomass = sqrt(transect_excavator_biomass)
  )

( excavator_log <- ggplot(excavator_transects, aes(x = log_excavator_biomass)) +
    geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
    labs(x = "log(biomass + 0.01)", y = "Frequency",
         title = "Log-transformed Large Excavator Biomass") +
    theme_bw() )

( excavator_sqrt <- ggplot(excavator_transects, aes(x = sqrt_excavator_biomass)) +
    geom_histogram(bins = 30, fill = "#d7191c", colour = "white") +
    labs(x = "sqrt(biomass)", y = "Frequency",
         title = "Sqrt-transformed Large Excavator Biomass") +
    theme_bw() )

jpeg("excavator_biomass_distributions.jpg", width = 33, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(excavator_raw, excavator_log, excavator_sqrt, ncol = 3)
dev.off()

# ── Box-Cox: what power transformation does the data suggest? ─
# Fit a simple lm on non-zero values only (boxcox requires y > 0)
excavator_nonzero <- excavator_transects %>% filter(transect_excavator_biomass > 0)

MASS::boxcox(
  lm(transect_excavator_biomass ~ 1, data = excavator_nonzero),
  lambda = seq(-2, 2, 0.1)
)
# lambda ~ 0  → log transformation appropriate

# ── Variation by site ─────────────────────────────────────────
ggplot(excavator_transects, aes(x = reorder(site, transect_excavator_biomass, median),
                                y = transect_excavator_biomass)) +
  geom_boxplot(fill = "#2c7bb6", alpha = 0.6,
               outlier.colour = "black",
               outlier.size = 1) +
  coord_flip() +
  labs(x = NULL, y = "Large excavator biomass (g)",
       title = "Biomass distribution by site") +
  theme_bw(base_size = 9)

# ── Overall zero inflation check ──────────────────────────────────────
zeros <- mean(excavator_transects$transect_excavator_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")

# ── Zeros by site ─────────────────────────────────────────────
excavator_transects %>%
  group_by(site) %>%
  summarise(
    n_transects  = n(),
    prop_zeros   = mean(transect_excavator_biomass == 0),
    mean_biomass = mean(transect_excavator_biomass),
    .groups = "drop"
  ) %>%
  arrange(desc(prop_zeros)) %>%
  print(n = Inf)

# ============================================================
#  PREDICTOR PREPARATION
# ============================================================

#  AGGREGATE ALL PREDICTORS TO SITE LEVEL --------------------

# ── Human gravity (new metrics) ──────────────────────────
gravity_sites <- gravity_2009 %>%
  st_drop_geometry() %>%
  group_by(site) %>%
  summarise(
    market_gravity  = mean(market_grav,           na.rm = TRUE),
    settlement_pop  = mean(settlement_tot_pop,    na.rm = TRUE),
    settlement_grav = mean(nearest_pop75_grav,    na.rm = TRUE),
    .groups = "drop"
  )

# ── Chlorophyll-a ─────────────────────────────────────────
chla_sites <- chla_2009 %>%
  group_by(site) %>%
  summarise(
    mean_annual_chla = mean(chla_annual_mean, na.rm = TRUE),
    .groups = "drop"
  )

# ── Sea surface temperature ───────────────────────────────
sst_sites <- sst_2009 %>%
  filter(!is.na(sst_annual_mean)) %>%
  group_by(site) %>%
  summarise(mean_annual_sst = mean(sst_annual_mean, na.rm = TRUE),
            .groups = "drop")

# ── Rugosity ──────────────────────────────────────────────
rugosity_sites <- rugosity_2009 %>%
  group_by(site) %>%
  summarise(
    rugosity = mean(rugosity, na.rm = TRUE),
    .groups = "drop"
  )

#  JOIN EVERYTHING INTO ONE DATAFRAME -----------------------

excavator_model_data <- excavator_transects %>%
  left_join(gravity_sites,  by = "site") %>%
  left_join(chla_sites,     by = "site") %>%
  left_join(sst_sites,      by = "site") %>%
  left_join(rugosity_sites, by = "site")

# Shared predictor order — raw names first, used to align both plots
predictor_order_raw  <- c("market_gravity", "settlement_grav", "settlement_pop",
                          "mean_annual_chla", "mean_annual_sst", "rugosity")
predictor_order_tran <- c("log_market_gravity", "log_settlement_grav", "log_settlement_pop",
                          "log_chla", "mean_annual_sst", "rugosity")

# Shared display labels (same order as above)
predictor_labels <- c("Market gravity", "Settlement gravity", "Settlement pop.",
                      "Chlorophyll-a", "SST", "Rugosity")

excavator_transects %>%
  left_join(gravity_sites,  by = "site") %>%
  left_join(chla_sites,     by = "site") %>%
  left_join(sst_sites,      by = "site") %>%
  left_join(rugosity_sites, by = "site") %>%
  dplyr::select(all_of(predictor_order_raw)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = predictor_order_raw,
                           labels = predictor_labels)) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
  facet_wrap(~ variable, scales = "free") +
  labs(title = "Raw predictors") +
  theme_bw()

# LOG-TRANSFORM & SCALE PREDICTORS --------------------------

excavator_model_data <- excavator_model_data %>%
  mutate(
    # Log-transformations (right-skewed variables)
    log_market_gravity  = log(market_gravity),
    log_settlement_grav = log(settlement_grav),
    log_settlement_pop  = log(settlement_pop),
    log_chla            = log(mean_annual_chla),
    
    # Scaled versions (mean = 0, SD = 1) — use these in models
    log_market_gravity_sc  = as.numeric(scale(log_market_gravity)),
    log_settlement_grav_sc = as.numeric(scale(log_settlement_grav)),
    log_settlement_pop_sc  = as.numeric(scale(log_settlement_pop)),
    log_chla_sc            = as.numeric(scale(log_chla)),
    sst_sc                 = as.numeric(scale(mean_annual_sst)),
    rugosity_sc            = as.numeric(scale(rugosity)),
    
    # Log-transform response for diagnostics / gaussian alternative
    log_excavator_biomass = log(transect_excavator_biomass + 0.01)
  )

# Check transformations
excavator_model_data %>%
  dplyr::select(all_of(predictor_order_tran)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = predictor_order_tran,
                           labels = predictor_labels)) %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
  facet_wrap(~ variable, scales = "free") +
  labs(title = "Transformed predictors") +
  theme_bw()

# ============================================================
# PREDICTOR CORRELATION MATRIX
# ============================================================

predictor_vars <- excavator_model_data %>%
  dplyr::select(
    log_market_gravity_sc, log_settlement_grav_sc, log_settlement_pop_sc,
    log_chla_sc, sst_sc, rugosity_sc
  ) %>%
  rename(
    "Market gravity"     = log_market_gravity_sc,
    "Settlement gravity" = log_settlement_grav_sc,
    "Settlement pop."    = log_settlement_pop_sc,
    "Chlorophyll-a"      = log_chla_sc,
    "SST"                = sst_sc,
    "Rugosity"           = rugosity_sc
  )

corr_matrix <- cor(predictor_vars, use = "complete.obs")
p_mat       <- cor_pmat(corr_matrix)

corrplot(corr_matrix,
         method      = "square",
         type        = "lower",
         tl.col      = "black",
         tl.srt      = 0,
         tl.offset   = 0.5,
         addCoef.col = "black",
         number.cex  = 0.8,
         col         = colorRampPalette(c("#d73027", "white", "#4575b4"))(200),
         mar         = c(0, 0, 4, 2))

# ── Decision guide ────────────────────────────────────────────
# Need to choose between gravity metrics but I will keep the rest

# ============================================================
#  FAMILY SELECTION
#  Work from simplest to most complex.
#  Each step fits the full predictor set so family differences
#  are not confounded with predictor differences.
#  Diagnose with DHARMa before moving to the next family.
#
#  Ladder:
#    1. Gaussian on log(y + 0.01)   — simplest baseline
#    2. Tweedie                      — handles zeros + skew natively
#    3. Zero-inflated Tweedie        — if ZI test flags excess zeros
# ============================================================

# ── STEP F1: Gaussian on log-transformed response ─────────────
mF1_gaussian <- glmmTMB(
  log_excavator_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = gaussian(),
  data   = excavator_model_data
)

resF1 <- simulateResiduals(mF1_gaussian, n = 1000)

jpeg("dharma_F1_gaussian_excavators.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF1, main = "DHARMa — Gaussian on log(y + 0.01)")
dev.off()

plot(resF1)
testZeroInflation(resF1)
testDispersion(resF1)

# ── STEP F2: Tweedie ──────────────────────────────────────────
mF2_tweedie <- glmmTMB(
  transect_excavator_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = excavator_model_data
)

resF2 <- simulateResiduals(mF2_tweedie, n = 1000)

jpeg("dharma_F2_tweedie_excavators.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF2, main = "DHARMa — Tweedie")
dev.off()

plot(resF2)
testZeroInflation(resF2)
testDispersion(resF2)

# ── STEP F3: Zero-inflated Tweedie ────────────────────────────
# Adds an explicit Bernoulli component (ziformula = ~1) that
# models the probability of a *structural* zero — a site/transect
# where piscivores are absent regardless of environmental conditions
# (e.g. chronically overfished). 
# This is separate from the sampling zeros already captured by Tweedie's compound structure.
# Use only if F2 diagnostics indicate remaining zero inflation.

mF3_tweedie_zi <- glmmTMB(
  transect_excavator_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = pisc_model_data
)

resF3 <- simulateResiduals(mF3_tweedie_zi, n = 1000)

jpeg("dharma_F3_tweedie_zi.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF3, main = "DHARMa — Zero-inflated Tweedie")
dev.off()

plot(resF3)
testZeroInflation(resF3)
testDispersion(resF3)

# Check AICs
# ── AICc comparison: F2 vs F3 ────────────────────────────────
aicc_family <- MuMIn::AICc(mF2_tweedie, mF3_tweedie_zi)

aicc_family$delta  <- aicc_family$AICc - min(aicc_family$AICc)
aicc_family$weight <- exp(-0.5 * aicc_family$delta) /
  sum(exp(-0.5 * aicc_family$delta))
aicc_family$weight <- round(aicc_family$weight, 4)
aicc_family$AICc   <- round(aicc_family$AICc, 2)
aicc_family$delta  <- round(aicc_family$delta, 2)

print(aicc_family)
# ── FAMILY SELECTION ──────────────────────────────────
# Update this section after inspecting DHARMa diagnostics and AICc.
# Follow the same decision logic as for piscivores:
#   - Prefer the simpler family unless diagnostics clearly favour the more complex one.
#   - ZI Tweedie warranted only if ZI test is significant AND AICc improves by > 2.

# ============================================================
#  CANDIDATE MODEL FITTING
#  Response:   transect_excavator_biomass (zero-inflated, continuous)
#  Family:     Tweedie (link = "log") — confirm via family selection above
#  Random fx:  (1 | site)
#  Predictors: sst_sc, log_chla_sc, rugosity_sc,
#              log_market_gravity_sc OR log_settlement_grav_sc
#              (not both — r = 0.56, compared separately via AICc)
# ============================================================

m0 <- glmmTMB(transect_excavator_biomass ~ 1 + (1 | site),
              family = tweedie(link = "log"), data = excavator_model_data)

m_env <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc +
                   (1 | site),
                 family = tweedie(link = "log"), data = excavator_model_data)

m_market <- glmmTMB(transect_excavator_biomass ~ log_market_gravity_sc + (1 | site),
                    family = tweedie(link = "log"), data = excavator_model_data)

m_settgrav <- glmmTMB(transect_excavator_biomass ~ log_settlement_grav_sc + (1 | site),
                      family = tweedie(link = "log"), data = excavator_model_data)

m_hab <- glmmTMB(transect_excavator_biomass ~ rugosity_sc + (1 | site),
                 family = tweedie(link = "log"), data = excavator_model_data)

m_env_market <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc +
                          log_market_gravity_sc + (1 | site),
                        family = tweedie(link = "log"), data = excavator_model_data)

m_env_settgrav <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc +
                            log_settlement_grav_sc + (1 | site),
                          family = tweedie(link = "log"), data = excavator_model_data)

m_habitat_market <- glmmTMB(transect_excavator_biomass ~ rugosity_sc + log_market_gravity_sc + (1 | site),
                            family = tweedie(link = "log"), data = excavator_model_data)

m_habitat_settgrav <- glmmTMB(transect_excavator_biomass ~ rugosity_sc + log_settlement_grav_sc + (1 | site),
                              family = tweedie(link = "log"), data = excavator_model_data)

m_full_market <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc +
                           log_market_gravity_sc + rugosity_sc +
                           (1 | site),
                         family = tweedie(link = "log"), data = excavator_model_data)

m_full_settgrav <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc +
                             log_settlement_grav_sc + rugosity_sc +
                             (1 | site),
                           family = tweedie(link = "log"), data = excavator_model_data)

# ── AICc TABLE ─────────────────────────────────────────────

model_list <- list(
  "Null"                     = m0,
  "Environment"              = m_env,
  "Market gravity"           = m_market,
  "Settlement gravity"       = m_settgrav,
  "Habitat"                  = m_hab,
  "Env + market gravity"     = m_env_market,
  "Env + settlement gravity" = m_env_settgrav,
  "Habitat + market gravity" = m_habitat_market,
  "Habitat + settlement gravity" = m_habitat_settgrav,
  "Full (market gravity)"    = m_full_market,
  "Full (settlement gravity)"= m_full_settgrav
)

aicc_vals <- sapply(model_list, AICc)
delta     <- aicc_vals - min(aicc_vals)
weights   <- exp(-0.5 * delta) / sum(exp(-0.5 * delta))

aicc_df <- data.frame(
  Model  = names(model_list),
  AICc   = round(aicc_vals, 2),
  Delta  = round(delta,     2),
  Weight = round(weights,   4),
  row.names = NULL
) %>% arrange(AICc)

print(aicc_df)

# ── DHARMA RESIDUAL DIAGNOSTICS ON BEST MODEL ──────────────

res <- simulateResiduals(m_hab, n = 1000)

jpeg("dharma_excavator_full_market_model.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res, main = "DHARMa diagnostics — large excavator best model")
dev.off()

testZeroInflation(res)
testDispersion(res)
testOutliers(res)
plot(res)

plotResiduals(res, excavator_model_data$rugosity_sc,            xlab = "Rugosity")

# ── SUMMARY ────────────────────────────────────────

summary(m_hab)

# ============================================================
#  MARGINAL EFFECT PLOTS
# ============================================================

plot_effect <- function(model, data, focal_var, x_label,
                        y_label = "Large excavator biomass (g)",
                        n = 200) {
  
  grid <- data %>%
    slice(1) %>%
    dplyr::select(all_of(names(data)[endsWith(names(data), "_sc")])) %>%
    mutate(across(everything(), ~ 0)) %>%
    slice(rep(1, n)) %>%
    mutate(
      !!sym(focal_var) := seq(
        min(data[[focal_var]], na.rm = TRUE),
        max(data[[focal_var]], na.rm = TRUE),
        length.out = n
      ),
      site = levels(data$site)[1]
    )
  
  pred     <- predict(model, newdata = grid, type = "response",
                      se.fit = TRUE, re.form = NA)
  grid$fit <- pred$fit
  grid$lwr <- pmax(pred$fit - 1.96 * pred$se.fit, 0)
  grid$upr <- pred$fit + 1.96 * pred$se.fit
  
  ggplot(grid, aes(x = !!sym(focal_var))) +
    geom_ribbon(aes(ymin = lwr, ymax = upr),
                fill = "#2c7bb6", alpha = 0.2) +
    geom_line(aes(y = fit), colour = "#2c7bb6", linewidth = 1.1) +
    labs(x = x_label, y = y_label) +
    theme_bw(base_size = 13) +
    theme(axis.title = element_text(face = "bold"))
}

p_rug  <- plot_effect(m_full_market, excavator_model_data,
                      "rugosity_sc",           "Rugosity (scaled)")

# ============================================================
#  COUNTRY EFFECTS
# ============================================================

mF2_tweedie_nested <- glmmTMB(
  transect_excavator_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | country/site),
  family = tweedie(link = "log"),
  data   = excavator_model_data
)

# ── Compare to original model ─────────────────────────────────
aicc_nested_vals <- MuMIn::AICc(mF2_tweedie, mF2_tweedie_nested)

aicc_nested_df <- data.frame(
  Model  = rownames(aicc_nested_vals),
  df     = aicc_nested_vals$df,
  AICc   = round(aicc_nested_vals$AICc, 2),
  row.names = NULL
) %>%
  mutate(
    Delta  = round(AICc - min(AICc), 2),
    Weight = round(exp(-0.5 * Delta) / sum(exp(-0.5 * Delta)), 4)
  ) %>%
  arrange(AICc)

print(aicc_nested_df)

# again non-nested wins