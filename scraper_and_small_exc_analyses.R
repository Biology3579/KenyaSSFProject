# ============================================================
#  SCRAPER / SMALL EXCAVATOR BIOMASS
#  Response:   Scraper + small excavator biomass (continuous, zero-inflated)
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

# ── PREPARE SCRAPER / SMALL EXCAVATOR MODEL DATA ─────────────────────────────────────────────────

fish_2009 <- readr::read_rds(here::here("processed_data", "clean_fish_2009.rds"))
location_2009 <- readr::read_rds(here::here("processed_data", "clean_location_2009.rds"))
gravity_2009 <- readr::read_rds(here::here("city_data", "locations_with_grav_combined.rds"))
old_grav_2009 <- read.csv(here::here("processed_data", "locations_with_gravity.csv"))
chla_2009 <- read.csv(here::here("processed_data", "locations_with_chla_2009.csv"))
sst_2009 <- read.csv(here::here("processed_data", "locations_with_sst_2009.csv"))
rugosity_2009 <- readr::read_rds(here::here("processed_data", "clean_dive_details_2009.rds"))

# ============================================================
#  SCRAPER / SMALL EXCAVATOR BIOMASS PER TRANSECT
# ============================================================

scraper_transects <- fish_2009 %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_scraper_biomass = sum(
      ifelse(trophic_group %in% c("scrapers", "small_excavators"), tot_wt_g, 0),
      na.rm = TRUE
    ),
    country = first(country),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(site = as.factor(site))

cat("Transects after filtering:", nrow(scraper_transects), "\n")
cat("Sites after filtering:",     n_distinct(scraper_transects$site), "\n")


# ============================================================
#  SCRAPER / SMALL EXCAVATOR BIOMASS EXPLORATION
# ============================================================

# ── Basic summary ─────────────────────────────────────────────
summary(scraper_transects$transect_scraper_biomass)

zeros <- mean(scraper_transects$transect_scraper_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")
# Check proportion — strong zero inflation may favour Tweedie or zero-inflated family

# ── Raw distribution ──────────────────────────────────────────
( scraper_raw <- ggplot(scraper_transects, aes(x = transect_scraper_biomass)) +
    geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
    labs(x = "Scraper / small excavator biomass per transect (g)", y = "Frequency",
         title = "Raw Scraper / Small Excavator Biomass") +
    theme_bw() )

# ── Log transformation ────────────────────────────────────────
# Small constant (+0.01) handles zeros; check sensitivity to constant choice
scraper_transects <- scraper_transects %>%
  mutate(
    log_scraper_biomass  = log(transect_scraper_biomass + 0.01),
    sqrt_scraper_biomass = sqrt(transect_scraper_biomass)
  )

( scraper_log <- ggplot(scraper_transects, aes(x = log_scraper_biomass)) +
    geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
    labs(x = "log(biomass + 0.01)", y = "Frequency",
         title = "Log-transformed Scraper / Small Excavator Biomass") +
    theme_bw() )

( scraper_sqrt <- ggplot(scraper_transects, aes(x = sqrt_scraper_biomass)) +
    geom_histogram(bins = 30, fill = "#d7191c", colour = "white") +
    labs(x = "sqrt(biomass)", y = "Frequency",
         title = "Sqrt-transformed Scraper / Small Excavator Biomass") +
    theme_bw() )

jpeg("scraper_biomass_distributions.jpg", width = 33, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(scraper_raw, scraper_log, scraper_sqrt, ncol = 3)
dev.off()

# ── Box-Cox: what power transformation does the data suggest? ─
# Fit a simple lm on non-zero values only (boxcox requires y > 0)
scraper_nonzero <- scraper_transects %>% filter(transect_scraper_biomass > 0)

MASS::boxcox(
  lm(transect_scraper_biomass ~ 1, data = scraper_nonzero),
  lambda = seq(-2, 2, 0.1)
)
# lambda ~ 0  → log transformation appropriate

# ── Variation by site ─────────────────────────────────────────
ggplot(scraper_transects, aes(x = reorder(site, transect_scraper_biomass, median),
                              y = transect_scraper_biomass)) +
  geom_boxplot(fill = "#2c7bb6", alpha = 0.6,
               outlier.colour = "black",
               outlier.size = 1) +
  coord_flip() +
  labs(x = NULL, y = "Scraper / small excavator biomass (g)",
       title = "Biomass distribution by site") +
  theme_bw(base_size = 9)

# ── Overall zero inflation check ──────────────────────────────────────
zeros <- mean(scraper_transects$transect_scraper_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")

# ── Zeros by site ─────────────────────────────────────────────
scraper_transects %>%
  group_by(site) %>%
  summarise(
    n_transects  = n(),
    prop_zeros   = mean(transect_scraper_biomass == 0),
    mean_biomass = mean(transect_scraper_biomass),
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

scraper_model_data <- scraper_transects %>%
  left_join(gravity_sites,  by = "site") %>%
  left_join(chla_sites,     by = "site") %>%
  left_join(sst_sites,      by = "site") %>%
  left_join(rugosity_sites, by = "site")

# Set predictor order and labels (for plot)
predictor_order_raw  <- c("market_gravity", "settlement_grav", "settlement_pop",
                          "mean_annual_chla", "mean_annual_sst", "rugosity")

predictor_labels <- c("Market gravity", "Settlement gravity", "Settlement pop.",
                      "Chlorophyll-a", "SST", "Rugosity")

scraper_transects %>%
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

scraper_model_data <- scraper_model_data %>%
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
    log_scraper_biomass = log(transect_scraper_biomass + 0.01)
  )

# Set predictor order (for plot)
predictor_order_tran <- c("log_market_gravity", "log_settlement_grav", "log_settlement_pop",
                          "log_chla", "mean_annual_sst", "rugosity")
# Plot transformed data
scraper_model_data %>%
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

predictor_vars <- scraper_model_data %>%
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
  log_scraper_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = gaussian(),
  data   = scraper_model_data
)

resF1 <- simulateResiduals(mF1_gaussian, n = 1000)

jpeg("dharma_F1_gaussian_scrapers.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF1, main = "DHARMa — Gaussian on log(y + 0.01)")
dev.off()

plot(resF1)
testZeroInflation(resF1)
testDispersion(resF1)

# ── STEP F2: Tweedie ──────────────────────────────────────────
mF2_tweedie <- glmmTMB(
  transect_scraper_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = scraper_model_data
)

resF2 <- simulateResiduals(mF2_tweedie, n = 1000)

jpeg("dharma_F2_tweedie_scrapers.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF2, main = "DHARMa — Tweedie")
dev.off()

plot(resF2)
testZeroInflation(resF2)
testDispersion(resF2)

# ============================================================
#  CANDIDATE MODEL FITTING
#  Response:   transect_scraper_biomass (zero-inflated, continuous)
#  Family:     Tweedie (link = "log") — confirm via family selection above
#  Random fx:  (1 | site)
#  Predictors: sst_sc, log_chla_sc, rugosity_sc,
#              log_market_gravity_sc OR log_settlement_grav_sc
#              (not both — r = 0.56, compared separately via AICc)
# ============================================================

m0 <- glmmTMB(transect_scraper_biomass ~ 1 + (1 | site),
              family = tweedie(link = "log"), data = scraper_model_data)

m_env <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc +
                   (1 | site),
                 family = tweedie(link = "log"), data = scraper_model_data)

m_market <- glmmTMB(transect_scraper_biomass ~ log_market_gravity_sc + (1 | site),
                    family = tweedie(link = "log"), data = scraper_model_data)

m_settgrav <- glmmTMB(transect_scraper_biomass ~ log_settlement_grav_sc + (1 | site),
                      family = tweedie(link = "log"), data = scraper_model_data)

m_hab <- glmmTMB(transect_scraper_biomass ~ rugosity_sc + (1 | site),
                 family = tweedie(link = "log"), data = scraper_model_data)

m_env_market <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc +
                          log_market_gravity_sc + (1 | site),
                        family = tweedie(link = "log"), data = scraper_model_data)

m_env_settgrav <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc +
                            log_settlement_grav_sc + (1 | site),
                          family = tweedie(link = "log"), data = scraper_model_data)

m_habitat_market <- glmmTMB(transect_scraper_biomass ~ rugosity_sc + log_market_gravity_sc + (1 | site),
                            family = tweedie(link = "log"), data = scraper_model_data)

m_habitat_settgrav <- glmmTMB(transect_scraper_biomass ~ rugosity_sc + log_settlement_grav_sc + (1 | site),
                              family = tweedie(link = "log"), data = scraper_model_data)

m_full_market <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc +
                           log_market_gravity_sc + rugosity_sc +
                           (1 | site),
                         family = tweedie(link = "log"), data = scraper_model_data)

m_full_settgrav <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc +
                             log_settlement_grav_sc + rugosity_sc +
                             (1 | site),
                           family = tweedie(link = "log"), data = scraper_model_data)

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

res <- simulateResiduals(m_full_settgrav, n = 1000)

jpeg("dharma_scraper_full_market_model.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res, main = "DHARMa diagnostics — scraper / small excavator best model")
dev.off()

testZeroInflation(res)
testDispersion(res)
testOutliers(res)
plot(res)

plotResiduals(res, scraper_model_data$sst_sc,                 xlab = "SST")
plotResiduals(res, scraper_model_data$log_chla_sc,            xlab = "Chl-a")
plotResiduals(res, scraper_model_data$rugosity_sc,            xlab = "Rugosity")
plotResiduals(res, scraper_model_data$log_settlement_grav_sc,  xlab = "Gravity")

# ── SUMMARY ────────────────────────────────────────

summary(m_full_settgrav)

# ============================================================
#  MARGINAL EFFECT PLOTS
# ============================================================

plot_effect <- function(model, data, focal_var, x_label,
                        y_label = "Scraper / small excavator biomass (g)",
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

p_sst  <- plot_effect(m_full_market, scraper_model_data,
                      "sst_sc",                "Sea surface temperature (scaled)")
p_chla <- plot_effect(m_full_market, scraper_model_data,
                      "log_chla_sc",           "Chlorophyll-a (scaled log)")
p_grav <- plot_effect(m_full_market, scraper_model_data,
                      "log_market_gravity_sc", "Market gravity (scaled log)")
p_rug  <- plot_effect(m_full_market, scraper_model_data,
                      "rugosity_sc",           "Rugosity (scaled)")

jpeg("effects_scrapers.jpg", width = 35, height = 14, units = "cm", res = 300)
gridExtra::grid.arrange(p_sst, p_chla, p_grav, p_rug, ncol = 2)
dev.off()

# ============================================================
#  COUNTRY EFFECTS
# ============================================================

mF2_tweedie_nested <- glmmTMB(
  transect_scraper_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | country/site),
  family = tweedie(link = "log"),
  data   = scraper_model_data
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

# Here the country nested model does not improve so we stay with the other model 