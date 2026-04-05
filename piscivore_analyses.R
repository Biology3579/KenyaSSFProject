# ============================================================
#  PISCIVORE BIOMASS 
#  Response:   Piscivore biomass (continuous, zero-inflated)
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

# ── PREPARE PISCIVORE MODEL DATA ─────────────────────────────────────────────────
# Replace these paths with your own files.
# Expected: one row per transect, columns for biomass per group + all predictors.

fish_2009 <- readr::read_rds(here::here("processed_data", "clean_fish_2009.rds"))
location_2009 <- readr::read_rds(here::here("processed_data", "clean_location_2009.rds"))
gravity_2009 <- readr::read_rds(here::here("city_data", "locations_with_grav_combined.rds"))
old_grav_2009 <- read.csv(here::here("processed_data", "locations_with_gravity.csv"))
chla_2009 <- read.csv(here::here("processed_data", "locations_with_chla_2009.csv"))
sst_2009 <- read.csv(here::here("processed_data", "locations_with_sst_2009.csv"))
rugosity_2009 <- readr::read_rds(here::here("processed_data", "clean_dive_details_2009.rds"))
# connectivity_pisc <- read.csv(here::here("connectivity", "primary_piscivores.csv"))
  
# ============================================================
#  PISCIVORE BIOMASS PER TRANSECT
# ============================================================

pisc_transects <- fish_2009 %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_pisc_biomass = sum(
      ifelse(trophic_group == "piscivores", tot_wt_g, 0),
      na.rm = TRUE
    ),
    country = first(country),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(site = as.factor(site)
  )

cat("Transects after filtering:", nrow(pisc_transects), "\n")
cat("Sites after filtering:",     n_distinct(pisc_transects$site), "\n")


# ============================================================
#  PISCIVORE BIOMASS EXPLORATION
# ============================================================

# ── Basic summary ─────────────────────────────────────────────
summary(pisc_transects$transect_pisc_biomass)

zeros <- mean(pisc_transects$transect_pisc_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")
# since 0.475 → strong case for Tweedie or zero-inflated family

# ── Raw distribution ──────────────────────────────────────────
( pisc_raw <- ggplot(pisc_transects, aes(x = transect_pisc_biomass)) +
  geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
  labs(x = "Piscivore biomass per transect (g)", y = "Frequency",
       title = "Raw Piscivore Biomass") +
  theme_bw() )

# ── Log transformation ────────────────────────────────────────
# Small constant (+0.01) handles zeros; check sensitivity to constant choice
pisc_transects <- pisc_transects %>%
  mutate(
    log_pisc_biomass  = log(transect_pisc_biomass + 0.01),
    sqrt_pisc_biomass = sqrt(transect_pisc_biomass)
  )

( pisc_log <- ggplot(pisc_transects, aes(x = log_pisc_biomass)) +
  geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
  labs(x = "log(biomass + 0.01)", y = "Frequency",
       title = "Log-transformed Piscivore Biomass") +
  theme_bw() )

( pisc_sqrt <- ggplot(pisc_transects, aes(x = sqrt_pisc_biomass)) +
  geom_histogram(bins = 30, fill = "#d7191c", colour = "white") +
  labs(x = "sqrt(biomass)", y = "Frequency",
       title = "Sqrt-transformed Piscivore Biomass") +
  theme_bw() )

jpeg("pisc_biomass_distributions.jpg", width = 33, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(pisc_raw, pisc_log, pisc_sqrt, ncol = 3)
dev.off()

# ── Box-Cox: what power transformation does the data suggest? ─
# Fit a simple lm on non-zero values only (boxcox requires y > 0)
pisc_nonzero <- pisc_transects %>% filter(transect_pisc_biomass > 0)

MASS::boxcox(
  lm(transect_pisc_biomass ~ 1, data = pisc_nonzero),
  lambda = seq(-2, 2, 0.1)
)
# lambda ~ 0  → log transformation appropriate

# ── Variation by site ─────────────────────────────────────────
ggplot(pisc_transects, aes(x = reorder(site, transect_pisc_biomass, median),
                           y = transect_pisc_biomass)) +
  geom_boxplot(fill = "#2c7bb6", alpha = 0.6,
               outlier.colour = "black",
               outlier.size = 1) +
  coord_flip() +
  labs(x = NULL, y = "Piscivore biomass (g)",
       title = "Biomass distribution by site") +
  theme_bw(base_size = 9)

# ── Overall zero inflation check ──────────────────────────────────────
zeros <- mean(pisc_transects$transect_pisc_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")

# 0.47% => will likely need a Tweedie model
# 
# ── Zeros by site ─────────────────────────────────────────────
pisc_transects %>%
  group_by(site) %>%
  summarise(
    n_transects  = n(),
    prop_zeros   = mean(transect_pisc_biomass == 0),
    mean_biomass = mean(transect_pisc_biomass),
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
    market_gravity = mean(market_grav, na.rm = TRUE),
    settlement_pop = mean(settlement_tot_pop, na.rm = TRUE),
    settlement_grav = mean(nearest_pop75_grav, na.rm = TRUE),
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

pisc_model_data <- pisc_transects %>%
  left_join(gravity_sites,  by = "site") %>%
  left_join(chla_sites,     by = "site") %>%
  left_join(sst_sites,      by = "site") %>%
  left_join(rugosity_sites, by = "site")

pisc_transects %>%
  left_join(gravity_sites,  by = "site") %>%
  left_join(chla_sites,     by = "site") %>%
  left_join(sst_sites,      by = "site") %>%
  left_join(rugosity_sites, by = "site") %>%
  dplyr::select(market_gravity, settlement_grav, settlement_pop,
                mean_annual_chla, mean_annual_sst, rugosity) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
  facet_wrap(~ variable, scales = "free") +
  theme_bw()

# LOG-TRANSFORM & SCALE PREDICTORS --------------------------

pisc_model_data <- pisc_model_data %>%
  mutate(
    # Log-transformations (right-skewed variables)
    log_market_gravity = log(market_gravity),
    log_settlement_grav = log(settlement_grav),
    log_settlement_pop = log(settlement_pop ),
    log_chla = log(mean_annual_chla),
    
    # Scaled versions (mean = 0, SD = 1) — use these in models
    log_market_gravity_sc = as.numeric(scale(log_market_gravity)),
    log_settlement_grav_sc = as.numeric(scale(log_settlement_grav)),
    log_settlement_pop_sc = as.numeric(scale(log_settlement_pop)),
    log_chla_sc = as.numeric(scale(log_chla)),
    sst_sc = as.numeric(scale(mean_annual_sst)),
    rugosity_sc = as.numeric(scale(rugosity)),
    
    # Log-transform response for diagnostics / gaussian alternative
    log_pisc_biomass        = log(transect_pisc_biomass + 0.01)
  )

# Check transformations
pisc_model_data %>%
  dplyr::select(log_market_gravity, log_settlement_grav, log_settlement_pop,
                log_chla, mean_annual_sst, rugosity) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
  facet_wrap(~ variable, scales = "free") +
  labs(title = "Transformed predictors") +
  theme_bw()

# ============================================================
# PREDICTOR CORRELATION MATRIX
# ============================================================

predictor_vars <- pisc_model_data %>%
  dplyr::select(
    log_market_gravity_sc, log_settlement_grav_sc, log_settlement_pop_sc,
    log_chla_sc, sst_sc, rugosity_sc
  ) %>%
  rename(
    "Market gravity"      = log_market_gravity_sc,
    "Settlement gravity"  = log_settlement_grav_sc,
    "Settlement pop."     = log_settlement_pop_sc,
    "Chlorophyll-a"       = log_chla_sc,
    "SST"                 = sst_sc,
    "Rugosity"            = rugosity_sc
  )

corr_matrix <- cor(predictor_vars, use = "complete.obs")
p_mat       <- cor_pmat(corr_matrix)

corrplot(corr_matrix,
         method      = "square",
         type        = "lower",
         tl.col      = "black",
         tl.srt      = 0,
         tl.offset = 0.5,    # ← pushes labels away from matrix
         addCoef.col = "black",
         number.cex  = 0.8,
         col         = colorRampPalette(c("#d73027", "white", "#4575b4"))(200),
         mar = c(0, 0, 4, 2))

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
# The simplest possible approach. Valid only if zeros are rare
# and the log-transformed response looks approximately normal.
# The constant (+0.01) is arbitrary — check sensitivity if used.

mF1_gaussian <- glmmTMB(
  log_pisc_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = gaussian(),
  data   = pisc_model_data
)

resF1 <- simulateResiduals(mF1_gaussian, n = 1000)

jpeg("dharma_F1_gaussian.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF1, main = "DHARMa — Gaussian on log(y + 0.01)")
dev.off()

plot(resF1)
testZeroInflation(resF1)
testDispersion(resF1)

# ── STEP F2: Tweedie ──────────────────────────────────────────
# Compound Poisson-Gamma distribution. 
# Natively produces exact # zeros alongside a continuous positive distribution. No need
# for a log-transformation of the response. The power parameter p
# (estimated internally, typically 1 < p < 2) controls the
# zero mass and the shape of positive values.

mF2_tweedie <- glmmTMB(
  transect_pisc_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = pisc_model_data
)

resF2 <- simulateResiduals(mF2_tweedie, n = 1000)

jpeg("dharma_F2_tweedie.jpg", width = 25, height = 15, units = "cm", res = 300)
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
  transect_pisc_biomass ~ sst_sc + log_chla_sc +
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
#
# F1 Gaussian on log(y + 0.01): REJECTED
#   KS p = 0.00019; systematic misfit from zero mass.
#   ZI test returns NaN — model has no zero-generating mechanism.
#
# F2 Plain Tweedie: CHOSEN
#   KS p = 0.44, dispersion p = 0.87, ZI p = 0.68 — all n.s.
#   AICc = 2861.87 (df = 8).
#
# F3 Zero-inflated Tweedie: NOT JUSTIFIED
#   Formal tests near-identical to F2: ZI p = 0.69,
#   dispersion p = 0.90. AICc = 2864.01 (df = 9),
#   delta AICc = +2.13 relative to F2. The ZI component
#   adds one parameter with no improvement in any formal
#   diagnostic and a worse AICc. Visual inspection of the
#   residuals vs predicted panel suggests marginally smoother
#   quantile tracking under F3, but this is not supported
#   by any test and does not outweigh the parsimony penalty.
#
# → Plain Tweedie (F2) is the chosen family.
#   All candidate models use tweedie(link = "log"),
#   no ziformula.

# ============================================================
#  CANDIDATE MODEL FITTING
#  Response:   transect_pisc_biomass (zero-inflated, continuous)
#  Family:     (plain) Tweedie (glmmTMB)
#  Random fx:  (1 | site)
#  Predictors: sst_sc, log_chla_sc, log_dhw_sc, rugosity_sc,
#              log_market_gravity_sc OR log_settlement_grav_sc
#              (not both — r = 0.56, compared separately via AICc)
# ============================================================

#  ── FIT CANDIDATE MODELS ─────────────────────────────────

m0 <- glmmTMB(transect_pisc_biomass ~ 1 + (1 | site),
              family = tweedie(link = "log"), data = pisc_model_data)

m_env <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc +
                   (1 | site), family = tweedie(link = "log"), data = pisc_model_data)

m_market <- glmmTMB(transect_pisc_biomass ~ log_market_gravity_sc + (1 | site),
                    family = tweedie(link = "log"), data = pisc_model_data)

m_settgrav <- glmmTMB(transect_pisc_biomass ~ log_settlement_grav_sc + (1 | site),
                      family = tweedie(link = "log"), data = pisc_model_data)

m_hab <- glmmTMB(transect_pisc_biomass ~ rugosity_sc + (1 | site),
                 family = tweedie(link = "log"), data = pisc_model_data)

m_env_market <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc +
                          log_market_gravity_sc + (1 | site),
                        family = tweedie(link = "log"), data = pisc_model_data)

m_env_settgrav <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc +
                            log_settlement_grav_sc + (1 | site),
                          family = tweedie(link = "log"), data = pisc_model_data)

m_habitat_market <- glmmTMB(transect_pisc_biomass ~ rugosity_sc + log_market_gravity_sc + (1 | site),
                            family = tweedie(link = "log"), data = pisc_model_data)

m_habitat_settgrav <- glmmTMB(transect_pisc_biomass ~ rugosity_sc + log_settlement_grav_sc + (1 | site),
                              family = tweedie(link = "log"), data = pisc_model_data)

m_full_market <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc +
                           log_market_gravity_sc + rugosity_sc +
                           (1 | site), family = tweedie(link = "log"), data = pisc_model_data)

m_full_settgrav <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc +
                             log_settlement_grav_sc + rugosity_sc +
                             (1 | site), family = tweedie(link = "log"), data = pisc_model_data)

#  ── AICc TABLE ─────────────────────────────────────────────

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
delta      <- aicc_vals - min(aicc_vals)
weights    <- exp(-0.5 * delta) / sum(exp(-0.5 * delta))

aicc_df <- data.frame(
  Model  = names(model_list),
  AICc   = round(aicc_vals, 2),
  Delta  = round(delta, 2),
  Weight = round(weights, 4),
  row.names = NULL
) %>% arrange(AICc)

print(aicc_df)
# Models within Delta < 2 have substantial support.
# Models within Delta < 7 have some support.
# The gravity metric (market vs settlement) with lower AICc
# is the one to carry forward into further analyses.

#  ── DHARMA RESIDUAL DIAGNOSTICS ON BEST MODELS ──────────────

res <- simulateResiduals(m_full_market, n = 1000)

jpeg("dharma_pisc__full_market_model.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res, main = "DHARMa diagnostics — piscivore best model")
dev.off()

# Specific tests to check
testZeroInflation(res)   # is the ZI component doing its job (fitting zeros correctly?)
testDispersion(res)      # overdispersion?
testOutliers(res)        # influential observations (e.g. dindini)?
plot(res)


# Plot residuals against each predictor separately
plotResiduals(res, pisc_model_data$sst_sc, xlab = "SST")
plotResiduals(res, pisc_model_data$log_chla_sc, xlab = "Chl-a")
plotResiduals(res, pisc_model_data$rugosity_sc, xlab = "Rugosity")
plotResiduals(res, pisc_model_data$log_market_gravity_sc, xlab = "Gravity")

# Check if site is driving the overall residual pattern
plotResiduals(res, pisc_model_data$site, 
              xlab = "Site")

# Check residuals against any variables NOT in the model
# e.g. depth, year, country — whatever you have available
plotResiduals(res, pisc_model_data$country,
              xlab = "Country")

#  ── SUMMARY ────────────────────────────────────────

summary(m_full_market)

# ============================================================
#  MARGINAL EFFECT PLOTS
#  Population-level predictions holding all other predictors
#  at their mean (= 0 on the scaled axis). i.e. he partial effect of one
#  variable while marginalising over the others.
#  RE set to zero (re.form = NA) so curves reflect the average site.
# ============================================================

plot_effect <- function(model, data, focal_var, x_label,
                        y_label = "Piscivore biomass (g)",
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

p_sst  <- plot_effect(m_full_market, pisc_model_data,
                      "sst_sc",                 "Sea surface temperature (scaled)")
p_chla <- plot_effect(m_full_market, pisc_model_data,
                      "log_chla_sc",            "Chlorophyll-a (scaled log)")
p_grav <- plot_effect(m_full_market, pisc_model_data,
                      "log_market_gravity_sc",  "Market gravity (scaled log)")
p_rug  <- plot_effect(m_full_market, pisc_model_data,
                      "rugosity_sc",            "Rugosity (scaled)")

jpeg("effects_piscivores.jpg", width = 35, height = 14, units = "cm", res = 300)
gridExtra::grid.arrange(p_sst, p_chla, p_grav, p_rug, ncol = 2)
dev.off()


# ============================================================
# COUNTRY EFFECTS
# Refit with nested random effects: site within country
# Adds a country-level random intercept to absorb the between-
# country variance flagged by the Levene test above.
# ============================================================

mF2_tweedie_nested <- glmmTMB(
  transect_pisc_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | country/site),
  family = tweedie(link = "log"),
  data   = pisc_model_data
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

# ── DHARMa diagnostics on nested model ───────────────────────
res_nested <- simulateResiduals(mF2_tweedie_nested, n = 1000)

jpeg("dharma_tweedie_nested.jpg", width = 25, height = 15,
     units = "cm", res = 300)
plot(res_nested, main = "DHARMa — Tweedie (1 | country/site)")
dev.off()

plot(res_nested)
testZeroInflation(res_nested)
testDispersion(res_nested)

# ── Check if country-level variance is resolved ───────────────
plotResiduals(res_nested, pisc_model_data$country,
              xlab = "Country")
plotResiduals(res_nested, pisc_model_data$sst_sc,
              xlab = "SST")
plotResiduals(res_nested, pisc_model_data$log_chla_sc,
              xlab = "Chl-a")
plotResiduals(res_nested, pisc_model_data$log_market_gravity_sc,
              xlab = "Market gravity")
plotResiduals(res_nested, pisc_model_data$rugosity_sc,
              xlab = "Rugosity")

# ============================================================
#  CANDIDATE MODEL FITTING
#  Response:   transect_pisc_biomass (continuous, zero-inflated)
#  Family:     Tweedie (link = "log") — justified via DHARMa
#              and AICc comparison against Gaussian and ZI Tweedie
#  Random fx:  (1 | country/site) — nested structure justified
#              by significant country-level residual variance
#              (country SD = 0.55, ~35% of combined RE variance)
#              and AICc improvement of 2.31 points over (1 | site)
#  Predictors: sst_sc, log_chla_sc, rugosity_sc,
#              log_market_gravity_sc OR log_settlement_grav_sc
#              (not both — r = 0.56, compared separately via AICc)
# ============================================================

m0_nested <- glmmTMB(transect_pisc_biomass ~ 1 + (1 | country/site),
              family = tweedie(link = "log"), data = pisc_model_data)

m_env_nested <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc +
                   (1 | country/site),
                 family = tweedie(link = "log"), data = pisc_model_data)

m_market_nested <- glmmTMB(transect_pisc_biomass ~ log_market_gravity_sc +
                      (1 | country/site),
                    family = tweedie(link = "log"), data = pisc_model_data)

m_settgrav_nested <- glmmTMB(transect_pisc_biomass ~ log_settlement_grav_sc +
                        (1 | country/site),
                      family = tweedie(link = "log"), data = pisc_model_data)

m_hab_nested <- glmmTMB(transect_pisc_biomass ~ rugosity_sc +
                   (1 | country/site),
                 family = tweedie(link = "log"), data = pisc_model_data)

m_env_market_nested <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc +
                          log_market_gravity_sc + (1 | country/site),
                        family = tweedie(link = "log"), data = pisc_model_data)

m_env_settgrav_nested <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc +
                            log_settlement_grav_sc + (1 | country/site),
                          family = tweedie(link = "log"), data = pisc_model_data)

m_full_market_nested <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc +
                           log_market_gravity_sc + rugosity_sc +
                           (1 | country/site),
                         family = tweedie(link = "log"), data = pisc_model_data)

m_full_settgrav_nested <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc +
                             log_settlement_grav_sc + rugosity_sc +
                             (1 | country/site),
                           family = tweedie(link = "log"), data = pisc_model_data)

# ── AICc TABLE ───────────────────────────────────────────────
model_list_nested <- list(
  "Null"                     = m0_nested,
  "Environment"              = m_env_nested,
  "Market gravity"           = m_market_nested,
  "Settlement gravity"       = m_settgrav_nested,
  "Habitat"                  = m_hab_nested,
  "Env + market gravity"     = m_env_market_nested,
  "Env + settlement gravity" = m_env_settgrav_nested,
  "Full (market gravity)"    = m_full_market_nested,
  "Full (settlement gravity)"= m_full_settgrav_nested
)

aicc_vals_nested <- sapply(model_list_nested, MuMIn::AICc)
delta_nested     <- aicc_vals - min(aicc_vals)
weights_nested   <- exp(-0.5 * delta) / sum(exp(-0.5 * delta))

aicc_df_nested <- data.frame(
  Model  = names(model_list_nested),
  AICc   = round(aicc_vals_nested, 2),
  Delta  = round(delta_nested,     2),
  Weight = round(weights_nested,   4),
  row.names = NULL
) %>% arrange(AICc)

print(aicc_df_nested)
# Models within Delta < 2 have substantial support.
# Models within Delta < 7 have some support.
# The gravity metric (market vs settlement) with lower AICc
# is the one to carry forward into further analyses.

# If AICc improves and country Levene test is no longer
# significant → use nested RE structure for all candidate
# models. If AICc is similar but diagnostics improve →
# still worth adopting on grounds of better model structure.

# ============================================================
#  MARGINAL EFFECT PLOTS
#  Population-level predictions holding all other predictors
#  at their mean (= 0 on the scaled axis). i.e. the partial effect of one
#  variable while marginalising over the others.
#  RE set to zero (re.form = NA) so curves reflect the average site
#  within an average country.
# ============================================================

p_sst_nested  <- plot_effect(mF2_tweedie_nested, pisc_model_data,
                             "sst_sc",                "Sea surface temperature (scaled)")
p_chla_nested <- plot_effect(mF2_tweedie_nested, pisc_model_data,
                             "log_chla_sc",           "Chlorophyll-a (scaled log)")
p_grav_nested <- plot_effect(mF2_tweedie_nested, pisc_model_data,
                             "log_market_gravity_sc", "Market gravity (scaled log)")
p_rug_nested  <- plot_effect(mF2_tweedie_nested, pisc_model_data,
                             "rugosity_sc",           "Rugosity (scaled)")

jpeg("effects_piscivores_nested.jpg", width = 35, height = 14, units = "cm", res = 300)
gridExtra::grid.arrange(p_sst_nested, p_chla_nested, 
                        p_grav_nested, p_rug_nested, ncol = 2)
dev.off()


# ... site levels with connectivity locations 
