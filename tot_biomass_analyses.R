# ============================================================
#  TOTAL FISH BIOMASS & ABUNDANCE — MIXED EFFECTS MODELS
#
#  Study design:
#    Transects nested within stations, stations within sites,
#    sites within locations, locations within countries.
#
#  Predictors measured at site level (averaged from station):
#    SST, Chl-a, Human gravity (market / settlement), Rugosity
#
#  Analytical structure:
#    PART 1 — Site-level analysis (PRIMARY)
#              Matches response resolution to predictor resolution.
#              Sites are the true unit of environmental inference.
#              Sections: no RE → country RE → location RE
#
#    PART 2 — Transect-level biomass (SENSITIVITY CHECK)
#              Retains within-site variation; (1 | site) accounts
#              for non-independence. Confirms site-level findings
#              are not an artefact of averaging.
#
#    PART 3 — Transect-level counts (COMPLEMENTARY ANALYSIS)
#              Models the raw data-generating process (discrete
#              counts) rather than derived biomass. Negative
#              Binomial family. Allows detection of whether
#              predictor effects operate through abundance,
#              body size, or both.
# ============================================================

options(scipen = 999)

# ── PACKAGES ─────────────────────────────────────────────────
library(tidyverse)
library(sf)
library(glmmTMB)
library(DHARMa)
library(MuMIn)
library(AICcmodavg)
library(ggcorrplot)
library(corrplot)
library(gridExtra)
library(MASS)
library(here)

# ──  LOAD ALL DATA ─────────────────────────────────────────────────
fish_2009     <- readr::read_rds(here::here("processed_data", "clean_fish_2009.rds"))
gravity_2009  <- readr::read_rds(here::here("city_data", "locations_with_grav_combined.rds"))
chla_2009     <- read.csv(here::here("processed_data", "locations_with_chla_2009.csv"))
sst_2009      <- read.csv(here::here("processed_data", "locations_with_sst_2009.csv"))
rugosity_2009 <- readr::read_rds(here::here("processed_data", "clean_dive_details_2009.rds"))

# ── FUNCTIONS ─────────────────────────────────────────────────

# Build an AICc comparison table from a named list of models
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

# Generate a marginal effect plot for a single focal predictor.
# All other scaled predictors (_sc suffix) are held at 0 (their mean).
# Returns a ggplot object with a fitted line and 95% confidence ribbon.
plot_effect <- function(model, data, focal_var,                         x_label,
                        y_label = "Fitted value",
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
  
  # glmmTMB requires a valid site value; lm() does not use it
  if (inherits(model, "glmmTMB") && "site" %in% names(data)) {
    grid$site <- levels(data$site)[1]
  }
  
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
                fill = "#2c7bb6", alpha = 0.2) +
    geom_line(aes(y = fit), colour = "#2c7bb6", linewidth = 1.1) +
    labs(x = x_label, y = y_label) +
    theme_bw(base_size = 13) +
    theme(axis.title = element_text(face = "bold"))
}

# ==============================================================================
#  AGGREGATE TRANSECT DATA
#  Minimum of 3 transects per site.
# ==============================================================================

total_transects <- fish_2009 %>%
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
    site     = as.factor(site),
    country  = as.factor(country))

cat("Number of transects:", nrow(total_transects), "\n")
cat("Number of sites:",     n_distinct(total_transects$site), "\n")
# cat("Number of connectivity clusters:", n_distinct(total_transects$connectivity_clusters), "\n")
cat("Number of countries:", n_distinct(total_transects$country), "\n")

# ==============================================================================
#  BIOMASS DATA EXPLORATION
# ==============================================================================

# ── Aggregate data at site level (mean ts_biomass per site) ───────────────────
site_data <- total_transects %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass = mean(transect_total_biomass, na.rm = TRUE),
    n_transects  = n(),
    .groups      = "drop"
  ) %>%
  mutate(
    site     = as.factor(site),
    country  = as.factor(country)
  )

# ── Basic summary ─────────────────────────────────────────────────────────────
summary(site_data$mean_biomass)

zeros <- mean(site_data$mean_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")
# Expected: 0 — mean total biomass per site should never be zero

# ── Raw distribution ──────────────────────────────────────────────────────────
( site_raw <- ggplot(site_data, aes(x = mean_biomass)) +
    geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
    labs(x = "Mean total biomass per site (g)", y = "Frequency",
         title = "Raw Site-Level Biomass") +
    theme_bw() )

# ── Box-Cox: what power transformation does the data suggest? ─
MASS::boxcox(
  lm(mean_biomass ~ 1, data = site_data),
  lambda = seq(-2, 2, 0.1)
)
# lambda ~ 0 → log transformation appropriate

# ── Apply transformations ─────────────────────────────────────────────────────
site_data <- site_data %>%
  mutate(
    log_mean_biomass  = log(mean_biomass),
    sqrt_mean_biomass = sqrt(mean_biomass)
  )

# ── Log and sqrt distributions ────────────────────────────────────────────────
( site_log <- ggplot(site_data, aes(x = log_mean_biomass)) +
    geom_histogram(bins = 25, fill = "#1a9641", colour = "white") +
    labs(x = "log(mean biomass)", y = "Frequency",
         title = "Log-transformed Site-Level Biomass") +
    theme_bw() )

( site_sqrt <- ggplot(site_data, aes(x = sqrt_mean_biomass)) +
    geom_histogram(bins = 25, fill = "#d7191c", colour = "white") +
    labs(x = "sqrt(mean biomass)", y = "Frequency",
         title = "Sqrt-transformed Site-Level Biomass") +
    theme_bw() )

jpeg("site_biomass_distributions.jpg", width = 33, height = 11,
     units = "cm", res = 300)
gridExtra::grid.arrange(site_raw, site_log, site_sqrt, ncol = 3)
dev.off()

# ── Normality checks on log-transformed response ──────────────────────────────
qqnorm(site_data$log_mean_biomass,
       main = "Q-Q plot: log(mean biomass per site)")
qqline(site_data$log_mean_biomass, col = "red")
shapiro.test(site_data$log_mean_biomass)

# ── Variation by site ─────────────────────────────────────────────────────────
ggplot(site_data, aes(x = reorder(site, mean_biomass, median),
                      y = mean_biomass)) +
  geom_col(fill = "#2c7bb6", alpha = 0.7) +
  coord_flip() +
  labs(x = NULL, y = "Mean total biomass (g)",
       title = "Mean biomass by site (raw)") +
  theme_bw(base_size = 9)

# ── Variation by country ──────────────────────────────────────────────────────
ggplot(site_data, aes(x = country, y = log_mean_biomass)) +
  geom_boxplot(outlier.shape = NA, fill = "grey92",
               colour = "grey40", linewidth = 0.4, width = 0.4) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7,
              colour = "#2c7bb6") +
  geom_hline(yintercept = mean(site_data$log_mean_biomass),
             linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  scale_x_discrete(labels = stringr::str_to_title) +
  labs(x = NULL, y = "log(mean biomass per site)") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.text.x        = element_text(face = "bold")
  )

# ── Summary table ─────────────────────────────────────────────────────────────
site_data %>%
  dplyr::select(site, country, n_transects, mean_biomass) %>%
  arrange(desc(mean_biomass)) %>%
  print(n = Inf)

# ============================================================
#  PREDICTOR PREPARATION
# ============================================================

# ── Human gravity metrics ──────────────────────────────
gravity_sites <- gravity_2009 %>%
  st_drop_geometry() %>%
  group_by(site) %>%
  summarise(
    market_gravity  = mean(market_grav, na.rm = TRUE),
    settlement_pop  = mean(settlement_tot_pop, na.rm = TRUE),
    settlement_grav = mean(nearest_pop75_grav, na.rm = TRUE),
    .groups = "drop"
  )

# ── Chlorophyll-a ─────────────────────────────────────────────
chla_sites <- chla_2009 %>%
  group_by(site) %>%
  summarise(mean_annual_chla = mean(chla_annual_mean, na.rm = TRUE),
    .groups = "drop")

# ── Sea surface temperature ───────────────────────────────────
sst_sites <- sst_2009 %>%
  filter(!is.na(sst_annual_mean)) %>%
  group_by(site) %>%
  summarise(mean_annual_sst = mean(sst_annual_mean, na.rm = TRUE),
    .groups = "drop")

# ── Rugosity ──────────────────────────────────────────────────
rugosity_sites <- rugosity_2009 %>%
  group_by(site) %>%
  summarise(rugosity = mean(rugosity, na.rm = TRUE),
    .groups = "drop")

# ============================================================
#  TRANSFORMATIONS AND CHECKS
# ============================================================

# ── Inspect raw predictors ───────────────────────────────────────────────
# Combine raw predictors for inspection
raw_predictors <- gravity_sites %>%
  left_join(chla_sites,     by = "site") %>%
  left_join(sst_sites,      by = "site") %>%
  left_join(rugosity_sites, by = "site")

predictor_labels <- c("Market gravity", "Settlement gravity", "Settlement pop.",
                      "Chlorophyll-a", "SST", "Rugosity")

predictor_order_raw  <- c("market_gravity", "settlement_grav", "settlement_pop",
                          "mean_annual_chla", "mean_annual_sst", "rugosity")
# Plot raw predictor histograms
( p_pred_raw <- raw_predictors %>%
    dplyr::select(all_of(predictor_order_raw)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, predictor_order_raw, predictor_labels)) %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
    facet_wrap(~ variable, scales = "free") +
    labs(title = "Raw predictors") + theme_bw() )

# ── Inspect transformed predictors ───────────────────────────────────────────────
# Transform predictors 
# - note SST and rugosity are already relatively normal so won't be transformed
transformed_predictors <- raw_predictors %>%
  transmute(
    site                = site,
    log_market_gravity  = log(market_gravity),
    log_settlement_grav = log(settlement_grav),
    log_settlement_pop  = log(settlement_pop),
    log_chla            = log(mean_annual_chla),
    mean_annual_sst     = mean_annual_sst,
    rugosity            = rugosity
  )
predictor_order_tran <- c("log_market_gravity", "log_settlement_grav",
                          "log_settlement_pop", "log_chla",
                          "mean_annual_sst", "rugosity")

# Plot transformed histograms
( p_pred_tran <- transformed_predictors %>%
    dplyr::select(all_of(predictor_order_tran)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, predictor_order_tran, predictor_labels)) %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
    facet_wrap(~ variable, scales = "free") +
    labs(title = "Transformed predictors") + theme_bw() )

# Combined figure
jpeg("predictor_distributions.jpg", width = 33, height = 22, units = "cm", res = 300)
gridExtra::grid.arrange(p_pred_raw, p_pred_tran, nrow = 2)
dev.off()


# Combine and scale predictors for analyses
scaled_predictors <- transformed_predictors %>%
  transmute(
    site = site,
    log_market_gravity_sc  = as.numeric(scale(log_market_gravity)),
    log_settlement_grav_sc = as.numeric(scale(log_settlement_grav)),
    log_settlement_pop_sc  = as.numeric(scale(log_settlement_pop)),
    log_chla_sc            = as.numeric(scale(log_chla)),
    sst_sc                 = as.numeric(scale(mean_annual_sst)),
    rugosity_sc            = as.numeric(scale(rugosity))
  )

# ============================================================
#  PREDICTOR CORRELATION MATRIX
# ============================================================

corr_matrix <- scaled_predictors %>%
  dplyr::select(ends_with("_sc")) %>%
  rename(
    "Market gravity"     = log_market_gravity_sc,
    "Settlement gravity" = log_settlement_grav_sc,
    "Settlement pop."    = log_settlement_pop_sc,
    "Chlorophyll-a"      = log_chla_sc,
    "SST"                = sst_sc,
    "Rugosity"           = rugosity_sc
  ) %>%
  cor(use = "complete.obs")

corrplot(abs(corr_matrix),
         method      = "square",
         type        = "lower",
         tl.col      = "black",
         tl.srt      = 0,
         tl.offset   = 0.5,
         addCoef.col = "black",
         number.cex  = 0.8,
         col         = colorRampPalette(c("white", "#d73027"))(200),
         is.corr     = FALSE,
         mar         = c(0, 0, 4, 2))

# ── Decision guide ────────────────────────────────────────────
# Need to choose between gravity metrics but keep the rest
# Maybe beware of chla?

# ============================================================
#  CHOOSING SETTLEMENT METRIC
# ============================================================
# Settlement gravity and settlement pop both proxy local human
# pressure. Select the better-performing metric via AICc before
# entering the main candidate set.

# Data
settlement_data <- total_transects %>%
  left_join(scaled_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    log_mean_biomass       = log(mean(transect_total_biomass, na.rm = TRUE)),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_settlement_pop_sc  = first(log_settlement_pop_sc),
    .groups = "drop"
  )

settgrav <- lm(log_mean_biomass ~ log_settlement_grav_sc, data = settlement_data)
settpop  <- lm(log_mean_biomass ~ log_settlement_pop_sc,  data = settlement_data)

make_aicc_df(list(
  "Settlement gravity" = settgrav,
  "Settlement pop."    = settpop
))
# Settlement gravity selected (AICc = 146.46 vs 146.93; delta = 0.46).
# Near-equivalent performance (weights 0.56 vs 0.44) confirms results
# robust to metric choice. Settlement gravity carried forward;
# log_settlement_pop_sc dropped from all subsequent candidate models.

rm(settlement_data) # remove settlement_data - not used after this

# ============================================================
#  ANALYSIS DATASETS
# ============================================================


# Drop log_settlement_pop_sc now that settlement gravity is confirmed
final_predictors <- scaled_predictors %>%
  dplyr::select(site, log_market_gravity_sc, log_settlement_grav_sc,
                log_chla_sc, sst_sc, rugosity_sc)

# ── Transect-level dataset  ─────────────────────
transect_model_data <- total_transects %>%
  left_join(final_predictors, by = "site") %>%
  mutate(log_transect_biomass = log(transect_total_biomass))

cat("\nTransect model data:", nrow(transect_model_data), "rows,",
    n_distinct(transect_model_data$site), "sites\n")
cat("Biomass zeros:", sum(transect_model_data$transect_total_biomass == 0), "\n")
cat("Count zeros:",  sum(transect_model_data$transect_total_count   == 0), "\n")

# ── Site-level dataset ───────────────────────────────
total_model_data <- transect_model_data %>%
  group_by(site, country) %>%
  summarise(
    log_mean_biomass       = log(mean(transect_total_biomass,  na.rm = TRUE)),
    mean_biomass           = mean(transect_total_biomass,   na.rm = TRUE),
    n_transects            = n(),
    log_market_gravity_sc  = first(log_market_gravity_sc),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_chla_sc            = first(log_chla_sc),
    sst_sc                 = first(sst_sc),
    rugosity_sc            = first(rugosity_sc),
    .groups = "drop"
  ) %>%
  mutate(
    site    = as.factor(site),
    country = as.factor(country)
  )

cat("\nSite model data:", nrow(total_model_data), "sites,",
    n_distinct(total_model_data$country), "countries\n")

# ============================================================
#  FAMILY SELECTION
# ============================================================

# ── F1: Gaussian on raw mean biomass ─────────────────────────
mS_F1 <- lm(mean_biomass ~ sst_sc + log_chla_sc +
              log_market_gravity_sc + rugosity_sc,
            data = total_model_data)

jpeg("diagnostics_site_F1_gaussian_raw.jpg",
     width = 25, height = 20, units = "cm", res = 300)
par(mfrow = c(2, 2))
plot(mS_F1)
par(mfrow = c(1, 1))
dev.off()

par(mfrow = c(2, 2)); plot(mS_F1); par(mfrow = c(1, 1))

# ── F2: Gaussian on log-transformed mean biomass ─────────────
mS_F2 <- lm(log_mean_biomass ~ sst_sc + log_chla_sc +
              log_market_gravity_sc + rugosity_sc,
            data = total_model_data)

jpeg("diagnostics_site_F2_gaussian_log.jpg",
     width = 25, height = 20, units = "cm", res = 300)
par(mfrow = c(2, 2))
plot(mS_F2)
par(mfrow = c(1, 1))
dev.off()

par(mfrow = c(2, 2)); plot(mS_F2); par(mfrow = c(1, 1))

# ── F3: Gamma (log link) on raw mean biomass ─────────────────
# Natural fit for continuous positive data; variance scales
# with mean squared.
mS_F3 <- glm(mean_biomass ~ sst_sc + log_chla_sc +
               log_market_gravity_sc + rugosity_sc,
             family = Gamma(link = "log"),
             data   = total_model_data)

jpeg("diagnostics_site_F3_gamma.jpg",
     width = 25, height = 20, units = "cm", res = 300)
par(mfrow = c(2, 2))
plot(mS_F3)
par(mfrow = c(1, 1))
dev.off()

par(mfrow = c(2, 2)); plot(mS_F3); par(mfrow = c(1, 1))

# ── Summary comparison ────────────────────────────────────────
# AICc comparable between F1 and F3 (same response: mean_biomass)
# AICc NOT comparable between F2 and F1/F3 (different response)
make_aicc_df(list(
  "Gaussian (raw)"  = mS_F1,
  "Gamma (log)"     = mS_F3
))
# Assess F2 on diagnostics only — look for:
#   Residuals vs Fitted: no pattern (homoscedasticity)
#   Q-Q: points on line (normality)
#   Scale-Location: flat (constant variance)
#   Residuals vs Leverage: no high-influence points


# ── Family selection decision ─────────────────────────────────
# F2 selected: Gaussian on log_mean_biomass.
# Log transformation fully stabilises variance (flat Scale-Location)
# and satisfies normality (Q-Q closely follows theoretical line).
# Gamma (F3) offers lower AICc on the raw scale but is not
# comparable to F2 (different response). F2 diagnostics are
# cleaner and the model is more interpretable. Proceed with lm().


# ============================================================
#  RANDOM EFFECT STRUCTURE SELECTION
# ============================================================

#glmmTMB structure for better comparisons
re_null <- glmmTMB(log_mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + rugosity_sc,
                   family = gaussian(), data = total_model_data) 

re_country <- glmmTMB(log_mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + rugosity_sc +
                         (1 | country), family = gaussian(), data = total_model_data)

re_comparison <- list(
  "No RE"                  = re_null,
  "(1 | country)"          = re_country
)

cat("\n--- RE structure comparison (full market gravity model) ---\n")
print(make_aicc_df(re_comparison))

# ── RE structure decision ─────────────────────────────────────
# No RE selected (AICc weight = 0.78 vs 0.22 for country RE;
# delta AICc = 2.53). Country-level clustering is not supported
# once environmental and human pressure predictors are included.
# All candidate models fitted as lm() — no random effects.

# ============================================================
#  CANDIDATE MODELS
# ============================================================

# ── No random effects  ────────────────────────────────────────────────────────
# Sites as independent observations.

s1_m0 <- lm(log_mean_biomass ~ 1,data = total_model_data) 
s1_m_env <- lm(log_mean_biomass ~ sst_sc + log_chla_sc, data = total_model_data)
s1_m_market <- lm(log_mean_biomass ~ log_market_gravity_sc, data = total_model_data)
s1_m_settgrav <- lm(log_mean_biomass ~ log_settlement_grav_sc, data = total_model_data)
s1_m_hab <- lm(log_mean_biomass ~ rugosity_sc, data = total_model_data)
s1_m_env_mkt <- lm(log_mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc, data = total_model_data)
s1_m_env_sett <- lm(log_mean_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc, data = total_model_data)
s1_m_habitat_market <- lm(log_mean_biomass ~ rugosity_sc + log_market_gravity_sc, data = total_model_data)
s1_m_habitat_settgrav <- lm(log_mean_biomass ~ rugosity_sc + log_settlement_grav_sc, data = total_model_data)
s1_m_full_mkt  <- lm(log_mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc, data = total_model_data)
s1_m_full_sett <- lm(log_mean_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc, data = total_model_data)

model_list_s1 <- list(
  "Null"                      = s1_m0,
  "Environment"               = s1_m_env,
  "Market gravity"            = s1_m_market,
  "Settlement gravity"        = s1_m_settgrav,
  "Habitat"                   = s1_m_hab,
  "Env + market gravity"      = s1_m_env_mkt,
  "Env + settlement gravity"  = s1_m_env_sett,
  "Habitat + market gravity" = s1_m_habitat_market,
  "Habitat + settlement gravity" = s1_m_habitat_settgrav,
  "Full (market gravity)"     = s1_m_full_mkt,
  "Full (settlement gravity)" = s1_m_full_sett
)

cat("\n--- AICc: Site-level candidate models ---\n")
aicc_site <- make_aicc_df(model_list_s1)
print(aicc_site)

# ── Residual diagnostics — top three models ───────────────────

# Model 1: Habitat
cat("\n--- Habitat ---\n")
par(mfrow = c(2, 2)); plot(s1_m_hab); par(mfrow = c(1, 1))
summary(s1_m_hab)

jpeg("diagnostics_site_habitat.jpg", width = 25, height = 20, units = "cm", res = 300)
par(mfrow = c(2, 2)); plot(s1_m_hab); par(mfrow = c(1, 1))
dev.off()

# Model 2: Habitat + settlement gravity
cat("\n--- Habitat + settlement gravity ---\n")
par(mfrow = c(2, 2)); plot(s1_m_habitat_settgrav); par(mfrow = c(1, 1))
summary(s1_m_habitat_settgrav)

jpeg("diagnostics_site_habitat_settgrav.jpg", width = 25, height = 20, units = "cm", res = 300)
par(mfrow = c(2, 2)); plot(s1_m_habitat_settgrav); par(mfrow = c(1, 1))
dev.off()


# Model 3: Habitat + market gravity
cat("\n--- Habitat + market gravity ---\n")
par(mfrow = c(2, 2)); plot(s1_m_habitat_market); par(mfrow = c(1, 1))
summary(s1_m_habitat_market)

jpeg("diagnostics_site_habitat_market.jpg", width = 25, height = 20, units = "cm", res = 300)
par(mfrow = c(2, 2)); plot(s1_m_habitat_market); par(mfrow = c(1, 1))
dev.off()

# ── Marginal effect plots — site level ───────────────────────
# Predictions holding all other predictors at their mean
# (= 0 on the scaled axis).

# Rugosity — present in all three top models
( p_site_rugosity <- plot_effect(s1_m_hab,
                                      total_model_data,
                                      "rugosity_sc",
                                      "Rugosity (scaled)") )

# Market gravity — from habitat + market gravity model
( p_site_market <- plot_effect(s1_m_habitat_market,
                                    total_model_data,
                                    "log_market_gravity_sc",
                                    "Market gravity (scaled)") )

# Settlement gravity — from habitat + settlement gravity model
( p_site_sett <- plot_effect(s1_m_habitat_settgrav,
                                  total_model_data,
                                  "log_settlement_grav_sc",
                                  "Settlement gravity (scaled)") )
# ============================================================
#  TRANSECT-LEVEL BIOMASS (SENSITIVITY CHECK)
#
#  Rationale: Retains within-site variation. (1 | site) accounts
#  for non-independence of transects. 
#  Confirms site-level findings are not an artefact of collapsing to site means.
#
#  Response:   log(transect_total_biomass) — continuous, no zeros
#  Family:     Gaussian (log-transformed response)
#              Gamma also tested (Step 11)
#  Random fx:  (1 | site), then (1 | country/site)
# ============================================================

# ── Explore transect-level response ─────────────────
summary(total_transects$transect_total_biomass)

zeros <- mean(total_transects$transect_total_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")

# Raw distribution
( total_raw <- ggplot(total_transects, aes(x = transect_total_biomass)) +
    geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
    labs(x = "Total biomass per transect (g)", y = "Frequency",
         title = "Raw Total Biomass") +
    theme_bw() )

( total_log <- ggplot(transect_model_data, aes(x = log_transect_biomass)) +
    geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
    labs(x = "log(biomass)", y = "Frequency",
         title = "Log-transformed Total Biomass") +
    theme_bw() )

( total_sqrt <- ggplot(transect_model_data, aes(x = sqrt(transect_total_biomass))) +
    geom_histogram(bins = 30, fill = "#d7191c", colour = "white") +
    labs(x = "sqrt(biomass)", y = "Frequency",
         title = "Sqrt-transformed Total Biomass") +
    theme_bw() )

jpeg("total_biomass_distributions.jpg", width = 33, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(total_raw, total_log, total_sqrt, ncol = 3)
dev.off()

# ── Box-Cox ───────────────────────────────────────────────────
MASS::boxcox(
  lm(transect_total_biomass ~ 1, data = transect_model_data),
  lambda = seq(-2, 2, 0.1)
)
# lambda ~ 0 → log transformation appropriate

# ── Variation by site ─────────────────────────────────────────
ggplot(transect_model_data,
       aes(x = reorder(site, transect_total_biomass, median),
           y = transect_total_biomass)) +
  geom_boxplot(fill = "#2c7bb6", alpha = 0.6,
               outlier.colour = "black", outlier.size = 1) +
  coord_flip() +
  labs(x = NULL, y = "Total biomass (g)",
       title = "Biomass distribution by site") +
  theme_bw(base_size = 9)

# ── Family selection ──────────────────────────────────────────
# AICc not comparable between F1 and F2 (different response).
# Select on DHARMa diagnostics alone.

# F1: Gaussian on log(y)
mF1_gaussian <- glmmTMB(
  log_transect_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = gaussian(),
  data   = transect_model_data
)

resF1 <- simulateResiduals(mF1_gaussian, n = 1000)

jpeg("dharma_F1_gaussian_total.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF1, main = "DHARMa — Gaussian on log(y)"); dev.off()

plot(resF1)
testDispersion(resF1)
testOutliers(resF1)

# F2: Gamma (log link)
mF2_gamma <- glmmTMB(
  transect_total_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = Gamma(link = "log"),
  data   = transect_model_data
)

resF2 <- simulateResiduals(mF2_gamma, n = 1000)

jpeg("dharma_F2_gamma_total.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF2, main = "DHARMa — Gamma"); dev.off()

plot(resF2)
testDispersion(resF2)
testOutliers(resF2)

# ── Family selection decision ─────────────────────────────────
# Gaussian on log(y) chosen:
# Dispersion p = 0.994, outlier p = 0.110 — all tests n.s.

# ── Random effect structure selection ────────────────────────
# Anchor: full market gravity model — same fixed effects across
# all RE structures.

re_t_null    <- glmmTMB(log_transect_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + rugosity_sc,
                        family = gaussian(), data = transect_model_data)

re_t_site    <- glmmTMB(log_transect_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + rugosity_sc +
                          (1 | site),
                        family = gaussian(), data = transect_model_data)

re_t_nested  <- glmmTMB(log_transect_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + rugosity_sc +
                          (1 | country/site),
                        family = gaussian(), data = transect_model_data)

re_comparison_t <- list(
  "No RE"              = re_t_null,
  "(1 | site)"         = re_t_site,
  "(1 | country/site)" = re_t_nested
)

cat("\n--- RE structure comparison (transect level) ---\n")
print(make_aicc_df(re_comparison_t))

# ── Candidate models ──────────────────────────────────────────
# Family: Gaussian on log(y); random effect: (1 | site).

m0             <- glmmTMB(log_transect_biomass ~ 1                                                                         + (1 | site), family = gaussian(), data = transect_model_data)
m_env          <- glmmTMB(log_transect_biomass ~ sst_sc + log_chla_sc                                                     + (1 | site), family = gaussian(), data = transect_model_data)
m_market       <- glmmTMB(log_transect_biomass ~ log_market_gravity_sc                                                    + (1 | site), family = gaussian(), data = transect_model_data)
m_settgrav     <- glmmTMB(log_transect_biomass ~ log_settlement_grav_sc                                                   + (1 | site), family = gaussian(), data = transect_model_data)
m_hab          <- glmmTMB(log_transect_biomass ~ rugosity_sc                                                              + (1 | site), family = gaussian(), data = transect_model_data)
m_env_market   <- glmmTMB(log_transect_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc                            + (1 | site), family = gaussian(), data = transect_model_data)
m_env_settgrav <- glmmTMB(log_transect_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc                           + (1 | site), family = gaussian(), data = transect_model_data)
m_hab_market   <- glmmTMB(log_transect_biomass ~ rugosity_sc + log_market_gravity_sc                                     + (1 | site), family = gaussian(), data = transect_model_data)
m_hab_settgrav <- glmmTMB(log_transect_biomass ~ rugosity_sc + log_settlement_grav_sc                                    + (1 | site), family = gaussian(), data = transect_model_data)
m_full_market  <- glmmTMB(log_transect_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc             + (1 | site), family = gaussian(), data = transect_model_data)
m_full_settgrav <- glmmTMB(log_transect_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc            + (1 | site), family = gaussian(), data = transect_model_data)

model_list_transect <- list(
  "Null"                         = m0,
  "Environment"                  = m_env,
  "Market gravity"               = m_market,
  "Settlement gravity"           = m_settgrav,
  "Habitat"                      = m_hab,
  "Env + market gravity"         = m_env_market,
  "Env + settlement gravity"     = m_env_settgrav,
  "Habitat + market gravity"     = m_hab_market,
  "Habitat + settlement gravity" = m_hab_settgrav,
  "Full (market gravity)"        = m_full_market,
  "Full (settlement gravity)"    = m_full_settgrav
)

cat("\n--- AICc: Transect-level biomass ---\n")
print(make_aicc_df(model_list_transect))

# ── Diagnostics on best model ─────────────────────────────────
res_t <- simulateResiduals(m_hab, n = 1000) 

jpeg("dharma_transect_best.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_t, main = "DHARMa — transect-level biomass best model")
dev.off()

testDispersion(res_t)
testOutliers(res_t)
plotResiduals(res_t, transect_model_data$rugosity_sc, xlab = "Rugosity")

summary(m_hab)

# ── Marginal effect plot ──────────────────────────────────────
( p_rug <- plot_effect(m_hab, transect_model_data, "rugosity_sc", "Rugosity (scaled)") )

# ============================================================
#  PART 3 — TRANSECT-LEVEL COUNTS (COMPLEMENTARY ANALYSIS)
#
#  Rationale: Biomass is derived from counts via length-weight
#  relationships. Modelling the raw count directly is more
#  appropriate because:
#    (1) The data-generating process is discrete (whole fish)
#    (2) Count models correctly handle the mean-variance
#        relationship inherent in ecological count data
#    (3) Biomass transformations discard information about the
#        original counting process
#
#  Response: Total fish count per transect (integer >= 0)
#  Family:   Poisson → NB2 → NB1; selected via AICc
#  Random fx: (1 | site)
# ============================================================

# ── Explore count distribution ────────────────────────────────
cat("Transects (count data):", nrow(transect_model_data), "\n")
cat("Zeros in count data:",    sum(transect_model_data$transect_total_count == 0), "\n")
cat("Proportion zeros:",       round(mean(transect_model_data$transect_total_count == 0), 3), "\n")

summary(transect_model_data$transect_total_count)

ggplot(transect_model_data, aes(x = transect_total_count)) +
  geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
  labs(x = "Total fish count per transect", y = "Frequency",
       title = "Raw count distribution") +
  theme_bw()

# ── Mean-variance relationship ────────────────────────────────
# Points above the Poisson line → overdispersion → Negative Binomial.
transect_model_data %>%
  group_by(site) %>%
  summarise(mean_count = mean(transect_total_count),
            var_count  = var(transect_total_count),
            .groups    = "drop") %>%
  ggplot(aes(x = mean_count, y = var_count)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "red") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Site mean count", y = "Site variance",
       title = "Mean-variance relationship (red = Poisson expectation)") +
  theme_bw()

# ── Family selection ──────────────────────────────────────────
# AICc comparable between C1-C3 (same response, same link, same data).

# C1: Poisson
mC1_poisson <- glmmTMB(
  transect_total_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = poisson(link = "log"),
  data   = transect_model_data
)

resC1 <- simulateResiduals(mC1_poisson, n = 1000)

jpeg("dharma_C1_poisson_counts.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC1, main = "DHARMa — Poisson"); dev.off()

plot(resC1)
testDispersion(resC1)
testZeroInflation(resC1)
testOutliers(resC1)

# C2: NB2 — quadratic variance (classic NB)
mC2_nb2 <- glmmTMB(
  transect_total_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = nbinom2(link = "log"),
  data   = transect_model_data
)

resC2 <- simulateResiduals(mC2_nb2, n = 1000)

jpeg("dharma_C2_nb2_counts.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC2, main = "DHARMa — NB2"); dev.off()

plot(resC2)
testDispersion(resC2)
testZeroInflation(resC2)
testOutliers(resC2)

# C3: NB1 — linear variance
mC3_nb1 <- glmmTMB(
  transect_total_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = nbinom1(link = "log"),
  data   = transect_model_data
)

resC3 <- simulateResiduals(mC3_nb1, n = 1000)

jpeg("dharma_C3_nb1_counts.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC3, main = "DHARMa — NB1"); dev.off()

plot(resC3)
testDispersion(resC3)
testZeroInflation(resC3)
testOutliers(resC3)

cat("\n--- Family selection: count models ---\n")
print(make_aicc_df(list(
  "Poisson" = mC1_poisson,
  "NB2"     = mC2_nb2,
  "NB1"     = mC3_nb1
)))
# NB2 selected: lowest AICc.

# ── Random effect structure selection ────────────────────────
# Anchor: full market gravity model — same fixed effects across
# all RE structures.

re_c_null   <- glmmTMB(transect_total_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + rugosity_sc,
                       family = nbinom2(link = "log"), data = transect_model_data)

re_c_site   <- glmmTMB(transect_total_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + rugosity_sc +
                         (1 | site),
                       family = nbinom2(link = "log"), data = transect_model_data)

re_c_nested <- glmmTMB(transect_total_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + rugosity_sc +
                         (1 | country/site),
                       family = nbinom2(link = "log"), data = transect_model_data)

re_comparison_c <- list(
  "No RE"              = re_c_null,
  "(1 | site)"         = re_c_site,
  "(1 | country/site)" = re_c_nested
)

cat("\n--- RE structure comparison (count level) ---\n")
print(make_aicc_df(re_comparison_c))

# ── Candidate models ──────────────────────────────────────────
# Family: NB2; random effect: (1 | site).

count_family <- nbinom2(link = "log")

cm0             <- glmmTMB(transect_total_count ~ 1                                                                        + (1 | site), family = count_family, data = transect_model_data)
cm_env          <- glmmTMB(transect_total_count ~ sst_sc + log_chla_sc                                                    + (1 | site), family = count_family, data = transect_model_data)
cm_market       <- glmmTMB(transect_total_count ~ log_market_gravity_sc                                                   + (1 | site), family = count_family, data = transect_model_data)
cm_settgrav     <- glmmTMB(transect_total_count ~ log_settlement_grav_sc                                                  + (1 | site), family = count_family, data = transect_model_data)
cm_hab          <- glmmTMB(transect_total_count ~ rugosity_sc                                                             + (1 | site), family = count_family, data = transect_model_data)
cm_env_mkt      <- glmmTMB(transect_total_count ~ sst_sc + log_chla_sc + log_market_gravity_sc                           + (1 | site), family = count_family, data = transect_model_data)
cm_env_sett     <- glmmTMB(transect_total_count ~ sst_sc + log_chla_sc + log_settlement_grav_sc                          + (1 | site), family = count_family, data = transect_model_data)
cm_hab_market   <- glmmTMB(transect_total_count ~ rugosity_sc + log_market_gravity_sc                                    + (1 | site), family = count_family, data = transect_model_data)
cm_hab_settgrav <- glmmTMB(transect_total_count ~ rugosity_sc + log_settlement_grav_sc                                   + (1 | site), family = count_family, data = transect_model_data)
cm_full_mkt     <- glmmTMB(transect_total_count ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc            + (1 | site), family = count_family, data = transect_model_data)
cm_full_sett    <- glmmTMB(transect_total_count ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc            + (1 | site), family = count_family, data = transect_model_data)

model_list_counts <- list(
  "Null"                         = cm0,
  "Environment"                  = cm_env,
  "Market gravity"               = cm_market,
  "Settlement gravity"           = cm_settgrav,
  "Habitat"                      = cm_hab,
  "Env + market gravity"         = cm_env_mkt,
  "Env + settlement gravity"     = cm_env_sett,
  "Habitat + market gravity"     = cm_hab_market,
  "Habitat + settlement gravity" = cm_hab_settgrav,
  "Full (market gravity)"        = cm_full_mkt,
  "Full (settlement gravity)"    = cm_full_sett
)

cat("\n--- AICc: count models ---\n")
print(make_aicc_df(model_list_counts))

# ── Diagnostics on best model ─────────────────────────────────
res_cm <- simulateResiduals(cm_full_sett, n = 1000)

jpeg("dharma_best_count_model.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_cm, main = "DHARMa — best count model"); dev.off()

plot(res_cm)
testDispersion(res_cm)
testZeroInflation(res_cm)
testOutliers(res_cm)
plotResiduals(res_cm, transect_model_data$rugosity_sc,            xlab = "Rugosity")
plotResiduals(res_cm, transect_model_data$sst_sc,                 xlab = "SST")
plotResiduals(res_cm, transect_model_data$log_chla_sc,            xlab = "Chl-a")
plotResiduals(res_cm, transect_model_data$log_settlement_grav_sc, xlab = "Settlement gravity")

summary(cm_full_sett)
exp(fixef(cm_full_sett)$cond)   # incidence rate ratios

# ── Marginal effect plots ─────────────────────────────────────
p_count_rug  <- plot_effect(cm_full_sett, transect_model_data, "rugosity_sc",
                            "Rugosity (scaled)",           "Expected fish count",
                            colour = "#d7191c")
p_count_sst  <- plot_effect(cm_full_sett, transect_model_data, "sst_sc",
                            "SST (scaled)",                "Expected fish count",
                            colour = "#d7191c")
p_count_chla <- plot_effect(cm_full_sett, transect_model_data, "log_chla_sc",
                            "Chl-a (scaled)",              "Expected fish count",
                            colour = "#d7191c")
p_count_grav <- plot_effect(cm_full_sett, transect_model_data, "log_settlement_grav_sc",
                            "Settlement gravity (scaled)", "Expected fish count",
                            colour = "#d7191c")

jpeg("count_marginal_effects.jpg", width = 33, height = 22, units = "cm", res = 300)
gridExtra::grid.arrange(p_count_rug, p_count_sst,
                        p_count_chla, p_count_grav, ncol = 2)
dev.off()

# ============================================================
#  SYNTHESIS: BIOMASS vs COUNT MODEL CONCLUSIONS
#
#  Part 1 best model:   s1_m_hab (lm, Gaussian on log-mean-biomass)
#                       — Habitat only (lowest AICc)
#  Part 2 best model:   m_hab (glmmTMB, Gaussian on log-transect-biomass)
#                       — Habitat only (lowest AICc)
#  Part 3 best model:   cm_full_sett (glmmTMB, NB2)
#                       — Full model with settlement gravity
#
#  Convergence across Parts 1 and 2 confirms site-level findings
#  are not an artefact of averaging.
#  Divergence between biomass and count models suggests human
#  pressure and chlorophyll operate through abundance rather
#  than body size.
# ============================================================

cat("\n=== PART 1 — Site-level biomass, best model ===\n")
summary(s1_m_hab)

cat("\n=== PART 2 — Transect-level biomass, best model ===\n")
summary(m_hab)

cat("\n=== PART 3 — Count model, best model ===\n")
summary(cm_full_sett)

# ── Rugosity: consistent across all three parts ───────────────
# Part 1: plain lm() — use coef()
# Parts 2 & 3: glmmTMB — use fixef()$cond
rug_site    <- coef(s1_m_hab)["rugosity_sc"]
rug_biomass <- fixef(m_hab)$cond["rugosity_sc"]
rug_count   <- fixef(cm_full_sett)$cond["rugosity_sc"]

cat("\n--- Rugosity effect ---\n")
cat("Part 1 — site biomass      beta:", round(rug_site,    3), "\n")
cat("Part 2 — transect biomass  beta:", round(rug_biomass, 3), "\n")
cat("Part 3 — count model       beta:", round(rug_count,   3), "\n")
cat("Part 3 — count model        IRR:", round(exp(rug_count), 3), "\n")
# Biomass beta > count beta: rugosity effect on biomass is
# disproportionately large relative to abundance, suggesting
# complex reefs support larger-bodied fish, not just more fish.

# ── Human pressure: appears in counts, absent from biomass ────
cat("\n--- Human pressure (settlement gravity) ---\n")
cat("Part 1 biomass: not in best model\n")
cat("Part 2 biomass: not in best model\n")
cat("Part 3 counts:  IRR =",
    round(exp(fixef(cm_full_sett)$cond["log_settlement_grav_sc"]), 3),
    "(p = 0.046)\n")
# Consistent non-detection in both biomass analyses but detection
# in counts is compatible with size-selective harvesting: large
# individuals removed, small fish largely remain intact.

# ── Chlorophyll-a: suppresses abundance, not biomass ──────────
cat("\n--- Chlorophyll-a ---\n")
cat("Part 3 counts  beta:", round(fixef(cm_full_sett)$cond["log_chla_sc"], 3),
    " IRR:", round(exp(fixef(cm_full_sett)$cond["log_chla_sc"]), 3), "\n")
cat("Parts 1 & 2 biomass: not in best model\n")
# High chl-a (turbidity/runoff) reduces fish numbers but does
# not detectably alter total biomass — individual fish in
# turbid sites may be larger on average, compensating.

# ── SST: weak positive trend in counts only ───────────────────
cat("\n--- SST ---\n")
cat("Part 3 counts  beta:", round(fixef(cm_full_sett)$cond["sst_sc"], 3),
    " IRR:", round(exp(fixef(cm_full_sett)$cond["sst_sc"]), 3),
    " (p = 0.094)\n")
cat("Parts 1 & 2 biomass: not in best model\n")
# Weak, non-significant positive trend — warmer sites may support
# marginally higher abundance but effect is not robust.

