# ============================================================
#  LARGE EXCAVATOR BIOMASS & ABUNDANCE — MIXED EFFECTS MODELS
#
#  Predictors measured at site level (averaged from station):
#    SST, Chl-a, Human gravity (market / settlement), Rugosity
#
#  Analytical structure:
#    PART 1 — Site-level analysis (PRIMARY)
#              Matches response resolution to predictor resolution.
#              Sites are the true unit of environmental inference.
#              Sections: no RE → country RE
#
#    PART 2 — Transect-level biomass (SENSITIVITY CHECK)
#              Retains within-site variation; (1 | site) accounts
#              for non-independence. Confirms site-level findings
#              are not an artefact of averaging.
#              NOTE: zero inflation may require Tweedie family.
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

# ── LOAD ALL DATA ─────────────────────────────────────────────
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
plot_effect <- function(model, data, focal_var,
                        x_label,
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
                fill = colour, alpha = 0.2) +
    geom_line(aes(y = fit), colour = colour, linewidth = 1.1) +
    labs(x = x_label, y = y_label) +
    theme_bw(base_size = 13) +
    theme(axis.title = element_text(face = "bold"))
}

# ==============================================================================
#  AGGREGATE TRANSECT DATA
#  Minimum of 3 transects per site.
# ==============================================================================

excavator_transects <- fish_2009 %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_excavator_biomass = sum(
      ifelse(trophic_group %in% c("large_excavators"), tot_wt_g, 0),
      na.rm = TRUE
    ),
    transect_excavator_count = sum(
      ifelse(trophic_group %in% c("large_excavators"), number, 0),
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

cat("Number of transects:", nrow(excavator_transects), "\n")
cat("Number of sites:",     n_distinct(excavator_transects$site), "\n")
cat("Number of countries:", n_distinct(excavator_transects$country), "\n")

# ==============================================================================
#  BIOMASS DATA EXPLORATION
# ==============================================================================

# ── Aggregate data at site level (mean large excavator biomass per site) ──────
site_data <- excavator_transects %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass = mean(transect_excavator_biomass, na.rm = TRUE),
    n_transects  = n(),
    .groups      = "drop"
  ) %>%
  mutate(
    site    = as.factor(site),
    country = as.factor(country)
  )

# ── Basic summary ─────────────────────────────────────────────────────────────
summary(site_data$mean_biomass)

zeros <- mean(site_data$mean_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")
# NOTE: site-level means may still be zero if all transects at a site
# recorded zero large excavators. Check whether these sites should be retained.

# ── Raw distribution ──────────────────────────────────────────────────────────
( site_raw <- ggplot(site_data, aes(x = mean_biomass)) +
    geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
    labs(x = "Mean large excavator biomass per site (g)", y = "Frequency",
         title = "Raw Site-Level Large Excavator Biomass") +
    theme_bw() )

# ── Box-Cox: what power transformation does the data suggest? ─────────────────
# Fit on non-zero values only (boxcox requires y > 0)
site_nonzero <- site_data %>% filter(mean_biomass > 0)

MASS::boxcox(
  lm(mean_biomass ~ 1, data = site_nonzero),
  lambda = seq(-2, 2, 0.1)
)
# lambda ~ 0 → log transformation appropriate

# ── Apply transformations ─────────────────────────────────────────────────────
# Small constant (+0.01) handles zero site means; check sensitivity if needed
site_data <- site_data %>%
  mutate(
    log_mean_biomass  = log(mean_biomass + 0.01),
    sqrt_mean_biomass = sqrt(mean_biomass)
  )

# ── Log and sqrt distributions ────────────────────────────────────────────────
( site_log <- ggplot(site_data, aes(x = log_mean_biomass)) +
    geom_histogram(bins = 25, fill = "#1a9641", colour = "white") +
    labs(x = "log(mean biomass + 0.01)", y = "Frequency",
         title = "Log-transformed Site-Level Large Excavator Biomass") +
    theme_bw() )

( site_sqrt <- ggplot(site_data, aes(x = sqrt_mean_biomass)) +
    geom_histogram(bins = 25, fill = "#d7191c", colour = "white") +
    labs(x = "sqrt(mean biomass)", y = "Frequency",
         title = "Sqrt-transformed Site-Level Large Excavator Biomass") +
    theme_bw() )

jpeg("site_excavator_biomass_distributions.jpg", width = 33, height = 11,
     units = "cm", res = 300)
gridExtra::grid.arrange(site_raw, site_log, site_sqrt, ncol = 3)
dev.off()

# ── Normality checks on log-transformed response ──────────────────────────────
qqnorm(site_data$log_mean_biomass,
       main = "Q-Q plot: log(mean large excavator biomass per site)")
qqline(site_data$log_mean_biomass, col = "red")
shapiro.test(site_data$log_mean_biomass)
# Shapiro-Wilk: W = 0.897, p < 0.001 — notable departure from
# normality on log scale, driven by the 10 zero-mean sites.
# Tweedie on raw scale (F2) is preferred precisely because it
# handles these zeros natively without transformation artefacts.

# ── Variation by site ─────────────────────────────────────────────────────────
ggplot(site_data, aes(x = reorder(site, mean_biomass, median),
                      y = mean_biomass)) +
  geom_col(fill = "#2c7bb6", alpha = 0.7) +
  coord_flip() +
  labs(x = NULL, y = "Mean large excavator biomass (g)",
       title = "Mean large excavator biomass by site (raw)") +
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
  labs(x = NULL, y = "log(mean large excavator biomass per site + 0.01)") +
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

# ── Human gravity metrics ─────────────────────────────────────
gravity_sites <- gravity_2009 %>%
  st_drop_geometry() %>%
  group_by(site) %>%
  summarise(
    market_gravity  = mean(market_grav,         na.rm = TRUE),
    settlement_pop  = mean(settlement_tot_pop,  na.rm = TRUE),
    settlement_grav = mean(nearest_pop75_grav,  na.rm = TRUE),
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

# ── Inspect raw predictors ────────────────────────────────────
raw_predictors <- gravity_sites %>%
  left_join(chla_sites,     by = "site") %>%
  left_join(sst_sites,      by = "site") %>%
  left_join(rugosity_sites, by = "site")

predictor_labels    <- c("Market gravity", "Settlement gravity", "Settlement pop.",
                         "Chlorophyll-a", "SST", "Rugosity")
predictor_order_raw <- c("market_gravity", "settlement_grav", "settlement_pop",
                         "mean_annual_chla", "mean_annual_sst", "rugosity")

( p_pred_raw <- raw_predictors %>%
    dplyr::select(all_of(predictor_order_raw)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, predictor_order_raw, predictor_labels)) %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
    facet_wrap(~ variable, scales = "free") +
    labs(title = "Raw predictors") + theme_bw() )

# ── Inspect transformed predictors ────────────────────────────
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

( p_pred_tran <- transformed_predictors %>%
    dplyr::select(all_of(predictor_order_tran)) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    mutate(variable = factor(variable, predictor_order_tran, predictor_labels)) %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
    facet_wrap(~ variable, scales = "free") +
    labs(title = "Transformed predictors") + theme_bw() )

jpeg("predictor_distributions_excavator.jpg", width = 33, height = 22, units = "cm", res = 300)
gridExtra::grid.arrange(p_pred_raw, p_pred_tran, nrow = 2)
dev.off()

# Combine and scale predictors for analyses
scaled_predictors <- transformed_predictors %>%
  transmute(
    site                   = site,
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

# ============================================================
#  CHOOSING SETTLEMENT METRIC
# ============================================================
# Settlement gravity and settlement pop both proxy local human
# pressure. Select the better-performing metric via AICc before
# entering the main candidate set.
# NOTE: zero-containing site means require glmmTMB Tweedie here.

settlement_data <- excavator_transects %>%
  left_join(scaled_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    mean_excavator_biomass  = mean(transect_excavator_biomass, na.rm = TRUE),
    log_settlement_grav_sc  = first(log_settlement_grav_sc),
    log_settlement_pop_sc   = first(log_settlement_pop_sc),
    .groups = "drop"
  )

settgrav <- glmmTMB(mean_excavator_biomass ~ log_settlement_grav_sc,
                    family = tweedie(link = "log"), data = settlement_data)
settpop  <- glmmTMB(mean_excavator_biomass ~ log_settlement_pop_sc,
                    family = tweedie(link = "log"), data = settlement_data)

make_aicc_df(list(
  "Settlement gravity" = settgrav,
  "Settlement pop."    = settpop
))

# Settlement metric selection: virtually identical performance
# (delta AICc = 0.11, weights 0.51 vs 0.49). Neither metric
# has any meaningful advantage over the other. Both carried
# forward as parallel candidate model sets throughout all
# subsequent analyses. Results compared across both metrics
# to assess robustness of human pressure effects.

rm(settlement_data)

# ============================================================
#  ANALYSIS DATASETS
# ============================================================

# Retain all human pressure metrics for parallel model sets
final_predictors <- scaled_predictors %>%
  dplyr::select(site, log_market_gravity_sc, log_settlement_grav_sc,
                log_settlement_pop_sc, log_chla_sc, sst_sc, rugosity_sc)

# ── Transect-level dataset ────────────────────────────────────
transect_model_data <- excavator_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect model data:", nrow(transect_model_data), "rows,",
    n_distinct(transect_model_data$site), "sites\n")
cat("Biomass zeros:", sum(transect_model_data$transect_excavator_biomass == 0), "\n")
cat("Count zeros:",  sum(transect_model_data$transect_excavator_count   == 0), "\n")

# ── Site-level dataset ────────────────────────────────────────
total_model_data <- transect_model_data %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_excavator_biomass, na.rm = TRUE),
    n_transects            = n(),
    log_market_gravity_sc  = first(log_market_gravity_sc),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_settlement_pop_sc  = first(log_settlement_pop_sc),
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
cat("Site-level zeros:", sum(total_model_data$mean_biomass == 0), "\n")

# ============================================================
#  PART 1 — SITE-LEVEL ANALYSIS (PRIMARY)
# ============================================================

# ── FAMILY SELECTION ──────────────────────────────────────────
# Work from simplest to most complex.
# AICc only comparable between models with the same response.
# Gaussian on log(y + 0.01) vs Tweedie vs ZI-Tweedie.

# ── F1: Gaussian on log-transformed mean biomass ──────────────
# Small constant handles zero site means.
total_model_data <- total_model_data %>%
  mutate(log_mean_biomass = log(mean_biomass + 0.01))

mS_F1 <- glmmTMB(log_mean_biomass ~ sst_sc + log_chla_sc +
                   log_market_gravity_sc + rugosity_sc,
                 family = gaussian(), data = total_model_data)

resS_F1 <- simulateResiduals(mS_F1, n = 1000)

jpeg("diagnostics_site_excavator_F1_gaussian_log.jpg",
     width = 25, height = 15, units = "cm", res = 300)
plot(resS_F1, main = "DHARMa — Gaussian on log(y + 0.01)"); dev.off()

plot(resS_F1)
testZeroInflation(resS_F1)
testDispersion(resS_F1)

# ── F2: Tweedie (log link) on raw mean biomass ────────────────
# Compound Poisson-Gamma: natively produces exact zeros alongside
# a continuous positive distribution. Preferred when zeros are
# present at the site level.
mS_F2 <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc +
                   log_market_gravity_sc + rugosity_sc,
                 family = tweedie(link = "log"), data = total_model_data)

resS_F2 <- simulateResiduals(mS_F2, n = 1000)

jpeg("diagnostics_site_excavator_F2_tweedie.jpg",
     width = 25, height = 15, units = "cm", res = 300)
plot(resS_F2, main = "DHARMa — Tweedie"); dev.off()

plot(resS_F2)
testZeroInflation(resS_F2)
testDispersion(resS_F2)


# ── Family selection decision ─────────────────────────────────
# Both F1 and F2 show adequate DHARMa diagnostics (64 sites, 10 zeros).
# Zero inflation test n.s. for both (F1: NaN p = 1.00; F2: p = 0.884).
# Dispersion test n.s. for both (F1: p = 0.914; F2: p = 0.804).
# Tweedie (F2) preferred: handles the 10 site-level zeros natively
# without a log + constant approximation and operates on the raw
# biomass scale. ZI Tweedie not attempted — no evidence of excess
# zeros in F2 and site-level n too small to support the ZI component.

# ── RANDOM EFFECT STRUCTURE SELECTION ────────────────────────
# Anchor: full market gravity model.
# glmmTMB allows consistent comparison across RE structures.

re_null <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc +
                     log_market_gravity_sc + rugosity_sc,
                   family = tweedie(link = "log"), data = total_model_data)

re_country <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc +
                        log_market_gravity_sc + rugosity_sc +
                        (1 | country),
                      family = tweedie(link = "log"), data = total_model_data)

cat("\n--- RE structure comparison (site-level large excavator) ---\n")
print(make_aicc_df(list(
  "No RE"         = re_null,
  "(1 | country)" = re_country
)))

# ── RE structure decision ─────────────────────────────────────
# No RE clearly preferred (AICc = 785.34 vs 787.96 for country RE;
# delta = 2.62, weight = 0.79 vs 0.21). Country-level variation is
# adequately absorbed by the included predictors. No RE carried
# forward for all site-level candidate models

# ── CANDIDATE MODELS — SITE LEVEL ────────────────────────────

# ── No random effects ─────────────────────────────────────────
# --- Null ---
s1_m0                 <- glmmTMB(mean_biomass ~ 1,
                                 family = tweedie(link = "log"), data = total_model_data)
# --- Single predictor ---
s1_m_env              <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
s1_m_market           <- glmmTMB(mean_biomass ~ log_market_gravity_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
s1_m_settgrav         <- glmmTMB(mean_biomass ~ log_settlement_grav_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
s1_m_settpop          <- glmmTMB(mean_biomass ~ log_settlement_pop_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
s1_m_hab              <- glmmTMB(mean_biomass ~ rugosity_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
# --- Environment + human pressure ---
s1_m_env_mkt          <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
s1_m_env_settgrav     <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
s1_m_env_settpop      <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
# --- Habitat + human pressure ---
s1_m_hab_market       <- glmmTMB(mean_biomass ~ rugosity_sc + log_market_gravity_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
s1_m_hab_settgrav     <- glmmTMB(mean_biomass ~ rugosity_sc + log_settlement_grav_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
s1_m_hab_settpop      <- glmmTMB(mean_biomass ~ rugosity_sc + log_settlement_pop_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
# --- Full (single human pressure metric) ---
s1_m_full_mkt         <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
s1_m_full_settgrav    <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
s1_m_full_settpop     <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc  + rugosity_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
# --- Combined gravity metrics ---
s1_m_both_grav        <- glmmTMB(mean_biomass ~ log_market_gravity_sc + log_settlement_grav_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
s1_m_hab_both_grav    <- glmmTMB(mean_biomass ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
s1_m_env_both_grav    <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc,
                                 family = tweedie(link = "log"), data = total_model_data)
s1_m_full_both_grav   <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc,
                                 family = tweedie(link = "log"), data = total_model_data)

model_list_s1 <- list(
  # Null
  "Null"                               = s1_m0,
  # Single predictor
  "Environment"                        = s1_m_env,
  "Market gravity"                     = s1_m_market,
  "Settlement gravity"                 = s1_m_settgrav,
  "Settlement pop."                    = s1_m_settpop,
  "Habitat"                            = s1_m_hab,
  # Environment + human pressure
  "Env + market gravity"               = s1_m_env_mkt,
  "Env + settlement gravity"           = s1_m_env_settgrav,
  "Env + settlement pop."              = s1_m_env_settpop,
  # Habitat + human pressure
  "Habitat + market gravity"           = s1_m_hab_market,
  "Habitat + settlement gravity"       = s1_m_hab_settgrav,
  "Habitat + settlement pop."          = s1_m_hab_settpop,
  # Full - single human pressure
  "Full (market gravity)"              = s1_m_full_mkt,
  "Full (settlement gravity)"          = s1_m_full_settgrav,
  "Full (settlement pop.)"             = s1_m_full_settpop,
  # Combined gravity
  "Both gravity"                       = s1_m_both_grav,
  "Habitat + both gravity"             = s1_m_hab_both_grav,
  "Env + both gravity"                 = s1_m_env_both_grav,
  "Full (both gravity)"                = s1_m_full_both_grav
)

cat("\n--- AICc: Site-level large excavator candidate models (no RE) ---\n")
aicc_site_s1 <- make_aicc_df(model_list_s1)
print(aicc_site_s1)


# ── Residual diagnostics — top models ────────────────────────
# Top model: Habitat only (weight = 0.3332).
# Habitat + settlement pop. is closest competitor (ΔAICc = 0.95).
# Habitat + market gravity also competitive (ΔAICc = 1.27).
# All three habitat-containing models within ΔAICc < 2.
# Both diagnostics models show clean residuals: no dispersion
# issues, no zero inflation, no outliers.
# Settlement pop. coefficient not significant (beta = 0.276, p = 0.244),
# confirming habitat only as the most parsimonious description.

cat("\n--- Diagnostics: Habitat only ---\n")
res_s1_hab <- simulateResiduals(s1_m_hab, n = 1000)

jpeg("diagnostics_site_excavator_hab.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_s1_hab, main = "DHARMa — Habitat only"); dev.off()

plot(res_s1_hab)
testZeroInflation(res_s1_hab)
testDispersion(res_s1_hab)
testOutliers(res_s1_hab)

cat("\n--- Diagnostics: Habitat + settlement pop.---\n")
res_s1_hab_settpop <- simulateResiduals(s1_m_hab_settpop, n = 1000)

jpeg("diagnostics_site_excavator_hab_settpop.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_s1_hab_settpop, main = "DHARMa — Habitat + settlement pop."); dev.off()

plot(res_s1_hab_settpop)
testZeroInflation(res_s1_hab_settpop)
testDispersion(res_s1_hab_settpop)
testOutliers(res_s1_hab_settpop)

# ── Summaries ─────────────────────────────────────────────────
cat("\n--- Summary: Habitat only ---\n")
summary(s1_m_hab)

cat("\n--- Summary: Habitat + settlement pop. ---\n")
summary(s1_m_hab_settpop)

# ── Rugosity coefficient stability across model structures ────
cat("\n--- Rugosity coefficient stability ---\n")
cat("Habitat only:           beta =", round(fixef(s1_m_hab)$cond["rugosity_sc"],        3), "\n")
cat("Habitat + market:       beta =", round(fixef(s1_m_hab_market)$cond["rugosity_sc"],  3), "\n")
cat("Habitat + settgrav:     beta =", round(fixef(s1_m_hab_settgrav)$cond["rugosity_sc"],3), "\n")
cat("Habitat + settpop:      beta =", round(fixef(s1_m_hab_settpop)$cond["rugosity_sc"], 3), "\n")

# ── Marginal effect plots ─────────────────────────────────────
# Rugosity from habitat-only model (top-ranked)
( p_site_rugosity <- plot_effect(s1_m_hab,
                                 total_model_data,
                                 "rugosity_sc",
                                 "Rugosity (scaled)",
                                 y_label = "Large excavator biomass (g)") )

# Settlement pop. from habitat + settlement pop. model
( p_site_settpop <- plot_effect(s1_m_hab_settpop,
                               total_model_data,
                               "log_settlement_pop_sc",
                               "Log Settlement pop (scaled)",
                               y_label = "Large excavator biomass (g)") )

jpeg("site_excavator_marginal_effects.jpg", width = 22, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(p_site_rugosity, p_site_settpop, ncol = 2)
dev.off()

# ============================================================
#  PART 2 — TRANSECT-LEVEL BIOMASS (SENSITIVITY CHECK)
#
#  Rationale: Retains within-site variation. (1 | site) accounts
#  for non-independence of transects.
#  Confirms site-level findings are not an artefact of collapsing
#  to site means.
#
#  Response:   transect_excavator_biomass — continuous, zero-inflated (~57.8% zeros at transect level)
#  Family:     Tweedie (log link) — confirm via family selection below
#  Random fx:  (1 | site), then (1 | country/site)
# ============================================================

# ── Explore transect-level response ──────────────────────────
summary(excavator_transects$transect_excavator_biomass)

zeros <- mean(excavator_transects$transect_excavator_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")

( excav_raw <- ggplot(excavator_transects, aes(x = transect_excavator_biomass)) +
    geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
    labs(x = "Large excavator biomass per transect (g)", y = "Frequency",
         title = "Raw Large Excavator Biomass") +
    theme_bw() )

excavator_transects <- excavator_transects %>%
  mutate(
    log_excavator_biomass  = log(transect_excavator_biomass + 0.01),
    sqrt_excavator_biomass = sqrt(transect_excavator_biomass)
  )

( excav_log <- ggplot(excavator_transects, aes(x = log_excavator_biomass)) +
    geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
    labs(x = "log(biomass + 0.01)", y = "Frequency",
         title = "Log-transformed Large Excavator Biomass") +
    theme_bw() )

( excav_sqrt <- ggplot(excavator_transects, aes(x = sqrt_excavator_biomass)) +
    geom_histogram(bins = 30, fill = "#d7191c", colour = "white") +
    labs(x = "sqrt(biomass)", y = "Frequency",
         title = "Sqrt-transformed Large Excavator Biomass") +
    theme_bw() )

jpeg("excavator_biomass_distributions.jpg", width = 33, height = 11,
     units = "cm", res = 300)
gridExtra::grid.arrange(excav_raw, excav_log, excav_sqrt, ncol = 3)
dev.off()

# Add log excavator biomass to transect model data
transect_model_data <- transect_model_data %>%
  mutate(log_excavator_biomass = log(transect_excavator_biomass + 0.01))

# ── Box-Cox on non-zero values ────────────────────────────────
excav_nonzero <- excavator_transects %>% filter(transect_excavator_biomass > 0)

MASS::boxcox(
  lm(transect_excavator_biomass ~ 1, data = excav_nonzero),
  lambda = seq(-2, 2, 0.1)
)
# lambda ~ 0 → log transformation appropriate for positive values

# ── Variation by site ─────────────────────────────────────────
ggplot(transect_model_data,
       aes(x = reorder(site, transect_excavator_biomass, median),
           y = transect_excavator_biomass)) +
  geom_boxplot(fill = "#2c7bb6", alpha = 0.6,
               outlier.colour = "black", outlier.size = 1) +
  coord_flip() +
  labs(x = NULL, y = "Large excavator biomass (g)",
       title = "Large excavator biomass distribution by site") +
  theme_bw(base_size = 9)

# ── Zeros by site ─────────────────────────────────────────────
transect_model_data %>%
  group_by(site) %>%
  summarise(
    n_transects  = n(),
    prop_zeros   = mean(transect_excavator_biomass == 0),
    mean_biomass = mean(transect_excavator_biomass),
    .groups = "drop"
  ) %>%
  arrange(desc(prop_zeros)) %>%
  print(n = Inf)

# ── Family selection ──────────────────────────────────────────
# AICc not comparable between F1 and F2/F3 (different response).
# Select on DHARMa diagnostics; use AICc only to compare F2 vs F3.

# F1: Gaussian on log(y + 0.01)
mF1_gaussian <- glmmTMB(
  log_excavator_biomass ~ sst_sc + log_chla_sc +
    log_settlement_pop_sc + rugosity_sc + (1 | site),
  family = gaussian(),
  data   = transect_model_data
)

resF1 <- simulateResiduals(mF1_gaussian, n = 1000)

jpeg("dharma_excavator_F1_gaussian.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF1, main = "DHARMa — Gaussian on log(y + 0.01)"); dev.off()

plot(resF1)
testZeroInflation(resF1)
testDispersion(resF1)

# F2: Plain Tweedie
mF2_tweedie <- glmmTMB(
  transect_excavator_biomass ~ sst_sc + log_chla_sc +
    log_settlement_pop_sc + rugosity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = transect_model_data
)

resF2 <- simulateResiduals(mF2_tweedie, n = 1000)

jpeg("dharma_excavator_F2_tweedie.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF2, main = "DHARMa — Tweedie"); dev.off()

plot(resF2)
testZeroInflation(resF2)
testDispersion(resF2)

# F3: Zero-inflated Tweedie
mF3_tweedie_zi <- glmmTMB(
  transect_excavator_biomass ~ sst_sc + log_chla_sc +
    log_settlement_pop_sc + rugosity_sc + (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = transect_model_data
)

resF3 <- simulateResiduals(mF3_tweedie_zi, n = 1000)

jpeg("dharma_excavator_F3_tweedie_zi.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF3, main = "DHARMa — Zero-inflated Tweedie"); dev.off()

plot(resF3)
testZeroInflation(resF3)
testDispersion(resF3)

cat("\n--- Family selection: transect-level large excavator biomass ---\n")
print(make_aicc_df(list(
  "Tweedie"    = mF2_tweedie,
  "ZI Tweedie" = mF3_tweedie_zi
)))

# ── Family selection decision ─────────────────────────────────
# Plain Tweedie (F2) retained.
# ZI Tweedie (F3) does not meet the prespecified adoption threshold:
#   - ΔAICc = 2.13 (ZI Tweedie vs plain Tweedie) — does not exceed
#     the >2 threshold
#   - Zero inflation test n.s. for both F2 (p = 0.662) and
#     F3 (p = 0.632) — no evidence of excess zeros
#
# Both conditions for adopting ZI Tweedie must be met:
#   (1) AICc must improve by > 2, AND
#   (2) Zero inflation test must be significant
# Neither condition is met.
#
# F1 (Gaussian on log scale) also shows clean diagnostics
# (dispersion p = 0.972) but requires an arbitrary +0.01 constant
# for the 57.8% transect-level zeros.
#
# Plain Tweedie carried forward — handles zeros natively and
# shows no dispersion or zero inflation issues.

# ── Random effect structure selection ────────────────────────
# Anchor: full market gravity model.

re_t_null   <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc +
                         log_settlement_pop_sc + rugosity_sc,
                       family = tweedie(link = "log"), data = transect_model_data)

re_t_site   <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc +
                         log_settlement_pop_sc + rugosity_sc +
                         (1 | site),
                       family = tweedie(link = "log"), data = transect_model_data)

re_t_nested <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc +
                         log_settlement_pop_sc + rugosity_sc +
                         (1 | country/site),
                       family = tweedie(link = "log"), data = transect_model_data)

cat("\n--- RE structure comparison (transect-level large excavator) ---\n")
print(make_aicc_df(list(
  "No RE"              = re_t_null,
  "(1 | site)"         = re_t_site,
  "(1 | country/site)" = re_t_nested
)))

# ── Random effect structure decision ─────────────────────────
# (1 | site) clearly preferred (AICc = 2087.63; delta = 2.13 over
# (1 | country/site), weight = 0.74 vs 0.26). No RE is strongly
# rejected (delta = 30.41). Site-level random intercept carried
# forward — the added country-level nesting is not supported.

# ── Candidate models ──────────────────────────────────────────
# Family: Tweedie (log link). RE: (1 | site).

excav_family <- tweedie(link = "log")

# --- Null ---
m0                 <- glmmTMB(transect_excavator_biomass ~ 1                                                                                        + (1 | site), family = excav_family, data = transect_model_data)

# --- Single predictor ---
m_env              <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc                                                                    + (1 | site), family = excav_family, data = transect_model_data)
m_market           <- glmmTMB(transect_excavator_biomass ~ log_market_gravity_sc                                                                   + (1 | site), family = excav_family, data = transect_model_data)
m_settgrav         <- glmmTMB(transect_excavator_biomass ~ log_settlement_grav_sc                                                                  + (1 | site), family = excav_family, data = transect_model_data)
m_settpop          <- glmmTMB(transect_excavator_biomass ~ log_settlement_pop_sc                                                                   + (1 | site), family = excav_family, data = transect_model_data)
m_hab              <- glmmTMB(transect_excavator_biomass ~ rugosity_sc                                                                             + (1 | site), family = excav_family, data = transect_model_data)

# --- Environment + human pressure ---
m_env_market       <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc                                           + (1 | site), family = excav_family, data = transect_model_data)
m_env_settgrav     <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc                                          + (1 | site), family = excav_family, data = transect_model_data)
m_env_settpop      <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc                                           + (1 | site), family = excav_family, data = transect_model_data)

# --- Habitat + human pressure ---
m_hab_market       <- glmmTMB(transect_excavator_biomass ~ rugosity_sc + log_market_gravity_sc                                                     + (1 | site), family = excav_family, data = transect_model_data)
m_hab_settgrav     <- glmmTMB(transect_excavator_biomass ~ rugosity_sc + log_settlement_grav_sc                                                    + (1 | site), family = excav_family, data = transect_model_data)
m_hab_settpop      <- glmmTMB(transect_excavator_biomass ~ rugosity_sc + log_settlement_pop_sc                                                     + (1 | site), family = excav_family, data = transect_model_data)

# --- Full (single human pressure metric) ---
m_full_market      <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc                            + (1 | site), family = excav_family, data = transect_model_data)
m_full_settgrav    <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc                            + (1 | site), family = excav_family, data = transect_model_data)
m_full_settpop     <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc  + rugosity_sc                            + (1 | site), family = excav_family, data = transect_model_data)

# --- Combined gravity metrics ---
m_both_grav        <- glmmTMB(transect_excavator_biomass ~ log_market_gravity_sc + log_settlement_grav_sc                                         + (1 | site), family = excav_family, data = transect_model_data)
m_hab_both_grav    <- glmmTMB(transect_excavator_biomass ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc                           + (1 | site), family = excav_family, data = transect_model_data)
m_env_both_grav    <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc                  + (1 | site), family = excav_family, data = transect_model_data)
m_full_both_grav   <- glmmTMB(transect_excavator_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc    + (1 | site), family = excav_family, data = transect_model_data)

model_list_transect <- list(
  # Null
  "Null"                               = m0,
  # Single predictor
  "Environment"                        = m_env,
  "Market gravity"                     = m_market,
  "Settlement gravity"                 = m_settgrav,
  "Settlement pop."                    = m_settpop,
  "Habitat"                            = m_hab,
  # Environment + human pressure
  "Env + market gravity"               = m_env_market,
  "Env + settlement gravity"           = m_env_settgrav,
  "Env + settlement pop."              = m_env_settpop,
  # Habitat + human pressure
  "Habitat + market gravity"           = m_hab_market,
  "Habitat + settlement gravity"       = m_hab_settgrav,
  "Habitat + settlement pop."          = m_hab_settpop,
  # Full - single human pressure
  "Full (market gravity)"              = m_full_market,
  "Full (settlement gravity)"          = m_full_settgrav,
  "Full (settlement pop.)"             = m_full_settpop,
  # Combined gravity
  "Both gravity"                       = m_both_grav,
  "Habitat + both gravity"             = m_hab_both_grav,
  "Env + both gravity"                 = m_env_both_grav,
  "Full (both gravity)"                = m_full_both_grav
)

cat("\n--- AICc: Transect-level large excavator biomass ---\n")
print(make_aicc_df(model_list_transect))

# ── Diagnostics on best model ─────────────────────────────────
# Habitat only: top model (weight = 0.2476).
# Habitat + settlement pop.: ΔAICc = 0.36 (weight = 0.2071).
# Both models show excellent diagnostics — no dispersion issues,
# no zero inflation (both p > 0.63), no outliers.
# Settlement pop. coefficient not significant (beta = 0.422,
# p = 0.188), confirming habitat only as the most parsimonious
# description. Rugosity robust: beta = 0.883 (habitat only)
# vs 0.930 (habitat + settlement pop.).

cat("\n--- Diagnostics: Habitat only ---\n")
res_t_hab <- simulateResiduals(m_hab, n = 1000)

jpeg("dharma_excavator_transect_hab.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_t_hab, main = "DHARMa — transect excavator: Habitat only"); dev.off()

plot(res_t_hab)
testZeroInflation(res_t_hab)
testDispersion(res_t_hab)
testOutliers(res_t_hab)
plotResiduals(res_t_hab, transect_model_data$rugosity_sc, xlab = "Rugosity")

cat("\n--- Diagnostics: Habitat + settlement pop. ---\n")
res_t_hab_settpop <- simulateResiduals(m_hab_settpop, n = 1000)

jpeg("dharma_excavator_transect_hab_settpop.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_t_hab_settpop, main = "DHARMa — transect excavator: Habitat + settlement pop."); dev.off()

plot(res_t_hab_settpop)
testZeroInflation(res_t_hab_settpop)
testDispersion(res_t_hab_settpop)
testOutliers(res_t_hab_settpop)
plotResiduals(res_t_hab_settpop, transect_model_data$rugosity_sc,          xlab = "Rugosity")
plotResiduals(res_t_hab_settpop, transect_model_data$log_settlement_pop_sc, xlab = "Settlement pop.")

cat("\n--- Summary: Habitat only ---\n")
summary(m_hab)

cat("\n--- Summary: Habitat + settlement pop. ---\n")
summary(m_hab_settpop)

# ── Convergence with site-level result ────────────────────────
# Strong convergence: habitat only is the top model at both
# site (weight = 0.33) and transect (weight = 0.25) levels.
# Rugosity coefficient highly stable across both levels
# (site beta = 0.856; transect beta = 0.883).
# Confirms site-level result is not an artefact of averaging.

# ── Marginal effect plots ─────────────────────────────────────
( p_t_rugosity <- plot_effect(m_hab,
                              transect_model_data,
                              "rugosity_sc",
                              "Rugosity (scaled)",
                              y_label = "Large excavator biomass (g)") )

( p_t_settpop <- plot_effect(m_hab_settpop,
                             transect_model_data,
                             "log_settlement_pop_sc",
                             "Settlement pop. (scaled)",
                             y_label = "Large excavator biomass (g)") )

jpeg("transect_excavator_marginal_effects.jpg", width = 22, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(p_t_rugosity, p_t_settpop, ncol = 2)
dev.off()

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
#  Response: Total large excavator count per transect (integer >= 0)
#  Family:   Poisson → NB2 → NB1; selected via AICc
#  Random fx: (1 | site)
# ============================================================

# ── Explore count distribution ────────────────────────────────
cat("Transects (count data):", nrow(transect_model_data), "\n")
cat("Zeros in count data:",    sum(transect_model_data$transect_excavator_count == 0), "\n")
cat("Proportion zeros:",       round(mean(transect_model_data$transect_excavator_count == 0), 3), "\n")

summary(transect_model_data$transect_excavator_count)

ggplot(transect_model_data, aes(x = transect_excavator_count)) +
  geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
  labs(x = "Total large excavator count per transect", y = "Frequency",
       title = "Raw large excavator count distribution") +
  theme_bw()

# ── Mean-variance relationship ────────────────────────────────
# Points above the Poisson line → overdispersion → Negative Binomial.
transect_model_data %>%
  group_by(site) %>%
  summarise(mean_count = mean(transect_excavator_count),
            var_count  = var(transect_excavator_count),
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
# C1: Poisson
mC1_poisson <- glmmTMB(
  transect_excavator_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = poisson(link = "log"),
  data   = transect_model_data
)

resC1 <- simulateResiduals(mC1_poisson, n = 1000)

jpeg("dharma_excavator_C1_poisson.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC1, main = "DHARMa — Poisson"); dev.off()

plot(resC1)
testDispersion(resC1)
testZeroInflation(resC1)
testOutliers(resC1)

# C2: NB2 — quadratic variance (classic NB)
mC2_nb2 <- glmmTMB(
  transect_excavator_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = nbinom2(link = "log"),
  data   = transect_model_data
)

resC2 <- simulateResiduals(mC2_nb2, n = 1000)

jpeg("dharma_excavator_C2_nb2.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC2, main = "DHARMa — NB2"); dev.off()

plot(resC2)
testDispersion(resC2)
testZeroInflation(resC2)
testOutliers(resC2)

# C3: NB1 — linear variance
mC3_nb1 <- glmmTMB(
  transect_excavator_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = nbinom1(link = "log"),
  data   = transect_model_data
)

resC3 <- simulateResiduals(mC3_nb1, n = 1000)

jpeg("dharma_excavator_C3_nb1.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC3, main = "DHARMa — NB1"); dev.off()

plot(resC3)
testDispersion(resC3)
testZeroInflation(resC3)
testOutliers(resC3)

cat("\n--- Family selection: large excavator count models ---\n")
print(make_aicc_df(list(
  "Poisson" = mC1_poisson,
  "NB2"     = mC2_nb2,
  "NB1"     = mC3_nb1
)))

# ── Family selection decision ─────────────────────────────────
# NB1 clearly preferred (AICc = 783.56; ΔAICc = 2.61 over NB2,
# weight = 0.79 vs 0.21). Poisson strongly rejected (ΔAICc = 69.96).
#
# DHARMa diagnostics:
#   NB1: no dispersion issues (p = 0.646), no zero inflation
#        (p = 0.982), no outliers — clean fit.
#   NB2: no dispersion issues (p = 0.970), no zero inflation
#        (p = 0.752), no outliers — also adequate.
#   Poisson: significant zero inflation (p = 0.012) and step
#        failure warning during DHARMa smoothing — inappropriate.
#
# NB1 (linear variance) carried forward for all candidate count models.

# ── Random effect structure selection ────────────────────────
count_family <- nbinom1(link = "log")

re_c_null   <- glmmTMB(transect_excavator_count ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc,
                       family = count_family, data = transect_model_data)

re_c_site   <- glmmTMB(transect_excavator_count ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | site),
                       family = count_family, data = transect_model_data)

re_c_nested <- glmmTMB(transect_excavator_count ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | country/site),
                       family = count_family, data = transect_model_data)

cat("\n--- RE structure comparison (large excavator counts) ---\n")
print(make_aicc_df(list(
  "No RE"              = re_c_null,
  "(1 | site)"         = re_c_site,
  "(1 | country/site)" = re_c_nested
)))

# (1 | site) clearly preferred (AICc = 783.56; delta = 2.12 over
# (1 | country/site), weight = 0.74 vs 0.26). No RE strongly
# rejected (delta = 12.60). Site-level random intercept carried
# forward — the added country-level nesting is not supported.

# ── Candidate count models ────────────────────────────────────
# Family: NB1 (log link). RE: (1 | site).

# --- Null ---
cm0 <- glmmTMB(
  transect_excavator_count ~ 1 + (1 | site),
  family = count_family,
  data   = transect_model_data
)

# --- Single predictor ---
cm_env      <- glmmTMB(transect_excavator_count ~ sst_sc + log_chla_sc + (1 | site), family = count_family, data = transect_model_data)
cm_market   <- glmmTMB(transect_excavator_count ~ log_market_gravity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_settgrav <- glmmTMB(transect_excavator_count ~ log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_settpop  <- glmmTMB(transect_excavator_count ~ log_settlement_pop_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab      <- glmmTMB(transect_excavator_count ~ rugosity_sc + (1 | site), family = count_family, data = transect_model_data)

# --- Environment + human pressure ---
cm_env_mkt      <- glmmTMB(transect_excavator_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_env_settgrav <- glmmTMB(transect_excavator_count ~ sst_sc + log_chla_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_env_settpop  <- glmmTMB(transect_excavator_count ~ sst_sc + log_chla_sc + log_settlement_pop_sc + (1 | site), family = count_family, data = transect_model_data)

# --- Habitat + human pressure ---
cm_hab_market   <- glmmTMB(transect_excavator_count ~ rugosity_sc + log_market_gravity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab_settgrav <- glmmTMB(transect_excavator_count ~ rugosity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab_settpop  <- glmmTMB(transect_excavator_count ~ rugosity_sc + log_settlement_pop_sc + (1 | site), family = count_family, data = transect_model_data)

# --- Full (single human pressure metric) ---
cm_full_mkt      <- glmmTMB(transect_excavator_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_full_settgrav <- glmmTMB(transect_excavator_count ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_full_settpop  <- glmmTMB(transect_excavator_count ~ sst_sc + log_chla_sc + log_settlement_pop_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)

# --- Combined gravity metrics ---
cm_both_grav      <- glmmTMB(transect_excavator_count ~ log_market_gravity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab_both_grav  <- glmmTMB(transect_excavator_count ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_env_both_grav  <- glmmTMB(transect_excavator_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_full_both_grav <- glmmTMB(transect_excavator_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)

model_list_counts <- list(
  # Null
  "Null"                               = cm0,
  # Single predictor
  "Environment"                        = cm_env,
  "Market gravity"                     = cm_market,
  "Settlement gravity"                 = cm_settgrav,
  "Settlement pop."                    = cm_settpop,
  "Habitat"                            = cm_hab,
  # Environment + human pressure
  "Env + market gravity"               = cm_env_mkt,
  "Env + settlement gravity"           = cm_env_settgrav,
  "Env + settlement pop."              = cm_env_settpop,
  # Habitat + human pressure
  "Habitat + market gravity"           = cm_hab_market,
  "Habitat + settlement gravity"       = cm_hab_settgrav,
  "Habitat + settlement pop."          = cm_hab_settpop,
  # Full - single human pressure
  "Full (market gravity)"              = cm_full_mkt,
  "Full (settlement gravity)"          = cm_full_settgrav,
  "Full (settlement pop.)"             = cm_full_settpop,
  # Combined gravity
  "Both gravity"                       = cm_both_grav,
  "Habitat + both gravity"             = cm_hab_both_grav,
  "Env + both gravity"                 = cm_env_both_grav,
  "Full (both gravity)"                = cm_full_both_grav
)

cat("\n--- AICc: large excavator count models ---\n")
print(make_aicc_df(model_list_counts))

# ── Diagnostics on best count model ──────────────────────────
# Habitat only: top model (weight = 0.2230).
# Null ranks third (ΔAICc = 1.38, weight = 0.1121) — very close,
# indicating the habitat signal is weak in abundance.
# DHARMa diagnostics on habitat only: no dispersion issues
# (p = 0.648), no zero inflation (p = 1.000), no outlier concerns
# (p = 0.56). Clean fit.
# Rugosity coefficient marginally non-significant (beta = 0.233,
# p = 0.060), consistent with the near-equivalence of habitat and
# null models.

cat("\n--- Diagnostics: Habitat only (count) ---\n")
res_cm_hab <- simulateResiduals(cm_hab, n = 1000)

jpeg("dharma_excavator_count_hab.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_cm_hab, main = "DHARMa — large excavator counts: Habitat only"); dev.off()

plot(res_cm_hab)
testDispersion(res_cm_hab)
testZeroInflation(res_cm_hab)
testOutliers(res_cm_hab)
plotResiduals(res_cm_hab, transect_model_data$rugosity_sc, xlab = "Rugosity")

cat("\n--- Summary: Habitat only (count) ---\n")
summary(cm_hab)

# ── SYNTHESIS ─────────────────────────────────────────────────
# Part 1 (site biomass):     Habitat only (weight = 0.33)
# Part 2 (transect biomass): Habitat only (weight = 0.25)
# Part 3 (transect counts): Habitat only (weight = 0.22) ≈ Null
#   (ΔAICc = 1.38). Rugosity marginally non-significant (p = 0.060).
#   The near-equivalence of habitat and null confirms that rugosity
#   influences excavator body size more than numerical abundance.
#
# Rugosity is the only consistent predictor across all three parts.
# Its effect weakens from biomass to counts, suggesting that habitat
# complexity influences excavator body size more than numerical
# abundance — complex reefs support larger individuals rather than
# more individuals.
#
# Human pressure metrics (market gravity, settlement gravity,
# settlement pop.) show no reliable signal at any level. Where
# they appear in competitive models (delta < 2 in Parts 1 and 2),
# the improvement is marginal and inconsistent across metrics,
# providing no basis for inference about fishing pressure effects.
#
# Environmental predictors (SST, Chl-a) are absent from all
# competitive models at every level of analysis (all delta > 6). 
# Large excavator communities appear insensitive to the productivity and 
# temperature gradients captured here.

## ── IRR: habitat only count model ────────────────────────────
# Only rugosity retained — IRR summarises its effect on abundance.
# Note: rugosity is marginally non-significant (p = 0.060) and the
# confidence interval crosses 1 (0.99–1.61). Treat cautiously —
# the near-equivalence of habitat and null models (ΔAICc = 1.38)
# means this effect is not robustly supported.
cat("\n--- IRR: Habitat only (count) ---\n")
irr_hab <- exp(fixef(cm_hab)$cond)
se_hab  <- summary(cm_hab)$coefficients$cond[, "Std. Error"]

cat("Rugosity: IRR =", round(irr_hab["rugosity_sc"], 2),
    " (95% CI:", round(exp(log(irr_hab["rugosity_sc"]) - 1.96 * se_hab["rugosity_sc"]), 2),
    "-",         round(exp(log(irr_hab["rugosity_sc"]) + 1.96 * se_hab["rugosity_sc"]), 2),
    ", p = 0.060)\n")
