# ============================================================
#  SCRAPER / SMALL EXCAVATOR BIOMASS & ABUNDANCE — MIXED EFFECTS MODELS
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
#              Sections: no RE → country RE
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

scraper_transects <- fish_2009 %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_scraper_biomass = sum(
      ifelse(trophic_group %in% c("scrapers", "small_excavators"), tot_wt_g, 0),
      na.rm = TRUE
    ),
    transect_scraper_count = sum(
      ifelse(trophic_group %in% c("scrapers", "small_excavators"), number, 0),
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

cat("Number of transects:", nrow(scraper_transects), "\n")
cat("Number of sites:",     n_distinct(scraper_transects$site), "\n")
cat("Number of countries:", n_distinct(scraper_transects$country), "\n")

# ==============================================================================
#  BIOMASS DATA EXPLORATION
# ==============================================================================

# ── Aggregate data at site level ─────────────────────────────────────────────
site_data <- scraper_transects %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass = mean(transect_scraper_biomass, na.rm = TRUE),
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

# ── Raw distribution ──────────────────────────────────────────────────────────
( site_raw <- ggplot(site_data, aes(x = mean_biomass)) +
    geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
    labs(x = "Mean scraper / small excavator biomass per site (g)", y = "Frequency",
         title = "Raw Site-Level Scraper / Small Excavator Biomass") +
    theme_bw() )

# ── Box-Cox ───────────────────────────────────────────────────────────────────
site_nonzero <- site_data %>% filter(mean_biomass > 0)

MASS::boxcox(
  lm(mean_biomass ~ 1, data = site_nonzero),
  lambda = seq(-2, 2, 0.1)
)
# lambda ~ 0 → log transformation appropriate

# ── Apply transformations ─────────────────────────────────────────────────────
site_data <- site_data %>%
  mutate(
    log_mean_biomass  = log(mean_biomass + 0.01),
    sqrt_mean_biomass = sqrt(mean_biomass)
  )

( site_log <- ggplot(site_data, aes(x = log_mean_biomass)) +
    geom_histogram(bins = 25, fill = "#1a9641", colour = "white") +
    labs(x = "log(mean biomass + 0.01)", y = "Frequency",
         title = "Log-transformed Site-Level Scraper / Small Excavator Biomass") +
    theme_bw() )

( site_sqrt <- ggplot(site_data, aes(x = sqrt_mean_biomass)) +
    geom_histogram(bins = 25, fill = "#d7191c", colour = "white") +
    labs(x = "sqrt(mean biomass)", y = "Frequency",
         title = "Sqrt-transformed Site-Level Scraper / Small Excavator Biomass") +
    theme_bw() )

jpeg("site_scraper_biomass_distributions.jpg", width = 33, height = 11,
     units = "cm", res = 300)
gridExtra::grid.arrange(site_raw, site_log, site_sqrt, ncol = 3)
dev.off()

# ── Normality checks ──────────────────────────────────────────────────────────
qqnorm(site_data$log_mean_biomass,
       main = "Q-Q plot: log(mean scraper / small excavator biomass per site)")
qqline(site_data$log_mean_biomass, col = "red")
shapiro.test(site_data$log_mean_biomass)
# Shapiro-Wilk: W = 0.537, p < 0.001 — extreme departure from
# normality. The one zero site (ch_rubu) collapses to log(0.01) = -4.6,
# creating a severe outlier on the log scale. Tweedie on raw scale
# is strongly preferred.

# ── Variation by site ─────────────────────────────────────────────────────────
ggplot(site_data, aes(x = reorder(site, mean_biomass, median),
                      y = mean_biomass)) +
  geom_col(fill = "#2c7bb6", alpha = 0.7) +
  coord_flip() +
  labs(x = NULL, y = "Mean scraper / small excavator biomass (g)",
       title = "Mean scraper / small excavator biomass by site (raw)") +
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
  labs(x = NULL, y = "log(mean scraper / small excavator biomass per site + 0.01)") +
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

jpeg("predictor_distributions_scraper.jpg", width = 33, height = 22, units = "cm", res = 300)
gridExtra::grid.arrange(p_pred_raw, p_pred_tran, nrow = 2)
dev.off()

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

settlement_data <- scraper_transects %>%
  left_join(scaled_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    mean_scraper_biomass   = mean(transect_scraper_biomass, na.rm = TRUE),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_settlement_pop_sc  = first(log_settlement_pop_sc),
    .groups = "drop"
  )

settgrav <- glmmTMB(mean_scraper_biomass ~ log_settlement_grav_sc,
                    family = tweedie(link = "log"), data = settlement_data)
settpop  <- glmmTMB(mean_scraper_biomass ~ log_settlement_pop_sc,
                    family = tweedie(link = "log"), data = settlement_data)

make_aicc_df(list(
  "Settlement gravity" = settgrav,
  "Settlement pop."    = settpop
))

# Settlement metric selection: settlement gravity clearly preferred
# (AICc = 1207.93 vs 1211.11; delta = 3.18, weight = 0.83 vs 0.17).
# Unlike piscivores and large excavators where metrics were
# near-equivalent, settlement gravity has meaningful support over
# settlement pop. for scrapers / small excavators. Both metrics
# carried forward in the full candidate set for completeness, but
# settlement gravity is the primary human pressure metric for
# this group.

rm(settlement_data)

# ============================================================
#  ANALYSIS DATASETS
# ============================================================

final_predictors <- scaled_predictors %>%
  dplyr::select(site, log_market_gravity_sc, log_settlement_grav_sc,
                log_settlement_pop_sc, log_chla_sc, sst_sc, rugosity_sc)

# ── Transect-level dataset ────────────────────────────────────
transect_model_data <- scraper_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect model data:", nrow(transect_model_data), "rows,",
    n_distinct(transect_model_data$site), "sites\n")
cat("Biomass zeros:", sum(transect_model_data$transect_scraper_biomass == 0), "\n")
cat("Count zeros:",  sum(transect_model_data$transect_scraper_count   == 0), "\n")

# ── Site-level dataset ────────────────────────────────────────
total_model_data <- transect_model_data %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_scraper_biomass, na.rm = TRUE),
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
# Gaussian on log(y + 0.01) vs Tweedie.
# ZI-Tweedie not attempted at site level — sample size is
# insufficient to support the additional ZI component and
# convergence failures are expected.

# ── F1: Gaussian on log-transformed mean biomass ──────────────
total_model_data <- total_model_data %>%
  mutate(log_mean_biomass = log(mean_biomass + 0.01))

mS_F1 <- glmmTMB(log_mean_biomass ~ sst_sc + log_chla_sc +
                   log_market_gravity_sc + rugosity_sc,
                 family = gaussian(), data = total_model_data)

resS_F1 <- simulateResiduals(mS_F1, n = 1000)

jpeg("diagnostics_site_scraper_F1_gaussian_log.jpg",
     width = 25, height = 15, units = "cm", res = 300)
plot(resS_F1, main = "DHARMa — Gaussian on log(y + 0.01)"); dev.off()

plot(resS_F1)
testZeroInflation(resS_F1)
testDispersion(resS_F1)

# ── F2: Tweedie (log link) on raw mean biomass ────────────────
mS_F2 <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc +
                   log_market_gravity_sc + rugosity_sc,
                 family = tweedie(link = "log"), data = total_model_data)

resS_F2 <- simulateResiduals(mS_F2, n = 1000)

jpeg("diagnostics_site_scraper_F2_tweedie.jpg",
     width = 25, height = 15, units = "cm", res = 300)
plot(resS_F2, main = "DHARMa — Tweedie"); dev.off()

plot(resS_F2)
testZeroInflation(resS_F2)
testDispersion(resS_F2)

# ── Family selection decision ─────────────────────────────────
# Both F1 and F2 show adequate DHARMa diagnostics (64 sites, 1 zero).
# F1 (Gaussian on log scale): no dispersion issues (p = 0.914),
# zero inflation NaN (expected — no zeros in residuals at this level).
# F2 (Tweedie): no dispersion issues (p = 0.168), zero inflation
# n.s. (p = 0.962).
#
# Tweedie (F2) preferred: handles the 1 site-level zero natively
# without a log + constant approximation and operates on the raw
# biomass scale. ZI Tweedie not attempted — site n too small to
# support the ZI component.

# ── RANDOM EFFECT STRUCTURE SELECTION ────────────────────────
re_null <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc +
                     log_market_gravity_sc + rugosity_sc,
                   family = tweedie(link = "log"), data = total_model_data)

re_country <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc +
                        log_market_gravity_sc + rugosity_sc +
                        (1 | country),
                      family = tweedie(link = "log"), data = total_model_data)

cat("\n--- RE structure comparison (site-level scraper) ---\n")
print(make_aicc_df(list(
  "No RE"         = re_null,
  "(1 | country)" = re_country
)))

# ── RE structure decision ─────────────────────────────────────
# No RE clearly preferred (ΔAICc = 2.62, weight = 0.79 vs 0.21).
# Country-level variation adequately absorbed by included predictors.
# No RE carried forward for all site-level candidate models.

# ── CANDIDATE MODELS — SITE LEVEL ────────────────────────────
# Family: Tweedie (log link). No RE — preferred over (1 | country)
# per RE structure comparison above (delta = 2.62).
# s2 (country RE) set not fitted — no RE decision is clear-cut
# (ΔAICc = 2.62, weight = 0.79 vs 0.21).

# ── No random effects ─────────────────────────────────────────
s1_m0              <- glmmTMB(mean_biomass ~ 1,                                                                                          family = tweedie(link = "log"), data = total_model_data)
s1_m_env           <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc,                                                                      family = tweedie(link = "log"), data = total_model_data)
s1_m_market        <- glmmTMB(mean_biomass ~ log_market_gravity_sc,                                                                     family = tweedie(link = "log"), data = total_model_data)
s1_m_settgrav      <- glmmTMB(mean_biomass ~ log_settlement_grav_sc,                                                                    family = tweedie(link = "log"), data = total_model_data)
s1_m_settpop       <- glmmTMB(mean_biomass ~ log_settlement_pop_sc,                                                                     family = tweedie(link = "log"), data = total_model_data)
s1_m_hab           <- glmmTMB(mean_biomass ~ rugosity_sc,                                                                               family = tweedie(link = "log"), data = total_model_data)
s1_m_env_mkt       <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc,                                             family = tweedie(link = "log"), data = total_model_data)
s1_m_env_settgrav  <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc,                                            family = tweedie(link = "log"), data = total_model_data)
s1_m_env_settpop   <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc,                                             family = tweedie(link = "log"), data = total_model_data)
s1_m_hab_market    <- glmmTMB(mean_biomass ~ rugosity_sc + log_market_gravity_sc,                                                       family = tweedie(link = "log"), data = total_model_data)
s1_m_hab_settgrav  <- glmmTMB(mean_biomass ~ rugosity_sc + log_settlement_grav_sc,                                                      family = tweedie(link = "log"), data = total_model_data)
s1_m_hab_settpop   <- glmmTMB(mean_biomass ~ rugosity_sc + log_settlement_pop_sc,                                                       family = tweedie(link = "log"), data = total_model_data)
s1_m_full_mkt      <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc,                              family = tweedie(link = "log"), data = total_model_data)
s1_m_full_settgrav <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc,                              family = tweedie(link = "log"), data = total_model_data)
s1_m_full_settpop  <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc  + rugosity_sc,                              family = tweedie(link = "log"), data = total_model_data)
s1_m_both_grav     <- glmmTMB(mean_biomass ~ log_market_gravity_sc + log_settlement_grav_sc,                                            family = tweedie(link = "log"), data = total_model_data)
s1_m_hab_both_grav <- glmmTMB(mean_biomass ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc,                              family = tweedie(link = "log"), data = total_model_data)
s1_m_env_both_grav <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc,                    family = tweedie(link = "log"), data = total_model_data)
s1_m_full_both_grav<- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc,      family = tweedie(link = "log"), data = total_model_data)

model_list_s1 <- list(
  "Null"                               = s1_m0,
  "Environment"                        = s1_m_env,
  "Market gravity"                     = s1_m_market,
  "Settlement gravity"                 = s1_m_settgrav,
  "Settlement pop."                    = s1_m_settpop,
  "Habitat"                            = s1_m_hab,
  "Env + market gravity"               = s1_m_env_mkt,
  "Env + settlement gravity"           = s1_m_env_settgrav,
  "Env + settlement pop."              = s1_m_env_settpop,
  "Habitat + market gravity"           = s1_m_hab_market,
  "Habitat + settlement gravity"       = s1_m_hab_settgrav,
  "Habitat + settlement pop."          = s1_m_hab_settpop,
  "Full (market gravity)"              = s1_m_full_mkt,
  "Full (settlement gravity)"          = s1_m_full_settgrav,
  "Full (settlement pop.)"             = s1_m_full_settpop,
  "Both gravity"                       = s1_m_both_grav,
  "Habitat + both gravity"             = s1_m_hab_both_grav,
  "Env + both gravity"                 = s1_m_env_both_grav,
  "Full (both gravity)"                = s1_m_full_both_grav
)

cat("\n--- AICc: Site-level scraper candidate models (no RE) ---\n")
print(make_aicc_df(model_list_s1))

# ── Residual diagnostics — top models ────────────────────────
# Full (settlement pop.): top model (weight = 0.2749).
# Env + settlement pop.: ΔAICc = 0.98 (weight = 0.1683).
# Both models show no dispersion issues and no outlier concerns.
# Settlement pop. significant in both (p = 0.040, p = 0.048).
# Chl-a significant in both (p = 0.002). SST and rugosity marginal
# (p = 0.059–0.089) — more consistent signal in Full model.

cat("\n--- Diagnostics: Full (settlement pop.) ---\n")
res_s1_full_settpop <- simulateResiduals(s1_m_full_settpop, n = 1000)

jpeg("diagnostics_site_scraper_full_settpop.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_s1_full_settpop, main = "DHARMa — Full (settlement pop.)"); dev.off()

plot(res_s1_full_settpop)
testZeroInflation(res_s1_full_settpop)
testDispersion(res_s1_full_settpop)
testOutliers(res_s1_full_settpop)

cat("\n--- Diagnostics: Env + settlement pop. ---\n")
res_s1_env_settpop <- simulateResiduals(s1_m_env_settpop, n = 1000)

jpeg("diagnostics_site_scraper_env_settpop.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_s1_env_settpop, main = "DHARMa — Env + settlement pop."); dev.off()

plot(res_s1_env_settpop)
testZeroInflation(res_s1_env_settpop)
testDispersion(res_s1_env_settpop)
testOutliers(res_s1_env_settpop)

# ── Summaries ─────────────────────────────────────────────────
cat("\n--- Summary: Full (settlement pop.) ---\n")
summary(s1_m_full_settpop)

cat("\n--- Summary: Env + settlement pop. ---\n")
summary(s1_m_env_settpop)

# ── Coefficient stability across human pressure metrics ───────
cat("\n--- SST coefficient stability ---\n")
cat("Full (settlement pop.):  beta =", round(fixef(s1_m_full_settpop)$cond["sst_sc"],      3), "\n")
cat("Env + settlement pop.:   beta =", round(fixef(s1_m_env_settpop)$cond["sst_sc"],       3), "\n")
cat("Full (settlement grav):  beta =", round(fixef(s1_m_full_settgrav)$cond["sst_sc"],     3), "\n")

cat("\n--- Chl-a coefficient stability ---\n")
cat("Full (settlement pop.):  beta =", round(fixef(s1_m_full_settpop)$cond["log_chla_sc"], 3), "\n")
cat("Env + settlement pop.:   beta =", round(fixef(s1_m_env_settpop)$cond["log_chla_sc"],  3), "\n")
cat("Full (settlement grav):  beta =", round(fixef(s1_m_full_settgrav)$cond["log_chla_sc"],3), "\n")

# ── Marginal effect plots ─────────────────────────────────────
# From Full (settlement pop.) — top-ranked model
( p_site_sst <- plot_effect(s1_m_full_settpop,
                            total_model_data,
                            "sst_sc",
                            "SST (scaled)",
                            y_label = "Scraper / small excavator biomass (g)") )

( p_site_chla <- plot_effect(s1_m_full_settpop,
                             total_model_data,
                             "log_chla_sc",
                             "Chlorophyll-a (scaled)",
                             y_label = "Scraper / small excavator biomass (g)") )

( p_site_settpop <- plot_effect(s1_m_full_settpop,
                                total_model_data,
                                "log_settlement_pop_sc",
                                "Settlement pop. (scaled)",
                                y_label = "Scraper / small excavator biomass (g)") )

jpeg("site_scraper_marginal_effects.jpg", width = 33, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(p_site_sst, p_site_chla, p_site_settpop, ncol = 3)
dev.off()

# ============================================================
#  PART 2 — TRANSECT-LEVEL BIOMASS (SENSITIVITY CHECK)
#
#  Rationale: Retains within-site variation. (1 | site) accounts
#  for non-independence of transects.
#  Confirms site-level findings are not an artefact of collapsing
#  to site means.
#
#  Response:   transect_scraper_biomass — continuous, ~8.5% zeros
#  Family:     Tweedie (log link) — confirm via family selection below
#  Random fx:  (1 | site), then (1 | country/site)
# ============================================================

# ── Explore transect-level response ──────────────────────────
summary(scraper_transects$transect_scraper_biomass)

zeros <- mean(scraper_transects$transect_scraper_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")

( scraper_raw <- ggplot(scraper_transects, aes(x = transect_scraper_biomass)) +
    geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
    labs(x = "Scraper / small excavator biomass per transect (g)", y = "Frequency",
         title = "Raw Scraper / Small Excavator Biomass") +
    theme_bw() )

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

jpeg("scraper_biomass_distributions.jpg", width = 33, height = 11,
     units = "cm", res = 300)
gridExtra::grid.arrange(scraper_raw, scraper_log, scraper_sqrt, ncol = 3)
dev.off()

# Add log scraper biomass to transect model data
transect_model_data <- transect_model_data %>%
  mutate(log_scraper_biomass = log(transect_scraper_biomass + 0.01))

# ── Box-Cox on non-zero values ────────────────────────────────
scraper_nonzero <- scraper_transects %>% filter(transect_scraper_biomass > 0)

MASS::boxcox(
  lm(transect_scraper_biomass ~ 1, data = scraper_nonzero),
  lambda = seq(-2, 2, 0.1)
)

# ── Variation by site ─────────────────────────────────────────
ggplot(transect_model_data,
       aes(x = reorder(site, transect_scraper_biomass, median),
           y = transect_scraper_biomass)) +
  geom_boxplot(fill = "#2c7bb6", alpha = 0.6,
               outlier.colour = "black", outlier.size = 1) +
  coord_flip() +
  labs(x = NULL, y = "Scraper / small excavator biomass (g)",
       title = "Scraper / small excavator biomass distribution by site") +
  theme_bw(base_size = 9)

# ── Zeros by site ─────────────────────────────────────────────
transect_model_data %>%
  group_by(site) %>%
  summarise(
    n_transects  = n(),
    prop_zeros   = mean(transect_scraper_biomass == 0),
    mean_biomass = mean(transect_scraper_biomass),
    .groups = "drop"
  ) %>%
  arrange(desc(prop_zeros)) %>%
  print(n = Inf)

# ── Family selection ──────────────────────────────────────────
# AICc not comparable between F1 and F2/F3 (different response).
# Select on DHARMa diagnostics; use AICc only to compare F2 vs F3.

# F1: Gaussian on log(y + 0.01)
mF1_gaussian <- glmmTMB(
  log_scraper_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = gaussian(),
  data   = transect_model_data
)

resF1 <- simulateResiduals(mF1_gaussian, n = 1000)

jpeg("dharma_scraper_F1_gaussian.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF1, main = "DHARMa — Gaussian on log(y + 0.01)"); dev.off()

plot(resF1)
testZeroInflation(resF1)
testDispersion(resF1)

# F2: Plain Tweedie
mF2_tweedie <- glmmTMB(
  transect_scraper_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = transect_model_data
)

resF2 <- simulateResiduals(mF2_tweedie, n = 1000)

jpeg("dharma_scraper_F2_tweedie.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF2, main = "DHARMa — Tweedie"); dev.off()

plot(resF2)
testZeroInflation(resF2)
testDispersion(resF2)

# F3: Zero-inflated Tweedie
mF3_tweedie_zi <- glmmTMB(
  transect_scraper_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = transect_model_data
)

resF3 <- simulateResiduals(mF3_tweedie_zi, n = 1000)

jpeg("dharma_scraper_F3_tweedie_zi.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF3, main = "DHARMa — Zero-inflated Tweedie"); dev.off()

plot(resF3)
testZeroInflation(resF3)
testDispersion(resF3)

cat("\n--- Family selection: transect-level scraper biomass ---\n")
print(make_aicc_df(list(
  "Tweedie"    = mF2_tweedie,
  "ZI Tweedie" = mF3_tweedie_zi
)))

# ── Family selection decision ─────────────────────────────────
# Plain Tweedie (F2) strongly preferred. ZI Tweedie is clearly
# not supported (delta = 6.95, weight = 0.03 vs 0.97). Plain
# Tweedie carried forward for all transect-level candidate models.

# ── Random effect structure selection ────────────────────────
re_t_null   <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc,
                       family = tweedie(link = "log"), data = transect_model_data)

re_t_site   <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | site),
                       family = tweedie(link = "log"), data = transect_model_data)

re_t_nested <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | country/site),
                       family = tweedie(link = "log"), data = transect_model_data)

cat("\n--- RE structure comparison (transect-level scraper) ---\n")
print(make_aicc_df(list(
  "No RE"              = re_t_null,
  "(1 | site)"         = re_t_site,
  "(1 | country/site)" = re_t_nested
)))

# ── Random effect structure decision ─────────────────────────
# (1 | site) clearly preferred (ΔAICc = 2.13, weight = 0.74 vs 0.26).
# No RE strongly rejected (ΔAICc = 37.82).
# Site-level random intercept carried forward.

# ── Candidate models ──────────────────────────────────────────
# Family: Tweedie (log link). RE: (1 | site).
scraper_family <- tweedie(link = "log")

# --- Null ---
m0                 <- glmmTMB(transect_scraper_biomass ~ 1                                                                                        + (1 | site), family = scraper_family, data = transect_model_data)

# --- Single predictor ---
m_env              <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc                                                                    + (1 | site), family = scraper_family, data = transect_model_data)
m_market           <- glmmTMB(transect_scraper_biomass ~ log_market_gravity_sc                                                                   + (1 | site), family = scraper_family, data = transect_model_data)
m_settgrav         <- glmmTMB(transect_scraper_biomass ~ log_settlement_grav_sc                                                                  + (1 | site), family = scraper_family, data = transect_model_data)
m_settpop          <- glmmTMB(transect_scraper_biomass ~ log_settlement_pop_sc                                                                   + (1 | site), family = scraper_family, data = transect_model_data)
m_hab              <- glmmTMB(transect_scraper_biomass ~ rugosity_sc                                                                             + (1 | site), family = scraper_family, data = transect_model_data)

# --- Environment + human pressure ---
m_env_market       <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc                                           + (1 | site), family = scraper_family, data = transect_model_data)
m_env_settgrav     <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc                                          + (1 | site), family = scraper_family, data = transect_model_data)
m_env_settpop      <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc                                           + (1 | site), family = scraper_family, data = transect_model_data)

# --- Habitat + human pressure ---
m_hab_market       <- glmmTMB(transect_scraper_biomass ~ rugosity_sc + log_market_gravity_sc                                                     + (1 | site), family = scraper_family, data = transect_model_data)
m_hab_settgrav     <- glmmTMB(transect_scraper_biomass ~ rugosity_sc + log_settlement_grav_sc                                                    + (1 | site), family = scraper_family, data = transect_model_data)
m_hab_settpop      <- glmmTMB(transect_scraper_biomass ~ rugosity_sc + log_settlement_pop_sc                                                     + (1 | site), family = scraper_family, data = transect_model_data)

# --- Full (single human pressure metric) ---
m_full_market      <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc                            + (1 | site), family = scraper_family, data = transect_model_data)
m_full_settgrav    <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc                            + (1 | site), family = scraper_family, data = transect_model_data)
m_full_settpop     <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc  + rugosity_sc                            + (1 | site), family = scraper_family, data = transect_model_data)

# --- Combined gravity metrics ---
m_both_grav        <- glmmTMB(transect_scraper_biomass ~ log_market_gravity_sc + log_settlement_grav_sc                                         + (1 | site), family = scraper_family, data = transect_model_data)
m_hab_both_grav    <- glmmTMB(transect_scraper_biomass ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc                           + (1 | site), family = scraper_family, data = transect_model_data)
m_env_both_grav    <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc                  + (1 | site), family = scraper_family, data = transect_model_data)
m_full_both_grav   <- glmmTMB(transect_scraper_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc    + (1 | site), family = scraper_family, data = transect_model_data)

model_list_transect <- list(
  "Null"                               = m0,
  "Environment"                        = m_env,
  "Market gravity"                     = m_market,
  "Settlement gravity"                 = m_settgrav,
  "Settlement pop."                    = m_settpop,
  "Habitat"                            = m_hab,
  "Env + market gravity"               = m_env_market,
  "Env + settlement gravity"           = m_env_settgrav,
  "Env + settlement pop."              = m_env_settpop,
  "Habitat + market gravity"           = m_hab_market,
  "Habitat + settlement gravity"       = m_hab_settgrav,
  "Habitat + settlement pop."          = m_hab_settpop,
  "Full (market gravity)"              = m_full_market,
  "Full (settlement gravity)"          = m_full_settgrav,
  "Full (settlement pop.)"             = m_full_settpop,
  "Both gravity"                       = m_both_grav,
  "Habitat + both gravity"             = m_hab_both_grav,
  "Env + both gravity"                 = m_env_both_grav,
  "Full (both gravity)"                = m_full_both_grav
)

cat("\n--- AICc: Transect-level scraper biomass ---\n")
print(make_aicc_df(model_list_transect))

# ── Transect-level AICc results ───────────────────────────────
# Broadly consistent with site-level: environment and human
# pressure dominate. No single model clearly dominant —
# Full (settlement pop.) top (weight = 0.2133).
#
# Key differences from site level: rugosity appears in competitive
# transect models (Habitat + settlement gravity ΔAICc = 1.31;
# Habitat + both gravity ΔAICc = 1.83) but was absent at site
# level (ΔAICc > 5). Likely reflects within-site habitat variation
# captured at transect level but averaged away at site level.
# Settlement gravity competitive here (ΔAICc = 0.49 for Full
# settlement gravity) consistent with settlement metric selection.

# ── Diagnostics on best model ─────────────────────────────────
# Full (settlement pop.): top model (weight = 0.2133).
# Full (settlement gravity): ΔAICc = 0.49 (weight = 0.1673).
# Both models show clean diagnostics — no dispersion issues
# (p = 0.238 and p = 0.214), no zero inflation, no outliers.
# Rugosity now significant at transect level (beta = 0.220,
# p = 0.016) — present in site model but marginal (p = 0.059).

cat("\n--- Diagnostics: Full (settlement pop.) ---\n")
res_t_full_settpop <- simulateResiduals(m_full_settpop, n = 1000)

jpeg("dharma_scraper_transect_full_settpop.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_t_full_settpop, main = "DHARMa — transect scraper: Full (settlement pop.)"); dev.off()

plot(res_t_full_settpop)
testZeroInflation(res_t_full_settpop)
testDispersion(res_t_full_settpop)
testOutliers(res_t_full_settpop)
plotResiduals(res_t_full_settpop, transect_model_data$sst_sc,                  xlab = "SST")
plotResiduals(res_t_full_settpop, transect_model_data$log_chla_sc,             xlab = "Chl-a")
plotResiduals(res_t_full_settpop, transect_model_data$log_settlement_pop_sc,   xlab = "Settlement pop.")
plotResiduals(res_t_full_settpop, transect_model_data$rugosity_sc,             xlab = "Rugosity")

cat("\n--- Diagnostics: Full (settlement gravity) ---\n")
res_t_full_settgrav <- simulateResiduals(m_full_settgrav, n = 1000)

jpeg("dharma_scraper_transect_full_settgrav.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_t_full_settgrav, main = "DHARMa — transect scraper: Full (settlement gravity)"); dev.off()

plot(res_t_full_settgrav)
testZeroInflation(res_t_full_settgrav)
testDispersion(res_t_full_settgrav)
testOutliers(res_t_full_settgrav)
plotResiduals(res_t_full_settgrav, transect_model_data$sst_sc,                 xlab = "SST")
plotResiduals(res_t_full_settgrav, transect_model_data$log_chla_sc,            xlab = "Chl-a")
plotResiduals(res_t_full_settgrav, transect_model_data$log_settlement_grav_sc, xlab = "Settlement gravity")
plotResiduals(res_t_full_settgrav, transect_model_data$rugosity_sc,            xlab = "Rugosity")

cat("\n--- Summary: Full (settlement pop.) ---\n")
summary(m_full_settpop)

cat("\n--- Summary: Full (settlement gravity) ---\n")
summary(m_full_settgrav)

# ── Convergence with site-level result ────────────────────────
# Part 1 (site):     Full (settlement pop.)     (weight = 0.27)
# Part 2 (transect): Full (settlement pop.)     (weight = 0.21)
# Identical top model structure across both levels confirms that
# environment + settlement pop. signal is not an artefact of
# collapsing to site means. The modest appearance of rugosity in
# competitive transect models (but not site models) suggests
# habitat structure captures within-site variation in biomass
# that averages out at the site level.

# ── Marginal effect plots ─────────────────────────────────────
( p_t_sst <- plot_effect(m_full_settpop,
                         transect_model_data,
                         "sst_sc",
                         "SST (scaled)",
                         y_label = "Scraper / small excavator biomass (g)") )

( p_t_chla <- plot_effect(m_full_settpop,
                          transect_model_data,
                          "log_chla_sc",
                          "Chlorophyll-a (scaled)",
                          y_label = "Scraper / small excavator biomass (g)") )

( p_t_settpop <- plot_effect(m_full_settpop,
                             transect_model_data,
                             "log_settlement_pop_sc",
                             "Settlement pop. (scaled)",
                             y_label = "Scraper / small excavator biomass (g)") )

( p_t_rugosity <- plot_effect(m_full_settpop,
                              transect_model_data,
                              "rugosity_sc",
                              "Rugosity (scaled)",
                              y_label = "Scraper / small excavator biomass (g)") )

jpeg("transect_scraper_marginal_effects.jpg", width = 33, height = 22, units = "cm", res = 300)
gridExtra::grid.arrange(p_t_sst, p_t_chla, p_t_settpop, p_t_rugosity, ncol = 2)
dev.off()

# ============================================================
#  PART 3 — TRANSECT-LEVEL COUNTS (COMPLEMENTARY ANALYSIS)
#
#  Response: Total scraper / small excavator count per transect
#  Family:   Poisson → NB2 → NB1; selected via AICc
#  Random fx: (1 | site)
# ============================================================

# ── Explore count distribution ────────────────────────────────
cat("Transects (count data):", nrow(transect_model_data), "\n")
cat("Zeros in count data:",    sum(transect_model_data$transect_scraper_count == 0), "\n")
cat("Proportion zeros:",       round(mean(transect_model_data$transect_scraper_count == 0), 3), "\n")

summary(transect_model_data$transect_scraper_count)

ggplot(transect_model_data, aes(x = transect_scraper_count)) +
  geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
  labs(x = "Total scraper / small excavator count per transect", y = "Frequency",
       title = "Raw scraper / small excavator count distribution") +
  theme_bw()

# ── Mean-variance relationship ────────────────────────────────
transect_model_data %>%
  group_by(site) %>%
  summarise(mean_count = mean(transect_scraper_count),
            var_count  = var(transect_scraper_count),
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
mC1_poisson <- glmmTMB(
  transect_scraper_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = poisson(link = "log"),
  data   = transect_model_data
)

resC1 <- simulateResiduals(mC1_poisson, n = 1000)

jpeg("dharma_scraper_C1_poisson.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC1, main = "DHARMa — Poisson"); dev.off()

plot(resC1)
testDispersion(resC1)
testZeroInflation(resC1)
testOutliers(resC1)

mC2_nb2 <- glmmTMB(
  transect_scraper_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = nbinom2(link = "log"),
  data   = transect_model_data
)

resC2 <- simulateResiduals(mC2_nb2, n = 1000)

jpeg("dharma_scraper_C2_nb2.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC2, main = "DHARMa — NB2"); dev.off()

plot(resC2)
testDispersion(resC2)
testZeroInflation(resC2)
testOutliers(resC2)

mC3_nb1 <- glmmTMB(
  transect_scraper_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = nbinom1(link = "log"),
  data   = transect_model_data
)

resC3 <- simulateResiduals(mC3_nb1, n = 1000)

jpeg("dharma_scraper_C3_nb1.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC3, main = "DHARMa — NB1"); dev.off()

plot(resC3)
testDispersion(resC3)
testZeroInflation(resC3)
testOutliers(resC3)

cat("\n--- Family selection: scraper count models ---\n")
print(make_aicc_df(list(
  "Poisson" = mC1_poisson,
  "NB2"     = mC2_nb2,
  "NB1"     = mC3_nb1
)))

# ── Family selection decision ─────────────────────────────────
# NB2 strongly preferred (AICc = 2050.41; ΔAICc = 6.73 over NB1,
# weight = 0.97 vs 0.03). Poisson catastrophically rejected
# (ΔAICc = 1037.42) — severe overdispersion confirmed.
#
# DHARMa diagnostics:
#   NB2: no dispersion issues (p = 0.738); zero inflation
#        borderline (p = 0.048) — minor concern; no outliers.
#   NB1: no dispersion issues (p = 0.248); no zero inflation
#        (p = 0.904); no outliers — but AICc strongly disfavoured.
#   Poisson: extreme zero inflation (p < 0.001) — inappropriate.
#
# NB2 carried forward. The borderline zero inflation in NB2 is
# noted as a minor caveat but does not change the model selection.

# ── Random effect structure selection ────────────────────────
count_family <- nbinom2(link = "log")  # NB2 selected

re_c_null   <- glmmTMB(transect_scraper_count ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc,
                       family = count_family, data = transect_model_data)

re_c_site   <- glmmTMB(transect_scraper_count ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | site),
                       family = count_family, data = transect_model_data)

re_c_nested <- glmmTMB(transect_scraper_count ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | country/site),
                       family = count_family, data = transect_model_data)

cat("\n--- RE structure comparison (scraper counts) ---\n")
print(make_aicc_df(list(
  "No RE"              = re_c_null,
  "(1 | site)"         = re_c_site,
  "(1 | country/site)" = re_c_nested
)))

# ── Random effect structure decision ─────────────────────────
# (1 | site) marginally preferred (AICc = 2050.41; delta = 0.69
# over (1 | country/site), weight = 0.59 vs 0.41). Decision is
# less clear-cut than for biomass models or large excavator counts
# — both RE structures are plausible. Simpler (1 | site) carried
# forward per the prespecified rule of preferring parsimony when
# delta < 2. No RE strongly rejected (delta = 44.33).

# ── Candidate count models ────────────────────────────────────
# Family: NB2 (log link). RE: (1 | site).

cm0            <- glmmTMB(transect_scraper_count ~ 1 + (1 | site), family = count_family, data = transect_model_data)
cm_env         <- glmmTMB(transect_scraper_count ~ sst_sc + log_chla_sc + (1 | site), family = count_family, data = transect_model_data)
cm_market      <- glmmTMB(transect_scraper_count ~ log_market_gravity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_settgrav    <- glmmTMB(transect_scraper_count ~ log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_settpop     <- glmmTMB(transect_scraper_count ~ log_settlement_pop_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab         <- glmmTMB(transect_scraper_count ~ rugosity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_env_mkt     <- glmmTMB(transect_scraper_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_env_settgrav<- glmmTMB(transect_scraper_count ~ sst_sc + log_chla_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_env_settpop <- glmmTMB(transect_scraper_count ~ sst_sc + log_chla_sc + log_settlement_pop_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab_market  <- glmmTMB(transect_scraper_count ~ rugosity_sc + log_market_gravity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab_settgrav<- glmmTMB(transect_scraper_count ~ rugosity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab_settpop <- glmmTMB(transect_scraper_count ~ rugosity_sc + log_settlement_pop_sc + (1 | site), family = count_family, data = transect_model_data)
cm_full_mkt    <- glmmTMB(transect_scraper_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_full_settgrav<-glmmTMB(transect_scraper_count ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_full_settpop<- glmmTMB(transect_scraper_count ~ sst_sc + log_chla_sc + log_settlement_pop_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_both_grav   <- glmmTMB(transect_scraper_count ~ log_market_gravity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab_both    <- glmmTMB(transect_scraper_count ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_env_both    <- glmmTMB(transect_scraper_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_full_both   <- glmmTMB(transect_scraper_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)

model_list_counts <- list(
  "Null"                               = cm0,
  "Environment"                        = cm_env,
  "Market gravity"                     = cm_market,
  "Settlement gravity"                 = cm_settgrav,
  "Settlement pop."                    = cm_settpop,
  "Habitat"                            = cm_hab,
  "Env + market gravity"               = cm_env_mkt,
  "Env + settlement gravity"           = cm_env_settgrav,
  "Env + settlement pop."              = cm_env_settpop,
  "Habitat + market gravity"           = cm_hab_market,
  "Habitat + settlement gravity"       = cm_hab_settgrav,
  "Habitat + settlement pop."          = cm_hab_settpop,
  "Full (market gravity)"              = cm_full_mkt,
  "Full (settlement gravity)"          = cm_full_settgrav,
  "Full (settlement pop.)"             = cm_full_settpop,
  "Both gravity"                       = cm_both_grav,
  "Habitat + both gravity"             = cm_hab_both,
  "Env + both gravity"                 = cm_env_both,
  "Full (both gravity)"                = cm_full_both
)

cat("\n--- AICc: scraper count models ---\n")
print(make_aicc_df(model_list_counts))

# ── Count model AICc results ──────────────────────────────────
# Settlement gravity dominates — a striking contrast with both the
# biomass models (where settlement pop. was top-ranked) and with
# large excavator counts (where habitat was top-ranked).
# Evidence is concentrated in settlement gravity and both-gravity
# models: top 8 models all contain settlement gravity, with market
# gravity and settlement pop. effectively absent from competitive
# models (all delta > 15).
#
# Environmental predictors add modest signal when combined with
# settlement gravity (Env + settlement gravity delta = 0.88) but
# the single-predictor settlement gravity model is top-ranked,
# suggesting that human pressure is the primary driver of scraper /
# small excavator abundance.
#
# Habitat (rugosity) appears in competitive models (delta < 2.5)
# when combined with settlement gravity but provides no independent
# signal (Habitat alone delta = 28.19).
#
# Key contrast with biomass models:
# Settlement gravity drives abundance but settlement pop. drives
# biomass — the two metrics capture different aspects of human
# pressure and their relative importance differs between response
# variables for this group.

# ── Diagnostics on best count model ──────────────────────────
# Settlement gravity: top model (weight = 0.2543).
# Both gravity: ΔAICc = 0.58 (weight = 0.1902).
# DHARMa diagnostics on Settlement gravity: no dispersion issues
# (p = 0.808); zero inflation borderline (p = 0.052); no outliers.
# Settlement gravity coefficient strongly significant
# (beta = -0.596, p < 0.001; IRR = 0.55, 95% CI: 0.45–0.67).

cat("\n--- Diagnostics: Settlement gravity (count) ---\n")
res_cm_settgrav <- simulateResiduals(cm_settgrav, n = 1000)

jpeg("dharma_scraper_count_settgrav.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_cm_settgrav, main = "DHARMa — scraper counts: Settlement gravity"); dev.off()

plot(res_cm_settgrav)
testDispersion(res_cm_settgrav)
testZeroInflation(res_cm_settgrav)
testOutliers(res_cm_settgrav)
plotResiduals(res_cm_settgrav, transect_model_data$log_settlement_grav_sc, xlab = "Settlement gravity")

cat("\n--- Summary: Settlement gravity (count) ---\n")
summary(cm_settgrav)


# ── SYNTHESIS ─────────────────────────────────────────────────
# Part 1 (site biomass):     Full (settlement pop.)     (weight = 0.27)
# Part 2 (transect biomass): Full (settlement pop.)     (weight = 0.21)
# Part 3 (transect counts):  Settlement gravity (weight = 0.25)
#   IRR = 0.55 (95% CI: 0.45–0.67, p < 0.001): each unit increase
#   in settlement gravity associated with a 45% reduction in
#   scraper / small excavator abundance.
#
# Environment (SST, Chl-a) consistently drives biomass across both
# levels — the strongest and most consistent environmental signal
# across all trophic groups analysed. This contrasts sharply with
# large excavators where environment was irrelevant.
#
# Human pressure metrics show a notable split by response variable:
# settlement pop. is the top human pressure predictor for biomass
# (Parts 1 & 2) while settlement gravity dominates abundance (Part 3).
# Market gravity is uncompetitive across all three parts (all delta > 15
# in counts; marginal in biomass), unlike piscivores where market
# gravity appeared consistently.
#
# Rugosity is absent from top biomass models (site level delta > 5)
# but appears in competitive transect biomass and count models when
# combined with settlement gravity (delta < 2.5). This suggests
# habitat structure captures fine-scale within-site variation but
# is not a site-level driver of scraper / small excavator biomass.
#
# The biomass vs. count divergence in settlement metric (pop. vs.
# gravity) is noteworthy — settlement pop. reflects total local
# population size while settlement gravity incorporates distance
# weighting. That gravity predicts abundance but pop. predicts
# biomass may indicate that proximity-weighted pressure affects
# how many individuals are present, while total population size
# shapes size structure and therefore biomass.

# ── IRR: Settlement gravity (top count model) ─────────────────
cat("\n--- IRR: Settlement gravity (count) ---\n")
irr_settgrav <- exp(fixef(cm_settgrav)$cond)
se_settgrav  <- summary(cm_settgrav)$coefficients$cond[, "Std. Error"]

cat("Settlement gravity: IRR =", round(irr_settgrav["log_settlement_grav_sc"], 2),
    " (95% CI:", round(exp(log(irr_settgrav["log_settlement_grav_sc"]) - 1.96 * se_settgrav["log_settlement_grav_sc"]), 2),
    "-",         round(exp(log(irr_settgrav["log_settlement_grav_sc"]) + 1.96 * se_settgrav["log_settlement_grav_sc"]), 2), ")\n")