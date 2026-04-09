# ============================================================
#  BROWSER BIOMASS & ABUNDANCE — MIXED EFFECTS MODELS
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
#              NOTE: check zero proportion — may require Tweedie.
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

browser_transects <- fish_2009 %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_browser_biomass = sum(
      ifelse(trophic_group == "browsers", tot_wt_g, 0),
      na.rm = TRUE
    ),
    transect_browser_count = sum(
      ifelse(trophic_group == "browsers", number, 0),
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

cat("Number of transects:", nrow(browser_transects), "\n")
cat("Number of sites:",     n_distinct(browser_transects$site), "\n")
cat("Number of countries:", n_distinct(browser_transects$country), "\n")

# ==============================================================================
#  BIOMASS DATA EXPLORATION
# ==============================================================================

# ── Aggregate data at site level ─────────────────────────────────────────────
site_data <- browser_transects %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass = mean(transect_browser_biomass, na.rm = TRUE),
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
    labs(x = "Mean browser biomass per site (g)", y = "Frequency",
         title = "Raw Site-Level Browser Biomass") +
    theme_bw() )

# ── Box-Cox on non-zero values ────────────────────────────────────────────────
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
         title = "Log-transformed Site-Level Browser Biomass") +
    theme_bw() )

( site_sqrt <- ggplot(site_data, aes(x = sqrt_mean_biomass)) +
    geom_histogram(bins = 25, fill = "#d7191c", colour = "white") +
    labs(x = "sqrt(mean biomass)", y = "Frequency",
         title = "Sqrt-transformed Site-Level Browser Biomass") +
    theme_bw() )

jpeg("site_browser_biomass_distributions.jpg", width = 33, height = 11,
     units = "cm", res = 300)
gridExtra::grid.arrange(site_raw, site_log, site_sqrt, ncol = 3)
dev.off()

# ── Normality checks ──────────────────────────────────────────────────────────
qqnorm(site_data$log_mean_biomass,
       main = "Q-Q plot: log(mean browser biomass per site)")
qqline(site_data$log_mean_biomass, col = "red")
shapiro.test(site_data$log_mean_biomass)

# ── Variation by site ─────────────────────────────────────────────────────────
ggplot(site_data, aes(x = reorder(site, mean_biomass, median),
                      y = mean_biomass)) +
  geom_col(fill = "#2c7bb6", alpha = 0.7) +
  coord_flip() +
  labs(x = NULL, y = "Mean browser biomass (g)",
       title = "Mean browser biomass by site (raw)") +
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
  labs(x = NULL, y = "log(mean browser biomass per site + 0.01)") +
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
    market_gravity  = mean(market_grav,        na.rm = TRUE),
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

jpeg("predictor_distributions_browser.jpg", width = 33, height = 22, units = "cm", res = 300)
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

# ============================================================
#  CHOOSING SETTLEMENT METRIC
# ============================================================

settlement_data <- browser_transects %>%
  left_join(scaled_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    mean_browser_biomass   = mean(transect_browser_biomass, na.rm = TRUE),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_settlement_pop_sc  = first(log_settlement_pop_sc),
    .groups = "drop"
  )

settgrav <- glmmTMB(mean_browser_biomass ~ log_settlement_grav_sc,
                    family = tweedie(link = "log"), data = settlement_data)
settpop  <- glmmTMB(mean_browser_biomass ~ log_settlement_pop_sc,
                    family = tweedie(link = "log"), data = settlement_data)

make_aicc_df(list(
  "Settlement gravity" = settgrav,
  "Settlement pop."    = settpop
))

# Settlement metric selection: virtually identical performance
# (delta AICc = 0.19, weights 0.52 vs 0.48). Neither metric
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
transect_model_data <- browser_transects %>%
  left_join(final_predictors, by = "site") %>%
  mutate(log_browser_biomass = log(transect_browser_biomass + 0.01))

cat("\nTransect model data:", nrow(transect_model_data), "rows,",
    n_distinct(transect_model_data$site), "sites\n")
cat("Biomass zeros:", sum(transect_model_data$transect_browser_biomass == 0), "\n")
cat("Count zeros:",  sum(transect_model_data$transect_browser_count   == 0), "\n")

# ── Site-level dataset ────────────────────────────────────────
total_model_data <- transect_model_data %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_browser_biomass, na.rm = TRUE),
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
    mean_biomass     = mean_biomass,
    log_mean_biomass = log(mean_biomass + 0.01),
    site             = as.factor(site),
    country          = as.factor(country)
  )

cat("\nSite model data:", nrow(total_model_data), "sites,",
    n_distinct(total_model_data$country), "countries\n")
cat("Site-level zeros:", sum(total_model_data$mean_biomass == 0), "\n")

# ============================================================
#  PART 1 — SITE-LEVEL ANALYSIS (PRIMARY)
# ============================================================

# ── FAMILY SELECTION ──────────────────────────────────────────
# Anchor: full settlement gravity model (marginally preferred
# in settlement metric selection; robust to anchor choice).

# F1: Gaussian on log(y + 0.01)
mS_F1 <- glmmTMB(log_mean_biomass ~ sst_sc + log_chla_sc +
                   log_settlement_grav_sc + rugosity_sc,
                 family = gaussian(), data = total_model_data)

resS_F1 <- simulateResiduals(mS_F1, n = 1000)

jpeg("diagnostics_site_browser_F1_gaussian_log.jpg",
     width = 25, height = 15, units = "cm", res = 300)
plot(resS_F1, main = "DHARMa — Gaussian on log(y + 0.01)"); dev.off()

plot(resS_F1)
testZeroInflation(resS_F1)
testDispersion(resS_F1)

# F2: Tweedie (log link) on raw mean biomass
mS_F2 <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc +
                   log_settlement_grav_sc + rugosity_sc,
                 family = tweedie(link = "log"), data = total_model_data)

resS_F2 <- simulateResiduals(mS_F2, n = 1000)

jpeg("diagnostics_site_browser_F2_tweedie.jpg",
     width = 25, height = 15, units = "cm", res = 300)
plot(resS_F2, main = "DHARMa — Tweedie"); dev.off()

plot(resS_F2)
testZeroInflation(resS_F2)
testDispersion(resS_F2)

# ── F3: Zero-inflated Tweedie — SKIPPED at site level ────────
# With only 64 sites and 11 zero means (17%), the ZI Bernoulli
# component is poorly identified — glmmTMB returns false convergence (8).
# The plain Tweedie (F2) already handles the zero mass via its compound
# Poisson structure. ZI Tweedie is not warranted here.
# Decision: proceed with plain Tweedie (F2) for all site-level models.

cat("\n--- Family selection: site-level browser ---\n")
cat("F2 (Tweedie): AICc =", round(AICc(mS_F2), 2), "\n")
cat("F3 (ZI Tweedie): convergence failure — excluded\n")
cat("Decision: plain Tweedie selected\n")

# ── RANDOM EFFECT STRUCTURE SELECTION ────────────────────────
re_null <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc +
                     log_settlement_grav_sc + rugosity_sc,
                   family = tweedie(link = "log"), data = total_model_data)

re_country <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc +
                        log_settlement_grav_sc + rugosity_sc +
                        (1 | country),
                      family = tweedie(link = "log"), data = total_model_data)

cat("\n--- RE structure comparison (site-level browser) ---\n")
print(make_aicc_df(list(
  "No RE"         = re_null,
  "(1 | country)" = re_country
)))

# ── RE structure decision — browser site level ────────────────
# Country RE clearly preferred (delta AICc = 1.91, weight = 0.72
# vs 0.28 for no RE).
#
# Browser biomass likely reflects country-level differences in
# grazing pressure, herbivore management, or reef condition
# that are not fully captured by the included predictors.
#
# All candidate models fitted with (1 | country).
# RE models (s1) retained for reference only — conclusions
# based on country RE models (s2).

# ── CANDIDATE MODELS — SITE LEVEL ────────────────────────────

# --- No country ------
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

cat("\n--- AICc: Site-level browser candidate models (no RE) ---\n")
print(make_aicc_df(model_list_s1))

summary(s1_m_hab_settpop)
summary(s1_m_hab)

# ── CANDIDATE MODELS — COUNTRY RE (PRIMARY) ──────────────────

# --- Null ---
s2_m0                 <- glmmTMB(mean_biomass ~ 1                                                                                        + (1 | country), family = tweedie(link = "log"), data = total_model_data)

# --- Single predictor ---
s2_m_env              <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc                                                                    + (1 | country), family = tweedie(link = "log"), data = total_model_data)
s2_m_market           <- glmmTMB(mean_biomass ~ log_market_gravity_sc                                                                   + (1 | country), family = tweedie(link = "log"), data = total_model_data)
s2_m_settgrav         <- glmmTMB(mean_biomass ~ log_settlement_grav_sc                                                                  + (1 | country), family = tweedie(link = "log"), data = total_model_data)
s2_m_settpop          <- glmmTMB(mean_biomass ~ log_settlement_pop_sc                                                                   + (1 | country), family = tweedie(link = "log"), data = total_model_data)
s2_m_hab              <- glmmTMB(mean_biomass ~ rugosity_sc                                                                             + (1 | country), family = tweedie(link = "log"), data = total_model_data)

# --- Environment + human pressure ---
s2_m_env_mkt          <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc                                           + (1 | country), family = tweedie(link = "log"), data = total_model_data)
s2_m_env_settgrav     <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc                                          + (1 | country), family = tweedie(link = "log"), data = total_model_data)
s2_m_env_settpop      <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc                                           + (1 | country), family = tweedie(link = "log"), data = total_model_data)

# --- Habitat + human pressure ---
s2_m_hab_market       <- glmmTMB(mean_biomass ~ rugosity_sc + log_market_gravity_sc                                                     + (1 | country), family = tweedie(link = "log"), data = total_model_data)
s2_m_hab_settgrav     <- glmmTMB(mean_biomass ~ rugosity_sc + log_settlement_grav_sc                                                    + (1 | country), family = tweedie(link = "log"), data = total_model_data)
s2_m_hab_settpop      <- glmmTMB(mean_biomass ~ rugosity_sc + log_settlement_pop_sc                                                     + (1 | country), family = tweedie(link = "log"), data = total_model_data)

# --- Full (single human pressure metric) ---
s2_m_full_mkt         <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc                            + (1 | country), family = tweedie(link = "log"), data = total_model_data)
s2_m_full_settgrav    <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc                            + (1 | country), family = tweedie(link = "log"), data = total_model_data)
s2_m_full_settpop     <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc  + rugosity_sc                            + (1 | country), family = tweedie(link = "log"), data = total_model_data)

# --- Combined gravity metrics ---
s2_m_both_grav        <- glmmTMB(mean_biomass ~ log_market_gravity_sc + log_settlement_grav_sc                                         + (1 | country), family = tweedie(link = "log"), data = total_model_data)
s2_m_hab_both_grav    <- glmmTMB(mean_biomass ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc                           + (1 | country), family = tweedie(link = "log"), data = total_model_data)
s2_m_env_both_grav    <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc                  + (1 | country), family = tweedie(link = "log"), data = total_model_data)
s2_m_full_both_grav   <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc    + (1 | country), family = tweedie(link = "log"), data = total_model_data)

model_list_s2 <- list(
  # Null
  "Null"                               = s2_m0,
  # Single predictor
  "Environment"                        = s2_m_env,
  "Market gravity"                     = s2_m_market,
  "Settlement gravity"                 = s2_m_settgrav,
  "Settlement pop."                    = s2_m_settpop,
  "Habitat"                            = s2_m_hab,
  # Environment + human pressure
  "Env + market gravity"               = s2_m_env_mkt,
  "Env + settlement gravity"           = s2_m_env_settgrav,
  "Env + settlement pop."              = s2_m_env_settpop,
  # Habitat + human pressure
  "Habitat + market gravity"           = s2_m_hab_market,
  "Habitat + settlement gravity"       = s2_m_hab_settgrav,
  "Habitat + settlement pop."          = s2_m_hab_settpop,
  # Full - single human pressure
  "Full (market gravity)"              = s2_m_full_mkt,
  "Full (settlement gravity)"          = s2_m_full_settgrav,
  "Full (settlement pop.)"             = s2_m_full_settpop,
  # Combined gravity
  "Both gravity"                       = s2_m_both_grav,
  "Habitat + both gravity"             = s2_m_hab_both_grav,
  "Env + both gravity"                 = s2_m_env_both_grav,
  "Full (both gravity)"                = s2_m_full_both_grav
)

cat("\n--- AICc: Site-level browser candidate models (1 | country) ---\n")
aicc_site_s2 <- make_aicc_df(model_list_s2)
print(aicc_site_s2)

# ── Cross-RE comparison ───────────────────────────────────────
cat("\n--- Top 5 models: No RE (s1) ---\n")
print(head(aicc_site_s1, 5))

cat("\n--- Top 5 models: (1 | country) (s2) ---\n")
print(head(aicc_site_s2, 5))

# ── Diagnostics — s1 (no RE) ─────────────────────────────────
res_s1_hab_settpop <- simulateResiduals(s1_m_hab_settpop, n = 1000)

jpeg("diagnostics_site_browser_s1_hab_settpop.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_s1_hab_settpop, main = "DHARMa — s1: Habitat + settlement pop."); dev.off()

plot(res_s1_hab_settpop)
testZeroInflation(res_s1_hab_settpop)
testDispersion(res_s1_hab_settpop)
testOutliers(res_s1_hab_settpop)

res_s1_hab <- simulateResiduals(s1_m_hab, n = 1000)

jpeg("diagnostics_site_browser_s1_hab.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_s1_hab, main = "DHARMa — s1: Habitat only"); dev.off()

plot(res_s1_hab)
testZeroInflation(res_s1_hab)
testDispersion(res_s1_hab)
testOutliers(res_s1_hab)

summary(s1_m_hab_settpop)
summary(s1_m_hab)

# ── Diagnostics — s2 (1 | country) ───────────────────
res_s2_hab_settpop <- simulateResiduals(s2_m_hab_settpop, n = 1000)

jpeg("diagnostics_site_browser_s2_hab_settpop.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_s2_hab_settpop, main = "DHARMa — s2: Habitat + settlement pop. (1 | country)"); dev.off()

plot(res_s2_hab_settpop)
testZeroInflation(res_s2_hab_settpop)
testDispersion(res_s2_hab_settpop)
testOutliers(res_s2_hab_settpop)

res_s2_hab <- simulateResiduals(s2_m_hab, n = 1000)

jpeg("diagnostics_site_browser_s2_hab.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_s2_hab, main = "DHARMa — s2: Habitat only (1 | country)"); dev.off()

plot(res_s2_hab)
testZeroInflation(res_s2_hab)
testDispersion(res_s2_hab)
testOutliers(res_s2_hab)

summary(s2_m_hab_settpop)
summary(s2_m_hab)

# ── Coefficient comparison across RE structures ───────────────
cat("\n--- Rugosity coefficient stability ---\n")
cat("s1 (no RE):      beta =", round(fixef(s1_m_hab)$cond["rugosity_sc"],         3), "\n")
cat("s2 (country RE): beta =", round(fixef(s2_m_hab)$cond["rugosity_sc"],         3), "\n")

cat("\n--- Settlement pop. coefficient stability ---\n")
cat("s1 (no RE):      beta =", round(fixef(s1_m_hab_settpop)$cond["log_settlement_pop_sc"], 3), "\n")
cat("s2 (country RE): beta =", round(fixef(s2_m_hab_settpop)$cond["log_settlement_pop_sc"], 3), "\n")


# ── Marginal effect plots — site level ────────────────────────
# Rugosity — from habitat-only model (most parsimonious rugosity estimate)
( p_site_rugosity <- plot_effect(s1_m_hab,
                                 total_model_data,
                                 "rugosity_sc",
                                 "Rugosity (scaled)",
                                 y_label = "Browser biomass (g)") )
# Settlement pop. — top-ranked human pressure metric
( p_site_settpop <- plot_effect(s1_m_hab_settpop,
                                total_model_data,
                                "log_settlement_pop_sc",
                                "Settlement pop. (scaled)",
                                y_label = "Browser biomass (g)") )
# Settlement gravity — near-equivalent alternative
( p_site_settgrav <- plot_effect(s1_m_hab_settgrav,
                                 total_model_data,
                                 "log_settlement_grav_sc",
                                 "Settlement gravity (scaled)",
                                 y_label = "Browser biomass (g)") )
# Market gravity — for comparison with other functional groups
( p_site_market <- plot_effect(s1_m_hab_market,
                               total_model_data,
                               "log_market_gravity_sc",
                               "Market gravity (scaled)",
                               y_label = "Browser biomass (g)") )

# s1 plots
jpeg("site_browser_marginal_effects_s1_noRE.jpg", width = 33, height = 22, units = "cm", res = 300)
gridExtra::grid.arrange(p_site_rugosity, p_site_settpop,
                        p_site_settgrav, p_site_market, ncol = 2)
dev.off()

# ── Marginal effect plots - country ─────────────────────────────────────
# Rugosity from habitat-only model
( p_site_rugosity <- plot_effect(s2_m_hab,
                                 total_model_data,
                                 "rugosity_sc",
                                 "Rugosity (scaled)",
                                 y_label = "Browser biomass (g)") )

# Settlement pop. — top-ranked human pressure metric
( p_site_settpop <- plot_effect(s2_m_hab_settpop,
                                total_model_data,
                                "log_settlement_pop_sc",
                                "Settlement pop. (scaled)",
                                y_label = "Browser biomass (g)") )

# Settlement gravity — near-equivalent alternative
( p_site_settgrav <- plot_effect(s2_m_hab_settgrav,
                                 total_model_data,
                                 "log_settlement_grav_sc",
                                 "Settlement gravity (scaled)",
                                 y_label = "Browser biomass (g)") )

# Market gravity — for comparison with other functional groups
( p_site_market <- plot_effect(s2_m_hab_market,
                               total_model_data,
                               "log_market_gravity_sc",
                               "Market gravity (scaled)",
                               y_label = "Browser biomass (g)") )

# s2 plots — primary
jpeg("site_browser_marginal_effects_s2_country.jpg", width = 33, height = 22, units = "cm", res = 300)
gridExtra::grid.arrange(p_site_rugosity, p_site_settpop,
                        p_site_settgrav, p_site_market, ncol = 2)
dev.off()

# ============================================================
#  PART 2 — TRANSECT-LEVEL BIOMASS (SENSITIVITY CHECK)
#
#  Rationale: Retains within-site variation. (1 | site) accounts
#  for non-independence of transects.
#  Confirms site-level findings are not an artefact of collapsing
#  to site means.
#
#  Response:   transect_browser_biomass — continuous, possibly zero-inflated
#  Family:     Tweedie (log link) — confirmed via DHARMa ladder below
#  Random fx:  (1 | site), then (1 | country/site)
# ============================================================

# ── Explore transect-level response ──────────────────────────
summary(browser_transects$transect_browser_biomass)

zeros <- mean(browser_transects$transect_browser_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")

( browser_raw <- ggplot(browser_transects, aes(x = transect_browser_biomass)) +
    geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
    labs(x = "Browser biomass per transect (g)", y = "Frequency",
         title = "Raw Browser Biomass") +
    theme_bw() )

( browser_log <- ggplot(transect_model_data, aes(x = log_browser_biomass)) +
    geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
    labs(x = "log(biomass + 0.01)", y = "Frequency",
         title = "Log-transformed Browser Biomass") +
    theme_bw() )

( browser_sqrt <- ggplot(transect_model_data, aes(x = sqrt(transect_browser_biomass))) +
    geom_histogram(bins = 30, fill = "#d7191c", colour = "white") +
    labs(x = "sqrt(biomass)", y = "Frequency",
         title = "Sqrt-transformed Browser Biomass") +
    theme_bw() )

jpeg("browser_biomass_distributions.jpg", width = 33, height = 11,
     units = "cm", res = 300)
gridExtra::grid.arrange(browser_raw, browser_log, browser_sqrt, ncol = 3)
dev.off()

# ── Box-Cox on non-zero values ────────────────────────────────
browser_nonzero <- browser_transects %>% filter(transect_browser_biomass > 0)

MASS::boxcox(
  lm(transect_browser_biomass ~ 1, data = browser_nonzero),
  lambda = seq(-2, 2, 0.1)
)

# ── Variation by site ─────────────────────────────────────────
ggplot(transect_model_data,
       aes(x = reorder(site, transect_browser_biomass, median),
           y = transect_browser_biomass)) +
  geom_boxplot(fill = "#2c7bb6", alpha = 0.6,
               outlier.colour = "black", outlier.size = 1) +
  coord_flip() +
  labs(x = NULL, y = "Browser biomass (g)",
       title = "Browser biomass distribution by site") +
  theme_bw(base_size = 9)

# ── Zeros by site ─────────────────────────────────────────────
transect_model_data %>%
  group_by(site) %>%
  summarise(
    n_transects  = n(),
    prop_zeros   = mean(transect_browser_biomass == 0),
    mean_biomass = mean(transect_browser_biomass),
    .groups = "drop"
  ) %>%
  arrange(desc(prop_zeros)) %>%
  print(n = Inf)

# ── Family selection ──────────────────────────────────────────

# F1: Gaussian on log(y + 0.01)
mF1_gaussian <- glmmTMB(
  log_browser_biomass ~ sst_sc + log_chla_sc +
    log_settlement_grav_sc  + rugosity_sc + (1 | site),
  family = gaussian(),
  data   = transect_model_data
)

resF1 <- simulateResiduals(mF1_gaussian, n = 1000)

jpeg("dharma_browser_F1_gaussian.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF1, main = "DHARMa — Gaussian on log(y + 0.01)"); dev.off()

plot(resF1)
testZeroInflation(resF1)
testDispersion(resF1)

# F2: Plain Tweedie
mF2_tweedie <- glmmTMB(
  transect_browser_biomass ~ sst_sc + log_chla_sc +
    log_settlement_grav_sc  + rugosity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = transect_model_data
)

resF2 <- simulateResiduals(mF2_tweedie, n = 1000)

jpeg("dharma_browser_F2_tweedie.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF2, main = "DHARMa — Tweedie"); dev.off()

plot(resF2)
testZeroInflation(resF2)
testDispersion(resF2)

# F3: Zero-inflated Tweedie
mF3_tweedie_zi <- glmmTMB(
  transect_browser_biomass ~ sst_sc + log_chla_sc +
    log_settlement_grav_sc  + rugosity_sc + (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = transect_model_data
)

resF3 <- simulateResiduals(mF3_tweedie_zi, n = 1000)

jpeg("dharma_browser_F3_tweedie_zi.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF3, main = "DHARMa — Zero-inflated Tweedie"); dev.off()

plot(resF3)
testZeroInflation(resF3)
testDispersion(resF3)

cat("\n--- Family selection: transect-level browser biomass ---\n")
print(make_aicc_df(list(
  "Tweedie"    = mF2_tweedie,
  "ZI Tweedie" = mF3_tweedie_zi
)))

# ── Family selection decision ─────────────────────────────────
# Transect-level browser biomass shows strong right-skew and a high
# proportion of zeros (~45%), consistent with a compound distribution.

# The plain Tweedie model (F2):
# - No evidence of zero-inflation (DHARMa p = 0.724)
# - No dispersion issues
# - Adequate residual behaviour

# The zero-inflated Tweedie (F3):
# - Also shows no evidence of excess zeros (p = 0.79)
# - Slightly lower AICc (ΔAICc = 1.1), but improvement is < 2

# Decision rule:
# Adopt ZI Tweedie only if:
#   (1) ZI test is significant, AND
#   (2) AICc improves by > 2

# Neither condition is met.

# Final decision:
# Retain plain Tweedie for consistency and parsimony.
# Note: ZI Tweedie has marginal AICc support (ΔAICc = 1.1),
# but this is insufficient to justify added complexity.

# ── Random effect structure selection ────────────────────────
re_t_null   <- glmmTMB(transect_browser_biomass ~ sst_sc + log_chla_sc +
                         log_settlement_grav_sc  + rugosity_sc,
                       family = tweedie(link = "log"), data = transect_model_data)

re_t_site   <- glmmTMB(transect_browser_biomass ~ sst_sc + log_chla_sc +
                         log_settlement_grav_sc  + rugosity_sc +
                         (1 | site),
                       family = tweedie(link = "log"), data = transect_model_data)

re_t_nested <- glmmTMB(transect_browser_biomass ~ sst_sc + log_chla_sc +
                         log_settlement_grav_sc  + rugosity_sc +
                         (1 | country/site),
                       family = tweedie(link = "log"), data = transect_model_data)

cat("\n--- RE structure comparison (transect-level browser) ---\n")
print(make_aicc_df(list(
  "No RE"              = re_t_null,
  "(1 | site)"         = re_t_site,
  "(1 | country/site)" = re_t_nested
)))

# ── RE structure decision — transect-level browser biomass ────
# (1 | site) has the lowest AICc and is preferred.
# (1 | country/site) is within ΔAICc = 0.98, indicating near-equivalence,
# but does not provide sufficient improvement to justify added complexity.
# No RE is strongly unsupported (ΔAICc = 64.97), confirming substantial
# within-site dependence among transects.
# Proceed with the simpler site-level random intercept (1 | site).
# ── Candidate models ──────────────────────────────────────────

browser_family <- tweedie(link = "log")
# --- Null ---
m0                 <- glmmTMB(transect_browser_biomass ~ 1                                                                                        + (1 | site), family = browser_family, data = transect_model_data)

# --- Single predictor ---
m_env              <- glmmTMB(transect_browser_biomass ~ sst_sc + log_chla_sc                                                                    + (1 | site), family = browser_family, data = transect_model_data)
m_market           <- glmmTMB(transect_browser_biomass ~ log_market_gravity_sc                                                                   + (1 | site), family = browser_family, data = transect_model_data)
m_settgrav         <- glmmTMB(transect_browser_biomass ~ log_settlement_grav_sc                                                                  + (1 | site), family = browser_family, data = transect_model_data)
m_settpop          <- glmmTMB(transect_browser_biomass ~ log_settlement_pop_sc                                                                   + (1 | site), family = browser_family, data = transect_model_data)
m_hab              <- glmmTMB(transect_browser_biomass ~ rugosity_sc                                                                             + (1 | site), family = browser_family, data = transect_model_data)

# --- Environment + human pressure ---
m_env_market       <- glmmTMB(transect_browser_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc                                           + (1 | site), family = browser_family, data = transect_model_data)
m_env_settgrav     <- glmmTMB(transect_browser_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc                                          + (1 | site), family = browser_family, data = transect_model_data)
m_env_settpop      <- glmmTMB(transect_browser_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc                                           + (1 | site), family = browser_family, data = transect_model_data)

# --- Habitat + human pressure ---
m_hab_market       <- glmmTMB(transect_browser_biomass ~ rugosity_sc + log_market_gravity_sc                                                     + (1 | site), family = browser_family, data = transect_model_data)
m_hab_settgrav     <- glmmTMB(transect_browser_biomass ~ rugosity_sc + log_settlement_grav_sc                                                    + (1 | site), family = browser_family, data = transect_model_data)
m_hab_settpop      <- glmmTMB(transect_browser_biomass ~ rugosity_sc + log_settlement_pop_sc                                                     + (1 | site), family = browser_family, data = transect_model_data)

# --- Full (single human pressure metric) ---
m_full_market      <- glmmTMB(transect_browser_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc                            + (1 | site), family = browser_family, data = transect_model_data)
m_full_settgrav    <- glmmTMB(transect_browser_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc                            + (1 | site), family = browser_family, data = transect_model_data)
m_full_settpop     <- glmmTMB(transect_browser_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc  + rugosity_sc                            + (1 | site), family = browser_family, data = transect_model_data)

# --- Combined gravity metrics ---
m_both_grav        <- glmmTMB(transect_browser_biomass ~ log_market_gravity_sc + log_settlement_grav_sc                                         + (1 | site), family = browser_family, data = transect_model_data)
m_hab_both_grav    <- glmmTMB(transect_browser_biomass ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc                           + (1 | site), family = browser_family, data = transect_model_data)
m_env_both_grav    <- glmmTMB(transect_browser_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc                  + (1 | site), family = browser_family, data = transect_model_data)
m_full_both_grav   <- glmmTMB(transect_browser_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc    + (1 | site), family = browser_family, data = transect_model_data)

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

cat("\n--- AICc: Transect-level browser biomass ---\n")
print(make_aicc_df(model_list_transect))

# ── Diagnostics on best model ─────────────────────────────────
# Best model confirmed: Habitat + market gravity (weight = 0.33)
# Also inspect habitat-only (delta 0.82) as robust parsimonious model

res_t_hab_market <- simulateResiduals(m_hab_market, n = 1000)

jpeg("dharma_browser_transect_hab_market.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_t_hab_market, main = "DHARMa — transect browser: Habitat + market gravity"); dev.off()

plot(res_t_hab_market)
testZeroInflation(res_t_hab_market)
testDispersion(res_t_hab_market)
testOutliers(res_t_hab_market)
plotResiduals(res_t_hab_market, transect_model_data$rugosity_sc,           xlab = "Rugosity")
plotResiduals(res_t_hab_market, transect_model_data$log_market_gravity_sc, xlab = "Market gravity")

res_t_hab <- simulateResiduals(m_hab, n = 1000)

jpeg("dharma_browser_transect_hab.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_t_hab, main = "DHARMa — transect browser: Habitat only"); dev.off()

plot(res_t_hab)
testZeroInflation(res_t_hab)
testDispersion(res_t_hab)
testOutliers(res_t_hab)

summary(m_hab_market)
summary(m_hab)

# ── Convergence with site-level result ────────────────────────
# Site level s1 (no RE):    Habitat + market gravity — weight = 0.33
# Site level s2 (country):  Habitat + settlement pop. — weight = 0.26
#                           (market gravity third, delta 0.49)
# Transect level:           Habitat + market gravity — weight = 0.33
#
# Transect result converges with s1 rather than s2 — consistent
# with country RE absorbing some human pressure variation at site
# level that manifests as settlement pop. preference. Market gravity
# remains the dominant human pressure signal at transect level.
# Rugosity robust across all three levels.
# Market gravity positive direction — interpret with caution

# ── Marginal effect plots ─────────────────────────────────────
( p_t_rugosity <- plot_effect(m_hab_market,
                              transect_model_data,
                              "rugosity_sc",
                              "Rugosity (scaled)",
                              y_label = "Browser biomass (g)") )

( p_t_market <- plot_effect(m_hab_market,
                            transect_model_data,
                            "log_market_gravity_sc",
                            "Market gravity (scaled)",
                            y_label = "Browser biomass (g)") )

jpeg("transect_browser_marginal_effects.jpg", width = 22, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(p_t_rugosity, p_t_market, ncol = 2)
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
#  Response: Total browser count per transect (integer >= 0)
#  Family:   Poisson → NB2 → NB1; selected via AICc
#  Random fx: (1 | site)
# ============================================================

# ── Explore count distribution ────────────────────────────────
cat("Transects (count data):", nrow(transect_model_data), "\n")
cat("Zeros in count data:",    sum(transect_model_data$transect_browser_count == 0), "\n")
cat("Proportion zeros:",       round(mean(transect_model_data$transect_browser_count == 0), 3), "\n")

summary(transect_model_data$transect_browser_count)

ggplot(transect_model_data, aes(x = transect_browser_count)) +
  geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
  labs(x = "Total browser count per transect", y = "Frequency",
       title = "Raw browser count distribution") +
  theme_bw()

# ── Mean-variance relationship ────────────────────────────────
transect_model_data %>%
  group_by(site) %>%
  summarise(mean_count = mean(transect_browser_count),
            var_count  = var(transect_browser_count),
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
  transect_browser_count ~ sst_sc + log_chla_sc +
    log_settlement_grav_sc + rugosity_sc + (1 | site),
  family = poisson(link = "log"),
  data   = transect_model_data
)

resC1 <- simulateResiduals(mC1_poisson, n = 1000)

jpeg("dharma_browser_C1_poisson.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC1, main = "DHARMa — Poisson"); dev.off()

plot(resC1)
testDispersion(resC1)
testZeroInflation(resC1)
testOutliers(resC1)

# C2: NB2 — quadratic variance
mC2_nb2 <- glmmTMB(
  transect_browser_count ~ sst_sc + log_chla_sc +
    log_settlement_grav_sc + rugosity_sc + (1 | site),
  family = nbinom2(link = "log"),
  data   = transect_model_data
)

resC2 <- simulateResiduals(mC2_nb2, n = 1000)

jpeg("dharma_browser_C2_nb2.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC2, main = "DHARMa — NB2"); dev.off()

plot(resC2)
testDispersion(resC2)
testZeroInflation(resC2)
testOutliers(resC2)

# C3: NB1 — linear variance
mC3_nb1 <- glmmTMB(
  transect_browser_count ~ sst_sc + log_chla_sc +
    log_settlement_grav_sc + rugosity_sc + (1 | site),
  family = nbinom1(link = "log"),
  data   = transect_model_data
)

resC3 <- simulateResiduals(mC3_nb1, n = 1000)

jpeg("dharma_browser_C3_nb1.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC3, main = "DHARMa — NB1"); dev.off()

plot(resC3)
testDispersion(resC3)
testZeroInflation(resC3)
testOutliers(resC3)

cat("\n--- Family selection: browser count models ---\n")
print(make_aicc_df(list(
  "Poisson" = mC1_poisson,
  "NB2"     = mC2_nb2,
  "NB1"     = mC3_nb1
)))

# ── Family selection ──────────────────────────────────────────
# NB1 is overwhelmingly supported (weight = 0.999).
# NB2 is substantially worse (ΔAICc = 13.32),
# and Poisson is strongly unsupported (ΔAICc = 283.73).
# Poisson also shows significant zero inflation (p < 0.001),
# confirming it is inappropriate. Residual diagnostics for NB1
# show no dispersion issues (p = 0.828), no zero inflation
# (p = 0.712), and no outlier concerns.

# ── Random effect structure selection ────────────────────────
# Using selected count family (NB1).
# Compare no RE, site-level RE, and nested country/site RE.

count_family <- nbinom1(link = "log")

re_c_null   <- glmmTMB(transect_browser_count ~ sst_sc + log_chla_sc +
                         log_settlement_grav_sc + rugosity_sc,
                       family = count_family, data = transect_model_data)

re_c_site   <- glmmTMB(transect_browser_count ~ sst_sc + log_chla_sc +
                         log_settlement_grav_sc + rugosity_sc +
                         (1 | site),
                       family = count_family, data = transect_model_data)

re_c_nested <- glmmTMB(transect_browser_count ~ sst_sc + log_chla_sc +
                         log_settlement_grav_sc + rugosity_sc +
                         (1 | country/site),
                       family = count_family, data = transect_model_data)

cat("\n--- RE structure comparison (browser counts) ---\n")
print(make_aicc_df(list(
  "No RE"              = re_c_null,
  "(1 | site)"         = re_c_site,
  "(1 | country/site)" = re_c_nested
)))

# ── RE structure decision — browser counts ────────────────────
# (1 | site) has the lowest AICc and is preferred.
# (1 | country/site) is within ΔAICc = 1.13, indicating near-equivalence,
# but does not justify additional complexity.
# No RE is worse (ΔAICc = 4.10), indicating some within-site dependence,
# though weaker than in biomass models.
# Proceed with the simpler site-level random intercept (1 | site).

# ── Candidate models ──────────────────────────────────────────

# --- Null ---
cm0                <- glmmTMB(transect_browser_count ~ 1                                                                                        + (1 | site), family = count_family, data = transect_model_data)

# --- Single predictor ---
cm_env             <- glmmTMB(transect_browser_count ~ sst_sc + log_chla_sc                                                                    + (1 | site), family = count_family, data = transect_model_data)
cm_market          <- glmmTMB(transect_browser_count ~ log_market_gravity_sc                                                                   + (1 | site), family = count_family, data = transect_model_data)
cm_settgrav        <- glmmTMB(transect_browser_count ~ log_settlement_grav_sc                                                                  + (1 | site), family = count_family, data = transect_model_data)
cm_settpop         <- glmmTMB(transect_browser_count ~ log_settlement_pop_sc                                                                   + (1 | site), family = count_family, data = transect_model_data)
cm_hab             <- glmmTMB(transect_browser_count ~ rugosity_sc                                                                             + (1 | site), family = count_family, data = transect_model_data)

# --- Environment + human pressure ---
cm_env_mkt         <- glmmTMB(transect_browser_count ~ sst_sc + log_chla_sc + log_market_gravity_sc                                           + (1 | site), family = count_family, data = transect_model_data)
cm_env_settgrav    <- glmmTMB(transect_browser_count ~ sst_sc + log_chla_sc + log_settlement_grav_sc                                          + (1 | site), family = count_family, data = transect_model_data)
cm_env_settpop     <- glmmTMB(transect_browser_count ~ sst_sc + log_chla_sc + log_settlement_pop_sc                                           + (1 | site), family = count_family, data = transect_model_data)

# --- Habitat + human pressure ---
cm_hab_market      <- glmmTMB(transect_browser_count ~ rugosity_sc + log_market_gravity_sc                                                     + (1 | site), family = count_family, data = transect_model_data)
cm_hab_settgrav    <- glmmTMB(transect_browser_count ~ rugosity_sc + log_settlement_grav_sc                                                    + (1 | site), family = count_family, data = transect_model_data)
cm_hab_settpop     <- glmmTMB(transect_browser_count ~ rugosity_sc + log_settlement_pop_sc                                                     + (1 | site), family = count_family, data = transect_model_data)

# --- Full (single human pressure metric) ---
cm_full_mkt        <- glmmTMB(transect_browser_count ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc                            + (1 | site), family = count_family, data = transect_model_data)
cm_full_settgrav   <- glmmTMB(transect_browser_count ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc                            + (1 | site), family = count_family, data = transect_model_data)
cm_full_settpop    <- glmmTMB(transect_browser_count ~ sst_sc + log_chla_sc + log_settlement_pop_sc  + rugosity_sc                            + (1 | site), family = count_family, data = transect_model_data)

# --- Combined gravity metrics ---
cm_both_grav       <- glmmTMB(transect_browser_count ~ log_market_gravity_sc + log_settlement_grav_sc                                         + (1 | site), family = count_family, data = transect_model_data)
cm_hab_both_grav   <- glmmTMB(transect_browser_count ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc                           + (1 | site), family = count_family, data = transect_model_data)
cm_env_both_grav   <- glmmTMB(transect_browser_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc                  + (1 | site), family = count_family, data = transect_model_data)
cm_full_both_grav  <- glmmTMB(transect_browser_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc    + (1 | site), family = count_family, data = transect_model_data)

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

cat("\n--- AICc: browser count models ---\n")
print(make_aicc_df(model_list_counts))

# ── Model selection: browser counts ───────────────────────────
# The best-supported model is Habitat + both gravity (AICc weight = 0.35).
# Habitat + market gravity is nearly equivalent (ΔAICc = 0.87),
# so the simpler model is acceptable for interpretation.
# Models without habitat are clearly worse, indicating that rugosity
# contributes meaningfully to browser abundance once human pressure
# is considered.
#
# For diagnostics and effect plots, use the near-equivalent simpler model:
# Habitat + market gravity.

# ── Diagnostics ───────────────────────────────────────────────
# DHARMa diagnostics on the near-equivalent simpler model
# (Habitat + market gravity) show no dispersion issues (p = 0.866),
# no zero inflation (p = 0.800), and no outliers.
# Rugosity (IRR = 1.24, p = 0.035) and market gravity
# (IRR = 1.42, p < 0.001) are both significant.

res_cm_hab_market <- simulateResiduals(cm_hab_market, n = 1000)

jpeg("dharma_browser_count_hab_market.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_cm_hab_market, main = "DHARMa — browser counts: Habitat + market gravity"); dev.off()

plot(res_cm_hab_market)
testDispersion(res_cm_hab_market)
testZeroInflation(res_cm_hab_market)
testOutliers(res_cm_hab_market)

plotResiduals(res_cm_hab_market, transect_model_data$rugosity_sc,           xlab = "Rugosity")
plotResiduals(res_cm_hab_market, transect_model_data$log_market_gravity_sc, xlab = "Market gravity")

summary(cm_hab_market)
exp(fixef(cm_hab_market)$cond)

# ── Marginal effect plots ─────────────────────────────────────
# # Use the simpler near-equivalent model for interpretation.
p_count_market <- plot_effect(cm_hab_market, transect_model_data,
                              "log_market_gravity_sc",
                              "Market gravity (scaled)",
                              "Expected browser count",
                              colour = "#d7191c")

p_count_rug    <- plot_effect(cm_hab_market, transect_model_data,
                              "rugosity_sc",
                              "Rugosity (scaled)",
                              "Expected browser count",
                              colour = "#d7191c")

jpeg("browser_count_marginal_effects.jpg", width = 22, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(p_count_market, p_count_rug, ncol = 2)
dev.off()

# ============================================================
#  SYNTHESIS: BIOMASS vs COUNT MODEL CONCLUSIONS
#
#  Part 1 best model:   s2_m_hab_settpop (Tweedie, 1 | country)
#                       — Habitat + settlement pop. (weight ≈ 0.26)
#                       — Habitat alone near-equivalent (ΔAICc ≈ 0.7)
#
#  Part 2 best model:   m_hab_market (Tweedie, 1 | site)
#                       — Habitat + market gravity (weight ≈ 0.33)
#                       — Habitat alone near-equivalent (ΔAICc ≈ 0.8)
#
#  Part 3 best model:   cm_hab_both_grav (NB1, 1 | site)
#                       — Habitat + both gravity (weight ≈ 0.35)
#                       — Habitat + market gravity near-equivalent (ΔAICc = 0.87)
#
#  Convergence across Parts 1 and 2:
#  Rugosity is the dominant predictor of browser biomass at both
#  site and transect scales. Human pressure metrics add weak,
#  model-dependent signal.
#
#  Part 3 nuance:
#  Habitat remains important for abundance when paired with human
#  pressure, but model uncertainty exists among top models.
#
#  Key divergence — biomass vs counts:
#  Rugosity strongly drives biomass (Parts 1 & 2), but its effect
#  on abundance is weaker and less consistent across top models.
#  This suggests habitat complexity primarily influences size
#  structure (biomass) rather than numerical abundance.
#
#  Human pressure metrics:
#  Market gravity appears in Parts 2 and 3, but with a positive
#  direction that is counterintuitive. This effect is not robust
#  (sensitive to metric choice) and likely reflects spatial
#  confounding rather than a causal fishing effect.
#
#  Settlement population appears only in Part 1 and is not retained
#  at finer scales, suggesting any signal operates at broader
#  spatial aggregation.
#
#  Environment (SST, chlorophyll-a):
#  Not competitive in any browser model (ΔAICc > 10 throughout).
#  No evidence for environmental control of browser biomass or
#  abundance in this system.
# ============================================================
# ============================================================

cat("\n=== PART 1 — Site-level browser biomass, best model ===\n")
summary(s2_m_hab_settpop)

cat("\n=== PART 1 — Site-level browser biomass, habitat only ===\n")
summary(s2_m_hab)

cat("\n=== PART 2 — Transect-level browser biomass, best model ===\n")
summary(m_hab_market)

cat("\n=== PART 2 — Transect-level browser biomass, habitat only ===\n")
summary(m_hab)

cat("\n=== PART 3 — Browser count model, near-equivalent simpler model (Habitat + market gravity) ===\n")
summary(cm_hab_market)

# ── Rugosity ─────────────────────────────────────────────────
# Part 1: glmmTMB country RE — use fixef()$cond
# Parts 2 & 3: glmmTMB site RE — use fixef()$cond
rug_site        <- fixef(s2_m_hab)$cond["rugosity_sc"]
rug_site_full   <- fixef(s2_m_hab_settpop)$cond["rugosity_sc"]
rug_biomass     <- fixef(m_hab)$cond["rugosity_sc"]
rug_count       <- fixef(cm_hab_market)$cond["rugosity_sc"]

cat("\n--- Rugosity effect ---\n")
cat("Part 1 — site biomass (habitat only)   beta:", round(rug_site,      3), "\n")
cat("Part 1 — site biomass (+ settpop)      beta:", round(rug_site_full, 3), "\n")
cat("Part 2 — transect biomass              beta:", round(rug_biomass,   3), "\n")
cat("Part 3 — count model                   beta:", round(rug_count,     3),
    " IRR:", round(exp(rug_count), 3), "\n")
# Rugosity present in biomass models (Parts 1 & 2) but weak in
# counts (Part 3, present only via habitat + market gravity model).
# Complex reefs support larger-bodied browsers not simply more
# individuals.

# ── Market gravity ────────────────────────────────────────────
mkt_biomass <- fixef(m_hab_market)$cond["log_market_gravity_sc"]
mkt_count   <- fixef(cm_hab_market)$cond["log_market_gravity_sc"]

cat("\n--- Market gravity effect ---\n")
cat("Part 1 — site biomass:     marginal (delta 0.49 vs habitat only)\n")
cat("Part 2 — transect biomass  beta:", round(mkt_biomass, 3), "\n")
cat("Part 3 — count model       beta:", round(mkt_count,   3),
    " IRR:", round(exp(mkt_count), 3), "\n")
# Positive direction across Parts 2 & 3 — counter intuitive.

# ── Settlement pop. ───────────────────────────────────────────
cat("\n--- Settlement pop. ---\n")
cat("Part 1 — site biomass      beta:",
    round(fixef(s2_m_hab_settpop)$cond["log_settlement_pop_sc"], 3), "\n")
cat("Part 2 — transect biomass: not in best model (delta 2.53)\n")
cat("Part 3 — counts:           not in best model (delta 13.40)\n")
# Settlement pop. only relevant at site level under country RE —
# country-level variation may be absorbing human pressure signal
# that appears as settlement pop. effect at site level.

# ── Chlorophyll-a ─────────────────────────────────────────────
cat("\n--- Chlorophyll-a ---\n")
cat("Parts 1, 2 & 3: not in best model in any analysis\n")
cat("Environment not competitive across all browser analyses (delta > 10)\n")
# No evidence for productivity or turbidity effects on browsers
# in this system.

# ── SST ───────────────────────────────────────────────────────
cat("\n--- SST ---\n")
cat("Parts 1, 2 & 3: not in best model in any analysis\n")
# No support for temperature effect on browser biomass or abundance.