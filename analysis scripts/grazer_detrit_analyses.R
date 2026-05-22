# ============================================================
#  GRAZER / DETRITIVORE BIOMASS & ABUNDANCE — MIXED EFFECTS MODELS
#
#  Trophic group definition:
#    Grazers, detritivores, and grazer-detritivores are treated as
#    a single functional guild. All three process benthic material
#    and contribute to the same ecosystem process (removal of algae
#    and detritus from reef surfaces). Their roles overlap
#    substantially and many species are classified differently across
#    studies depending on gut content methodology. Combining them
#    reduces classification uncertainty and produces a more
#    ecologically coherent and statistically tractable response
#    variable. This approach is consistent with standard practice
#    in Indo-Pacific reef fish ecology.
#
#  Study design:
#    Transects nested within stations, stations within sites,
#    sites within locations, locations within countries.
#
#  Predictors measured at site level (averaged from station):
#    SST, Chl-a, Human gravity (market / settlement), Rugosity
#
#  Zero structure note:
#    No site-level zeros (0/64) and only 2.1% transect-level zeros
#    (6/282), reflecting the ubiquity of this combined functional
#    guild. Family selection differs from other trophic groups:
#    Gaussian on log(y) is preferred at site level (no zeros, no
#    constant needed); Tweedie adopted at transect level to handle
#    the small number of zeros natively.
#
#  Analytical structure:
#    PART 1 — Site-level analysis (PRIMARY)
#              Matches response resolution to predictor resolution.
#              Sites are the true unit of environmental inference.
#              Family ladder: Gaussian raw → Gaussian log → Gamma
#
#    PART 2 — Transect-level biomass (SENSITIVITY CHECK)
#              Retains within-site variation; (1 | site) accounts
#              for non-independence. Confirms site-level findings
#              are not an artefact of averaging.
#              Family ladder: Gaussian log → Tweedie → ZI Tweedie
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

grazer_transects <- fish_2009 %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_grazer_biomass = sum(
      ifelse(trophic_group %in% c("grazers", "detritivores", "grazer-detritivores"),
             tot_wt_g, 0),
      na.rm = TRUE
    ),
    transect_grazer_count = sum(
      ifelse(trophic_group %in% c("grazers", "detritivores", "grazer-detritivores"),
             number, 0),
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

cat("Number of transects:", nrow(grazer_transects), "\n")
cat("Number of sites:",     n_distinct(grazer_transects$site), "\n")
cat("Number of countries:", n_distinct(grazer_transects$country), "\n")

# ==============================================================================
#  BIOMASS DATA EXPLORATION
# ==============================================================================

# ── Aggregate data at site level ─────────────────────────────────────────────
site_data <- grazer_transects %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass = mean(transect_grazer_biomass, na.rm = TRUE),
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
# Expected: 0 — grazers/detritivores are ubiquitous; site means
# will not be zero when averaging across 3+ transects per site.

# ── Raw distribution ──────────────────────────────────────────────────────────
( site_raw <- ggplot(site_data, aes(x = mean_biomass)) +
    geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
    labs(x = "Mean grazer / detritivore biomass per site (g)", y = "Frequency",
         title = "Raw Site-Level Grazer / Detritivore Biomass") +
    theme_bw() )

# ── Box-Cox ───────────────────────────────────────────────────────────────────
# No zeros — boxcox can be applied directly to site means
MASS::boxcox(
  lm(mean_biomass ~ 1, data = site_data),
  lambda = seq(-2, 2, 0.1)
)
# lambda ~ 0 → log transformation appropriate

# ── Apply transformations ─────────────────────────────────────────────────────
# No constant needed — no zeros at site level
site_data <- site_data %>%
  mutate(
    log_mean_biomass  = log(mean_biomass),
    sqrt_mean_biomass = sqrt(mean_biomass)
  )

( site_log <- ggplot(site_data, aes(x = log_mean_biomass)) +
    geom_histogram(bins = 25, fill = "#1a9641", colour = "white") +
    labs(x = "log(mean biomass)", y = "Frequency",
         title = "Log-transformed Site-Level Grazer / Detritivore Biomass") +
    theme_bw() )

( site_sqrt <- ggplot(site_data, aes(x = sqrt_mean_biomass)) +
    geom_histogram(bins = 25, fill = "#d7191c", colour = "white") +
    labs(x = "sqrt(mean biomass)", y = "Frequency",
         title = "Sqrt-transformed Site-Level Grazer / Detritivore Biomass") +
    theme_bw() )

jpeg("site_grazer_biomass_distributions.jpg", width = 33, height = 11,
     units = "cm", res = 300)
gridExtra::grid.arrange(site_raw, site_log, site_sqrt, ncol = 3)
dev.off()

# ── Normality checks ──────────────────────────────────────────────────────────
qqnorm(site_data$log_mean_biomass,
       main = "Q-Q plot: log(mean grazer / detritivore biomass per site)")
qqline(site_data$log_mean_biomass, col = "red")
shapiro.test(site_data$log_mean_biomass)

# Shapiro-Wilk: W = 0.944, p = 0.006 — marginal departure from
# normality on log scale. 
# Inspect Q-Q plot; moderate tail
# behaviour expected given site-level sample size (n = 64).
# Proceed with lm() — OLS is robust to mild non-normality
# at this sample size.

# ── Variation by site ─────────────────────────────────────────────────────────
ggplot(site_data, aes(x = reorder(site, mean_biomass, median),
                      y = mean_biomass)) +
  geom_col(fill = "#2c7bb6", alpha = 0.7) +
  coord_flip() +
  labs(x = NULL, y = "Mean grazer / detritivore biomass (g)",
       title = "Mean grazer / detritivore biomass by site (raw)") +
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
  labs(x = NULL, y = "log(mean grazer / detritivore biomass per site)") +
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

jpeg("predictor_distributions_grazer.jpg", width = 33, height = 22, units = "cm", res = 300)
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
# Settlement gravity and settlement pop. both proxy local human
# pressure. Select the better-performing metric via AICc before
# entering the main candidate set.
# NOTE: use lm() here — no zeros at site level, Gaussian appropriate.
# Settlement gravity produced a convergence warning in glmmTMB
# Tweedie — lm() on log(mean_biomass) avoids this issue entirely.

settlement_data <- grazer_transects %>%
  left_join(scaled_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    log_mean_biomass       = log(mean(transect_grazer_biomass, na.rm = TRUE)),
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

# Settlement metric selection: settlement pop. is preferred by AICc.
# Delta AICc = 1.64, so settlement gravity is still a plausible
# near-equivalent alternative. Because the difference is < 2, both
# metrics can be retained as parallel candidate predictors in the
# main model set if you want to keep the comparison explicit.
# If you want a single primary human-pressure metric, use
# settlement pop. going forward.

rm(settlement_data)

# ============================================================
#  ANALYSIS DATASETS
# ============================================================

final_predictors <- scaled_predictors %>%
  dplyr::select(site, log_market_gravity_sc, log_settlement_grav_sc,
                log_settlement_pop_sc, log_chla_sc, sst_sc, rugosity_sc)

# ── Transect-level dataset ────────────────────────────────────
transect_model_data <- grazer_transects %>%
  left_join(final_predictors, by = "site") %>%
  mutate(
    log_grazer_biomass = log(transect_grazer_biomass + 0.01)
    # Small constant for the 2.1% transect-level zeros
  )

cat("\nTransect model data:", nrow(transect_model_data), "rows,",
    n_distinct(transect_model_data$site), "sites\n")
cat("Biomass zeros:", sum(transect_model_data$transect_grazer_biomass == 0), "\n")
cat("Count zeros:",  sum(transect_model_data$transect_grazer_count   == 0), "\n")
# Transect-level zeros: ~2.1% (6/282) — very low, reflecting the
# ubiquity of this combined functional guild. Site-level means are
# all positive because the few zero transects are averaged with
# non-zero transects within the same site.

# ── Site-level dataset ────────────────────────────────────────
total_model_data <- transect_model_data %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_grazer_biomass, na.rm = TRUE),
    log_mean_biomass       = log(mean(transect_grazer_biomass, na.rm = TRUE)),
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
# No site-level zeros — Gaussian on log(y) is the preferred
# starting point.
# Tweedie not needed: zero-handling machinery is unnecessary
# for a strictly positive response. ZI Tweedie not attempted.
#
# Ladder:
#   F1: Gaussian on raw mean biomass   — simplest baseline
#   F2: Gaussian on log(mean biomass)  — preferred; no constant needed
#   F3: Gamma (log link)               — comparable to F1 via AICc;
#                                        assess F2 on diagnostics only

# ── F1: Gaussian on raw mean biomass ─────────────────────────
mS_F1 <- lm(mean_biomass ~ sst_sc + log_chla_sc +
              log_market_gravity_sc + rugosity_sc,
            data = total_model_data)

jpeg("diagnostics_site_grazer_F1_gaussian_raw.jpg",
     width = 25, height = 20, units = "cm", res = 300)
par(mfrow = c(2, 2)); plot(mS_F1); par(mfrow = c(1, 1))
dev.off()

par(mfrow = c(2, 2)); plot(mS_F1); par(mfrow = c(1, 1))

# ── F2: Gaussian on log(mean biomass) ────────────────────────
# No constant needed — no zeros at site level.
mS_F2 <- lm(log_mean_biomass ~ sst_sc + log_chla_sc +
              log_market_gravity_sc + rugosity_sc,
            data = total_model_data)

jpeg("diagnostics_site_grazer_F2_gaussian_log.jpg",
     width = 25, height = 20, units = "cm", res = 300)
par(mfrow = c(2, 2)); plot(mS_F2); par(mfrow = c(1, 1))
dev.off()

par(mfrow = c(2, 2)); plot(mS_F2); par(mfrow = c(1, 1))

# ── F3: Gamma (log link) ──────────────────────────────────────
# Natural fit for continuous positive data.
mS_F3 <- glm(mean_biomass ~ sst_sc + log_chla_sc +
               log_market_gravity_sc + rugosity_sc,
             family = Gamma(link = "log"),
             data   = total_model_data)

jpeg("diagnostics_site_grazer_F3_gamma.jpg",
     width = 25, height = 20, units = "cm", res = 300)
par(mfrow = c(2, 2)); plot(mS_F3); par(mfrow = c(1, 1))
dev.off()

par(mfrow = c(2, 2)); plot(mS_F3); par(mfrow = c(1, 1))

# AICc comparison: F1 vs F3 only (same response: mean_biomass)
# AICc NOT comparable between F2 and F1/F3 (different response)
cat("\n--- Family selection: F1 vs F3 (same response scale) ---\n")
print(make_aicc_df(list(
  "Gaussian (raw)" = mS_F1,
  "Gamma (log)"    = mS_F3
)))
# Assess F2 on diagnostics alone — look for:
#   Residuals vs Fitted: no pattern
#   Q-Q: points on line
#   Scale-Location: flat (constant variance)
#   Residuals vs Leverage: no high-influence points

# ── Family selection decision ─────────────────────────────────
# Gamma (log link) strongly outperforms Gaussian on the raw scale
# (ΔAICc = 119.32), indicating clear non-normality and strong
# right-skew in the untransformed response.
#
# F2 (Gaussian on log scale) is still the preferred working model
# provided diagnostics are satisfactory, as it addresses skewness
# directly.
#
# Decision:
#   - Reject F1 (Gaussian raw)
#   - Retain F2 (Gaussian log) as primary model IF diagnostics are clean
#   - Use F3 (Gamma) as a supporting / sensitivity model if needed
#
# Proceed with F2 

# ── RANDOM EFFECT STRUCTURE SELECTION ────────────────────────
# Using glmmTMB for consistent RE comparison.
# Anchor: full market gravity model.

re_null <- glmmTMB(log_mean_biomass ~ sst_sc + log_chla_sc +
                     log_market_gravity_sc + rugosity_sc,
                   family = gaussian(), data = total_model_data)

re_country <- glmmTMB(log_mean_biomass ~ sst_sc + log_chla_sc +
                        log_market_gravity_sc + rugosity_sc +
                        (1 | country),
                      family = gaussian(), data = total_model_data)

cat("\n--- RE structure comparison (site-level grazer) ---\n")
print(make_aicc_df(list(
  "No RE"         = re_null,
  "(1 | country)" = re_country
)))

# ── RE structure decision ─────────────────────────────────────
# No random effect is marginally preferred (ΔAICc = 0.74,
# weight = 0.59 vs 0.41 for country RE). The difference is
# small and does not clearly favour either structure, but
# parsimony supports no RE.
# Proceed with lm() for candidate models.

# ── CANDIDATE MODELS — SITE LEVEL ────────────────────────────
# Family: Gamma (log link) — strongly preferred over Gaussian raw.
# F2 (Gaussian on log(mean_biomass)) remains a useful working
# scale for diagnostics and interpretation, but the raw-scale AICc
# clearly favours Gamma.
# RE structure: no random effect preferred, so all candidate models
# are fitted with lm() below.

# --- Null ---
s1_m0              <- lm(log_mean_biomass ~ 1,                                                                                     data = total_model_data)

# --- Single predictor ---
s1_m_env           <- lm(log_mean_biomass ~ sst_sc + log_chla_sc,                                                                 data = total_model_data)
s1_m_market        <- lm(log_mean_biomass ~ log_market_gravity_sc,                                                                data = total_model_data)
s1_m_settgrav      <- lm(log_mean_biomass ~ log_settlement_grav_sc,                                                               data = total_model_data)
s1_m_settpop       <- lm(log_mean_biomass ~ log_settlement_pop_sc,                                                                data = total_model_data)
s1_m_hab           <- lm(log_mean_biomass ~ rugosity_sc,                                                                          data = total_model_data)

# --- Environment + human pressure ---
s1_m_env_mkt       <- lm(log_mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc,                                        data = total_model_data)
s1_m_env_settgrav  <- lm(log_mean_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc,                                       data = total_model_data)
s1_m_env_settpop   <- lm(log_mean_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc,                                        data = total_model_data)

# --- Habitat + human pressure ---
s1_m_hab_market    <- lm(log_mean_biomass ~ rugosity_sc + log_market_gravity_sc,                                                  data = total_model_data)
s1_m_hab_settgrav  <- lm(log_mean_biomass ~ rugosity_sc + log_settlement_grav_sc,                                                 data = total_model_data)
s1_m_hab_settpop   <- lm(log_mean_biomass ~ rugosity_sc + log_settlement_pop_sc,                                                  data = total_model_data)

# --- Full (single human pressure metric) ---
s1_m_full_mkt      <- lm(log_mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc,                         data = total_model_data)
s1_m_full_settgrav <- lm(log_mean_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc,                         data = total_model_data)
s1_m_full_settpop  <- lm(log_mean_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc  + rugosity_sc,                         data = total_model_data)

# --- Combined gravity metrics ---
s1_m_both_grav     <- lm(log_mean_biomass ~ log_market_gravity_sc + log_settlement_grav_sc,                                       data = total_model_data)
s1_m_hab_both_grav <- lm(log_mean_biomass ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc,                         data = total_model_data)
s1_m_env_both_grav <- lm(log_mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc,               data = total_model_data)
s1_m_full_both_grav<- lm(log_mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc,  data = total_model_data)

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

cat("\n--- AICc: Site-level grazer candidate models ---\n")
print(make_aicc_df(model_list_s1))

# ── Residual diagnostics — top models ────────────────────────
# Based on AICc results:
#   - Full (settlement gravity) is the top model (weight = 0.2759)
#   - Habitat is the closest competitor (ΔAICc = 1.18)
#   - Habitat + settlement pop. is also close (ΔAICc = 1.46)
# Diagnostics are therefore run on the top model and the two
# near-equivalent alternatives.

# Top model: Full (settlement gravity)
cat("\n--- Diagnostics: Full (settlement gravity) ---\n")
par(mfrow = c(2, 2))
plot(s1_m_full_settgrav)
par(mfrow = c(1, 1))

# Close competitor: Habitat
cat("\n--- Diagnostics: Habitat ---\n")
par(mfrow = c(2, 2))
plot(s1_m_hab)
par(mfrow = c(1, 1))

# Close competitor: Habitat + settlement pop.
cat("\n--- Diagnostics: Habitat + settlement pop. ---\n")
par(mfrow = c(2, 2))
plot(s1_m_hab_settpop)
par(mfrow = c(1, 1))

# ── Summaries ─────────────────────────────────────────────────
cat("\n--- Summary: Full (settlement gravity) ---\n")
summary(s1_m_full_settgrav)

cat("\n--- Summary: Habitat ---\n")
summary(s1_m_hab)

cat("\n--- Summary: Habitat + settlement pop. ---\n")
summary(s1_m_hab_settpop)

# ── Coefficient comparison across top models ──────────────────
cat("\n--- Rugosity coefficient stability ---\n")
cat("Full (settlement gravity): beta =",
    round(coef(s1_m_full_settgrav)["rugosity_sc"], 3), "\n")
cat("Habitat:                  beta =",
    round(coef(s1_m_hab)["rugosity_sc"], 3), "\n")
cat("Habitat + settlement pop.: beta =",
    round(coef(s1_m_hab_settpop)["rugosity_sc"], 3), "\n")

cat("\n--- Settlement gravity coefficient stability ---\n")
cat("Full (settlement gravity): beta =",
    round(coef(s1_m_full_settgrav)["log_settlement_grav_sc"], 3), "\n")
cat("Habitat + settlement gravity: beta =",
    round(coef(s1_m_hab_settgrav)["log_settlement_grav_sc"], 3), "\n")

# ── Marginal effect plots ─────────────────────────────────────
# Full model effects (primary interpretation)
( p_site_rugosity <- plot_effect(
  s1_m_full_settgrav,
  total_model_data,
  "rugosity_sc",
  "Rugosity (scaled)",
  y_label = "log(mean grazer / detritivore biomass)"
) )

( p_site_settgrav <- plot_effect(
  s1_m_full_settgrav,
  total_model_data,
  "log_settlement_grav_sc",
  "Settlement gravity (scaled)",
  y_label = "log(mean grazer / detritivore biomass)"
) )

# Optional: show the habitat-only competitor for rugosity robustness
( p_site_hab_rugosity <- plot_effect(
  s1_m_hab,
  total_model_data,
  "rugosity_sc",
  "Rugosity (scaled)",
  y_label = "log(mean grazer / detritivore biomass)"
) )

jpeg("site_grazer_marginal_effects.jpg", width = 28, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(p_site_rugosity, p_site_settgrav, p_site_hab_rugosity, ncol = 3)
dev.off()

# ============================================================
#  PART 2 — TRANSECT-LEVEL BIOMASS (SENSITIVITY CHECK)
#
#  Rationale: Retains within-site variation. (1 | site) accounts
#  for non-independence of transects.
#  Confirms site-level findings are not an artefact of collapsing
#  to site means.
#
#  Response:   transect_grazer_biomass — continuous, 2.1% zeros
#  Family ladder:
#    F1: Gaussian on log(y + 0.01) — small constant for zeros
#    F2: Tweedie (log link)        — handles zeros natively
#    F3: ZI Tweedie                — only if F2 flags excess zeros
#  Random fx:  (1 | site), then (1 | country/site)
# ============================================================

# ── Explore transect-level response ──────────────────────────
summary(grazer_transects$transect_grazer_biomass)

zeros <- mean(grazer_transects$transect_grazer_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")
# ~2.1% zeros — very low; ZI Tweedie unlikely to be warranted

( grazer_raw <- ggplot(grazer_transects, aes(x = transect_grazer_biomass)) +
    geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
    labs(x = "Grazer / detritivore biomass per transect (g)", y = "Frequency",
         title = "Raw Grazer / Detritivore Biomass") +
    theme_bw() )

grazer_transects <- grazer_transects %>%
  mutate(
    log_grazer_biomass  = log(transect_grazer_biomass + 0.01),
    sqrt_grazer_biomass = sqrt(transect_grazer_biomass)
  )

( grazer_log <- ggplot(grazer_transects, aes(x = log_grazer_biomass)) +
    geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
    labs(x = "log(biomass + 0.01)", y = "Frequency",
         title = "Log-transformed Grazer / Detritivore Biomass") +
    theme_bw() )

( grazer_sqrt <- ggplot(grazer_transects, aes(x = sqrt_grazer_biomass)) +
    geom_histogram(bins = 30, fill = "#d7191c", colour = "white") +
    labs(x = "sqrt(biomass)", y = "Frequency",
         title = "Sqrt-transformed Grazer / Detritivore Biomass") +
    theme_bw() )

jpeg("grazer_biomass_distributions.jpg", width = 33, height = 11,
     units = "cm", res = 300)
gridExtra::grid.arrange(grazer_raw, grazer_log, grazer_sqrt, ncol = 3)
dev.off()

# ── Box-Cox on non-zero values ────────────────────────────────
grazer_nonzero <- grazer_transects %>% filter(transect_grazer_biomass > 0)

MASS::boxcox(
  lm(transect_grazer_biomass ~ 1, data = grazer_nonzero),
  lambda = seq(-2, 2, 0.1)
)

# ── Variation by site ─────────────────────────────────────────
ggplot(transect_model_data,
       aes(x = reorder(site, transect_grazer_biomass, median),
           y = transect_grazer_biomass)) +
  geom_boxplot(fill = "#2c7bb6", alpha = 0.6,
               outlier.colour = "black", outlier.size = 1) +
  coord_flip() +
  labs(x = NULL, y = "Grazer / detritivore biomass (g)",
       title = "Grazer / detritivore biomass distribution by site") +
  theme_bw(base_size = 9)

# ── Zeros by site ─────────────────────────────────────────────
transect_model_data %>%
  group_by(site) %>%
  summarise(
    n_transects  = n(),
    prop_zeros   = mean(transect_grazer_biomass == 0),
    mean_biomass = mean(transect_grazer_biomass),
    .groups = "drop"
  ) %>%
  arrange(desc(prop_zeros)) %>%
  print(n = Inf)

# ── Family selection ──────────────────────────────────────────
# AICc not comparable between F1 and F2/F3 (different response).
# Select on DHARMa diagnostics; use AICc only to compare F2 vs F3.

# F1: Gaussian on log(y + 0.01)
mF1_gaussian <- glmmTMB(
  log_grazer_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = gaussian(),
  data   = transect_model_data
)

resF1 <- simulateResiduals(mF1_gaussian, n = 1000)

jpeg("dharma_grazer_F1_gaussian.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF1, main = "DHARMa — Gaussian on log(y + 0.01)"); dev.off()

plot(resF1)
testZeroInflation(resF1)
testDispersion(resF1)

# F2: Plain Tweedie
mF2_tweedie <- glmmTMB(
  transect_grazer_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = transect_model_data
)

resF2 <- simulateResiduals(mF2_tweedie, n = 1000)

jpeg("dharma_grazer_F2_tweedie.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF2, main = "DHARMa — Tweedie"); dev.off()

plot(resF2)
testZeroInflation(resF2)
testDispersion(resF2)

# F3: Zero-inflated Tweedie
# Not fitted.
# With only 2.1% zeros, there is no empirical justification for a
# zero-inflation component. Any ZI structure would be weakly
# identifiable and adds unnecessary computational cost.

# ── Family selection decision ─────────────────────────────────
# F1 (Gaussian on log scale) shows clean diagnostics (no
# dispersion or zero inflation issues), but requires an arbitrary
# +0.01 constant to handle the 2.1% transect-level zeros.
#
# F2 (Plain Tweedie) shows significant overdispersion
# (dispersion = 3.78, p = 0.014) but handles zeros natively
# without a constant. Zero inflation is not significant (p = 0.10).
#
# Despite the Tweedie dispersion flag, Tweedie is retained because:
#   (1) It handles zeros natively — no arbitrary constant needed
#   (2) The dispersion issue is consistent across models and will
#       be noted as a caveat in interpretation
#   (3) AICc-based model comparison remains valid within family
#
# ZI Tweedie not considered — zero inflation not significant.

# ── Random effect structure selection ────────────────────────
re_t_null   <- glmmTMB(transect_grazer_biomass ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc,
                       family = tweedie(link = "log"), data = transect_model_data)

re_t_site   <- glmmTMB(transect_grazer_biomass ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | site),
                       family = tweedie(link = "log"), data = transect_model_data)

re_t_nested <- glmmTMB(transect_grazer_biomass ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | country/site),
                       family = tweedie(link = "log"), data = transect_model_data)

cat("\n--- RE structure comparison (transect-level grazer) ---\n")
print(make_aicc_df(list(
  "No RE"              = re_t_null,
  "(1 | site)"         = re_t_site,
  "(1 | country/site)" = re_t_nested
)))

# ── Random effect structure decision ─────────────────────────
# (1 | site) is preferred (weight = 0.60 vs 0.40 for nested).
# No RE is strongly unsupported (ΔAICc = 86.53).
# Retain (1 | site).
#
# Both random-effect structures are plausible, but the nested
# country effect adds little explanatory power relative to its
# additional complexity.
#
# Decision:
#   - Retain (1 | site) as the primary random effect
#   - Do not include country-level nesting
#
# Interpretation:
#   - Strong within-site clustering of transects
#   - Minimal additional variation attributable to country once
#     site-level structure is accounted for

# ── Candidate models ──────────────────────────────────────────
# Family: Tweedie (log link) selected.
# RE structure: (1 | site) retained as the primary random effect.
# The nested country effect was not retained because it improved
# AICc only marginally relative to the added complexity.

grazer_family <- tweedie(link = "log")

# --- Null ---
m0                 <- glmmTMB(transect_grazer_biomass ~ 1                                                                                        + (1 | site), family = grazer_family, data = transect_model_data)

# --- Single predictor ---
m_env              <- glmmTMB(transect_grazer_biomass ~ sst_sc + log_chla_sc                                                                    + (1 | site), family = grazer_family, data = transect_model_data)
m_market           <- glmmTMB(transect_grazer_biomass ~ log_market_gravity_sc                                                                   + (1 | site), family = grazer_family, data = transect_model_data)
m_settgrav         <- glmmTMB(transect_grazer_biomass ~ log_settlement_grav_sc                                                                  + (1 | site), family = grazer_family, data = transect_model_data)
m_settpop          <- glmmTMB(transect_grazer_biomass ~ log_settlement_pop_sc                                                                   + (1 | site), family = grazer_family, data = transect_model_data)
m_hab              <- glmmTMB(transect_grazer_biomass ~ rugosity_sc                                                                             + (1 | site), family = grazer_family, data = transect_model_data)

# --- Environment + human pressure ---
m_env_market       <- glmmTMB(transect_grazer_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc                                           + (1 | site), family = grazer_family, data = transect_model_data)
m_env_settgrav     <- glmmTMB(transect_grazer_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc                                          + (1 | site), family = grazer_family, data = transect_model_data)
m_env_settpop      <- glmmTMB(transect_grazer_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc                                           + (1 | site), family = grazer_family, data = transect_model_data)

# --- Habitat + human pressure ---
m_hab_market       <- glmmTMB(transect_grazer_biomass ~ rugosity_sc + log_market_gravity_sc                                                     + (1 | site), family = grazer_family, data = transect_model_data)
m_hab_settgrav     <- glmmTMB(transect_grazer_biomass ~ rugosity_sc + log_settlement_grav_sc                                                    + (1 | site), family = grazer_family, data = transect_model_data)
m_hab_settpop      <- glmmTMB(transect_grazer_biomass ~ rugosity_sc + log_settlement_pop_sc                                                     + (1 | site), family = grazer_family, data = transect_model_data)

# --- Full (single human pressure metric) ---
m_full_market      <- glmmTMB(transect_grazer_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc                            + (1 | site), family = grazer_family, data = transect_model_data)
m_full_settgrav    <- glmmTMB(transect_grazer_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc                            + (1 | site), family = grazer_family, data = transect_model_data)
m_full_settpop     <- glmmTMB(transect_grazer_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc  + rugosity_sc                            + (1 | site), family = grazer_family, data = transect_model_data)

# --- Combined gravity metrics ---
m_both_grav        <- glmmTMB(transect_grazer_biomass ~ log_market_gravity_sc + log_settlement_grav_sc                                         + (1 | site), family = grazer_family, data = transect_model_data)
m_hab_both_grav    <- glmmTMB(transect_grazer_biomass ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc                           + (1 | site), family = grazer_family, data = transect_model_data)
m_env_both_grav    <- glmmTMB(transect_grazer_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc                  + (1 | site), family = grazer_family, data = transect_model_data)
m_full_both_grav   <- glmmTMB(transect_grazer_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc    + (1 | site), family = grazer_family, data = transect_model_data)

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

cat("\n--- AICc: Transect-level grazer biomass ---\n")
print(make_aicc_df(model_list_transect))

## ── Diagnostics on best model ─────────────────────────────────
# AICc results confirmed:
#   - Full (settlement gravity): top model (weight = 0.3511)
#   - Full (both gravity): ΔAICc = 1.98
#   - Full (settlement pop.): ΔAICc = 2.23
#
# Both top models show significant overdispersion
# (Full settgrav: dispersion = 3.91, p = 0.010;
#  Full both grav: dispersion = 3.96, p = 0.012).
# Outlier test significant for Full (both gravity) (p = 0.0003,
# 5 outliers) but borderline for Full (settlement gravity)
# (p = 0.020, 3 outliers).
# Zero inflation not significant in either model.
#
# Overdispersion is a known limitation of Tweedie at the
# transect level for this guild. Results are interpreted
# cautiously. Convergence with Part 1 (same top model) is
# reassuring, but coefficient uncertainty may be understated.

# Top model: Full (settlement gravity)
cat("\n--- Diagnostics: Full (settlement gravity) ---\n")
res_best <- simulateResiduals(m_full_settgrav, n = 1000)

jpeg("dharma_grazer_transect_full_settgrav.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_best, main = "DHARMa — Full (settlement gravity)"); dev.off()

plot(res_best)
testZeroInflation(res_best)
testDispersion(res_best)
testOutliers(res_best)

# Close competitor: Full (both gravity)
cat("\n--- Diagnostics: Full (both gravity) ---\n")
res_both <- simulateResiduals(m_full_both_grav, n = 1000)

jpeg("dharma_grazer_transect_full_bothgrav.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_both, main = "DHARMa — Full (both gravity)"); dev.off()

plot(res_both)
testZeroInflation(res_both)
testDispersion(res_both)
testOutliers(res_both)

# ── Convergence with site-level result ────────────────────────
# Strong convergence with site-level analysis:
#   - The same top model structure is selected: Full (settlement gravity)
#   - Confirms that site-level results are not an artefact of averaging
#
# Interpretation:
#   - Settlement gravity is consistently the dominant human-pressure
#     predictor across both spatial resolutions
#   - Rugosity remains an important co-predictor in top models

# ── Marginal effect plots ─────────────────────────────────────
# Effects from the top transect-level model

( p_transect_rugosity <- plot_effect(
  m_full_settgrav,
  transect_model_data,
  "rugosity_sc",
  "Rugosity (scaled)",
  y_label = "Grazer / detritivore biomass (g)"
) )

( p_transect_settgrav <- plot_effect(
  m_full_settgrav,
  transect_model_data,
  "log_settlement_grav_sc",
  "Settlement gravity (scaled)",
  y_label = "Grazer / detritivore biomass (g)"
) )

jpeg("transect_grazer_marginal_effects.jpg", width = 22, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(p_transect_rugosity, p_transect_settgrav, ncol = 2)
dev.off()

# ============================================================
#  PART 3 — TRANSECT-LEVEL COUNTS (COMPLEMENTARY ANALYSIS)
#
#  Response: Total grazer / detritivore count per transect
#  Family:   Poisson → NB2 → NB1; selected via AICc
#  Random fx: (1 | site)
# ============================================================

# ── Explore count distribution ────────────────────────────────
cat("Transects (count data):", nrow(transect_model_data), "\n")
cat("Zeros in count data:",    sum(transect_model_data$transect_grazer_count == 0), "\n")
cat("Proportion zeros:",       round(mean(transect_model_data$transect_grazer_count == 0), 3), "\n")

summary(transect_model_data$transect_grazer_count)

ggplot(transect_model_data, aes(x = transect_grazer_count)) +
  geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
  labs(x = "Total grazer / detritivore count per transect", y = "Frequency",
       title = "Raw grazer / detritivore count distribution") +
  theme_bw()

# ── Mean-variance relationship ────────────────────────────────
transect_model_data %>%
  group_by(site) %>%
  summarise(mean_count = mean(transect_grazer_count),
            var_count  = var(transect_grazer_count),
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
  transect_grazer_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = poisson(link = "log"),
  data   = transect_model_data
)

resC1 <- simulateResiduals(mC1_poisson, n = 1000)

jpeg("dharma_grazer_C1_poisson.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC1, main = "DHARMa — Poisson"); dev.off()

plot(resC1)
testDispersion(resC1)
testZeroInflation(resC1)
testOutliers(resC1)

mC2_nb2 <- glmmTMB(
  transect_grazer_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = nbinom2(link = "log"),
  data   = transect_model_data
)

resC2 <- simulateResiduals(mC2_nb2, n = 1000)

jpeg("dharma_grazer_C2_nb2.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC2, main = "DHARMa — NB2"); dev.off()

plot(resC2)
testDispersion(resC2)
testZeroInflation(resC2)
testOutliers(resC2)

mC3_nb1 <- glmmTMB(
  transect_grazer_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = nbinom1(link = "log"),
  data   = transect_model_data
)

resC3 <- simulateResiduals(mC3_nb1, n = 1000)

jpeg("dharma_grazer_C3_nb1.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC3, main = "DHARMa — NB1"); dev.off()

plot(resC3)
testDispersion(resC3)
testZeroInflation(resC3)
testOutliers(resC3)

cat("\n--- Family selection: grazer count models ---\n")
print(make_aicc_df(list(
  "Poisson" = mC1_poisson,
  "NB2"     = mC2_nb2,
  "NB1"     = mC3_nb1
)))

# ── Family selection decision ─────────────────────────────────
# NB1 is overwhelmingly supported by AICc (weight = 0.9998),
# with NB2 far worse (ΔAICc = 17.21) and Poisson completely
# unsupported (ΔAICc = 868.94).
#
# However, NB1 shows significant zero inflation (p = 0.012),
# suggesting the model slightly underestimates the frequency
# of zero counts. NB2 also shows significant zero inflation
# (p < 0.001), but NB1 handles it considerably better.
# NB1 also produced a step failure warning during DHARMa residual
# smoothing (from mgcv, not from glmmTMB fitting) — check
# diagnostic plots carefully.
#
# NB2 shows highly significant zero inflation (p < 0.001),
# confirming NB1 is the better choice despite the caveat.
# Poisson shows extreme zero inflation (ratioObsSim = 375,
# p < 0.001) and significant outliers — clearly inappropriate.
#
# Decision: retain NB1. The residual zero inflation is noted
# as a caveat but does not invalidate the model comparison.

# ── Random effect structure selection ────────────────────────
count_family <- nbinom1(link = "log")

re_c_null   <- glmmTMB(transect_grazer_count ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc,
                       family = count_family, data = transect_model_data)

re_c_site   <- glmmTMB(transect_grazer_count ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | site),
                       family = count_family, data = transect_model_data)

re_c_nested <- glmmTMB(transect_grazer_count ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | country/site),
                       family = count_family, data = transect_model_data)

cat("\n--- RE structure comparison (grazer counts) ---\n")
print(make_aicc_df(list(
  "No RE"              = re_c_null,
  "(1 | site)"         = re_c_site,
  "(1 | country/site)" = re_c_nested
)))

# ── Random effect structure decision ─────────────────────────
# (1 | site) is clearly preferred (lowest AICc; weight = 0.7426).
# The nested structure (1 | country/site) performs worse
# (ΔAICc = 2.12), indicating limited additional variation at the
# country level once site is accounted for.
# The no-RE model is strongly unsupported (ΔAICc = 39.31).
#
# Decision:
#   - Retain (1 | site) as the random effect structure
#   - Do not include country-level nesting
#
# Interpretation:
#   - Strong clustering of counts within sites
#   - Minimal additional country-level variation beyond site effects


# ── Candidate count models ────────────────────────────────────
# Family: NB1 (log link) — overwhelmingly supported by AICc
# (weight ≈ 1), indicating strong overdispersion with variance
# scaling linearly with the mean.
#
# RE structure: (1 | site) retained — clearly preferred over both
# no RE and nested country structure.
#
# Model set mirrors biomass analyses:
#   - Environment (SST, Chl-a)
#   - Human pressure (market, settlement gravity, settlement pop.)
#   - Habitat (rugosity)
#   - Full combinations and alternative human-pressure formulations
#
# Purpose:
#   - Identify whether predictor effects observed in biomass models
#     are driven by abundance (counts), body size, or both

cm0            <- glmmTMB(transect_grazer_count ~ 1 + (1 | site), family = count_family, data = transect_model_data)
cm_env         <- glmmTMB(transect_grazer_count ~ sst_sc + log_chla_sc + (1 | site), family = count_family, data = transect_model_data)
cm_market      <- glmmTMB(transect_grazer_count ~ log_market_gravity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_settgrav    <- glmmTMB(transect_grazer_count ~ log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_settpop     <- glmmTMB(transect_grazer_count ~ log_settlement_pop_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab         <- glmmTMB(transect_grazer_count ~ rugosity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_env_mkt     <- glmmTMB(transect_grazer_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_env_settgrav<- glmmTMB(transect_grazer_count ~ sst_sc + log_chla_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_env_settpop <- glmmTMB(transect_grazer_count ~ sst_sc + log_chla_sc + log_settlement_pop_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab_market  <- glmmTMB(transect_grazer_count ~ rugosity_sc + log_market_gravity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab_settgrav<- glmmTMB(transect_grazer_count ~ rugosity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab_settpop <- glmmTMB(transect_grazer_count ~ rugosity_sc + log_settlement_pop_sc + (1 | site), family = count_family, data = transect_model_data)
cm_full_mkt    <- glmmTMB(transect_grazer_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_full_settgrav<-glmmTMB(transect_grazer_count ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_full_settpop<- glmmTMB(transect_grazer_count ~ sst_sc + log_chla_sc + log_settlement_pop_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_both_grav   <- glmmTMB(transect_grazer_count ~ log_market_gravity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab_both_grav    <- glmmTMB(transect_grazer_count ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_env_both_grav    <- glmmTMB(transect_grazer_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_full_both_grav   <- glmmTMB(transect_grazer_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)

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
  "Habitat + both gravity"             = cm_hab_both_grav,
  "Env + both gravity"                 = cm_env_both_grav,
  "Full (both gravity)"                = cm_full_both_grav
)

cat("\n--- AICc: grazer count models ---\n")
print(make_aicc_df(model_list_counts))

# ── Diagnostics on best count model ──────────────────────────
# Full (settlement pop.) overwhelmingly supported
# (weight = 0.8469). Next best model (Full settlement gravity)
# is ΔAICc = 4.87 away — clear single-model dominance.
#
# DHARMa diagnostics on Full (settlement pop.):
#   - No dispersion issues (dispersion = 1.026, p = 0.764)
#   - Significant zero inflation (ratioObsSim = 4.31, p = 0.006)
#     — a minor concern; likely reflects the same structural
#     zeros seen in Parts 1 and 2. Interpret cautiously.
#   - No outlier concerns (p = 0.22)

# Top model: Full (settlement pop.)
cat("\n--- Diagnostics: Full (settlement pop.) ---\n")
res_count_best <- simulateResiduals(cm_full_settpop, n = 1000)

jpeg("dharma_grazer_count_full_settpop.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_count_best, main = "DHARMa — Full (settlement pop.)"); dev.off()

plot(res_count_best)
testDispersion(res_count_best)
testZeroInflation(res_count_best)
testOutliers(res_count_best)

# ── Summary ──────────────────────────────────────────────────
cat("\n--- Summary: Full (settlement pop.) ---\n")
summary(cm_full_settpop)

# ── Key coefficient checks ────────────────────────────────────
cat("\n--- Rugosity effect ---\n")
cat("beta =", round(fixef(cm_full_settpop)$cond["rugosity_sc"], 3), "\n")

cat("\n--- Settlement population effect ---\n")
cat("beta =", round(fixef(cm_full_settpop)$cond["log_settlement_pop_sc"], 3), "\n")

# ── Marginal effect plots ─────────────────────────────────────
( p_count_rugosity <- plot_effect(
  cm_full_settpop,
  transect_model_data,
  "rugosity_sc",
  "Rugosity (scaled)",
  y_label = "Grazer / detritivore count"
) )

( p_count_settpop <- plot_effect(
  cm_full_settpop,
  transect_model_data,
  "log_settlement_pop_sc",
  "Settlement population (scaled)",
  y_label = "Grazer / detritivore count"
) )

jpeg("grazer_count_marginal_effects.jpg", width = 22, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(p_count_rugosity, p_count_settpop, ncol = 2)
dev.off()

# ── SYNTHESIS ─────────────────────────────────────────────────
# Part 1 (site biomass):     Full (settlement gravity)
# Part 2 (transect biomass): Full (settlement gravity)
# Part 3 (transect counts):  Full (settlement pop.)
#   Settlement pop. is negative (IRR = 0.87, p = 0.011) —
#   higher settlement pressure associated with fewer grazers.
#   Chl-a also negative (IRR = 0.78, p < 0.001).
#   Rugosity positive (IRR = 1.22, p < 0.001).
#   SST not significant.
#
# Key contrasts to report:
#   - Biomass models and count models do not select the same
#     human-pressure proxy.
#   - Settlement gravity is preferred for biomass, whereas
#     settlement population is preferred for counts.
#   - Rugosity is retained in the top models for both biomass and
#     counts, indicating a robust habitat signal across response
#     types.
#   - Site-level and transect-level biomass analyses are highly
#     consistent: both identify the same top model structure.
#   - The count model suggests that part of the biomass signal
#     likely reflects abundance differences, but the biomass
#     models imply an additional size-structure or biomass-per-
#     individual component tied to settlement gravity.

# ── IRR summary table ─────────────────────────────────────────
cat("\n--- IRR: Full (settlement pop.) ---\n")

irr <- exp(fixef(cm_full_settpop)$cond)
se  <- summary(cm_full_settpop)$coefficients$cond[, "Std. Error"]

cat("SST:               IRR =", round(irr["sst_sc"], 2),
    " (95% CI:", round(exp(log(irr["sst_sc"]) - 1.96*se["sst_sc"]), 2),
    "-",         round(exp(log(irr["sst_sc"]) + 1.96*se["sst_sc"]), 2), ")\n")

cat("Chla:              IRR =", round(irr["log_chla_sc"], 2),
    " (95% CI:", round(exp(log(irr["log_chla_sc"]) - 1.96*se["log_chla_sc"]), 2),
    "-",         round(exp(log(irr["log_chla_sc"]) + 1.96*se["log_chla_sc"]), 2), ")\n")

cat("Settlement pop.:   IRR =", round(irr["log_settlement_pop_sc"], 2),
    " (95% CI:", round(exp(log(irr["log_settlement_pop_sc"]) - 1.96*se["log_settlement_pop_sc"]), 2),
    "-",         round(exp(log(irr["log_settlement_pop_sc"]) + 1.96*se["log_settlement_pop_sc"]), 2), ")\n")

cat("Rugosity:          IRR =", round(irr["rugosity_sc"], 2),
    " (95% CI:", round(exp(log(irr["rugosity_sc"]) - 1.96*se["rugosity_sc"]), 2),
    "-",         round(exp(log(irr["rugosity_sc"]) + 1.96*se["rugosity_sc"]), 2), ")\n")