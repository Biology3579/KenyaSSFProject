# ============================================================
#  CORALLIVORE BIOMASS
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
#              NOTE: zero inflation requires Tweedie family.
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
  
  if (inherits(model, "glmmTMB")) {
    if ("site" %in% names(data))    grid$site <- levels(as.factor(data$site))[1]
    if ("country" %in% names(data)) grid$country <- levels(as.factor(data$country))[1]
  }
  
  pred <- predict(model, newdata = grid, type = "response",
                  se.fit = TRUE, re.form = NA)
  
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

coralliv_transects <- fish_2009 %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_coralliv_biomass = sum(
      ifelse(trophic_group == "corallivores", tot_wt_g, 0),
      na.rm = TRUE
    ),
    transect_coralliv_count = sum(
      ifelse(trophic_group == "corallivores", number, 0),
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

cat("Number of transects:", nrow(coralliv_transects), "\n")
cat("Number of sites:",     n_distinct(coralliv_transects$site), "\n")
cat("Number of countries:", n_distinct(coralliv_transects$country), "\n")

# ==============================================================================
#  BIOMASS DATA EXPLORATION
# ==============================================================================

site_data <- coralliv_transects %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass = mean(transect_coralliv_biomass, na.rm = TRUE),
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
# recorded zero corallivores. Check whether these sites should be retained.

# ── Raw distribution ──────────────────────────────────────────────────────────
( site_raw <- ggplot(site_data, aes(x = mean_biomass)) +
    geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
    labs(x = "Mean corallivore biomass per site (g)", y = "Frequency",
         title = "Raw Site-Level corallivore Biomass") +
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
         title = "Log-transformed Site-Level corallivore Biomass") +
    theme_bw() )

( site_sqrt <- ggplot(site_data, aes(x = sqrt_mean_biomass)) +
    geom_histogram(bins = 25, fill = "#d7191c", colour = "white") +
    labs(x = "sqrt(mean biomass)", y = "Frequency",
         title = "Sqrt-transformed Site-Level corallivore Biomass") +
    theme_bw() )

jpeg("site_coralliv_biomass_distributions.jpg", width = 33, height = 11,
     units = "cm", res = 300)
gridExtra::grid.arrange(site_raw, site_log, site_sqrt, ncol = 3)
dev.off()

# ── Normality checks on log-transformed response ──────────────────────────────
qqnorm(site_data$log_mean_biomass,
       main = "Q-Q plot: log(mean corallivore biomass per site)")
qqline(site_data$log_mean_biomass, col = "red")
shapiro.test(site_data$log_mean_biomass)

# ── Variation by site ─────────────────────────────────────────────────────────
ggplot(site_data, aes(x = reorder(site, mean_biomass, median),
                      y = mean_biomass)) +
  geom_col(fill = "#2c7bb6", alpha = 0.7) +
  coord_flip() +
  labs(x = NULL, y = "Mean corallivore biomass (g)",
       title = "Mean corallivore biomass by site (raw)") +
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
  labs(x = NULL, y = "log(mean corallivore biomass per site + 0.01)") +
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

jpeg("predictor_distributions_coralliv.jpg", width = 33, height = 22, units = "cm", res = 300)
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


# this version oincldues directoion of correlation that might be more useful...
# but do need to change the colorus (swap negative = blue and psotive is red)
# corrplot(corr_matrix,
#          method = "square", type = "lower",
#          tl.col = "black", tl.srt = 0, tl.offset = 0.5,
#          addCoef.col = "black", number.cex = 0.8,
#          col = colorRampPalette(c("#d73027", "white", "#4575b4"))(200),
#          mar = c(0, 0, 4, 2))


# ── Decision guide ────────────────────────────────────────────
# ... (fill out)

# ============================================================
#  CHOOSING SETTLEMENT
# Settlement gravity and settlement pop both proxy local human
# pressure. Select the better-performing metric via AICc before
# entering the main candidate set.
# NOTE: zero-containing site means require glmmTMB Tweedie here.

settlement_data <- coralliv_transects %>%
  left_join(scaled_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    mean_coralliv_biomass     = mean(transect_coralliv_biomass, na.rm = TRUE),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_settlement_pop_sc  = first(log_settlement_pop_sc),
    .groups = "drop"
  )

settgrav <- glmmTMB(mean_coralliv_biomass ~ log_settlement_grav_sc,
                    family = tweedie(link = "log"), data = settlement_data)
settpop  <- glmmTMB(mean_coralliv_biomass ~ log_settlement_pop_sc,
                    family = tweedie(link = "log"), data = settlement_data)

make_aicc_df(list(
  "Settlement gravity" = settgrav,
  "Settlement pop."    = settpop
))

# Settlement metric selection: virtually identical performance
# (delta AICc = 0.48, weights 0.42 vs 0.43). Neither metric
# has any meaningful advantage over the other. Both carried
# forward as parallel candidate model sets throughout all
# subsequent analyses. Results compared across both metrics
# to assess robustness of human pressure effects.)

rm(settlement_data)

# ============================================================
#  ANALYSIS DATASETS
# ============================================================

# Retain all human pressure metrics for parallel model sets
final_predictors <- scaled_predictors %>%
  dplyr::select(site, log_market_gravity_sc, log_settlement_grav_sc, 
                log_settlement_pop_sc, log_chla_sc, sst_sc, rugosity_sc)

# ── Transect-level dataset ────────────────────────────────────
coralliv_model_data <- coralliv_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect model data:", nrow(coralliv_model_data), "rows,",
    n_distinct(coralliv_model_data$site), "sites\n")
cat("Biomass zeros:", sum(coralliv_model_data$transect_coralliv_biomass == 0), "\n")
cat("Count zeros:",  sum(coralliv_model_data$transect_coralliv_count   == 0), "\n")

# ── Site-level dataset ────────────────────────────────────────
total_model_data <- coralliv_model_data %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_coralliv_biomass, na.rm = TRUE),
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

jpeg("diagnostics_site_coralliv_F1_gaussian_log.jpg",
     width = 25, height = 15, units = "cm", res = 300)
plot(resS_F1, main = "DHARMa — Gaussian on log(y + 0.01)")
dev.off()

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

jpeg("diagnostics_site_coralliv_F2_tweedie.jpg",
     width = 25, height = 15, units = "cm", res = 300)
plot(resS_F2, main = "DHARMa — Tweedie")
dev.off()

plot(resS_F2)
testZeroInflation(resS_F2)
testDispersion(resS_F2)

# ── F3: Zero-inflated Tweedie ─────────────────────────────────
# Adds an explicit Bernoulli component for structural zeros
# (sites where corallivores are chronically absent). Only justified
# if F2 diagnostics show remaining zero inflation.
mS_F3 <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc +
                   log_market_gravity_sc + rugosity_sc,
                 family    = tweedie(link = "log"),
                 ziformula = ~1,
                 data      = total_model_data)

resS_F3 <- simulateResiduals(mS_F3, n = 1000)

jpeg("diagnostics_site_coralliv_F3_tweedie_zi.jpg",
     width = 25, height = 15, units = "cm", res = 300)
plot(resS_F3, main = "DHARMa — Zero-inflated Tweedie")
dev.off()

plot(resS_F3)
testZeroInflation(resS_F3)
testDispersion(resS_F3)

# AICc comparison: F2 vs F3 only (same response, same link)
cat("\n--- Family selection: site-level corallivore ---\n")
print(make_aicc_df(list(
  "Tweedie"    = mS_F2,
  "ZI Tweedie" = mS_F3
)))

# ── Family selection decision ─────────────────────────────────
# Plain Tweedie is preferred.
# ZI Tweedie shows slightly higher AICc (ΔAICc = 1.09),
# indicating no meaningful improvement in model fit.
#
# Following the decision rule:
# adopt ZI only if AICc improves by > 2 and zero-inflation is supported.
# Neither condition is met.
#
# Final decision: retain plain Tweedie for parsimony.

# ── RANDOM EFFECT STRUCTURE SELECTION ────────────────────────
# Anchor: full market gravity model.
# Compare no RE vs country RE.

re_null <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc +
                     log_market_gravity_sc + rugosity_sc,
                   family = tweedie(link = "log"),
                   data = total_model_data)

re_country <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc +
                        log_market_gravity_sc + rugosity_sc +
                        (1 | country),
                      family = tweedie(link = "log"),
                      data = total_model_data)

cat("\n--- RE structure comparison (site-level corallivore) ---\n")
print(make_aicc_df(list(
  "No RE"         = re_null,
  "(1 | country)" = re_country
)))

# ── RE structure decision — site-level corallivore ────────────
# No random effect is preferred (lowest AICc).
# Adding (1 | country) worsens model fit (ΔAICc = 2.62),
# indicating no meaningful country-level structure.
#
# Final decision: retain fixed-effects-only model (no RE).

# ── CANDIDATE MODELS — SITE LEVEL ────────────────────────────

# --- Null ---
s1_m0 <- glmmTMB(mean_biomass ~ 1,
                 family = tweedie(link = "log"),
                 data = total_model_data)

# --- Single predictor ---
s1_m_env      <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc,
                         family = tweedie(link = "log"),
                         data = total_model_data)
s1_m_market   <- glmmTMB(mean_biomass ~ log_market_gravity_sc,
                         family = tweedie(link = "log"),
                         data = total_model_data)
s1_m_settgrav <- glmmTMB(mean_biomass ~ log_settlement_grav_sc,
                         family = tweedie(link = "log"),
                         data = total_model_data)
s1_m_settpop  <- glmmTMB(mean_biomass ~ log_settlement_pop_sc,
                         family = tweedie(link = "log"),
                         data = total_model_data)
s1_m_hab      <- glmmTMB(mean_biomass ~ rugosity_sc,
                         family = tweedie(link = "log"),
                         data = total_model_data)

# --- Environment + human pressure ---
s1_m_env_mkt      <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc,
                             family = tweedie(link = "log"),
                             data = total_model_data)
s1_m_env_settgrav <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc,
                             family = tweedie(link = "log"),
                             data = total_model_data)
s1_m_env_settpop  <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc,
                             family = tweedie(link = "log"),
                             data = total_model_data)

# --- Habitat + human pressure ---
s1_m_hab_market   <- glmmTMB(mean_biomass ~ rugosity_sc + log_market_gravity_sc,
                             family = tweedie(link = "log"),
                             data = total_model_data)
s1_m_hab_settgrav <- glmmTMB(mean_biomass ~ rugosity_sc + log_settlement_grav_sc,
                             family = tweedie(link = "log"),
                             data = total_model_data)
s1_m_hab_settpop  <- glmmTMB(mean_biomass ~ rugosity_sc + log_settlement_pop_sc,
                             family = tweedie(link = "log"),
                             data = total_model_data)

# --- Full (single human pressure metric) ---
s1_m_full_mkt      <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc,
                              family = tweedie(link = "log"),
                              data = total_model_data)
s1_m_full_settgrav <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc,
                              family = tweedie(link = "log"),
                              data = total_model_data)
s1_m_full_settpop  <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc  + rugosity_sc,
                              family = tweedie(link = "log"),
                              data = total_model_data)

# --- Combined gravity metrics ---
s1_m_both_grav     <- glmmTMB(mean_biomass ~ log_market_gravity_sc + log_settlement_grav_sc,
                              family = tweedie(link = "log"),
                              data = total_model_data)
s1_m_hab_both_grav <- glmmTMB(mean_biomass ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc,
                              family = tweedie(link = "log"),
                              data = total_model_data)
s1_m_env_both_grav <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc,
                              family = tweedie(link = "log"),
                              data = total_model_data)
s1_m_full_both_grav <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc,
                               family = tweedie(link = "log"),
                               data = total_model_data)

model_list_s1 <- list(
  "Null"                     = s1_m0,
  "Environment"              = s1_m_env,
  "Market gravity"           = s1_m_market,
  "Settlement gravity"       = s1_m_settgrav,
  "Settlement pop."          = s1_m_settpop,
  "Habitat"                  = s1_m_hab,
  "Env + market gravity"     = s1_m_env_mkt,
  "Env + settlement gravity" = s1_m_env_settgrav,
  "Env + settlement pop."    = s1_m_env_settpop,
  "Habitat + market gravity" = s1_m_hab_market,
  "Habitat + settlement gravity" = s1_m_hab_settgrav,
  "Habitat + settlement pop."    = s1_m_hab_settpop,
  "Full (market gravity)"    = s1_m_full_mkt,
  "Full (settlement gravity)"= s1_m_full_settgrav,
  "Full (settlement pop.)"   = s1_m_full_settpop,
  "Both gravity"             = s1_m_both_grav,
  "Habitat + both gravity"   = s1_m_hab_both_grav,
  "Env + both gravity"       = s1_m_env_both_grav,
  "Full (both gravity)"      = s1_m_full_both_grav
)

cat("\n--- AICc: Site-level corallivore candidate models (no RE) ---\n")
aicc_site_s1 <- make_aicc_df(model_list_s1)
print(aicc_site_s1)

# ── Residual diagnostics — top site-only models ──────────────
# Based on AICc results:
# s1 (no RE):      Null
# Habitat is the closest alternative site-only model.
# Diagnostics run on the top site-only model, plus habitat as a
# simple sensitivity check.

# s1 best: Null
cat("\n--- Diagnostics: s1 Null ---\n")
res_s1_null <- simulateResiduals(s1_m0, n = 1000)

jpeg("diagnostics_site_coralliv_s1_null.jpg",
     width = 25, height = 15, units = "cm", res = 300)
plot(res_s1_null, main = "DHARMa — s1: Null")
dev.off()

plot(res_s1_null)
testZeroInflation(res_s1_null)
testDispersion(res_s1_null)
testOutliers(res_s1_null)

# s1 habitat: closest alternative
cat("\n--- Diagnostics: s1 Habitat ---\n")
res_s1_hab <- simulateResiduals(s1_m_hab, n = 1000)

jpeg("diagnostics_site_coralliv_s1_hab.jpg",
     width = 25, height = 15, units = "cm", res = 300)
plot(res_s1_hab, main = "DHARMa — s1: Habitat")
dev.off()

plot(res_s1_hab)
testZeroInflation(res_s1_hab)
testDispersion(res_s1_hab)
testOutliers(res_s1_hab)

# ── Summaries ─────────────────────────────────────────────────
cat("\n--- Summary: s1 Null ---\n")
summary(s1_m0)

cat("\n--- Summary: s1 Habitat ---\n")
summary(s1_m_hab)

# ── Marginal effect plots ─────────────────────────────────────
# No marginal effects for the null model.
# Habitat is retained only as a visual sensitivity check.

( p_site_rugosity <- plot_effect(s1_m_hab,
                                 total_model_data,
                                 "rugosity_sc",
                                 "Rugosity (scaled)",
                                 y_label = "Corallivore biomass (g)") )

jpeg("site_coralliv_marginal_effects.jpg",
     width = 11, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(p_site_rugosity, ncol = 1)
dev.off()

# ============================================================
#  PART 2 — TRANSECT-LEVEL BIOMASS (SENSITIVITY CHECK)
#
#  Rationale: Retains within-site variation. (1 | site) accounts
#  for non-independence of transects.
#  Confirms site-level findings are not an artefact of collapsing
#  to site means.
#
#  Response:   transect_coralliv_biomass — continuous, zero-inflated
#              (~23% zeros at transect level require Tweedie)
#  Family:     Tweedie (log link) — justified above
#  Random fx:  (1 | site), then (1 | country/site)
# ============================================================

# ── Explore transect-level response ──────────────────────────
summary(coralliv_transects$transect_coralliv_biomass)

zeros <- mean(coralliv_transects$transect_coralliv_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")
# ~0.23 → strong case for Tweedie

( coralliv_raw <- ggplot(coralliv_transects, aes(x = transect_coralliv_biomass)) +
    geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
    labs(x = "corallivore biomass per transect (g)", y = "Frequency",
         title = "Raw corallivore Biomass") +
    theme_bw() )

coralliv_transects <- coralliv_transects %>%
  mutate(
    log_coralliv_biomass  = log(transect_coralliv_biomass + 0.01),
    sqrt_coralliv_biomass = sqrt(transect_coralliv_biomass)
  )

( coralliv_log <- ggplot(coralliv_transects, aes(x = log_coralliv_biomass)) +
    geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
    labs(x = "log(biomass + 0.01)", y = "Frequency",
         title = "Log-transformed corallivore Biomass") +
    theme_bw() )

( coralliv_sqrt <- ggplot(coralliv_transects, aes(x = sqrt_coralliv_biomass)) +
    geom_histogram(bins = 30, fill = "#d7191c", colour = "white") +
    labs(x = "sqrt(biomass)", y = "Frequency",
         title = "Sqrt-transformed corallivore Biomass") +
    theme_bw() )

jpeg("coralliv_biomass_distributions.jpg", width = 33, height = 11,
     units = "cm", res = 300)
gridExtra::grid.arrange(coralliv_raw, coralliv_log, coralliv_sqrt, ncol = 3)
dev.off()

# Add log coralliv biomass 
transect_model_data <- coralliv_transects %>%
  left_join(final_predictors, by = "site") %>%
  mutate(log_coralliv_biomass = log(transect_coralliv_biomass + 0.01)) 

# ── Box-Cox on non-zero values ────────────────────────────────
coralliv_nonzero <- coralliv_transects %>% filter(transect_coralliv_biomass > 0)

MASS::boxcox(
  lm(transect_coralliv_biomass ~ 1, data = coralliv_nonzero),
  lambda = seq(-2, 2, 0.1)
)
# lambda ~ 0 → log transformation appropriate for positive values

# ── Variation by site ─────────────────────────────────────────
ggplot(transect_model_data,
       aes(x = reorder(site, transect_coralliv_biomass, median),
           y = transect_coralliv_biomass)) +
  geom_boxplot(fill = "#2c7bb6", alpha = 0.6,
               outlier.colour = "black", outlier.size = 1) +
  coord_flip() +
  labs(x = NULL, y = "corallivore biomass (g)",
       title = "corallivore biomass distribution by site") +
  theme_bw(base_size = 9)

# ── Zeros by site ─────────────────────────────────────────────
transect_model_data %>%
  group_by(site) %>%
  summarise(
    n_transects  = n(),
    prop_zeros   = mean(transect_coralliv_biomass == 0),
    mean_biomass = mean(transect_coralliv_biomass),
    .groups = "drop"
  ) %>%
  arrange(desc(prop_zeros)) %>%
  print(n = Inf)

# ── Family selection ──────────────────────────────────────────
# AICc not comparable between F1 and F2/F3 (different response).
# Select on DHARMa diagnostics; use AICc only to compare F2 vs F3.

# F1: Gaussian on log(y + 0.01)
mF1_gaussian <- glmmTMB(
  log_coralliv_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = gaussian(),
  data   = transect_model_data
)

resF1 <- simulateResiduals(mF1_gaussian, n = 1000)

jpeg("dharma_coralliv_F1_gaussian.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF1, main = "DHARMa — Gaussian on log(y + 0.01)"); dev.off()

plot(resF1)
testZeroInflation(resF1)
testDispersion(resF1)

# F2: Plain Tweedie
mF2_tweedie <- glmmTMB(
  transect_coralliv_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = transect_model_data
)

resF2 <- simulateResiduals(mF2_tweedie, n = 1000)

jpeg("dharma_coralliv_F2_tweedie.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF2, main = "DHARMa — Tweedie"); dev.off()

plot(resF2)
testZeroInflation(resF2)
testDispersion(resF2)

# F3: Zero-inflated Tweedie
mF3_tweedie_zi <- glmmTMB(
  transect_coralliv_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = transect_model_data
)

resF3 <- simulateResiduals(mF3_tweedie_zi, n = 1000)

jpeg("dharma_coralliv_F3_tweedie_zi.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF3, main = "DHARMa — Zero-inflated Tweedie"); dev.off()

plot(resF3)
testZeroInflation(resF3)
testDispersion(resF3)

cat("\n--- Family selection: transect-level corallivore biomass ---\n")
print(make_aicc_df(list(
  "Tweedie"    = mF2_tweedie,
  "ZI Tweedie" = mF3_tweedie_zi
)))

# ── Family selection decision ─────────────────────────────────
# Plain Tweedie is preferred.
# ZI Tweedie has higher AICc (ΔAICc = 2.13), indicating worse fit.
#
# Following the decision rule:
# adopt ZI only if AICc improves by > 2 and zero-inflation is supported.
# Neither condition is met.
#
# Final decision: retain plain Tweedie.

# ── Random effect structure selection ────────────────────────
# Anchor: full market gravity model.

re_t_null   <- glmmTMB(transect_coralliv_biomass ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc,
                       family = tweedie(link = "log"), data = transect_model_data)

re_t_site   <- glmmTMB(transect_coralliv_biomass ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | site),
                       family = tweedie(link = "log"), data = transect_model_data)

re_t_nested <- glmmTMB(transect_coralliv_biomass ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | country/site),
                       family = tweedie(link = "log"), data = transect_model_data)

cat("\n--- RE structure comparison (transect-level corallivore) ---\n")
print(make_aicc_df(list(
  "No RE"              = re_t_null,
  "(1 | site)"         = re_t_site,
  "(1 | country/site)" = re_t_nested
)))

# ── RE structure decision — transect-level corallivore biomass ──
# (1 | site) has the lowest AICc and is preferred.
# (1 | country/site) is within ΔAICc = 1.71, indicating near-equivalence,
# but does not justify added complexity.
# No RE is strongly unsupported (ΔAICc = 48.12), confirming substantial
# within-site dependence among transects.
#
# Final decision: retain (1 | site).

# ── Candidate models ──────────────────────────────────────────
# Family: Tweedie 
# RE: (1 | site)

coralliv_family <- tweedie(link = "log")

# --- Null ---
m0                 <- glmmTMB(transect_coralliv_biomass ~ 1                                                                                        + (1 | site), family = coralliv_family, data = transect_model_data)

# --- Single predictor ---
m_env              <- glmmTMB(transect_coralliv_biomass ~ sst_sc + log_chla_sc                                                                    + (1 | site), family = coralliv_family, data = transect_model_data)
m_market           <- glmmTMB(transect_coralliv_biomass ~ log_market_gravity_sc                                                                   + (1 | site), family = coralliv_family, data = transect_model_data)
m_settgrav         <- glmmTMB(transect_coralliv_biomass ~ log_settlement_grav_sc                                                                  + (1 | site), family = coralliv_family, data = transect_model_data)
m_settpop          <- glmmTMB(transect_coralliv_biomass ~ log_settlement_pop_sc                                                                   + (1 | site), family = coralliv_family, data = transect_model_data)
m_hab              <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc                                                                             + (1 | site), family = coralliv_family, data = transect_model_data)

# --- Environment + human pressure ---
m_env_market       <- glmmTMB(transect_coralliv_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc                                           + (1 | site), family = coralliv_family, data = transect_model_data)
m_env_settgrav     <- glmmTMB(transect_coralliv_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc                                          + (1 | site), family = coralliv_family, data = transect_model_data)
m_env_settpop      <- glmmTMB(transect_coralliv_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc                                           + (1 | site), family = coralliv_family, data = transect_model_data)

# --- Habitat + human pressure ---
m_hab_market       <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc + log_market_gravity_sc                                                     + (1 | site), family = coralliv_family, data = transect_model_data)
m_hab_settgrav     <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc + log_settlement_grav_sc                                                    + (1 | site), family = coralliv_family, data = transect_model_data)
m_hab_settpop      <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc + log_settlement_pop_sc                                                     + (1 | site), family = coralliv_family, data = transect_model_data)

# --- Full (single human pressure metric) ---
m_full_market      <- glmmTMB(transect_coralliv_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc                            + (1 | site), family = coralliv_family, data = transect_model_data)
m_full_settgrav    <- glmmTMB(transect_coralliv_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc                            + (1 | site), family = coralliv_family, data = transect_model_data)
m_full_settpop     <- glmmTMB(transect_coralliv_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc  + rugosity_sc                            + (1 | site), family = coralliv_family, data = transect_model_data)

# --- Combined gravity metrics ---
m_both_grav        <- glmmTMB(transect_coralliv_biomass ~ log_market_gravity_sc + log_settlement_grav_sc                                         + (1 | site), family = coralliv_family, data = transect_model_data)
m_hab_both_grav    <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc                           + (1 | site), family = coralliv_family, data = transect_model_data)
m_env_both_grav    <- glmmTMB(transect_coralliv_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc                  + (1 | site), family = coralliv_family, data = transect_model_data)
m_full_both_grav   <- glmmTMB(transect_coralliv_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc    + (1 | site), family = coralliv_family, data = transect_model_data)

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

cat("\n--- AICc: Transect-level corallivore biomass ---\n")
print(make_aicc_df(model_list_transect))

# ── Diagnostics on best models ───────────────────────────────
# Best-supported models (ΔAICc < 2):
#   1. Null (weight = 0.17)
#   2. Habitat (weight = 0.15)

# ── Null model ───────────────────────────────────────────────
cat("\n--- Diagnostics: Null model ---\n")
res_t_null <- simulateResiduals(m0, n = 1000)

jpeg("dharma_coralliv_transect_null.jpg",
     width = 25, height = 15, units = "cm", res = 300)
plot(res_t_null, main = "DHARMa — corallivore: Null model")
dev.off()

plot(res_t_null)
testZeroInflation(res_t_null)
testDispersion(res_t_null)
testOutliers(res_t_null)

summary(m0)

# ── Habitat model ────────────────────────────────────────────
cat("\n--- Diagnostics: Habitat model ---\n")
res_t_hab <- simulateResiduals(m_hab, n = 1000)

jpeg("dharma_coralliv_transect_habitat.jpg",
     width = 25, height = 15, units = "cm", res = 300)
plot(res_t_hab, main = "DHARMa — corallivore: Habitat model")
dev.off()

plot(res_t_hab)
testZeroInflation(res_t_hab)
testDispersion(res_t_hab)
testOutliers(res_t_hab)

plotResiduals(res_t_hab, transect_model_data$rugosity_sc,
              xlab = "Rugosity")

summary(m_hab)

# ── Interpretation ───────────────────────────────────────────
# Transect-level corallivore biomass shows no strong support
# for any predictor set. The null model is best-supported,
# with habitat providing only marginal improvement (ΔAICc = 0.26).
# This indicates weak or highly variable relationships at the
# transect scale.

# ── Marginal effect plots ─────────────────────────────────────
# Only meaningful for habitat model (null has no predictors)

( p_t_rugosity <- plot_effect(m_hab,
                              transect_model_data,
                              "rugosity_sc",
                              "Rugosity (scaled)",
                              y_label = "Corallivore biomass (g)") )

jpeg("transect_coralliv_marginal_effects.jpg",
     width = 11, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(p_t_rugosity, ncol = 1)
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
#  Response: Total corallivore count per transect (integer >= 0)
#  Family:   Poisson → NB2 → NB1; selected via AICc
#  Random fx: (1 | site)
# ============================================================

# ── Explore count distribution ────────────────────────────────
cat("Transects (count data):", nrow(transect_model_data), "\n")
cat("Zeros in count data:",    sum(transect_model_data$transect_coralliv_count == 0), "\n")
cat("Proportion zeros:",       round(mean(transect_model_data$transect_coralliv_count == 0), 3), "\n")

summary(transect_model_data$transect_coralliv_count)

ggplot(transect_model_data, aes(x = transect_coralliv_count)) +
  geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
  labs(x = "Total corallivore count per transect", y = "Frequency",
       title = "Raw corallivore count distribution") +
  theme_bw()

# ── Mean-variance relationship ────────────────────────────────
# Points above the Poisson line → overdispersion → Negative Binomial.
transect_model_data %>%
  group_by(site) %>%
  summarise(mean_count = mean(transect_coralliv_count),
            var_count  = var(transect_coralliv_count),
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
  transect_coralliv_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = poisson(link = "log"),
  data   = transect_model_data
)

resC1 <- simulateResiduals(mC1_poisson, n = 1000)

jpeg("dharma_coralliv_C1_poisson.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC1, main = "DHARMa — Poisson"); dev.off()

plot(resC1)
testDispersion(resC1)
testZeroInflation(resC1)
testOutliers(resC1)

# C2: NB2 — quadratic variance (classic NB)
mC2_nb2 <- glmmTMB(
  transect_coralliv_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = nbinom2(link = "log"),
  data   = transect_model_data
)

resC2 <- simulateResiduals(mC2_nb2, n = 1000)

jpeg("dharma_coralliv_C2_nb2.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC2, main = "DHARMa — NB2"); dev.off()

plot(resC2)
testDispersion(resC2)
testZeroInflation(resC2)
testOutliers(resC2)

# C3: NB1 — linear variance
mC3_nb1 <- glmmTMB(
  transect_coralliv_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = nbinom1(link = "log"),
  data   = transect_model_data
)

resC3 <- simulateResiduals(mC3_nb1, n = 1000)

jpeg("dharma_coralliv_C3_nb1.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC3, main = "DHARMa — NB1"); dev.off()

plot(resC3)
testDispersion(resC3)
testZeroInflation(resC3)
testOutliers(resC3)

cat("\n--- Family selection: corallivore count models ---\n")
print(make_aicc_df(list(
  "Poisson" = mC1_poisson,
  "NB2"     = mC2_nb2,
  "NB1"     = mC3_nb1
)))

## NB1 strongly supported (ΔAICc > 20 vs NB2; weight = 1.00).
# All subsequent count models fitted using NB1.

# ── Random effect structure selection ────────────────────────
# Family selected above: NB1.
# Compare no RE, site-level RE, and nested country/site RE.

re_c_null   <- glmmTMB(transect_coralliv_count ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc,
                       family = nbinom1(link = "log"),
                       data = transect_model_data)

re_c_site   <- glmmTMB(transect_coralliv_count ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | site),
                       family = nbinom1(link = "log"),
                       data = transect_model_data)

re_c_nested <- glmmTMB(transect_coralliv_count ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | country/site),
                       family = nbinom1(link = "log"),
                       data = transect_model_data)

print(make_aicc_df(list(
  "No RE"              = re_c_null,
  "(1 | site)"         = re_c_site,
  "(1 | country/site)" = re_c_nested
)))

# ── RE structure decision — corallivore counts ────────────────
# (1 | site) has the lowest AICc and is preferred.
# (1 | country/site) is within ΔAICc = 1.73, indicating near-equivalence,
# but does not justify additional complexity.
# No RE is strongly unsupported (ΔAICc = 55.64), confirming strong
# within-site dependence among transects.
#
# Final decision: retain (1 | site).

# ── Candidate models ──────────────────────────────────────────
# Family: NB1; RE: (1 | site).

count_family <- nbinom1(link = "log")

# --- Null ---
cm0                <- glmmTMB(transect_coralliv_count ~ 1                                                                                        + (1 | site), family = count_family, data = transect_model_data)

# --- Single predictor ---
cm_env             <- glmmTMB(transect_coralliv_count ~ sst_sc + log_chla_sc                                                                    + (1 | site), family = count_family, data = transect_model_data)
cm_market          <- glmmTMB(transect_coralliv_count ~ log_market_gravity_sc                                                                   + (1 | site), family = count_family, data = transect_model_data)
cm_settgrav        <- glmmTMB(transect_coralliv_count ~ log_settlement_grav_sc                                                                  + (1 | site), family = count_family, data = transect_model_data)
cm_settpop         <- glmmTMB(transect_coralliv_count ~ log_settlement_pop_sc                                                                   + (1 | site), family = count_family, data = transect_model_data)
cm_hab             <- glmmTMB(transect_coralliv_count ~ rugosity_sc                                                                             + (1 | site), family = count_family, data = transect_model_data)

# --- Environment + human pressure ---
cm_env_mkt         <- glmmTMB(transect_coralliv_count ~ sst_sc + log_chla_sc + log_market_gravity_sc                                           + (1 | site), family = count_family, data = transect_model_data)
cm_env_settgrav    <- glmmTMB(transect_coralliv_count ~ sst_sc + log_chla_sc + log_settlement_grav_sc                                          + (1 | site), family = count_family, data = transect_model_data)
cm_env_settpop     <- glmmTMB(transect_coralliv_count ~ sst_sc + log_chla_sc + log_settlement_pop_sc                                           + (1 | site), family = count_family, data = transect_model_data)

# --- Habitat + human pressure ---
cm_hab_market      <- glmmTMB(transect_coralliv_count ~ rugosity_sc + log_market_gravity_sc                                                     + (1 | site), family = count_family, data = transect_model_data)
cm_hab_settgrav    <- glmmTMB(transect_coralliv_count ~ rugosity_sc + log_settlement_grav_sc                                                    + (1 | site), family = count_family, data = transect_model_data)
cm_hab_settpop     <- glmmTMB(transect_coralliv_count ~ rugosity_sc + log_settlement_pop_sc                                                     + (1 | site), family = count_family, data = transect_model_data)

# --- Full (single human pressure metric) ---
cm_full_mkt        <- glmmTMB(transect_coralliv_count ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc                            + (1 | site), family = count_family, data = transect_model_data)
cm_full_settgrav   <- glmmTMB(transect_coralliv_count ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc                            + (1 | site), family = count_family, data = transect_model_data)
cm_full_settpop    <- glmmTMB(transect_coralliv_count ~ sst_sc + log_chla_sc + log_settlement_pop_sc  + rugosity_sc                            + (1 | site), family = count_family, data = transect_model_data)

# --- Combined gravity metrics ---
cm_both_grav       <- glmmTMB(transect_coralliv_count ~ log_market_gravity_sc + log_settlement_grav_sc                                         + (1 | site), family = count_family, data = transect_model_data)
cm_hab_both_grav   <- glmmTMB(transect_coralliv_count ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc                           + (1 | site), family = count_family, data = transect_model_data)
cm_env_both_grav   <- glmmTMB(transect_coralliv_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc                  + (1 | site), family = count_family, data = transect_model_data)
cm_full_both_grav  <- glmmTMB(transect_coralliv_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc    + (1 | site), family = count_family, data = transect_model_data)

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

cat("\n--- AICc: corallivore count models ---\n")
print(make_aicc_df(model_list_counts))
# ── Diagnostics on best model ─────────────────────────────────
# Null is best-supported (AICc weight = 0.17).
# Six models fall within ΔAICc < 2 (Null, Habitat, Both gravity,
# Market gravity, Settlement gravity, Habitat + market gravity),
# indicating no meaningful signal — model uncertainty is spread
# across intercept-only and weakly-supported single-predictor models.
# DHARMa diagnostics on the null model show no dispersion issues
# (p = 0.244), no zero inflation (p = 0.406), and no outliers.
# Habitat model diagnostics are near-identical, confirming rugosity
# (IRR = 1.13, p = 0.252) adds no meaningful explanatory power.

res_cm_null <- simulateResiduals(cm0, n = 1000)

jpeg("dharma_coralliv_count_null.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_cm_null, main = "DHARMa — corallivore counts: Null"); dev.off()

plot(res_cm_null)
testDispersion(res_cm_null)
testZeroInflation(res_cm_null)
testOutliers(res_cm_null)

# Optional comparison model for interpretation
res_cm_hab <- simulateResiduals(cm_hab, n = 1000)

jpeg("dharma_coralliv_count_hab.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_cm_hab, main = "DHARMa — corallivore counts: Habitat"); dev.off()

plot(res_cm_hab)
testDispersion(res_cm_hab)
testZeroInflation(res_cm_hab)
testOutliers(res_cm_hab)

# Residual plots against predictors — only for the comparison model
# if you want to inspect possible weak structure.
plotResiduals(res_cm_hab, transect_model_data$rugosity_sc,           xlab = "Rugosity")
plotResiduals(res_cm_hab, transect_model_data$log_market_gravity_sc, xlab = "Market gravity")

summary(cm0)
summary(cm_hab)
exp(fixef(cm_hab)$cond)

# ── Marginal effect plots ─────────────────────────────────────
# Use the comparison model cautiously; effects are weak and do not
# outperform the null model in AICc.
p_count_hab <- plot_effect(cm_hab, transect_model_data, "rugosity_sc",
                           "Rugosity (scaled)", "Expected corallivore count",
                           colour = "#d7191c")

jpeg("coralliv_count_marginal_effects.jpg",
     width = 11, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(p_count_hab, ncol = 1)
dev.off()

# ── Key divergence from biomass models ───────────────────────
# Biomass (Parts 1 & 2) and counts (Part 3) are consistent:
#   the null model is best-supported across all three analyses.
#   No predictor set clearly outperforms the intercept-only model.
#   Corallivore biomass and abundance appear weakly structured by
#   SST, Chl-a, human gravity, and rugosity at this scale.

# ── SYNTHESIS ─────────────────────────────────────────
# Part 1 (site biomass):      Null model best-supported.
#                             No predictor set meaningfully improves
#                             on the intercept-only model. Habitat
#                             (rugosity) is the closest alternative
#                             but does not meet the threshold for
#                             clear support.
#
# Part 2 (transect biomass):  Null model best-supported (weight = 0.17).
#                             Habitat provides only marginal improvement
#                             (ΔAICc = 0.26), confirming the site-level
#                             finding is not an artefact of averaging.
#                             No predictor is robustly supported at the
#                             transect scale.
#
# Part 3 (transect counts):   Null model best-supported (weight = 0.17).
#                             Six models fall within ΔAICc < 2, but none
#                             meaningfully improves on the null. Results
#                             are consistent with Parts 1 & 2 —
#                             corallivore biomass and abundance are both
#                             weakly structured by the tested predictors
#                             at this scale.
#
# Overall: No predictor is robustly supported across any of the three
#   analytical approaches. The null model is best or near-best in all
#   parts. SST, Chl-a, human gravity metrics, and rugosity show no
#   clear or consistent signal in corallivore biomass or abundance.
#   Habitat complexity (rugosity) is the most recurrent near-null
#   alternative but never clears the ΔAICc > 2 threshold.

# ── IRR summary table ─────────────────────────────────────────
# Reported from the habitat model (closest alternative to null;
# ΔAICc = 0.76). Treat as indicative only — rugosity is not
# significant (p = 0.252) and this model does not outperform the null.
cat("\n--- IRR: Habitat model (closest alternative to null) ---\n")
irr_hab <- exp(fixef(cm_hab)$cond)
se_hab  <- summary(cm_hab)$coefficients$cond[, "Std. Error"]

cat("Rugosity:  IRR =", round(irr_hab["rugosity_sc"], 2),
    " (95% CI:", round(exp(log(irr_hab["rugosity_sc"]) - 1.96 * se_hab["rugosity_sc"]), 2),
    "-",         round(exp(log(irr_hab["rugosity_sc"]) + 1.96 * se_hab["rugosity_sc"]), 2),
    ", p = 0.252)\n")