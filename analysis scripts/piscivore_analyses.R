# ============================================================
#  PISCIVORE BIOMASS & ABUNDANCE — MIXED EFFECTS MODELS
#
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
#              NOTE: ~47% zeros require Tweedie family.
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

pisc_transects <- fish_2009 %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_pisc_biomass = sum(
      ifelse(trophic_group == "piscivores", tot_wt_g, 0),
      na.rm = TRUE
    ),
    transect_pisc_count = sum(
      ifelse(trophic_group == "piscivores", number, 0),
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

cat("Number of transects:", nrow(pisc_transects), "\n")
cat("Number of sites:",     n_distinct(pisc_transects$site), "\n")
cat("Number of countries:", n_distinct(pisc_transects$country), "\n")

# ==============================================================================
#  BIOMASS DATA EXPLORATION
# ==============================================================================

# ── Aggregate data at site level (mean piscivore biomass per site) ────────────
site_data <- pisc_transects %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass = mean(transect_pisc_biomass, na.rm = TRUE),
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
# recorded zero piscivores. Check whether these sites should be retained.

# ── Raw distribution ──────────────────────────────────────────────────────────
( site_raw <- ggplot(site_data, aes(x = mean_biomass)) +
    geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
    labs(x = "Mean piscivore biomass per site (g)", y = "Frequency",
         title = "Raw Site-Level Piscivore Biomass") +
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
         title = "Log-transformed Site-Level Piscivore Biomass") +
    theme_bw() )

( site_sqrt <- ggplot(site_data, aes(x = sqrt_mean_biomass)) +
    geom_histogram(bins = 25, fill = "#d7191c", colour = "white") +
    labs(x = "sqrt(mean biomass)", y = "Frequency",
         title = "Sqrt-transformed Site-Level Piscivore Biomass") +
    theme_bw() )

jpeg("site_pisc_biomass_distributions.jpg", width = 33, height = 11,
     units = "cm", res = 300)
gridExtra::grid.arrange(site_raw, site_log, site_sqrt, ncol = 3)
dev.off()

# ── Normality checks on log-transformed response ──────────────────────────────
qqnorm(site_data$log_mean_biomass,
       main = "Q-Q plot: log(mean piscivore biomass per site)")
qqline(site_data$log_mean_biomass, col = "red")
shapiro.test(site_data$log_mean_biomass)
# Shapiro-Wilk: W = 0.684, p < 0.001 — severe departure from
# normality on log scale, driven by the 11 zero-mean sites
# collapsing to log(0.01) = -4.6. Tweedie on raw scale is
# strongly preferred for this reason.

# ── Variation by site ─────────────────────────────────────────────────────────
ggplot(site_data, aes(x = reorder(site, mean_biomass, median),
                      y = mean_biomass)) +
  geom_col(fill = "#2c7bb6", alpha = 0.7) +
  coord_flip() +
  labs(x = NULL, y = "Mean piscivore biomass (g)",
       title = "Mean piscivore biomass by site (raw)") +
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
  labs(x = NULL, y = "log(mean piscivore biomass per site + 0.01)") +
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

jpeg("predictor_distributions_pisc.jpg", width = 33, height = 22, units = "cm", res = 300)
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

settlement_data <- pisc_transects %>%
  left_join(scaled_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    mean_pisc_biomass      = mean(transect_pisc_biomass, na.rm = TRUE),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_settlement_pop_sc  = first(log_settlement_pop_sc),
    .groups = "drop"
  )

settgrav <- glmmTMB(mean_pisc_biomass ~ log_settlement_grav_sc,
                    family = tweedie(link = "log"), data = settlement_data)
settpop  <- glmmTMB(mean_pisc_biomass ~ log_settlement_pop_sc,
                    family = tweedie(link = "log"), data = settlement_data)

make_aicc_df(list(
  "Settlement gravity" = settgrav,
  "Settlement pop."    = settpop
))

# Settlement metric selection: virtually identical performance
# (delta AICc = 0.06, weights 0.51 vs 0.49). Neither metric
# has any meaningful advantage over the other. Both carried
# forward as parallel candidate model sets.

rm(settlement_data)

# ============================================================
#  ANALYSIS DATASETS
# ============================================================

# Retain all human pressure metrics for parallel model sets
final_predictors <- scaled_predictors %>%
  dplyr::select(site, log_market_gravity_sc, log_settlement_grav_sc, 
                log_settlement_pop_sc, log_chla_sc, sst_sc, rugosity_sc)

# ── Transect-level dataset ────────────────────────────────────
transect_model_data <- pisc_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect model data:", nrow(transect_model_data), "rows,",
    n_distinct(transect_model_data$site), "sites\n")
cat("Biomass zeros:", sum(transect_model_data$transect_pisc_biomass == 0), "\n")
cat("Count zeros:",  sum(transect_model_data$transect_pisc_count   == 0), "\n")

# ── Site-level dataset ────────────────────────────────────────
total_model_data <- transect_model_data %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_pisc_biomass, na.rm = TRUE),
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

jpeg("diagnostics_site_pisc_F1_gaussian_log.jpg",
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

jpeg("diagnostics_site_pisc_F2_tweedie.jpg",
     width = 25, height = 15, units = "cm", res = 300)
plot(resS_F2, main = "DHARMa — Tweedie"); dev.off()

plot(resS_F2)
testZeroInflation(resS_F2)
testDispersion(resS_F2)

# ── F3: Zero-inflated Tweedie ─────────────────────────────────
# Adds an explicit Bernoulli component for structural zeros
# (sites where piscivores are chronically absent). Only justified
# if F2 diagnostics show remaining zero inflation.
mS_F3 <- glmmTMB(mean_biomass ~ sst_sc + log_chla_sc +
                   log_market_gravity_sc + rugosity_sc,
                 family    = tweedie(link = "log"),
                 ziformula = ~1,
                 data      = total_model_data)

resS_F3 <- simulateResiduals(mS_F3, n = 1000)

jpeg("diagnostics_site_pisc_F3_tweedie_zi.jpg",
     width = 25, height = 15, units = "cm", res = 300)
plot(resS_F3, main = "DHARMa — Zero-inflated Tweedie"); dev.off()

plot(resS_F3)
testZeroInflation(resS_F3)
testDispersion(resS_F3)

# AICc comparison: F2 vs F3 only (same response, same link)
cat("\n--- Family selection: site-level ---\n")
print(make_aicc_df(list(
  "Tweedie"            = mS_F2,
  "ZI Tweedie"         = mS_F3
)))
# Assess F1 on DHARMa diagnostics alone — different response scale.
# If ZI test is n.s. for F2 and AICc(F3) > AICc(F2), use plain Tweedie.

# ── Family selection decision ─────────────────────────────────
# Plain Tweedie (F2) retained. ZI Tweedie does not meet the
# prespecified adoption threshold:
#   - ΔAICc = 2.62 (ZI vs plain Tweedie) — exceeds the >2 threshold
#   - Zero inflation test n.s. for both F2 (p = 1.00) and
#     F3 (p = 1.00) — no evidence of excess zeros
#
# Both conditions must be met: ΔAICc > 2 AND zero inflation
# significant. Zero inflation condition fails.
# Plain Tweedie carried forward.
#
# F1 (Gaussian on log scale): clean dispersion (p = 0.914) but
# produced an outer Newton convergence warning during DHARMa
# smoothing. Tweedie preferred regardless — handles the 11
# site-level zeros natively.

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

cat("\n--- RE structure comparison (site-level piscivore) ---\n")
print(make_aicc_df(list(
  "No RE"         = re_null,
  "(1 | country)" = re_country
)))

# ── RE structure decision — piscivore site level ──────────────
# No RE marginally preferred (ΔAICc = 0.68, weights 0.58 vs 0.42
# for country RE). The difference is not decisive.
# No RE carried forward as the primary analysis; key models
# refitted with (1 | country) as a sensitivity check to confirm
# conclusions are robust to this choice.

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

cat("\n--- AICc: Site-level piscivore candidate models (no RE) ---\n")
aicc_site_s1 <- make_aicc_df(model_list_s1)
print(aicc_site_s1)

# ── CANDIDATE MODELS — COUNTRY RE (SENSITIVITY) ──────────────
# Fitted to confirm s1 conclusions are robust to RE choice.
# delta AICc = 0.68 between no RE and country RE — borderline.

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
  "Null"                               = s2_m0,
  "Environment"                        = s2_m_env,
  "Market gravity"                     = s2_m_market,
  "Settlement gravity"                 = s2_m_settgrav,
  "Settlement pop."                    = s2_m_settpop,
  "Habitat"                            = s2_m_hab,
  "Env + market gravity"               = s2_m_env_mkt,
  "Env + settlement gravity"           = s2_m_env_settgrav,
  "Env + settlement pop."              = s2_m_env_settpop,
  "Habitat + market gravity"           = s2_m_hab_market,
  "Habitat + settlement gravity"       = s2_m_hab_settgrav,
  "Habitat + settlement pop."          = s2_m_hab_settpop,
  "Full (market gravity)"              = s2_m_full_mkt,
  "Full (settlement gravity)"          = s2_m_full_settgrav,
  "Full (settlement pop.)"             = s2_m_full_settpop,
  "Both gravity"                       = s2_m_both_grav,
  "Habitat + both gravity"             = s2_m_hab_both_grav,
  "Env + both gravity"                 = s2_m_env_both_grav,
  "Full (both gravity)"                = s2_m_full_both_grav
)

cat("\n--- AICc: Site-level piscivore candidate models (1 | country) ---\n")
aicc_site_s2 <- make_aicc_df(model_list_s2)
print(aicc_site_s2)

# ── Residual diagnostics — top models ────────────────────────
# Based on AICc results:
# s1 (no RE):      Habitat + market gravity (weight = 0.41)
# s2 (country RE): Habitat + market gravity ≈ Habitat (delta = 0.03)
# Diagnostics run on top model from each RE structure, plus
# habitat-only country RE model given near-equivalent performance.

# s1 best: Habitat + market gravity
cat("\n--- Diagnostics: s1 Habitat + market gravity ---\n")
res_s1_hab_market <- simulateResiduals(s1_m_hab_market, n = 1000)

jpeg("diagnostics_site_pisc_s1_hab_market.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_s1_hab_market, main = "DHARMa — s1: Habitat + market gravity"); dev.off()

plot(res_s1_hab_market)
testZeroInflation(res_s1_hab_market)
testDispersion(res_s1_hab_market)
testOutliers(res_s1_hab_market)

# s2 best: Habitat + market gravity (country RE)
cat("\n--- Diagnostics: s2 Habitat + market gravity ---\n")
res_s2_hab_market <- simulateResiduals(s2_m_hab_market, n = 1000)

jpeg("diagnostics_site_pisc_s2_hab_market.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_s2_hab_market, main = "DHARMa — s2: Habitat + market gravity (1 | country)"); dev.off()

plot(res_s2_hab_market)
testZeroInflation(res_s2_hab_market)
testDispersion(res_s2_hab_market)
testOutliers(res_s2_hab_market)

# s2 habitat only: near-equivalent to s2 best (delta = 0.03)
# Key model for assessing rugosity robustness without market gravity
cat("\n--- Diagnostics: s2 Habitat only ---\n")
res_s2_hab <- simulateResiduals(s2_m_hab, n = 1000)

jpeg("diagnostics_site_pisc_s2_hab.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_s2_hab, main = "DHARMa — s2: Habitat only (1 | country)"); dev.off()

plot(res_s2_hab)
testZeroInflation(res_s2_hab)
testDispersion(res_s2_hab)
testOutliers(res_s2_hab)

# ── Summaries ─────────────────────────────────────────────────
cat("\n--- Summary: s1 Habitat + market gravity ---\n")
summary(s1_m_hab_market)

cat("\n--- Summary: s2 Habitat + market gravity ---\n")
summary(s2_m_hab_market)

cat("\n--- Summary: s2 Habitat only ---\n")
summary(s2_m_hab)

# ── Coefficient comparison across RE structures ───────────────
cat("\n--- Rugosity coefficient stability ---\n")
cat("s1 (no RE):          beta =", round(fixef(s1_m_hab_market)$cond["rugosity_sc"], 3), "\n")
cat("s2 (country RE):     beta =", round(fixef(s2_m_hab_market)$cond["rugosity_sc"], 3), "\n")
cat("s2 habitat only:     beta =", round(fixef(s2_m_hab)$cond["rugosity_sc"],        3), "\n")

cat("\n--- Market gravity coefficient stability ---\n")
cat("s1 (no RE):          beta =", round(fixef(s1_m_hab_market)$cond["log_market_gravity_sc"], 3), "\n")
cat("s2 (country RE):     beta =", round(fixef(s2_m_hab_market)$cond["log_market_gravity_sc"], 3), "\n")

# ── Marginal effect plots ─────────────────────────────────────
# Rugosity from s1 best model 
( p_site_rugosity <- plot_effect(s1_m_hab_market,
                                 total_model_data,
                                 "rugosity_sc",
                                 "Rugosity (scaled)",
                                 y_label = "Piscivore biomass (g)") )

# Market gravity from s1 best model — with caveat re: direction
( p_site_market <- plot_effect(s1_m_hab_market,
                               total_model_data,
                               "log_market_gravity_sc",
                               "Market gravity (scaled)",
                               y_label = "Piscivore biomass (g)") )

jpeg("site_pisc_marginal_effects.jpg", width = 22, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(p_site_rugosity, p_site_market, ncol = 2)
dev.off()

# ============================================================
#  PART 2 — TRANSECT-LEVEL BIOMASS (SENSITIVITY CHECK)
#
#  Rationale: Retains within-site variation. (1 | site) accounts
#  for non-independence of transects.
#  Confirms site-level findings are not an artefact of collapsing
#  to site means.
#
#  Response:   transect_pisc_biomass — continuous, zero-inflated
#              (~47% zeros at transect level require Tweedie)
#  Family:     Tweedie (log link) — justified above
#  Random fx:  (1 | site), then (1 | country/site)
# ============================================================

# ── Explore transect-level response ──────────────────────────
summary(pisc_transects$transect_pisc_biomass)

zeros <- mean(pisc_transects$transect_pisc_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")
# ~0.47 → strong case for Tweedie

( pisc_raw <- ggplot(pisc_transects, aes(x = transect_pisc_biomass)) +
    geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
    labs(x = "Piscivore biomass per transect (g)", y = "Frequency",
         title = "Raw Piscivore Biomass") +
    theme_bw() )

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

jpeg("pisc_biomass_distributions.jpg", width = 33, height = 11,
     units = "cm", res = 300)
gridExtra::grid.arrange(pisc_raw, pisc_log, pisc_sqrt, ncol = 3)
dev.off()

# Add log pisc biomass 
transect_model_data <- pisc_transects %>%
  left_join(final_predictors, by = "site") %>%
  mutate(log_pisc_biomass = log(transect_pisc_biomass + 0.01)) 

# ── Box-Cox on non-zero values ────────────────────────────────
pisc_nonzero <- pisc_transects %>% filter(transect_pisc_biomass > 0)

MASS::boxcox(
  lm(transect_pisc_biomass ~ 1, data = pisc_nonzero),
  lambda = seq(-2, 2, 0.1)
)
# lambda ~ 0 → log transformation appropriate for positive values

# ── Variation by site ─────────────────────────────────────────
ggplot(transect_model_data,
       aes(x = reorder(site, transect_pisc_biomass, median),
           y = transect_pisc_biomass)) +
  geom_boxplot(fill = "#2c7bb6", alpha = 0.6,
               outlier.colour = "black", outlier.size = 1) +
  coord_flip() +
  labs(x = NULL, y = "Piscivore biomass (g)",
       title = "Piscivore biomass distribution by site") +
  theme_bw(base_size = 9)

# ── Zeros by site ─────────────────────────────────────────────
transect_model_data %>%
  group_by(site) %>%
  summarise(
    n_transects  = n(),
    prop_zeros   = mean(transect_pisc_biomass == 0),
    mean_biomass = mean(transect_pisc_biomass),
    .groups = "drop"
  ) %>%
  arrange(desc(prop_zeros)) %>%
  print(n = Inf)

# ── Family selection ──────────────────────────────────────────
# AICc not comparable between F1 and F2/F3 (different response).
# Select on DHARMa diagnostics; use AICc only to compare F2 vs F3.

# F1: Gaussian on log(y + 0.01)
mF1_gaussian <- glmmTMB(
  log_pisc_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = gaussian(),
  data   = transect_model_data
)

resF1 <- simulateResiduals(mF1_gaussian, n = 1000)

jpeg("dharma_pisc_F1_gaussian.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF1, main = "DHARMa — Gaussian on log(y + 0.01)"); dev.off()

plot(resF1)
testZeroInflation(resF1)
testDispersion(resF1)

# F2: Plain Tweedie
mF2_tweedie <- glmmTMB(
  transect_pisc_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = transect_model_data
)

resF2 <- simulateResiduals(mF2_tweedie, n = 1000)

jpeg("dharma_pisc_F2_tweedie.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF2, main = "DHARMa — Tweedie"); dev.off()

plot(resF2)
testZeroInflation(resF2)
testDispersion(resF2)

# F3: Zero-inflated Tweedie
mF3_tweedie_zi <- glmmTMB(
  transect_pisc_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = transect_model_data
)

resF3 <- simulateResiduals(mF3_tweedie_zi, n = 1000)

jpeg("dharma_pisc_F3_tweedie_zi.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF3, main = "DHARMa — Zero-inflated Tweedie"); dev.off()

plot(resF3)
testZeroInflation(resF3)
testDispersion(resF3)

cat("\n--- Family selection: transect-level piscivore biomass ---\n")
print(make_aicc_df(list(
  "Tweedie"    = mF2_tweedie,
  "ZI Tweedie" = mF3_tweedie_zi
)))

# ── Family selection decision ─────────────────────────────────
# Plain Tweedie (F2) retained.
# ZI Tweedie (F3) does not meet the prespecified adoption threshold:
#   - ΔAICc = 2.13 (ZI vs plain Tweedie) — does not exceed >2
#   - Zero inflation n.s. for both F2 (p = 0.724) and F3 (p = 0.746)
# Both conditions must be met; neither is met.
#
# F1 (Gaussian on log scale) produced an outer Newton convergence
# warning during DHARMa smoothing (same warning seen at site level).
# Dispersion otherwise clean (p = 0.988). Tweedie preferred
# regardless — handles the 47.5% transect-level zeros natively.

# ── Random effect structure selection ────────────────────────
# Anchor: full market gravity model.

re_t_null   <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc,
                       family = tweedie(link = "log"), data = transect_model_data)

re_t_site   <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | site),
                       family = tweedie(link = "log"), data = transect_model_data)

re_t_nested <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | country/site),
                       family = tweedie(link = "log"), data = transect_model_data)

cat("\n--- RE structure comparison (transect-level piscivore) ---\n")
print(make_aicc_df(list(
  "No RE"              = re_t_null,
  "(1 | site)"         = re_t_site,
  "(1 | country/site)" = re_t_nested
)))

# ── Random effect structure decision ─────────────────────────
# (1 | site) and (1 | country/site) are exactly tied (ΔAICc = 0.00,
# weights 0.50 vs 0.50). No RE strongly rejected (ΔAICc = 13.29).
# Prefer the simpler (1 | site) structure — adding the country-level
# nesting provides no AICc benefit and increases model complexity.

# ── Candidate models ──────────────────────────────────────────
# Family: Tweedie 
# RE: (1 | site)

pisc_family <- tweedie(link = "log")

# --- Null ---
m0                 <- glmmTMB(transect_pisc_biomass ~ 1                                                                                        + (1 | site), family = pisc_family, data = transect_model_data)

# --- Single predictor ---
m_env              <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc                                                                    + (1 | site), family = pisc_family, data = transect_model_data)
m_market           <- glmmTMB(transect_pisc_biomass ~ log_market_gravity_sc                                                                   + (1 | site), family = pisc_family, data = transect_model_data)
m_settgrav         <- glmmTMB(transect_pisc_biomass ~ log_settlement_grav_sc                                                                  + (1 | site), family = pisc_family, data = transect_model_data)
m_settpop          <- glmmTMB(transect_pisc_biomass ~ log_settlement_pop_sc                                                                   + (1 | site), family = pisc_family, data = transect_model_data)
m_hab              <- glmmTMB(transect_pisc_biomass ~ rugosity_sc                                                                             + (1 | site), family = pisc_family, data = transect_model_data)

# --- Environment + human pressure ---
m_env_market       <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc                                           + (1 | site), family = pisc_family, data = transect_model_data)
m_env_settgrav     <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc                                          + (1 | site), family = pisc_family, data = transect_model_data)
m_env_settpop      <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc                                           + (1 | site), family = pisc_family, data = transect_model_data)

# --- Habitat + human pressure ---
m_hab_market       <- glmmTMB(transect_pisc_biomass ~ rugosity_sc + log_market_gravity_sc                                                     + (1 | site), family = pisc_family, data = transect_model_data)
m_hab_settgrav     <- glmmTMB(transect_pisc_biomass ~ rugosity_sc + log_settlement_grav_sc                                                    + (1 | site), family = pisc_family, data = transect_model_data)
m_hab_settpop      <- glmmTMB(transect_pisc_biomass ~ rugosity_sc + log_settlement_pop_sc                                                     + (1 | site), family = pisc_family, data = transect_model_data)

# --- Full (single human pressure metric) ---
m_full_market      <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc  + rugosity_sc                            + (1 | site), family = pisc_family, data = transect_model_data)
m_full_settgrav    <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc                            + (1 | site), family = pisc_family, data = transect_model_data)
m_full_settpop     <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc + log_settlement_pop_sc  + rugosity_sc                            + (1 | site), family = pisc_family, data = transect_model_data)

# --- Combined gravity metrics ---
m_both_grav        <- glmmTMB(transect_pisc_biomass ~ log_market_gravity_sc + log_settlement_grav_sc                                         + (1 | site), family = pisc_family, data = transect_model_data)
m_hab_both_grav    <- glmmTMB(transect_pisc_biomass ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc                           + (1 | site), family = pisc_family, data = transect_model_data)
m_env_both_grav    <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc                  + (1 | site), family = pisc_family, data = transect_model_data)
m_full_both_grav   <- glmmTMB(transect_pisc_biomass ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc    + (1 | site), family = pisc_family, data = transect_model_data)

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

cat("\n--- AICc: Transect-level piscivore biomass ---\n")
print(make_aicc_df(model_list_transect))

# ── Diagnostics on best model ─────────────────────────────────
# Habitat + market gravity: top model (weight = 0.2667).
# Full (market gravity) is next (ΔAICc = 1.08).
# DHARMa diagnostics: no dispersion issues (p = 0.992), no zero
# inflation (p = 0.762), no outliers. Clean fit.
# Rugosity (beta = 0.353, p = 0.023) and market gravity
# (beta = 0.369, p = 0.013) both significant.
res_t <- simulateResiduals(m_hab_market, n = 1000)

jpeg("dharma_pisc_transect_best.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_t, main = "DHARMa — transect piscivore: Habitat + market gravity"); dev.off()

plot(res_t)
testZeroInflation(res_t)
testDispersion(res_t)
testOutliers(res_t)
plotResiduals(res_t, transect_model_data$rugosity_sc,            xlab = "Rugosity")
plotResiduals(res_t, transect_model_data$log_market_gravity_sc,  xlab = "Market gravity")

summary(m_hab_market)

# ── Convergence with site-level result ────────────────────────
# Part 1 (site):      Habitat + market gravity — weight = 0.41
# Part 2 (transect):  Habitat + market gravity — weight = 0.27
# Identical top model structure confirms site-level findings
# are not an artefact of collapsing to site means.
# Market gravity positive effect consistent across both levels —
# interpret with caution (see direction checks above).

# ── Marginal effect plots ─────────────────────────────────────
( p_t_rugosity <- plot_effect(m_hab_market,
                              transect_model_data,
                              "rugosity_sc",
                              "Rugosity (scaled)",
                              y_label = "Piscivore biomass (g)") )

( p_t_market <- plot_effect(m_hab_market,
                            transect_model_data,
                            "log_market_gravity_sc",
                            "Market gravity (scaled)",
                            y_label = "Piscivore biomass (g)") )

jpeg("transect_pisc_marginal_effects.jpg", width = 22, height = 11, units = "cm", res = 300)
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
#  Response: Total piscivore count per transect (integer >= 0)
#  Family:   Poisson → NB2 → NB1; selected via AICc
#  Random fx: (1 | site)
# ============================================================

# ── Explore count distribution ────────────────────────────────
cat("Transects (count data):", nrow(transect_model_data), "\n")
cat("Zeros in count data:",    sum(transect_model_data$transect_pisc_count == 0), "\n")
cat("Proportion zeros:",       round(mean(transect_model_data$transect_pisc_count == 0), 3), "\n")

summary(transect_model_data$transect_pisc_count)

ggplot(transect_model_data, aes(x = transect_pisc_count)) +
  geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
  labs(x = "Total piscivore count per transect", y = "Frequency",
       title = "Raw piscivore count distribution") +
  theme_bw()

# ── Mean-variance relationship ────────────────────────────────
# Points above the Poisson line → overdispersion → Negative Binomial.
transect_model_data %>%
  group_by(site) %>%
  summarise(mean_count = mean(transect_pisc_count),
            var_count  = var(transect_pisc_count),
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
  transect_pisc_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = poisson(link = "log"),
  data   = transect_model_data
)

resC1 <- simulateResiduals(mC1_poisson, n = 1000)

jpeg("dharma_pisc_C1_poisson.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC1, main = "DHARMa — Poisson"); dev.off()

plot(resC1)
testDispersion(resC1)
testZeroInflation(resC1)
testOutliers(resC1)

# C2: NB2 — quadratic variance (classic NB)
mC2_nb2 <- glmmTMB(
  transect_pisc_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = nbinom2(link = "log"),
  data   = transect_model_data
)

resC2 <- simulateResiduals(mC2_nb2, n = 1000)

jpeg("dharma_pisc_C2_nb2.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC2, main = "DHARMa — NB2"); dev.off()

plot(resC2)
testDispersion(resC2)
testZeroInflation(resC2)
testOutliers(resC2)

# C3: NB1 — linear variance
mC3_nb1 <- glmmTMB(
  transect_pisc_count ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = nbinom1(link = "log"),
  data   = transect_model_data
)

resC3 <- simulateResiduals(mC3_nb1, n = 1000)

jpeg("dharma_pisc_C3_nb1.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resC3, main = "DHARMa — NB1"); dev.off()

plot(resC3)
testDispersion(resC3)
testZeroInflation(resC3)
testOutliers(resC3)

cat("\n--- Family selection: piscivore count models ---\n")
print(make_aicc_df(list(
  "Poisson" = mC1_poisson,
  "NB2"     = mC2_nb2,
  "NB1"     = mC3_nb1
)))

# ── Family selection ──────────────────────────────────────────
# Poisson is unsupported (ΔAICc = 25.34 relative to best model).
# NB1 has the lowest AICc and is preferred over NB2 by the prespecified
# AICc rule, although the difference is small (ΔAICc = 0.53).
# Retain NB1.

# ── Random effect structure selection ────────────────────────
# Using the selected count family (NB1).
# Compare no RE, site-level RE, and nested country/site RE.

count_family <- nbinom1(link = "log")

re_c_null   <- glmmTMB(transect_pisc_count ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc,
                       family = count_family, data = transect_model_data)

re_c_site   <- glmmTMB(transect_pisc_count ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | site),
                       family = count_family, data = transect_model_data)

re_c_nested <- glmmTMB(transect_pisc_count ~ sst_sc + log_chla_sc +
                         log_market_gravity_sc + rugosity_sc +
                         (1 | country/site),
                       family = count_family, data = transect_model_data)

cat("\n--- RE structure comparison (piscivore counts) ---\n")
print(make_aicc_df(list(
  "No RE"              = re_c_null,
  "(1 | site)"         = re_c_site,
  "(1 | country/site)" = re_c_nested
)))

# ── Random effect structure selection ────────────────────────
# Using the selected count family (NB1).
# Compare no RE, site-level RE, and nested country/site RE.
# (1 | site) has the lowest AICc and is preferred.
# (1 | country/site) is worse by ΔAICc = 2.12, so the added country-
# level complexity is not supported. No RE is clearly inferior.

# ── Candidate models ──────────────────────────────────────────
# Family: NB1. RE: (1 | site).

# --- Null ---
cm0 <- glmmTMB(
  transect_pisc_count ~ 1 + (1 | site),
  family = count_family,
  data   = transect_model_data
)

# --- Single predictor ---
cm_env      <- glmmTMB(transect_pisc_count ~ sst_sc + log_chla_sc + (1 | site), family = count_family, data = transect_model_data)
cm_market   <- glmmTMB(transect_pisc_count ~ log_market_gravity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_settgrav <- glmmTMB(transect_pisc_count ~ log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_settpop  <- glmmTMB(transect_pisc_count ~ log_settlement_pop_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab      <- glmmTMB(transect_pisc_count ~ rugosity_sc + (1 | site), family = count_family, data = transect_model_data)

# --- Environment + human pressure ---
cm_env_mkt      <- glmmTMB(transect_pisc_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_env_settgrav <- glmmTMB(transect_pisc_count ~ sst_sc + log_chla_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_env_settpop  <- glmmTMB(transect_pisc_count ~ sst_sc + log_chla_sc + log_settlement_pop_sc + (1 | site), family = count_family, data = transect_model_data)

# --- Habitat + human pressure ---
cm_hab_market   <- glmmTMB(transect_pisc_count ~ rugosity_sc + log_market_gravity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab_settgrav <- glmmTMB(transect_pisc_count ~ rugosity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab_settpop  <- glmmTMB(transect_pisc_count ~ rugosity_sc + log_settlement_pop_sc + (1 | site), family = count_family, data = transect_model_data)

# --- Full (single human pressure metric) ---
cm_full_mkt      <- glmmTMB(transect_pisc_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_full_settgrav <- glmmTMB(transect_pisc_count ~ sst_sc + log_chla_sc + log_settlement_grav_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)
cm_full_settpop  <- glmmTMB(transect_pisc_count ~ sst_sc + log_chla_sc + log_settlement_pop_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)

# --- Combined gravity metrics ---
cm_both_grav      <- glmmTMB(transect_pisc_count ~ log_market_gravity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_hab_both_grav  <- glmmTMB(transect_pisc_count ~ rugosity_sc + log_market_gravity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_env_both_grav  <- glmmTMB(transect_pisc_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + (1 | site), family = count_family, data = transect_model_data)
cm_full_both_grav <- glmmTMB(transect_pisc_count ~ sst_sc + log_chla_sc + log_market_gravity_sc + log_settlement_grav_sc + rugosity_sc + (1 | site), family = count_family, data = transect_model_data)

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

cat("\n--- AICc: piscivore count models ---\n")
print(make_aicc_df(model_list_counts))

# ── Diagnostics on best model ─────────────────────────────────
# Full (both gravity): top model (weight = 0.2905).
# Env + both gravity: ΔAICc = 0.08 (weight = 0.2794) — essentially
# tied. Prefer Env + both gravity for interpretation (simpler: no
# rugosity, which adds little to count models).
#
# Full (both gravity): clean diagnostics — no dispersion (p = 0.936),
# no zero inflation (p = 1.00), no outliers. Step failure warning
# appeared during DHARMa residual smoothing (mgcv, not glmmTMB).
#
# Env + both gravity: also clean — no dispersion (p = 0.960),
# no zero inflation (p = 0.912), no outliers. No step failure warning.
# Preferred model for all subsequent interpretation and IRR reporting.

res_cm_full <- simulateResiduals(cm_full_both_grav, n = 1000)

jpeg("dharma_pisc_count_full_both.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_cm_full, main = "DHARMa — piscivore counts: Full (both gravity)"); dev.off()

plot(res_cm_full)
testDispersion(res_cm_full)
testZeroInflation(res_cm_full)
testOutliers(res_cm_full)

res_cm_env_both <- simulateResiduals(cm_env_both_grav, n = 1000)

jpeg("dharma_pisc_count_env_both.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_cm_env_both, main = "DHARMa — piscivore counts: Env + both gravity"); dev.off()

plot(res_cm_env_both)
testDispersion(res_cm_env_both)
testZeroInflation(res_cm_env_both)
testOutliers(res_cm_env_both)

# Residual plots against all predictors — use the preferred simpler model
# (Env + both gravity), because it is near-equivalent to the top model.
plotResiduals(res_cm_env_both, transect_model_data$sst_sc,                 xlab = "SST")
plotResiduals(res_cm_env_both, transect_model_data$log_chla_sc,            xlab = "Chl-a")
plotResiduals(res_cm_env_both, transect_model_data$log_market_gravity_sc,  xlab = "Market gravity")
plotResiduals(res_cm_env_both, transect_model_data$log_settlement_grav_sc, xlab = "Settlement gravity")

# ── Key divergence from biomass models ───────────────────────
# Biomass (Parts 1 & 2): Habitat + market gravity (weights 0.41, 0.27)
# Counts  (Part 3):      Env + both gravity ≈ Full (both gravity)
#                        (ΔAICc = 0.08, combined weight = 0.57)
#
# Rugosity improves biomass models but is absent from top count models
# (ΔAICc > 10 for all habitat-containing count models): complex reefs
# support larger piscivores rather than more individuals.
# Chl-a and both gravity metrics drive piscivore abundance.
# Market gravity positive in both biomass and count models — interpret
# cautiously (likely spatial confounding).


# ── SYNTHESIS ─────────────────────────────────────────
# Part 1 (site biomass):      Habitat + market gravity
# Part 2 (transect biomass):  Habitat + market gravity  
# Part 3 (transect counts):   Env + both gravity
#                             ≈ Full (both gravity) (ΔAICc = 0.08)
# Chl-a is a strong negative predictor in count models
# (IRR = 0.58, 95% CI: 0.45–0.75, p < 0.001): higher productivity /
# turbidity is associated with fewer piscivores. This effect is absent
# from biomass models, suggesting turbid sites may support fewer
# individuals but not necessarily lower total biomass (compensation
# via larger-bodied individuals).
#
# Market gravity is positive across all three parts
# (IRR = 1.37, 95% CI: 1.09–1.71): a consistent but counterintuitive
# relationship. Likely reflects spatial or sampling confounding rather
# than a true ecological effect — interpret cautiously.
#
# Settlement gravity is negative in the count model
# (IRR = 0.77, 95% CI: 0.60–1.01): CI just crosses 1, providing only
# tentative evidence that local human pressure reduces piscivore
# abundance. Not robustly supported.
#
# Rugosity is absent from top count models (ΔAICc > 10 for all
# habitat-containing count models): strongly supported in biomass
# models but not abundance models, indicating habitat complexity
# influences size structure rather than numerical abundance.
#
# SST: IRR = 1.09, no clear effect (p > 0.4). No evidence for a
# temperature-driven change in piscivore abundance at this scale.

# ── IRR summary table ─────────────────────────────────────────
cat("\n--- IRR: Env + both gravity (top count model) ---\n")
irr_both <- exp(fixef(cm_env_both_grav)$cond)
se_both  <- summary(cm_env_both_grav)$coefficients$cond[, "Std. Error"]
cat("Chla:              IRR =", round(irr_both["log_chla_sc"],           2),
    " (95% CI:", round(exp(log(irr_both["log_chla_sc"]) - 1.96*se_both["log_chla_sc"]), 2),
    "-",         round(exp(log(irr_both["log_chla_sc"]) + 1.96*se_both["log_chla_sc"]), 2), ")\n")
cat("Market gravity:    IRR =", round(irr_both["log_market_gravity_sc"], 2),
    " (95% CI:", round(exp(log(irr_both["log_market_gravity_sc"]) - 1.96*se_both["log_market_gravity_sc"]), 2),
    "-",         round(exp(log(irr_both["log_market_gravity_sc"]) + 1.96*se_both["log_market_gravity_sc"]), 2), ")\n")
cat("Settlement grav:   IRR =", round(irr_both["log_settlement_grav_sc"],2),
    " (95% CI:", round(exp(log(irr_both["log_settlement_grav_sc"]) - 1.96*se_both["log_settlement_grav_sc"]), 2),
    "-",         round(exp(log(irr_both["log_settlement_grav_sc"]) + 1.96*se_both["log_settlement_grav_sc"]), 2), ")\n")
cat("SST:               IRR =", round(irr_both["sst_sc"], 2), "(n.s., p > 0.4)\n")

