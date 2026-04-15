# ============================================================
#  CORALLIVORE BIOMASS & ABUNDANCE — MIXED EFFECTS MODELS
#
#  Study design:
#    Transects nested within stations, stations within sites,
#    sites within locations, locations within countries.
#
#  Predictors measured at site level (averaged from station):
#    Chl-a, DHW, Human gravity (market / settlement),
#    Rugosity, Connectivity, MPA status
#
#  Analytical structure:
#    PART 1 — Site-level analysis (PRIMARY)
#              Matches response resolution to predictor resolution.
#              Sites are the true unit of environmental inference.
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
fish_data <- readr::read_rds(here::here("processed_data", "clean_fish_connectivity.rds"))
gravity_data <- readr::read_rds(here::here("city_data", "locations_with_grav_combined.rds"))
chla_data <- read.csv(here::here("processed_data", "locations_with_chla_2009.csv"))
dhw_data <- readr::read_rds(here::here("processed_data", "locations_with_dhw_2009.rds"))
rugosity_data <- readr::read_rds(here::here("processed_data", "clean_dive_details_connectivity.rds"))
location_data <- readr::read_rds(here::here("processed_data", "clean_location_connectivity.rds"))

# ── FUNCTIONS ─────────────────────────────────────────────────

make_aicc_df <- function(model_list) {
  aicc_v <- sapply(model_list, AICc)
  delta_v <- aicc_v - min(aicc_v)
  wt_v <- exp(-0.5 * delta_v) / sum(exp(-0.5 * delta_v))
  data.frame(
    Model = names(model_list),
    AICc = round(aicc_v, 2),
    Delta = round(delta_v, 2),
    Weight = round(wt_v, 4),
    row.names = NULL
  ) %>% arrange(AICc)
}

plot_effect <- function(model, data, focal_var,
                        x_label,
                        y_label = "Fitted value",
                        colour = "#2c7bb6",
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
  pred <- if (is_lm) {
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

coralliv_transects <- fish_data %>%
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
    site = as.factor(site),
    country = as.factor(country)
  )

cat("Number of transects:", nrow(coralliv_transects), "\n")
cat("Number of sites:", n_distinct(coralliv_transects$site), "\n")
cat("Number of countries:", n_distinct(coralliv_transects$country), "\n")

# ==============================================================================
#  BIOMASS DATA EXPLORATION
# ==============================================================================

site_data <- coralliv_transects %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass = mean(transect_coralliv_biomass, na.rm = TRUE),
    n_transects = n(),
    .groups = "drop"
  ) %>%
  mutate(
    site = as.factor(site),
    country = as.factor(country)
  )

summary(site_data$mean_biomass)

zeros <- mean(site_data$mean_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")

(site_raw <- ggplot(site_data, aes(x = mean_biomass)) +
    geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
    labs(x = "Mean corallivore biomass per site (g)", y = "Frequency",
         title = "Raw Site-Level Corallivore Biomass") +
    theme_bw())

site_nonzero <- site_data %>% filter(mean_biomass > 0)

MASS::boxcox(
  lm(mean_biomass ~ 1, data = site_nonzero),
  lambda = seq(-2, 2, 0.1)
)

site_data <- site_data %>%
  mutate(
    log_mean_biomass = log(mean_biomass + 0.01),
    sqrt_mean_biomass = sqrt(mean_biomass)
  )

(site_log <- ggplot(site_data, aes(x = log_mean_biomass)) +
    geom_histogram(bins = 25, fill = "#1a9641", colour = "white") +
    labs(x = "log(mean biomass + 0.01)", y = "Frequency",
         title = "Log-transformed") +
    theme_bw())

(site_sqrt <- ggplot(site_data, aes(x = sqrt_mean_biomass)) +
    geom_histogram(bins = 25, fill = "#d7191c", colour = "white") +
    labs(x = "sqrt(mean biomass)", y = "Frequency",
         title = "Sqrt-transformed") +
    theme_bw())

jpeg("site_coralliv_biomass_distributions.jpg", width = 33, height = 11,
     units = "cm", res = 300)
gridExtra::grid.arrange(site_raw, site_log, site_sqrt, ncol = 3)
dev.off()

qqnorm(site_data$log_mean_biomass,
       main = "Q-Q plot: log(mean corallivore biomass per site)")
qqline(site_data$log_mean_biomass, col = "red")
shapiro.test(site_data$log_mean_biomass)

ggplot(site_data, aes(x = reorder(site, mean_biomass, median),
                      y = mean_biomass)) +
  geom_col(fill = "#2c7bb6", alpha = 0.7) +
  coord_flip() +
  labs(x = NULL, y = "Mean corallivore biomass (g)",
       title = "Mean corallivore biomass by site") +
  theme_bw(base_size = 9)

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
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold")
  )

site_data %>%
  dplyr::select(site, country, n_transects, mean_biomass) %>%
  arrange(desc(mean_biomass)) %>%
  print(n = Inf)

# ============================================================
#  PREDICTOR PREPARATION
# ============================================================

gravity_sites <- gravity_data %>%
  st_drop_geometry() %>%
  group_by(site) %>%
  summarise(
    market_gravity = mean(market_grav, na.rm = TRUE),
    settlement_pop = mean(settlement_tot_pop, na.rm = TRUE),
    settlement_grav = mean(nearest_pop75_grav, na.rm = TRUE),
    .groups = "drop"
  )

chla_sites <- chla_data %>%
  group_by(site) %>%
  summarise(mean_annual_chla = mean(chla_annual_mean, na.rm = TRUE),
            .groups = "drop")

dhw_sites <- dhw_data %>%
  filter(!is.na(max_dhw)) %>%
  group_by(site) %>%
  summarise(max_annual_dhw = max(max_dhw, na.rm = TRUE),
            .groups = "drop")

rugosity_sites <- rugosity_data %>%
  group_by(site) %>%
  summarise(rugosity = mean(rugosity, na.rm = TRUE),
            .groups = "drop")

location_sites <- location_data %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(
    mpa_status = first(mpa_status),
    connectivity = mean(prop_connectivity, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
#  TRANSFORMATIONS AND CHECKS
# ============================================================

raw_predictors <- location_sites %>%
  left_join(chla_sites, by = "site") %>%
  left_join(dhw_sites, by = "site") %>%
  left_join(rugosity_sites, by = "site") %>%
  left_join(gravity_sites, by = "site")

transformed_predictors <- raw_predictors %>%
  transmute(
    site = site,
    log_market_gravity = log(market_gravity),
    log_settlement_grav = log(settlement_grav),
    log_settlement_pop = log(settlement_pop),
    log_chla = log(mean_annual_chla),
    log_max_dhw = log(max_annual_dhw + 1),
    rugosity = rugosity,
    connectivity = connectivity,
    mpa_status = mpa_status
  )

scaled_predictors <- transformed_predictors %>%
  transmute(
    site = site,
    rugosity_sc = as.numeric(scale(rugosity)),
    log_settlement_grav_sc = as.numeric(scale(log_settlement_grav)),
    connectivity_sc = as.numeric(scale(connectivity)),
    log_chla_sc = as.numeric(scale(log_chla)),
    log_max_dhw_sc = as.numeric(scale(log_max_dhw)),
    log_market_gravity_sc = as.numeric(scale(log_market_gravity)),
    log_settlement_pop_sc = as.numeric(scale(log_settlement_pop)),
    mpa_status = mpa_status
  )

# ============================================================
#  CHOOSING SETTLEMENT METRIC
# ============================================================

settlement_data <- coralliv_transects %>%
  left_join(scaled_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    mean_coralliv_biomass = mean(transect_coralliv_biomass, na.rm = TRUE),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_settlement_pop_sc = first(log_settlement_pop_sc),
    log_market_gravity_sc = first(log_market_gravity_sc),
    .groups = "drop"
  )

make_aicc_df(list(
  "Settlement gravity" = glmmTMB(mean_coralliv_biomass ~ log_settlement_grav_sc,
                                 family = tweedie(link = "log"), data = settlement_data),
  "Settlement pop." = glmmTMB(mean_coralliv_biomass ~ log_settlement_pop_sc,
                              family = tweedie(link = "log"), data = settlement_data),
  "Market gravity" = glmmTMB(mean_coralliv_biomass ~ log_market_gravity_sc,
                             family = tweedie(link = "log"), data = settlement_data)
))

# ── Decision ──────────────────────────────────────────────────
# Settlement gravity is the preferred human pressure metric:
#   Settlement gravity: AICc = 631.17 , weight = 0.9737(best)
#   Settlement pop.:    delta = 7.23, weight = 0.0262
#   Market gravity:     delta = 17.91 , weight = 0.0001
#
# Settlement gravity is clearly preferred. 
# Delta > 2 for both # alternatives indicates meaningful separation.
#
# Primary metric:           settlement gravity (main analyses)
# Sensitivity only:         market gravity
#                           settlement pop.   

rm(settlement_data)

# ============================================================
#  ANALYSIS DATASETS
# ============================================================

final_predictors <- scaled_predictors %>%
  dplyr::select(
    site,
    rugosity_sc,
    log_settlement_grav_sc,
    connectivity_sc,
    mpa_status,
    log_chla_sc,
    log_max_dhw_sc,
    log_settlement_pop_sc,
    log_market_gravity_sc
  )

transect_model_data <- coralliv_transects %>%
  left_join(final_predictors, by = "site") %>%
  mutate(log_coralliv_biomass = log(transect_coralliv_biomass + 0.01))

cat("\nTransect model data:", nrow(transect_model_data), "rows,",
    n_distinct(transect_model_data$site), "sites\n")
cat("Biomass zeros:", sum(transect_model_data$transect_coralliv_biomass == 0), "\n")
cat("Count zeros:", sum(transect_model_data$transect_coralliv_count == 0), "\n")

total_model_data <- transect_model_data %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass = mean(transect_coralliv_biomass, na.rm = TRUE),
    log_mean_biomass = log(mean(transect_coralliv_biomass, na.rm = TRUE) + 0.01),
    n_transects = n(),
    rugosity_sc = first(rugosity_sc),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    connectivity_sc = first(connectivity_sc),
    mpa_status = first(mpa_status),
    log_chla_sc = first(log_chla_sc),
    log_max_dhw_sc = first(log_max_dhw_sc),
    log_settlement_pop_sc = first(log_settlement_pop_sc),
    log_market_gravity_sc = first(log_market_gravity_sc),
    .groups = "drop"
  ) %>%
  mutate(
    site = as.factor(site),
    country = as.factor(country)
  )

cat("\nSite model data:", nrow(total_model_data), "sites,",
    n_distinct(total_model_data$country), "countries\n")
cat("Site-level zeros:", sum(total_model_data$mean_biomass == 0), "\n")

total_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                connectivity_sc, mpa_status, log_chla_sc,
                log_max_dhw_sc) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

# ============================================================
#  PART 1 — SITE-LEVEL ANALYSIS (PRIMARY)
# ============================================================

# ── Family selection ──────────────────────────────────────────
# Tweedie appropriate given continuous positive response with zeros.
# ZI Tweedie tested — adopt only if ZI test significant AND
# AICc improves by > 2.

m_tweedie <- glmmTMB(mean_biomass ~ rugosity_sc +
                       log_settlement_grav_sc +
                       log_chla_sc +
                       log_max_dhw_sc,
                     family = tweedie(link = "log"), data = total_model_data)

m_tweedie_zi <- glmmTMB(mean_biomass ~ rugosity_sc +
                          log_settlement_grav_sc +
                          log_chla_sc +
                          log_max_dhw_sc,
                        family = tweedie(link = "log"),
                        ziformula = ~1,
                        data = total_model_data)

res_tweedie <- simulateResiduals(m_tweedie, n = 1000)
res_tweedie_zi <- simulateResiduals(m_tweedie_zi, n = 1000)

jpeg("diagnostics_site_coralliv_tweedie.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_tweedie, main = "DHARMa — Tweedie"); dev.off()

jpeg("diagnostics_site_coralliv_tweedie_zi.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_tweedie_zi, main = "DHARMa — ZI Tweedie"); dev.off()

plot(res_tweedie); testZeroInflation(res_tweedie); testDispersion(res_tweedie)
plot(res_tweedie_zi); testZeroInflation(res_tweedie_zi); testDispersion(res_tweedie_zi)

cat("\n--- Family selection: site-level corallivore ---\n")
print(make_aicc_df(list(
  "Tweedie" = m_tweedie,
  "ZI Tweedie" = m_tweedie_zi
)))
# Update decision after running.

# ── Random effect structure ───────────────────────────────────
# Country not included — purposive sampling, no biological meaning.

re_null <- glmmTMB(mean_biomass ~ rugosity_sc +
                     log_settlement_grav_sc +
                     log_chla_sc +
                     log_max_dhw_sc,
                   family = tweedie(link = "log"), data = total_model_data)

re_country <- glmmTMB(mean_biomass ~ rugosity_sc +
                        log_settlement_grav_sc +
                        log_chla_sc +
                        log_max_dhw_sc +
                        (1 | country),
                      family = tweedie(link = "log"), data = total_model_data)

cat("\n--- RE structure comparison (site-level corallivore) ---\n")
print(make_aicc_df(list(
  "No RE" = re_null,
  "(1 | country)" = re_country
)))
# Update decision after running.

# ── Candidate models ──────────────────────────────────────────
# Same ladder as browser analyses for comparability.
# Family: update after family selection decision above.

site_family <- tweedie(link = "log") # update if ZI selected

m1_hab <- glmmTMB(mean_biomass ~ rugosity_sc,
                  family = site_family, data = total_model_data)

m2_hab_press <- glmmTMB(mean_biomass ~ rugosity_sc +
                          log_settlement_grav_sc,
                        family = site_family, data = total_model_data)

m3_hab_press_mpa <- glmmTMB(mean_biomass ~ rugosity_sc +
                              log_settlement_grav_sc +
                              mpa_status,
                            family = site_family, data = total_model_data)

m4_conn <- glmmTMB(mean_biomass ~ rugosity_sc +
                     log_settlement_grav_sc +
                     mpa_status +
                     connectivity_sc,
                   family = site_family, data = total_model_data)

m5_chla <- glmmTMB(mean_biomass ~ rugosity_sc +
                     log_settlement_grav_sc +
                     mpa_status +
                     connectivity_sc +
                     log_chla_sc,
                   family = site_family, data = total_model_data)

m6_dhw <- glmmTMB(mean_biomass ~ rugosity_sc +
                    log_settlement_grav_sc +
                    mpa_status +
                    connectivity_sc +
                    log_max_dhw_sc,
                  family = site_family, data = total_model_data)

m7_mpa_conn <- glmmTMB(mean_biomass ~ rugosity_sc +
                         log_settlement_grav_sc +
                         mpa_status * connectivity_sc,
                       family = site_family, data = total_model_data)

m8_mpa_press <- glmmTMB(mean_biomass ~ rugosity_sc +
                          mpa_status * log_settlement_grav_sc +
                          connectivity_sc,
                        family = site_family, data = total_model_data)

m9_conn_press <- glmmTMB(mean_biomass ~ rugosity_sc +
                           mpa_status +
                           connectivity_sc * log_settlement_grav_sc,
                         family = site_family, data = total_model_data)

sens_settpop <- glmmTMB(mean_biomass ~ rugosity_sc +
                          log_settlement_pop_sc +
                          mpa_status +
                          connectivity_sc,
                        family = site_family, data = total_model_data)

sens_mktgrav <- glmmTMB(mean_biomass ~ rugosity_sc +
                          log_market_gravity_sc +
                          mpa_status +
                          connectivity_sc,
                        family = site_family, data = total_model_data)

model_list_s1 <- list(
  "Habitat" = m1_hab,
  "Habitat + pressure" = m2_hab_press,
  "Habitat + pressure + MPA" = m3_hab_press_mpa,
  "Above + connectivity" = m4_conn,
  "Above + chla" = m5_chla,
  "Above + DHW" = m6_dhw,
  "MPA x connectivity" = m7_mpa_conn,
  "MPA x pressure" = m8_mpa_press,
  "Connectivity x pressure" = m9_conn_press,
  "Settlement pop. (sensitivity)" = sens_settpop,
  "Market gravity (sensitivity)" = sens_mktgrav
)

cat("\n--- AICc: Site-level corallivore candidate models ---\n")
print(make_aicc_df(model_list_s1))

# ── Diagnostics — models within delta < 2 ────────────────────
# Update after running — run DHARMa on top models

# ── Marginal effect plots ─────────────────────────────────────
# Update after running — use best model for continuous predictors
# Use interaction models for MPA x pressure / MPA x connectivity

# ============================================================
#  PART 2 — TRANSECT-LEVEL BIOMASS (SENSITIVITY CHECK)
# ============================================================

summary(coralliv_transects$transect_coralliv_biomass)

zeros <- mean(coralliv_transects$transect_coralliv_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")

(coralliv_raw <- ggplot(coralliv_transects, aes(x = transect_coralliv_biomass)) +
    geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
    labs(x = "Corallivore biomass per transect (g)", y = "Frequency",
         title = "Raw") + theme_bw())

(coralliv_log <- ggplot(transect_model_data, aes(x = log_coralliv_biomass)) +
    geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
    labs(x = "log(biomass + 0.01)", y = "Frequency",
         title = "Log-transformed") + theme_bw())

jpeg("coralliv_biomass_distributions.jpg", width = 22, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(coralliv_raw, coralliv_log, ncol = 2)
dev.off()

coralliv_nonzero <- coralliv_transects %>% filter(transect_coralliv_biomass > 0)

MASS::boxcox(
  lm(transect_coralliv_biomass ~ 1, data = coralliv_nonzero),
  lambda = seq(-2, 2, 0.1)
)

ggplot(transect_model_data,
       aes(x = reorder(site, transect_coralliv_biomass, median),
           y = transect_coralliv_biomass)) +
  geom_boxplot(fill = "#2c7bb6", alpha = 0.6,
               outlier.colour = "black", outlier.size = 1) +
  coord_flip() +
  labs(x = NULL, y = "Corallivore biomass (g)",
       title = "Corallivore biomass by site") +
  theme_bw(base_size = 9)

transect_model_data %>%
  group_by(site) %>%
  summarise(
    n_transects = n(),
    prop_zeros = mean(transect_coralliv_biomass == 0),
    mean_biomass = mean(transect_coralliv_biomass),
    .groups = "drop"
  ) %>%
  arrange(desc(prop_zeros)) %>%
  print(n = Inf)

# ── Family selection ──────────────────────────────────────────
# Zero proportion determines transformation strategy.
# Update comment after checking zeros above.
# ZI Tweedie adopted only if: (1) ZI test significant AND
# (2) AICc improves by > 2.

transect_m_tweedie <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    log_max_dhw_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data = transect_model_data
)

transect_m_tweedie_zi <- glmmTMB(
  transect_coralliv_biomass ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    log_max_dhw_sc +
    (1 | site),
  family = tweedie(link = "log"),
  ziformula = ~1,
  data = transect_model_data
)

transect_res_tweedie <- simulateResiduals(transect_m_tweedie, n = 1000)
transect_res_tweedie_zi <- simulateResiduals(transect_m_tweedie_zi, n = 1000)

jpeg("dharma_coralliv_transect_tweedie.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(transect_res_tweedie, main = "DHARMa — Tweedie"); dev.off()

jpeg("dharma_coralliv_transect_tweedie_zi.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(transect_res_tweedie_zi, main = "DHARMa — ZI Tweedie"); dev.off()

plot(transect_res_tweedie)
testZeroInflation(transect_res_tweedie)
testDispersion(transect_res_tweedie)

plot(transect_res_tweedie_zi)
testZeroInflation(transect_res_tweedie_zi)
testDispersion(transect_res_tweedie_zi)

cat("\n--- Family selection: transect-level corallivore biomass ---\n")
print(make_aicc_df(list(
  "Tweedie" = transect_m_tweedie,
  "ZI Tweedie" = transect_m_tweedie_zi
)))
# Update decision after running.

# ── Random effect structure ───────────────────────────────────

transect_re_null <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                              log_settlement_grav_sc +
                              log_chla_sc +
                              log_max_dhw_sc,
                            family = tweedie(link = "log"),
                            data = transect_model_data)

transect_re_site <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                              log_settlement_grav_sc +
                              log_chla_sc +
                              log_max_dhw_sc +
                              (1 | site),
                            family = tweedie(link = "log"),
                            data = transect_model_data)

cat("\n--- RE structure comparison (transect-level corallivore) ---\n")
print(make_aicc_df(list(
  "No RE" = transect_re_null,
  "(1 | site)" = transect_re_site
)))
# Update decision after running.

# ── Candidate models ──────────────────────────────────────────

coralliv_family <- tweedie(link = "log") # update if ZI selected

transect_m1_hab <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                             (1 | site),
                           family = coralliv_family, data = transect_model_data)

transect_m2_hab_press <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                                   log_settlement_grav_sc +
                                   (1 | site),
                                 family = coralliv_family, data = transect_model_data)

transect_m3_hab_press_mpa <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                                       log_settlement_grav_sc +
                                       mpa_status +
                                       (1 | site),
                                     family = coralliv_family, data = transect_model_data)

transect_m4_conn <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                              log_settlement_grav_sc +
                              mpa_status +
                              connectivity_sc +
                              (1 | site),
                            family = coralliv_family, data = transect_model_data)

transect_m5_chla <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                              log_settlement_grav_sc +
                              mpa_status +
                              connectivity_sc +
                              log_chla_sc +
                              (1 | site),
                            family = coralliv_family, data = transect_model_data)

transect_m6_dhw <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                             log_settlement_grav_sc +
                             mpa_status +
                             connectivity_sc +
                             log_max_dhw_sc +
                             (1 | site),
                           family = coralliv_family, data = transect_model_data)

transect_m7_mpa_conn <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                                  log_settlement_grav_sc +
                                  mpa_status * connectivity_sc +
                                  (1 | site),
                                family = coralliv_family, data = transect_model_data)

transect_m8_mpa_press <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                                   mpa_status * log_settlement_grav_sc +
                                   connectivity_sc +
                                   (1 | site),
                                 family = coralliv_family, data = transect_model_data)

transect_m9_conn_press <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                                    mpa_status +
                                    connectivity_sc * log_settlement_grav_sc +
                                    (1 | site),
                                  family = coralliv_family, data = transect_model_data)

transect_sens_settpop <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                                   log_settlement_pop_sc +
                                   mpa_status +
                                   connectivity_sc +
                                   (1 | site),
                                 family = coralliv_family, data = transect_model_data)

transect_sens_mktgrav <- glmmTMB(transect_coralliv_biomass ~ rugosity_sc +
                                   log_market_gravity_sc +
                                   mpa_status +
                                   connectivity_sc +
                                   (1 | site),
                                 family = coralliv_family, data = transect_model_data)

model_list_transect <- list(
  "Habitat" = transect_m1_hab,
  "Habitat + pressure" = transect_m2_hab_press,
  "Habitat + pressure + MPA" = transect_m3_hab_press_mpa,
  "Above + connectivity" = transect_m4_conn,
  "Above + chla" = transect_m5_chla,
  "Above + DHW" = transect_m6_dhw,
  "MPA x connectivity" = transect_m7_mpa_conn,
  "MPA x pressure" = transect_m8_mpa_press,
  "Connectivity x pressure" = transect_m9_conn_press,
  "Settlement pop. (sensitivity)" = transect_sens_settpop,
  "Market gravity (sensitivity)" = transect_sens_mktgrav
)

cat("\n--- AICc: Transect-level corallivore biomass ---\n")
print(make_aicc_df(model_list_transect))

# ── Results — update after running ───────────────────────────
# Key question: does transect-level analysis converge with site-level?

# ============================================================
#  PART 3 — TRANSECT-LEVEL COUNTS (COMPLEMENTARY ANALYSIS)
# ============================================================

cat("Transects:", nrow(transect_model_data), "\n")
cat("Zero counts:", sum(transect_model_data$transect_coralliv_count == 0), "\n")
cat("Proportion zeros:", round(mean(transect_model_data$transect_coralliv_count == 0), 3), "\n")

summary(transect_model_data$transect_coralliv_count)

ggplot(transect_model_data, aes(x = transect_coralliv_count)) +
  geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
  labs(x = "Total corallivore count per transect", y = "Frequency") +
  theme_bw()

transect_model_data %>%
  group_by(site) %>%
  summarise(mean_count = mean(transect_coralliv_count),
            var_count = var(transect_coralliv_count),
            .groups = "drop") %>%
  ggplot(aes(x = mean_count, y = var_count)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Site mean count", y = "Site variance",
       title = "Mean-variance (red = Poisson expectation)") +
  theme_bw()

# ── Family selection ──────────────────────────────────────────

m_count_poisson <- glmmTMB(
  transect_coralliv_count ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    log_max_dhw_sc +
    (1 | site),
  family = poisson(link = "log"),
  data = transect_model_data
)

m_count_nb2 <- glmmTMB(
  transect_coralliv_count ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    log_max_dhw_sc +
    (1 | site),
  family = nbinom2(link = "log"),
  data = transect_model_data
)

m_count_nb1 <- glmmTMB(
  transect_coralliv_count ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    log_max_dhw_sc +
    (1 | site),
  family = nbinom1(link = "log"),
  data = transect_model_data
)

res_poisson <- simulateResiduals(m_count_poisson, n = 1000)
res_nb2 <- simulateResiduals(m_count_nb2, n = 1000)
res_nb1 <- simulateResiduals(m_count_nb1, n = 1000)

jpeg("dharma_coralliv_count_poisson.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_poisson, main = "DHARMa — Poisson"); dev.off()

jpeg("dharma_coralliv_count_nb2.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_nb2, main = "DHARMa — NB2"); dev.off()

jpeg("dharma_coralliv_count_nb1.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_nb1, main = "DHARMa — NB1"); dev.off()

plot(res_poisson); testDispersion(res_poisson); testZeroInflation(res_poisson)
plot(res_nb2); testDispersion(res_nb2); testZeroInflation(res_nb2)
plot(res_nb1); testDispersion(res_nb1); testZeroInflation(res_nb1)

cat("\n--- Family selection: corallivore count models ---\n")
print(make_aicc_df(list(
  "Poisson" = m_count_poisson,
  "NB2" = m_count_nb2,
  "NB1" = m_count_nb1
)))
# Update decision after running.

# ── Random effect structure ───────────────────────────────────

count_family <- nbinom2(link = "log") # update after family selection

re_c_null <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                       log_settlement_grav_sc +
                       log_chla_sc +
                       log_max_dhw_sc,
                     family = count_family, data = transect_model_data)

re_c_site <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                       log_settlement_grav_sc +
                       log_chla_sc +
                       log_max_dhw_sc +
                       (1 | site),
                     family = count_family, data = transect_model_data)

cat("\n--- RE structure: corallivore count models ---\n")
print(make_aicc_df(list(
  "No RE" = re_c_null,
  "(1 | site)" = re_c_site
)))
# (1 | site) expected to be supported — retain for consistency
# with all other count analyses regardless of delta.

# ── Candidate models ──────────────────────────────────────────

c_m1_hab <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                      (1 | site),
                    family = count_family, data = transect_model_data)

c_m2_hab_press <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                            log_settlement_grav_sc +
                            (1 | site),
                          family = count_family, data = transect_model_data)

c_m3_hab_press_mpa <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                                log_settlement_grav_sc +
                                mpa_status +
                                (1 | site),
                              family = count_family, data = transect_model_data)

c_m4_conn <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                       log_settlement_grav_sc +
                       mpa_status +
                       connectivity_sc +
                       (1 | site),
                     family = count_family, data = transect_model_data)

c_m5_chla <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                       log_settlement_grav_sc +
                       mpa_status +
                       connectivity_sc +
                       log_chla_sc +
                       (1 | site),
                     family = count_family, data = transect_model_data)

c_m6_dhw <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                      log_settlement_grav_sc +
                      mpa_status +
                      connectivity_sc +
                      log_max_dhw_sc +
                      (1 | site),
                    family = count_family, data = transect_model_data)

c_m7_mpa_conn <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                           log_settlement_grav_sc +
                           mpa_status * connectivity_sc +
                           (1 | site),
                         family = count_family, data = transect_model_data)

c_m8_mpa_press <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                            mpa_status * log_settlement_grav_sc +
                            connectivity_sc +
                            (1 | site),
                          family = count_family, data = transect_model_data)

c_m9_conn_press <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                             mpa_status +
                             connectivity_sc * log_settlement_grav_sc +
                             (1 | site),
                           family = count_family, data = transect_model_data)

c_sens_settpop <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                            log_settlement_pop_sc +
                            mpa_status +
                            connectivity_sc +
                            (1 | site),
                          family = count_family, data = transect_model_data)

c_sens_mktgrav <- glmmTMB(transect_coralliv_count ~ rugosity_sc +
                            log_market_gravity_sc +
                            mpa_status +
                            connectivity_sc +
                            (1 | site),
                          family = count_family, data = transect_model_data)

model_list_counts <- list(
  "Habitat" = c_m1_hab,
  "Habitat + pressure" = c_m2_hab_press,
  "Habitat + pressure + MPA" = c_m3_hab_press_mpa,
  "Above + connectivity" = c_m4_conn,
  "Above + chla" = c_m5_chla,
  "Above + DHW" = c_m6_dhw,
  "MPA x connectivity" = c_m7_mpa_conn,
  "MPA x pressure" = c_m8_mpa_press,
  "Connectivity x pressure" = c_m9_conn_press,
  "Settlement pop. (sensitivity)" = c_sens_settpop,
  "Market gravity (sensitivity)" = c_sens_mktgrav
)

cat("\n--- AICc: Corallivore count models ---\n")
print(make_aicc_df(model_list_counts))

# ── Results — update after running ───────────────────────────
# Key comparison: do count results converge with biomass?
# Divergence indicates whether predictors operate through
# fish abundance, body size, or both.

# ============================================================
#  SYNTHESIS: CORALLIVORE BIOMASS vs COUNT CONCLUSIONS
# ============================================================
# Update after running all three parts.
# Key questions:
#   1. Is the null model dominant as in the old analysis?
#   2. Do MPA interactions appear for corallivores?
#   3. Does DHW affect corallivore counts?
#   4. How do corallivores differ from browsers?

cat("\n=== CORALLIVORE — Site-level biomass, best model ===\n")
# summary(best model — update after running)

cat("\n=== CORALLIVORE — Transect-level biomass, best model ===\n")
# summary(best model — update after running)

cat("\n=== CORALLIVORE — Count models, top models ===\n")
# summary(best model — update after running)