# ============================================================
#  CORALLIVORE BIOMASS 
#  Response:   Corallivore biomass (continuous, zero-inflated)
#  Predictors: SST, Chl-a, Human gravity, Rugosity
#  Random fx:  (1 | country/site)
#  Family:     Tweedie (glmmTMB)
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
library(MASS)
library(here)

# ── DATA ─────────────────────────────────────────────────────
fish_2009      <- readr::read_rds(here("processed_data", "clean_fish_2009.rds"))
location_2009  <- readr::read_rds(here("processed_data", "clean_location_2009.rds"))
gravity_2009   <- readr::read_rds(here("city_data", "locations_with_grav_combined.rds"))
chla_2009      <- read.csv(here("processed_data", "locations_with_chla_2009.csv"))
sst_2009       <- read.csv(here("processed_data", "locations_with_sst_2009.csv"))
rugosity_2009  <- readr::read_rds(here("processed_data", "clean_dive_details_2009.rds"))

# ============================================================
#  CORALLIVORE BIOMASS PER TRANSECT
# ============================================================

coral_transects <- fish_2009 %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_coral_biomass = sum(
      ifelse(trophic_group == "corallivores", tot_wt_g, 0),
      na.rm = TRUE
    ),
    country = first(country),
    .groups = "drop"
  ) %>%
  group_by(site) %>%
  filter(n() >= 3) %>%
  ungroup() %>%
  mutate(site = as.factor(site))

cat("Transects after filtering:", nrow(coral_transects), "\n")
cat("Sites after filtering:",     n_distinct(coral_transects$site), "\n")


# ============================================================
#  CORALLIVORE BIOMASS EXPLORATION
# ============================================================

summary(coral_transects$transect_coral_biomass)

zeros <- mean(coral_transects$transect_coral_biomass == 0, na.rm = TRUE)
cat("Proportion of zeros:", round(zeros, 3), "\n")

coral_transects <- coral_transects %>%
  mutate(
    log_coral_biomass  = log(transect_coral_biomass + 0.01),
    sqrt_coral_biomass = sqrt(transect_coral_biomass)
  )

( coral_raw <- ggplot(coral_transects, aes(x = transect_coral_biomass)) +
  geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
  labs(x = "Corallivore biomass per transect (g)", y = "Frequency",
       title = "Raw Corallivore Biomass") +
  theme_bw() )

( coral_log <- ggplot(coral_transects, aes(x = log_coral_biomass)) +
  geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
  labs(x = "log(biomass + 0.01)", y = "Frequency",
       title = "Log-transformed Corallivore Biomass") +
  theme_bw() )

( coral_sqrt <- ggplot(coral_transects, aes(x = sqrt_coral_biomass)) +
  geom_histogram(bins = 30, fill = "#d7191c", colour = "white") +
  labs(x = "sqrt(biomass)", y = "Frequency",
       title = "Sqrt-transformed Corallivore Biomass") +
  theme_bw() )

jpeg("coral_biomass_distributions.jpg", width = 33, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(coral_raw, coral_log, coral_sqrt, ncol = 3)
dev.off()

# Box-Cox on non-zero values only
coral_nonzero <- coral_transects %>% filter(transect_coral_biomass > 0)
MASS::boxcox(lm(transect_coral_biomass ~ 1, data = coral_nonzero),
             lambda = seq(-2, 2, 0.1))


# ── Zeros by site ─────────────────────────────────────────────
# Are zeros spread evenly or concentrated in a few sites?
# prop_zeros = 1 suggests structural absence — biological/
# fishing signal, strengthens case for zero-inflated model.

coral_transects %>%
  group_by(site) %>%
  summarise(
    n_transects  = n(),
    prop_zeros   = mean(transect_coral_biomass == 0),
    mean_biomass = mean(transect_coral_biomass),
    .groups = "drop"
  ) %>%
  arrange(desc(prop_zeros)) %>%
  print(n = Inf)

# ── Variation by site ─────────────────────────────────────────
ggplot(coral_transects,
       aes(x = reorder(site, transect_coral_biomass, median),
           y = transect_coral_biomass)) +
  geom_boxplot(fill = "#2c7bb6", alpha = 0.6,
               outlier.colour = "black", outlier.size = 1) +
  coord_flip() +
  labs(x = NULL, y = "Corallivore biomass (g)",
       title = "Biomass distribution by site") +
  theme_bw(base_size = 9)

# ============================================================
#  PREDICTOR PREPARATION
# ============================================================

gravity_sites <- gravity_2009 %>%
  st_drop_geometry() %>%
  group_by(site) %>%
  summarise(
    market_gravity  = mean(market_grav,         na.rm = TRUE),
    settlement_pop  = mean(settlement_tot_pop,  na.rm = TRUE),
    settlement_grav = mean(nearest_pop75_grav,  na.rm = TRUE),
    .groups = "drop"
  )

chla_sites <- chla_2009 %>%
  group_by(site) %>%
  summarise(mean_annual_chla = mean(chla_annual_mean, na.rm = TRUE),
            .groups = "drop")

sst_sites <- sst_2009 %>%
  group_by(site) %>%
  summarise(mean_annual_sst = mean(sst_annual_mean, na.rm = TRUE),
            .groups = "drop")

rugosity_sites <- rugosity_2009 %>%
  group_by(site) %>%
  summarise(rugosity = mean(rugosity, na.rm = TRUE),
            .groups = "drop")

coral_model_data <- coral_transects %>%
  left_join(gravity_sites,  by = "site") %>%
  left_join(chla_sites,     by = "site") %>%
  left_join(sst_sites,      by = "site") %>%
  left_join(rugosity_sites, by = "site")

# ── Raw predictor distributions ───────────────────────────────
coral_model_data %>%
  dplyr::select(market_gravity, settlement_grav, settlement_pop,
                mean_annual_chla, mean_annual_sst, rugosity) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "#2c7bb6", colour = "white") +
  facet_wrap(~ variable, scales = "free") +
  theme_bw()

# ── Log-transform and scale ───────────────────────────────────
coral_model_data <- coral_model_data %>%
  mutate(
    log_market_gravity     = log(market_gravity),
    log_settlement_grav    = log(settlement_grav),
    log_settlement_pop     = log(settlement_pop),
    log_chla               = log(mean_annual_chla),
    
    log_market_gravity_sc  = as.numeric(scale(log_market_gravity)),
    log_settlement_grav_sc = as.numeric(scale(log_settlement_grav)),
    log_settlement_pop_sc  = as.numeric(scale(log_settlement_pop)),
    log_chla_sc            = as.numeric(scale(log_chla)),
    sst_sc                 = as.numeric(scale(mean_annual_sst)),
    rugosity_sc            = as.numeric(scale(rugosity)),
    
    log_coral_biomass      = log(transect_coral_biomass + 0.01)
  )

# Check transformations (need to reorder these to match the others...)
coral_model_data %>%
  dplyr::select(log_market_gravity, log_settlement_grav, log_settlement_pop,
                log_chla, mean_annual_sst, rugosity) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "#1a9641", colour = "white") +
  facet_wrap(~ variable, scales = "free") +
  labs(title = "Transformed predictors") +
  theme_bw()

# ── Predictor correlation matrix ─────────────────────────────
predictor_vars <- coral_model_data %>%
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

corrplot(corr_matrix,
         method = "square", type = "lower",
         tl.col = "black", tl.srt = 0, tl.offset = 0.5,
         addCoef.col = "black", number.cex = 0.8,
         col = colorRampPalette(c("#d73027", "white", "#4575b4"))(200),
         mar = c(0, 0, 4, 2))


# ============================================================
#  FAMILY SELECTION
#  Ladder: Gaussian → Tweedie → ZI Tweedie
#  Full predictor set used throughout so family differences
#  are not confounded with predictor differences.
# ============================================================

# ── F1: Gaussian on log(y + 0.01) ────────────────────────────
mF1_gaussian <- glmmTMB(
  log_coral_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = gaussian(),
  data   = coral_model_data
)

resF1 <- simulateResiduals(mF1_gaussian, n = 1000)

jpeg("dharma_F1_gaussian_coral.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF1, main = "DHARMa — Gaussian on log(y + 0.01) — corallivores")
dev.off()

plot(resF1)
testZeroInflation(resF1)
testDispersion(resF1)

# ── F1 RESULTS: Gaussian on log(y + 0.01) — REJECTED ─────────
# KS p ≈ 0; highly significant deviation despite a relatively
# low zero proportion (0.227 — lower than piscivores at 0.475).
# The severe misfit is driven primarily by extreme right-skew
# in the positive values rather than zero mass — corallivore
# biomass is highly aggregated, with most transects near zero
# and a few with very high values.
# QQ plot shows the characteristic two-cloud pattern plus
# strong curvature in the upper tail.
# ZI test returns NaN (p = 1) — Gaussian has no zero-generating
# mechanism. Dispersion fine (p = 0.946).
# → The combination of ~23% zeros and extreme positive skew
#   is sufficient to reject Gaussian. Move to Tweedie (F2).

# ── F2: Plain Tweedie ─────────────────────────────────────────
mF2_tweedie <- glmmTMB(
  transect_coral_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family = tweedie(link = "log"),
  data   = coral_model_data
)

resF2 <- simulateResiduals(mF2_tweedie, n = 1000)

jpeg("dharma_F2_tweedie_coral.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF2, main = "DHARMa — Tweedie — corallivores")
dev.off()

plot(resF2)
testZeroInflation(resF2)
testDispersion(resF2)

plotResiduals(resF2, coral_model_data$sst_sc,       xlab = "SST")
plotResiduals(resF2, coral_model_data$log_chla_sc,  xlab = "Chl-a")
plotResiduals(resF2, coral_model_data$log_market_gravity_sc, xlab = "Market gravity")
plotResiduals(resF2, coral_model_data$rugosity_sc,  xlab = "Rugosity")


# If ZI test p < 0.05 → proceed to F3.
# If diagnostics clean → Tweedie is the chosen family.

# ── F3: Zero-inflated Tweedie ─────────────────────────────────
mF3_tweedie_zi <- glmmTMB(
  transect_coral_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = coral_model_data
)

resF3 <- simulateResiduals(mF3_tweedie_zi, n = 1000)

jpeg("dharma_F3_tweedie_zi_coral.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(resF3, main = "DHARMa — ZI Tweedie — corallivores")
dev.off()

plot(resF3)
testZeroInflation(resF3)
testDispersion(resF3)

# ── AICc: F2 vs F3 ───────────────────────────────────────────
aicc_family <- MuMIn::AICc(mF2_tweedie, mF3_tweedie_zi)
aicc_family$delta  <- aicc_family$AICc - min(aicc_family$AICc)
aicc_family$weight <- round(exp(-0.5 * aicc_family$delta) /
                              sum(exp(-0.5 * aicc_family$delta)), 4)
aicc_family$AICc   <- round(aicc_family$AICc, 2)
aicc_family$delta  <- round(aicc_family$delta, 2)
print(aicc_family)
# Delta AICc > 2 in favour of F3 → ZI component justified.
# Otherwise retain plain Tweedie on parsimony grounds.

# The plain tweedie is better but still great amount of deviation! likely due to Kenya!
# ============================================================
#  RANDOM EFFECTS STRUCTURE
#  Check whether country-level variance justifies
#  nested RE (1 | country/site) over (1 | site).
# ============================================================

mF2_tweedie_nested <- glmmTMB(
  transect_coral_biomass ~ sst_sc + log_chla_sc +
    log_market_gravity_sc + rugosity_sc + (1 | country/site),
  family = tweedie(link = "log"),
  data   = coral_model_data
)

aicc_nested_vals <- MuMIn::AICc(mF2_tweedie, mF2_tweedie_nested)
aicc_nested_df <- data.frame(
  Model = rownames(aicc_nested_vals),
  df    = aicc_nested_vals$df,
  AICc  = round(aicc_nested_vals$AICc, 2),
  row.names = NULL
) %>%
  mutate(
    Delta  = round(AICc - min(AICc), 2),
    Weight = round(exp(-0.5 * Delta) / sum(exp(-0.5 * Delta)), 4)
  ) %>%
  arrange(AICc)
print(aicc_nested_df)

# Here the non-nested model is a lot better!

# Check variance components ...
VarCorr(mF2_tweedie_nested)
# If country SD is substantial relative to site SD → keep nested.
# If country SD ≈ 0 → revert to (1 | site).

# So we keep plain tweedie

# ============================================================
#  CANDIDATE MODEL FITTING
#  Update family and RE structure below based on
#  family selection and RE comparison above.
#  Placeholder: Tweedie + (1 | country/site)
# ============================================================

m0_coral <- glmmTMB(transect_coral_biomass ~ 1 + (1 | site),
                    family = tweedie(link = "log"), data = coral_model_data)

m_env_coral <- glmmTMB(transect_coral_biomass ~ sst_sc + log_chla_sc +
                         (1 | site),
                       family = tweedie(link = "log"), data = coral_model_data)

m_market_coral <- glmmTMB(transect_coral_biomass ~ log_market_gravity_sc +
                            (1 | site),
                          family = tweedie(link = "log"), data = coral_model_data)

m_settgrav_coral <- glmmTMB(transect_coral_biomass ~ log_settlement_grav_sc +
                              (1 | site),
                            family = tweedie(link = "log"), data = coral_model_data)

m_hab_coral <- glmmTMB(transect_coral_biomass ~ rugosity_sc +
                         (1 | site),
                       family = tweedie(link = "log"), data = coral_model_data)

m_env_market_coral <- glmmTMB(transect_coral_biomass ~ sst_sc + log_chla_sc +
                                log_market_gravity_sc + (1 | site),
                              family = tweedie(link = "log"), data = coral_model_data)

m_env_settgrav_coral <- glmmTMB(transect_coral_biomass ~ sst_sc + log_chla_sc +
                                  log_settlement_grav_sc + (1 | site),
                                family = tweedie(link = "log"), data = coral_model_data)

m_habitat_market <- glmmTMB(transect_coral_biomass ~ rugosity_sc + log_market_gravity_sc + (1 | site),
                            family = tweedie(link = "log"), data = coral_model_data)

m_habitat_settgrav <- glmmTMB(transect_coral_biomass ~ rugosity_sc + log_settlement_grav_sc +(1 | site),
                              family = tweedie(link = "log"), data = coral_model_data)

m_full_market_coral <- glmmTMB(transect_coral_biomass ~ sst_sc + log_chla_sc +
                                 log_market_gravity_sc + rugosity_sc +
                                 (1 | site),
                               family = tweedie(link = "log"), data = coral_model_data)

m_full_settgrav_coral <- glmmTMB(transect_coral_biomass ~ sst_sc + log_chla_sc +
                                   log_settlement_grav_sc + rugosity_sc +
                                   (1 | site),
                                 family = tweedie(link = "log"), data = coral_model_data)

# ── AICc TABLE ───────────────────────────────────────────────
model_list_coral <- list(
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

aicc_vals_coral <- sapply(model_list_coral, MuMIn::AICc)
delta_coral     <- aicc_vals_coral - min(aicc_vals_coral)
weights_coral   <- exp(-0.5 * delta_coral) / sum(exp(-0.5 * delta_coral))

aicc_df_coral <- data.frame(
  Model  = names(model_list_coral),
  AICc   = round(aicc_vals_coral, 2),
  Delta  = round(delta_coral,     2),
  Weight = round(weights_coral,   4),
  row.names = NULL
) %>% arrange(AICc)

print(aicc_df_coral)


res_m_market_coral <- simulateResiduals(m_market_coral, n = 1000)

testZeroInflation(res_m_market_coral)
testDispersion(res_m_market_coral)
testOutliers(res_m_market_coral)

plotResiduals(res_m_market_coral, coral_model_data$sst_sc,               xlab = "SST")
plotResiduals(res_m_market_coral, coral_model_data$log_chla_sc,          xlab = "Chl-a")
plotResiduals(res_m_market_coral, coral_model_data$rugosity_sc,          xlab = "Rugosity")
plotResiduals(res_m_market_coral, coral_model_data$log_market_gravity_sc,xlab = "Market gravity")
plotResiduals(res_m_market_coral, coral_model_data$country,              xlab = "Country")

summary(m_market_coral)

# ============================================================
#  MARGINAL EFFECT PLOTS
# ============================================================

plot_effect <- function(model, data, focal_var, x_label,
                        y_label = "Corallivore biomass (g)",
                        n = 200) {
  
  sc_vars <- names(data)[endsWith(names(data), "_sc")]
  
  grid <- data %>%
    slice(1) %>%
    dplyr::select(all_of(sc_vars)) %>%
    mutate(across(everything(), ~ 0)) %>%
    slice(rep(1, n)) %>%
    mutate(
      !!sym(focal_var) := seq(
        min(data[[focal_var]], na.rm = TRUE),
        max(data[[focal_var]], na.rm = TRUE),
        length.out = n
      ),
      site    = levels(data$site)[1],
      country = levels(as.factor(data$country))[1]
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


(p_grav_coral <- plot_effect(m_market_coral, coral_model_data,
                            "log_market_gravity_sc", "Market gravity (scaled log)"))
