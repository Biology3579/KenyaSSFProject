# ============================================================
#  DRIVERS OF REEF FISH BIOMASS
#  Chapter 1 — Site-level Analysis
#
#  Scientific questions:
#
#  Q1 — Pressure metric validity
#       Does settlement gravity outperform market gravity and
#       settlement population as a proxy for small-scale
#       fisheries pressure on reef fish biomass?
#       Tested via multivariate AICc comparison — three models
#       identical in structure, differing only in pressure
#       metric. Robustness confirmed with and without
#       connectivity as a control.
#
#  Q2 — Are connectivity and MPA baseline drivers?
#       Does connectivity explain biomass variation
#       independently of human pressure and ecological context,
#       extending Warmuth et al. (2024) by including fishing
#       pressure as a covariate and using biomass rather than
#       abundance as the response?
#       Does MPA status explain additional variance beyond
#       the correctly-specified local baseline?
#       Tested via hierarchical model comparison — baseline
#       model extended with connectivity and MPA separately
#       and together.
#
#  Q3 — Do management and connectivity modify pressure effects?
#       Do MPA status and connectivity change the relationship
#       between human pressure and biomass?
#       Three a priori interaction hypotheses tested separately
#       against best Q2 model. All run regardless of Q2
#       main-effect results — absence of main effect does not
#       preclude interaction. Interpretation conditional on Q2.
#
#  Baseline model (fixed a priori for all analyses):
#       log(biomass) ~ rugosity + settlement_gravity + chla
#       — rugosity: sets local carrying capacity
#       — settlement gravity: SSF exploitation pressure
#       — chla: background productivity (Warmuth et al. 2024,
#               Samoilys et al. 2025)
#
#  Sensitivity analysis:
#       (a) Alternative pressure metrics — confirms Q1/Q2
#           conclusions not metric-dependent
#       (b) Transect-level GLMM — confirms site-level findings
#           not an artefact of aggregation
#
#  Study design:
#       Transects nested within stations, stations within
#       sites, sites within locations, locations within
#       countries.
#
#  Response variable:
#       log(mean total fish biomass per site) [primary]
#       log(transect biomass) [sensitivity, lmer]
#       Biomass used throughout — integrates abundance and
#       body size, more sensitive to exploitation than counts
#       alone (fishing reduces mean body size before
#       detectably reducing abundance).
# ============================================================

options(scipen = 999)

# ── PACKAGES ──────────────────────────────────────────────────
library(tidyverse)
library(sf)
library(MuMIn)
library(gridExtra)
library(here)
library(spdep)
library(ggplot2)
library(lme4)
library(car)

# ── FUNCTIONS ─────────────────────────────────────────────────

# AICc comparison table from a named list of models
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

# Marginal effect plot for a single focal predictor
# All other scaled predictors held at 0 (their mean)
# Observed data points overlaid for transparency
plot_effect <- function(model, data, focal_var,
                        x_label,
                        y_label = "Fitted log(biomass)",
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
  is_lm <- inherits(model, "lm")
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
                fill = colour, alpha = 0.15) +
    geom_line(aes(y = fit), colour = colour, linewidth = 1.1) +
    geom_point(data = data,
               aes(x = .data[[focal_var]],
                   y = log_mean_biomass),
               colour = "grey40", size = 1.5,
               alpha  = 0.5, inherit.aes = FALSE) +
    labs(x = x_label, y = y_label) +
    theme_bw(base_size = 13) +
    theme(axis.title      = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}


# ============================================================
#  DATA LOADING
# ============================================================

fish_data     <- readr::read_rds(here("processed_data",
                                      "clean_fish_connectivity.rds"))
gravity_data  <- readr::read_rds(here("city_data",
                                      "locations_with_grav_combined.rds"))
chla_data     <- read.csv(here("processed_data",
                               "locations_with_chla_2009.csv"))
rugosity_data <- readr::read_rds(here("processed_data",
                                      "clean_dive_details_connectivity.rds"))
location_data <- readr::read_rds(here("processed_data",
                                      "clean_location_connectivity.rds"))


# ============================================================
#  DATA AGGREGATION
#  Minimum 3 transects per site retained.
# ============================================================

total_transects <- fish_data %>%
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
    site    = as.factor(site),
    country = as.factor(country)
  )

cat("Transects:", nrow(total_transects), "\n")
cat("Sites:",     n_distinct(total_transects$site), "\n")
cat("Countries:", n_distinct(total_transects$country), "\n")


# ============================================================
#  PREDICTOR PREPARATION
# ============================================================

# ── Human gravity metrics ─────────────────────────────────────
# Three metrics proxy small-scale fisheries pressure.
# All retained here — primary metric selected in Q1.
gravity_sites <- gravity_data %>%
  st_drop_geometry() %>%
  group_by(site) %>%
  summarise(
    market_gravity  = mean(market_grav,        na.rm = TRUE),
    settlement_pop  = mean(settlement_tot_pop, na.rm = TRUE),
    settlement_grav = mean(nearest_pop75_grav, na.rm = TRUE),
    .groups = "drop"
  )

# ── Chlorophyll-a ─────────────────────────────────────────────
# Baseline covariate: background productivity.
# Standard environmental control in reef fish models
# (Warmuth et al. 2024, Samoilys et al. 2025).
chla_sites <- chla_data %>%
  group_by(site) %>%
  summarise(mean_annual_chla = mean(chla_annual_mean, na.rm = TRUE),
            .groups = "drop")

# ── Rugosity ──────────────────────────────────────────────────
# Baseline covariate: habitat structural complexity.
# Sets local carrying capacity for reef fish
# (Samoilys et al. 2025, Darling et al. 2017).
rugosity_sites <- rugosity_data %>%
  group_by(site) %>%
  summarise(rugosity = mean(rugosity, na.rm = TRUE),
            .groups = "drop")

# ── MPA status and connectivity ───────────────────────────────
# MPA: governance modifier — tested in Q2 and Q3 only,
#   not part of baseline (operates through pressure reduction,
#   not as an independent biomass driver).
# Connectivity: candidate driver tested in Q2 (main effect)
#   and Q3 (interaction with pressure and MPA).
location_sites <- location_data %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(
    mpa_status   = first(mpa_status),
    connectivity = mean(prop_connectivity, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
#  TRANSFORMATIONS AND SCALING
# ============================================================

raw_predictors <- location_sites %>%
  left_join(chla_sites,     by = "site") %>%
  left_join(rugosity_sites, by = "site") %>%
  left_join(gravity_sites,  by = "site")

# ── Transformations ───────────────────────────────────────────
# Rugosity:        no transformation (approximately normal)
# Gravity metrics: log (right-skewed)
# Chla:            log (right-skewed)
# Connectivity:    no transformation
# MPA status:      unordered factor (none / low / medium)

transformed_predictors <- raw_predictors %>%
  transmute(
    site                = site,
    rugosity            = rugosity,
    log_settlement_grav = log(settlement_grav),
    log_settlement_pop  = log(settlement_pop),
    log_market_gravity  = log(market_gravity),
    log_chla            = log(mean_annual_chla),
    connectivity        = connectivity,
    mpa_status          = factor(mpa_status,
                                 levels  = c("none", "low", "medium"),
                                 ordered = FALSE)
  )

# ── Standardise continuous predictors ────────────────────────
# z-score scaling (mean = 0, SD = 1) enables direct comparison
# of effect sizes across predictors with different units.

scaled_predictors <- transformed_predictors %>%
  transmute(
    site                   = site,
    rugosity_sc            = as.numeric(scale(rugosity)),
    log_settlement_grav_sc = as.numeric(scale(log_settlement_grav)),
    log_settlement_pop_sc  = as.numeric(scale(log_settlement_pop)),
    log_market_gravity_sc  = as.numeric(scale(log_market_gravity)),
    log_chla_sc            = as.numeric(scale(log_chla)),
    connectivity_sc        = as.numeric(scale(connectivity)),
    mpa_status             = mpa_status
  )


# ============================================================
#  PREDICTOR CORRELATION CHECK
#  Check for blocking collinearity before modelling.
#  Rule of thumb: |r| > 0.70 warrants caution.
# ============================================================

corr_matrix <- scaled_predictors %>%
  mutate(mpa_numeric = as.numeric(mpa_status)) %>%
  dplyr::select(ends_with("_sc"), mpa_numeric) %>%
  rename(
    "Market gravity"     = log_market_gravity_sc,
    "Settlement gravity" = log_settlement_grav_sc,
    "Settlement pop."    = log_settlement_pop_sc,
    "Chlorophyll-a"      = log_chla_sc,
    "Rugosity"           = rugosity_sc,
    "Connectivity"       = connectivity_sc,
    "MPA status"         = mpa_numeric
  ) %>%
  cor(use = "complete.obs")

vars <- colnames(corr_matrix)

corr_long <- as.data.frame(as.table(corr_matrix)) %>%
  rename(x = Var1, y = Var2, corr = Freq) %>%
  mutate(
    x = factor(x, levels = vars),
    y = factor(y, levels = rev(vars))
  ) %>%
  filter(as.numeric(x) < length(vars) - as.numeric(y) + 1)

ggplot(corr_long, aes(x = x, y = y, fill = corr)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", corr)),
            size = 3.2, color = "black") +
  scale_fill_gradient2(
    low    = "#2c7bb6",
    mid    = "white",
    high   = "#d7191c",
    limits = c(-1, 1),
    name   = "Correlation"
  ) +
  scale_x_discrete(limits = vars, drop = FALSE,
                   position = "bottom") +
  scale_y_discrete(limits = rev(vars), drop = FALSE,
                   position = "left") +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid        = element_blank(),
    axis.title        = element_blank(),
    axis.text.x       = element_text(angle = 45, hjust = 1,
                                     vjust = 1, size = 10),
    axis.text.y       = element_text(hjust = 1, size = 9.5),
    legend.position   = c(1.015, 0.48),
    legend.key.height = unit(1.75, "cm"),
    legend.key.width  = unit(0.5,  "cm"),
    legend.title      = element_text(size = 9.75, hjust = 1.5,
                                     vjust = 1.5, colour = "grey40"),
    legend.text       = element_text(size = 9, vjust = 0.3,
                                     colour = "grey40"),
    plot.margin       = margin(10, 10, 10, 10)
  )

# ── Collinearity notes ────────────────────────────────────────
#
# GRAVITY METRICS (settlement gravity / pop / market): r = 0.53–0.54
#   Moderate positive — all proxy the same construct.
#   One per model only. Primary metric selected in Q1.
#
# CHLA vs settlement gravity: r = -0.57
#   Strongest pairwise correlation in the dataset.
#   Negative — productive sites tend to be more remote and
#   less fished. Geographic covariation, not collinearity
#   that blocks inference. VIFs confirm no meaningful
#   variance inflation (chla VIF = 1.43, settlement
#   gravity VIF = 1.44). Monitor chla coefficient in
#   models including settlement gravity simultaneously.
#
# MPA vs settlement gravity: r = -0.30
#   Negative — MPAs preferentially placed in lower-pressure
#   areas. Expected and acceptable. Noted when interpreting
#   MPA effects in Q2 and Q3.
#
# MPA vs settlement pop.: r = -0.41
#   Stronger negative association with population-based
#   metric than with settlement gravity — consistent with
#   MPAs placed away from densely populated areas.
#
# CONNECTIVITY vs market gravity: r = 0.30
#   Weak positive — more connected sites tend to be near
#   larger markets. Acceptable; monitor in models including
#   both terms.
#
# CONNECTIVITY vs chla: r = -0.30
#   Weak negative — more connected sites tend to be less
#   productive. No inference concern.
#
# RUGOSITY: all |r| ≤ 0.19 — orthogonal to all predictors.
# CONNECTIVITY: all |r| ≤ 0.34 — no blocking collinearity.
#
# No blocking collinearity. All pairwise |r| < 0.60.
# Gravity metrics handled by single-metric-per-model rule.


# ============================================================
#  ANALYSIS DATASETS
# ============================================================

final_predictors <- scaled_predictors %>%
  dplyr::select(
    site,
    rugosity_sc,                # BASELINE    — habitat complexity
    log_settlement_grav_sc,     # Q1/Q2/Q3   — primary pressure metric
    log_chla_sc,                # BASELINE    — productivity
    connectivity_sc,            # Q2/Q3      — larval connectivity
    mpa_status,                 # Q2/Q3      — governance modifier
    log_settlement_pop_sc,      # SENSITIVITY
    log_market_gravity_sc       # SENSITIVITY
  )

# ── Transect-level dataset ────────────────────────────────────
transect_model_data <- total_transects %>%
  left_join(final_predictors, by = "site") %>%
  mutate(log_transect_biomass = log(transect_total_biomass))

cat("\nTransect data:", nrow(transect_model_data), "rows,",
    n_distinct(transect_model_data$site), "sites\n")
cat("Biomass zeros:",
    sum(transect_model_data$transect_total_biomass == 0), "\n")

# ── Site-level dataset ────────────────────────────────────────
total_model_data <- transect_model_data %>%
  group_by(site, country) %>%
  summarise(
    log_mean_biomass       = log(mean(transect_total_biomass,
                                      na.rm = TRUE)),
    mean_biomass           = mean(transect_total_biomass,
                                  na.rm = TRUE),
    n_transects            = n(),
    rugosity_sc            = first(rugosity_sc),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_chla_sc            = first(log_chla_sc),
    connectivity_sc        = first(connectivity_sc),
    mpa_status             = first(mpa_status),
    log_settlement_pop_sc  = first(log_settlement_pop_sc),
    log_market_gravity_sc  = first(log_market_gravity_sc),
    .groups = "drop"
  ) %>%
  mutate(
    site    = as.factor(site),
    country = as.factor(country)
  )

cat("\nSite data:", nrow(total_model_data), "sites,",
    n_distinct(total_model_data$country), "countries\n")

# ── Data checks ───────────────────────────────────────────────
total_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                log_chla_sc, connectivity_sc, mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

cat("\nZeros in mean_biomass:",
    sum(total_model_data$mean_biomass == 0), "\n")
cat("-Inf in log_mean_biomass:",
    sum(is.infinite(total_model_data$log_mean_biomass)), "\n")
cat("NAs in log_mean_biomass:",
    sum(is.na(total_model_data$log_mean_biomass)), "\n")
cat("\nResponse variable summary:\n")
print(summary(total_model_data$log_mean_biomass))

# ── MPA classification check ──────────────────────────────────
# Verify MPA counts after any raw data corrections.
cat("\nMPA status counts:\n")
print(table(total_model_data$mpa_status))


# ── Variance inflation factors ────────────────────────────────
# Checked on baseline model — most demanding collinearity
# case for the primary predictor set. Predictors are
# identical across all response variables (total biomass
# and functional groups), so one check covers all analyses.
# Rule of thumb: VIF > 5 warrants concern, > 10 problematic.

cat("\n--- VIF: baseline model ---\n")
vif(lm(log_mean_biomass ~ rugosity_sc +
         log_settlement_grav_sc +
         log_chla_sc,
       data = total_model_data))

# rugosity:           VIF = 1.01
# settlement gravity: VIF = 1.44
# chla:               VIF = 1.43
# All VIFs < 2 — no multicollinearity concern.
# The moderate negative correlation between chla and
# settlement gravity (r = -0.57) does not translate into
# meaningful variance inflation in the multivariate model.


# ============================================================
#  MODEL FAMILY SELECTION
#  Run once on global model — if distribution holds with all
#  predictors, it holds for all reduced models.
#
#  Gaussian (raw):   rejected — severe heteroscedasticity,
#                    heavy upper Q-Q tail, sites 6 and 48
#                    influential.
#
#  Gamma (log link): rejected — systematic Q-Q deviation,
#                    sites 6, 27, 35 pulling from half-normal
#                    line. Diagnostic failure outweighs better
#                    raw-scale AICc (1153.73 vs 1197.84).
#
#  Gaussian (log):   SELECTED — flat residuals, Q-Q closely
#                    follows theoretical line, no influential
#                    sites. Not AICc-comparable (transformed
#                    response); selected on diagnostics.
#
#  Proceed: lm() on log_mean_biomass throughout.
# ============================================================

lm_gaussian_raw <- lm(mean_biomass ~ rugosity_sc +
                        log_settlement_grav_sc +
                        log_chla_sc +
                        connectivity_sc +
                        mpa_status,
                      data = total_model_data)

lm_gaussian_log <- lm(log_mean_biomass ~ rugosity_sc +
                        log_settlement_grav_sc +
                        log_chla_sc +
                        connectivity_sc +
                        mpa_status,
                      data = total_model_data)

glm_gamma <- glm(mean_biomass ~ rugosity_sc +
                   log_settlement_grav_sc +
                   log_chla_sc +
                   connectivity_sc +
                   mpa_status,
                 family = Gamma(link = "log"),
                 data   = total_model_data)

par(mfrow = c(2, 2))
plot(lm_gaussian_raw, main = "Gaussian raw")
plot(lm_gaussian_log, main = "Gaussian log")
plot(glm_gamma,       main = "Gamma log-link")
par(mfrow = c(1, 1))

cat("\n--- Family selection: AICc (raw scale only) ---\n")
print(make_aicc_df(list(
  "Gaussian (raw)" = lm_gaussian_raw,
  "Gamma"          = glm_gamma
)))


# ============================================================
#  RANDOM EFFECT STRUCTURE
#  Tests whether country-level clustering requires a random
#  intercept once all predictors are included.
#  Global predictor set used deliberately — the most
#  demanding test. If country RE not supported here,
#  it will not be supported in any reduced model.
#  Both models fitted by ML (REML = FALSE) for AICc
#  comparison.
# ============================================================

re_null <- lm(log_mean_biomass ~ rugosity_sc +
                log_settlement_grav_sc +
                log_chla_sc +
                connectivity_sc +
                mpa_status,
              data = total_model_data)

re_country <- lmer(log_mean_biomass ~ rugosity_sc +
                     log_settlement_grav_sc +
                     log_chla_sc +
                     connectivity_sc +
                     mpa_status +
                     (1 | country),
                   data = total_model_data,
                   REML = FALSE)

cat("\n--- Random effect structure: country-level ---\n")
print(make_aicc_df(list(
  "No RE"         = re_null,
  "(1 | country)" = re_country
)))

# No RE: AICc = 105.43, weight = 0.807
# (1 | country): ΔAICc = 2.86, weight = 0.193
# Country clustering not supported once predictors included.
# All models fitted as lm() throughout.


# ============================================================
#  Q1 — PRESSURE METRIC SELECTION
#
#  Scientific question:
#  Does settlement gravity outperform market gravity and
#  settlement population as a proxy for SSF pressure?
#
#  Three models identical in structure — only pressure metric
#  differs. Baseline controls (rugosity + chla) held constant.
#  Connectivity included as additional control to ensure
#  pressure ranking is robust to spatial structure.
#  Both versions compared — with and without connectivity —
#  to confirm ranking is independent of this choice.
#
#  Settlement gravity selected if AICc weight substantially
#  higher than alternatives and result consistent across
#  both model versions.
# ============================================================

# ── Without connectivity control ─────────────────────────────
q1_settgrav <- lm(log_mean_biomass ~ rugosity_sc +
                    log_settlement_grav_sc +
                    log_chla_sc,
                  data = total_model_data)

q1_mktgrav  <- lm(log_mean_biomass ~ rugosity_sc +
                    log_market_gravity_sc +
                    log_chla_sc,
                  data = total_model_data)

q1_settpop  <- lm(log_mean_biomass ~ rugosity_sc +
                    log_settlement_pop_sc +
                    log_chla_sc,
                  data = total_model_data)

cat("\n--- Q1: Pressure metric selection (without connectivity) ---\n")
print(make_aicc_df(list(
  "Settlement gravity" = q1_settgrav,
  "Market gravity"     = q1_mktgrav,
  "Settlement pop."    = q1_settpop
)))

# ── With connectivity control ─────────────────────────────────
# Confirms ranking is robust to inclusion of spatial structure
q1_settgrav_conn <- lm(log_mean_biomass ~ rugosity_sc +
                         log_settlement_grav_sc +
                         log_chla_sc +
                         connectivity_sc,
                       data = total_model_data)

q1_mktgrav_conn  <- lm(log_mean_biomass ~ rugosity_sc +
                         log_market_gravity_sc +
                         log_chla_sc +
                         connectivity_sc,
                       data = total_model_data)

q1_settpop_conn  <- lm(log_mean_biomass ~ rugosity_sc +
                         log_settlement_pop_sc +
                         log_chla_sc +
                         connectivity_sc,
                       data = total_model_data)

cat("\n--- Q1: Pressure metric selection (with connectivity) ---\n")
print(make_aicc_df(list(
  "Settlement gravity" = q1_settgrav_conn,
  "Market gravity"     = q1_mktgrav_conn,
  "Settlement pop."    = q1_settpop_conn
)))

# ── Q1 results ────────────────────────────────────────────────
#
# Without connectivity control:
#   Settlement gravity: AICc = 101.36, weight = 0.878 (SELECTED)
#   Settlement pop.:    ΔAICc = 4.92,  weight = 0.075
#   Market gravity:     ΔAICc = 5.86,  weight = 0.047
#
# With connectivity control:
#   Settlement gravity: AICc = 103.77, weight = 0.881 (SELECTED)
#   Settlement pop.:    ΔAICc = 5.02,  weight = 0.071
#   Market gravity:     ΔAICc = 5.84,  weight = 0.048
#
# Settlement gravity is the best-supported pressure metric
# regardless of whether connectivity is included as a control.
# Ranking consistent and weights near-identical across both
# comparisons — selection robust to spatial structure.
# Weight increases slightly with connectivity included
# (0.878 → 0.881), confirming no confounding between
# pressure metric ranking and larval connectivity.
#
# Market gravity — validated globally by Cinner et al.
# (2016) and used by Samoilys et al. (2025) — is the
# weakest performer in this SSF-dominated system
# (ΔAICc = 5.86 and 5.84), supporting the hypothesis
# that settlement proximity better captures the spatial
# footprint of subsistence-oriented exploitation than
# market access.
#
# Settlement gravity used as primary pressure metric
# throughout Q2, Q3, and sensitivity analyses.
# Market gravity and settlement population retained for
# sensitivity analysis only.

cat("\n--- Q1: Coefficient summary — settlement gravity model ---\n")
summary(q1_settgrav)

# lm(log_mean_biomass ~ rugosity_sc + log_settlement_grav_sc +
#    log_chla_sc)
# n = 54 sites, adj. R² = 0.194, F(3,50) = 5.25, p = 0.003
#
# Rugosity:           β = +0.216, p = 0.009 **
# Settlement gravity: β = -0.251, p = 0.012 *
# Chla:               β = -0.125, p = 0.208 ns
#   Non-significant. Retained as baseline control —
#   negative direction reflects collinearity with
#   settlement gravity (r = -0.57): productive sites
#   tend to be more remote and less fished.


# ============================================================
#  Q2 — ARE CONNECTIVITY AND MPA BASELINE DRIVERS?
#
#  Baseline model (fixed a priori):
#    log(biomass) ~ rugosity + settlement_gravity + chla
#
#  Sequence:
#    Null         → intercept only
#    Baseline     → rugosity + settlement gravity + chla
#    + conn       → adds connectivity as main effect
#    + MPA        → adds MPA as main effect
#    + conn + MPA → global additive upper bound
#
#  Both ΔAICc and ΔR² reported.
#  Best-supported model from Q2 becomes reference for Q3.
# ============================================================

m_null <- lm(log_mean_biomass ~ 1,
             data = total_model_data)

# Baseline (fixed a priori)
# rugosity: habitat carrying capacity (Samoilys et al. 2025)
# settlement gravity: SSF pressure (selected in Q1)
# chla: background productivity (Warmuth et al. 2024)
m_baseline <- lm(log_mean_biomass ~ rugosity_sc +
                   log_settlement_grav_sc +
                   log_chla_sc,
                 data = total_model_data)

# Baseline + connectivity
# Tests Q2: does connectivity add as main effect?
# Extends Warmuth et al. (2024) — includes fishing pressure
# as covariate and uses biomass not abundance.
m_baseline_conn <- lm(log_mean_biomass ~ rugosity_sc +
                        log_settlement_grav_sc +
                        log_chla_sc +
                        connectivity_sc,
                      data = total_model_data)

# Baseline + MPA
# Tests Q2: does formal protection add beyond local baseline?
# Partial overlap with settlement gravity expected (r = -0.30).
m_baseline_mpa <- lm(log_mean_biomass ~ rugosity_sc +
                       log_settlement_grav_sc +
                       log_chla_sc +
                       mpa_status,
                     data = total_model_data)

# Baseline + connectivity + MPA (global additive)
# Upper bound of additive explanatory power.
m_global_additive <- lm(log_mean_biomass ~ rugosity_sc +
                          log_settlement_grav_sc +
                          log_chla_sc +
                          connectivity_sc +
                          mpa_status,
                        data = total_model_data)

q2_models <- list(
  "Null"                  = m_null,
  "Baseline"              = m_baseline,
  "Baseline + conn"       = m_baseline_conn,
  "Baseline + MPA"        = m_baseline_mpa,
  "Baseline + conn + MPA" = m_global_additive
)

cat("\n--- Q2: Hierarchical model comparison (AICc ranked) ---\n")
print(make_aicc_df(q2_models))

cat("\n--- Q2: Variance explained relative to baseline ---\n")
baseline_r2 <- summary(m_baseline)$adj.r.squared

q2_models %>%
  imap_dfr(~ tibble(
    Model  = .y,
    Adj_R2 = round(summary(.x)$adj.r.squared, 3)
  )) %>%
  mutate(
    Delta_R2 = round(Adj_R2 - baseline_r2, 3),
    Delta_R2 = ifelse(Model %in% c("Null", "Baseline"),
                      NA, Delta_R2)
  ) %>%
  print()

cat("\n--- Q2: Baseline model coefficients ---\n")
summary(m_baseline)

cat("\n--- Q2: Baseline + connectivity coefficients ---\n")
summary(m_baseline_conn)

cat("\n--- Q2: Baseline + MPA coefficients ---\n")
summary(m_baseline_mpa)

# ── Q2 results ────────────────────────────────────────────────
#
# AICc comparison:
#   Baseline:              AICc = 101.36, weight = 0.605 (BEST)
#   Baseline + conn:       ΔAICc = 2.41,  weight = 0.181
#   Baseline + MPA:        ΔAICc = 3.01,  weight = 0.134
#   Baseline + conn + MPA: ΔAICc = 4.39,  weight = 0.067
#   Null:                  ΔAICc = 7.78,  weight = 0.012
#
# Variance explained (adj. R² relative to baseline = 0.194):
#   Baseline + conn:       ΔR² = -0.015
#   Baseline + MPA:        ΔR² =  0.000
#   Baseline + conn + MPA: ΔR² = +0.003
#
# Connectivity (β = +0.028, p = 0.739, ΔAICc = 2.41):
#   No independent effect on biomass once human pressure
#   and ecological context are controlled. Directly extends
#   Warmuth et al. (2024): regional connectivity signal on
#   herbivore abundance does not translate to total biomass
#   once fishing pressure is included as a covariate.
#
# MPA status (low: β = -0.342, p = 0.171;
#             medium: β = -0.084, p = 0.669; ΔAICc = 3.01):
#   Neither protection level significantly improves biomass.
#   Negative coefficients reflect preferential MPA placement
#   in lower-pressure areas (r = -0.30 with settlement
#   gravity). MPA adds zero adj. R² beyond baseline.
#
# Rugosity (β = +0.216, p = 0.009) and settlement gravity
# (β = -0.251, p = 0.012) stable across all Q2 models.
#
# Best-supported model: Baseline (weight = 0.605).
# Baseline used as reference model in Q3.


# ============================================================
#  SPATIAL AUTOCORRELATION CHECK
#  Residuals from best Q2 model (baseline) tested for
#  spatial structure. Reported as diagnostic.
# ============================================================

site_coords <- location_data %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(lon = first(longitude),
            lat = first(latitude),
            .groups = "drop")

total_model_data_coords <- total_model_data %>%
  left_join(site_coords, by = "site")

coords_mat <- cbind(total_model_data_coords$lon,
                    total_model_data_coords$lat)
listw5 <- nb2listw(knn2nb(knearneigh(coords_mat, k = 5)),
                   style = "W")
cat("\n--- Spatial autocorrelation: best Q2 model residuals ---\n")
print(moran.test(residuals(m_baseline), listw5))

# Moran's I = 0.140, p = 0.015 — weak but significant.
#
# Spatial error modelling not pursued for two reasons:
# (1) Country-level random effects tested and not
#     supported (ΔAICc = 2.86 vs no RE — see random
#     effect structure section). Between-country
#     clustering is already absorbed by the fixed
#     predictors, meaning the four geographically
#     separate clusters do not represent unmodelled
#     spatial structure requiring correction.
# (2) Discontinuous sampling design — k-NN spatial
#     weights bridge across isolated country clusters,
#     producing a weights matrix that does not reflect
#     true within-region spatial covariance structure
#     (Kissling & Carl 2008, Dormann et al. 2007).
#
# Weak residual autocorrelation acknowledged as a
# limitation. May slightly inflate type I error rates
# for pressure and rugosity coefficients.

# ============================================================
#  Q3 — DO MANAGEMENT AND CONNECTIVITY MODIFY
#       PRESSURE EFFECTS?
#
#  Three a priori interaction hypotheses:
#
#  H1 — MPA × pressure:
#       Protection only detectable where external pressure
#       is low enough for recovery.
#       (Cinner et al. 2016)
#
#  H2 — Connectivity × pressure:
#       Well-connected sites sustain higher biomass under
#       pressure through larval replenishment.
#       (Warmuth et al. 2024)
#
#  H3 — MPA × connectivity:
#       Protected sites recover faster where connectivity
#       is high enough to subsidise recruitment.
#
#  Reference model: best Q2 model (baseline).
#  Interaction preferred if ΔAICc > 2 vs reference AND
#  pattern consistent with a priori hypothesis.
#  All three run regardless of Q2 main-effect results.
#  Interpretation conditional on Q2 support.
# ============================================================

# Baseline confirmed as best Q2 model (weight = 0.605)
q3_reference <- m_baseline

# H1: MPA effectiveness depends on external fishing intensity
m_int_mpa_press <- lm(log_mean_biomass ~ rugosity_sc +
                        log_chla_sc +
                        mpa_status * log_settlement_grav_sc,
                      data = total_model_data)

# H2: Connectivity buffers exploitation effects
m_int_conn_press <- lm(log_mean_biomass ~ rugosity_sc +
                         log_chla_sc +
                         connectivity_sc * log_settlement_grav_sc,
                       data = total_model_data)

# H3: MPA effectiveness depends on larval supply
m_int_mpa_conn <- lm(log_mean_biomass ~ rugosity_sc +
                       log_chla_sc +
                       log_settlement_grav_sc +
                       mpa_status * connectivity_sc,
                     data = total_model_data)

q3_models <- list(
  "Reference (best Q2)"    = q3_reference,
  "H1: MPA × pressure"     = m_int_mpa_press,
  "H2: Conn × pressure"    = m_int_conn_press,
  "H3: MPA × connectivity" = m_int_mpa_conn
)

cat("\n--- Q3: Interaction model comparison ---\n")
print(make_aicc_df(q3_models))

cat("\n--- Q3: H1 MPA × pressure coefficients ---\n")
summary(m_int_mpa_press)

cat("\n--- Q3: H2 Connectivity × pressure coefficients ---\n")
summary(m_int_conn_press)

cat("\n--- Q3: H3 MPA × connectivity coefficients ---\n")
summary(m_int_mpa_conn)

# ── Q3 results ────────────────────────────────────────────────
#
# AICc comparison:
#   H1: MPA × pressure:     AICc = 99.40,  weight = 0.587 (BEST)
#   Reference (baseline):   ΔAICc = 1.95,  weight = 0.221
#   H2: Conn × pressure:    ΔAICc = 3.01,  weight = 0.130
#   H3: MPA × connectivity: ΔAICc = 4.50,  weight = 0.062
#
# H1: MPA × pressure (ΔAICc = 1.95 vs reference)
#   Statistically competitive but not decisively better
#   than baseline. Medium MPA × pressure significant
#   (β = +1.071, p = 0.003). See exploratory tests below.
#
# H2: Connectivity × pressure (ΔAICc = 3.01)
#   Not supported. Interaction marginal (β = -0.225,
#   p = 0.060) and in wrong direction — negative sign
#   suggests connectivity amplifies rather than buffers
#   pressure effects, contrary to hypothesis. Pressure
#   main effect weakens (β = -0.176, p = 0.095).
#
# H3: MPA × connectivity (ΔAICc = 4.50)
#   Not supported. Medium MPA × connectivity significant
#   (β = +0.443, p = 0.021) but low MPA × connectivity
#   not (β = -1.165, p = 0.476) — inconsistent pattern,
#   not interpretable as a coherent biological mechanism.
#
# Rugosity stable across all Q3 models (β = 0.21–0.26,
# p < 0.01) — local habitat signal robust throughout.

# ── H1 exploratory tests ──────────────────────────────────────
# Two targeted tests assess whether the significant medium
# MPA × pressure interaction reflects a genuine ecological
# pattern or a data structure artefact.

# ── VIF check on H1 interaction model ────────────────────────
cat("\n--- H1: VIF check ---\n")
vif(m_int_mpa_press)
# All GVIF^(1/(2*Df)) < 1.5 — no multicollinearity concern.
# Interaction instability is not a collinearity artefact.

# ── Cook's distance on H1 interaction model ───────────────────
par(mfrow = c(1, 2))
plot(m_int_mpa_press, which = c(4, 5),
     main = "H1: Cook's distance and leverage")
par(mfrow = c(1, 1))

# Identify the high leverage site
cat("\n--- H1: Top 5 leverage sites (post-correction) ---\n")
hatvalues(m_int_mpa_press) %>%
  sort(decreasing = TRUE) %>%
  head(5) %>%
  print()

# Check which site this is
high_lev_idx <- which.max(hatvalues(m_int_mpa_press))
total_model_data[high_lev_idx, ] %>%
  dplyr::select(site, mpa_status,
                log_settlement_grav_sc,
                log_mean_biomass)

# Leverage diagnostics on H1 interaction model:
#   Site 28 (mirereni): leverage = 0.880 — highest in model.
#     Low MPA site at extreme high pressure (z = 2.03),
#     the only low MPA observation at the far right of the
#     pressure gradient. Anchors the low MPA slope at high
#     pressure. Standardised residual near zero — influential
#     on slope estimates but well-fitted by the model.
#     Low MPA × pressure interaction not significant
#     (β = +0.235, p = 0.451) so mirereni is not driving
#     a spurious result.
#   Sites 32, 9, 47, 13: leverage = 0.29–0.37, acceptable.
#   No sites exceed Cook's distance threshold of 0.5.
#   Overall influence structure acceptable for n = 54.

# Check medium MPA site distribution on pressure gradient
cat("\n--- H1: Medium MPA sites on pressure gradient ---\n")
total_model_data %>%
  filter(mpa_status == "medium") %>%
  dplyr::select(site, log_settlement_grav_sc,
                log_mean_biomass) %>%
  arrange(log_settlement_grav_sc) %>%
  print()

# Test 1 — raw data: pressure-biomass by MPA category
ggplot(total_model_data,
       aes(x = log_settlement_grav_sc,
           y = log_mean_biomass,
           colour = mpa_status,
           shape  = mpa_status)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  scale_colour_manual(
    values = c("none"   = "#636363",
               "low"    = "#74a9cf",
               "medium" = "#0570b0"),
    labels = c("No MPA", "Low protection", "Medium protection")
  ) +
  scale_shape_manual(
    values = c("none" = 16, "low" = 17, "medium" = 15),
    labels = c("No MPA", "Low protection", "Medium protection")
  ) +
  labs(x      = "log(Settlement gravity) (standardised)",
       y      = "log(Mean biomass)",
       colour = "MPA status",
       shape  = "MPA status") +
  theme_bw(base_size = 12) +
  theme(axis.title        = element_text(face = "bold"),
        legend.position   = c(0.82, 0.85),
        legend.background = element_rect(fill      = "white",
                                         colour    = "grey80",
                                         linewidth = 0.3),
        panel.grid.minor  = element_blank())

# Test 2 — biomass comparison within observed pressure range
cat("\n--- H1: Biomass by MPA within observed range (z < 0.20) ---\n")
total_model_data %>%
  filter(log_settlement_grav_sc < 0.20) %>%
  group_by(mpa_status) %>%
  summarise(
    n            = n(),
    mean_biomass = round(mean(log_mean_biomass), 3),
    sd_biomass   = round(sd(log_mean_biomass),   3),
    min_pressure = round(min(log_settlement_grav_sc), 3),
    max_pressure = round(max(log_settlement_grav_sc), 3),
    .groups = "drop"
  ) %>%
  print()

# ── H1 exploratory results ────────────────────────────────────
#
# Test 1 — raw data plot:
#   Medium MPA sites (squares) show high scatter within
#   a narrow low-pressure range (z = -0.95 to 0.17).
#   Positive fitted slope chases noise — no genuine
#   positive trend visible in the raw data points.
#   No MPA sites span the full pressure range with a
#   clean negative slope — this is the genuine signal.
#
# Test 2 — biomass at equivalent pressure (z < 0.20):
#   No MPA:  9.92 ± 0.620 (n = 18)
#   Low MPA: 8.96 ± 0.068 (n = 3)
#   Medium:  9.77 ± 0.729 (n = 17)
#   Medium MPA sites show no biomass advantage over
#   unprotected sites at equivalent pressure values
#   (mean difference = -0.15). All 17 medium MPA sites
#   cluster between z = -0.95 and z = 0.17 — no medium
#   MPA observations exist at high pressure. The positive
#   slope at high pressure in the interaction model is
#   therefore extrapolation beyond the observed data range,
#   not an empirical pattern.
#   High within-group variance at medium MPA sites
#   (SD = 0.729 vs 0.620 for none) likely reflects
#   heterogeneity in enforcement quality across sites
#   sharing the same nominal protection category —
#   consistent with Cinner et al. (2016) showing
#   compliance rather than MPA status determines outcomes.
#
# ── Overall Q3 conclusion ─────────────────────────────────────
# No interaction model provides a convincing, ecologically
# interpretable improvement over the baseline. H1
# statistically competitive (ΔAICc = 1.95) but exploratory
# tests confirm a data structure artefact — MPA placement
# confounded with pressure gradient, no biomass advantage
# at medium MPA sites within observed data range. H2 and
# H3 not supported. Reef fish biomass is not detectably
# modified by connectivity or formal protection at this
# spatial scale.


# ============================================================
#  SENSITIVITY ANALYSIS
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
# Mirrors baseline structure — only pressure metric swapped.
# Purpose: robustness check, not model selection.
# AICc comparison not reported here — see Q1 for metric
# selection. 
# Coefficients confirm conclusions are not metric-dependent.

sens_mktgrav <- lm(log_mean_biomass ~ rugosity_sc +
                     log_market_gravity_sc +
                     log_chla_sc,
                   data = total_model_data)

sens_settpop  <- lm(log_mean_biomass ~ rugosity_sc +
                      log_settlement_pop_sc +
                      log_chla_sc,
                    data = total_model_data)

cat("\n--- Sensitivity (a): alternative pressure metrics ---\n")
cat("Market gravity:\n")
print(summary(sens_mktgrav)$coefficients)
cat("\nSettlement population:\n")
print(summary(sens_settpop)$coefficients)

# ── Sensitivity (a) results ───────────────────────────────────
#
# Alternative pressure metrics — baseline structure retained,
# only pressure metric substituted. Purpose: confirm Q1 and
# Q2 conclusions are not metric-dependent.
#
# Market gravity (β = -0.086, p = 0.346):
#   Not significant. Direction consistent with settlement
#   gravity but effect size approximately one third as
#   large (β = -0.086 vs -0.251). Rugosity remains
#   significant (β = +0.220, p = 0.013). Chla near-zero
#   and non-significant (β = -0.013, p = 0.884) —
#   weaker than in primary model, consistent with market
#   gravity capturing less of the pressure-productivity
#   geographic confound than settlement gravity.
#
# Settlement population (β = -0.116, p = 0.186):
#   Not significant. Direction consistent but effect
#   size less than half that of settlement gravity.
#   Rugosity remains significant (β = +0.217, p = 0.013).
#   Chla near-zero and non-significant (β = -0.013,
#   p = 0.880).
#
# Both alternative metrics fail to reach significance,
# confirming settlement gravity as the strongest pressure
# proxy in this SSF-dominated system. The human pressure
# signal is consistent in direction across all three
# metrics but only detectable with settlement gravity —
# supporting the hypothesis that residential proximity
# better captures exploitation intensity than market
# access or aggregate population size in this context.
# Rugosity significant and stable across all three
# metrics (β = 0.217–0.220, p < 0.05) — habitat
# complexity signal robust to pressure metric choice.

# ── (b) Transect-level replication ───────────────────────────
# Replicates Q2 sequence at transect level using lmer —
# conventional choice for Gaussian mixed models with simple
# random intercept structure.
# (1 | site) accounts for non-independence of transects
# within sites.
# Models fitted with REML = TRUE for final inference.
# Refitted with REML = FALSE for AICc comparison —
# REML likelihoods not valid for comparing models with
# different fixed effect structures.

sens_t_null <- lmer(log_transect_biomass ~ 1 +
                      (1 | site),
                    data = transect_model_data,
                    REML = TRUE)

sens_t_baseline <- lmer(log_transect_biomass ~ rugosity_sc +
                          log_settlement_grav_sc +
                          log_chla_sc +
                          (1 | site),
                        data = transect_model_data,
                        REML = TRUE)

sens_t_conn <- lmer(log_transect_biomass ~ rugosity_sc +
                      log_settlement_grav_sc +
                      log_chla_sc +
                      connectivity_sc +
                      (1 | site),
                    data = transect_model_data,
                    REML = TRUE)

sens_t_mpa <- lmer(log_transect_biomass ~ rugosity_sc +
                     log_settlement_grav_sc +
                     log_chla_sc +
                     mpa_status +
                     (1 | site),
                   data = transect_model_data,
                   REML = TRUE)

sens_t_global <- lmer(log_transect_biomass ~ rugosity_sc +
                        log_settlement_grav_sc +
                        log_chla_sc +
                        connectivity_sc +
                        mpa_status +
                        (1 | site),
                      data = transect_model_data,
                      REML = TRUE)

# ── Diagnostics on best transect model ────────────────────────
# Run on REML-fitted baseline before AICc comparison.
# Standard residual diagnostics sufficient for Gaussian lmer.

par(mfrow = c(1, 2))
plot(fitted(sens_t_baseline), residuals(sens_t_baseline),
     xlab = "Fitted values", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2, col = "grey60")
lines(lowess(fitted(sens_t_baseline),
             residuals(sens_t_baseline)), col = "red")
qqnorm(residuals(sens_t_baseline),
       main = "Q-Q Residuals")
qqline(residuals(sens_t_baseline), col = "red")
par(mfrow = c(1, 1))

# ── AICc comparison — ML refits ───────────────────────────────
sens_t_null_ml     <- update(sens_t_null,     REML = FALSE)
sens_t_baseline_ml <- update(sens_t_baseline, REML = FALSE)
sens_t_conn_ml     <- update(sens_t_conn,     REML = FALSE)
sens_t_mpa_ml      <- update(sens_t_mpa,      REML = FALSE)
sens_t_global_ml   <- update(sens_t_global,   REML = FALSE)

cat("\n--- Sensitivity (b): transect-level model comparison ---\n")
print(make_aicc_df(list(
  "Null"                  = sens_t_null_ml,
  "Baseline"              = sens_t_baseline_ml,
  "Baseline + conn"       = sens_t_conn_ml,
  "Baseline + MPA"        = sens_t_mpa_ml,
  "Baseline + conn + MPA" = sens_t_global_ml
)))

# ── Coefficient summary — REML-fitted baseline ────────────────
cat("\n--- Sensitivity (b): baseline coefficients (REML) ---\n")
summary(sens_t_baseline)

# ── Transect model diagnostics ────────────────────────────────
# Residuals vs Fitted: flat with minor upward trend at high
#   fitted values — acceptable for n = 243 transects.
# Q-Q: excellent — points follow theoretical line closely
#   across full range, minor deviations at lower tail
#   (2-3 points) within acceptable range.
# Gaussian lmer structure confirmed appropriate at
# transect level.

# ── Sensitivity (b) results ───────────────────────────────────
#
# AICc comparison (ML):
#   Baseline:              AICc = 634.93, weight = 0.601 (BEST)
#   Baseline + conn:       ΔAICc = 2.10,  weight = 0.210
#   Baseline + MPA:        ΔAICc = 3.08,  weight = 0.129
#   Baseline + conn + MPA: ΔAICc = 4.64,  weight = 0.059
#   Null:                  ΔAICc = 13.05, weight = 0.001
#
# Model ordering identical to site-level Q2 — baseline
# best-supported at both levels of analysis.
#
# Baseline coefficients (REML, n = 243 transects,
# 54 sites):
#   Rugosity:           β = +0.240, t = 3.100 **
#   Settlement gravity: β = -0.285, t = -3.006 **
#   Chla:               β = -0.205, t = -2.187 *
#
# Random effects:
#   Site-level variance:    0.168 (SD = 0.410)
#   Residual variance:      0.651 (SD = 0.807)
#   ICC = 0.168 / (0.168 + 0.651) = 0.205 — approximately
#   20% of total variance attributable to between-site
#   differences beyond what fixed predictors explain.
#   Confirms (1 | site) random intercept is justified.
#
# Notable difference from site-level: chla reaches
# significance at transect level (β = -0.205, t = -2.187)
# but not at site level (β = -0.125, p = 0.208). Greater
# statistical power at transect level (n = 243 vs 54)
# allows detection of a weaker signal. Direction
# consistent across both levels.
#
# Rugosity and settlement gravity significant and stable
# in direction and magnitude at both levels — primary
# conclusions robust to aggregation.
#
# Connectivity and MPA not supported at transect level
# (ΔAICc = 2.10 and 3.08 respectively) — consistent
# with site-level Q2. Site-level aggregation does not
# alter qualitative inference.

# ============================================================
#  MARGINAL EFFECT PLOTS
#  All predictors in baseline model.
#  Observed data overlaid on fitted lines.
#  Chla non-significant at site level (p = 0.208) but
#  retained as productivity control and reaches significance
#  at transect level (β = -0.205, t = -2.187).
# ============================================================

best_model <- m_baseline

p_rugosity <- plot_effect(best_model, total_model_data,
                          "rugosity_sc",
                          "Rugosity (standardised)")

p_pressure <- plot_effect(best_model, total_model_data,
                          "log_settlement_grav_sc",
                          "log(Settlement gravity) (standardised)")

p_chla     <- plot_effect(best_model, total_model_data,
                          "log_chla_sc",
                          "log(Chlorophyll-a) (standardised)")

gridExtra::grid.arrange(p_rugosity, p_pressure, p_chla,
                        ncol = 3)

# jpeg("marginal_effects_baseline.jpg",
#      width = 33, height = 11, units = "cm", res = 300)
# gridExtra::grid.arrange(p_rugosity, p_pressure, p_chla,
#                         ncol = 3)
# dev.off()

# ── MPA marginal means ────────────────────────────────────────
# From m_baseline_mpa — descriptive only.
# MPA not supported in Q2 (ΔAICc = 3.01).
# Shows direction of effect; not a supported relationship.

mpa_grid <- data.frame(
  mpa_status             = factor(c("none", "low", "medium"),
                                  levels = c("none", "low",
                                             "medium"),
                                  labels = c("None", "Low",
                                             "Medium")),
  rugosity_sc            = 0,
  log_settlement_grav_sc = 0,
  log_chla_sc            = 0
)

mpa_pred     <- predict(m_baseline_mpa,
                        newdata = mpa_grid, se.fit = TRUE)
mpa_grid$fit <- mpa_pred$fit
mpa_grid$lwr <- mpa_pred$fit - 1.96 * mpa_pred$se.fit
mpa_grid$upr <- mpa_pred$fit + 1.96 * mpa_pred$se.fit

p_mpa <- ggplot(mpa_grid, aes(x = mpa_status, y = fit)) +
  geom_hline(yintercept = mpa_grid$fit[1],
             linetype   = "dashed",
             colour     = "grey70",
             linewidth  = 0.4) +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  colour    = "#0570b0",
                  linewidth = 0.7,
                  size      = 0.6) +
  labs(x       = "MPA status",
       y       = "Fitted log(biomass)",
       caption = paste("From baseline + MPA model —",
                       "MPA not supported in Q2 (ΔAICc = 3.01).",
                       "\nDescriptive only; dashed line =",
                       "no MPA reference.")) +
  theme_bw(base_size = 12) +
  theme(
    axis.title         = element_text(face = "bold"),
    plot.caption       = element_text(colour = "grey50", size = 8),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

p_mpa

# jpeg("marginal_effect_mpa.jpg",
#      width = 12, height = 12, units = "cm", res = 300)
# print(p_mpa)
# dev.off()

# ============================================================
#  RESULTS SUMMARY
#  Quick reference for writing — not part of analysis.
#  Verify all values match reported results before writing.
# ============================================================

results_summary <- tribble(
  ~Question,  ~Best_model,      ~Key_finding,
  "Q1",       "Sett. gravity",  "weight = 0.878; ΔAICc = 4.92 vs next best",
  "Q2 conn",  "Baseline",       "β = +0.028, p = 0.739, ΔAICc = 2.41",
  "Q2 MPA",   "Baseline",       "ns, ΔAICc = 3.01, ΔR² = 0.000",
  "Q3 H1",    "H1 (marginal)",  "ΔAICc = 1.95, data structure artefact",
  "Q3 H2",    "Reference",      "ΔAICc = 3.01, wrong direction",
  "Q3 H3",    "Reference",      "ΔAICc = 4.50, inconsistent pattern"
)

cat("\n--- Results summary ---\n")
print(results_summary)

# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()