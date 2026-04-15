# ============================================================
#  DRIVERS OF REEF FISH BIOMASS
#  Chapter 1 — Site-level Analysis
#
#  Analytical framework (four sequential stages):
#
#  STAGE 1 — Variance partitioning
#             Quantifies unique and shared variance attributable
#             to local ecological, spatial, and environmental
#             process groups. MPA status excluded (governance
#             variable — tested separately in Stage 2).
#
#  STAGE 2 — Hierarchical model comparison
#             Nested model sequence adds process groups
#             progressively. Tests incremental explanatory value
#             of each group beyond the local baseline.
#             Reports both ΔAICc and ΔR² at each step.
#
#  STAGE 3 — Interaction testing (conditional on Stage 2)
#             Three a priori interactions test whether spatial
#             management and connectivity modify local-scale
#             relationships. Run only if Global outperforms
#             Local baseline by ΔAICc > 2.
#
#  STAGE 4 — Sensitivity analysis
#             (a) Alternative pressure metrics — best Stage 2
#                 model refitted with settlement pop. and market
#                 gravity substituted for settlement gravity.
#             (b) Transect-level mixed model — confirms
#                 site-level findings not an artefact of
#                 aggregation.
#
#  Study design:
#    Transects nested within stations, stations within sites,
#    sites within locations, locations within countries.
#
#  Response:
#    log(mean total fish biomass per site)  [site-level LM]
#    log(transect total biomass)            [sensitivity, GLMM]
#
#  Primary pressure metric: settlement gravity (selected via
#    pre-analysis univariate AICc comparison — see Section 2.3)
# ============================================================

options(scipen = 999)

# ── PACKAGES ─────────────────────────────────────────────────
library(tidyverse)
library(sf)
library(glmmTMB)
library(DHARMa)
library(MuMIn)
library(AICcmodavg)
library(vegan)       # varpart()
library(ggcorrplot)
library(corrplot)
library(gridExtra)
library(MASS)
library(here)

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


# ============================================================
#  DATA LOADING
# ============================================================

fish_data      <- readr::read_rds(here("processed_data", "clean_fish_connectivity.rds"))
gravity_data   <- readr::read_rds(here("city_data", "locations_with_grav_combined.rds"))
chla_data      <- read.csv(here("processed_data", "locations_with_chla_2009.csv"))
rugosity_data  <- readr::read_rds(here("processed_data", "clean_dive_details_connectivity.rds"))
location_data  <- readr::read_rds(here("processed_data", "clean_location_connectivity.rds"))
dhw_data       <- readr::read_rds(here("processed_data", "locations_with_dhw_2009.rds"))


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
gravity_sites <- gravity_data %>%
  st_drop_geometry() %>%
  group_by(site) %>%
  summarise(
    market_gravity  = mean(market_grav,          na.rm = TRUE),
    settlement_pop  = mean(settlement_tot_pop,   na.rm = TRUE),
    settlement_grav = mean(nearest_pop75_grav,   na.rm = TRUE),
    .groups = "drop"
  )

# ── Chlorophyll-a ─────────────────────────────────────────────
chla_sites <- chla_data %>%
  group_by(site) %>%
  summarise(mean_annual_chla = mean(chla_annual_mean, na.rm = TRUE),
            .groups = "drop")

# ── Degree heating weeks ──────────────────────────────────────
dhw_sites <- dhw_data %>%
  filter(!is.na(max_dhw)) %>%
  group_by(site) %>%
  summarise(mean_annual_dhw = mean(max_dhw, na.rm = TRUE),
            .groups = "drop")

# ── Rugosity ──────────────────────────────────────────────────
rugosity_sites <- rugosity_data %>%
  group_by(site) %>%
  summarise(rugosity = mean(rugosity, na.rm = TRUE),
            .groups = "drop")

# ── MPA status and connectivity ───────────────────────────────
location_sites <- location_data %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(
    mpa_status   = first(mpa_status),
    connectivity = mean(prop_connectivity, na.rm = TRUE),
    .groups = "drop"
  )


# ============================================================
#  TRANSFORMATIONS, SCALING, AND PREDICTOR CHECKS
# ============================================================

raw_predictors <- location_sites %>%
  left_join(chla_sites,     by = "site") %>%
  left_join(dhw_sites,      by = "site") %>%
  left_join(rugosity_sites, by = "site") %>%
  left_join(gravity_sites,  by = "site")

# ── Apply transformations ─────────────────────────────────────
# Rugosity:         no transformation (approximately normal)
# Gravity metrics:  log (right-skewed)
# Chla:             log (right-skewed)
# DHW:              log(x + 1) (right-skewed with zeros)
# MPA status:       ordered factor (governance gradient)
# Connectivity:     no transformation (inspect distribution)

transformed_predictors <- raw_predictors %>%
  transmute(
    site                = site,
    rugosity            = rugosity,
    log_settlement_grav = log(settlement_grav),
    log_settlement_pop  = log(settlement_pop),
    log_market_gravity  = log(market_gravity),
    log_chla            = log(mean_annual_chla),
    log_max_dhw         = log(mean_annual_dhw + 1),
    connectivity        = connectivity,
    mpa_status          = factor(mpa_status,
                                 levels  = c("none", "low", "medium"),
                                 ordered = FALSE)   # ← change here
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
    connectivity_sc        = as.numeric(scale(connectivity)),
    log_chla_sc            = as.numeric(scale(log_chla)),
    log_max_dhw_sc         = as.numeric(scale(log_max_dhw)),
    mpa_status             = mpa_status   # ordered factor — not scaled
  )

# ── Predictor correlation matrix ──────────────────────────────
# Check for blocking collinearity before modelling.
# Rule of thumb: |r| > 0.70 warrants caution; > 0.80 is problematic.

corr_matrix <- scaled_predictors %>%
  mutate(mpa_numeric = as.numeric(mpa_status)) %>%
  dplyr::select(ends_with("_sc"), mpa_numeric) %>%
  rename(
    "Market gravity"     = log_market_gravity_sc,
    "Settlement gravity" = log_settlement_grav_sc,
    "Settlement pop."    = log_settlement_pop_sc,
    "Chlorophyll-a"      = log_chla_sc,
    "Max DHW"            = log_max_dhw_sc,
    "Rugosity"           = rugosity_sc,
    "Connectivity"       = connectivity_sc,
    "MPA status"         = mpa_numeric
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

# ── Collinearity summary ──────────────────────────────────────
# Gravity metrics (settlement grav / pop / market): r = 0.53–0.54
#   → Use one per model only. Settlement gravity selected as
#     primary metric (see pre-analysis selection below).
# Chla vs settlement gravity: r = 0.57 — monitor.
# MPA vs chla: r = 0.42 — acceptable.
# All other pairs: |r| < 0.40 — no concerns.


# ============================================================
#  PRE-ANALYSIS: HUMAN PRESSURE METRIC SELECTION
#
#  Three gravity metrics all proxy local human pressure.
#  Select the best-performing metric via univariate AICc
#  before entering the main candidate set. This is a
#  pre-analysis selection step, not part of main inference.
# ============================================================

# ── Build site-level dataset for selection ───────────────────
pressure_selection_data <- total_transects %>%
  left_join(scaled_predictors, by = "site") %>%
  group_by(site) %>%
  summarise(
    log_mean_biomass       = log(mean(transect_total_biomass, na.rm = TRUE)),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_settlement_pop_sc  = first(log_settlement_pop_sc),
    log_market_gravity_sc  = first(log_market_gravity_sc),
    .groups = "drop"
  )

# ── Univariate models ────────────────────────────────────────
press_settgrav <- lm(log_mean_biomass ~ log_settlement_grav_sc,
                     data = pressure_selection_data)
press_settpop  <- lm(log_mean_biomass ~ log_settlement_pop_sc,
                     data = pressure_selection_data)
press_mktgrav  <- lm(log_mean_biomass ~ log_market_gravity_sc,
                     data = pressure_selection_data)

cat("\n--- Pre-analysis: human pressure metric selection ---\n")
print(make_aicc_df(list(
  "Settlement gravity" = press_settgrav,
  "Settlement pop."    = press_settpop,
  "Market gravity"     = press_mktgrav
)))

# ── Decision ─────────────────────────────────────────────────
# Settlement gravity: AICc = 105.87, weight = 0.67 (SELECTED)
# Settlement pop.:    ΔAICc = 2.41,  weight = 0.20
# Market gravity:     ΔAICc = 3.34,  weight = 0.12
#
# Settlement gravity is the primary pressure metric throughout.
# Settlement pop. and market gravity retained for sensitivity
# analysis only (Stage 4).

rm(pressure_selection_data)


# ============================================================
#  ANALYSIS DATASETS
# ============================================================

# ── Combine all predictors ────────────────────────────────────
final_predictors <- scaled_predictors %>%
  dplyr::select(
    site,
    rugosity_sc,                # LOCAL      — habitat complexity
    log_settlement_grav_sc,     # LOCAL      — human pressure (primary)
    connectivity_sc,            # SPATIAL    — larval network position
    mpa_status,                 # GOVERNANCE — protection gradient
    log_chla_sc,                # ENV        — primary productivity
    log_max_dhw_sc,             # ENV        — thermal stress history
    log_settlement_pop_sc,      # SENSITIVITY — alternative pressure
    log_market_gravity_sc       # SENSITIVITY — alternative pressure
  )

# ── Transect-level dataset (sensitivity check) ───────────────
transect_model_data <- total_transects %>%
  left_join(final_predictors, by = "site") %>%
  mutate(log_transect_biomass = log(transect_total_biomass))

cat("\nTransect model data:", nrow(transect_model_data), "rows,",
    n_distinct(transect_model_data$site), "sites\n")
cat("Biomass zeros:", sum(transect_model_data$transect_total_biomass == 0), "\n")

# ── Site-level dataset (primary analysis) ────────────────────
total_model_data <- transect_model_data %>%
  group_by(site, country) %>%
  summarise(
    log_mean_biomass       = log(mean(transect_total_biomass, na.rm = TRUE)),
    mean_biomass           = mean(transect_total_biomass, na.rm = TRUE),
    n_transects            = n(),
    rugosity_sc            = first(rugosity_sc),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    connectivity_sc        = first(connectivity_sc),
    mpa_status             = first(mpa_status),
    log_chla_sc            = first(log_chla_sc),
    log_max_dhw_sc         = first(log_max_dhw_sc),
    log_settlement_pop_sc  = first(log_settlement_pop_sc),
    log_market_gravity_sc  = first(log_market_gravity_sc),
    .groups = "drop"
  ) %>%
  mutate(
    site    = as.factor(site),
    country = as.factor(country)
  )

cat("\nSite model data:", nrow(total_model_data), "sites,",
    n_distinct(total_model_data$country), "countries\n")

# ── Verify no NAs in primary predictors ──────────────────────
total_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                connectivity_sc, mpa_status,
                log_chla_sc, log_max_dhw_sc) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)


# ============================================================
#  MODEL FAMILY SELECTION (documented, not re-run each time)
#
#  Three candidate families tested on an anchor model
#  (rugosity + settlement gravity + chla + DHW):
#
#  Gaussian (raw):   rejected — heteroscedasticity, non-normality,
#                    curved residual pattern.
#  Gamma (log link): rejected — Q-Q shows systematic deviation
#                    from half-normal line.
#  Gaussian (log):   SELECTED — flat residuals vs fitted,
#                    Q-Q closely follows theoretical line,
#                    homoscedasticity met. Minor upper tail
#                    deviation at 2 sites, within acceptable range.
#
#  AICc (raw scale only — log not comparable):
#    Gamma:         AICc = 1146.93, weight = 1.00
#    Gaussian raw:  AICc = 1191.14, ΔAICc = 44.21
#  Gamma outperforms raw Gaussian, confirming transformation
#  is required. Gaussian log selected on diagnostic quality.
#
#  Proceed: lm() on log_mean_biomass throughout.
# ============================================================


# ============================================================
#  RANDOM EFFECT STRUCTURE (documented, not re-run each time)
#
#  Country-level RE tested using glmmTMB (Gaussian):
#    No RE:          AICc weight = 0.78
#    (1 | country):  AICc weight = 0.22, ΔAICc = 2.53
#
#  Country-level clustering not supported once environmental
#  and human pressure predictors are included.
#  All site-level models fitted as lm() — no random effects.
# ============================================================


# ============================================================
#  STAGE 1 — VARIANCE PARTITIONING
#
#  Quantifies unique and shared variance attributable to three
#  a priori ecological process groups.
#
#  Groups defined before analysis based on mechanistic role:
#    Local:       rugosity + settlement gravity
#    Spatial:     connectivity
#    Environment: chla + DHW
#
#  MPA status excluded — governance variable, not an ecological
#  filter. Tested independently in Stage 2.
#
#  Uses vegan::varpart() with adjusted R² throughout.
# ============================================================

# ── Build matrix inputs for varpart() ────────────────────────
# varpart() requires predictor matrices, not a formula.
# MPA status is excluded here (see rationale above).

vp_local <- total_model_data %>%
  dplyr::select(rugosity_sc, log_settlement_grav_sc) %>%
  as.data.frame()

vp_spatial <- total_model_data %>%
  dplyr::select(connectivity_sc) %>%
  as.data.frame()

vp_environ <- total_model_data %>%
  dplyr::select(log_chla_sc, log_max_dhw_sc) %>%
  as.data.frame()

y_biomass <- total_model_data$log_mean_biomass

# ── Run variance partition ────────────────────────────────────
vp_result <- varpart(y_biomass,
                     vp_local,
                     vp_spatial,
                     vp_environ)

cat("\n--- Stage 1: Variance partitioning ---\n")
print(vp_result)

# ── Venn diagram ──────────────────────────────────────────────
jpeg("varpart_venn.jpg", width = 18, height = 16,
     units = "cm", res = 300)
plot(vp_result,
     Xnames = c("Local", "Spatial", "Environment"),
     bg     = c("#4dac26", "#2c7bb6", "#d7191c"),
     alpha  = 80,
     digits = 2,
     cex    = 1.1)
title("Variance partitioning — reef fish biomass",
      cex.main = 1.1)
dev.off()

# ── Significance tests for individual fractions ───────────────
# Tests the unique fraction of each group (conditioned on others)
# using permutation-based RDA.
cat("\n--- Significance of unique fractions ---\n")
cat("Local unique fraction:\n")
print(anova(rda(y_biomass ~ rugosity_sc + log_settlement_grav_sc +
                  Condition(connectivity_sc) +
                  Condition(log_chla_sc) +
                  Condition(log_max_dhw_sc),
                data = total_model_data)))

cat("Spatial unique fraction:\n")
print(anova(rda(y_biomass ~ connectivity_sc +
                  Condition(rugosity_sc) +
                  Condition(log_settlement_grav_sc) +
                  Condition(log_chla_sc) +
                  Condition(log_max_dhw_sc),
                data = total_model_data)))

cat("Environmental unique fraction:\n")
print(anova(rda(y_biomass ~ log_chla_sc + log_max_dhw_sc +
                  Condition(rugosity_sc) +
                  Condition(log_settlement_grav_sc) +
                  Condition(connectivity_sc),
                data = total_model_data)))

# ── Variance partition results and interpretation ─────────────
#
# Unique local fraction    [a] =  0.242  p = 0.001 ***
# Unique spatial fraction  [b] = -0.016  p = 0.904
# Unique environ fraction  [c] =  0.054  p = 0.077 .
# Total explained (adj R²)     =  0.227
# Residual                     =  0.773
#
# Local processes are the only statistically significant
# contributor to reef fish biomass variance. Rugosity and
# human pressure together explain 24.2% of variance uniquely —
# i.e. independently of connectivity and environmental context.
#
# Connectivity explains no unique variance (adj R² < 0,
# reported as ~0; F = 0.017, p = 0.904). This is not a
# borderline result — the F-statistic indicates a near-complete
# absence of independent spatial signal at this scale.
#
# Environmental context (chla + DHW) explains a modest
# independent fraction (5.4%) that does not reach conventional
# significance (p = 0.077). Reported as a trend only.
#
# All shared fractions [d–g] are near zero or negative,
# indicating that the three process groups are largely
# orthogonal in this dataset. The unique local fraction is
# therefore a reliable estimate, not a conservative lower
# bound inflated by collinearity.
#
# Negative adjusted R² values arise when the penalty for
# degrees of freedom exceeds the raw explained variance.
# They indicate the fraction is indistinguishable from zero
# and are reported as ~0 throughout.
#
# CAUTION: The absence of a spatial signal should be
# interpreted as scale-dependent, not as evidence that
# connectivity is ecologically irrelevant. At the spatial
# grain of this study, biomass is more strongly constrained
# by local carrying capacity and exploitation intensity than
# by larval network position.

# ============================================================
#  STAGE 2 — HIERARCHICAL MODEL COMPARISON
#
#  Nested sequence adds process groups progressively.
#  Tests incremental explanatory value of each group beyond
#  the local ecological baseline.
#
#  Sequence:
#    Null            → intercept only
#    Local           → habitat + pressure
#    Local + env     → adds environmental context
#    Local + spatial → adds larval connectivity
#    Local + MPA     → adds governance gradient (independent test)
#    Global          → all predictors
#
#  Both ΔAICc and ΔR² reported at each step.
# ============================================================

# ── Null model ────────────────────────────────────────────────
m_null <- lm(log_mean_biomass ~ 1,
             data = total_model_data)

# ── Local ecological baseline (Tier 1) ───────────────────────
# Habitat complexity + human pressure only.
# Represents the hypothesis that local processes are sufficient
# to explain biomass without spatial or environmental context.
m_local <- lm(log_mean_biomass ~ rugosity_sc +
                log_settlement_grav_sc,
              data = total_model_data)

# ── Local + environmental context ────────────────────────────
# Adds productivity and thermal stress to the local baseline.
# Tests whether background abiotic conditions explain additional
# variance beyond local habitat and pressure.
m_local_env <- lm(log_mean_biomass ~ rugosity_sc +
                    log_settlement_grav_sc +
                    log_chla_sc +
                    log_max_dhw_sc,
                  data = total_model_data)

# And separate
m_local_chla <- lm(log_mean_biomass ~ rugosity_sc +
                     log_settlement_grav_sc +
                     log_chla_sc,
                   data = total_model_data)

m_local_dhw  <- lm(log_mean_biomass ~ rugosity_sc +
                     log_settlement_grav_sc +
                     log_max_dhw_sc,
                   data = total_model_data)

# ── Local + spatial processes ────────────────────────────────
# Adds larval connectivity to the local baseline.
# Tests whether network position explains additional variance
# beyond local habitat and pressure.
m_local_spatial <- lm(log_mean_biomass ~ rugosity_sc +
                        log_settlement_grav_sc +
                        connectivity_sc,
                      data = total_model_data)

# ── Local + MPA status (governance test) ─────────────────────
# Adds the protection gradient to the local baseline.
# Key test: does management regime explain variance beyond
# local ecological processes alone?
# MPA treated as ordered factor: none < low < medium.
m_local_mpa <- lm(log_mean_biomass ~ rugosity_sc +
                    log_settlement_grav_sc +
                    mpa_status,
                  data = total_model_data)

# ── Global model (all processes) ─────────────────────────────
# Full additive model including all hypothesised drivers.
# Upper bound on explained variance for this predictor set.
m_global <- lm(log_mean_biomass ~ rugosity_sc +
                 log_settlement_grav_sc +
                 connectivity_sc +
                 mpa_status +
                 log_chla_sc +
                 log_max_dhw_sc,
               data = total_model_data)

# ── Hierarchical comparison table ────────────────────────────
# Listed in intended hierarchical order — ΔR² is the increment
# relative to the previous model in the sequence.
model_list <- list(
  "Null"            = m_null,
  "Local"           = m_local,
  "Local + chla"    = m_local_chla,
  "Local + DHW"     = m_local_dhw,
  "Local + env"     = m_local_env,
  "Local + MPA"     = m_local_mpa,
  "Local + spatial" = m_local_spatial,
  "Global"          = m_global
)

# AICc-ranked table — for model selection
cat("\n--- Stage 2: Model comparison (AICc ranked) ---\n")
print(make_aicc_df(model_list))

# Hierarchical table — for R² increments, in sequence order
cat("\n--- Stage 2: Variance explained (hierarchical sequence) ---\n")
local_r2 <- summary(m_local)$adj.r.squared

model_list %>%
  imap_dfr(~ tibble(
    Model  = .y,
    Adj_R2 = round(summary(.x)$adj.r.squared, 3)
  )) %>%
  mutate(Delta_R2 = round(Adj_R2 - local_r2, 3),
         Delta_R2 = ifelse(Model %in% c("Null", "Local"), NA, Delta_R2)) %>%
  arrange(Model == "Null",        # Null always last
          desc(Model == "Local"), # Local always first
          desc(Adj_R2)) %>%       # everything else by R²
  print()

# ── Key contrasts ─────────────────────────────────────────────
# Local vs Null:           Local baseline explains 18.4% of variance
#                          (ΔAICc = 8.49) — local processes clearly
#                          supported over intercept-only model.
#
# Local + env vs Local:    Adding environmental context (chla + DHW)
#                          improves fit by ΔR² = 0.059 — best supported
#                          model (ΔAICc = 0.00, weight = 0.30).
#
# Local + DHW vs Local:    DHW alone adds ΔR² = 0.038 — thermal stress
#                          is the stronger environmental driver.
#
# Local + chla vs Local:   Chla alone adds ΔR² = 0.010 — productivity
#                          contributes modestly and independently.
#
# Local + MPA vs Local:    MPA status adds ΔR² = 0.008 (ΔAICc = 3.52)
#                          — governance gradient not supported beyond
#                          local ecological baseline.
#
# Local + spatial vs Local: Connectivity adds ΔR² = -0.011 (ΔAICc = 3.30)
#                          — spatial processes add no explanatory value.
#
# Global vs Local:         Full model adds ΔR² = 0.056 over Local but
#                          ΔAICc = 5.11 — gains driven by environmental
#                          terms only; MPA and connectivity contribute
#                          dead-weight parameters. Global is the worst
#                          supported model after Null.
#
# Best-supported model: Local + env (AICc = 99.43, weight = 0.30)
# Model selection uncertainty: Local + DHW competitive (ΔAICc = 0.02)
# Combined weight of Local + env and Local + DHW = 0.59
# ΔAICc > 2 threshold: only Local, Local + env, Local + DHW supported

# ── Coefficient summary: best-supported model (Local + env) ───
summary(m_local_env)   # already done
summary(m_local_dhw)   
summary(m_local)

# ── Results: best-supported model (Local + env) ───────────────
# lm(log_mean_biomass ~ rugosity_sc + log_settlement_grav_sc +
#    log_chla_sc + log_max_dhw_sc)
# n = 54 sites, 4 countries
#
# Coefficients:
#                        Estimate Std. Error t value Pr(>|t|)
# (Intercept)             9.6539     0.0765  126.23  < 0.001 ***
# rugosity_sc             0.2095     0.0773    2.71   0.0092 **
# log_settlement_grav_sc -0.2794     0.0947   -2.95   0.0049 **
# log_chla_sc            -0.1476     0.0958   -1.54   0.1297
# log_max_dhw_sc          0.1611     0.0783    2.06   0.0451 *
#
# Residual standard error: 0.561 on 49 df
# Adjusted R² = 0.243
# F(4, 49) = 5.25, p = 0.001
#
# Coefficient stability across competitive models (ΔAICc < 2):
#
#                      Local      Local+DHW   Local+env
# rugosity_sc         +0.221**   +0.217**    +0.209**
# settlement_grav_sc  -0.184*    -0.198*     -0.279**
# log_max_dhw_sc          —      +0.147.     +0.161*
# log_chla_sc             —          —       -0.148 ns
#
# Rugosity and settlement gravity are stable in direction and
# magnitude across all three competitive models, confirming
# that the local signal is robust and does not depend on which
# environmental terms are included.
#
# Interpretation:
#
# Rugosity (β = +0.209, p = 0.009):
#   Positive effect — greater structural complexity associated
#   with higher biomass, consistent with rugosity setting local
#   carrying capacity through provision of refuge and foraging
#   habitat. Effect is significant and stable across all
#   competitive models.
#
# Settlement gravity (β = -0.279, p = 0.005):
#   Negative effect — higher human pressure associated with
#   lower biomass, consistent with exploitation reducing
#   standing stock. Effect strengthens slightly when
#   environmental terms are included, suggesting modest
#   suppression by collinearity with chla (r = 0.57).
#   Significant and stable across all competitive models.
#
# Chla (β = -0.148, p = 0.130):
#   Non-significant. Negative direction is counterintuitive
#   and likely reflects moderate collinearity with settlement
#   gravity (r = 0.57), which destabilises the chla coefficient
#   when both are included. Interpret with caution.
#
# DHW (β = +0.161, p = 0.045):
#   Positive and significant — counterintuitive given expected
#   negative effect of thermal stress on coral and fish biomass.
#   Not driven by zeros (15 zero-DHW sites have near-identical
#   mean biomass to non-zero sites: 9.62 vs 9.67) and not
#   confounded with fishing pressure (r = 0.10 with settlement
#   gravity). Most likely reflects geographic covariation —
#   higher DHW sites tend to be in warmer, more tropical
#   environments that support greater reef fish biomass
#   independently of bleaching history. Interpret as a
#   geographic context variable rather than a direct
#   thermal stress effect. Address in discussion.

# ============================================================
#  STAGE 3 — INTERACTION TESTING
#
#  Tests whether spatial management (MPA) and connectivity
#  modify local-scale relationships — the mechanistic
#  hypothesis underlying the Tier 2 structure.
#
#  Run only if Global outperforms Local by ΔAICc > 2.
#  If the Tier 2 variables add little, interactions are
#  untestable and this stage is reported as not supported.
#
#  Three a priori interactions, each representing a distinct
#  mechanistic hypothesis:
#    MPA × connectivity:  effectiveness depends on larval supply
#    MPA × pressure:      effectiveness depends on fishing intensity
#    Connectivity × pressure: connectivity buffers fishing impact
#
#  Compared against the additive Global model as reference.
#  Interaction preferred only if ΔAICc > 2 vs Global additive.
# ============================================================

# ── Check Stage 2 before proceeding ──────────────────────────
delta_global_vs_local <- AICc(m_local) - AICc(m_global)
cat("\nΔAICc (Local vs Global):", round(delta_global_vs_local, 2), "\n")

# ── Gate check result ─────────────────────────────────────────
# ΔAICc (Local vs Global) = -3.88
# Global is WORSE than Local by 3.88 AICc units.
# MPA and connectivity add no explanatory value as main effects.
# Interaction models are fitted for completeness and to satisfy
# the a priori analytical framework, but any significant
# interaction should be interpreted with considerable caution —
# it would be difficult to justify a context-dependent effect
# of variables whose main effects are not supported.
# Results reported in supplementary material if no interaction
# is supported; in main text only if ΔAICc > 2 vs Local + env.

# ── Hypothesis 1: MPA effectiveness depends on larval supply ──
# Protected sites with high connectivity should recover faster
# and maintain higher biomass than isolated protected sites
m_int_mpa_conn <- lm(log_mean_biomass ~ rugosity_sc +
                       log_settlement_grav_sc +
                       log_chla_sc +
                       log_max_dhw_sc +
                       mpa_status * connectivity_sc,
                     data = total_model_data)

# ── Hypothesis 2: MPA effectiveness depends on fishing intensity
# MPAs should only be detectable where external pressure is low
# enough that protection translates into a biomass difference
m_int_mpa_press <- lm(log_mean_biomass ~ rugosity_sc +
                        log_chla_sc +
                        log_max_dhw_sc +
                        mpa_status * log_settlement_grav_sc,
                      data = total_model_data)

# ── Hypothesis 3: connectivity buffers fishing pressure ───────
# Well-connected sites should be more resilient to exploitation
# through sustained larval replenishment offsetting mortality
m_int_conn_press <- lm(log_mean_biomass ~ rugosity_sc +
                         log_chla_sc +
                         log_max_dhw_sc +
                         connectivity_sc * log_settlement_grav_sc,
                       data = total_model_data)

# ── Interaction candidate set ─────────────────────────────────
model_list_interactions <- list(
  "Local + env (additive)"  = m_local_env,
  "MPA × connectivity"      = m_int_mpa_conn,
  "MPA × pressure"          = m_int_mpa_press,
  "Connectivity × pressure" = m_int_conn_press
)

cat("\n--- Stage 3: Interaction model comparison ---\n")
print(make_aicc_df(model_list_interactions))


# ── Interaction candidate set ─────────────────────────────────
# Global additive model is the reference — we ask whether any
# interaction improves on the additive structure.
model_list_interactions <- list(
  "Global (additive)"       = m_global,
  "MPA × connectivity"      = m_int_mpa_conn,
  "MPA × pressure"          = m_int_mpa_press,
  "Connectivity × pressure" = m_int_conn_press
)

cat("\n--- Stage 3: Interaction model comparison ---\n")
print(make_aicc_df(model_list_interactions))

# ── Interpretation ────────────────────────────────────────────
# An interaction model is preferred only if ΔAICc > 2 vs Global
# additive. If no interaction is supported, report this as
# evidence that Tier 2 effects are additive rather than
# context-dependent at this spatial scale.

summary(m_int_mpa_press)
summary(m_int_conn_press)

# ── Stage 3 conclusion ────────────────────────────────────────
# No interaction model provides a convincing, ecologically
# interpretable improvement over Local + env.
#
# MPA × pressure (ΔAICc = 1.29) is statistically competitive
# but the interaction pattern is not consistent with the a
# priori hypothesis. The main effect of settlement gravity
# collapses entirely within the interaction model, suggesting
# the signal reflects geographic covariation rather than a
# genuine management-by-pressure effect.
#
# Connectivity × pressure (ΔAICc = 1.81, p = 0.079) is
# marginal and not supported at conventional thresholds.
#
# MPA × connectivity (ΔAICc = 5.13) is not supported.
#
# Primary inference is unchanged: reef fish biomass is
# structured by local ecological processes. Spatial management
# and connectivity do not modify this relationship in a
# statistically or ecologically robust way at this scale.

# ============================================================
#  STAGE 4 — SENSITIVITY ANALYSIS
#
#  (a) Alternative pressure metrics
#      Best Stage 2 model refitted substituting settlement pop.
#      and market gravity for settlement gravity. All other
#      predictors retained to ensure comparability.
#      Results reported in Supplementary Table S2.
#
#  (b) Transect-level replication
#      Site-level hierarchical sequence repeated at transect
#      level using GLMM with (1 | site). Confirms site-level
#      findings are not an artefact of averaging.
#      Results reported in Supplementary Table S1.
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
# Mirrors the Global model structure — only pressure metric swapped.
# NOTE: earlier sensitivity models dropped env. variables,
# making them non-comparable. These are corrected here.

sens_settpop <- lm(log_mean_biomass ~ rugosity_sc +
                     log_settlement_pop_sc +      # swapped
                     connectivity_sc +
                     mpa_status +
                     log_chla_sc +
                     log_max_dhw_sc,
                   data = total_model_data)

sens_mktgrav <- lm(log_mean_biomass ~ rugosity_sc +
                     log_market_gravity_sc +      # swapped
                     connectivity_sc +
                     mpa_status +
                     log_chla_sc +
                     log_max_dhw_sc,
                   data = total_model_data)

cat("\n--- Stage 4a: Sensitivity — alternative pressure metrics ---\n")
cat("(Coefficients only — these are not compared via AICc\n",
    "because they test robustness, not process importance)\n\n")
cat("Settlement population:\n"); print(summary(sens_settpop)$coefficients)
cat("\nMarket gravity:\n");      print(summary(sens_mktgrav)$coefficients)

# ── Stage 4a results ──────────────────────────────────────────
# Robustness check: primary conclusions hold across alternative
# pressure metrics.
#
# Rugosity coefficient stability:
#   Primary (settlement gravity): β = +0.209, p = 0.009
#   Settlement population:        β = +0.210, p = 0.014
#   Market gravity:               β = +0.192, p = 0.024
#
# Rugosity is significant and stable in direction and magnitude
# across all three pressure metrics. The habitat complexity
# signal is robust to pressure proxy choice.
#
# Pressure metric performance:
#   Settlement population: β = -0.110, p = 0.258 — not significant
#   Market gravity:        β = -0.175, p = 0.081 — marginal
#
# Neither alternative reaches significance, confirming
# settlement gravity as the strongest pressure proxy.
#
# All other predictors (connectivity, MPA, chla) non-significant
# across both models — consistent with Stage 2.
#
# DHW positive in both models (β = +0.144 and +0.209) —
# consistent with primary model interpretation.
#
# Conclusion: sensitivity analysis supports primary inference.
# Rugosity is the most robust predictor regardless of pressure
# metric. Human pressure direction is consistent but weakens
# with alternative metrics, further justifying settlement
# gravity selection.

# ── (b) Transect-level sensitivity ───────────────────────────
# Replicates the Stage 2 hierarchical sequence at transect level.
# (1 | site) accounts for non-independence of transects within sites.

sens_t_null <- glmmTMB(log_transect_biomass ~ 1 +
                         (1 | site),
                       family = gaussian(),
                       data   = transect_model_data)

sens_t_local <- glmmTMB(log_transect_biomass ~ rugosity_sc +
                          log_settlement_grav_sc +
                          (1 | site),
                        family = gaussian(),
                        data   = transect_model_data)

sens_t_local_env <- glmmTMB(log_transect_biomass ~ rugosity_sc +
                              log_settlement_grav_sc +
                              log_chla_sc +
                              log_max_dhw_sc +
                              (1 | site),
                            family = gaussian(),
                            data   = transect_model_data)

sens_t_local_spatial <- glmmTMB(log_transect_biomass ~ rugosity_sc +
                                  log_settlement_grav_sc +
                                  connectivity_sc +
                                  (1 | site),
                                family = gaussian(),
                                data   = transect_model_data)

sens_t_local_mpa <- glmmTMB(log_transect_biomass ~ rugosity_sc +
                              log_settlement_grav_sc +
                              mpa_status +
                              (1 | site),
                            family = gaussian(),
                            data   = transect_model_data)

sens_t_global <- glmmTMB(log_transect_biomass ~ rugosity_sc +
                           log_settlement_grav_sc +
                           connectivity_sc +
                           mpa_status +
                           log_chla_sc +
                           log_max_dhw_sc +
                           (1 | site),
                         family = gaussian(),
                         data   = transect_model_data)

model_list_transect <- list(
  "Null"            = sens_t_null,
  "Local"           = sens_t_local,
  "Local + env"     = sens_t_local_env,
  "Local + spatial" = sens_t_local_spatial,
  "Local + MPA"     = sens_t_local_mpa,
  "Global"          = sens_t_global
)

cat("\n--- Stage 4b: Sensitivity — transect-level ---\n")
print(make_aicc_df(model_list_transect))

# ── Stage 4b results ──────────────────────────────────────────
# Transect-level mixed model (1 | site) replicates site-level
# findings — conclusions are not an artefact of aggregation.
#
# Model ranking (transect-level):
#   Local + env:     AICc = 634.92, weight = 0.663 — best supported
#   Local:           ΔAICc = 2.83,  weight = 0.161 — competitive
#   Local + spatial: ΔAICc = 4.63  — not supported
#   Global:          ΔAICc = 4.82  — not supported
#   Local + MPA:     ΔAICc = 5.16  — not supported
#   Null:            ΔAICc = 13.06 — worst supported
#
# Comparison with site-level ranking:
#   Site-level:      Local + env best (weight = 0.297)
#   Transect-level:  Local + env best (weight = 0.663)
#
# Model ordering is identical across both levels of analysis.
# Local + env is best supported at both levels. Local + spatial,
# Local + MPA, and Global are all unsupported at both levels.
# The stronger weight at transect level reflects greater
# statistical power from 243 observations vs 54 sites.
#
# Conclusion: site-level aggregation does not alter qualitative
# inference. Local ecological processes (habitat complexity and
# human pressure) with environmental context (chla + DHW)
# provide the best explanation of reef fish biomass at both
# the transect and site level. Spatial and governance variables
# add no explanatory value at either level of analysis.

# ============================================================
#  MARGINAL EFFECT PLOTS
#  Generated for primary predictors from best supported model:
#  Local + env (AICc = 99.43, weight = 0.297)
#  All other standardised predictors held at their mean (= 0).
# ============================================================

best_model <- m_local_env

# ── Continuous predictors ─────────────────────────────────────
p_rugosity <- plot_effect(best_model, total_model_data,
                          "rugosity_sc",
                          "Rugosity (standardised)")

p_pressure <- plot_effect(best_model, total_model_data,
                          "log_settlement_grav_sc",
                          "log(Settlement gravity) (standardised)")

p_chla <- plot_effect(best_model, total_model_data,
                      "log_chla_sc",
                      "log(Chlorophyll-a) (standardised)")

p_dhw <- plot_effect(best_model, total_model_data,
                     "log_max_dhw_sc",
                     "log(Max DHW + 1) (standardised)")

# Connectivity excluded from main marginal effect plots —
# not included in best supported model and adds no explanatory
# value (Stage 2: ΔAICc = 3.30, Stage 1: unique fraction ≈ 0)

jpeg("marginal_effects_main.jpg",
     width = 33, height = 11, units = "cm", res = 300)
gridExtra::grid.arrange(p_rugosity, p_pressure,
                        p_chla,    p_dhw,
                        ncol = 4)
dev.off()

# ── MPA marginal means ────────────────────────────────────────
# MPA status not included in best supported model (Local + env).
# Marginal means generated from m_local_mpa for completeness
# and to report the direction of the protection gradient.
# Interpret with caution — MPA not supported in Stage 2
# (ΔAICc = 3.52 vs Local).
# Unordered factor — coefficients are dummy contrasts:
#   low vs none, medium vs none.

mpa_grid <- data.frame(
  mpa_status             = factor(c("none", "low", "medium"),
                                  levels  = c("none", "low", "medium"),
                                  ordered = FALSE),
  rugosity_sc            = 0,
  log_settlement_grav_sc = 0,
  log_chla_sc            = 0,
  log_max_dhw_sc         = 0
)

mpa_pred <- predict(m_local_mpa, newdata = mpa_grid, se.fit = TRUE)
mpa_grid$fit <- mpa_pred$fit
mpa_grid$lwr <- mpa_pred$fit - 1.96 * mpa_pred$se.fit
mpa_grid$upr <- mpa_pred$fit + 1.96 * mpa_pred$se.fit

p_mpa <- ggplot(mpa_grid, aes(x = mpa_status, y = fit)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  colour = "#2c7bb6", size = 0.8) +
  labs(x = "MPA status",
       y = "Fitted log(biomass)",
       caption = "From m_local_mpa — MPA not supported in Stage 2") +
  theme_bw(base_size = 13) +
  theme(axis.title   = element_text(face = "bold"),
        plot.caption = element_text(colour = "grey50", size = 9))

jpeg("marginal_effect_mpa.jpg",
     width = 12, height = 12, units = "cm", res = 300)
print(p_mpa)
dev.off()

# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()

#-------------------------------------------------------------------------------
# ============================================================
#  TRANSECT-LEVEL COUNTS (COMPLEMENTARY ANALYSIS)
#
#  Rationale: Biomass integrates both fish abundance and body
#  size via length-weight relationships. Modelling counts
#  directly allows decomposition of biomass patterns into
#  abundance vs body size components:
#
#    Predictor in both biomass and count models
#    → effect operates through fish abundance
#
#    Predictor in biomass but not count models
#    → effect operates through individual body size
#
#    Predictor in count but not biomass models
#    → effect on abundance is masked in biomass by
#      compensatory changes in body size
#
#  This decomposition directly addresses whether MPA protection
#  and connectivity operate through recovering fish numbers,
#  recovering large individuals, or both.
#
#  Response: Total fish count per transect (integer >= 0)
#  Family:   Poisson → NB2 → NB1 — selected via AICc + DHARMa
#  Random fx: (1 | site)
#  Model ladder: identical to biomass analyses for comparability
# ============================================================

# ── Explore count distribution ────────────────────────────────
cat("Transects:", nrow(transect_model_data), "\n")
cat("Zero counts:", sum(transect_model_data$transect_total_count == 0), "\n")
cat("Proportion zeros:", round(mean(transect_model_data$transect_total_count == 0), 3), "\n")

summary(transect_model_data$transect_total_count)

ggplot(transect_model_data, aes(x = transect_total_count)) +
  geom_histogram(bins = 50, fill = "#2c7bb6", colour = "white") +
  labs(x = "Total fish count per transect", y = "Frequency") +
  theme_bw()

# ── Mean-variance relationship ────────────────────────────────
transect_model_data %>%
  group_by(site) %>%
  summarise(mean_count = mean(transect_total_count),
            var_count = var(transect_total_count),
            .groups = "drop") %>%
  ggplot(aes(x = mean_count, y = var_count)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Site mean count", y = "Site variance",
       title = "Mean-variance (red = Poisson expectation)") +
  theme_bw()

# ── Family selection ──────────────────────────────────────────
# AICc directly comparable across Poisson, NB1, NB2 (same response)

m_count_poisson <- glmmTMB(
  transect_total_count ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    log_max_dhw_sc +
    (1 | site),
  family = poisson(link = "log"),
  data = transect_model_data
)

m_count_nb2 <- glmmTMB(
  transect_total_count ~ rugosity_sc +
    log_settlement_grav_sc +
    log_chla_sc +
    log_max_dhw_sc +
    (1 | site),
  family = nbinom2(link = "log"),
  data = transect_model_data
)

m_count_nb1 <- glmmTMB(
  transect_total_count ~ rugosity_sc +
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

jpeg("dharma_count_poisson.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_poisson, main = "DHARMa — Poisson"); dev.off()

jpeg("dharma_count_nb2.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_nb2, main = "DHARMa — NB2"); dev.off()

jpeg("dharma_count_nb1.jpg", width = 25, height = 15, units = "cm", res = 300)
plot(res_nb1, main = "DHARMa — NB1"); dev.off()

plot(res_poisson); testDispersion(res_poisson); testZeroInflation(res_poisson)
plot(res_nb2);     testDispersion(res_nb2);     testZeroInflation(res_nb2)
plot(res_nb1);     testDispersion(res_nb1);     testZeroInflation(res_nb1)

cat("\n--- Family selection: count models ---\n")
print(make_aicc_df(list(
  "Poisson" = m_count_poisson,
  "NB2" = m_count_nb2,
  "NB1" = m_count_nb1
)))

# ── Family selection decision ───────────────────────────────── update specifics...
# NB2 selected — confirmed by both AICc and DHARMa diagnostics.
#
# Poisson (image 1): rejected decisively
#   KS test p = 0.0002, dispersion p = 0.006, outlier p = 0
#   Dispersion = 2.20 — severe overdispersion
#   QQ plot shows systematic deviation throughout
#   Residuals vs predicted shows strong quantile violations
#
# NB2 (image 2): selected
#   KS p = 0.952, dispersion p = 0.848, outlier p = 0.98
#   All tests n.s. — no significant problems detected
#   Dispersion = 1.005 — near perfect
#   QQ plot closely follows theoretical line
#   Residuals vs predicted flat with no systematic pattern
#   Zero inflation: NaN p = 1 — no zero inflation present
#
# NB1 (image 3): adequate but inferior to NB2
#   KS p = 0.659, dispersion p = 0.174, outlier p = 0.32
#   Formal tests pass but residuals vs predicted shows
#   quantile deviations flagged by DHARMa — upper and lower
#   quantile curves show systematic pattern not present in NB2
#   AICc delta = 4.48 vs NB2 — not competitive
#
# Proceed with nbinom2(link = "log") + (1 | site) throughout
# count analyses.

# ── RE structure ──────────────────────────────────────────────
# (1 | site) tested against no random effect structure
# Test confirms this is strongly supported.
re_c_null <- glmmTMB(transect_total_count ~ rugosity_sc +
                       log_settlement_grav_sc +
                       log_chla_sc +
                       log_max_dhw_sc,
                     family = nbinom2(link = "log"),
                     data = transect_model_data)

re_c_site <- glmmTMB(transect_total_count ~ rugosity_sc +
                       log_settlement_grav_sc +
                       log_chla_sc +
                       log_max_dhw_sc +
                       (1 | site),
                     family = nbinom2(link = "log"),
                     data = transect_model_data)

cat("\n--- RE structure: count models ---\n")
print(make_aicc_df(list(
  "No RE" = re_c_null,
  "(1 | site)" = re_c_site
)))

# ── RE structure decision ─────────────────────────────────────
# (1 | site) strongly supported (difference in AICc = 27.73, weight = 1.00).
# Even stronger site-level clustering in counts than in biomass
# (delta 27.73 vs 9.89) — fish counts are more variable within
# sites than biomass, making the site random effect essential.
# Proceed with (1 | site) throughout count analyses.

# ── Candidate models — identical ladder to biomass ────────────
count_family <- nbinom2(link = "log")

# Model 1: Habitat only
count_m1_hab <- glmmTMB(transect_total_count ~ rugosity_sc +
                      (1 | site),
                    family = count_family, data = transect_model_data)

# Model 2: Habitat + pressure
count_m2_hab_press <- glmmTMB(transect_total_count ~ rugosity_sc +
                            log_settlement_grav_sc +
                            (1 | site),
                          family = count_family, data = transect_model_data)

# Model 3: Habitat + pressure + MPA
count_m3_hab_press_mpa <- glmmTMB(transect_total_count ~ rugosity_sc +
                                log_settlement_grav_sc +
                                mpa_status +
                                (1 | site),
                              family = count_family, data = transect_model_data)

# Model 4: Above + connectivity
count_m4_conn <- glmmTMB(transect_total_count ~ rugosity_sc +
                       log_settlement_grav_sc +
                       mpa_status +
                       connectivity_sc +
                       (1 | site),
                     family = count_family, data = transect_model_data)

# Model 5: Above + chla
count_m5_chla <- glmmTMB(transect_total_count ~ rugosity_sc +
                       log_settlement_grav_sc +
                       mpa_status +
                       connectivity_sc +
                       log_chla_sc +
                       (1 | site),
                     family = count_family, data = transect_model_data)

# Model 6: Above + DHW
count_m6_dhw <- glmmTMB(transect_total_count ~ rugosity_sc +
                      log_settlement_grav_sc +
                      mpa_status +
                      connectivity_sc +
                      log_max_dhw_sc +
                      (1 | site),
                    family = count_family, data = transect_model_data)

# Model 7: MPA x connectivity
count_m7_mpa_conn <- glmmTMB(transect_total_count ~ rugosity_sc +
                           log_settlement_grav_sc +
                           mpa_status * connectivity_sc +
                           (1 | site),
                         family = count_family, data = transect_model_data)

# Model 8: MPA x pressure
count_m8_mpa_press <- glmmTMB(transect_total_count ~ rugosity_sc +
                            mpa_status * log_settlement_grav_sc +
                            connectivity_sc +
                            (1 | site),
                          family = count_family, data = transect_model_data)

# Model 9: Connectivity x pressure
count_m9_conn_press <- glmmTMB(transect_total_count ~ rugosity_sc +
                             mpa_status +
                             connectivity_sc * log_settlement_grav_sc +
                             (1 | site),
                           family = count_family, data = transect_model_data)

# Sensitivity
c_sens_settpop <- glmmTMB(transect_total_count ~ rugosity_sc +
                            log_settlement_pop_sc +
                            mpa_status +
                            connectivity_sc +
                            (1 | site),
                          family = count_family, data = transect_model_data)

c_sens_mktgrav <- glmmTMB(transect_total_count ~ rugosity_sc +
                            log_market_gravity_sc +
                            mpa_status +
                            connectivity_sc +
                            (1 | site),
                          family = count_family, data = transect_model_data)

# ── Model list ────────────────────────────────────────────────
model_list_counts <- list(
  "Habitat" = count_m1_hab,
  "Habitat + pressure" = count_m2_hab_press,
  "Habitat + pressure + MPA" = count_m3_hab_press_mpa,
  "Above + connectivity" = count_m4_conn,
  "Above + chla" = count_m5_chla,
  "Above + DHW" = count_m6_dhw,
  "MPA x connectivity" = count_m7_mpa_conn,
  "MPA x pressure" = count_m8_mpa_press,
  "Connectivity x pressure" = count_m9_conn_press,
  "Settlement pop. (sensitivity)" = c_sens_settpop,
  "Market gravity (sensitivity)" = c_sens_mktgrav
)

cat("\n--- AICc: Count models ---\n")
print(make_aicc_df(model_list_counts))

# ── Count model results ───────────────────────────────────────
#
#                         Model   AICc Delta Weight
#                   Above + DHW 2481.03  0.00  0.21
#  Market gravity (sensitivity) 2481.08  0.05  0.21
# Settlement pop. (sensitivity) 2481.71  0.68  0.15
#      Habitat + pressure + MPA 2481.82  0.79  0.14
#          Above + connectivity 2482.44  1.40  0.11
#       Connectivity x pressure 2484.08  3.05  0.05
#                       Habitat 2484.39  3.35  0.04
#                  Above + chla 2484.57  3.54  0.04
#            MPA x connectivity 2485.18  4.14  0.03
#              Habitat + pressure 2486.32  5.29  0.02
#                  MPA x pressure 2486.36  5.33  0.01
#
# STRIKING DIVERGENCE FROM BIOMASS ANALYSES:
#
# 1. HABITAT + PRESSURE is no longer the best model — it sits
#    at delta 5.29, essentially unsupported for counts.
#    This is the opposite of biomass where it was best by far.
#
# 2. DHW is the best-supported predictor for counts (delta 0.00)
#    but was not supported for biomass (delta 4.27 at site level).
#    Thermal stress appears to operate primarily through fish
#    abundance rather than body size — bleaching events reduce
#    fish numbers without necessarily reducing the size of
#    remaining individuals.
#
# 3. SENSITIVITY METRICS (market gravity, settlement pop.) are
#    within delta < 2 for counts but not for biomass — the
#    choice of pressure metric matters more for abundance than
#    for biomass. This may reflect that different pressure types
#    (market access vs local settlement) target different
#    components of the fish community.
#
# 4. MPA STATUS has marginal support as a main effect (delta 0.79)
#    — some evidence protection increases fish numbers.
#
# 5. MPA INTERACTIONS not supported for counts (delta > 4) —
#    the MPA x pressure interaction that was prominent in
#    biomass analyses does not appear in counts. This suggests
#    the interaction operates primarily through body size
#    recovery rather than abundance recovery at protected sites.
#
# 6. RUGOSITY is not driving counts as strongly as biomass —
#    habitat complexity appears to support larger-bodied fish
#    rather than simply more fish, consistent with the idea
#    that structural complexity provides refuge for large
#    individuals rather than increasing total recruitment.
#
# KEY BIOLOGICAL INTERPRETATION:
#   Biomass and abundance are structured by fundamentally
#   different processes in this system:
#   - Biomass: driven by habitat complexity and human pressure,
#     with MPA x pressure interaction suggesting protection
#     modifies the pressure-body size relationship
#   - Abundance: driven by thermal stress and MPA status,
#     suggesting fish numbers respond to disturbance history
#     and protection independently of habitat structure
#   This decomposition suggests MPAs recover fish communities
#   primarily through protecting individual body size rather
#   than increasing total fish numbers — consistent with
#   size-selective fishing pressure targeting large individuals.

# ============================================================
#  SYNTHESIS: BIOMASS vs COUNT MODEL CONCLUSIONS - update...
#
#  Site-level biomass best models (delta < 2):
#    1. Habitat + pressure (weight 0.31)
#    2. MPA x pressure (weight 0.18)
#
#  Transect-level biomass best models (delta < 2):
#    1. Habitat + pressure (weight 0.33)
#    2. MPA x connectivity (weight 0.14)
#    3. Above + chla (weight 0.13)
#
#  Count models best models (delta < 2):
#    1. Above + DHW (weight 0.21)
#    2. Market gravity sensitivity (weight 0.21)
#    3. Settlement pop. sensitivity (weight 0.15)
#    4. Habitat + pressure + MPA (weight 0.14)
#    5. Above + connectivity (weight 0.11)
#
#  Key finding: biomass and abundance are structured by
#  fundamentally different processes — the dominant predictors
#  for each response do not overlap.
# ============================================================

# ── Rugosity ──────────────────────────────────────────────────
# Extract rugosity coefficients across all analyses
rug_site     <- coef(s1_m2_hab_press)["rugosity_sc"]
rug_transect <- fixef(t_m2_hab_press)$cond["rugosity_sc"]
rug_count    <- fixef(c_m4_conn)$cond["rugosity_sc"]

cat("\n--- Rugosity effect across analyses ---\n")
cat("Site biomass:      beta =", round(rug_site, 3), "\n")
cat("Transect biomass:  beta =", round(rug_transect, 3), "\n")
cat("Counts:            beta =", round(rug_count, 3),
    " IRR =", round(exp(rug_count), 3), "\n")

# Rugosity beta for biomass > counts at both levels — complex
# reefs support larger-bodied fish, not just more fish.
# Structural complexity provides refuge for large individuals
# rather than simply increasing total recruitment.

# ── Human pressure ────────────────────────────────────────────
press_site     <- coef(s1_m2_hab_press)["log_settlement_grav_sc"]
press_transect <- fixef(t_m2_hab_press)$cond["log_settlement_grav_sc"]
press_count    <- fixef(c_m4_conn)$cond["log_settlement_grav_sc"]

cat("\n--- Human pressure (settlement gravity) ---\n")
cat("Site biomass:      beta =", round(press_site, 3), "\n")
cat("Transect biomass:  beta =", round(press_transect, 3), "\n")
cat("Counts:            beta =", round(press_count, 3),
    " IRR =", round(exp(press_count), 3), "\n")

# Habitat + pressure is the best biomass model at both levels
# but sits at delta 5.29 for counts — pressure operates
# primarily through body size reduction (size-selective
# harvesting removes large individuals) rather than through
# reducing total fish numbers.

# ── MPA effects ───────────────────────────────────────────────
cat("\n--- MPA effects ---\n")
cat("Site biomass — MPA x pressure: delta AICc = 1.10 (within delta 2)\n")
cat("Transect biomass — MPA x connectivity: delta AICc = 1.74 (within delta 2)\n")
cat("Counts — MPA main effect: delta AICc = 0.79 (within delta 2)\n")
cat("Counts — MPA x pressure: delta AICc = 5.33 (not supported)\n")
cat("Counts — MPA x connectivity: delta AICc = 4.14 (not supported)\n")

# MPA interactions supported in biomass but not counts —
# protection appears to operate primarily through recovering
# individual body size rather than fish numbers. The MPA x
# pressure interaction in biomass suggests medium-protection
# sites buffer the negative body size effects of fishing
# pressure. The absence of this interaction in counts confirms
# the effect is on size structure not abundance.

# ── DHW ───────────────────────────────────────────────────────
cat("\n--- Thermal stress (DHW) ---\n")
cat("Site biomass:      delta AICc = 2.47 (marginal)\n")
cat("Transect biomass:  delta AICc = 4.27 (not supported)\n")
cat("Counts:            delta AICc = 0.00 (BEST MODEL)\n")

# DHW is the best predictor for counts but not for biomass —
# thermal stress reduces fish numbers without proportionally
# reducing total biomass. Bleaching events may remove small
# and large fish equally in terms of numbers, but surviving
# fish (or recruits) maintain community biomass through
# compensatory growth or community restructuring.

# ── Connectivity ──────────────────────────────────────────────
cat("\n--- Connectivity ---\n")
cat("Site biomass — MPA x connectivity: delta AICc = 2.05 (marginal)\n")
cat("Transect biomass — MPA x connectivity: delta AICc = 1.74 (within delta 2)\n")
cat("Counts — connectivity main effect: delta AICc = 1.40 (within delta 2)\n")

# Connectivity appears in the top model sets for both biomass
# and counts but in different forms — as an interaction with
# MPA status for biomass (medium-protection sites benefit more
# from connectivity) and as a main effect for counts (more
# connected sites have higher fish numbers regardless of
# protection). This suggests connectivity supports fish
# abundance generally via larval supply, while its biomass
# effect is contingent on protection status.

# ── Overall synthesis ─────────────────────────────────────────
cat("\n--- Overall synthesis ---\n")
cat("
Biomass is primarily structured by local habitat complexity
and human pressure, with MPA protection modifying the
pressure-biomass relationship — particularly at medium-
protection sites where the negative pressure effect on
body size is attenuated or reversed.

Fish abundance is primarily structured by thermal stress
history and MPA status, with connectivity providing
additional larval supply benefits. The absence of MPA
interactions in count models suggests protection operates
through recovering individual body size rather than
increasing total fish numbers.

Together these results indicate that reef fish community
recovery under protection is a body-size mediated process:
MPAs protect large individuals from size-selective fishing
pressure, generating biomass recovery that is not detectable
as an increase in total fish numbers. Connectivity amplifies
this effect at medium-protection sites by maintaining larval
supply to recovering reefs.
")