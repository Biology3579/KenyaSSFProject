# ============================================================
#  DRIVERS OF PISCIVORE BIOMASS
#  Chapter 1 — Functional Group Analysis: Piscivores
#
#  Piscivores are fish-eating predators at high trophic
#  positions. Among the most heavily targeted functional
#  groups by artisanal and recreational fishing.
#
#  Analytical framework mirrors all other functional groups.
#  Structure: Q1 (pressure) → Q2 (connectivity) → Q3 (MPA)
#
#  Key note: prior analysis found market gravity preferred
#  in Q1 — first group besides total biomass with decisive
#  metric differentiation. Confirm with current data.
#
#  Key differences from other groups:
#    ~11% zeros at site level — Tweedie required.
#    ~47% zeros at transect level — ZI Tweedie tested.
#
#  Sensitivity analyses:
#    (a) Alternative pressure metrics
#    (b) Transect-level replication (Tweedie GLMM)
# ============================================================

source(here::here("data_preparation.R"))


# ============================================================
#  DATA AGGREGATION
# ============================================================

pisc_transects <- fish_data %>%
  group_by(site, station, ts_no, date) %>%
  summarise(
    transect_pisc_biomass = sum(
      ifelse(trophic_group == "piscivores", tot_wt_g, 0),
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

cat("Piscivore transects:", nrow(pisc_transects), "\n")
cat("Sites:",               n_distinct(pisc_transects$site), "\n")

# ── Site-level dataset ────────────────────────────────────────
pisc_model_data <- pisc_transects %>%
  left_join(final_predictors, by = "site") %>%
  group_by(site, country) %>%
  summarise(
    mean_biomass           = mean(transect_pisc_biomass,
                                  na.rm = TRUE),
    n_transects            = n(),
    rugosity_sc            = first(rugosity_sc),
    log_settlement_grav_sc = first(log_settlement_grav_sc),
    log_market_gravity_sc  = first(log_market_gravity_sc),
    log_chla_sc            = first(log_chla_sc),
    connectivity_sc        = first(connectivity_sc),
    mpa_status             = first(mpa_status),
    ecoregion              = first(ecoregion),
    log_settlement_pop_sc  = first(log_settlement_pop_sc),
    .groups = "drop"
  ) %>%
  mutate(
    site       = as.factor(site),
    country    = as.factor(country),
    ecoregion  = as.factor(ecoregion),
    mpa_status = factor(mpa_status,
                        levels  = c("none", "low", "medium"),
                        ordered = FALSE)
  )

cat("\nPiscivore model data:", nrow(pisc_model_data), "sites\n")

# ── Data checks ───────────────────────────────────────────────
pisc_model_data %>%
  dplyr::select(site, rugosity_sc, log_settlement_grav_sc,
                log_chla_sc, connectivity_sc, mpa_status) %>%
  filter(if_any(everything(), is.na)) %>%
  print(n = Inf)

cat("\nZeros in mean_biomass:",
    sum(pisc_model_data$mean_biomass == 0), "\n")
cat("Site-level zero proportion:",
    round(mean(pisc_model_data$mean_biomass == 0), 3), "\n")
cat("\nResponse summary:\n")
print(summary(pisc_model_data$mean_biomass))

cat("\nMPA status counts:\n")
print(table(pisc_model_data$mpa_status))

# ── Transect-level dataset ────────────────────────────────────
pisc_transect_data <- pisc_transects %>%
  left_join(final_predictors, by = "site")

cat("\nTransect zeros:",
    sum(pisc_transect_data$transect_pisc_biomass == 0),
    "/", nrow(pisc_transect_data),
    "(", round(mean(pisc_transect_data$transect_pisc_biomass == 0),
               3), ")\n")


# ============================================================
#  MODEL FAMILY SELECTION
#  Run on baseline model. Zeros present — Tweedie required.
#  ZI Tweedie tested given zero proportion.
# ============================================================
# ── Gaussian log — baseline ───────────────────────────────────
pisc_lm_baseline <- lm(
  log(mean_biomass + 0.01) ~ rugosity_sc +
    log_chla_sc,
  data = pisc_model_data
)

par(mfrow = c(2, 2))
plot(pisc_lm_baseline, main = "Gaussian log — baseline")
par(mfrow = c(1, 1))

# ── Tweedie — baseline ───────────────────────────────────────
pisc_tw_base <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data
)

pisc_tw_res <- simulateResiduals(pisc_tw_base, n = 1000)
plot(pisc_tw_res)
testZeroInflation(pisc_tw_res)
testDispersion(pisc_tw_res)

# ── Family selection decision ─────────────────────────────────
#
# Gaussian log: REJECTED
#   Zero sites (n = ~6) pull residuals to -10 at low
#   fitted values — severe lower tail Q-Q deviation
#   (sites 4, 41, 5). Scale-location strong downward
#   trend. Sites 9, 42 approach Cook's distance
#   threshold. Identical problem to browsers and
#   excavators — bimodal log distribution from zero
#   sites cannot be resolved with an offset constant.
#
# Tweedie (log link): SELECTED
#   DHARMa diagnostics (n = 1000):
#     KS test:        p = 0.987 — excellent fit
#     Dispersion:     p = 0.600, ratio = 1.117 — acceptable
#     Zero inflation: p = 1.000, ratio = 0.970 — not needed
#     Outlier test:   p = 1.000 — no outliers
#   Residuals vs predicted: slight upward trend in
#     smoothing line but well within confidence band.
#     No significant problems detected.
#
#  Proceed: glmmTMB(family = tweedie(link = "log")) on
#  raw mean_biomass throughout all piscivore analyses.


# ============================================================
#  RANDOM EFFECT STRUCTURE
#  Tested on baseline model. Tweedie family throughout.
# ============================================================

pisc_re_null <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data
)

pisc_re_ecoregion <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | ecoregion),
  family = tweedie(link = "log"),
  data   = pisc_model_data
)

cat("\n--- Piscivore RE structure: ecoregion ---\n")
print(make_aicc_df(list(
  "No RE"           = pisc_re_null,
  "(1 | ecoregion)" = pisc_re_ecoregion
)))

# Random effect structure: ecoregion
# Tested on baseline model to avoid circularity.
# No RE:           AICc = 800.68, weight = 0.802 (BEST)
# (1 | ecoregion): DAICc = 2.80,  weight = 0.198
# Ecoregion RE not supported — consistent with total
# biomass (DAICc = 2.25), corallivores (DAICc = 2.54),
# excavators (DAICc = 2.54), grazer-detritivores
# (DAICc = 2.93). Remarkably consistent pattern across
# all functional groups.
# Not pursued — only 4 ecoregions with severely uneven
# group sizes (n = 2, 8, 9, 35; Gelman & Hill 2007).
# All piscivore models fitted without RE throughout.

# ── Variance inflation factors ────────────────────────────────
# Sequence: baseline → baseline + market gravity →
#           connectivity x market gravity →
#           connectivity x market gravity + MPA

cat("\n--- VIF: baseline ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data))

cat("\n--- VIF: baseline + market gravity ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data))

cat("\n--- VIF: baseline + connectivity x market gravity ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + connectivity_sc * log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data))

cat("\n--- VIF: baseline + connectivity x market gravity + MPA ---\n")
check_collinearity(glmmTMB(
  mean_biomass ~ rugosity_sc + log_chla_sc + connectivity_sc * log_market_gravity_sc + mpa_status,
  family = tweedie(link = "log"),
  data   = pisc_model_data))

# ============================================================
#  Q1 — HUMAN PRESSURE
#
#  Does human pressure explain variation in piscivore
#  biomass beyond local ecological context?
#  Which metric best captures SSF exploitation intensity?
#
#  Prior analysis found market gravity preferred —
#  ecologically motivated as piscivores are large-bodied
#  commercial species targeted by market-oriented fishing.
#  Confirm with current data.
#
#  If pressure supported: carry best metric forward.
#  If baseline best: no pressure carried forward.
# ============================================================

p_baseline <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data
)

p_q1_settgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data
)

p_q1_mktgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data
)

p_q1_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_pop_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data
)

cat("\n--- Q1 Step 1: Piscivore metric comparison ---\n")
print(make_aicc_df(list(
  "Baseline"                      = p_baseline,
  "Baseline + settlement gravity" = p_q1_settgrav,
  "Baseline + market gravity"     = p_q1_mktgrav,
  "Baseline + settlement pop."    = p_q1_settpop
)))

cat("\n--- Q1 Step 2: Baseline coefficients ---\n")
summary(p_baseline)

summary(p_q1_mktgrav)

# ── Q1: Market gravity range and fold difference ─────────────
cat("\n--- Piscivore Q1: market gravity range ---\n")
mkt_range_p <- range(pisc_model_data$log_market_gravity_sc,
                     na.rm = TRUE)
cat("Range:", mkt_range_p, "\n")
mkt_span_p  <- diff(mkt_range_p)
cat(sprintf("Span: %.3f SD units\n", mkt_span_p))

cat("\n--- Piscivore Q1: market gravity model confint ---\n")
print(confint(p_q1_mktgrav))

cat("\n--- Q1: Pressure metric direction checks ---\n")
cat("Settlement gravity:\n")
print(summary(p_q1_settgrav)$coefficients$cond[
  "log_settlement_grav_sc", ])
cat("\nMarket gravity:\n")
print(summary(p_q1_mktgrav)$coefficients$cond[
  "log_market_gravity_sc", ])
cat("\nSettlement population:\n")
print(summary(p_q1_settpop)$coefficients$cond[
  "log_settlement_pop_sc", ])

# Q1 results:
#   Market gravity:     AICc = 799.72, weight = 0.488 (BEST)
#   Baseline:           DAICc = 0.96,  weight = 0.302
#   Settlement pop.:    DAICc = 2.96,  weight = 0.111
#   Settlement gravity: DAICc = 3.19,  weight = 0.099
#
#   Market gravity preferred — weight = 0.488, nearly
#   twice the baseline weight. However genuine model
#   selection uncertainty (DAICc = 0.96 vs baseline) —
#   not decisive support. Contrast with total biomass
#   where settlement gravity was unambiguous (weight =
#   0.826, DAICc = 4.39 vs baseline).
#
# Baseline coefficients:
#   Rugosity: b = +0.015, p = 0.897 ns
#   Chla:     b = -0.130, p = 0.409 ns
#   Neither significant — piscivore biomass poorly
#   explained by ecological baseline alone. High
#   dispersion (18.1) reflects patchy distribution
#   of large predators.
#
# Pressure metric direction checks:
#   Market gravity:     b = +0.286, p = 0.058 . marginal
#     Positive direction — counterintuitive for
#     exploitation pressure. Accessible high-market
#     sites may be historically productive or coincide
#     with urban centres where enforcement is higher.
#     Positive direction likely conditional on
#     connectivity — see Q2 interaction.
#   Settlement gravity: b = -0.102, p = 0.578 ns
#     Negative but not significant.
#   Settlement population: b = +0.108, p = 0.460 ns
#     Positive, not significant.
#
#   Market gravity carried forward as best supported
#   metric despite marginal support — ecologically
#   motivated (piscivores are commercially targeted)
#   and consistent with prior analysis result.
#
# Best Q1 model: p_q1_mktgrav
# (rugosity + chla + market gravity)
p_best_q1 <- p_q1_mktgrav

# ============================================================
#  Q2 — LARVAL CONNECTIVITY
#
#  Does connectivity explain additional variation beyond
#  the best Q1 model?
#
#  Pressure supported in Q1 (market gravity, weight = 0.488)
#  — connectivity x pressure interaction testable.
#
#  Two steps:
#  Step 1 — connectivity as main effect
#  Step 2 — connectivity x pressure interaction
#  Tested against best Q1 throughout.
#
#  Tweedie family throughout.
# ============================================================

# ── Step 1: Connectivity main effect ─────────────────────────
p_q2_conn <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_market_gravity_sc +
    connectivity_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data
)

cat("\n--- Q2 Step 1: Piscivore connectivity main effect ---\n")
print(make_aicc_df(list(
  "Best Q1"        = p_best_q1,
  "Best Q1 + conn" = p_q2_conn
)))

cat("\n--- Q2 Step 1: Connectivity coefficients ---\n")
summary(p_q2_conn)

# Q2 Step 1 results:
#   Best Q1:        AICc = 799.72, weight = 0.549
#   Best Q1 + conn: DAICc = 0.40,  weight = 0.451
#   Genuine uncertainty — neither clearly preferred.
#   Connectivity: b = +0.214, p = 0.130 ns


# ── Step 2: Connectivity x pressure interaction ───────────────
# A priori hypothesis: connectivity moderates the
# market gravity-biomass relationship — well-connected
# sites with high market access experience stronger
# exploitation pressure through broader fishing networks.
# Tested against best Q1, not the main effect model,
# since main effect not clearly supported in Step 1.

p_q2_conn_int <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc * log_market_gravity_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data
)

cat("\n--- Q2 Step 2: Connectivity x pressure interaction ---\n")
print(make_aicc_df(list(
  "Best Q1"                   = p_best_q1,
  "Best Q1 + conn"            = p_q2_conn,
  "Best Q1 + conn x pressure" = p_q2_conn_int
)))

cat("\n--- Q2 Step 2: Interaction coefficients ---\n")
summary(p_q2_conn_int)

# ── Q2: Interaction confint ───────────────────────────────────
cat("\n--- Piscivore Q2: interaction model confint ---\n")
print(confint(p_q2_conn_int))

# ── Q2: Connectivity range ────────────────────────────────────
cat("\n--- Piscivore Q2: connectivity range ---\n")
conn_range_p <- range(pisc_model_data$connectivity_sc,
                      na.rm = TRUE)
cat("Range:", conn_range_p, "\n")
cat(sprintf("Span: %.3f SD units\n", diff(conn_range_p)))

# Q2 Step 2 results:
#   Conn x pressure: AICc = 795.27, weight = 0.836 (BEST)
#   Best Q1:         DAICc = 4.45,  weight = 0.090
#   Best Q1 + conn:  DAICc = 4.85,  weight = 0.074
#   Interaction strongly supported — resolves Q1 and
#   Q2 main effect uncertainty.
#
#   Conn x market gravity: b = -0.511, p = 0.004 **
#     Negative — connectivity moderates market gravity.
#     At low connectivity, market access associates
#     positively with piscivore biomass. At high
#     connectivity, market access depletes biomass —
#     consistent with connectivity facilitating
#     commercial fishing access to well-connected reefs.
#   Market gravity: b = +0.407, p = 0.008 **
#     Positive at mean connectivity — resolves
#     counterintuitive Q1 direction.
#   Connectivity:   b = +0.265, p = 0.048 *
#     Positive at mean market gravity.
#   Rugosity:       b = +0.240, p = 0.067 .
#     Marginal positive — consistent throughout.

# ── Best Q2 model ─────────────────────────────────────────────
# Connectivity x pressure interaction strongly supported.
p_best_q2 <- p_q2_conn_int

# ============================================================
#  Q3 — FORMAL PROTECTION
#
#  Does MPA status explain additional variation beyond
#  the connectivity x pressure model?
#
#  MPA x connectivity interaction explored if MPA
#  supported as main effect (DAICc > 2).
#
#  Tweedie family throughout.
# ============================================================

p_q3_mpa <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc * log_market_gravity_sc +
    mpa_status,
  family = tweedie(link = "log"),
  data   = pisc_model_data
)

cat("\n--- Q3: Piscivore MPA main effect ---\n")
print(make_aicc_df(list(
  "Best Q2"       = p_best_q2,
  "Best Q2 + MPA" = p_q3_mpa
)))

cat("\n--- Q3: MPA coefficients ---\n")
summary(p_q3_mpa)

# ── Q3: MPA confint and fold difference ──────────────────────
cat("\n--- Piscivore Q3: MPA model confint ---\n")
print(confint(p_q3_mpa))

cat("\n--- Piscivore Q3: medium MPA fold difference ---\n")
b_medium_mpa_p <- 0.868
fold_mpa_p     <- exp(b_medium_mpa_p)
cat(sprintf("Medium MPA fold difference: %.2fx\n", fold_mpa_p))

# ── Q3: MPA placement check ──────────────────────────────────
cat("\n--- Piscivore Q3: pressure and connectivity by MPA ---\n")
pisc_model_data %>%
  group_by(mpa_status) %>%
  summarise(
    n                    = n(),
    mean_biomass         = round(mean(mean_biomass), 1),
    mean_settlement_grav = round(mean(log_settlement_grav_sc),
                                 3),
    mean_market_grav     = round(mean(log_market_gravity_sc),
                                 3),
    mean_connectivity    = round(mean(connectivity_sc), 3),
    .groups = "drop"
  ) %>%
  print()

# Q3 results:
#   Best Q2 + MPA: AICc = 792.22, weight = 0.821 (BEST)
#   Best Q2:       DAICc = 3.05,  weight = 0.179
#   MPA strongly supported beyond connectivity x
#   pressure model.
#
#   Medium MPA: b = +0.868, p = 0.003 **
#     Highly significant — medium protection sites
#     have substantially higher piscivore biomass
#     at equivalent habitat, pressure, and connectivity
#     conditions. Consistent with piscivores being
#     heavily targeted by fishing and benefiting
#     directly from harvest exclusion.
#   Low MPA:    b = +0.197, p = 0.649 ns
#     Not significant — low protection insufficient
#     to benefit large predators.
#   Rugosity:   b = +0.299, p = 0.015 * — now
#     significant once MPA and interaction included.
#     Habitat complexity is a genuine driver of
#     piscivore biomass when other factors controlled.
#   Interaction: b = -0.578, p = 0.001 ***
#     Stable and strengthened (was -0.511, p = 0.004).
#     Connectivity x pressure dynamic robust to
#     inclusion of MPA.
#
#   This contrasts with total biomass (MPA not
#   supported, DAICc = 3.01) — piscivores are
#   demonstrably more sensitive to formal protection
#   than the total fish community, consistent with
#   their status as the most heavily targeted group.

# ── MPA data structure check ──────────────────────────────────
cat("\n--- Q3: MPA site distribution ---\n")
pisc_model_data %>%
  group_by(mpa_status) %>%
  summarise(
    n                 = n(),
    mean_biomass      = round(mean(mean_biomass), 1),
    mean_connectivity = round(mean(connectivity_sc), 3),
    min_connectivity  = round(min(connectivity_sc),  3),
    max_connectivity  = round(max(connectivity_sc),  3),
    .groups = "drop"
  ) %>%
  print()

# ── Q3 MPA × connectivity interaction ────────────────────────
# Run only if MPA supported as main effect (DAICc > 2).

p_q3_mpa_conn_int <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_market_gravity_sc +  
    mpa_status * connectivity_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data
)

cat("\n--- Q3: MPA x connectivity interaction ---\n")
print(make_aicc_df(list(
  "Best Q2"              = p_best_q2,
  "Best Q2 + MPA"        = p_q3_mpa,
  "Best Q2 + MPA x conn" = p_q3_mpa_conn_int
)))

cat("\n--- Q3: MPA x connectivity coefficients ---\n")
summary(p_q3_mpa_conn_int)

# Q3 MPA x connectivity interaction:
#   Best Q2 + MPA:         AICc = 792.22, weight = 0.815
#   Best Q2:               DAICc = 3.05,  weight = 0.178
#   Best Q2 + MPA x conn:  DAICc = 9.49,  weight = 0.007
#
#   Interaction clearly not supported (DAICc = 9.49) —
#   substantially worse than MPA main effect model.
#   MPA effectiveness does not vary with connectivity
#   for piscivores.
#   Note: interaction model does not include the
#   connectivity x pressure term from best Q2 —
#   coefficients not directly comparable. MPA main
#   effect model (p_best_q3) is the correct reference.
#
# Final model sequence:
#   Q1: market gravity preferred (weight = 0.488,
#       DAICc = 0.96 vs baseline) — marginal but
#       ecologically motivated
#   Q2: connectivity x pressure interaction strongly
#       supported (weight = 0.836, DAICc = 4.45)
#   Q3: MPA main effect supported (weight = 0.821,
#       DAICc = 3.05); medium MPA b = +0.868, p = 0.003
#   MPA x connectivity not supported (DAICc = 9.49)
#
# Best Q3 model: rugosity + chla +
#   connectivity x market gravity + MPA status
p_best_q3 <- p_q3_mpa

cat("\n--- Piscivore: predicted vs observed ---\n")
pred_p <- predict(p_best_q3, type = "response")
obs_p  <- pisc_model_data$mean_biomass
cat(sprintf("Pearson r: %.3f\n", cor(pred_p, obs_p)))

cat("\n--- Best model: DHARMa diagnostics ---\n")
pisc_best_sim <- simulateResiduals(p_best_q3, n = 1000)
testZeroInflation(pisc_best_sim)
testDispersion(pisc_best_sim)
testOutliers(pisc_best_sim)

# ── Q3 exploratory checks ─────────────────
cat("\n--- Q3: Biomass by MPA status ---\n")
pisc_model_data %>%
  group_by(mpa_status) %>%
  summarise(
    n            = n(),
    mean_biomass = round(mean(mean_biomass), 1),
    sd_biomass   = round(sd(mean_biomass),   1),
    mean_conn    = round(mean(connectivity_sc), 3),
    .groups = "drop"
  ) %>%
  print()

pisc_mpa_sim <- simulateResiduals(p_q3_mpa, n = 1000)
plot(pisc_mpa_sim)
testOutliers(pisc_mpa_sim)

ggplot(pisc_model_data,
       aes(x = mpa_status, y = mean_biomass,
           fill = mpa_status)) +
  geom_boxplot(outlier.shape = 16, outlier.size = 1.5,
               alpha = 0.7, width = 0.5) +
  geom_jitter(width = 0.1, size = 1.8,
              alpha = 0.6, colour = "grey30") +
  scale_fill_manual(values = c("none"   = "#bdbdbd",
                               "low"    = "#74a9cf",
                               "medium" = "#0570b0")) +
  scale_x_discrete(labels = c("none"   = "No MPA",
                               "low"    = "Low",
                               "medium" = "Medium")) +
  labs(x = "MPA status", y = "Piscivore biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title         = element_text(face = "bold"),
        legend.position    = "none",
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank())

# ── Q3 exploratory checks ────────────────────────────────────
#
# Biomass by MPA status:
#   none:   n = 30, mean = 579g,  conn = -0.345
#   low:    n = 7,  mean = 754g,  conn = +0.879
#   medium: n = 17, mean = 1284g, conn = +0.303
#
#   Clear biomass gradient across protection levels —
#   medium MPA sites have 2.2x higher mean biomass than
#   unprotected sites. Low MPA sites also elevated
#   (754g vs 579g) but not significant in model —
#   consistent with low protection being insufficient
#   for large predators.
#   Note: low MPA sites have highest mean connectivity
#   (z = +0.879) — same pattern as other functional
#   groups. Low MPA coefficient not artefactual here
#   since it is not significant (p = 0.649) and the
#   raw biomass advantage is modest (754g vs 579g).
#   Medium MPA advantage genuine and ecologically
#   interpretable — harvest exclusion benefits large
#   predators more than any other functional group
#   except browsers.
#
# DHARMa diagnostics (n = 1000):
#   Outlier test: p = 1.000 — no outliers.
#   (Full KS and dispersion from plot — update if
#   significant problems detected.)
#
# MPA x connectivity: DAICc = 9.49, not supported.
# MPA main effect model retained as best Q3.

# ============================================================
#  SPATIAL AUTOCORRELATION
# ============================================================

site_coords <- location_data %>%
  mutate(site = as.character(site)) %>%
  group_by(site) %>%
  summarise(lon = first(longitude),
            lat = first(latitude),
            .groups = "drop")

pisc_model_data_coords <- pisc_model_data %>%
  left_join(site_coords, by = "site")

coords_mat_p <- cbind(pisc_model_data_coords$lon,
                      pisc_model_data_coords$lat)
listw5_p <- nb2listw(knn2nb(knearneigh(coords_mat_p, k = 5)),
                     style = "W")

cat("\n--- Spatial autocorrelation: piscivore best model ---\n")
print(moran.test(residuals(p_best_q3, type = "pearson"),
                 listw5_p))


# Spatial autocorrelation: piscivore best model
# (rugosity + chla + conn x market gravity + MPA)
# Moran's I = -0.043, p = 0.630 — no significant
# spatial autocorrelation in residuals.
#
# Clean — consistent with corallivores (I = -0.021)
# and excavators (I = -0.015). Contrasts with total
# biomass (I = 0.140, p = 0.015) and grazer-
# detritivores (I = 0.210, p = 0.001).
# The connectivity x pressure interaction and MPA
# adequately capture spatial variation in piscivore
# biomass without leaving a residual geographic signal.
# No spatial error modelling required.

# ============================================================
#  SENSITIVITY ANALYSIS
# ============================================================

# ── (a) Alternative pressure metrics ─────────────────────────
p_sens_settgrav <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_grav_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data
)

p_sens_settpop <- glmmTMB(
  mean_biomass ~ rugosity_sc +
    log_chla_sc +
    log_settlement_pop_sc,
  family = tweedie(link = "log"),
  data   = pisc_model_data
)

cat("\n--- Sensitivity (a): piscivore alternative metrics ---\n")
cat("Settlement gravity:\n")
print(summary(p_sens_settgrav)$coefficients$cond)
cat("\nSettlement population:\n")
print(summary(p_sens_settpop)$coefficients$cond)

# Sensitivity (a) results:
#   Settlement gravity:    b = -0.102, p = 0.578 ns
#   Settlement population: b = +0.108, p = 0.460 ns
#
#   Both non-significant — consistent with Q1 where
#   neither outperformed the baseline. Directions
#   inconsistent (settlement gravity negative,
#   settlement population positive) — no coherent
#   pressure signal with alternative metrics.
#   Market gravity remains the only metric showing
#   a meaningful signal for piscivores (Q1 weight =
#   0.488, marginal positive in Q2 interaction).
#   Q1 metric choice robust — alternative metrics
#   do not detect the pressure signal that market
#   gravity captures through the interaction.


# ── (b) Transect-level replication ───────────────────────────

# ── Family selection ──────────────────────────────────────────
p_trans_tw_base <- glmmTMB(
  transect_pisc_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = pisc_transect_data
)

p_trans_tw_zi_base <- glmmTMB(
  transect_pisc_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | site),
  family    = tweedie(link = "log"),
  ziformula = ~1,
  data      = pisc_transect_data
)

p_trans_res    <- simulateResiduals(p_trans_tw_base,    n = 500)
p_trans_res_zi <- simulateResiduals(p_trans_tw_zi_base, n = 500)

plot(p_trans_res);    testZeroInflation(p_trans_res)
plot(p_trans_res_zi); testZeroInflation(p_trans_res_zi)

cat("\n--- Sensitivity (b): transect family selection ---\n")
print(make_aicc_df(list(
  "Tweedie"    = p_trans_tw_base,
  "ZI Tweedie" = p_trans_tw_zi_base
)))

# Tweedie:    AICc = 2557.28, weight = 0.743 — SELECTED
# ZI Tweedie: DAICc = 2.12,   weight = 0.257
# ZI not supported — consistent with site-level.

# ── Transect models mirroring Q1-Q3 sequence ─────────────────
p_trans_null <- glmmTMB(
  transect_pisc_biomass ~ 1 + (1 | site),
  family = tweedie(link = "log"),
  data   = pisc_transect_data
)

p_trans_baseline <- glmmTMB(
  transect_pisc_biomass ~ rugosity_sc +
    log_chla_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = pisc_transect_data
)

# Mirrors Q1 — pressure baseline
p_trans_pressure <- glmmTMB(
  transect_pisc_biomass ~ rugosity_sc +
    log_chla_sc +
    log_market_gravity_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = pisc_transect_data
)

# Mirrors Q2 Step 1 — connectivity main effect
# Note: pressure included since Q1 supported
p_trans_conn <- glmmTMB(
  transect_pisc_biomass ~ rugosity_sc +
    log_chla_sc +
    log_market_gravity_sc +
    connectivity_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = pisc_transect_data
)

# Mirrors Q2 Step 2 — connectivity x pressure interaction
p_trans_conn_int <- glmmTMB(
  transect_pisc_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc * log_market_gravity_sc +
    (1 | site),
  family = tweedie(link = "log"),
  data   = pisc_transect_data
)

# Mirrors Q3 — interaction + MPA
p_trans_mpa_int <- glmmTMB(
  transect_pisc_biomass ~ rugosity_sc +
    log_chla_sc +
    connectivity_sc * log_market_gravity_sc +
    mpa_status +
    (1 | site),
  family = tweedie(link = "log"),
  data   = pisc_transect_data
)

cat("\n--- Sensitivity (b): piscivore transect comparison ---\n")
print(make_aicc_df(list(
  "Null"               = p_trans_null,
  "Baseline"           = p_trans_baseline,
  "Baseline + press"   = p_trans_pressure,
  "Baseline + conn"    = p_trans_conn,
  "Conn x press"       = p_trans_conn_int,
  "Conn x press + MPA" = p_trans_mpa_int
)))

cat("\n--- Sensitivity (b): interaction coefficients ---\n")
summary(p_trans_conn_int)

# ── ICC calculation ───────────────────────────────────────────
vc_p <- VarCorr(p_trans_baseline)
site_sd_p <- sqrt(as.numeric(vc_p$cond$site))
cat(sprintf("\nSite random intercept SD = %.3f\n", site_sd_p))

# Sensitivity (b) results:
#
# AICc comparison (n = 243 transects, 54 sites):
#   Conn x press + MPA: AICc = 2546.70, weight = 0.893
#   Conn x press:       DAICc = 4.92,   weight = 0.076
#   All other models:   DAICc > 9 — not competitive
#
#   Conn x press + MPA best supported at transect level
#   — stronger than site level (site weight = 0.821).
#   Both primary Q2 and Q3 findings fully replicated.
#
# Interaction coefficients (REML):
#   Conn x market gravity: b = -0.529, p = 0.007 **
#     Consistent with site-level (b = -0.511, p = 0.004)
#     — direction and magnitude stable across scales.
#   Market gravity:        b = +0.429, p = 0.012 *
#     Stable (site: b = +0.407, p = 0.008).
#   Connectivity:          b = +0.292, p = 0.047 *
#     Stable (site: b = +0.265, p = 0.048).
#   Rugosity:              b = +0.246, t = 1.65 .
#     Direction consistent — marginal at both scales.
#   Site variance: 0.417 (SD = 0.646)
#   Dispersion:    59 — high within-site variance,
#     consistent with patchily distributed predators.
#
# ICC = 0.000 — between-site variance effectively zero
#   in the baseline model. This is because piscivore
#   biomass variability is dominated by within-site
#   transect-to-transect variation (dispersion = 58.9)
#   rather than between-site differences. Once the
#   interaction structure is included, between-site
#   clustering increases (site variance = 0.417 in
#   interaction model). Confirms that fixed effects
#   (not random structure) drive the piscivore signal.
#
# Overall: primary findings fully replicated at
#   transect level and with stronger support.
#   Interaction coefficient stable in direction and
#   magnitude across both analytical scales — not
#   an aggregation artefact.

# ============================================================
#  MARGINAL EFFECT PLOTS
#  Best model: rugosity + chla + connectivity x market
#  gravity + MPA status (p_best_q3)
#  Three plots:
#  (1) Interaction surface — connectivity x market gravity
#  (2) MPA marginal means
#  (3) Rugosity effect (marginal, p = 0.067 at site level)
#  Predictions on response scale (raw biomass).
# ============================================================

best_model_p <- p_best_q3  # conn x market gravity + MPA

# ── Interaction plot: connectivity x market gravity ───────────
# Primary result — show at low/medium/high connectivity.
# MPA held at "none", rugosity and chla at mean (z = 0).

conn_quantiles_p <- quantile(pisc_model_data$connectivity_sc,
                             c(0.10, 0.50, 0.90))

pred_grid_p <- expand.grid(
  log_market_gravity_sc = seq(
    min(pisc_model_data$log_market_gravity_sc),
    max(pisc_model_data$log_market_gravity_sc),
    length.out = 200),
  connectivity_sc = conn_quantiles_p,
  rugosity_sc     = 0,
  log_chla_sc     = 0,
  mpa_status      = factor("none",
                           levels = c("none", "low", "medium"))
)

pred_grid_p$fit <- predict(best_model_p,
                           newdata = pred_grid_p,
                           type    = "response",
                           re.form = NA)

pred_grid_p$conn_label <- factor(
  round(pred_grid_p$connectivity_sc, 2),
  labels = c("Low connectivity (10th percentile)",
             "Medium connectivity (50th percentile)",
             "High connectivity (90th percentile)"))

p_p_interaction <- ggplot(pred_grid_p,
                          aes(x      = log_market_gravity_sc,
                              y      = fit,
                              colour = conn_label,
                              group  = conn_label)) +
  geom_line(linewidth = 1.1) +
  geom_rug(data = pisc_model_data,
           aes(x = log_market_gravity_sc),
           inherit.aes = FALSE,
           alpha = 0.4, sides = "b") +
  scale_colour_manual(
    values = c("Low connectivity (10th percentile)"    = "#d7191c",
               "Medium connectivity (50th percentile)" = "#fdae61",
               "High connectivity (90th percentile)"   = "#2c7bb6")
  ) +
  labs(x      = "log(Market gravity) (standardised)",
       y      = "Piscivore biomass (g)",
       colour = NULL) +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        legend.position  = "top",
        panel.grid.minor = element_blank())


# ── MPA marginal means ────────────────────────────────────────
# MPA supported in Q3 (weight = 0.821).
# Medium protection: substantial biomass advantage.
# All continuous predictors held at mean (z = 0).

mpa_grid_p <- data.frame(
  mpa_status            = factor(c("none", "low", "medium"),
                                 levels = c("none", "low",
                                            "medium")),
  rugosity_sc           = 0,
  log_chla_sc           = 0,
  log_market_gravity_sc = 0,
  connectivity_sc       = 0
)

mpa_pred_p     <- predict(best_model_p,
                          newdata = mpa_grid_p,
                          se.fit  = TRUE,
                          type    = "response",
                          re.form = NA)
mpa_grid_p$fit <- mpa_pred_p$fit
mpa_grid_p$lwr <- mpa_pred_p$fit - 1.96 * mpa_pred_p$se.fit
mpa_grid_p$upr <- mpa_pred_p$fit + 1.96 * mpa_pred_p$se.fit

p_p_mpa <- ggplot(mpa_grid_p,
                  aes(x = mpa_status, y = fit)) +
  geom_hline(yintercept = mpa_grid_p$fit[1],
             linetype   = "dashed",
             colour     = "grey70",
             linewidth  = 0.4) +
  geom_pointrange(aes(ymin = lwr, ymax = upr),
                  colour    = "#0570b0",
                  linewidth = 0.7,
                  size      = 0.6) +
  scale_x_discrete(labels = c("none"   = "No MPA",
                              "low"    = "Low",
                              "medium" = "Medium")) +
  labs(x = "MPA status", y = "Piscivore biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title         = element_text(face = "bold"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank())


# ── Rugosity effect ───────────────────────────────────────────
# Marginal at site level (p = 0.067) — significant at
# transect level (t = 1.65). Included for completeness.
# All other predictors held at mean (z = 0), MPA = none.

rug_grid_p <- data.frame(
  rugosity_sc           = seq(
    min(pisc_model_data$rugosity_sc),
    max(pisc_model_data$rugosity_sc),
    length.out = 200),
  log_chla_sc           = 0,
  log_market_gravity_sc = 0,
  connectivity_sc       = 0,
  mpa_status            = factor("none",
                                 levels = c("none", "low",
                                            "medium"))
)

rug_pred_p     <- predict(best_model_p,
                          newdata = rug_grid_p,
                          se.fit  = TRUE,
                          type    = "response",
                          re.form = NA)
rug_grid_p$fit <- rug_pred_p$fit
rug_grid_p$lwr <- rug_pred_p$fit - 1.96 * rug_pred_p$se.fit
rug_grid_p$upr <- rug_pred_p$fit + 1.96 * rug_pred_p$se.fit

p_p_rugosity <- ggplot(rug_grid_p,
                       aes(x = rugosity_sc, y = fit)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr),
              fill = "#2c7bb6", alpha = 0.15) +
  geom_line(colour = "#2c7bb6", linewidth = 1.1) +
  geom_point(data = pisc_model_data,
             aes(x = rugosity_sc, y = mean_biomass),
             colour = "grey40", size = 1.5,
             alpha  = 0.5, inherit.aes = FALSE) +
  labs(x = "Rugosity (standardised)",
       y = "Piscivore biomass (g)") +
  theme_bw(base_size = 12) +
  theme(axis.title       = element_text(face = "bold"),
        panel.grid.minor = element_blank())


# ── Arrange plots ─────────────────────────────────────────────
gridExtra::grid.arrange(p_p_interaction, p_p_mpa,
                        p_p_rugosity,
                        ncol = 3)

# jpeg("piscivore_marginal_effects.jpg",
#      width = 33, height = 11, units = "cm", res = 300)
# gridExtra::grid.arrange(p_p_interaction, p_p_mpa,
#                         p_p_rugosity, ncol = 3)
# dev.off()


# ============================================================
#  RESULTS SUMMARY
# ============================================================

cat("\n--- Piscivore results summary ---\n")
tribble(
  ~Question,   ~Result,              ~Key_finding,
  "Q1",        "Market gravity",     "weight = 0.488, DAICc = 0.96 vs baseline — marginal, ecologically motivated",
  "Q2 conn",   "Interaction",        "weight = 0.836, DAICc = 4.45, conn x mkt gravity b = -0.511, p = 0.004",
  "Q3 MPA",    "Supported",          "weight = 0.821, DAICc = 3.05, medium b = +0.868, p = 0.003",
  "Q3 int",    "Not supported",      "MPA x conn DAICc = 9.49",
  "Spatial",   "Clean",              "Moran's I = -0.043, p = 0.630",
  "Sens (a)",  "Consistent",         "alternative metrics ns — market gravity signal metric-specific",
  "Sens (b)",  "Consistent",         "conn x press + MPA weight = 0.893, interaction b = -0.529, p = 0.007"
) %>% print()


# ============================================================
#  SESSION INFO
# ============================================================
cat("\n--- Session info ---\n")
sessionInfo()