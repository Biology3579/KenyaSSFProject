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