# Sensitivity analyses 
# 
# 
# Generalised additive models (GAMs) were used to assess potential non-linear 
# relationships between predictors and biomass. Smooth terms indicated largely 
# linear relationships (EDF ≈ 1), and model fit was not improved relative to 
# linear models (ΔAIC < 2). Linear models were therefore retained for all 
# subsequent analyses.
# 
# For total biomass:
library(mgcv)

gam_global <- gam(
  log_mean_biomass ~
    s(rugosity_sc, k = 4) +
    s(log_settlement_grav_sc, k = 4) +
    mpa_status +
    s(connectivity_sc, k = 4) +
    s(log_chla_sc, k = 4) +
    s(log_max_dhw_sc, k = 4),
  data = total_model_data,
  method = "REML"
)

summary(gam_global)
gam.check(gam_global)
plot(gam_global, pages = 1)
AIC(m_global, gam_global)

# | Variable           | EDF  | Interpretation              |
# | ------------------ | ---- | --------------------------- |
# | rugosity           | 1.00 | linear                      |
# | settlement gravity | 1.54 | *slightly nonlinear (weak)* |
# | connectivity       | 1.34 | *very weak nonlinearity*    |
# | chl-a              | 1.00 | linear                      |
# | DHW                | 1.00 | linear                      |
# 
# #Almost everything is essentially linear
# Settlement + connectivity show minor curvature, but:
# not statistically significant
# not strong enough to improve fit
# 
library(MASS)

m_robust <- rlm(log_mean_biomass ~ rugosity_sc +
                  log_settlement_grav_sc +
                  mpa_status +
                  connectivity_sc +
                  log_chla_sc +
                  log_max_dhw_sc,
                data = total_model_data)
  