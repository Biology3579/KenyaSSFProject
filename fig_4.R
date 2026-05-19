library(ggplot2)
library(dplyr)
library(scales)

mod <- p_best_q2

grid <- expand.grid(
  connectivity_sc = seq(
    min(pisc_model_data$connectivity_sc, na.rm = TRUE),
    max(pisc_model_data$connectivity_sc, na.rm = TRUE),
    length.out = 800
  ),
  log_market_gravity_sc = c(-1, 1),
  rugosity_sc = 0,
  log_chla_sc = 0
)

pred <- predict(mod, newdata = grid, type = "link", se.fit = TRUE)

grid <- grid %>%
  mutate(
    fit    = exp(pred$fit),
    lo     = exp(pred$fit - 1.96 * pred$se.fit),
    hi     = exp(pred$fit + 1.96 * pred$se.fit),
    mg_lab = factor(
      log_market_gravity_sc,
      levels = c(-1, 1),
      labels = c("Far from markets", "Near markets")
    )
  )

mg_cols <- c(
  "Far from markets" = "#C97C3B",
  "Near markets"     = "#2F4858"
)

lab_dat <- grid %>%
  group_by(mg_lab) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  mutate(
    label = ifelse(
      mg_lab == "Far from markets",
      "far from markets\nmore fish",
      "near markets\nfewer fish"
    ),
    x = connectivity_sc + 0.03
  )

ggplot() +
  geom_ribbon(
    data = grid,
    aes(x = connectivity_sc, ymin = lo, ymax = hi, fill = mg_lab),
    alpha = 0.10,
    colour = NA
  ) +
  geom_point(
    data = pisc_model_data,
    aes(x = connectivity_sc, y = mean_biomass),
    colour = "grey50",
    alpha = 0.25,
    size = 1.7
  ) +
  geom_line(
    data = grid,
    aes(x = connectivity_sc, y = fit, colour = mg_lab),
    linewidth = 1.05,
    lineend = "round"
  ) +
  geom_text(
    data = lab_dat,
    aes(x = x, y = fit, label = label, colour = mg_lab),
    hjust = 0,
    size = 3.0,
    fontface = "italic",
    lineheight = 1.0,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = mg_cols, guide = "none") +
  scale_fill_manual(values = mg_cols, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  scale_y_log10(
    labels = comma,
    breaks = c(10, 50, 100, 500, 1000, 2000, 5000)
  ) +
  labs(
    x = "Larval connectivity (standardised)",
    y = "Mean piscivore biomass (g per site, log scale)"
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 11, base_family = "serif") +
  theme(
    axis.title      = element_text(face = "bold", size = 11, colour = "black"),
    axis.text       = element_text(size = 10, colour = "black"),
    axis.line       = element_line(colour = "black", linewidth = 0.5),
    axis.ticks      = element_line(colour = "black", linewidth = 0.5),
    legend.position = "none",
    plot.margin     = margin(5.5, 35, 5.5, 5.5)
  )