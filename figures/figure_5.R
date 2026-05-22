library(ggplot2)
library(dplyr)
library(tibble)

font <- "serif"
border_col <- "grey75"

mpa_cols <- c(
  "none"   = "#8F8F8F",
  "low"    = "#58C3A6",
  "medium" = "#118C7A"
)

mpa_levels <- c("none", "low", "medium")

# ------------------------------------------------------------
# Corrected results with proper 95% CIs
# ------------------------------------------------------------
plot_df <- tribble(
  ~outcome,                  ~mpa_status, ~mean,   ~lwr,    ~upr,
  
  "Settlement gravity",      "none",       0.242,  -0.165,  0.650,
  "Settlement gravity",      "low",        0.473,  -0.077,  1.023,
  "Settlement gravity",      "medium",    -0.532,  -0.725, -0.338,
  
  "Settlement population",   "none",       0.220,  -0.087,  0.528,
  "Settlement population",   "low",        0.747,   0.549,  0.945,
  "Settlement population",   "medium",    -0.666,  -1.180, -0.149,
  
  "Connectivity",            "none",      -0.345,  -0.687, -0.002,
  "Connectivity",            "low",        0.879,   0.769,  0.988,
  "Connectivity",            "medium",     0.303,  -0.170,  0.775
)

plot_df$mpa_status <- factor(plot_df$mpa_status,
                             levels = mpa_levels)

plot_df$outcome <- factor(
  plot_df$outcome,
  levels = c(
    "Settlement gravity",
    "Settlement population",
    "Connectivity"
  )
)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
fig_mpa <- ggplot(
  plot_df,
  aes(x = mpa_status,
      y = mean,
      color = mpa_status)
) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.5,
    linetype = "dashed",
    color = "grey80"
  ) +
  
  geom_errorbar(
    aes(ymin = lwr, ymax = upr),
    width = 0.12,
    linewidth = 0.9
  ) +
  
  geom_point(
    shape = 16,
    size = 3.5
  ) +
  
  facet_wrap(
    ~ outcome,
    nrow = 1
  ) +
  
  scale_color_manual(values = mpa_cols) +
  
  scale_x_discrete(
    labels = c(
      "none"   = "None",
      "low"    = "Low",
      "medium" = "Medium"
    )
  ) +
  
  scale_y_continuous(
    limits = c(-1.5, 1.3),
    breaks = c(-1.0, -0.5, 0, 0.5, 1.0),
    expand = expansion(mult = c(0.02, 0.03))
  ) +
  
  labs(
    x = "Protection level",
    y = "Standardised\nvalue (z-score)"
  ) +
  
  theme_classic(base_family = font) +
  
  theme(
    legend.position = "none",
    
    axis.title.x = element_text(
      size = 13,
      face = "bold",
      margin = margin(t = 8)
    ),
    
    axis.title.y = element_text(
      size = 13,
      face = "bold",
      angle = 0,
      lineheight = 0.95,
      hjust = 0.5,
      vjust = 0.5,
      margin = margin(r = 10)
    ),
    
    axis.text.x = element_text(
      size = 12,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 11,
      colour = "black"
    ),
    
    axis.line = element_line(
      colour = border_col,
      linewidth = 0.45
    ),
    
    axis.ticks = element_line(
      colour = border_col,
      linewidth = 0.45
    ),
    
    strip.background = element_rect(
      fill = "grey92",
      colour = border_col,
      linewidth = 0.4
    ),
    
    strip.text = element_text(
      size = 12,
      face = "bold",
      family = font
    ),
    
    panel.border = element_rect(
      color = border_col,
      fill = NA,
      linewidth = 0.4
    ),
    
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    panel.spacing.x = unit(0.18, "lines"),
    
    aspect.ratio = 1.15,
    
    plot.margin = margin(
      t = 5,
      r = 75,
      b = 5,
      l = 5
    )
  )

fig_mpa

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
ggsave(
  "mpa_placement_final.pdf",
  fig_mpa,
  width = 8.0,
  height = 5.6,
  units = "in",
  bg = "white",
  useDingbats = FALSE
)

ggsave(
  "mpa_placement_final.png",
  fig_mpa,
  width = 8.0,
  height = 5.6,
  units = "in",
  dpi = 300,
  bg = "white"
)