library(ggplot2)
library(dplyr)
library(scales)
library(ggbeeswarm)

font       <- "serif"
border_col <- "grey75"

mpa_cols <- c(
  "None"   = "#8F8F8F",
  "Low"    = "#58C3A6",
  "Medium" = "#118C7A"
)

# ------------------------------------------------------------
# Data
# ------------------------------------------------------------
mpa_fig_data <- bind_rows(
  pisc_model_data %>%
    dplyr::select(site, mpa_status, mean_biomass) %>%
    mutate(group = "Piscivores"),
  browser_model_data %>%
    dplyr::select(site, mpa_status, mean_biomass) %>%
    mutate(group = "Browsers")
) %>%
  mutate(
    mpa_status = factor(mpa_status,
                        levels = c("none", "low", "medium"),
                        labels = c("None", "Low", "Medium")),
    group = factor(group, levels = c("Piscivores", "Browsers"))
  )

# ------------------------------------------------------------
# Annotations — p-values only
# ------------------------------------------------------------
pval_dat <- data.frame(
  group = factor(c("Piscivores", "Browsers"),
                 levels = c("Piscivores", "Browsers")),
  x     = c("Medium", "Medium"),
  y     = c(5200, 5200),
  label = c("p = 0.001", "p = 0.002")
)

sig_dat <- data.frame(
  group = factor(c("Piscivores", "Browsers"),
                 levels = c("Piscivores", "Browsers")),
  x     = c(1, 1),
  xend  = c(3, 3),
  y     = c(4800, 4800)
)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
fig_protection <- ggplot(
  mpa_fig_data,
  aes(x = mpa_status, y = mean_biomass, colour = mpa_status)
) +
  
  geom_boxplot(
    aes(fill = mpa_status),
    alpha         = 0.10,
    colour        = "grey60",
    outlier.shape = NA,
    width         = 0.45,
    linewidth     = 0.35
  ) +
  
  geom_beeswarm(
    size  = 2.0,
    alpha = 0.80,
    cex   = 2.5
  ) +
  
  # Significance brackets
  geom_segment(
    data = sig_dat,
    aes(x = x, xend = xend, y = y, yend = y),
    colour      = "grey55",
    linewidth   = 0.35,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = sig_dat,
    aes(x = x, xend = x, y = y, yend = y - 150),
    colour      = "grey55",
    linewidth   = 0.35,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = sig_dat,
    aes(x = xend, xend = xend, y = y, yend = y - 150),
    colour      = "grey55",
    linewidth   = 0.35,
    inherit.aes = FALSE
  ) +
  
  # p-value annotations
  geom_text(
    data = pval_dat,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust    = 0.5,
    vjust    = -0.3,
    size     = 3.4,
    colour   = "grey30",
    fontface = "plain",
    family   = font
  ) +
  
  scale_colour_manual(values = mpa_cols, guide = "none") +
  scale_fill_manual(values   = mpa_cols, guide = "none") +
  
  scale_y_continuous(
    labels = comma,
    limits = c(0, 5800),
    expand = c(0, 0)
  ) +
  
  facet_wrap(~ group, nrow = 1) +
  
  labs(
    x = "Protection level",
    y = "Mean biomass\n(g per site)"
  ) +
  
  theme_classic(base_family = font) +
  
  theme(
    axis.title.x = element_text(
      size   = 13,
      face   = "bold",
      margin = margin(t = 8)
    ),
    
    axis.title.y = element_text(
      size       = 13,
      face       = "bold",
      angle      = 0,
      lineheight = 0.95,
      hjust      = 0.5,
      vjust      = 0.5,
      margin     = margin(r = 10)
    ),
    
    axis.text.x = element_text(
      size   = 12,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size   = 11,
      colour = "black"
    ),
    
    axis.line = element_line(
      colour    = border_col,
      linewidth = 0.45
    ),
    
    axis.ticks = element_line(
      colour    = border_col,
      linewidth = 0.45
    ),
    
    strip.background = element_rect(
      fill      = "grey92",
      colour    = border_col,
      linewidth = 0.4
    ),
    
    strip.text = element_text(
      size   = 12,
      face   = "bold",
      family = font
    ),
    
    panel.border = element_rect(
      color     = border_col,
      fill      = NA,
      linewidth = 0.4
    ),
    
    panel.background = element_rect(
      fill  = "white",
      color = NA
    ),
    
    panel.spacing.x = unit(0.18, "lines"),
    
    aspect.ratio = 1.15,
    
    plot.margin = margin(
      t = 5,
      r = 55,
      b = 5,
      l = 5
    )
  )

fig_protection

# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
ggsave(
  "mpa_protection_final.pdf",
  fig_protection,
  width       = 8.0,
  height      = 5.6,
  units       = "in",
  bg          = "white",
  useDingbats = FALSE
)

ggsave(
  "mpa_protection_final.png",
  fig_protection,
  width  = 8.0,
  height = 5.6,
  units  = "in",
  dpi    = 300,
  bg     = "white"
)