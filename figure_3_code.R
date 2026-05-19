# Packages
library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)
library(cowplot)

font <- "serif"

pred_cols <- c(
  "Rugosity"           = "#A0522D",
  "Chlorophyll-a"      = "#4A7C59",
  "Settlement gravity" = "#6A5ACD",
  "Human pressure"     = "#6A5ACD",
  "Connectivity"       = "#4A7FA5"
)

x_scale <- scale_x_continuous(
  limits = c(-0.925, 0.925),
  breaks = c(-0.8, -0.4, 0, 0.4, 0.8),
  expand = expansion(mult = c(0, 0))
)

base_panel_theme <- theme_classic(base_family = font) +
  theme(
    axis.title.x     = element_text(size = 13),
    axis.text.y      = element_text(size = 12),
    axis.text.x      = element_text(size = 12),
    axis.line        = element_blank(),
    axis.ticks       = element_line(linewidth = 0.45, color = "grey30"),
    panel.border     = element_rect(color = "grey75", fill = NA, linewidth = 0.4),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.position  = "none",
    plot.margin      = margin(4, 4, 4, 4),
    strip.background = element_rect(fill = "grey92", color = "grey75", linewidth = 0.4),
    strip.text       = element_text(size = 11, face = "bold", family = font,
                                    margin = margin(3, 0, 3, 0))
  )

# 1. Total Biomass ----
total_biomass <- tribble(
  ~term,               ~beta,  ~lower,  ~upper,  ~status,
  "Rugosity",           0.216,  0.056,   0.376,  "supported",
  "Chlorophyll-a",     -0.125, -0.322,   0.072,  "supported",
  "Settlement gravity",-0.251, -0.445,  -0.057,  "supported",
  "Connectivity",       0.000,     NA,      NA,   "unsupported"
)

plot_total <- ggplot(total_biomass, aes(x = beta, y = term)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed", color = "grey80") +
  geom_errorbarh(
    data = filter(total_biomass, status == "supported"),
    aes(xmin = lower, xmax = upper, color = term),
    height = 0.12, linewidth = 0.95
  ) +
  geom_point(
    data = filter(total_biomass, status == "supported"),
    aes(color = term), shape = 16, size = 2.5
  ) +
  geom_point(
    data = filter(total_biomass, status == "unsupported"),
    aes(x = 0), shape = 16, color = "grey55", size = 2.5
  ) +
  scale_color_manual(values = pred_cols) +
  scale_y_discrete(limits = c("Connectivity", "Settlement gravity", "Chlorophyll-a", "Rugosity")) +
  x_scale +
  labs(x = expression("Standardised effect size (" * beta * ")"), y = NULL) +
  base_panel_theme   # no facet_wrap

# 2. Scrapers ----
scrapers <- tribble(
  ~term,             ~beta,  ~lower,  ~upper,  ~status,
  "Rugosity",         0.143, -0.024,   0.351,  "supported",
  "Chlorophyll-a",    0.204,  0.013,   0.396,  "supported",
  "Human pressure",  -0.261, -0.475,  -0.047,  "supported",
  "Connectivity",     0.000,     NA,      NA,   "unsupported"
) |> mutate(panel = "Scrapers/small excavators")

plot_scrapers <- ggplot(scrapers, aes(x = beta, y = term)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed", color = "grey80") +
  geom_errorbarh(
    data = filter(scrapers, status == "supported"),
    aes(xmin = lower, xmax = upper, color = term),
    height = 0.12, linewidth = 0.95
  ) +
  geom_point(
    data = filter(scrapers, status == "supported"),
    aes(color = term), shape = 16, size = 2.5
  ) +
  geom_point(
    data = filter(scrapers, status == "unsupported"),
    aes(x = 0), shape = 16, color = "grey55", size = 2.5
  ) +
  scale_color_manual(values = pred_cols) +
  scale_y_discrete(limits = c("Connectivity", "Human pressure", "Chlorophyll-a", "Rugosity")) +
  x_scale +
  facet_wrap(~panel) +
  labs(x = expression("Standardised effect size (" * beta * ")"), y = NULL) +
  base_panel_theme

# 3. Piscivores ----
piscivores <- tribble(
  ~term,             ~beta,  ~lower,  ~upper,  ~status,
  "Rugosity",         0.392,  0.103,   0.682,  "supported",
  "Chlorophyll-a",    0.018, -0.274,   0.311,  "supported",
  "Human pressure",   0.218, -0.135,   0.571,  "supported",
  "Connectivity",     0.353,  0.084,   0.622,  "supported"
) |> mutate(panel = "Piscivores")

plot_piscivores <- ggplot(piscivores, aes(x = beta, y = term)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed", color = "grey80") +
  geom_errorbarh(
    data = filter(piscivores, status == "supported"),
    aes(xmin = lower, xmax = upper, color = term),
    height = 0.12, linewidth = 0.95
  ) +
  geom_point(
    data = filter(piscivores, status == "supported"),
    aes(color = term), shape = 16, size = 2.5
  ) +
  scale_color_manual(values = pred_cols) +
  scale_y_discrete(limits = c("Connectivity", "Human pressure", "Chlorophyll-a", "Rugosity")) +
  x_scale +
  facet_wrap(~panel) +
  labs(x = expression("Standardised effect size (" * beta * ")"), y = NULL) +
  base_panel_theme

# 4. Browsers ----
browsers <- tribble(
  ~term,             ~beta,  ~lower,  ~upper,  ~status,
  "Rugosity",         0.479,  0.101,   0.854,  "supported",
  "Chlorophyll-a",   -0.106, -0.481,   0.269,  "supported",
  "Human pressure",   0.000,     NA,      NA,   "unsupported",
  "Connectivity",     0.421,  0.108,   0.734,  "supported"
) |> mutate(panel = "Browsers")

plot_browsers <- ggplot(browsers, aes(x = beta, y = term)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed", color = "grey80") +
  geom_errorbarh(
    data = filter(browsers, status == "supported"),
    aes(xmin = lower, xmax = upper, color = term),
    height = 0.12, linewidth = 0.95
  ) +
  geom_point(
    data = filter(browsers, status == "supported"),
    aes(color = term), shape = 16, size = 2.5
  ) +
  geom_point(
    data = filter(browsers, status == "unsupported"),
    aes(x = 0), shape = 16, color = "grey55", size = 2.5
  ) +
  scale_color_manual(values = pred_cols) +
  scale_y_discrete(limits = c("Connectivity", "Human pressure", "Chlorophyll-a", "Rugosity")) +
  x_scale +
  facet_wrap(~panel) +
  labs(x = expression("Standardised effect size (" * beta * ")"), y = NULL) +
  base_panel_theme

# 5. Corallivores ----
corallivores <- tribble(
  ~term,             ~beta,  ~lower,  ~upper,  ~status,
  "Rugosity",         0.108, -0.163,   0.379,  "supported",
  "Chlorophyll-a",    0.121, -0.063,   0.305,  "supported",
  "Human pressure",   0.000,     NA,      NA,   "unsupported",
  "Connectivity",    -0.221, -0.399,  -0.043,  "supported"
) |> mutate(panel = "Corallivores")

plot_corallivores <- ggplot(corallivores, aes(x = beta, y = term)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed", color = "grey80") +
  geom_errorbarh(
    data = filter(corallivores, status == "supported"),
    aes(xmin = lower, xmax = upper, color = term),
    height = 0.12, linewidth = 0.95
  ) +
  geom_point(
    data = filter(corallivores, status == "supported"),
    aes(color = term), shape = 16, size = 2.5
  ) +
  geom_point(
    data = filter(corallivores, status == "unsupported"),
    aes(x = 0), shape = 16, color = "grey55", size = 2.5
  ) +
  scale_color_manual(values = pred_cols) +
  scale_y_discrete(limits = c("Connectivity", "Human pressure", "Chlorophyll-a", "Rugosity")) +
  x_scale +
  facet_wrap(~panel) +
  labs(x = expression("Standardised effect size (" * beta * ")"), y = NULL) +
  base_panel_theme

# 6. Grazer-detritivores ----
grazer_detritivores <- tribble(
  ~term,             ~beta,  ~lower,  ~upper,  ~status,
  "Rugosity",         0.241,  0.047,   0.435,  "supported",
  "Chlorophyll-a",   -0.125, -0.335,   0.084,  "supported",
  "Human pressure",   0.000,     NA,      NA,   "unsupported",
  "Connectivity",     0.000,     NA,      NA,   "unsupported"
) |> mutate(panel = "Grazers/detritivores")

plot_grazers <- ggplot(grazer_detritivores, aes(x = beta, y = term)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed", color = "grey80") +
  geom_errorbarh(
    data = filter(grazer_detritivores, status == "supported"),
    aes(xmin = lower, xmax = upper, color = term),
    height = 0.12, linewidth = 0.95
  ) +
  geom_point(
    data = filter(grazer_detritivores, status == "supported"),
    aes(color = term), shape = 16, size = 2.5
  ) +
  geom_point(
    data = filter(grazer_detritivores, status == "unsupported"),
    aes(x = 0), shape = 16, color = "grey55", size = 2.5
  ) +
  scale_color_manual(values = pred_cols) +
  scale_y_discrete(limits = c("Connectivity", "Human pressure", "Chlorophyll-a", "Rugosity")) +
  x_scale +
  facet_wrap(~panel) +
  labs(x = expression("Standardised effect size (" * beta * ")"), y = NULL) +
  base_panel_theme

# 7. Large excavators ----
large_excavators <- tribble(
  ~term,             ~beta,  ~lower,  ~upper,  ~status,
  "Rugosity",         0.506,  0.128,   0.884,  "supported",
  "Chlorophyll-a",   -0.300, -0.715,   0.116,  "supported",
  "Human pressure",   0.000,     NA,      NA,   "unsupported",
  "Connectivity",     0.000,     NA,      NA,   "unsupported"
) |> mutate(panel = "Large excavators")

plot_excavators <- ggplot(large_excavators, aes(x = beta, y = term)) +
  geom_vline(xintercept = 0, linewidth = 0.5, linetype = "dashed", color = "grey80") +
  geom_errorbarh(
    data = filter(large_excavators, status == "supported"),
    aes(xmin = lower, xmax = upper, color = term),
    height = 0.12, linewidth = 0.95
  ) +
  geom_point(
    data = filter(large_excavators, status == "supported"),
    aes(color = term), shape = 16, size = 2.5
  ) +
  geom_point(
    data = filter(large_excavators, status == "unsupported"),
    aes(x = 0), shape = 16, color = "grey55", size = 2.5
  ) +
  scale_color_manual(values = pred_cols) +
  scale_y_discrete(limits = c("Connectivity", "Human pressure", "Chlorophyll-a", "Rugosity")) +
  x_scale +
  facet_wrap(~panel) +
  labs(x = expression("Standardised effect size ( " * beta * " )"), y = NULL) +
  base_panel_theme

# Legend
legend_data <- tibble(
  label  = c("Rugosity", "Chlorophyll-a", "Human pressure", "Connectivity", "Unsupported"),
  colour = c("#A0522D", "#4A7C59", "#6A5ACD", "#4A7FA5", "grey55"),
  y      = 5:1
)

legend_plot <- ggplot(legend_data, aes(x = 0.2, y = y, color = colour)) +
  geom_point(shape = 16, size = 3.5) +  
  geom_text(
    aes(x = 0.38, label = label),
    color = "grey15", hjust = 0, size = 3.6, family = font
  ) +
  annotate(
    "text", x = 0.2, y = 5.9,
    label = "Legend", hjust = 0,
    fontface = "bold", size = 4.5, family = font, color = "grey10"
  ) +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.1, 2.0)) +
  scale_y_continuous(limits = c(0.5, 6.0)) +
  theme_void(base_family = font) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(6, 10, 6, 6)
  )

# Assemble ----
hide_y <- theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
hide_x <- theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                axis.ticks.x = element_blank())

row1 <- (plot_scrapers   + hide_x + theme(plot.margin = margin(4,4,8,4))) |
  (plot_piscivores + hide_y + hide_x + theme(plot.margin = margin(4,4,8,4))) |
  (plot_browsers   + hide_y + hide_x + theme(plot.margin = margin(4,4,8,4)))

row2 <- (plot_corallivores + theme(axis.title.x = element_blank())) |
  (plot_grazers      + hide_y) |
  (plot_excavators   + hide_y + theme(axis.title.x = element_blank()))

panel_b <- (row1 / row2) +
  plot_layout(heights = c(1, 1))

# Assemble without tags
top_row <- wrap_elements(full =
                           (plot_total | legend_plot) + plot_layout(widths = c(3.4, 1))
)

fig3_notags <- (top_row / wrap_elements(full = panel_b)) +
  plot_layout(heights = c(1, 1.7)) &
  theme(plot.margin = margin(t = 15, r = 5, b = 5, l = 5))

fig3 <- ggdraw(fig3_notags) +
  draw_label("A.  Total reef fish biomass",
             x = 0.02, y = 0.985,
             hjust = 0, vjust = 1,
             fontface = "bold", size = 13, fontfamily = font) +
  draw_label("B.  Functional group responses",
             x = 0.02, y = 0.625,
             hjust = 0, vjust = 1,
             fontface = "bold", size = 13, fontfamily = font)

ggsave("fig3.pdf", fig3,
       width = 10, height = 7.5,
       units = "in",
       bg = "white",
       useDingbats = FALSE)

ggsave("fig3.png", fig3,
       width = 10, height = 7.5,
       units = "in",
       dpi = 96,
       bg = "white")


