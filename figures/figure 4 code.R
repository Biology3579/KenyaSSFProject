# FUNCTIONAL GROUPS FIGURE

library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

# ── DATA ─────────────────────────────────────────────────────

biomass_effects <- data.frame(
  group = rep(c("Grazers", "Scrapers", "Lg. excav.",
                "Browsers", "Piscivores", "Corallivores"), each = 6),
  predictor = rep(c("Rugosity", "SST", "Chl-a",
                    "Sett. pop.", "Sett. grav.", "Mkt grav."), times = 6),
  effect = c(
    "pos",  "neg",  "none", "none", "none",  "none",
    "none", "none", "pos",  "neg",  "neg",   "none",
    "pos",  "none", "none", "none", "none",  "none",
    "pos",  "none", "none", "none", "none",  "none",
    "pos",  "none", "none", "none", "none",  "pos_c",
    "null", "null", "null", "null", "null",  "null"
  )
)

abundance_effects <- data.frame(
  group = rep(c("Grazers", "Scrapers", "Lg. excav.",
                "Browsers", "Piscivores", "Corallivores"), each = 6),
  predictor = rep(c("Rugosity", "SST", "Chl-a",
                    "Sett. pop.", "Sett. grav.", "Mkt grav."), times = 6),
  effect = c(
    "pos",  "none", "neg",  "neg",  "none", "none",
    "none", "none", "none", "none", "neg",  "none",
    "none", "none", "none", "none", "none", "none",
    "pos",  "none", "none", "none", "none", "pos",
    "none", "none", "neg",  "none", "none", "pos",
    "null", "null", "null", "null", "null", "null"
  )
)

# ── LOOKUPS ──────────────────────────────────────────────────

group_levels <- c("Corallivores", "Piscivores", "Browsers",
                  "Lg. excav.", "Scrapers", "Grazers")

pred_levels  <- c("Rugosity", "SST", "Chl-a",
                  "Sett. pop.", "Sett. grav.", "Mkt grav.")

fill_values <- c(
  "pos"   = "#D55E00",
  "neg"   = "#0072B2",
  "pos_c" = "#F0A070",
  "none"  = "#f0f0f0",
  "null"  = "#c8c6c0"
)

label_values <- c(
  "pos"   = "+",
  "neg"   = "\u2212",
  "pos_c" = "\u2020",   # † dagger
  "none"  = "",
  "null"  = ""

text_cols <- c(
  "pos"   = "white",
  "neg"   = "white",
  "pos_c" = "#7a2e00",
  "none"  = NA,
  "null"  = NA
)

# ── SHARED THEME ─────────────────────────────────────────────

theme_heat <- theme_minimal(base_size = 11) +
  theme(
    panel.grid      = element_blank(),
    axis.ticks      = element_blank(),
    axis.title      = element_blank(),
    axis.text.y     = element_blank(),
    axis.text.x     = element_text(size = 8.5, colour = "grey30",
                                   angle = 40, hjust = 0, vjust = 0),
    plot.title      = element_text(hjust = 0.5, size = 12,
                                   face = "bold", colour = "grey10",
                                   margin = margin(b = 2)),
    plot.subtitle   = element_text(hjust = 0.5, size = 8,
                                   colour = "grey55",
                                   margin = margin(b = 10)),
    legend.position = "none",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin     = margin(t = 6, r = 20, b = 6, l = 6)
  )

# ── BUILD PANELS ─────────────────────────────────────────────

panels <- list()

for (i in 1:2) {
  
  df  <- if (i == 1) biomass_effects else abundance_effects
  ttl <- if (i == 1) "A  Biomass" else "B  Abundance"
  sub <- if (i == 1) "site-level and transect-level models" else "transect-level count models"
  
  df <- df %>%
    mutate(
      group     = factor(group, levels = group_levels),
      predictor = factor(predictor, levels = pred_levels),
      label     = label_values[effect]
    )
  
  outlined <- df %>% dplyr::filter(effect == "pos_c")
  
  panels[[i]] <- ggplot(df, aes(x = predictor, y = group)) +
    geom_tile(aes(fill = effect), colour = "white",
              linewidth = 0.8, width = 0.88, height = 0.82) +
    geom_tile(data = outlined, fill = NA, colour = "#a04000",
              linewidth = 1.2, width = 0.88, height = 0.82) +
    geom_text(aes(label = label, colour = effect),
              size = 4, fontface = "bold", na.rm = TRUE) +
    scale_fill_manual(values = fill_values, na.value = "#f0f0f0") +
    scale_colour_manual(values = text_cols, na.value = NA) +
    scale_x_discrete(position = "top", expand = c(0.04, 0.04)) +
    scale_y_discrete(expand = c(0.06, 0.06), labels = NULL) +
    labs(title = ttl, subtitle = sub) +
    theme_heat
}

# ── ROW LABEL GROB ───────────────────────────────────────────
# Calculate y positions analytically.
#
# The full grob height = panel area + legend area (heights c(8,1)).
# Panel area fraction = 8/9 of total grob height.
#
# Within the panel grob:
#   - title + subtitle + top margin consume the top portion
#   - bottom margin at the bottom
#   - the actual grid (plot area) sits in between
#
# With theme_minimal, title (12pt) + subtitle (8pt) + margins:
#   top consumed  ~ 0.38 of panel grob height
#   bottom margin ~ 0.04 of panel grob height
#   plot area     ~ 0.58 of panel grob height
#
# Within the plot area, scale_y_discrete with expand=c(0.06,0.06)
# places 6 levels at positions that account for the expand padding.
# The expand adds 0.06 * total_range padding on each side.
# With 6 levels: total data range = 1 to 6, expanded range = 0.7 to 6.3
# Each level i sits at position i within [0.7, 6.3]
# As fraction of expanded range: (i - 0.7) / (6.3 - 0.7)
#
# Then map to grob npc (note: grob y=1 is TOP, data y=6 is TOP):
# grob_y = panel_top_in_grob + plot_area_height * (1 - data_fraction)
# But we also need to account for legend row:
# label_grob sits in top row which is 8/9 of total height
# so npc coordinates within label_grob are already 0-1 for that row

n_groups    <- 6
expand_pad  <- 0.06
data_range  <- c(1 - expand_pad * n_groups,
                 n_groups + expand_pad * n_groups)  # ~0.64 to 6.36

# fraction of plot area from bottom for each level (1=bottom, 6=top)
level_frac  <- (seq(1, n_groups) - data_range[1]) /
  (data_range[2] - data_range[1])

# within label_grob npc (0=bottom, 1=top of label_grob row)
# plot area sits between plot_bottom_npc and plot_top_npc
plot_top_npc    <- 0.78   # was 0.80, shift down
plot_bottom_npc <- 0.015   # was 0.22, wider spacing
plot_area_h     <- plot_top_npc - plot_bottom_npc

row_ys_npc <- plot_bottom_npc + plot_area_h * level_frac

cat("row_ys_npc:", round(row_ys_npc, 3), "\n")

legend_cols   <- c("#D55E00", "#0072B2", "#f0f0f0", "#c8c6c0", "#F0A070")
legend_labels <- c("Positive (p < 0.05)", "Negative (p < 0.05)",
                   "Not supported", "Null result",
                   "Supported\u2020")
legend_xs     <- seq(0.10, 0.90, length.out = 5)

legend_grob <- gTree(children = do.call(gList, c(
  list(rectGrob(gp = gpar(fill = "white", col = NA))),
  lapply(seq_along(legend_cols), function(i) {
    gList(
      rectGrob(
        x      = legend_xs[i],
        y      = 0.68,
        width  = unit(0.09, "npc"),
        height = unit(0.28, "npc"),
        gp     = gpar(
          fill = legend_cols[i],
          col  = if (i == 5) "#a04000" else NA,
          lwd  = if (i == 5) 2 else 0
        )
      ),
      textGrob(
        label = legend_labels[i],
        x     = legend_xs[i],
        y     = 0.25,
        just  = "centre",
        gp    = gpar(fontsize = 7.5, col = "grey30")
      )
    )
  })
)))

fig4 <- arrangeGrob(
  arrangeGrob(
    label_grob,
    panels[[1]],
    panels[[2]],
    ncol   = 3,
    widths = c(1.2, 1, 1)
  ),
  legend_grob,
  nrow    = 2,
  heights = c(8, 1)
)

ggsave("figure4_functional_heatmap.pdf", plot = fig4,
       width = 190, height = 135, units = "mm", device = cairo_pdf)

ggsave("figure4_functional_heatmap.png", plot = fig4,
       width = 190, height = 135, units = "mm", dpi = 300)

cat("Figure 4 saved.\n")

