library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)
library(cowplot)
library(grid)
library(scales)

font <- "serif"

col_conn <- "#4A7FA5"
col_near <- "#6A5ACD"
col_far  <- "#4A7C59"
col_obs  <- "grey45"

border_col <- "grey75"

base_panel_theme <- theme_classic(base_family = font) +
  theme(
    axis.title.x     = element_text(size = 13, face = "bold"),
    axis.title.y     = element_text(size = 13, face = "bold"),
    axis.text.x      = element_text(size = 12, colour = "black"),
    axis.text.y      = element_text(size = 12, colour = "black"),
    axis.line        = element_line(colour = border_col, linewidth = 0.45),
    axis.ticks       = element_line(linewidth = 0.45, colour = border_col),
    panel.border     = element_rect(color = border_col, fill = NA, linewidth = 0.4),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.position  = "none",
    plot.margin      = margin(8, 8, 8, 8),
    strip.background = element_rect(fill = "grey92", color = border_col, linewidth = 0.4),
    strip.text       = element_text(
      size = 11, face = "bold", family = font,
      margin = margin(3, 0, 3, 0)
    )
  )

make_shared_label <- function(txt) {
  textGrob(
    txt,
    x = 0.5, y = 0.5,
    hjust = 0.5, vjust = 0.5,
    gp = gpar(
      fontfamily = font,
      fontface   = "bold",
      fontsize   = 13
    )
  )
}

make_shared_x <- function(txt) {
  textGrob(
    txt,
    x = 0.55, y = 0.60,   # fixed: was 0.55, now centred
    hjust = 0.5, vjust = 0.5,
    gp = gpar(
      fontfamily = font,
      fontface   = "bold",
      fontsize   = 13
    )
  )
}

add_title_box <- function(
    p,
    txt,
    strip_h = 0.085,
    gap_h   = 0.005,
    inset   = 0
) {
  title_strip <- ggdraw() +
    draw_grob(
      rectGrob(
        gp = gpar(fill = "grey92", col = border_col, lwd = 0.4)
      ),
      x = inset, y = 0,
      width  = 1 - 2 * inset,
      height = 1,
      hjust = 0, vjust = 0
    ) +
    draw_label(
      txt,
      x = 0.5, y = 0.5,
      hjust = 0.5, vjust = 0.5,
      fontface = "bold", fontfamily = font, size = 13
    ) +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  plot_grid(
    title_strip,
    NULL,
    p + theme(plot.margin = margin(0, 8, 8, 8)),
    ncol = 1,
    rel_heights = c(strip_h, gap_h, 1),
    align = "v", axis = "lr"
  )
}

# ── Panel A: Browsers ─────────────────────────────────────────
plot_browser_raw <- ggplot() +
  geom_ribbon(
    data = b_conn_grid,
    aes(x = connectivity_sc, ymin = lwr, ymax = upr),
    fill = col_conn, alpha = 0.10
  ) +
  geom_line(
    data = b_conn_grid,
    aes(x = connectivity_sc, y = fit),
    colour = col_conn, linewidth = 1.1
  ) +
  geom_point(
    data = b_obs,
    aes(x = connectivity_sc, y = mean_biomass),
    colour = col_obs, size = 1.5, alpha = 0.5
  ) +
  scale_y_log10(
    breaks = c(100, 200, 500, 1000, 2000, 5000, 10000),
    labels = scales::comma,
    expand = expansion(mult = c(0.03, 0.02))
  ) +
  labs(x = NULL, y = NULL) +
  base_panel_theme

plot_browser <- add_title_box(plot_browser_raw, "Browsers")

# ── Panel A: Corallivores ─────────────────────────────────────
plot_corallivore_raw <- ggplot() +
  geom_ribbon(
    data = c_conn_grid,
    aes(x = connectivity_sc, ymin = lwr, ymax = upr),
    fill = col_conn, alpha = 0.10
  ) +
  geom_line(
    data = c_conn_grid,
    aes(x = connectivity_sc, y = fit),
    colour = col_conn, linewidth = 1.1
  ) +
  geom_point(
    data = c_obs,
    aes(x = connectivity_sc, y = mean_biomass),
    colour = col_obs, size = 1.5, alpha = 0.5
  ) +
  scale_y_log10(
    breaks = c(10, 20, 50, 100, 200),
    labels = scales::comma,
    expand = expansion(mult = c(0.03, 0.10))
  ) +
  labs(x = NULL, y = NULL) +
  base_panel_theme

plot_corallivore <- add_title_box(plot_corallivore_raw, "Corallivores")

# ── Panel B: Piscivores ───────────────────────────────────────
p_grid <- rbind(p_grid_high, p_grid_low)

p_grid$market <- ifelse(
  p_grid$market == "Near markets (+1 SD)",
  "Near markets",
  "Far from markets"
)

p_grid$market <- factor(
  p_grid$market,
  levels = c("Near markets", "Far from markets")
)

plot_pisc_raw <- ggplot() +
  geom_ribbon(
    data = p_grid,
    aes(x = connectivity_sc, ymin = lwr, ymax = upr, fill = market),
    alpha = 0.10
  ) +
  geom_line(
    data = p_grid,
    aes(x = connectivity_sc, y = fit, colour = market),
    linewidth = 1.1
  ) +
  scale_colour_manual(
    values = c("Near markets" = col_near, "Far from markets" = col_far),
    name = NULL
  ) +
  scale_fill_manual(
    values = c("Near markets" = col_near, "Far from markets" = col_far),
    name = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  scale_y_log10(
    breaks = c(100, 200, 500, 1000, 2000),
    labels = scales::comma,
    expand = expansion(mult = c(0.03, 0.06))
  ) +
  labs(x = NULL, y = NULL) +
  # Use base_panel_theme then override legend settings only
  base_panel_theme +
  theme(
    legend.position      = c(0.95, 0.35),
    legend.justification = c(1, 1),
    legend.background    = element_rect(
      fill = "white", colour = border_col, linewidth = 0.25
    ),
    legend.text     = element_text(size = 12, family = font),
    legend.key.size = unit(0.7, "lines"),
    legend.margin   = margin(2, 3, 2, 3)
  )

plot_pisc <- add_title_box(plot_pisc_raw, "Piscivores")

# ── Shared labels ─────────────────────────────────────────────
# Single y-label grob reused for both rows
y_lab <- make_shared_label("Biomass\nper site (g)")
x_lab <- make_shared_x("Larval connectivity (standardized)")

# ── Assemble top row ──────────────────────────────────────────
top_row_plots <- plot_grid(
  plot_browser,
  plot_corallivore,
  ncol = 2, rel_widths = c(1, 1), align = "h"
)

top_block <- plot_grid(
  y_lab, top_row_plots,
  ncol = 2, rel_widths = c(0.12, 0.88), align = "h"
)

top_block <- plot_grid(
  top_block, x_lab,
  ncol = 1, rel_heights = c(1, 0.065)
)

# ── Bottom row ────────────────────────────────────────────────
bottom_block <- plot_grid(
  y_lab, plot_pisc,
  ncol = 2, rel_widths = c(0.12, 0.88), align = "h"
)

bottom_block <- plot_grid(
  bottom_block, x_lab,
  ncol = 1, rel_heights = c(1, 0.065)
)

# ── Final figure ──────────────────────────────────────────────
fig4_core <- plot_grid(
  top_block,
  NULL,
  bottom_block,
  ncol = 1,
  rel_heights = c(1, 0.10, 1),   # fixed: equal row heights (was 1, 0.10, 1.10)
  align = "v"
)

fig4 <- ggdraw(fig4_core) +
  theme(plot.margin = margin(12, 12, 12, 12)) +
  draw_label(
    "A.",
    x = 0.025, y = 0.998,
    hjust = 0, vjust = 1,
    fontface = "bold", size = 15, fontfamily = font
  ) +
  draw_label(
    "B.",
    x = 0.025, y = 0.52,   # nudged slightly up from 0.50 to clear the bottom block
    hjust = 0, vjust = 1,
    fontface = "bold", size = 15, fontfamily = font
  )

fig4

ggsave(
  "fig4_v5.pdf", fig4,
  width = 9, height = 7,
  units = "in", bg = "white",
  useDingbats = FALSE
)

ggsave(
  "fig4_v5.png", fig4,
  width = 9, height = 7,
  units = "in", dpi = 300,
  bg = "white"
)

cat("Figure 4 v5 saved\n")