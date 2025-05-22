setwd("D:/college/year_3/CMML/ICA2/scenario1")

library(ggplot2)
library(patchwork)
library(readr)
library(colorspace)
library(RColorBrewer)
library(Polychrome)
library(scales)

# === Step 1: Load UMAP data
umap_files <- c(
  raw       = "umap_raw.rds",
  seurat    = "umap_seurat_integrated.rds",
  scvi      = "umap_scvi.rds",
  fastmnn   = "umap_fastmnn.rds",
  harmony   = "umap_harmony.rds"
)
dfs <- lapply(umap_files, read_rds)

# === Step 2: Rename cell_type labels before setting factor levels
for (i in seq_along(dfs)) {
  dfs[[i]]$cell_type <- recode(
    dfs[[i]]$cell_type,
    "Type 2" = "Alveolar Type 2"
  )
}

# Step 3: Build palettes based on ORDER, not label content
all_celltypes <- sort(Reduce(intersect, lapply(dfs, \(d) unique(na.omit(d$cell_type)))))
all_batches   <- sort(Reduce(intersect, lapply(dfs, \(d) unique(na.omit(d$batch)))))

# Enforce consistent factor levels for all UMAPs
for (i in seq_along(dfs)) {
  dfs[[i]]$cell_type <- factor(dfs[[i]]$cell_type, levels = all_celltypes)
  dfs[[i]]$batch     <- factor(dfs[[i]]$batch, levels = all_batches)
}

# === Step 4: Define palettes
dimplot_hue_palette <- function(values) {
  pal <- scales::hue_pal(h = c(0, 360), c = 100, l = 65)(length(values))
  setNames(pal, values)
}
celltype_cols <- dimplot_hue_palette(all_celltypes)
batch_cols    <- dimplot_hue_palette(all_batches)

# === Step 5: Plotting functions with auto-scaling square axes
expand_equal_range <- function(x, y, pad = 0.02) {
  x_span <- diff(range(x, na.rm = TRUE))
  y_span <- diff(range(y, na.rm = TRUE))
  max_span <- max(x_span, y_span) * (1 + pad)
  
  x_center <- mean(range(x, na.rm = TRUE))
  y_center <- mean(range(y, na.rm = TRUE))
  
  list(
    x = c(x_center - max_span / 2, x_center + max_span / 2),
    y = c(y_center - max_span / 2, y_center + max_span / 2)
  )
}

plot_umap <- function(df, colour_col, title, palette) {
  lims <- expand_equal_range(df$UMAP_1, df$UMAP_2)
  ggplot(df, aes(UMAP_1, UMAP_2, colour = .data[[colour_col]])) +
    geom_point(size = 0.6, alpha = 0.85) +
    coord_fixed(xlim = lims$x, ylim = lims$y) +
    theme_void(base_size = 14) +
    scale_colour_manual(
      values = palette,
      guide  = guide_legend(override.aes = list(size = 4), ncol = 2)
    ) +
    labs(title = title, colour = colour_col) +
    theme(
      plot.title     = element_text(hjust = 0.5, size = 18, face = "bold"),
      legend.title   = element_text(size = 14),
      legend.text    = element_text(size = 12),
      legend.key.size = unit(0.7, "cm"),
      panel.border   = element_rect(colour = "grey70", linetype = "dashed", linewidth = 0.5, fill = NA)
    )
}

plot_umap_batch <- function(df, title, palette) {
  lims <- expand_equal_range(df$UMAP_1, df$UMAP_2)
  ggplot(df, aes(UMAP_1, UMAP_2, colour = batch)) +
    geom_point(size = 0.6, alpha = 0.85) +
    coord_fixed(xlim = lims$x, ylim = lims$y) +
    theme_void(base_size = 14) +
    scale_colour_manual(
      values = palette,
      guide  = guide_legend(override.aes = list(size = 4), ncol = 1)
    ) +
    labs(title = title, colour = "batch") +
    theme(
      plot.title     = element_text(hjust = 0.5, size = 18, face = "bold"),
      legend.title   = element_text(size = 14),
      legend.text    = element_text(size = 12),
      legend.key.size = unit(0.7, "cm"),
      panel.border   = element_rect(colour = "grey70", linetype = "dashed", linewidth = 0.5, fill = NA)
    )
}

# === Step 6: Create cell type plots
plots_type <- list(
  plot_umap(dfs$raw,     "cell_type", "Raw",     celltype_cols),
  plot_umap(dfs$seurat,  "cell_type", "Seurat",  celltype_cols),
  plot_umap(dfs$scvi,    "cell_type", "scVI",    celltype_cols),
  plot_umap(dfs$fastmnn, "cell_type", "fastMNN", celltype_cols),
  plot_umap(dfs$harmony, "cell_type", "Harmony", celltype_cols)
)

fig_type <- wrap_plots(plots_type, ncol = 3) +
  plot_layout(guides = "collect", widths = rep(1, 3), heights = rep(1, 2)) &
  theme(
    legend.position = "right",
    plot.margin = margin(2, 2, 2, 2)
  )

ggsave("UMAP_all_celltype.png", fig_type, width = 16, height = 9, dpi = 300)
print(fig_type)

# === Step 7: Create batch plots
plots_batch <- list(
  plot_umap_batch(dfs$raw,     "Raw",     batch_cols),
  plot_umap_batch(dfs$seurat,  "Seurat",  batch_cols),
  plot_umap_batch(dfs$scvi,    "scVI",    batch_cols),
  plot_umap_batch(dfs$fastmnn, "fastMNN", batch_cols),
  plot_umap_batch(dfs$harmony, "Harmony", batch_cols)
)

fig_batch <- wrap_plots(plots_batch, ncol = 3) +
  plot_layout(guides = "collect", widths = rep(1, 3), heights = rep(1, 2)) &
  theme(
    legend.position = "right",
    plot.margin = margin(2, 2, 2, 2)
  )

ggsave("UMAP_all_batch.png", fig_batch, width = 13, height = 9, dpi = 300)
print(fig_batch)
