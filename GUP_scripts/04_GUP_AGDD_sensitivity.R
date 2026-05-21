library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(data.table)
library(patchwork)

# ── Configuration ─────────────────────────────────────────────────────────────
ecoregs  <- c('2.1.3','2.1.5','2.1.6','2.1.7','2.2.1','2.2.2','2.2.3','2.2.4','2.3.1','2.4.1','2.4.2','2.4.3','2.4.4')
grid_shp        <- '/projectnb/modislc/users/seamorez/HLS_Pheno/shapefiles/Arctic_grid_5km.shp'
ecoreg_shp      <- '/projectnb/modislc/users/seamorez/HLS_Pheno/shapefiles/NAA_ecoregions_final.shp'
fig_dir         <- '/projectnb/modislc/users/seamorez/HLS_Pheno/GUP_climate_sensitivity/figures'
output_dir_base <- '/projectnb/modislc/users/seamorez/HLS_Pheno/GUP_climate_sensitivity/jags_output'
input_dir_base  <- '/projectnb/modislc/users/seamorez/HLS_Pheno/GUP_climate_sensitivity/jags_input'
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

hmax     <- 100
h_thresh <- 0.01 * hmax  # same as 03_analyze_all_ecoregions.R

# ── Load spatial layers ────────────────────────────────────────────────────────
grid_sf    <- st_read(grid_shp)
ecoreg_sf  <- st_read(ecoreg_shp)
world      <- ne_countries(scale = 'medium', returnclass = 'sf')
crs_use    <- st_crs("ESRI:102001")

world_proj  <- st_transform(world, crs_use)
ecoreg_proj <- st_transform(
  ecoreg_sf %>% filter(AT_L3CODE %in% ecoregs),
  crs_use
)

# ── Loop over ecoregions — collect sensitivity results ─────────────────────────
all_sens_list     <- list()
scatter_data_list <- list()   # stores pxyr_dt for every 10th successful cell
scatter_n         <- 0L       # counts successfully processed cells

for (ecoreg in ecoregs) {

  cat(sprintf('\n══ Processing ecoregion %s ══\n', ecoreg))

  output_dir <- file.path(output_dir_base, ecoreg)
  input_dir  <- file.path(input_dir_base,  ecoreg)

  rda_files <- list.files(
    output_dir,
    pattern    = paste0('^ma_', ecoreg, '_.*\\.rda$'),
    full.names = TRUE
  )
  cat(sprintf('  Found %d rda files\n', length(rda_files)))
  if (length(rda_files) == 0) {
    cat(sprintf('  WARNING: no rda files for %s — skipping\n', ecoreg))
    next
  }

  sens_list <- vector('list', length(rda_files))

  for (k in seq_along(rda_files)) {
    f    <- rda_files[k]
    base <- basename(f)

    cell_id_str <- sub(paste0('^ma_', ecoreg, '_'), '', sub('\\.rda$', '', base))
    cell_id     <- as.integer(cell_id_str)
    if (is.na(cell_id)) next

    # ── Load model output rda ──────────────────────────────────────────────────
    env <- new.env(parent = emptyenv())
    tryCatch(load(f, envir = env),
             error = function(e) cat(sprintf('  ERROR loading %s\n', base)))

    if (!exists('gibbs', envir = env)) next
    gibbs <- env$gibbs

    if (!all(c('dat', 'data') %in% ls(env)) || is.null(gibbs$h)) next
    dat_out  <- env$dat
    data_out <- env$data

    # ── Derive preseason windows per (px_location, year) ──────────────────────
    # doy_force: first DOY the model puts hitting probability above threshold
    # doy_onset: first observed DOY with Y = 1 (actual GUP date)
    h_med <- apply(gibbs$h, 2, median)

    compare <- data.table(
      px_locs = dat_out$px_location,
      year    = dat_out$year,
      doy     = dat_out$doy,
      y       = data_out$Y,
      h       = h_med
    )

    preseason_wins <- compare[, {
      force_idx <- which(h > h_thresh)
      doy_force <- if (length(force_idx) > 0) doy[force_idx[1]] else NA_integer_
      onset_idx <- which(y == 1)
      doy_onset <- if (length(onset_idx) > 0) doy[onset_idx[1]] else NA_integer_
      list(doy_force = doy_force, doy_onset = doy_onset)
    }, by = .(px_locs, year)]

    preseason_wins <- preseason_wins[!is.na(doy_force) & !is.na(doy_onset) & doy_onset > doy_force]
    if (nrow(preseason_wins) < 5) next

    # ── Load matching input rda ────────────────────────────────────────────────
    input_file <- file.path(input_dir, paste0(ecoreg, '_', cell_id, '.rda'))
    if (!file.exists(input_file)) {
      cat(sprintf('  WARNING: no input rda for cell %d — skipping\n', cell_id))
      next
    }
    env2 <- new.env(parent = emptyenv())
    tryCatch(load(input_file, envir = env2),
             error = function(e) cat(sprintf('  ERROR loading input for cell %d\n', cell_id)))
    if (!exists('dat', envir = env2)) next
    dat_in <- as.data.table(env2$dat)

    # ── Sum GDD (°C above 0) over the preseason for each pixel-year ───────────
    # Tavg is stored in Kelvin; 0°C = 273.15 K
    # Merge preseason windows into daily data then filter to preseason DOY range
    dat_merged <- merge(
      dat_in,
      preseason_wins,
      by.x = c('px_location', 'year'),
      by.y = c('px_locs',     'year')
    )

    agdd_dt <- dat_merged[
      doy >= doy_force & doy <= doy_onset,
      .(AGDD = sum(pmax(0, Tavg - 273.15), na.rm = TRUE)),
      by = .(px_location, year)
    ]

    # Attach GUP DOY (doy_onset) to each pixel-year AGDD estimate
    pxyr_dt <- merge(
      agdd_dt,
      preseason_wins[, .(px_locs, year, doy_onset)],
      by.x = c('px_location', 'year'),
      by.y = c('px_locs',     'year')
    )
    if (nrow(pxyr_dt) < 5) next

    # ── Compute anomalies at the grid-cell level ──────────────────────────────
    # Treat all pixel-years in the cell as one population: subtract the
    # cell-wide mean from every observation. This avoids the instability of
    # per-pixel means when pixels have very few sampled years (~2 on average).
    pxyr_dt[, GUP_anom  := doy_onset - mean(doy_onset, na.rm = TRUE)]
    pxyr_dt[, AGDD_anom := AGDD      - mean(AGDD,      na.rm = TRUE)]

    # ── Linear regression: GUP_anom ~ AGDD_anom ───────────────────────────────
    # Slope units: days per degree-day (days / [°C·day])
    # Negative slope → more preseason heat → earlier GUP
    fit  <- lm(GUP_anom ~ AGDD_anom, data = pxyr_dt)
    sfit <- summary(fit)

    sens_list[[k]] <- data.frame(
      id     = cell_id,
      ecoreg = ecoreg,
      slope  = coef(fit)[['AGDD_anom']],
      r_sq   = sfit$r.squared,
      p_val  = sfit$coefficients['AGDD_anom', 'Pr(>|t|)'],
      n_obs  = nrow(pxyr_dt)
    )

    # ── Collect scatter data for every 10th successful cell ────────────────────
    scatter_n <- scatter_n + 1L
    if (scatter_n %% 10 == 1) {
      scatter_data_list[[length(scatter_data_list) + 1]] <- data.frame(
        pxyr_dt[, .(px_location, year, AGDD, doy_onset, AGDD_anom, GUP_anom)],
        cell_id     = cell_id,
        ecoreg      = ecoreg,
        panel_label = sprintf('Cell %d | %s\nR²=%.3f  β=%.4f  n=%d',
                              cell_id, ecoreg,
                              sfit$r.squared,
                              coef(fit)[['AGDD_anom']],
                              nrow(pxyr_dt))
      )
    }
  }

  all_sens_list <- c(all_sens_list, sens_list)

  cat(sprintf('  Ecoregion %s done: %d sensitivity results\n',
              ecoreg, sum(!sapply(sens_list, is.null))))
}

# ── Combine across ecoregions ──────────────────────────────────────────────────
sens_df <- bind_rows(all_sens_list)
cat(sprintf('\nTotal grid cells with sensitivity estimates: %d\n', nrow(sens_df)))

# ── Spatial join ───────────────────────────────────────────────────────────────
grid_sub <- grid_sf %>%
  filter(ecoregion %in% unique(sens_df$ecoreg)) %>%
  inner_join(sens_df, by = 'id')
sens_centroids <- st_centroid(grid_sub)
sens_proj      <- st_transform(sens_centroids, crs_use)

cat(sprintf('Spatial join: %d points projected\n', nrow(sens_proj)))

# ── Shared bounding box ────────────────────────────────────────────────────────
bbox     <- st_bbox(sens_proj)
pad      <- 400000
xlim_map <- c(bbox['xmin'] - pad, bbox['xmax'] + pad)
ylim_map <- c(bbox['ymin'] - pad, bbox['ymax'] + pad)

# ── Shared base layers ────────────────────────────────────────────────────────
base_layers <- list(
  geom_sf(data = world_proj,  fill = 'grey55',    color = 'grey40', linewidth = 0.2),
  geom_sf(data = ecoreg_proj, fill = 'lightgrey', color = 'grey20', linewidth = 0.5, linetype = 'solid'),
  theme_void(base_size = 12),
  theme(
    panel.background = element_rect(fill = '#e0edf5', color = NA),
    legend.position  = 'bottom',
    legend.title     = element_text(size = 10, face = 'bold'),
    plot.title       = element_text(size = 12, face = 'bold', hjust = 0.5)
  )
)

# ── Map: temperature sensitivity (regression slope) ────────────────────────────
# Flag statistically significant cells (p < 0.05) with larger points
sens_proj <- sens_proj %>% mutate(sig = p_val < 0.05)

slope_lim <- quantile(abs(sens_proj$slope), 0.98, na.rm = TRUE)  # robust colour limits

p_slope <- ggplot() +
  base_layers +
  geom_sf(data  = sens_proj,
          aes(color = slope, size = sig),
          alpha = 0.9, shape = 19) +
  scale_color_gradient2(
    low      = 'steelblue3',
    mid      = 'white',
    high     = 'red3',
    midpoint = 0,
    limits   = c(-slope_lim, slope_lim),
    oob      = scales::squish,
    name     = 'Sensitivity\n(days / °C·day)',
    guide    = guide_colorbar(barwidth = 8, barheight = 0.8)
  ) +
  scale_size_manual(
    values = c('FALSE' = 0.8, 'TRUE' = 2.0),
    guide  = 'none'
  ) +
  coord_sf(crs = crs_use, xlim = xlim_map, ylim = ylim_map) +
  labs(title = 'GUP anomaly sensitivity to preseason AGDD anomaly\n(larger points: p < 0.05)')

# ── Map: regression R² ────────────────────────────────────────────────────────
p_rsq <- ggplot() +
  base_layers +
  geom_sf(data  = sens_proj,
          aes(color = r_sq),
          size  = 1.5, alpha = 0.9, shape = 19) +
  scale_color_distiller(
    palette   = 'YlOrRd',
    direction = 1,
    limits    = c(0, 1),
    name      = 'R²',
    guide     = guide_colorbar(barwidth = 8, barheight = 0.8)
  ) +
  coord_sf(crs = crs_use, xlim = xlim_map, ylim = ylim_map) +
  labs(title = 'Regression R²\n(GUP anomaly ~ AGDD anomaly)')

# ── Save ───────────────────────────────────────────────────────────────────────
ecoreg_tag <- paste(ecoregs, collapse = '_')

ggsave(file.path(fig_dir, paste0('AGDD_GUP_slope_', ecoreg_tag, '.png')),
       plot = p_slope, width = 10, height = 8, dpi = 300, bg = 'white')

p_combined <- p_slope + p_rsq +
  plot_annotation(
    title = 'Temperature sensitivity of GUP timing via preseason AGDD anomalies',
    theme = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5))
  )

ggsave(file.path(fig_dir, paste0('AGDD_GUP_sensitivity_', ecoreg_tag, '.png')),
       plot = p_combined, width = 18, height = 8, dpi = 300, bg = 'white')

# ── Diagnostic scatter plots: every 10th cell ─────────────────────────────────
# Each panel shows the GUP_anom ~ AGDD_anom cloud with the fitted lm line,
# R², slope, and n labelled in the strip. Use these to visually confirm the
# relationship and check whether the cell-level anomaly approach is working.
if (length(scatter_data_list) > 0) {
  scatter_df <- bind_rows(scatter_data_list)

  # Fix panel order so facets appear in processing sequence
  scatter_df$panel_label <- factor(scatter_df$panel_label,
                                   levels = unique(scatter_df$panel_label))

  p_scatter <- ggplot(scatter_df, aes(x = AGDD_anom, y = GUP_anom)) +
    geom_point(alpha = 0.4, size = 0.9, color = 'steelblue4') +
    geom_smooth(method = 'lm', se = TRUE, color = 'red3',
                linewidth = 0.8, fill = 'red3', alpha = 0.15) +
    geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.3, color = 'grey50') +
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.3, color = 'grey50') +
    facet_wrap(~ panel_label, scales = 'free') +
    labs(
      x        = 'AGDD anomaly (°C·days)',
      y        = 'GUP anomaly (days)',
      title    = 'GUP anomaly vs. preseason AGDD anomaly — sample of grid cells (every 10th)',
      subtitle = 'Each panel: one grid cell. Points = pixel-years. Anomalies relative to cell-wide mean.'
    ) +
    theme_bw(base_size = 9) +
    theme(
      strip.text       = element_text(size = 7),
      strip.background = element_rect(fill = 'grey92'),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = 'bold', size = 11),
      plot.subtitle    = element_text(size = 9, color = 'grey40')
    )

  # Scale figure height to number of panels (4 columns)
  n_panels   <- length(unique(scatter_df$panel_label))
  n_rows_fig <- ceiling(n_panels / 4)
  fig_h      <- max(4, n_rows_fig * 3.2)

  ggsave(file.path(fig_dir, paste0('scatter_GUP_AGDD_sample_', ecoreg_tag, '.png')),
         plot   = p_scatter,
         width  = 16,
         height = fig_h,
         dpi    = 300,
         bg     = 'white')
  cat(sprintf('Scatter diagnostics saved (%d panels).\n', n_panels))
}

cat('Done. Figures saved to', fig_dir, '\n')
