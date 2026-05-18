library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(data.table)

# ── Configuration ─────────────────────────────────────────────────────────────
ecoreg     <- '2.2.2'
output_dir <- paste0('/projectnb/modislc/users/seamorez/HLS_Pheno/GUP_climate_sensitivity/jags_output/', ecoreg)
grid_shp   <- '/projectnb/modislc/users/seamorez/HLS_Pheno/shapefiles/Arctic_grid_5km.shp'
fig_dir    <- '/projectnb/modislc/users/seamorez/HLS_Pheno/GUP_climate_sensitivity/figures'
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# beta column indices (in X matrix order: intercept, Tavg, Pprd, SWin)
idx_tavg <- 2   # beta[,2] = Tavg
idx_pprd <- 3   # beta[,3] = Pprd
idx_swin <- 4   # beta[,4] = SWin

# ── Load grid shapefile ────────────────────────────────────────────────────────
grid_sf <- st_read(grid_shp) %>%
  filter(ecoregion == ecoreg)
cat(sprintf('Grid cells in shapefile for %s: %d\n', ecoreg, nrow(grid_sf)))

# ── Loop through rda output files ─────────────────────────────────────────────
rda_files <- list.files(
  output_dir,
  pattern   = paste0('^ma_', ecoreg, '_.*\\.rda$'),
  full.names = TRUE
)
cat(sprintf('Found %d rda output files\n', length(rda_files)))

if (length(rda_files) == 0) stop('No rda files found — check output_dir and filename pattern')

results_list <- vector('list', length(rda_files))

for (k in seq_along(rda_files)) {
  f    <- rda_files[k]
  base <- basename(f)
  
  # Parse cell_id from filename: ma_2.2.2_668891.rda → 668891
  cell_id_str <- sub(paste0('^ma_', ecoreg, '_'), '', sub('\\.rda$', '', base))
  cell_id     <- as.integer(cell_id_str)
  
  if (is.na(cell_id)) {
    cat(sprintf('  WARNING: could not parse cell_id from %s — skipping\n', base))
    next
  }
  
  # Load into isolated environment to avoid overwriting global objects
  env <- new.env(parent = emptyenv())
  tryCatch(
    load(f, envir = env),
    error = function(e) cat(sprintf('  ERROR loading %s: %s\n', base, e$message))
  )
  
  if (!exists('gibbs', envir = env)) {
    cat(sprintf('  WARNING: no gibbs object in %s — skipping\n', base))
    next
  }
  
  gibbs <- env$gibbs
  
  # Relative thermal forcing per MCMC sample (excludes intercept)
  # Each row = one MCMC iteration
  b_tavg <- abs(gibbs$beta[, idx_tavg])
  b_pprd <- abs(gibbs$beta[, idx_pprd])
  b_swin <- abs(gibbs$beta[, idx_swin])
  b_sum  <- b_tavg + b_pprd + b_swin
  
  rel_thermal <- b_tavg / b_sum   # vector of length n_sims
  
  # Also store relative Pprd and SWin for potential later use
  rel_pprd <- b_pprd / b_sum
  rel_swin <- b_swin / b_sum
  
  results_list[[k]] <- data.frame(
    id              = cell_id,
    # Thermal forcing
    rel_thermal     = median(rel_thermal),
    rel_thermal_lo  = quantile(rel_thermal, 0.025),
    rel_thermal_hi  = quantile(rel_thermal, 0.975),
    # Photoperiod forcing
    rel_pprd        = median(rel_pprd),
    # SWin forcing
    rel_swin        = median(rel_swin),
    # Raw beta medians for reference
    beta_tavg       = median(gibbs$beta[, idx_tavg]),
    beta_pprd       = median(gibbs$beta[, idx_pprd]),
    beta_swin       = median(gibbs$beta[, idx_swin]),
    beta_kappa      = median(gibbs$kappa)
  )
  
  cat(sprintf('  [%d/%d] cell %d | rel_thermal = %.3f [%.3f, %.3f]\n',
              k, length(rda_files), cell_id,
              results_list[[k]]$rel_thermal,
              results_list[[k]]$rel_thermal_lo,
              results_list[[k]]$rel_thermal_hi))
}

results_df <- bind_rows(results_list)
cat(sprintf('\nProcessed %d / %d grid cells successfully\n',
            nrow(results_df), length(rda_files)))

# ── Join with spatial grid (use centroids for dot plot) ───────────────────────
grid_joined <- grid_sf %>%
  inner_join(results_df, by = 'id')

# Check join
n_matched <- nrow(grid_joined)
n_missing <- nrow(results_df) - n_matched
cat(sprintf('Spatial join: %d matched, %d unmatched (id not in shapefile)\n',
            n_matched, n_missing))

# Convert polygons to centroids for dot-style plot
grid_centroids <- st_centroid(grid_joined)

# ── Base map ───────────────────────────────────────────────────────────────────
world <- ne_countries(scale = 'medium', returnclass = 'sf')

# Arctic Polar Stereographic (EPSG:3995) — clean for pan-Arctic view
# crs_use <- st_crs(3995)
crs_use <- st_crs("ESRI:102001")   # Canada Albers Equal Area Conic

world_proj     <- st_transform(world, crs_use)
points_proj    <- st_transform(grid_centroids, crs_use)

# Bounding box: derive from data with padding
bbox     <- st_bbox(points_proj)
pad      <- 500000   # 500 km padding around data extent
xlim_map <- c(bbox['xmin'] - pad, bbox['xmax'] + pad)
ylim_map <- c(bbox['ymin'] - pad, bbox['ymax'] + pad)

# ── Plot 1: Relative thermal forcing ──────────────────────────────────────────
p_thermal <- ggplot() +
  geom_sf(data  = world_proj,
          fill  = 'grey55',
          color = 'grey40',
          linewidth = 0.2) +
  geom_sf(data  = points_proj,
          aes(color = rel_thermal),
          size  = 1.5,
          alpha = 0.9) +
  scale_color_gradient(
    low    = 'white',
    high   = 'red3',
    limits = c(0, 1),
    breaks = c(0.1, 0.5, 0.9),
    labels = c('< 0.1', '0.5', '> 0.9'),
    name   = 'Thermal\nforcing'
  ) +
  coord_sf(
    crs  = crs_use,
    xlim = xlim_map,
    ylim = ylim_map
  ) +
  theme_void(base_size = 12) +
  theme(
    panel.background  = element_rect(fill = 'grey80', color = NA),
    legend.position   = 'left',
    legend.title      = element_text(size = 10, face = 'bold'),
    legend.text       = element_text(size = 9),
    plot.title        = element_text(size = 12, face = 'bold', hjust = 0.5),
    plot.margin       = margin(5, 5, 5, 5)
  ) +
  labs(title = paste0('GUP thermal forcing dependence — ecoregion ', ecoreg))

print(p_thermal)
ggsave(
  file.path(fig_dir, paste0('rel_thermal_', ecoreg, '.png')),
  plot   = p_thermal,
  width  = 10, height = 8, dpi = 300, bg = 'white'
)

# ── Plot 2: Three-panel all drivers side by side ───────────────────────────────
# Helper to build one panel
make_panel <- function(var, title, color_high) {
  ggplot() +
    geom_sf(data = world_proj, fill = 'grey55', color = 'grey40', linewidth = 0.2) +
    geom_sf(data  = points_proj,
            aes(color = .data[[var]]),
            size  = 1, alpha = 0.9) +
    scale_color_gradient(
      low    = 'white',
      high   = color_high,
      limits = c(0, 1),
      breaks = c(0.1, 0.5, 0.9),
      labels = c('< 0.1', '0.5', '> 0.9'),
      name   = 'Rel.\ndependence'
    ) +
    coord_sf(crs = crs_use, xlim = xlim_map, ylim = ylim_map) +
    theme_void(base_size = 11) +
    theme(
      panel.background = element_rect(fill = 'grey80', color = NA),
      legend.position  = 'bottom',
      plot.title       = element_text(size = 11, face = 'bold', hjust = 0.5)
    ) +
    labs(title = title)
}

p1 <- make_panel('rel_thermal', 'Thermal forcing (Tavg)',    'red3')
p2 <- make_panel('rel_pprd',    'Photoperiod (Pprd)',        'steelblue3')
p3 <- make_panel('rel_swin',    'Shortwave radiation (SWin)', 'darkorange2')

library(patchwork)
p_all <- p1 + p2 + p3 +
  plot_annotation(
    title   = paste0('Relative climate driver dependence of GUP — ecoregion ', ecoreg),
    theme   = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5))
  )

print(p_all)
ggsave(
  file.path(fig_dir, paste0('rel_all_drivers_', ecoreg, '.png')),
  plot   = p_all,
  width  = 18, height = 7, dpi = 300, bg = 'white'
)

# ── Save results table ─────────────────────────────────────────────────────────
# write.csv(
#   results_df,
#   file.path(fig_dir, paste0('rel_forcing_', ecoreg, '.csv')),
#   row.names = FALSE
# )
# cat('Done. Figures and CSV saved to', fig_dir, '\n')

## Plotting pre-season length
hmax      <- 100
h_thresh  <- 0.01 * hmax   # h must exceed this to count as "forcing started" (1% of hmax)

# ── Load grid shapefile ────────────────────────────────────────────────────────
grid_sf <- st_read(grid_shp) %>% filter(ecoregion == ecoreg)

# ── Find rda files ─────────────────────────────────────────────────────────────
rda_files <- list.files(
  output_dir,
  pattern    = paste0('^ma_', ecoreg, '_.*\\.rda$'),
  full.names = TRUE
)
cat(sprintf('Found %d rda files\n', length(rda_files)))
if (length(rda_files) == 0) stop('No rda files found — check output_dir and pattern')

# ── Loop: extract preseason stats per grid cell ────────────────────────────────
results_list <- vector('list', length(rda_files))

for (k in seq_along(rda_files)) {
  f    <- rda_files[k]
  base <- basename(f)
  
  # Parse cell_id
  cell_id_str <- sub(paste0('^ma_', ecoreg, '_'), '', sub('\\.rda$', '', base))
  cell_id     <- as.integer(cell_id_str)
  if (is.na(cell_id)) {
    cat(sprintf('  WARNING: could not parse cell_id from %s — skipping\n', base))
    next
  }
  
  # Load into isolated environment
  env <- new.env(parent = emptyenv())
  tryCatch(
    load(f, envir = env),
    error = function(e) cat(sprintf('  ERROR loading %s: %s\n', base, e$message))
  )
  
  # Check required objects exist
  if (!all(c('gibbs', 'dat', 'data') %in% ls(env))) {
    cat(sprintf('  WARNING: missing objects in %s — skipping\n', base))
    next
  }
  
  gibbs     <- env$gibbs
  dat_cell  <- env$dat
  data_cell <- env$data
  
  # Check h was saved
  if (is.null(gibbs$h)) {
    cat(sprintf('  WARNING: h not in gibbs for cell %d — was it saved? skipping\n', cell_id))
    next
  }
  
  # ── Summarise h and Y into per-block preseason stats ────────────────────────
  h_med <- apply(gibbs$h, 2, median)   # posterior median h, length = n observations
  
  compare <- data.table(
    px_locs = dat_cell$px_location,
    year    = dat_cell$year,
    doy     = dat_cell$doy,
    y       = data_cell$Y,             # observed binary PS
    h       = h_med
  )
  
  # Per block (px_locs × year): forcing start DOY and observed onset DOY
  block_stats <- compare[, {
    # First DOY where h meaningfully lifts off zero
    force_idx <- which(h > h_thresh)
    doy_force <- if (length(force_idx) > 0) doy[force_idx[1]] else NA_integer_
    
    # First DOY where observed green-up occurred
    onset_idx <- which(y == 1)
    doy_onset <- if (length(onset_idx) > 0) doy[onset_idx[1]] else NA_integer_
    
    list(
      doy_force  = doy_force,
      doy_onset  = doy_onset,
      preseason  = doy_onset - doy_force   # NA if either is missing
    )
  }, by = .(px_locs, year)]
  
  # Cell-level summaries (median across all px × year blocks)
  n_blocks_total  <- nrow(block_stats)
  n_blocks_valid  <- sum(!is.na(block_stats$preseason))
  
  results_list[[k]] <- data.frame(
    id                  = cell_id,
    preseason_median    = median(block_stats$preseason,  na.rm = TRUE),
    preseason_mean      = mean(block_stats$preseason,    na.rm = TRUE),
    preseason_sd        = sd(block_stats$preseason,      na.rm = TRUE),
    onset_median        = median(block_stats$doy_onset,  na.rm = TRUE),
    force_start_median  = median(block_stats$doy_force,  na.rm = TRUE),
    n_blocks_valid      = n_blocks_valid,
    n_blocks_total      = n_blocks_total
  )
  
  cat(sprintf(
    '  [%d/%d] cell %d | preseason = %.1f days (median) | onset = DOY %.1f | force start = DOY %.1f | valid blocks: %d/%d\n',
    k, length(rda_files), cell_id,
    results_list[[k]]$preseason_median,
    results_list[[k]]$onset_median,
    results_list[[k]]$force_start_median,
    n_blocks_valid, n_blocks_total
  ))
}

results_df <- bind_rows(results_list)
cat(sprintf('\nProcessed %d / %d cells\n', nrow(results_df), length(rda_files)))

# ── Spatial join ───────────────────────────────────────────────────────────────
grid_joined <- grid_sf %>%
  inner_join(results_df, by = 'id')

cat(sprintf('Matched %d / %d cells to shapefile\n', nrow(grid_joined), nrow(results_df)))

# Centroids for dot plot
grid_centroids <- st_centroid(grid_joined)

# ── Project ────────────────────────────────────────────────────────────────────
crs_use        <- st_crs("ESRI:102001")
world          <- ne_countries(scale = 'medium', returnclass = 'sf')
world_proj     <- st_transform(world, crs_use)
points_proj    <- st_transform(grid_centroids, crs_use)

bbox     <- st_bbox(points_proj)
pad      <- 300000
xlim_map <- c(bbox['xmin'] - pad, bbox['xmax'] + pad)
ylim_map <- c(bbox['ymin'] - pad, bbox['ymax'] + pad)

# ── Plot 1: Preseason length ───────────────────────────────────────────────────
p_preseason <- ggplot() +
  geom_sf(data = world_proj, fill = 'grey55', color = 'grey40', linewidth = 0.2) +
  geom_sf(data  = points_proj,
          aes(color = preseason_median),
          size  = 3, alpha = 0.9, shape = 19) +
  scale_color_distiller(
    palette  = 'RdYlBu',
    direction = -1,           # red = long preseason, blue = short
    limits   = c(20, 50),
    name     = 'Preseason\nlength (days)',
    guide    = guide_colorbar(barwidth = 8, barheight = 0.8)
  ) +
  coord_sf(crs = crs_use, xlim = xlim_map, ylim = ylim_map) +
  theme_void(base_size = 12) +
  theme(
    panel.background = element_rect(fill = 'grey80', color = NA),
    legend.position  = 'bottom',
    legend.title     = element_text(size = 10, face = 'bold'),
    plot.title       = element_text(size = 12, face = 'bold', hjust = 0.5)
  ) +
  labs(title = paste0('Preseason length — ecoregion ', ecoreg))

print(p_preseason)
ggsave(
  file.path(fig_dir, paste0('preseason_', ecoreg, '.png')),
  plot  = p_preseason,
  width = 10, height = 8, dpi = 300, bg = 'white'
)