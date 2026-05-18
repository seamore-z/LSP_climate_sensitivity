library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(data.table)
library(patchwork)

# ── Configuration ─────────────────────────────────────────────────────────────
ecoregs  <- c('2.1.3','2.1.6','2.2.1','2.2.2','2.2.3','2.2.4','2.3.1','2.4.1', '2.4.2', '2.4.3', '2.4.4')   # all ecoregions to include
grid_shp <- '/projectnb/modislc/users/seamorez/HLS_Pheno/shapefiles/Arctic_grid_5km.shp'
ecoreg_shp <- '/projectnb/modislc/users/seamorez/HLS_Pheno/shapefiles/NAA_ecoregions_final.shp'
fig_dir  <- '/projectnb/modislc/users/seamorez/HLS_Pheno/GUP_climate_sensitivity/figures'
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

idx_tavg <- 2
idx_pprd <- 3
idx_swin <- 4
hmax     <- 100
h_thresh <- 0.01 * hmax

# ── Load spatial layers ────────────────────────────────────────────────────────
grid_sf    <- st_read(grid_shp)
ecoreg_sf  <- st_read(ecoreg_shp)
world      <- ne_countries(scale = 'medium', returnclass = 'sf')
crs_use    <- st_crs("ESRI:102001")

world_proj   <- st_transform(world, crs_use)
ecoreg_proj <- st_transform(
  ecoreg_sf %>% filter(AT_L3CODE %in% ecoregs),
  crs_use
)

# ── Loop over ecoregions — collect all results ─────────────────────────────────
all_beta_list     <- list()
all_preseason_list <- list()

for (ecoreg in ecoregs) {
  
  cat(sprintf('\n══ Processing ecoregion %s ══\n', ecoreg))
  
  output_dir <- paste0('/projectnb/modislc/users/seamorez/HLS_Pheno/GUP_climate_sensitivity/jags_output/', ecoreg)
  
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
  
  beta_list     <- vector('list', length(rda_files))
  preseason_list <- vector('list', length(rda_files))
  
  for (k in seq_along(rda_files)) {
    f    <- rda_files[k]
    base <- basename(f)
    
    cell_id_str <- sub(paste0('^ma_', ecoreg, '_'), '', sub('\\.rda$', '', base))
    cell_id     <- as.integer(cell_id_str)
    if (is.na(cell_id)) next
    
    env <- new.env(parent = emptyenv())
    tryCatch(load(f, envir = env),
             error = function(e) cat(sprintf('  ERROR loading %s\n', base)))
    
    if (!exists('gibbs', envir = env)) next
    gibbs <- env$gibbs
    
    # ── Beta / relative forcing ──────────────────────────────────────────────
    b_tavg      <- abs(gibbs$beta[, idx_tavg])
    b_pprd      <- abs(gibbs$beta[, idx_pprd])
    b_swin      <- abs(gibbs$beta[, idx_swin])
    b_sum       <- b_tavg + b_pprd + b_swin
    rel_thermal <- b_tavg / b_sum
    rel_pprd    <- b_pprd / b_sum
    rel_swin    <- b_swin / b_sum
    
    beta_list[[k]] <- data.frame(
      id          = cell_id,
      ecoreg      = ecoreg,
      rel_thermal = median(rel_thermal),
      rel_thermal_lo = quantile(rel_thermal, 0.025),
      rel_thermal_hi = quantile(rel_thermal, 0.975),
      rel_pprd    = median(rel_pprd),
      rel_swin    = median(rel_swin),
      beta_tavg   = median(gibbs$beta[, idx_tavg]),
      beta_pprd   = median(gibbs$beta[, idx_pprd]),
      beta_swin   = median(gibbs$beta[, idx_swin]),
      beta_kappa  = median(gibbs$kappa)
    )
    
    # ── Preseason ────────────────────────────────────────────────────────────
    if (!all(c('dat', 'data') %in% ls(env)) || is.null(gibbs$h)) next
    
    dat_cell  <- env$dat
    data_cell <- env$data
    h_med     <- apply(gibbs$h, 2, median)
    
    compare <- data.table(
      px_locs = dat_cell$px_location,
      year    = dat_cell$year,
      doy     = dat_cell$doy,
      y       = data_cell$Y,
      h       = h_med
    )
    
    block_stats <- compare[, {
      force_idx <- which(h > h_thresh)
      doy_force <- if (length(force_idx) > 0) doy[force_idx[1]] else NA_integer_
      onset_idx <- which(y == 1)
      doy_onset <- if (length(onset_idx) > 0) doy[onset_idx[1]] else NA_integer_
      list(doy_force = doy_force,
           doy_onset = doy_onset,
           preseason = doy_onset - doy_force)
    }, by = .(px_locs, year)]
    
    preseason_list[[k]] <- data.frame(
      id               = cell_id,
      ecoreg           = ecoreg,
      preseason_median = median(block_stats$preseason,  na.rm = TRUE),
      onset_median     = median(block_stats$doy_onset,  na.rm = TRUE),
      force_start_median = median(block_stats$doy_force, na.rm = TRUE),
      n_blocks_valid   = sum(!is.na(block_stats$preseason)),
      n_blocks_total   = nrow(block_stats)
    )
  }
  
  all_beta_list     <- c(all_beta_list,     beta_list)
  all_preseason_list <- c(all_preseason_list, preseason_list)
  
  cat(sprintf('  Ecoregion %s done: %d beta results, %d preseason results\n',
              ecoreg,
              sum(!sapply(beta_list, is.null)),
              sum(!sapply(preseason_list, is.null))))
}

# ── Combine across ecoregions ──────────────────────────────────────────────────
beta_df     <- bind_rows(all_beta_list)
preseason_df <- bind_rows(all_preseason_list)

cat(sprintf('\nTotal grid cells: %d beta, %d preseason\n',
            nrow(beta_df), nrow(preseason_df)))

# ── Spatial join — all ecoregions ─────────────────────────────────────────────
join_to_centroids <- function(results_df, grid_sf) {
  grid_sub <- grid_sf %>%
    filter(ecoregion %in% unique(results_df$ecoreg)) %>%
    inner_join(results_df, by = 'id')
  st_centroid(grid_sub)
}

beta_centroids     <- join_to_centroids(beta_df,     grid_sf)
preseason_centroids <- join_to_centroids(preseason_df, grid_sf)

beta_proj     <- st_transform(beta_centroids,     crs_use)
preseason_proj <- st_transform(preseason_centroids, crs_use)

cat(sprintf('Spatial join: %d beta points, %d preseason points projected\n',
            nrow(beta_proj), nrow(preseason_proj)))

# ── Shared bounding box across all points ─────────────────────────────────────
all_geom <- c(st_geometry(beta_proj), st_geometry(preseason_proj))
bbox     <- st_bbox(all_geom)
pad      <- 400000
xlim_map <- c(bbox['xmin'] - pad, bbox['xmax'] + pad)
ylim_map <- c(bbox['ymin'] - pad, bbox['ymax'] + pad)

# ── Shared base layers (reused in every plot) ─────────────────────────────────
base_layers <- list(
  geom_sf(data = world_proj,  fill = 'grey55', color = 'grey40', linewidth = 0.2),
  geom_sf(data = ecoreg_proj, fill = 'lightgrey',        color = 'grey20', linewidth = 0.5, linetype = 'solid'),
  theme_void(base_size = 12),
  theme(
    panel.background = element_rect(fill = '#e0edf5', color = NA),
    legend.position  = 'bottom',
    legend.title     = element_text(size = 10, face = 'bold'),
    plot.title       = element_text(size = 12, face = 'bold', hjust = 0.5)
  )
)


# ── Extract latitude from centroids ───────────────────────────────────────────
beta_wgs84 <- st_transform(beta_centroids, 4326)
beta_df$lat <- st_coordinates(beta_wgs84)[, 2]
beta_df$lon <- st_coordinates(beta_wgs84)[, 1]

# ── Reshape for plotting ───────────────────────────────────────────────────────
beta_lat_long <- beta_df %>%
  select(id, ecoreg, lat, rel_thermal, rel_pprd, rel_swin) %>%
  pivot_longer(
    cols      = c(rel_thermal, rel_pprd, rel_swin),
    names_to  = "driver",
    values_to = "rel_dependence"
  ) %>%
  mutate(
    driver = recode(driver,
                    rel_thermal = "Therm.",
                    rel_pprd    = "Photo.",
                    rel_swin    = "SWin."),
    driver = factor(driver, levels = c("Therm.", "Photo.", "SWin."))
  )

# ── Plot: Relative forcings across lat ────────────────────────────────────────
p_lat <- ggplot(beta_lat_long, aes(x = lat, y = rel_dependence, color = driver)) +
  geom_point(alpha = 0.2, size = 0.7) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.2) +
  scale_color_manual(
    values = c("Therm." = "red3", "Photo." = "steelblue3", "SWin." = "yellow3"),
    name   = NULL
  ) +
  labs(
    x     = "Latitude (°N)",
    y     = "Relative dependence",
    title = "GUP climate driver dependence vs. latitude"
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

ggsave(file.path(fig_dir, paste0("rel_dependence_latitude_", ecoreg_tag, ".png")),
       plot = p_lat, width = 8, height = 5, dpi = 300, bg = "white")

# ── Plot: Relative thermal forcing ────────────────────────────────────────────
p_thermal <- ggplot() +
  base_layers +
  geom_sf(data  = beta_proj,
          aes(color = rel_thermal),
          size  = 1.5, alpha = 0.9, shape = 19) +
  # scale_color_gradient(
  #   low    = 'white', high = 'red3',
  #   limits = c(0, 1),
  #   breaks = c(0.1, 0.5, 0.9),
  #   labels = c('< 0.1', '0.5', '> 0.9'),
  #   name   = 'Thermal\nforcing',
  #   guide  = guide_colorbar(barwidth = 8, barheight = 0.8)
  # ) +
  scale_color_distiller(
    palette   = 'Sunset-Sunrise Diverging',
    limits    = c(0, 1),
    name      = 'SWin\nforcing',
    guide     = guide_colorbar(barwidth = 8, barheight = 0.8)
  ) +
  coord_sf(crs = crs_use, xlim = xlim_map, ylim = ylim_map) +
  labs(title = 'GUP thermal forcing dependence')

# ── Plot: Relative photoperiod forcing ────────────────────────────────────────
p_pprd <- ggplot() +
  base_layers +
  geom_sf(data  = beta_proj,
          aes(color = rel_pprd),
          size  = 1.5, alpha = 0.9, shape = 19) +
  # scale_color_gradient(
  #   low    = 'white', high = 'steelblue3',
  #   limits = c(0, 1),
  #   breaks = c(0.1, 0.5, 0.9),
  #   labels = c('< 0.1', '0.5', '> 0.9'),
  #   name   = 'Photoperiod\nforcing',
  #   guide  = guide_colorbar(barwidth = 8, barheight = 0.8)
  # ) +
  scale_color_distiller(
    palette   = 'Cividis',
    limits    = c(0, 1),
    name      = 'SWin\nforcing',
    guide     = guide_colorbar(barwidth = 8, barheight = 0.8)
  ) +
  coord_sf(crs = crs_use, xlim = xlim_map, ylim = ylim_map) + 
  labs(title = 'GUP photoperiod forcing dependence')

# ── Plot: Relative SWin forcing ───────────────────────────────────────────────
p_swin <- ggplot() +
  base_layers +
  geom_sf(data  = beta_proj,
          aes(color = rel_swin),
          size  = 1.5, alpha = 0.9, shape = 19) +
  # scale_color_gradient(
  #   low    = 'white', high = 'darkorange2',
  #   limits = c(0, 1),
  #   breaks = c(0.1, 0.5, 0.9),
  #   labels = c('< 0.1', '0.5', '> 0.9'),
  #   name   = 'SWin\nforcing',
  #   guide  = guide_colorbar(barwidth = 8, barheight = 0.8)
  # ) +
  scale_color_distiller(
    palette   = 'ag_Sunset',
    limits    = c(0, 1),
    name      = 'SWin\nforcing',
    guide     = guide_colorbar(barwidth = 8, barheight = 0.8)
  ) +
  coord_sf(crs = crs_use, xlim = xlim_map, ylim = ylim_map) + 
  labs(title = 'GUP shortwave radiation forcing dependence')

# ── Plot: Preseason length ────────────────────────────────────────────────────
p_preseason <- ggplot() +
  base_layers +
  geom_sf(data  = preseason_proj,
          aes(color = preseason_median),
          size  = 2.5, alpha = 0.9, shape = 19) +
  scale_color_distiller(
    palette   = 'YlGn',
    direction = -1,
    limits    = c(20, 50),
    name      = 'Preseason\nlength (days)',
    guide     = guide_colorbar(barwidth = 8, barheight = 0.8)
  ) +
  coord_sf(crs = crs_use, xlim = xlim_map, ylim = ylim_map) + 
  labs(title = 'GUP preseason length')

# ── Save individual plots ─────────────────────────────────────────────────────
ecoreg_tag <- paste(ecoregs, collapse = '_')

ggsave(file.path(fig_dir, paste0('rel_thermal_',  ecoreg_tag, '.png')),
       plot = p_thermal,  width = 10, height = 8, dpi = 300, bg = 'white')
ggsave(file.path(fig_dir, paste0('preseason_',    ecoreg_tag, '.png')),
       plot = p_preseason, width = 10, height = 8, dpi = 300, bg = 'white')

# ── Save combined three-panel forcing plot ────────────────────────────────────
p_all <- p_thermal + p_pprd + p_swin +
  plot_annotation(
    title = 'Relative climate driver dependence of GUP',
    theme = theme(plot.title = element_text(size = 13, face = 'bold', hjust = 0.5))
  )

ggsave(file.path(fig_dir, paste0('rel_all_drivers_', ecoreg_tag, '.png')),
       plot = p_all, width = 18, height = 7, dpi = 300, bg = 'white')

# ── Save relative importance plot ────────────────────────────────────
# beta_proj <- beta_proj %>%
#   mutate(
#     dominance = (abs(beta_tavg) - abs(beta_pprd)) / (abs(beta_tavg) + abs(beta_pprd)),
#     dom_cat  = cut(dominance, breaks = c(-1, -0.5, -0.2, 0.2, 0.5, 1),
#                    labels = c('Strong Photo.', 'Mod. Photo.', 'Balanced', 'Mod. Therm.', 'Strong Therm.')),
#     size_cat = cut(abs(dominance), breaks = c(0, 0.2, 0.5, 1),
#                    labels = c('Weak', 'Moderate', 'Strong'))
#   )
# 
# p_dominance <- ggplot() +
#   base_layers +
#   geom_sf(data  = beta_proj,
#           aes(color = dom_cat, size = size_cat),
#           alpha = 0.35, shape = 19) +
#   scale_color_manual(
#     values = c(
#       'Strong Photo.' = 'blue',
#       'Mod. Photo.'   = rgb(0.05, 0, 0.95),
#       'Balanced'      = rgb(0.5,  0, 0.5),
#       'Mod. Therm.'   = rgb(0.95, 0, 0.05),
#       'Strong Therm.' = 'red'
#     ),
#     name = NULL
#   ) +
#   scale_size_manual(
#     values = c('Weak' = 2, 'Moderate' = 4, 'Strong' = 6),
#     name   = NULL
#   ) +
#   coord_sf(crs = crs_use, xlim = xlim_map, ylim = ylim_map)
beta_proj <- beta_proj %>%
  mutate(
    dominance     = (abs(beta_tavg) - abs(beta_pprd)) / (abs(beta_tavg) + abs(beta_pprd)),
    dom_magnitude = abs(dominance)
  )

p_dominance <- ggplot() +
  base_layers +
  geom_sf(data  = beta_proj,
          aes(color = dominance, size = dom_magnitude),
          alpha = 0.7, shape = 19) +
  scale_color_gradient2(
    low      = 'blue',
    mid      = 'mediumorchid2',
    high     = 'red',
    midpoint = 0,
    limits   = c(-1, 1),
    name     = 'Dominance\n← Photo.  Therm. →'
  ) +
  scale_size_continuous(
    range  = c(0.5, 5),
    name   = 'Magnitude'
  ) +
  coord_sf(crs = crs_use, xlim = xlim_map, ylim = ylim_map) +
  labs(title = 'Thermal vs. Photoperiod dominance of GUP')


# BOXPLOTS
library(tidyr)

beta_long <- beta_df %>%
  select(id, ecoreg, rel_thermal, rel_pprd, rel_swin) %>%
  mutate(ecoreg = factor(ecoreg, levels = ecoregs)) %>%  # ← your order here
  pivot_longer(
    cols      = c(rel_thermal, rel_pprd, rel_swin),
    names_to  = "driver",
    values_to = "rel_dependence"
  ) %>%
  mutate(
    driver = recode(driver,
                    rel_thermal = "Therm.",
                    rel_pprd    = "Photo.",
                    rel_swin    = "SWin."
    ),
    driver = factor(driver, levels = c("Therm.", "Photo.", "SWin."))
  )

p_box <- ggplot(beta_long, aes(x = driver, y = rel_dependence, fill = ecoreg)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4, width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(
    values = c("2.1.3" = "lavenderblush2", "2.1.6" = "thistle3", "2.2.1" = "slategray2", "2.2.2" = "goldenrod2", "2.2.3" = "turquoise4", "2.2.4" = "salmon2", "2.3.1" = "mistyrose4", "2.4.1" = "navajowhite3", "2.4.2" = "gold1", "2.4.3" = "darkslategray1", "2.4.4" = "palegreen4"),
    # values = c("2.2.1" = "slategray2", "2.2.2" = "goldenrod2", "2.3.1" = "darkseagreen4", "2.2.4" = "salmon2", "2.2.3" = "turquoise4"),
    name   = NULL
  ) +
  labs(x = NULL, y = "Relative dependence") +
  theme_classic(base_size = 24) +
  theme(
    legend.position = "top",
    legend.text      = element_text(size = 20),  # ← legend labels
    axis.text.x      = element_text(size = 24),  # ← Therm. Photo. SWin.
    axis.text.y      = element_text(size = 20),  # ← y-axis numbers
    axis.title.y     = element_text(size = 24)   # ← "Relative dependence"
  )

ggsave(file.path(fig_dir, paste0("rel_dependence_boxplot_", ecoreg_tag, ".png")),
       plot = p_box, width = 6, height = 5, dpi = 300, bg = "white")

# ── Save results tables ───────────────────────────────────────────────────────
# write.csv(beta_df,      file.path(fig_dir, paste0('rel_forcing_',  ecoreg_tag, '.csv')), row.names = FALSE)
# write.csv(preseason_df, file.path(fig_dir, paste0('preseason_',    ecoreg_tag, '.csv')), row.names = FALSE)
# 
# cat('Done. All figures and CSVs saved to', fig_dir, '\n')