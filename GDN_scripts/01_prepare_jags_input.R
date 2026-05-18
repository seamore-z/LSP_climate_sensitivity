library(terra)
library(sf)
library(dplyr)
library(data.table)

######################################################################
## Script for preparing JAGS inputs for GDN model                   ##
## Accumulation window: DOY 150-320 for all variables               ##
## GUP_doy derived from full time series before trimming and        ##
## stored as a constant per pixel-year for use in JAGS logit link   ##
######################################################################

# 0) Args
ecoreg         <- '2.2.2'
total_nsamples <- 15000
years          <- 2001:2024
gdn_start      <- 151L   # post-solstice accumulation start
ndoys_full     <- 320L   # extract through DOY 320 to capture GUP timing
ndoys_model    <- ndoys_full - gdn_start + 1L   # 170 DOYs saved per block

Tavg_path <- paste0('/projectnb/modislc/data/climate/CHELSA/arctic/MODIS_gridded/', ecoreg, '/', 'tas_',         ecoreg, '_2001.tif')
SWin_path <- paste0('/projectnb/modislc/data/climate/CHELSA/arctic/MODIS_gridded/', ecoreg, '/', 'rsds_',        ecoreg, '_2001.tif')
Pprd_path <- paste0('/projectnb/modislc/data/climate/CHELSA/arctic/MODIS_gridded/', ecoreg, '/', 'photoperiod_', ecoreg, '_2013.tif')
LSP_path  <- paste0('/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/2001/mod_', ecoreg, '_2001_pheno.tif')

source(file = '/projectnb/modislc/users/seamorez/HLS_Pheno/scripts/LSP_clim_sensitivity/LSPCS_Functions.R')

# 1) Sample grid cells
panArctic_grid_5km <- st_read("/projectnb/modislc/users/seamorez/HLS_Pheno/shapefiles/Arctic_grid_5km.shp")
ecoreg_grid_5km    <- panArctic_grid_5km %>% filter(ecoregion == ecoreg)

nsamples <- round(total_nsamples * (nrow(ecoreg_grid_5km) / nrow(panArctic_grid_5km)))
set.seed(7)
ecoreg_samples_5km <- ecoreg_grid_5km %>% slice_sample(n = nsamples)

# 2) Loop over grid cells
for (i in 1:nrow(ecoreg_samples_5km)) {
  cell_id      <- ecoreg_samples_5km[i, ]$id
  px_locations <- samplePixelsInGrid(LSP_path, ecoreg_samples_5km[i, ], 50)
  
  if (is.null(px_locations)) {
    cat(sprintf('Grid cell %d: no valid pixels — skipping\n', cell_id))
    next
  }
  
  # ── Build valid combinations (px-yr where GDN is not NA) ──────────────
  # Check GDN band directly in pheno TIFs before sampling — avoids the
  # old approach of sampling blindly and then discarding sparse cells
  valid_list <- vector('list', length(years))
  for (yr_idx in seq_along(years)) {
    yr         <- years[yr_idx]
    pheno_path <- paste0('/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/',
                         yr, '/mod_', ecoreg, '_', yr, '_pheno.tif')
    if (!file.exists(pheno_path)) next
    gup_vals <- values(rast(pheno_path)[[3]])[px_locations]  # band 3 = GDN DOY
    valid_px <- px_locations[!is.na(gup_vals)]
    if (length(valid_px) == 0) next
    valid_list[[yr_idx]] <- data.frame(Var1 = valid_px, Var2 = yr)
  }
  valid_combinations <- bind_rows(valid_list)
  
  # ── Skip guards ────────────────────────────────────────────────────────
  combinations <- expand.grid(px_locations, years)
  if (nrow(combinations) < 600) {
    cat(sprintf('  Grid cell %d: only %d total combinations — skipping\n',
                cell_id, nrow(combinations)))
    next
  }
  if (nrow(valid_combinations) < 200) {
    cat(sprintf('  Grid cell %d: only %d valid (non-NA GUP) combinations — skipping\n',
                cell_id, nrow(valid_combinations)))
    next
  }
  
  # ── Sample 100 from valid combinations only ────────────────────────────
  n <- 100
  sample_indices <- sample(1:nrow(valid_combinations), n, replace = FALSE)
  sampled_tuples <- valid_combinations[sample_indices, ]
  
  unique_yrs <- sort(unique(sampled_tuples$Var2))
  cat(sprintf('Grid cell %d: %d unique years across %d sampled tuples\n',
              cell_id, length(unique_yrs), nrow(sampled_tuples)))
  
  yr_dat_list <- vector('list', length(unique_yrs))
  
  for (yr_idx in seq_along(unique_yrs)) {
    yr         <- unique_yrs[yr_idx]
    yr_rows    <- sampled_tuples$Var2 == yr
    yr_px_locs <- sampled_tuples$Var1[yr_rows]
    
    cat(sprintf('  Year %d | %d pixels\n', yr, length(yr_px_locs)))
    
    yr_dat_list[[yr_idx]] <- extract_year_timeseries(
      ecoreg       = ecoreg,
      grid_cell_id = cell_id,
      px_locs      = yr_px_locs,
      year         = yr,
      gup_band     = 1,
      gdn_band     = 3,
      do_gdn       = TRUE
    )
  }
  
  dat <- bind_rows(yr_dat_list)
  
  # Filter to extraction window and sort
  dat <- dat %>%
    filter(doy <= ndoys_full) %>%
    arrange(px_location, year, doy)
  
  # ── Derive GUP_doy BEFORE trimming ───────────────────────────────────
  # GUP occurs pre-DOY 150; after trimming the binary column loses timing
  # information (it's all 1s). Derive here and store as constant column.
  dat <- as.data.table(dat)
  setorder(dat, px_location, year, doy)
  
  # GUP is already the scalar DOY from extract_year_timeseries
  setnames(dat, 'GUP', 'gup_doy')
  
  # Fill any pixel-years with no detected GUP with the cell mean
  if (any(is.na(dat$gup_doy))) {
    gup_mean <- mean(unique(dat$gup_doy), na.rm = TRUE)
    cat(sprintf('  WARNING: %d pixel-years with no GUP — filling with mean (DOY %.0f)\n',
                sum(is.na(unique(dat[, .(px_location, year)])$gup_doy)), gup_mean))
    dat$gup_doy[is.na(dat$gup_doy)] <- gup_mean
  }
  
  # ── Trim to post-solstice accumulation window ─────────────────────────
  dat <- dat[doy >= gdn_start]
  stopifnot(nrow(dat) == n * ndoys_model)
  
  # ── Implausible value removal ─────────────────────────────────────────
  setorder(dat, px_location, year, doy)
  
  cat(sprintf('  Implausible Tavg values: %d\n', sum(dat$Tavg < 173 | dat$Tavg > 400, na.rm = TRUE)))
  cat(sprintf('  Implausible SWin values: %d\n', sum(dat$SWin < 0   | dat$SWin > 900, na.rm = TRUE)))
  dat$Tavg[dat$Tavg < 173 | dat$Tavg > 400] <- NA
  dat$SWin[dat$SWin < 0   | dat$SWin > 900] <- NA
  
  # ── NA rate check ─────────────────────────────────────────────────────
  na_thresh  <- 0.01
  total_rows <- nrow(dat)
  
  na_rates <- c(
    Tavg    = sum(is.na(dat$Tavg))    / total_rows,
    SWin    = sum(is.na(dat$SWin))    / total_rows,
    Pprd    = sum(is.na(dat$Pprd))    / total_rows,
    PS      = sum(is.na(dat$PS))      / total_rows,
    gup_doy = sum(is.na(dat$gup_doy)) / total_rows
  )
  
  cat(sprintf('  NA rates — Tavg: %.2f%%  SWin: %.2f%%  Pprd: %.2f%%  PS: %.2f%%  gup_doy: %.2f%%\n',
              na_rates['Tavg'] * 100, na_rates['SWin'] * 100,
              na_rates['Pprd'] * 100, na_rates['PS'] * 100,
              na_rates['gup_doy'] * 100))
  
  if (any(na_rates >= na_thresh)) {
    cat(sprintf('  Grid cell %d: NA rate exceeds %.0f%% in [%s] — skipping\n',
                cell_id, na_thresh * 100,
                paste(names(na_rates)[na_rates >= na_thresh], collapse = ', ')))
    next
  }
  
  print(summary(dat[, c('Tavg', 'SWin', 'Pprd', 'PS', 'gup_doy')]))
  
  # ── Save ──────────────────────────────────────────────────────────────
  out_path <- paste0('/projectnb/modislc/users/seamorez/HLS_Pheno/GDN_climate_sensitivity/jags_input/', ecoreg)
  if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
  save(dat, file = paste0(out_path, '/', ecoreg, '_', cell_id, '.rda'))
  
  cat(sprintf('Grid cell %d saved: %d rows (%d pixel-years x %d DOYs, DOY %d-%d)\n',
              cell_id, nrow(dat), n, ndoys_model, gdn_start, ndoys_full))
}