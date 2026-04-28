library(terra)
library(sf)
library(dplyr)

######################################################################
## Script for preparing the JAGS inputs for each ecoregion          ##
######################################################################

# 0) Args
ecoreg <- '2.2.2'
total_nsamples <- 10000
years <- 2001:2024
Tavg_path <- paste0('/projectnb/modislc/data/climate/CHELSA/arctic/MODIS_gridded/', ecoreg, '/', 'tas_', ecoreg, '_2001.tif')
SWin_path <- paste0('/projectnb/modislc/data/climate/CHELSA/arctic/MODIS_gridded/', ecoreg, '/', 'rsds_', ecoreg, '_2001.tif')
Pprd_path <- paste0('/projectnb/modislc/data/climate/CHELSA/arctic/MODIS_gridded/', ecoreg, '/', 'photoperiod_', ecoreg, '_2013.tif')
LSP_path <-  paste0('/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/2001/mod_', ecoreg, '_2001_pheno.tif')
# Load functions
source(file='/projectnb/modislc/users/seamorez/HLS_Pheno/scripts/LSP_clim_sensitivity/LSPCS_Functions.R')

# 1) Sample grid cells for the ecoregion
# Read in the pan-Arctic 5km grid and filter to the ecoregion
panArctic_grid_5km <- st_read("/projectnb/modislc/users/seamorez/HLS_Pheno/shapefiles/Arctic_grid_5km.shp")
ecoreg_grid_5km <- panArctic_grid_5km %>% filter(ecoregion == ecoreg)

# Calculate the n samples for the ecoregion and take the random sample of grid cells
nsamples <- round(total_nsamples * (nrow(ecoreg_grid_5km)/nrow(panArctic_grid_5km)))
set.seed(7)  # Seed for repeatable results (same sample every time)
ecoreg_samples_5km <- ecoreg_grid_5km %>% slice_sample(n = nsamples)

# 2) For each grid cell, clip LSP, climate rasters to the grid cell, randomly sample 50 MODIS pixel locations in each, and then randomly sample 100 px-yrs from the time series data
# Sample 50 MODIS pixel locations, using Tavg data as a basemap
for (i in 1:2) {#nrow(ecoreg_samples_5km)) {  # UNCOMMENT THIS WHEN STEP 2 IS FULLY READY AND OPERATIONAL!
  cell_id <- ecoreg_samples_5km[i, ]$id
  px_locations <- samplePixelsInGrid(Tavg_path, ecoreg_samples_5km[i, ], 50)
  print(px_locations)
  
  if (is.null(px_locations)) {
    cat(sprintf('Grid cell %d: no valid pixels — skipping\n', cell_id))
    next
  }
  
  # Sample 100 px-yrs from the time series data
  # Create a dataframe of all pixel-years
  combinations <- expand.grid(px_locations, years)
  # Randomly sample 100 indices without replacement
  n <- 100
  sample_indices <- sample(1:nrow(combinations), n, replace = FALSE)
  # Retrieve the unique tuples
  sampled_tuples <- combinations[sample_indices, ]
  print(sampled_tuples)
  
  # For each sampled year, retrieve the daily time series for Tavg, SWin, Pprd, and GUP
  unique_yrs <- sort(unique(sampled_tuples$Var2))
  cat(sprintf('Grid cell %d: %d unique years across %d sampled tuples\n',
              cell_id, length(unique_yrs), nrow(sampled_tuples)))
  # Collect per-year data.frames in a list (much faster than rbind in loop)
  yr_dat_list <- vector('list', length(unique_yrs))
  for (yr_idx in seq_along(unique_yrs)) {
    yr <- unique_yrs[yr_idx]
    
    # Find all pixel locations sampled in this year
    # sampled_tuples columns: px_locations (cell numbers), years (int)
    yr_rows    <- sampled_tuples$Var2 == yr
    yr_px_locs <- sampled_tuples$Var1[yr_rows]
    
    cat(sprintf('  Year %d | %d pixels | %d days each | expected %d rows\n',
                yr,
                length(yr_px_locs),
                if (is_leap_year(yr)) 366L else 365L,
                length(yr_px_locs) * (if (is_leap_year(yr)) 366L else 365L)))
    
    # Extract daily Tavg, SWin, Pprd, and PS for all pixels in this year
    yr_dat_list[[yr_idx]] <- extract_year_timeseries(
      ecoreg       = ecoreg,
      grid_cell_id = cell_id,
      px_locs      = yr_px_locs,
      year         = yr,
      gup_band     = 1    # band 1 = GUP DOY in mod_{ecoreg}_{year}_pheno.tif
    )
  }

  # Combine all years into one data.frame for this grid cell
  # bind_rows handles NULL entries from missing files gracefully
  dat <- bind_rows(yr_dat_list)
  
  # For GUP model, we only need DOYs 1-250
  ndoys_model <- 250L   # consistent window for JAGS
  dat <- dat %>%
    filter(doy <= ndoys_model) %>%   # truncate to model window
    arrange(px_location, doy)        # ensure sequential block order
  # Verify structure is exactly n_pix × ndoys_model rows
  stopifnot(nrow(dat) == n * ndoys_model)
  
  # # Drop any px_location that has NA in ANY climate variable across ANY of its rows
  # bad_px <- dat %>%
  #   group_by(px_location) %>%
  #   filter(any(is.na(Tavg) | is.na(SWin) | is.na(Pprd))) %>%
  #   pull(px_location) %>%
  #   unique()
  # if (length(bad_px) > 0) {
  #   cat(sprintf('  Dropping %d px_locations with NA climate data\n', length(bad_px)))
  #   dat <- dat %>% filter(!px_location %in% bad_px)
  # }
  
  cat(sprintf('Grid cell %d complete: %d rows (%d pixel-years)\n',
              cell_id, nrow(dat), nrow(sampled_tuples)))
  
  print(summary(dat[, c('Tavg', 'SWin', 'Pprd', 'PS')]))
  
}

#jags_input/[ecoreg]/[ecoreg]_[gridcell].rda
