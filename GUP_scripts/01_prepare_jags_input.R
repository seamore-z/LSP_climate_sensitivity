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
SWin_path <- paste0('/projectnb/modislc/data/climate/CHELSA/arctic/MODIS_gridded/', ecoreg, '/', 'rsps_', ecoreg, '_2001.tif')
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
for (i in 1:1) {#nrow(ecoreg_samples_5km)) {  # UNCOMMENT THIS WHEN STEP 2 IS FULLY READY AND OPERATIONAL!
  px_locations <- samplePixelsInGrid(Tavg_path, ecoreg_samples_5km[i, ], 50)
  print(px_locations)
  
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
}

#jags_input/filt/[tile]/[tile]_[panno].rda
