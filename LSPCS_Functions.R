#---------------------------------------------------------------------
#Returns the indices for the sampled pixels
#---------------------------------------------------------------------
samplePixelsInGrid <-function(rastPath, gridCell, numPix) {
  r <- rast(rastPath)

  if (st_crs(r) == st_crs(gridCell)) {
    message("CRS already matches.")
  } else {
    message("CRS mismatch. Reprojecting vector to match raster...")
    # 1. Convert sf to terra's SpatVector
    v <- vect(gridCell)

    # 2. Use terra::project to ALIGN to the raster's specific Sinusoidal sphere
    # This avoids the 'blown up' polygon issue common with st_transform
    v_proj <- project(v, r)
    
    # 3. Double-check the size before extraction
    # For 5km grid and 500m pixels, ext_area should be ~25,000,000
    print(paste("Projected Polygon Area:", expanse(v_proj)))
  }
  
  # Extract raster cells (returns index in 'id' column)
  extracted_data <- terra::extract(r[[1]], v_proj, cells=TRUE)
  
  # Get unique cell indices
  overlapping_indices <- unique(extracted_data$cell)
  
  # Drop NA pixels and get valid cell indices in one step
  valid_cells <- extracted_data$cell[!is.na(extracted_data[, 2])]
  if (length(valid_cells) == 0) return(NULL)
  
  # Sample numPix locations
  sampled_indices <- sample(overlapping_indices, min(length(overlapping_indices), numPix))
  
  # # Plotting for verification
  # Create a small buffer so you can see the surroundings
  plot_ext <- ext(v_proj) + 5000 # Add 5km padding for context
  # 2. Crop the raster to that area (fast and light)
  r_crop <- crop(r, plot_ext)
  # 3. Plot the results
  plot(r_crop, y=1, main = paste("Verification: Grid Cell", i))
  plot(v_proj, add = TRUE, border = "red", lwd = 2) # The 5km grid cell
  # 4. Plot the sampled pixels to see their locations
  # xyFromCell converts the indices back to coordinates
  sampled_coords <- xyFromCell(r, sampled_indices)
  points(sampled_coords, col = "black", pch = 3, cex = 0.5)

  return(sampled_indices)
}


#---------------------------------------------------------------------
#Check if a year is a leap year
#---------------------------------------------------------------------
#
is_leap_year <- function(year) {
  (year %% 4 == 0) & ((year %% 100 != 0) | (year %% 400 == 0))
}


#---------------------------------------------------------------------
#Build file paths for all four rasters for a given ecoregion and year.
#Photoperiod uses 2013 (non-leap) or 2012 (leap) as representative years
#since daylength is identical across all years of the same type.
#---------------------------------------------------------------------
build_raster_paths <- function(ecoreg, year) {
  base_clim <- paste0('/projectnb/modislc/data/climate/CHELSA/arctic/MODIS_gridded/',
                      ecoreg, '/')
  base_lsp  <- paste0('/projectnb/modislc/users/seamorez/HLS_Pheno/modis_data/',
                      year, '/')
  pprd_year <- if (is_leap_year(year)) 2012 else 2013
  
  list(
    Tavg = paste0(base_clim, 'tas_',         ecoreg, '_', year,      '.tif'),
    SWin = paste0(base_clim, 'rsds_',        ecoreg, '_', year,      '.tif'),
    Pprd = paste0(base_clim, 'photoperiod_', ecoreg, '_', pprd_year, '.tif'),
    LSP  = paste0(base_lsp,  'mod_',         ecoreg, '_', year,      '_pheno.tif')
  )
}


#---------------------------------------------------------------------
#Extract full daily time series for a set of sampled px in one year.
#
#Loads all four rasters ONCE per year, extracts all sampled px in a
#single vectorised pass, and returns a long-format dataframe with
#365 or 366 rows per pixel.
#
#' @return data.frame with n_pixels x ndoys rows and columns:
#'         ecoregion | grid_cell | px_location | year | doy |
#'         Tavg (K) | SWin (W m-2) | Pprd (s) | PS (0/1/NA)
#---------------------------------------------------------------------
extract_year_timeseries <- function(ecoreg, grid_cell_id, px_locs, year,
                                    gup_band = 1) {
  
  ndoys <- if (is_leap_year(year)) 366L else 365L
  n_pix <- length(px_locs)
  paths <- build_raster_paths(ecoreg, year)
  
  # ── Validate all four paths before loading to avoid partial reads ─────────
  missing_files <- names(paths)[!file.exists(unlist(paths))]
  if (length(missing_files) > 0) {
    warning(sprintf('Year %d ecoregion %s: missing rasters [%s] — skipping year',
                    year, ecoreg, paste(missing_files, collapse = ', ')))
    return(NULL)
  }
  
  # ── Load rasters ──────────────────────────────────────────────────────────
  r_tavg <- rast(paths$Tavg)   # 365/366 bands, stored as int16 Kelvin × 10
  r_swin <- rast(paths$SWin)   # 365/366 bands, stored as int16 W m-2 × 10
  r_pprd <- rast(paths$Pprd)   # 365/366 bands, stored as int32 seconds (no scaling)
  r_lsp  <- rast(paths$LSP)    # multi-band phenology; band gup_band = GUP DOY
  
  # ── Extract values at all sampled pixel locations (one call per raster) ───
  # terra::extract returns data.frame: rows = pixels, cols = bands
  # ID = FALSE drops the auto-generated ID column terra adds by default
  gup_df  <- terra::extract(r_lsp[[gup_band]], px_locs, ID = FALSE)  # n_pix × 1
  tavg_df <- terra::extract(r_tavg,            px_locs, ID = FALSE)  # n_pix × ndoys
  swin_df <- terra::extract(r_swin,            px_locs, ID = FALSE)  # n_pix × ndoys
  pprd_df <- terra::extract(r_pprd,            px_locs, ID = FALSE)  # n_pix × ndoys

              
  # ── Convert stored int16 values back to physical units ────────────────────
  # tas:  stored as Kelvin × 10 → divide by 10 to recover Kelvin
  # rsds: stored as W m-2 × 10  → divide by 10 to recover W m-2
  # pprd: stored as raw seconds (int32) → no conversion needed
  gup_vec  <- as.integer(gup_df[, 1])       # GUP DOY for each pixel
  tavg_mat <- as.matrix(tavg_df) / 10       # n_pix × ndoys, Kelvin
  swin_mat <- as.matrix(swin_df) / 10       # n_pix × ndoys, W m-2
  pprd_mat <- as.matrix(pprd_df)            # n_pix × ndoys, seconds

  # Mask implausible GUP fill values (e.g. 32767, 0) to NA
  gup_vec[gup_vec < 1 | gup_vec > ndoys] <- NA_integer_
  
  # ── Build long-format data.frame (fully vectorised — no loops) ────────────
  #
  # rep(px_locs, each=ndoys):  px1 px1 ... px1 px2 px2 ... px2 ...
  #                             |---ndoys--|  |---ndoys--|
  # rep(1:ndoys, times=n_pix): 1 2 .. ndoys 1 2 .. ndoys ...
  #
  # as.vector(t(matrix)): flattens matrix row-by-row (pixel 1 all DOYs,
  #                        then pixel 2 all DOYs, etc.) to match the
  #                        pixel ordering from rep(..., each=ndoys)
  px_rep  <- rep(px_locs,          each  = ndoys)
  doy_rep <- rep(seq_len(ndoys),   times = n_pix)
  gup_rep <- rep(gup_vec,          each  = ndoys)  # broadcast GUP across DOYs
  
  tavg_vec <- as.vector(t(tavg_mat))
  swin_vec <- as.vector(t(swin_mat))
  pprd_vec <- as.vector(t(pprd_mat))
  
  # PS: 0 = before GUP (growing season not yet started)
  #     1 = on or after GUP (growing season underway)
  #    NA = GUP missing for this pixel in this year
  ps_vec              <- as.integer(doy_rep >= gup_rep)
  ps_vec[is.na(gup_rep)] <- NA_integer_
  
  data.frame(
    ecoregion   = ecoreg,
    grid_cell   = grid_cell_id,
    px_location = px_rep,
    year        = year,
    doy         = doy_rep,
    Tavg        = tavg_vec,
    SWin        = swin_vec,
    Pprd        = pprd_vec,
    PS          = ps_vec,
    stringsAsFactors = FALSE
  )
}