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
  extracted_data <- terra::extract(r, v_proj, cells=TRUE)
  
  # Get unique cell indices
  overlapping_indices <- unique(extracted_data$cell)
  
  # Sample numPix locations
  sampled_indices <- sample(overlapping_indices, min(length(overlapping_indices), numPix))
  
  # # Plotting for verification
  # # Create a small buffer so you can see the surroundings
  # plot_ext <- ext(v_proj) + 5000 # Add 5km padding for context
  # # 2. Crop the raster to that area (fast and light)
  # r_crop <- crop(r, plot_ext)
  # # 3. Plot the results
  # plot(r_crop, y=1, main = paste("Verification: Grid Cell", i))
  # plot(v_proj, add = TRUE, border = "red", lwd = 2) # The 5km grid cell
  # # 4. Plot the sampled pixels to see their locations
  # # xyFromCell converts the indices back to coordinates
  # sampled_coords <- xyFromCell(r, sampled_indices)
  # points(sampled_coords, col = "black", pch = 3, cex = 0.5)

  return(sampled_indices)
}