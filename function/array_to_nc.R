#' =============================================================================
#' @name array_to_nc
#' @description Write CEPHALOPOD projection arrays to NetCDF.
#'
#'
#' =============================================================================


write_projection_netcdf <- function(array_in,
                                    fn_out,
                                    lon,
                                    lat,
                                    month,
                                    bootstrap,
                                    depth = NULL,
                                    layer_depth_range = NULL,
                                    varname = "proj_yhat",
                                    units = "",
                                    longname = "CEPHALOPOD projected target",
                                    compression = 4) {
  
  # --- 1. Basic checks
  
  if (is.null(dim(array_in))) {
    stop("'array_in' must be an array with dimensions.")
  }
  
  if (length(lon) < 1 || length(lat) < 1) {
    stop("'lon' and 'lat' must be non-empty.")
  }
  
  if (length(month) < 1) {
    stop("'month' must be non-empty.")
  }
  
  if (length(bootstrap) < 1) {
    stop("'bootstrap' must be non-empty.")
  }
  
  
  # --- 2. Determine whether this is a 2-D or 3-D projection
  
  d <- dim(array_in)
  
  if (length(d) == 3) {
    
    # -------------------------------------------------------------------------
    # CEPHALOPOD 2-D format:
    #   cell x bootstrap x month
    #
    # IMPORTANT:
    #
    # The first dimension is NOT an arbitrary spatial ordering.
    # It is the cell ordering inherited from terra::values(r0).
    #
    # For a raster with nlon columns and nlat rows, terra numbers cells:
    #
    #   cell 1 ... nlon       = northernmost row, west -> east
    #   cell nlon+1 ...       = next row south, west -> east
    #   ...
    #
    # Therefore the input vector must first be reconstructed as:
    #
    #   [row north -> south, longitude west -> east]
    #
    # and then flipped vertically and transposed to obtain:
    #
    #   [longitude west -> east, latitude south -> north]
    #
    # Output:
    #   lon x lat x layer x month x bootstrap
    # -------------------------------------------------------------------------
    
    if (!is.null(depth)) {
      stop(
        "A 3-dimensional array with 3 dimensions is interpreted as a ",
        "2-D projection [cell, bootstrap, month]. Do not provide 'depth'."
      )
    }
    
    if (is.null(layer_depth_range) ||
        length(layer_depth_range) != 2) {
      stop(
        "For a 2-D projection, 'layer_depth_range' must contain ",
        "c(zmin, zmax)."
      )
    }
    
    ncell <- d[1]
    nboot <- d[2]
    nmonth <- d[3]
    
    nlon <- length(lon)
    nlat <- length(lat)
    
    if (ncell != nlon * nlat) {
      stop(
        "Spatial dimension mismatch: array contains ",
        ncell,
        " cells, but lon x lat contains ",
        nlon * nlat,
        " cells."
      )
    }
    
    if (nboot != length(bootstrap)) {
      stop(
        "Bootstrap dimension mismatch: array contains ",
        nboot,
        " bootstrap replicates, but ",
        length(bootstrap),
        " bootstrap coordinates were supplied."
      )
    }
    
    if (nmonth != length(month)) {
      stop(
        "Month dimension mismatch: array contains ",
        nmonth,
        " months, but ",
        length(month),
        " month coordinates were supplied."
      )
    }
    
    # Check coordinate orientation.
    #
    # The conversion below assumes that the supplied coordinates are:
    #
    #   lon: west -> east
    #   lat: south -> north
    
    if (nlon > 1 && any(diff(lon) <= 0)) {
      stop("'lon' must be strictly increasing (west -> east).")
    }
    
    if (nlat > 1 && any(diff(lat) <= 0)) {
      stop(
        "'lat' must be strictly increasing (south -> north)."
      )
    }
    
    # Output dimensions:
    #
    #   lon x lat x layer x month x bootstrap
    #
    # This is the LOGICAL order of the array we construct.
    
    out <- array(
      NA_real_,
      dim = c(
        nlon,
        nlat,
        1,
        nmonth,
        nboot
      )
    )
    
    # The CEPHALOPOD raster cell ordering is inherited from:
    #
    #   terra::values(r0)
    #
    # and therefore follows terra's raster cell numbering.
    #
    # Reconstruct each cell vector as:
    #
    #   [latitude north -> south, longitude west -> east]
    
    for (b in seq_len(nboot)) {
      
      for (m in seq_len(nmonth)) {
        
        cell_values <- array_in[, b, m]
        
        # terra::values() returns a ROW-MAJOR vector: the first nlon
        # entries are raster row 1 (north), the next nlon entries are
        # row 2, and so on. To rebuild a [lat, lon] matrix whose rows
        # correspond to raster rows, we must fill it row-by-row
        # (byrow = TRUE). Using byrow = FALSE here would interleave
        # values from different raster rows into the same column
        # whenever nlon != nlat, silently scrambling the grid.
        mat <- matrix(
          cell_values,
          nrow = nlat,
          ncol = nlon,
          byrow = TRUE
        )
        
        # terra raster ordering is north -> south.
        #
        # The supplied latitude coordinate is south -> north,
        # therefore reverse the row order.
        
        mat <- mat[nrow(mat):1, , drop = FALSE]
        
        # mat is now:
        #
        #   [lat south -> north, lon west -> east]
        #
        # Transpose to obtain:
        #
        #   [lon west -> east, lat south -> north]
        
        out[, , 1, m, b] <- t(mat)
      }
    }
    
    # Use the midpoint as the numeric coordinate of the singleton layer.
    # The exact depth interval is stored as layer_bounds below.
    
    layer <- mean(layer_depth_range)
    
    
    # --- 3. Define NetCDF dimensions
    
    dim_lon <- ncdf4::ncdim_def(
      name = "lon",
      units = "degrees_east",
      vals = lon,
      longname = "Longitude"
    )
    
    dim_lat <- ncdf4::ncdim_def(
      name = "lat",
      units = "degrees_north",
      vals = lat,
      longname = "Latitude"
    )
    
    dim_depth <- ncdf4::ncdim_def(
      name = "depth",
      units = "m",
      vals = layer,
      longname = "Mid-depth of the layer"
    )
    
    dim_month <- ncdf4::ncdim_def(
      name = "month",
      units = "1",
      vals = month,
      longname = "Calendar month"
    )
    
    dim_bootstrap <- ncdf4::ncdim_def(
      name = "bootstrap",
      units = "1",
      vals = bootstrap,
      longname = "Bootstrap replicate"
    )
    
    
    # --- 4. Define projection variable
    
    # IMPORTANT:
    #
    # ncdf4 has a different internal storage convention from R arrays.
    #
    # We therefore define the NetCDF variable using the dimension order
    # corresponding to the desired NetCDF representation.
    #
    # The actual spatial reconstruction above is independent of this
    # R/NetCDF storage convention.
    
    var <- ncdf4::ncvar_def(
      name = varname,
      units = units,
      dim = list(
        dim_lon,
        dim_lat,
        dim_depth,
        dim_month,
        dim_bootstrap
      ),
      missval = -9999,
      longname = longname,
      prec = "double",
      compression = compression
    )
    
    
    # --- 5. Create or open file and write projection
    
    nc <- ncdf4::nc_create(
      filename = fn_out,
      vars = var
    )
    
    values_to_write <- out
    values_to_write[is.na(values_to_write)] <- -9999
    
    ncdf4::ncvar_put(
      nc,
      var,
      values_to_write
    )
    
    
    # --- 6. Metadata
    
    ncdf4::ncatt_put(
      nc,
      0,
      "Conventions",
      "CF-1.10"
    )
    
    ncdf4::ncatt_put(
      nc,
      0,
      "projection_dimensions",
      "lon, lat, depth, month, bootstrap"
    )
    
    ncdf4::ncatt_put(
      nc,
      "depth",
      "depth_range",
      paste0(
        layer_depth_range[1],
        "-",
        layer_depth_range[2],
        " m"
      )
    )
    
    ncdf4::ncatt_put(
      nc,
      "depth",
      "bounds",
      paste(
        layer_depth_range,
        collapse = ", "
      )
    )
    
    ncdf4::nc_close(nc)
    
    return(invisible(NULL))
  }
  
  
  # --- 7. Future 3-D projection
  
  if (length(d) == 4) {
    
    # -------------------------------------------------------------------------
    # CEPHALOPOD 3-D format:
    #
    #   cell x depth x bootstrap x month
    #
    # The first dimension is again assumed to follow terra cell ordering:
    #
    #   north -> south, west -> east.
    #
    # Output:
    #
    #   lon x lat x depth x month x bootstrap
    # -------------------------------------------------------------------------
    
    if (is.null(depth)) {
      stop(
        "A 3-D projection requires a 'depth' coordinate."
      )
    }
    
    ncell <- d[1]
    ndepth <- d[2]
    nboot <- d[3]
    nmonth <- d[4]
    
    nlon <- length(lon)
    nlat <- length(lat)
    
    if (ncell != nlon * nlat) {
      stop(
        "Spatial dimension mismatch: array contains ",
        ncell,
        " cells, but lon x lat contains ",
        nlon * nlat,
        " cells."
      )
    }
    
    if (ndepth != length(depth)) {
      stop(
        "Depth dimension mismatch: array contains ",
        ndepth,
        " depth levels, but ",
        length(depth),
        " depth coordinates were supplied."
      )
    }
    
    if (nboot != length(bootstrap)) {
      stop(
        "Bootstrap dimension mismatch."
      )
    }
    
    if (nmonth != length(month)) {
      stop(
        "Month dimension mismatch."
      )
    }
    
    if (nlon > 1 && any(diff(lon) <= 0)) {
      stop("'lon' must be strictly increasing (west -> east).")
    }
    
    if (nlat > 1 && any(diff(lat) <= 0)) {
      stop(
        "'lat' must be strictly increasing (south -> north)."
      )
    }
    
    # Output:
    #
    #   lon x lat x depth x month x bootstrap
    
    out <- array(
      NA_real_,
      dim = c(
        nlon,
        nlat,
        ndepth,
        nmonth,
        nboot
      )
    )
    
    for (b in seq_len(nboot)) {
      
      for (m in seq_len(nmonth)) {
        
        for (z in seq_len(ndepth)) {
          
          cell_values <- array_in[, z, b, m]
          
          # Reconstruct terra cell ordering:
          #
          #   [north -> south, west -> east]
          #
          # terra::values() is row-major, so the matrix must be filled
          # row-by-row (byrow = TRUE) to correctly align each raster
          # row with a matrix row. See explanation in the 2-D branch
          # above.
          
          mat <- matrix(
            cell_values,
            nrow = nlat,
            ncol = nlon,
            byrow = TRUE
          )
          
          # Convert north -> south to south -> north.
          
          mat <- mat[
            nrow(mat):1,
            ,
            drop = FALSE
          ]
          
          # Convert [lat, lon] -> [lon, lat].
          
          out[, , z, m, b] <- t(mat)
        }
      }
    }
    
    
    # --- 8. Define NetCDF dimensions
    
    dim_lon <- ncdf4::ncdim_def(
      name = "lon",
      units = "degrees_east",
      vals = lon,
      longname = "Longitude"
    )
    
    dim_lat <- ncdf4::ncdim_def(
      name = "lat",
      units = "degrees_north",
      vals = lat,
      longname = "Latitude"
    )
    
    dim_depth <- ncdf4::ncdim_def(
      name = "depth",
      units = "m",
      vals = depth,
      longname = "Depth"
    )
    
    dim_month <- ncdf4::ncdim_def(
      name = "month",
      units = "1",
      vals = month,
      longname = "Calendar month"
    )
    
    dim_bootstrap <- ncdf4::ncdim_def(
      name = "bootstrap",
      units = "1",
      vals = bootstrap,
      longname = "Bootstrap replicate"
    )
    
    
    # --- 9. Define variable
    
    var <- ncdf4::ncvar_def(
      name = varname,
      units = units,
      dim = list(
        dim_lon,
        dim_lat,
        dim_depth,
        dim_month,
        dim_bootstrap
      ),
      missval = -9999,
      longname = longname,
      prec = "double",
      compression = compression
    )
    
    
    # --- 10. Write
    
    nc <- ncdf4::nc_create(
      filename = fn_out,
      vars = var
    )
    
    values_to_write <- out
    values_to_write[is.na(values_to_write)] <- -9999
    
    ncdf4::ncvar_put(
      nc,
      var,
      values_to_write
    )
    
    ncdf4::ncatt_put(
      nc,
      0,
      "Conventions",
      "CF-1.10"
    )
    
    ncdf4::ncatt_put(
      nc,
      0,
      "projection_dimensions",
      "lon, lat, depth, month, bootstrap"
    )
    
    ncdf4::nc_close(nc)
    
    return(invisible(NULL))
  }
  
  
  stop(
    "Unsupported projection array. Expected either ",
    "[cell, bootstrap, month] for 2-D or ",
    "[cell, depth, bootstrap, month] for 3-D."
  )
}