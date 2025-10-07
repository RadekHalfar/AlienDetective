fetch_gbif_data <- function(species,
                          hasCoordinate = TRUE,
                          continent = "europe",
                          basisOfRecord = c("OBSERVATION", "MACHINE_OBSERVATION", "HUMAN_OBSERVATION", 
                                          "MATERIAL_SAMPLE", "LIVING_SPECIMEN", "OCCURRENCE"),
                          fields = c("decimalLatitude", "decimalLongitude", "year", "month", "country"),
                          limit = 10000,
                          output_dir) {
  
  # Input validation
  if (!is.character(species) || length(species) != 1) {
    stop("species must be a single character string")
  }
  
  # Get GBIF data
  data_list <- rgbif::occ_search(
    scientificName = species,
    hasCoordinate = hasCoordinate,
    continent = continent,
    basisOfRecord = basisOfRecord,
    fields = fields,
    limit = limit
  )

  # Collapse all returned data.frames in one step – rbindlist handles conversion
  res <- rbindlist(
    lapply(data_list, `[[`, "data"),   # pull the $data element
    idcol = "basisOfRecord",           # names(data_list) becomes the id values
    use.names = TRUE,
    fill = TRUE
  )
  if (nrow(res) == 0) {
    message(sprintf("No GBIF records found for species '%s'", species))
    return(NULL)
  }
  
  # Ensure all required fields exist
  missing_cols <- setdiff(fields, names(res))
  if (length(missing_cols) > 0) {
    res[, (missing_cols) := NA]
  }
  
  # Rename latitude/longitude columns once
  rename_map <- c(decimalLatitude = "latitude", decimalLongitude = "longitude")
  setnames(res, old = names(rename_map), new = unname(rename_map), skip_absent = TRUE)
  
  # Prepare output column order (after renaming)
  out_cols <- c("latitude", "longitude", setdiff(fields, names(rename_map)), "basisOfRecord")
  out_cols <- unique(out_cols[out_cols %in% names(res)])
  
  # Filter rows with valid coords and keep relevant columns
  res <- res[!is.na(latitude) & !is.na(longitude), ..out_cols]
  
  # Key by coordinates for faster joins later
  setkey(res, latitude, longitude)
  
  return(res[])
}


# Function to check if point is on land (TRUE = land, FALSE = sea)
is_on_land <- function(lat, lon, r) {
  point <- sp::SpatialPoints(cbind(lon, lat), proj4string = sp::CRS(proj4string(r)))
  return(is.na(raster::extract(r, point)))
}

# -----------------------------------------------------------------------------
# Fast nearest-sea search ------------------------------------------------------
# -----------------------------------------------------------------------------
# Move a land point to the nearest sea raster cell.  The expensive extraction of
# *all* sea-cell coordinates and the KD-tree index are done only once and cached
# in the global environment, so repeated calls are cheap.
# Returns list(coords = c(lon, lat), dist = distance_m).
move_to_sea <- function(lat, lon, r, cost_matrix, sea_cell_coords = NULL) {
  # -------------------------------------------------------------------
  # Build sea-cell coordinate matrix once, or reuse what the caller sent
  # -------------------------------------------------------------------
  if (is.null(sea_cell_coords)) {
    trans_matrix <- gdistance::transitionMatrix(cost_matrix)
    connected_cells <- which(rowSums(trans_matrix != 0) > 0)
    coords <- raster::xyFromCell(r, connected_cells)
    sea_idx <- which(raster::extract(r, coords) == 1)
    if (!length(sea_idx)) return(NULL)
    sea_cell_coords <- coords[sea_idx, , drop = FALSE]
  }
  
  # Filter to only retain sea cells (safety check in case caller supplied
  # a matrix that is out-of-date for this raster)
  is_sea <- raster::extract(r, sea_cell_coords) == 1
  sea_coords <- sea_cell_coords[is_sea, , drop = FALSE]
  
  if (nrow(sea_coords) == 0) {
    return(NULL)  # failure signal
  } 
  
  point_coords <- c(lon, lat)
  
  for (radius_km in seq(5, 100, 5)) {
    radius_deg <- radius_km / 111 # rough conversion km to degrees
    
    # Borders filter box
    lon_min <- lon - radius_deg
    lon_max <- lon + radius_deg
    lat_min <- lat - radius_deg
    lat_max <- lat + radius_deg
    
    # Define sea cells within radius
    sea_in_radius <- which (sea_coords[,1] >= lon_min & sea_coords[,1] <= lon_max &
                            sea_coords[,2] >= lat_min & sea_coords[,2] <= lat_max)
    
    if (length(sea_in_radius) > 0) {
      sea_coords_radius <- sea_coords[sea_in_radius, , drop = FALSE]
      dists <- geosphere::distVincentySphere(point_coords, sea_coords_radius)
      
      # return nearest seapoint
      nearest_idx <- which.min(dists)
      dist <- dists[nearest_idx]
      new_coords <- sea_coords_radius[nearest_idx, , drop = FALSE]
      return(list(
        coords = as.vector(new_coords),
        dist = dist,
        sea_cell_coords = sea_cell_coords
      ))
    }
    # Otherwise, continues with next larger radius
  }
  return(NULL) # when no sea point found
}

# Main function: calculates both sea route and geodesic distances from every downloaded GBIF occurrence to the species occurrence in question
calculate.distances <- function(data, raster_map, cost_matrix, chunk_size = 5000) {
  # Validate input -----------------------------------------------------------
  if (is.null(data) || nrow(data) == 0) {
    return(list(sea_distances = NULL,
                geodesic_distances = NULL,
                error_messages = "Input table is NULL or empty"))
  }
  required_cols <- c("Latitude_missing_locs", "Longitude_missing_locs",
                     "latitude", "longitude",
                     "latitude_moved", "longitude_moved")
  missing <- setdiff(required_cols, names(data))
  if (length(missing)) {
    stop(sprintf("data is missing required columns: %s", paste(missing, collapse = ", ")))
  }

  # Build origin (missing_locs) & destination (occurrence) coordinates
  origin_lon <- data$Longitude_missing_locs
  origin_lat <- data$Latitude_missing_locs
  dest_lon <- ifelse(is.na(data$longitude_moved), data$longitude, data$longitude_moved)
  dest_lat <- ifelse(is.na(data$latitude_moved),  data$latitude,  data$latitude_moved)

  n <- length(origin_lon)
  sea_dist <- rep(NA_real_, n)
  # Only compute sea-route where destination is at sea
  at_sea <- raster::extract(raster_map, cbind(dest_lon, dest_lat)) == 1L
  at_sea[is.na(at_sea)] <- FALSE # avoid propagating NAs

  if (any(at_sea)) {
    crs_wgs84 <- sp::CRS("+proj=longlat +datum=WGS84")
    idx <- which(at_sea)
    # Chunk indices
    chunk_idxs <- split(idx, ceiling(seq_along(idx) / chunk_size))
    for (chunk in chunk_idxs) {
      sea_from <- sp::SpatialPoints(cbind(origin_lon[chunk], origin_lat[chunk]), proj4string = crs_wgs84)
      sea_to   <- sp::SpatialPoints(cbind(dest_lon[chunk],   dest_lat[chunk]),   proj4string = crs_wgs84)
      # Try-catch to avoid errors for unreachable points
#      res <- tryCatch(
#        as.numeric(diag(gdistance::costDistance(cost_matrix, sea_from, sea_to))) / 1000,
#        error = function(e) rep(NA_real_, length(chunk))
#      )
#      sea_dist[chunk] <- res

      sea_dist[chunk] <- as.numeric(diag(gdistance::costDistance(cost_matrix, sea_from, sea_to))) / 1000

    }
  }

  # Row-wise geodesic distances (vectorized, not memory intensive)
  geo_dist <- as.numeric(
    geodist::geodist(
      data.table::data.table(lon = origin_lon, lat = origin_lat),
      data.table::data.table(lon = dest_lon,   lat = dest_lat),
      paired = TRUE,
      measure = "geodesic"
    )
  ) / 1000

  list(sea_distances       = round(sea_dist, 0),
       geodesic_distances = round(geo_dist, 0),
       error_messages     = NULL)
}

# Process a data.table of coordinates and, if necessary, move points on land to the nearest sea cell.
# Expects a data.table with at least the columns: "species", "latitude", "longitude".
# Returns the same table plus: latitude_moved, longitude_moved, dist_moved (km).
process_coords <- function(coords_dt, r, cost_matrix, chunk_size = 5000) {
  # -------------------------
  # Input validation
  # -------------------------
  
  if (!data.table::is.data.table(coords_dt)) {
    stop("coords_dt must be a data.table")
  }
  required_cols <- c("species", "latitude", "longitude")
  missing <- setdiff(required_cols, names(coords_dt))
  if (length(missing)) {
    stop(sprintf("coords_dt is missing required columns: %s", paste(missing, collapse = ", ")))
  }
  
  # Make a copy so we don't modify the original
  result <- data.table::copy(coords_dt)[,
    `:=`(latitude_moved = as.numeric(NA),
         longitude_moved = as.numeric(NA),
         dist_moved = 0.0)]

  # Work on unique coordinate rows ACROSS SPECIES to avoid redundant processing
  uniq <- unique(result[, .(latitude, longitude)])
  uniq[, `:=`(latitude_moved = as.numeric(NA),
              longitude_moved = as.numeric(NA),
              dist_moved = 0.0)]

  # Process coordinates in chunks to avoid memory issues (helps with very large tables)
  n_chunks <- ceiling(nrow(uniq) / chunk_size)

  # Reusable sea-cell coordinate cache local to this call
  sea_cell_coords <- NULL
  # Memoization cache for move_to_sea() results keyed by raster cell id
  cache <- new.env(parent = emptyenv())
  get_move <- function(cell, lat, lon) {
    key <- as.character(cell)
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache))
    }
    res <- move_to_sea(lat, lon, r, cost_matrix, sea_cell_coords)
    assign(key, res, envir = cache)
    res
  }

  for (i in seq_len(n_chunks)) {
    idx_start <- (i - 1) * chunk_size + 1
    idx_end <- min(i * chunk_size, nrow(uniq))
    chunk <- uniq[idx_start:idx_end]

    # Determine whether each point is on land (vectorized raster extract)
    # r has sea == 1 and land == NA; land if extracted value is NA
    xy <- cbind(chunk$longitude, chunk$latitude)
    vals <- raster::extract(r, xy)
    cells <- raster::cellFromXY(r, xy)
    chunk[, is_land := is.na(vals)]

    # Handle points on land
    land_idx <- which(chunk$is_land)
    moved <- vector("list", length(land_idx))
    if (length(land_idx) > 0) {
      for (k in seq_along(land_idx)) {
        idx_pt <- land_idx[k]
        res <- get_move(cells[idx_pt], chunk$latitude[idx_pt], chunk$longitude[idx_pt])
        # update local cache for subsequent iterations (simple assignment)
        if (!is.null(res$sea_cell_coords)) {
          sea_cell_coords <- res$sea_cell_coords
        }
        moved[[k]] <- res
      }
    }

    # Vectorized assignments: build output vectors once and assign
    n <- nrow(chunk)
    lat_mov  <- chunk$latitude
    lon_mov  <- chunk$longitude
    dist_mov <- numeric(n)

    # Only update those that were on land and got moved (if Cannot Be Moved to Sea return NA)
    if (length(land_idx) > 0) {
      lat_mov[land_idx]  <- vapply(moved, function(x) if (is.null(x) || is.null(x$coords)) NA_real_ else x$coords[2], numeric(1))
      lon_mov[land_idx]  <- vapply(moved, function(x) if (is.null(x) || is.null(x$coords)) NA_real_ else x$coords[1], numeric(1))
      dist_mov[land_idx] <- vapply(moved, function(x) if (is.null(x) || is.null(x$coords)) NA_real_ else if (is.null(x$dist)) 0 else x$dist, numeric(1))
    }
    chunk[, `:=`(
      latitude_moved  = lat_mov,
      longitude_moved = lon_mov,
      dist_moved      = dist_mov
    )]

    # Remove temporary column
    chunk[, is_land := NULL]

    # Update unique table
    uniq[idx_start:idx_end] <- chunk
  }

  # Join moved results back to the full table (including duplicates across species)
  data.table::setkeyv(uniq, c("latitude", "longitude"))
  data.table::setkeyv(result, c("latitude", "longitude"))
  result[uniq, on = .(latitude, longitude), `:=`(
    latitude_moved  = i.latitude_moved,
    longitude_moved = i.longitude_moved,
    dist_moved      = i.dist_moved
  )]
  
  # -------------------------
  # Final housekeeping & return
  # -------------------------
  setcolorder(result, c("species", "latitude", "longitude",
                        "latitude_moved", "longitude_moved", "dist_moved"))
  data.table::setkeyv(result, c("species", "latitude", "longitude"))
  return(result[])
}

# compute missing distance locations
#' locations_to_process: optional character vector of location column names to compute. If NULL, process all detected locations.
process_species_locations <- function(species, species_location, location_coordinates,
                                      unique_coords, r, cost_matrix,
                                      locations_to_process = NULL) {
  # Input validation
  stopifnot(
    is.character(species) && length(species) == 1,
    data.table::is.data.table(species_location) || is.data.frame(species_location),
    data.table::is.data.table(location_coordinates) || is.data.frame(location_coordinates),
    data.table::is.data.table(unique_coords) || is.data.frame(unique_coords)
  )
  
  # Convert inputs to data.table if needed
  if (!data.table::is.data.table(species_location)) {
    species_location <- as.data.table(species_location)
  }
  if (!data.table::is.data.table(location_coordinates)) {
    location_coordinates <- as.data.table(location_coordinates)
  }
  if (!data.table::is.data.table(unique_coords)) {
    unique_coords <- as.data.table(unique_coords)
  }
  
  # Make a copy to avoid modifying the original
  result_dt <- data.table::copy(unique_coords)
  
  # Set keys for faster joins
  species_col <- names(species_location)[1]
  data.table::setkeyv(species_location, species_col)
  data.table::setkeyv(location_coordinates, "Observatory.ID")
  # Get locations where species was detected (using data.table's fast subset)
  loc_presence <- species_location[.(species), .SD, .SDcols = -1]
  
  #  if (nrow(loc_presence) == 0) {
  #    message(sprintf("Species '%s' not found in species_location table", species))
  #    return(result_dt)
  #  }
  
  presence_vals <- as.numeric(loc_presence[1])
  detected_locations <- names(loc_presence)[!is.na(presence_vals) & presence_vals >= 1]
  
  if (length(detected_locations) == 0) {
    message(sprintf("No detections found for species: %s", species))
    return(result_dt)
  }
  
  # If caller supplied subset of locations, filter
  if (!is.null(locations_to_process)) {
    detected_locations <- intersect(detected_locations, locations_to_process)
  }
  
  # Process locations in parallel if possible
  results <- lapply(detected_locations, function(loc) {
    # Get coordinates for this location
    coords <- location_coordinates[.(loc), nomatch = NULL]
    if (nrow(coords) != 1) {
      message(sprintf("Could not retrieve coordinates for location: %s", loc))
      return(NULL)
    }
    
    # Convert coordinates to numeric (handling comma as decimal separator)
    coords[, `:=`(
      lat = as.numeric(gsub(",", ".", Latitude)),
      lon = as.numeric(gsub(",", ".", Longitude))
    )]
    # Calculate distances
    message(sprintf("Calculating distances to %s occurrences from %s", species, loc))
    
    tryCatch({
      result <- calculate.distances(
        data = result_dt,
        latitude = coords$lat[1],
        longitude = coords$lon[1],
        raster_map = r,
        cost_matrix = cost_matrix
      )
      # Return results with location prefix
      if (!is.null(result$sea_distances) && !is.null(result$geodesic_distances)) {
        data.table(
          id = seq_len(nrow(result_dt)),
          sea = result$sea_distances,
          geo = result$geodesic_distances
        )
      } else {
        NULL
      }
    }, error = function(e) {
      message(sprintf("Error processing location %s: %s", loc, e$message))
      NULL
    })
  })
  
  # Combine results
  valid_results <- Filter(Negate(is.null), results)
  if (length(valid_results) > 0) {
    # Combine all results into one long table with an index per location
    combined <- rbindlist(valid_results, idcol = "idx")
    combined[, location := detected_locations[idx]]
    
    # Reshape to wide format in one step (creates columns like "sea_<loc>", "geo_<loc>")
    wide <- data.table::dcast(
      combined,
      id ~ location,
      value.var = c("sea", "geo"),
      fill = NA_real_
    )
    
    # Rename columns to "<loc>_seaway" / "<loc>_geodesic"
    new_names <- names(wide)
    new_names <- gsub("^sea_(.*)$", "\\1_seaway", new_names)
    new_names <- gsub("^geo_(.*)$", "\\1_geodesic", new_names)
    data.table::setnames(wide, new_names)
    
    # Append reshaped columns to result_dt (drop the "id" column first)
    result_dt <- cbind(result_dt, wide[, -1, with = FALSE])
    
    # Remove duplicated columns if any (can occur after reruns)
    dup_cols <- duplicated(names(result_dt))
    if (any(dup_cols)) {
      #      warning(sprintf(
      #        "Removing duplicated columns from result_dt: %s",
      #        paste(names(result_dt)[dup_cols], collapse = ", ")
      #      ))
      result_dt <- result_dt[, !dup_cols, with = FALSE]
    }
  }
  
  # Optimize memory usage
  data.table::setalloccol(result_dt)
  return(result_dt[])
}