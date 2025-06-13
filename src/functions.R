#' @title Alien Detective Functions
#' @description This script contains the custom functions used by the other scripts
#' @import data.table
#' @importFrom utils globalVariables
NULL

# This file contains various functions for the Alien Detective project
# including data fetching, processing, and visualization utilities

## Declare global variables to avoid R CMD check notes
#globalVariables(
#  c(
#    # data.table special symbols
#    ".", 
#    ":=", 
#    ".SD", 
#    ".N", 
#    ".I", 
#    ".GRP", 
#    ".BY", 
#    ".EACHI",
#    
#    # Our variables
#    "longitude_moved",
#    "latitude_moved",
#    "x_coord",
#    "y_coord",
#    "Observatory.ID",
#    "Latitude",
#    "Longitude",
#    "location",
#    "count",
#    "sea_distances",
#    "geodesic_distances"
#  )
#)


#' Fetch GBIF occurrence data for a species
#' 
#' @param species Character. The species name to search for.
#' @param hasCoordinate Logical. Return only records with coordinates? Default TRUE.
#' @param continent Character. Continent to search in. Default "europe".
#' @param basisOfRecord Character vector. Types of records to include.
#' @param fields Character vector. Fields to retrieve from GBIF.
#' @param limit Numeric. Maximum number of records to return.
#' @param output_dir Character. Directory to save results (unused).
#' @return A data.table with occurrence data or NULL if no data found.
#' @import data.table
#' @export
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
is_on_land <- function(lat, lon) {
  point <- sp::SpatialPoints(cbind(lon, lat), proj4string = sp::CRS(proj4string(r)))
  return(is.na(raster::extract(r, point)))
}

# Function to move point on land to sea
move_to_sea <- function(lat, lon) {
  # Get transition matrix & all connected cells
  trans_matrix <- gdistance::transitionMatrix(cost_matrix)
  connected_cells <- which(rowSums(trans_matrix != 0) > 0)
  connected_coords <- raster::xyFromCell(r, connected_cells)
  
  # Filter to only retain sea cells
  is_sea <- raster::extract(r, connected_coords) == 1
  sea_coords <- connected_coords[is_sea, , drop = FALSE]
  
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
        dist = dist
      ))
    }
    # Otherwise, continues with next larger radius
  }
  return(NULL) # when no sea point found
}


# Main function: calculates both sea route and geodesic distances from every downloaded GBIF occurrence to the species occurrence in question
calculate.distances <- function(data, latitude, longitude, raster_map, cost_matrix){
  
  if (is.null(data)) return(list(sea_distances = NULL, geodesic_distances = NULL, error_messages = "Input table is NULL"))
  if (nrow(data) < 1) return(list(sea_distances = NULL, geodesic_distances = NULL, error_messages = "Input table is has no entries"))
  
  tryCatch({
    # Specify the PROJ4 string for WGS84
    proj4_crs <- sp::CRS("+init=EPSG:4326")
    
    # Create SpatialPoints objects from the coordinates
    query_point <- sp::SpatialPoints(cbind(longitude, latitude), proj4string = proj4_crs)
    ref_points <- sp::SpatialPoints(cbind(ifelse(is.na(data$longitude_moved), data$longitude, data$longitude_moved),
                                          ifelse(is.na(data$latitude_moved), data$latitude, data$latitude_moved)),
                                    proj4string = proj4_crs)
    
    
    # Get raster cell values of the GBIF occurrence points (1 for sea, Inf for land)
    cell_values <- raster::extract(raster_map, ref_points)
    # Initialize result vectors
    sea_distances <- rep(NA_real_, length(cell_values))
    geodesic_distances <- rep(NA_real_, length(cell_values))
    # Get indexes of the points that are in the sea
    indexes <- which(cell_values == 1L)
    if (length(indexes) > 0) {
      # Subset points that are in the sea
      ref_points_sea <- ref_points[indexes,]
      # Vectorized sea distance calculation to all GBIF occurrences in the sea
      sea_distances[indexes] <- as.numeric(gdistance::costDistance(cost_matrix, query_point, ref_points_sea)[1,])
      # Convert points to simple table format for use with geodist
      query_point_table <- data.frame(lon = sp::coordinates(query_point)[,1],
                                      lat = sp::coordinates(query_point)[,2])
      ref_points_sea_table <- data.frame(lon = sp::coordinates(ref_points_sea)[,1],
                                         lat = sp::coordinates(ref_points_sea)[,2])
      # Vecotrized geodesic distance calculation to all GBIF occurrences in the sea
      geodesic_distances[indexes] <- as.numeric(geodist::geodist(query_point_table, ref_points_sea_table, measure = "geodesic"))
      # Convert distances to kilometres
      sea_distances <- round(sea_distances / 1000, 0)
      geodesic_distances <- round(geodesic_distances / 1000, 0)
    }
    # Return result
    return(list(sea_distances = sea_distances, geodesic_distances = geodesic_distances, error_messages = NULL))
  }, error = function(e) {
    error_messages <- paste0("An error occurred during distance calculation for ", species, " in ", location, ": ", e$message)
    return(list(sea_distances = NULL, geodesic_distances = NULL, error_messages = error_messages))
  })
}



##########################
### PLOTTING FUNCTIONS ###
##########################

# Generic helper --------------------------------------------------------------
make_hist_plot <- function(data,
                           x_col,
                           fill_col = NULL,
                           title = "",
                           binwidth = 50,
                           breaks_by = 250,
                           x_label = "Distance (km)",
                           y_label = "Frequency") {
  stopifnot(is.data.frame(data), x_col %in% names(data))

  max_x <- max(data[[x_col]], na.rm = TRUE)

  # Dynamically build aesthetic mapping
  mapping <- ggplot2::aes_string(x = x_col)
  if (!is.null(fill_col)) {
    mapping <- ggplot2::aes_string(x = x_col, fill = fill_col)
  }

  ggplot2::ggplot(data, mapping) +
    ggplot2::geom_histogram(binwidth = binwidth, boundary = 0, alpha = 0.7,
                            position = "stack", colour = "#e9ecef") +
    ggplot2::labs(title = title, x = x_label, y = y_label, fill = fill_col) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::scale_x_continuous(breaks = seq(0, max_x * 1.1, by = breaks_by),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(0, max_x * 1.1))
}

# -----------------------------------------------------------------------------
# Plot wrappers (all rely on make_hist_plot) ----------------------------------
# -----------------------------------------------------------------------------

country.final <- function(species, data, output_dir) {
  plot <- make_hist_plot(
    data      = data,
    x_col     = "x",
    fill_col  = "location",
    title     = sprintf("Sea-route distances for %s (all ARMS locations)", species)
  )
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir, sprintf("%s_from_all_ARMS_locations.png", gsub(" ", "_", species))),
                  plot, width = 2400, height = 1200, units = "px", dpi = 300)
  invisible(plot)
}

plot.dist.sea <- function(species, location, data, output_dir) {
  plot <- make_hist_plot(
    data      = data,
    x_col     = "x",
    title     = sprintf("Sea-route distances for %s from %s", species, location)
  )
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir, sprintf("%s_from_%s.png", gsub(" ", "_", species), location)),
                  plot, width = 2000, height = 1200, units = "px", dpi = 300)
  invisible(plot)
}

plot.dist.both <- function(species, location, data, output_dir) {
  plot <- make_hist_plot(
    data      = data,
    x_col     = "x",
    fill_col  = "DistanceType",
    title     = sprintf("Sea vs geodesic distances for %s (%s)", species, location)
  ) + ggplot2::labs(fill = "DistanceType")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir, sprintf("%s_from_%s_seadist&geodesic.png", gsub(" ", "_", species), location)),
                  plot, width = 2400, height = 1200, units = "px", dpi = 300)
  invisible(plot)
}

plot.dist.by.country <- function(species, location, data, output_dir) {
  plot <- make_hist_plot(
    data      = data,
    x_col     = "x",
    fill_col  = "country",
    title     = sprintf("Sea-route distances by country for %s (%s)", species, location)
  )
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir, sprintf("%s_from_%s_by_country.png", gsub(" ", "_", species), location)),
                  plot, width = 2400, height = 1200, units = "px", dpi = 300)
  invisible(plot)
}

plot.dist.by.year <- function(species, location, data, output_dir) {
  plot <- make_hist_plot(
    data      = data,
    x_col     = "x",
    fill_col  = "year_category",
    title     = sprintf("Sea-route distances by year for %s (%s)", species, location)
  ) +
    ggplot2::scale_fill_brewer(palette = "YlOrRd", na.value = "black") +
    ggplot2::labs(fill = "Year")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir, sprintf("%s_from_%s_by_year.png", gsub(" ", "_", species), location)),
                  plot, width = 2400, height = 1200, units = "px", dpi = 300)
  invisible(plot)
}

# Process each coordinate pair
#' Process coordinates to ensure they are at sea
#' 
#' @param lat Numeric vector of latitudes
#' @param lon Numeric vector of longitudes
#' @return A data.table with processed coordinates and distances moved
#' @import data.table
#' @export
process_coords <- function(lat, lon) {
  # Input validation
  if (length(lat) != length(lon)) {
    stop("lat and lon must be of the same length")
  }
  
  # Initialize result data.table
  result <- data.table(
    latitude = lat,
    longitude = lon,
    latitude_moved = as.numeric(NA),
    longitude_moved = as.numeric(NA),
    dist_moved = 0.0
  )
  
  # Process coordinates in chunks to avoid memory issues
  chunk_size <- 1000
  n_chunks <- ceiling(nrow(result) / chunk_size)
  
  for (i in seq_len(n_chunks)) {
    idx_start <- (i - 1) * chunk_size + 1
    idx_end <- min(i * chunk_size, nrow(result))
    chunk <- result[idx_start:idx_end]
    
    # Process coordinates in parallel if possible
    chunk[, is_land := mapply(is_on_land, latitude, longitude)]
    
    # Handle points on land
    land_idx <- which(chunk$is_land)
    if (length(land_idx) > 0) {
      moved <- lapply(land_idx, function(i) {
        move_to_sea(chunk$latitude[i], chunk$longitude[i])
      })
      
      # Update moved coordinates
      for (j in seq_along(land_idx)) {
        idx <- land_idx[j]
        if (!is.null(moved[[j]])) {
          set(chunk, i = idx, j = "latitude_moved", value = moved[[j]]$coords[2])
          set(chunk, i = idx, j = "longitude_moved", value = moved[[j]]$coords[1])
          set(chunk, i = idx, j = "dist_moved", value = round(moved[[j]]$dist/1000, 2))
        } else {
          # If move_to_sea failed, keep original coords with NA for moved columns
          set(chunk, i = idx, j = "latitude_moved", value = chunk$latitude[idx])
          set(chunk, i = idx, j = "longitude_moved", value = chunk$longitude[idx])
        }
      }
    }
    
    # For points already at sea, just copy the coordinates
    sea_idx <- which(!chunk$is_land)
    if (length(sea_idx) > 0) {
      set(chunk, i = sea_idx, j = "latitude_moved", value = chunk$latitude[sea_idx])
      set(chunk, i = sea_idx, j = "longitude_moved", value = chunk$longitude[sea_idx])
    }
    
    # Remove temporary column
    chunk[, is_land := NULL]
    
    # Update result
    result[idx_start:idx_end] <- chunk
  }
  
  # Set column order and return
  setcolorder(result, c("latitude", "longitude", "latitude_moved", "longitude_moved", "dist_moved"))
  setkeyv(result, c("latitude", "longitude"))
  
  return(result[])
}

# Process species locations and calculate distances
#' Process species locations and calculate distances
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

# In src/functions.R

# plot_leaflet_map <- function(df, species_name, output_dir = NULL) {
#   if (!requireNamespace("leaflet", quietly = TRUE)) {
#     install.packages("leaflet")
#   }
#   library(leaflet)
#   
#   # Basic validation
#   if (!("latitude" %in% tolower(names(df))) || !("longitude" %in% tolower(names(df)))) {
#     stop("Dataframe must contain 'latitude' and 'longitude' columns")
#   }
#   
#   lat_col <- grep("latitude", names(df), ignore.case = TRUE, value = TRUE)
#   lon_col <- grep("longitude", names(df), ignore.case = TRUE, value = TRUE)
#   
#   map <- leaflet(df) |>
#     addTiles() |>
#     addCircleMarkers(
#       lng = ~get(lon_col), lat = ~get(lat_col),
#       popup = ~paste("Year:", df$year),
#       color = "blue", radius = 3, stroke = FALSE, fillOpacity = 0.6
#     ) |>
#     addLegend("bottomright", colors = "blue", labels = "Observations",
#               title = species_name)
#   
#   # Optional: save HTML map
#   if (!is.null(output_dir)) {
#     html_file <- file.path(output_dir, paste0(gsub(" ", "_", species_name), "_map.html"))
#     htmlwidgets::saveWidget(map, file = html_file, selfcontained = TRUE)
#   }
#   
#   return(map)
# }