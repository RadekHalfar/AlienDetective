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
calculate.distances <- function(data, latitude, longitude, raster_map, cost_matrix) {
  # Input validation
  if (is.null(data)) {
    return(list(sea_distances = NULL, 
                geodesic_distances = NULL, 
                error_messages = "Input table is NULL"))
  }
  if (nrow(data) < 1) {
    return(list(sea_distances = NULL, 
                geodesic_distances = NULL, 
                error_messages = "Input table has no entries"))
  }
  
  # Ensure data is a data.table (no copy)
  data.table::setDT(data)
  
  tryCatch({
    # Specify the PROJ4 string for WGS84
    proj4_crs <- sp::CRS("+init=EPSG:4326")
    
    # Work directly on the incoming table
    coords_dt <- data
    
    # Check if the required columns exist
    if (!all(c("longitude", "latitude") %in% names(coords_dt))) {
      return(list(sea_distances = NULL, 
                 geodesic_distances = NULL, 
                 error_messages = "Missing required longitude/latitude columns"))
    }
    
    # Handle missing moved coordinates
    if (!"longitude_moved" %in% names(coords_dt)) {
      coords_dt[, longitude_moved := NA_real_]
    }
    if (!"latitude_moved" %in% names(coords_dt)) {
      coords_dt[, latitude_moved := NA_real_]
    }
    
    # Update coordinates
    coords_dt[, `:=`(
      x_coord = data.table::fifelse(is.na(longitude_moved), longitude, longitude_moved),
      y_coord = data.table::fifelse(is.na(latitude_moved), latitude, latitude_moved)
    )]
    
    # Create SpatialPoints objects
    query_point <- sp::SpatialPoints(
      cbind(longitude, latitude), 
      proj4string = proj4_crs
    )
    
    ref_points <- sp::SpatialPoints(
      coords_dt[, .(x_coord, y_coord)],
      proj4string = proj4_crs
    )
    
    # Get raster cell values (1 for sea, NA for land)
    cell_values <- raster::extract(raster_map, ref_points)
    
    # Initialize result vectors with NAs
    n_points <- nrow(coords_dt)
    sea_distances <- rep(NA_real_, n_points)
    geodesic_distances <- rep(NA_real_, n_points)
    
    # Process points that are in the sea
    sea_indexes <- which(cell_values == 1L)
    
    if (length(sea_indexes) > 0) {
      ref_points_sea <- ref_points[sea_indexes, ]
      
      # Calculate sea distances (vectorized)
      sea_distances[sea_indexes] <- as.numeric(
        gdistance::costDistance(cost_matrix, query_point, ref_points_sea)[1,]
      )
      
      # Calculate geodesic distances (vectorized)
      query_coords <- sp::coordinates(query_point)
      sea_coords <- sp::coordinates(ref_points_sea)
      
      geodesic_distances[sea_indexes] <- as.numeric(
        geodist::geodist(
          x1 = query_coords[1, 1], y1 = query_coords[1, 2],
          x2 = sea_coords[, 1], y2 = sea_coords[, 2],
          measure = "geodesic"
        )
      )
      
      # Convert distances to kilometers and round
      sea_distances <- round(sea_distances / 1000, 0)
      geodesic_distances <- round(geodesic_distances / 1000, 0)
    }
    
    # Return results
    return(list(
      sea_distances = sea_distances,
      geodesic_distances = geodesic_distances,
      error_messages = NULL
    ))
    
  }, error = function(e) {
    # Error handling
    error_msg <- paste0(
      "Error in calculate.distances: ",
      conditionMessage(e), "\n",
      "Call: ", deparse(conditionCall(e))
    )
    
    return(list(
      sea_distances = NULL,
      geodesic_distances = NULL,
      error_messages = error_msg
    ))
  })
}



##########################
### PLOTTING FUNCTIONS ###
##########################

country.final <- function(species, distances, output_dir) {
  #hist_info <- hist(long_sea$x, plot = FALSE)
  #max_count <- max(hist_info$counts)
  max_x <- max(long_sea$x)
  plot <- ggplot(long_sea, aes(x = x, fill = location)) +
    geom_histogram(binwidth = 50, boundary = 0, position = "stack") + # adjust the binwidth to personal preference
    labs(title = paste0("Frequencies of Sea distances for ", species," from all ARMS locations" ),
         x = "Sea distance in km", y = "Frequency of species") +
    theme_bw() +
    #scale_fill_brewer(palette = "Set1") +  # You can choose a different palette if you like
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 10),           # Set font size for axis numbers
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 14),   # Increase legend title size
          legend.text = element_text(size = 12),    # Increase legend text size
          legend.key.size = unit(1.5, "lines")) +   # Increase legend key size
    scale_x_continuous(breaks = seq(0, max_x*1.1, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, max_x*1.1)) # Use coord_cartesian for setting limits
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(output_dir, paste0(gsub(" ", "_", species), "_from_all_ARMS_locations", ".png")), 
         plot = plot, width = 2400, height = 1200, units = "px", dpi = 300)
  return(plot)
}

# Make histogram of sea distances
plot.dist.sea <- function(species, location, distances, output_dir) {
  hist_info <- hist(sea_loc_data$x, plot = FALSE)
  max_count <- max(hist_info$counts)
  max_x <- max(sea_loc_data$x)
  # make histograms of distances per species, with filtering on distance limit 40000
  plot <- ggplot(sea_loc_data, aes(x = x, fill = location)) +
    geom_histogram(binwidth = 50, boundary = 0, position = "stack", alpha = 0.7) +  # default is position = "stack"
    labs(title = paste("Distances for", species), x = "Distance (km)", y = "Count") +
    theme_minimal()+
    #scale_fill_brewer(palette = "Set1") +
    ggtitle(paste0("Distribution of ", species, " from", loc)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement with hjust
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 10),           # Set font size for axis numbers
          axis.title = element_text(size = 20)) +         # Set font size for axis titles
    scale_x_continuous(breaks = seq(0, max_x*1.1, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, max_x*1.1), ylim = c(0,max_count*1.1))
  
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(output_dir, paste0(gsub(" ", "_", species), "_from_", location, ".png")), 
         plot = plot, width = 2000, height = 1200, units = "px", dpi = 300)
  return(plot)
}

# Make combined histogram of sea distances and fly distances
plot.dist.both <- function(species, location, distances, output_dir) {
  combined_distances <- rbind(sea_loc_data, geo_loc_data)
  hist_info <- hist(combined_distances$x, plot = FALSE)
  max_count <- max(hist_info$counts)
  max_x <- max(sea_loc_data$x)
  plot <- ggplot(combined_distances, aes(x = x, fill = DistanceType)) +
    geom_histogram(binwidth = 50, color="#e9ecef", alpha=0.6, position = 'identity') +
    theme_bw() +
    #scale_fill_brewer(palette = "Set1") +
    labs(x = "Distance in km", y = "Frequency") +
    ggtitle(paste0("Distribution of ", species, " from sea and fly distances")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
      plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
      axis.text = element_text(size = 10),  # Set font size for axis numbers
      axis.title = element_text(size = 16), # Set font size for title
      legend.title = element_text(size = 18, face="bold"), # Settings for legend title
      legend.text = element_text(size = 16)) +  # settings for legend text
    scale_x_continuous(breaks = seq(0, max_x*1.1, by = 250), expand = c(0, 0)) +  # settings for x axis
    scale_y_continuous(expand = c(0, 0)) +
    # used expand to make sure the axes are on the lines of the axes and not above them floating
    coord_cartesian(xlim = c(0, max_x*1.1), ylim = c(0,max_count*1.1)) + # Use coord_cartesian for setting limits
    # set legend title and labels
    labs(x = "Distance in km", y = "Frequency", fill = "DistanceType")
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(output_dir, paste0(gsub(" ", "_", species), "_from_", location, "_seadist&geodesic.png")), 
         plot = plot, width = 2400, height = 1200, units = "px", dpi = 300)
  return(plot)
}


# Make histograms of locations
plot.dist.by.country <- function(species, location, distances, output_dir) {
  hist_info <- hist(sea_loc_data$x, plot = FALSE)
  max_count <- max(hist_info$counts)
  max_x <- max(sea_loc_data$x)
  plot <- ggplot(sea_loc_data, aes(x = x, fill = country)) +
    geom_histogram(binwidth = 50, boundary = 0, position = "stack") +  # adjust the binwidth to personal preference
    labs(title = paste0("Frequencies of Sea distances/country for ", species," in ", location),
         x = "Sea distance in km", y = "Frequency of species") +
    theme_bw() +
    #scale_fill_brewer(palette = "Set1") +  # You can choose a different palette if you like
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 10),           # Set font size for axis numbers
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 14),   # Increase legend title size
          legend.text = element_text(size = 12),    # Increase legend text size
          legend.key.size = unit(1.5, "lines")) +   # Increase legend key size
    scale_x_continuous(breaks = seq(0, max_x*1.1, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, max_x*1.1), ylim = c(0,max_count*1.1)) # Use coord_cartesian for setting limits
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(output_dir, paste0(gsub(" ", "_", species), "_from_", location, "_by_country.png")), 
         plot = plot, width = 2400, height = 1200, units = "px", dpi = 300)
  return(plot)
}


# Make year categories (Used by plot.dist.by.year function)
assign_year_category <- function(year) {
  if (is.na(year)) {
    return(NA)   # return NA when year is not present
  }
  for (category in year_categories) {
    range <- as.numeric(unlist(strsplit(category, "-"))) # save years as numeric without "-"
    if (year >= range[1] & year < range[2]) {  # if the year falls into this category
      return(category) # return this category
    }
  }
  return(NA) # If year doesn't fall into any category, return NA
}


# Make histograms of year categories
plot.dist.by.year <- function(species, location, distances, output_dir) {
  hist_info <- hist(sea_loc_data$x, plot = FALSE)
  max_count <- max(hist_info$counts)
  plot <- ggplot(sea_loc_data, aes(x = x, fill = year_category)) +
    geom_histogram(binwidth = 50, boundary = 0, position = "stack") +  # adjust the binwidth to personal preference
    labs(title = paste0("Frequencies of Sea distances/year for ", species," in ", location),
         x = "Sea distance in km", y = "Frequency of species") +
    theme_bw() +
    scale_fill_brewer(palette = "YlOrRd", na.value = "black") + # You can choose a different palette if you like
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"), # set title font size, placement
          plot.margin = margin(0.3, 0.3, 0.4, 0.4, "cm"),
          axis.text = element_text(size = 10),           # Set font size for axis numbers
          axis.title = element_text(size = 20),
          legend.title = element_text(size = 14),   # set legend title size
          legend.text = element_text(size = 12),    # set legend text size
          legend.key.size = unit(1.5, "lines")) +   # set legend key size
    scale_x_continuous(breaks = seq(0, 7500, by = 250), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(0, 7500), ylim = c(0,max_count*1.1)) # Use coord_cartesian for setting limits
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  ggsave(filename = file.path(output_dir, paste0(gsub(" ", "_", species), "_from_", location, "_by_year.png")), 
         plot = plot, width = 2400, height = 1200, units = "px", dpi = 300)
  return(plot)
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
process_species_locations <- function(species, species_location, location_coordinates, 
                                    unique_coords, r, cost_matrix) {
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
  detected_locations <- names(which(
    species_location[.(species), .SD, .SDcols = -1] >= 1
  ))
  
  if (length(detected_locations) == 0) {
    message(sprintf("No detections found for species: %s", species))
    return(result_dt)
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
    # Combine all results
    combined <- rbindlist(valid_results, idcol = "location_idx")
    
    # Reshape and add to result_dt
    for (i in seq_along(valid_results)) {
      loc <- detected_locations[i]
      loc_data <- combined[location_idx == i]
      
      if (nrow(loc_data) > 0) {
        # Ensure we have the right number of rows
        if (nrow(loc_data) == nrow(result_dt)) {
          set(result_dt, j = paste0(loc, "_seaway"), value = loc_data$sea)
          set(result_dt, j = paste0(loc, "_geodesic"), value = loc_data$geo)
        } else {
          warning(sprintf("Mismatch in row counts for location: %s", loc))
        }
      }
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