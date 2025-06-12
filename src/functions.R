#' @title Alien Detective Functions
#' @description This script contains the custom functions used by the other scripts
#' @import data.table
#' @importFrom utils globalVariables
NULL

# This file contains various functions for the Alien Detective project
# including data fetching, processing, and visualization utilities

# Declare global variables to avoid R CMD check notes
globalVariables(
  c(
    # data.table special symbols
    ".", 
    ":=", 
    ".SD", 
    ".N", 
    ".I", 
    ".GRP", 
    ".BY", 
    ".EACHI",
    
    # Our variables
    "longitude_moved",
    "latitude_moved",
    "x_coord",
    "y_coord",
    "Observatory.ID",
    "Latitude",
    "Longitude",
    "location",
    "count",
    "sea_distances",
    "geodesic_distances"
  )
)


fetch_gbif_data <- function(species,
                            hasCoordinate = TRUE,
                            continent = "europe",
                            basisOfRecord = c("OBSERVATION", "MACHINE_OBSERVATION", "HUMAN_OBSERVATION", "MATERIAL_SAMPLE", "LIVING_SPECIMEN", "OCCURRENCE"),
                            fields = c("decimalLatitude", "decimalLongitude", "year", "month", "country"),
                            limit = 10000,
                            output_dir) {
  
  # Get GBIF data
  data_list <- rgbif::occ_search(
    scientificName = species,
    hasCoordinate = hasCoordinate,
    continent = continent,
    basisOfRecord = basisOfRecord,
    fields = fields,
    limit = limit
  )
  
  # Initialize result as data.table
  res <- data.table::data.table()
  
  # Process each data frame in the list
  for (i in seq_along(data_list)) {
    if (is.null(data_list[[i]]$data) || nrow(data_list[[i]]$data) == 0) {
      next
    }
    
    # Convert to data.table
    dt <- data.table::as.data.table(data_list[[i]]$data)
    
    # Add missing columns
    missing_cols <- setdiff(fields, names(dt))
    if (length(missing_cols) > 0) {
      dt[, (missing_cols) := NA]
    }
    
    # Add basisOfRecord
    dt[, basisOfRecord := names(data_list)[i]]
    
    # Combine results
    res <- data.table::rbindlist(list(res, dt), use.names = TRUE, fill = TRUE)
  }
  
  if (nrow(res) > 0) {
    # Reorder columns
    res <- res[, c(fields, "basisOfRecord"), with = FALSE]
    
    # Rename lat/long columns
    data.table::setnames(res,
                        old = c("decimalLatitude", "decimalLongitude"),
                        new = c("latitude", "longitude"),
                        skip_absent = TRUE)
    
    # Remove NA coordinates
    res <- res[!is.na(latitude) & !is.na(longitude)]
    
    return(res)
  } else {
    error_message <- paste0("No GBIF records found for species \"", species, "\"")
    message(error_message)
    return(NULL)
  }
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
  
  # Ensure data is a data.table
  if (!data.table::is.data.table(data)) {
    data <- data.table::as.data.table(data)
  }
  
  tryCatch({
    # Specify the PROJ4 string for WGS84
    proj4_crs <- sp::CRS("+init=EPSG:4326")
    
    # Create a data.table with coordinates (using := for in-place modification)
    coords_dt <- data.table::copy(data)
    
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
      x_coord = ifelse(is.na(longitude_moved), longitude, longitude_moved),
      y_coord = ifelse(is.na(latitude_moved), latitude, latitude_moved)
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
process_coords <- function(lat, lon) {
  if (!is_on_land(lat, lon)) {
    return(data.table::data.table(
      latitude_moved = lat, 
      longitude_moved = lon, 
      dist_moved = 0
    ))
  }
    
  moved <- move_to_sea(lat, lon)
  if (is.null(moved)) return(NULL)
    
  data.table::data.table(
    latitude_moved = moved$coords[2],
    longitude_moved = moved$coords[1],
    dist_moved = round(moved$dist/1000, 2)
  )
}

# Process species locations and calculate distances
#' Process Species Locations
#' 
#' @param species Character string of the species name to process
#' @param species_location data.table with species detection data
#' @param location_coordinates data.table with location coordinates
#' @param unique_coords data.table with unique coordinates to calculate distances to
#' @param r Raster layer for spatial operations
#' @param cost_matrix Cost matrix for distance calculations
#' @return Updated unique_coords data.table with distance columns added
#' @import data.table
#' @importFrom methods is
#' @export
process_species_locations <- function(species, species_location, location_coordinates, 
                                    unique_coords, r, cost_matrix) {
  # Input validation
  if (!is.character(species) || length(species) != 1) {
    stop("species must be a single character string")
  }
  
  # Convert inputs to data.table if needed
  if (!data.table::is.data.table(species_location)) {
    species_location <- data.table::as.data.table(species_location)
  }
  if (!data.table::is.data.table(location_coordinates)) {
    location_coordinates <- data.table::as.data.table(location_coordinates)
  }
  if (!data.table::is.data.table(unique_coords)) {
    unique_coords <- data.table::as.data.table(unique_coords)
  }
  
  # Set keys for faster joins
  species_col <- names(species_location)[1]
  data.table::setkeyv(species_location, species_col)
  data.table::setkeyv(location_coordinates, "Observatory.ID")
  
  # Get locations where species was detected
  detected_locations <- names(
    which(species_location[get(species_col) == species, 
                         -1, with = FALSE] >= 1)
  )
  
  # Process each detected location
  for (location in detected_locations) {
    # Get coordinates for this location
    coords <- location_coordinates[.(location), on = "Observatory.ID", nomatch = NULL]
    if (nrow(coords) != 1) {
      message("Could not retrieve coordinates for \"", location, "\"")
      next
    }
    
    # Convert coordinates to numeric (handling comma as decimal separator)
    coords[, `:=`(
      lat = as.numeric(gsub(",", ".", Latitude)),
      lon = as.numeric(gsub(",", ".", Longitude))
    )]
    
    # Calculate distances
    cat(">>> [DIST] Calculating distances to", species, "occurrences from", location, "\n")
    result <- calculate.distances(
      data = unique_coords,
      latitude = coords$lat[1],
      longitude = coords$lon[1],
      raster_map = r,
      cost_matrix = cost_matrix
    )
    
    # Add distance columns to unique_coords
    col_names <- paste0(location, c("_seaway", "_geodesic"))
    if (!is.null(result$sea_distances) || !is.null(result$geodesic_distances)) {
      unique_coords[, (col_names) := .(
        sea_distances = result$sea_distances,
        geodesic_distances = result$geodesic_distances
      )]
    } else {
      unique_coords[, (col_names) := .(
        sea_distances = NA_real_,
        geodesic_distances = NA_real_
      )]
    }
  }
  
  return(unique_coords)
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
