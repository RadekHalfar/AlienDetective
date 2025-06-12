# Clear workspace
rm(list = ls())

# AlienDetective.R
# Main script

# Load required packages
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)

#############
### SETUP ###
#############
setup_start <- Sys.time()
# Define number of cores  (place in comments for use on Windows OS)
# num_cores <- 4
# if (!is.numeric(num_cores) || num_cores <= 0 || num_cores != floor(num_cores)) {
#   stop("Number of cores must be a whole number!")
# }

# Set CRAN mirror for downloading packages on server
# options(repos = c(CRAN = "https://cloud.r-project.org"))

#setwd("~/AlienDetective")
source("src/functions.R")

# Reset graphics settings
graphics.off()

# NB! No packages loaded here, only installed if missing. Better to use explicit namespaces instead [e.g. raster::extract() rather than just extract()].
# That way it's easier to maintain the code and see which packages are actually required as development progresses, and you also avoid clashes between
# package namespaces, making sure that the correct function is always used regardless of which other packages the user has installed and loaded.
cat(">>> [INIT] Checking for required packages...\n")
# Core packages currently required
packages <- c(
  "rgbif", "sf", "sp", "gdistance", "geodist", "raster", "fasterize",
  "ggplot2", "rnaturalearth", "rnaturalearthdata", "geosphere"
)

# Optional extras if you re-enable parallelism or leaflet maps later:
# extras <- c("dplyr", "foreach", "doParallel", "leaflet", "htmlwidgets", "tidyr")
# packages <- c(packages, extras)

# Install missing packages
missing_pkgs <- setdiff(packages, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs)
}

# Load all packages
invisible(lapply(packages, library, character.only = TRUE, quietly = TRUE))

# Parse command-line arguments with sensible defaults
arg_defaults <- c(
  file.path("data", "Species_Location_NIS.csv"),      # species_location_path
  file.path("data", "Coordinates_NIS.csv"),            # location_coordinates_path
  file.path("data", "rasterized_land_polygons.rds"),   # rasterized_path
  file.path("data", "cost_matrix.rds"),                # cost_matrix_path
  "output"                                              # output_dir
)

args <- commandArgs(trailingOnly = TRUE)
# Use defaults if no args supplied
if (length(args) == 0) args <- arg_defaults
# Pad with defaults if fewer than expected
if (length(args) < length(arg_defaults)) {
  args <- c(args, arg_defaults[(length(args) + 1):length(arg_defaults)])
}

# Assign to named variables in the current environment
list2env(setNames(as.list(args), c("species_location_path", "location_coordinates_path", 
                          "rasterized_path", "cost_matrix_path", "output_dir")),
         envir = environment())

# Read species-location presence/absence matrix using data.table
species_location      <- data.table::fread(species_location_path,  sep = ";")
# If there are more than one row per species, keep only the first row for each species
species_location      <- species_location[!duplicated(species_location, by = names(species_location)[1])]
location_coordinates  <- data.table::fread(location_coordinates_path, sep = ";")
# Explicitly ensure both are data.table objects (one-time conversion)
data.table::setDT(species_location)
data.table::setDT(location_coordinates)

# Set keys for faster lookups
data.table::setkeyv(species_location, names(species_location)[1])
data.table::setkey(location_coordinates, "Observatory.ID")

# INSERT LIST OF NATIVE SPECIES TO REMOVE NATIVE SPECIES FROM DF LIST

# Subselect species to run the script for (optional). Can also be used to exclude species, e.g. known natives, by negating the which function
species_subset <- c("Aurelia solida")
species_location <- species_location[which(species_location$Specieslist %in% species_subset),]
#species_location <- species_location[c(2, 10, 57),] # Or subset a few species to try at random

# Create a simple character vector of species names for easy iteration
species_vec <- as.character(species_location[[1]])

required_columns <- c("decimalLatitude", "decimalLongitude", "year", "month", "country")

#########################
### MAP CONFIGURATION ###
#########################

# Load rasterized world map if it exists, otherwise load custom vector shapefile and rasterize it
cat(">>> [MAP] Loading world map...\n")

if(file.exists(rasterized_path)) {
  r <- readRDS(rasterized_path)
} else {
  # Read vector map as sf object
  #land_polygons <- sf::st_read(land_polygons_path)
  land_polygons <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")
  cat(">>> [MAP] Rasterizing land polygons...\n")
  # Create raster
  r <- raster::raster(raster::extent(-180, 180, -90, 90), crs = sp::CRS("+init=EPSG:4326"), resolution = 0.1)
  # Rasterize vector map using fasterize
  r <- fasterize::fasterize(land_polygons, r, field = NULL, fun = "max")
  # Set sea cells to value 1 and land cells to NA (Opposite of what fasterize outputs)
  r <- raster::calc(r, function(x) ifelse(is.na(x), 1, NA))
  saveRDS(r, rasterized_path)
  rm(land_polygons)
  cat(">>> [MAP] Rasterization done. Saved raster to \"", file.path(getwd(), rasterized_path), "\"\n")
}

if (file.exists(cost_matrix_path)) {
  cat(">>> [MAP] Loading cost matrix...\n")
  cost_matrix <- readRDS(cost_matrix_path)
} else {
  cat(">>> [MAP] Generating cost matrix...\n")
  # Create a transition object for adjacent cells
  cost_matrix <- gdistance::transition(r, transitionFunction = mean, directions = 16)
  # Set infinite costs to NA to prevent travel through these cells
  cost_matrix <- gdistance::geoCorrection(cost_matrix, type = "c", scl = FALSE)
  # Save transition matrix
  saveRDS(cost_matrix, file = cost_matrix_path)
  cat(">>> [MAP] Saved cost matrix to \"", file.path(getwd(), cost_matrix_path), "\"\n")
}

####################################
### Check input coordinates file ###
####################################
# Check if input coordinates are in sea, if not, move them to sea
cat(">>> [COORD] Checking if input coordinates are in sea ...\n")
for (i in 1:nrow(location_coordinates)) {
  loc_name <- location_coordinates$Observatory.ID[i]
  longitude <- as.numeric(gsub(",", ".", location_coordinates$Longitude[i]))
  latitude <- as.numeric(gsub(",", ".", location_coordinates$Latitude[i]))  
  #cat("Checking", loc_name,": latitude", latitude, ", longitude", longitude, "\n")
  
  if (is_on_land(latitude, longitude)) {
    #cat(loc_name, "is on land, searching nearest sea coordinates...\n")
    moved <- move_to_sea(latitude, longitude)
    
    if (is.null(moved)) {
      #cat("No valid sea coordinates found\n")
      message(loc_name, " is on land, no valid sea coordinates found")
    } else {
      # Update df with coordinates moved point
      location_coordinates$Longitude[i] <- moved$coords[1]
      location_coordinates$Latitude[i] <- moved$coords[2]
      dist <- round((moved$dist/1000), 2)
      #cat("Updated", loc_name, "to", location_coordinates$Latitude[i], ", ", location_coordinates$Longitude[i], "; moved", dist, "km.\n")
    }
  } else {
    #cat(loc_name, "is already in sea\n")
  }
  #cat("\n")
}
#cat(">>> [DONE] All coordinates updated to nearest sea point\n")

#############################
### DISTANCES CALCULATION ###
#############################
dist_start <- Sys.time()
# For non-parallel execution -> use "for" loop
# For parallel execution -> use "foreach" loop + parallel setup

# Setup parallelisation (place in comments for use on Windows OS)
# cluster <- makeCluster(num_cores)
# registerDoParallel(cluster)

# Process species one by one
for (species in species_vec) {
  # skip empty / NA entries
  if (is.na(species) || nchar(trimws(species)) == 0) next
   
  # Process GBIF data for the species
  safe_name <- gsub(" ", "_", species)
  species_dir <- file.path(output_dir, safe_name)
  
  # Ensure species_dir exists
  if (!dir.exists(species_dir)) {
    dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  gbif_file <- file.path(species_dir, paste0(safe_name, ".csv"))
  
  # Load or fetch GBIF data
  if (file.exists(gbif_file)) {
    #cat(">>> [GBIF] Loading GBIF data for", species, "\n")
    gbif_occurrences <- data.table::fread(gbif_file)
  } else {
    #cat(">>> [GBIF] Fetching GBIF data for", species, "\n")
    gbif_occurrences <- fetch_gbif_data(species, fields = required_columns)
    if (is.null(gbif_occurrences)) return(NULL)
    
    dir.create(species_dir, recursive = TRUE, showWarnings = FALSE)
    data.table::fwrite(gbif_occurrences, file = gbif_file)
  }
  
  # Ensure gbif_occurrences is a data.table (single call, no repeated checks)
  data.table::setDT(gbif_occurrences)
  
  # Process coordinates to ensure they're at sea
  #cat(">>> [GBIF] Ensuring GBIF occurrence coordinates are at sea\n")
  
  # Get unique coordinates using data.table
  unique_coords <- unique(gbif_occurrences[, .(latitude, longitude)])
  
  # Apply processing to all coordinates using data.table's := operator
  processed_list <- lapply(1:nrow(unique_coords), function(i) {
    process_coords(unique_coords$latitude[i], unique_coords$longitude[i])
  })
  
  # Combine results using rbindlist
  results <- data.table::rbindlist(processed_list, fill = TRUE)
  unique_coords <- cbind(unique_coords, results)
  
  # Report statistics
  moved_count <- sum(!is.na(unique_coords$dist_moved) & unique_coords$dist_moved > 0)
  failed_count <- sum(is.na(unique_coords$dist_moved))
  
  #cat(moved_count, "of", nrow(unique_coords), "coordinate pairs were moved to sea.\n")
  if (failed_count > 0) {
    #cat("Moving to sea failed for", failed_count, "coordinate pairs\n")
    message("Species ", species, ": moving to sea failed for ", failed_count, " coordinate pairs.")
  }
  
  # Process locations and calculate distances
  unique_coords <- process_species_locations(
    species = species,
    species_location = species_location,
    location_coordinates = location_coordinates,
    unique_coords = unique_coords,
    r = r,
    cost_matrix = cost_matrix
  )
  # Convert to data.table if not already
  if (!data.table::is.data.table(unique_coords)) {
    data.table::setDT(unique_coords)
  }
  
  # Joining using data.table merge
  gbif_occurrences <- unique_coords[gbif_occurrences, 
                                  on = c("latitude", "longitude"),
                                  nomatch = NA]
  
  # Save to csv file using fwrite
  data.table::fwrite(gbif_occurrences, file = gbif_file)
  cat("\n")
  return(TRUE)  # To avoid printing NULL in stdout
}

cat(">>> [DONE] Finished calculating distances for all species. \n")

dist_end <- Sys.time()
dist_time <- as.numeric(difftime(dist_end, dist_start, units = "secs"))

# Close the cluster   (place in comments for use on Windows OS)
# stopCluster(cluster)

##############
## Plotting ##
##############

plot_start <- Sys.time()

# Iterate over species names for plotting
for (species in species_vec) {
  if (is.na(species) || nchar(trimws(species)) == 0) next
  
  safe_name <- gsub(" ", "_", species)  # Change spaces to underscores for filenames
  species_dir <- file.path(output_dir, safe_name)
  gbif_file <- file.path(species_dir, paste0(safe_name, ".csv"))
  
  # If the file exists, execute following lines
  if (file.exists(gbif_file)) {
    # Read the species csv using the correct filename
    distance_df <- read.csv(gbif_file)
    #Change the format to a long format with pivot function
    long_df <- distance_df |>
      pivot_longer(
        cols = contains("_seaway") | contains("_geodesic"), #Select all columns with _seaway and _geodesic
        names_to = "location", #Change name of original columns (cols) to location column
        values_to = "x" #Put values (distances) in a column called x
      )|>
      separate(location, into = c("location", "DistanceType"), sep = "_") #separate the location column into location where the location represents a ARMS location and the DistanceType the type of distance
    long_sea <- long_df |>
      filter(grepl("seaway", DistanceType)) #Put all data of seaway into a dataframe
    long_geo <- long_df |>
      filter(grepl("geodesic", DistanceType)) #Put all data of geodesic into dataframe
    
  } else {
    warning("No output directory found for species \"", species, "\". Skipping plotting.")
    next
  }
  
  # Assign year categories
  year_categories <- c("1965-1985", "1985-1990", "1990-1995",
                       "1995-2000", "2000-2005", "2005-2010",
                       "2010-2015", "2015-2020", "2020-2025")
  long_sea$year_category <- sapply(long_sea$year, assign_year_category)
  long_geo$year_category <- sapply(long_sea$year, assign_year_category)
  
  # clean dataframe from rows with Inf and NA in them
  long_sea <- long_sea |>
    filter(!is.na(x), is.finite(x))
  
  #Make a graph of all locations where that species is found
  country_final_plot <- country.final(
    species = species,
    distances = long_sea$x,
    output_dir = species_dir)
  
  #Make a for loop that goes over every occurence location to make seperate graphs
  for (loc in unique(long_sea$location)) {
    sea_loc_data <- long_sea[long_sea$location == loc, ] #filter data on that specific location
    geo_loc_data <- long_geo[long_geo$location == loc, ]
    
    # Plot functions by location
    plot_dist_sea <- plot.dist.sea(
      species = species,
      location = loc,
      distances = sea_loc_data$x,
      output_dir = species_dir
    )
    
    plot_both <- plot.dist.both(
      species = species,
      location = loc,
      distances = combined_distances$x,
      output_dir = species_dir
    )
    
    plot_country <- plot.dist.by.country(
      species = species,
      location = loc,
      distances = sea_loc_data$x,
      output_dir = species_dir
    )
    
    plot_year <- plot.dist.by.year(
      species = species,
      location = loc,
      distances = sea_loc_data$x,
      output_dir = species_dir
    )
  }
}

plot_end <- Sys.time()
plot_time <- as.numeric(difftime(plot_end, plot_start, units = "secs"))
cat(">>> [DONE] Finished plotting for all species.\n")

end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, setup_start, units = "secs"))

cat(">>> [TIMING] Distance calculations completed in", round(dist_time, 2), "seconds.\n")
cat(">>> [TIMING] Plotting completed in", round(plot_time, 2), "seconds.\n")
cat(">>> [TIMING] Total runtime: ", round(total_time, 2), "seconds.\n")