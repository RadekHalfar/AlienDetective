setwd("C:/Users/radek/Documents/IT4I_projects/BioFlow/AlienDetective")

# Clear workspace
rm(list = ls())

# load functions
source("src/functions_setup.R")
source("src/functions_computation.R")
source("src/functions_plotting.R")

# setup workspace
paths <- setup_workspace()

# Get species to analyze
species <- get_species(paths)

# create safe name
species$safe_name <- gsub(" ", "_", species$species_vec)

# create output directories
species$directory <- file.path(paths$output_dir, species$safe_name)
sapply(species$directory, function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
})


# create gbif file name
species$gbif_file <- file.path(species$directory, paste0(species$safe_name, ".csv"))

# get world map
r <- get_world_map(paths)

# get cost matrix
cost_matrix <- get_cost_matrix(paths)

# Check input coordinates file
species$location_coordinates <- check_coordinates(species$location_coordinates)

# get gbif data
gbif_data <- gbif_data(species)

# get missing locations for distance computation
missing_locs <- get_location(species, gbif_data)

unique_coords <- process_gbif_coords(gbif_data)

# add information about lattitude and longitude to missing locations
setkey(species$location_coordinates, Observatory.ID) # Set keys for fast join

# Perform the join (left join)
missing_locs <- species$location_coordinates[
  missing_locs,
  on = .(Observatory.ID = missing_locs)
]

# create data table for row wise calculation of calculate.distances function 

setnames(missing_locs, c("Longitude", "Latitude"),
         c("Longitude_missing_locs", "Latitude_missing_locs")) # Rename Longitude and Latitude in missing_locs
setkey(unique_coords, species)
setkey(missing_locs, species)
merged_dt <- missing_locs[unique_coords, allow.cartesian = TRUE] # Perform the join (many-to-many by species)

# compute distance for missing locations
dists <- calculate.distances_v2(
  data = merged_dt,
  raster_map = r,
  cost_matrix = cost_matrix
)

# merge results
merged_dt[, `:=`(dist_seaway    = dists$sea_distances,
                 dist_geodesic  = dists$geodesic_distances)]

