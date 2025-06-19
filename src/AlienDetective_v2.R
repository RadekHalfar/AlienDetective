setwd("C:/Users/radek/Documents/IT4I_projects/BioFlow/AlienDetective")

# Clear workspace
rm(list = ls())

library(profvis)

# Start profiler: captures time & memory usage per function/line
pv <- profvis({

# load functions
lapply(c("setup", "computation", "plotting", "data_manipulation"),
       function(f) source(file.path("src", paste0("functions_", f, ".R"))))

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
species$location_coordinates <- check_coordinates(species$location_coordinates, r, cost_matrix)

# get gbif data
gbif_data <- gbif_data(species)

# get missing locations for distance computation
missing_locs <- get_location(species, gbif_data)

unique_coords <- process_gbif_coords(gbif_data, r, cost_matrix)

# create data table for row wise calculation of calculate.distances function
distances_dt <- add_missing_dist(species, missing_locs, unique_coords)

# compute distance for missing locations
dists <- calculate.distances(
  data = distances_dt,
  raster_map = r,
  cost_matrix = cost_matrix
)

# merge results
distances_dt[, `:=`(dist_seaway    = dists$sea_distances,
                    dist_geodesic  = dists$geodesic_distances)]

# create new gbif occurences file
distances_gbif_list <- create_gbif_occurrences_file(species, gbif_data, distances_dt)

# Plotting
plot_data(distances_gbif_list)

unlink("output", recursive = TRUE, force = TRUE)

}) # end profvis

print(pv)
