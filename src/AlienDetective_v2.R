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

# Determine which locations still need distance computation
# species_row <- species_location[Specieslist == species]
# presence_vals <- as.numeric(species_row[, -1])
# detected_locs <- names(species_row)[-1][!is.na(presence_vals) & presence_vals >= 1]
# 
# existing_dist_cols <- grep("(_seaway|_geodesic)$", names(gbif_occurrences), value = TRUE)
# processed_locs <- unique(sub("_(seaway|geodesic)$", "", existing_dist_cols))
# 
# missing_locs <- setdiff(detected_locs, processed_locs)



