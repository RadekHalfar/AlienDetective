setwd("C:/Users/radek/Documents/IT4I_projects/BioFlow/AlienDetective")

# Clear workspace
rm(list = ls())

# Load profiling functions
library("profiling")

# Initialize profiling
.init_profiling(
  script_name = "AlienDetective.R",
  workers = 4,
  data_source = "GBIF",
  species = "Aurelia solida",
  version = "1.0.0"
)

.start_profiling_step("Script initialization")

#library(profvis)

# Start profiler: captures time & memory usage per function/line
#pv <- profvis({

# load functions
lapply(c("setup", "computation", "plotting", "data_manipulation"),
       function(f) source(file.path("src", paste0("functions_", f, ".R"))))

# setup workspace
paths <- setup_workspace()

# End initialization profiling
.end_profiling_step("Script initialization")

# Read input data
profile_code("Read input data", {
  # Get species to analyze
  species <- get_species(paths)
})

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
profile_code("Get world map", {
  r <- get_world_map(paths)
})

# get cost matrix
profile_code("Get cost matrix", {
  cost_matrix <- get_cost_matrix(paths, r)
})

# Check input coordinates file
profile_code("Check input coordinates file", {
  species$location_coordinates <- check_coordinates(species$location_coordinates, r, cost_matrix)
})

# get gbif data
profile_code("gbif data", {
  gbif_data <- gbif_data(species)
})

# get missing locations for distance computation
profile_code("get location", {
  missing_locs <- get_location(species, gbif_data)
})

if(dim(missing_locs)[1] == 0){
  print("No missing locations")
} else {

  profile_code("process gbif coords", {
    unique_coords <- process_gbif_coords(gbif_data, r, cost_matrix)
  })

  # create data table for row wise calculation of calculate.distances function
  profile_code("add_missing_dist", {
    distances_dt <- add_missing_dist(species, missing_locs, unique_coords)
  })

  # compute distance for missing locations
  profile_code("calculate distances", {
    dists <- calculate.distances(
      data = distances_dt,
      raster_map = r,
      cost_matrix = cost_matrix
    )
  })

  # merge results
  profile_code("merge results", {
    distances_dt[, `:=`(dist_seaway    = dists$sea_distances,
                        dist_geodesic  = dists$geodesic_distances)]
  })

  # create new gbif occurences file
  profile_code("create gbif occurrences file", {
    distances_gbif_list <- create_gbif_occurrences_file(species, gbif_data, distances_dt)
  })


  # Plotting
  profile_code("plot data", {
    plot_data(distances_gbif_list)
  })
}


#unlink("output", recursive = TRUE, force = TRUE)

# Generate final profiling report without saving
generate_profiling_report(save_report = FALSE,
                          show_report = FALSE)

#print(pv)
