# _targets.R - Pipeline definition for AlienDetective
# --------------------------------------------------
# This script defines the reproducible workflow with the {targets} package.
# Each tar_target below corresponds to one logical step that existed in the
# previous sequential script.  Running `targets::tar_make()` will execute the
# pipeline, building only the steps whose inputs (code, data, or upstream
# targets) have changed.

library(targets)

# Packages needed by every target ------------------------------------------------
# Add/remove packages here as your analysis evolves.
tar_option_set(packages = c(
  "data.table", "sf", "sp", "gdistance", "geodist", "raster", "fasterize",
  "ggplot2", "rnaturalearth", "rnaturalearthdata", "geosphere", "fs"
))

# Automatically source helper files in src/ that start with "functions_" ---------
# This keeps your helper code modular while letting {targets} track changes.
tar_source(files = list.files("src", pattern = "^functions_.*\\.R$", full.names = TRUE))

# Pipeline -----------------------------------------------------------------------
list(
  # 1. Prepare workspace paths and parameters
  tar_target(paths, setup_workspace()),

  # 2. Read species/location tables and metadata
  tar_target(species_raw, get_species(paths)),

  # 3. Augment species object with paths and safe filenames (adds safe_name, directory, gbif_file)
  tar_target(species, {
    species <- species_raw
    species$safe_name <- gsub(" ", "_", species$species_vec)
    species$directory <- file.path(paths$output_dir, species$safe_name)
    if (!dir.exists(species$directory)) dir.create(species$directory, recursive = TRUE)
    species$gbif_file <- file.path(species$directory, paste0(species$safe_name, ".csv"))
    species
  }),

  # 3. Spatial data prerequisites (world map raster + cost matrix)
  tar_target(raster_map, get_world_map(paths)),
  tar_target(cost_matrix, get_cost_matrix(paths)),

  # 4. GBIF occurrences download/cache per species
  tar_target(gbif_occ, gbif_data(species)),

  # 5. Coordinate processing & land-to-sea correction
  tar_target(coords, process_gbif_coords(gbif_occ, raster_map, cost_matrix)),

  # 6. Determine missing locations that need distance computation
  tar_target(missing_locs, get_location(species, gbif_occ)),

  # 7. Prepare distance data table
  tar_target(distances_dt, add_missing_dist(species, missing_locs, coords)),

  # 8. Compute seaway & geodesic distances (heavy; parallelisable)
  tar_target(dists, calculate.distances(
    data        = distances_dt,
    raster_map  = raster_map,
    cost_matrix = cost_matrix
  )),

  # 9. Merge distances back to data table
  tar_target(distances_merged,
             {
               distances_dt[, `:=`(dist_seaway   = dists$sea_distances,
                                     dist_geodesic = dists$geodesic_distances)]
               distances_dt
             }),

  # 10. Write species-specific CSV of occurrences with distances
  tar_target(distances_gbif_list,
             create_gbif_occurrences_file(species, gbif_occ, distances_merged)),

  # 11. Generate plots and track all files in output directory
  tar_target(output_files,
             {
               dir.create(paths$output_dir, showWarnings = FALSE, recursive = TRUE)
               plot_data(distances_gbif_list)
               fs::dir_ls(paths$output_dir)   # return vector of file paths
             },
             format = "file")
)
