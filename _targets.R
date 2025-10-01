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
tar_option_set(
  packages = c(
  "data.table", "sf", "sp", "gdistance", "geodist", "raster", "fasterize",
  "ggplot2", "rnaturalearth", "rnaturalearthdata", "geosphere"
  ),
  resources = tar_resources(
    future = tar_resources_future(plan = "sequential")
  )
)

# Automatically source helper files in src/ that start with "functions_" ---------
# This keeps your helper code modular while letting {targets} track changes.
tar_source(files = list.files("src", pattern = "^functions_.*\\.R$", full.names = TRUE))

# ---- Parallel backend selection ----------------------------------------------
#
# How to run the pipeline
# ----------------------
#   • Sequential run (single R process)
#       targets::tar_make()
#
#   • Local workstation (default multisession via {future})
#       targets::tar_make_future()
#
#   • HPC/Slurm cluster using {clustermq}
#       Sys.setenv(ALIEN_PARALLEL = "hpc")   # or export ALIEN_PARALLEL=hpc in the job script
#       targets::tar_make_clustermq()
#
# The `ALIEN_PARALLEL` variable only tells this file which backend to configure;
# you must still call the matching tar_make_*() function from your R session or
# batch script. Adjust `workers` (local) or `n_jobs` (Slurm) in the section
# below to match your resources.
#
# Set the environment variable `ALIEN_PARALLEL=hpc` (e.g. in .Renviron or the
# job script) to dispatch targets as separate Slurm jobs via {clustermq}.
# Otherwise the pipeline falls back to a local multisession pool via {future}.
if (tolower(Sys.getenv("ALIEN_PARALLEL")) == "hpc") {
  suppressPackageStartupMessages(library(clustermq))
  # Configure clustermq for the scheduler (writes ~/.clustermq_slurm.tmpl once).
  try(cmq_enable_hpc(template = "slurm"), silent = TRUE)
  tar_option_set(resources = list(clustermq = list(n_jobs = 200))) # adjust as needed
} else {
  suppressPackageStartupMessages(library(future))
  plan(multisession, workers = max(1, parallel::detectCores() - 1))
}


# Pipeline -----------------------------------------------------------------------
list(
  # 1. Prepare workspace paths and parameters
  tar_target(paths, setup_workspace()),

  # 2. Read species/location tables and metadata
  tar_target(species_raw, get_species(paths, species_select = "Aurelia solida")),

  # 3. Augment species object with paths and safe filenames (adds safe_name, directory, gbif_file)
  tar_target(species, {
    species <- species_raw
    species$safe_name <- gsub(" ", "_", species$species_vec)
    species$directory <- file.path(paths$output_dir, species$safe_name)
    # Create any directories that do not yet exist (vector-safe)
    missing_dirs <- species$directory[!dir.exists(species$directory)]
    if (length(missing_dirs)) {
      for (d in missing_dirs) {
        if (!is.na(d) && nzchar(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
      }
    }
    species$gbif_file <- file.path(species$directory, paste0(species$safe_name, ".csv"))
    species
  }),

  # 3. Spatial data prerequisites (world map raster + cost matrix)
  tar_target(raster_map, get_world_map(paths)),
  tar_target(cost_matrix, get_cost_matrix(paths, raster_map)),

  # 4. Check input coordinates file (validate coordinates)
  tar_target(species_checked, {
    species_checked <- species
    species_checked$location_coordinates <- check_coordinates(species$location_coordinates, raster_map, cost_matrix)
    species_checked
  }),

  # 5. GBIF occurrences download/cache per species
  tar_target(gbif_occ, gbif_data(species_checked)),

  # 6. Coordinate processing & land-to-sea correction
  tar_target(coords, process_gbif_coords(gbif_occ, raster_map, cost_matrix),
    resources = tar_resources(
      future = tar_resources_future(plan = "multisession")
    )
  ),

  # 7. Determine missing locations that need distance computation
  tar_target(missing_locs, get_location(species_checked, gbif_occ)),

  # 8. Prepare distance data table
  tar_target(distances_dt,
             if (nrow(missing_locs) == 0) {
               NULL   # no missing locations, nothing to prepare
             } else {
               add_missing_dist(species_checked, missing_locs, coords)
             }),

  # 9. Compute seaway & geodesic distances (heavy; parallelisable)
  tar_target(dists,
             if (is.null(distances_dt)) {
               list(sea_distances = NULL, geodesic_distances = NULL)
             } else {
              calculate.distances(
                 data        = distances_dt,
                 raster_map  = raster_map,
                 cost_matrix = cost_matrix
               )
             }#,
             #resources = tar_resources(
             #  future = tar_resources_future(plan = "multisession")
             #)
 ),

  # 10. Merge distances back to data table
  tar_target(distances_merged,
             if (is.null(distances_dt)) {
               NULL
             } else {
               distances_dt[, dist_seaway   := if (!is.null(dists$sea_distances)) dists$sea_distances else NA_real_]
               distances_dt[, dist_geodesic := if (!is.null(dists$geodesic_distances)) dists$geodesic_distances else NA_real_]
               distances_dt
             }),

  # 11. Write species-specific CSV of occurrences with distances
  tar_target(distances_gbif_list,
             if (is.null(distances_merged)) {
               NULL
             } else {
               create_gbif_occurrences_file(species_checked, gbif_occ, distances_merged)
             }),

  # 12. Final plots (side-effect)
  tar_target(plotting,
             {
               if (!is.null(distances_gbif_list)) {
                 plot_data(distances_gbif_list)
               }
               NULL    # targets should return an object; NULL is fine for side-effects
             })
)