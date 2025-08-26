
# Load required packages --------------------------------------------------------
library(data.table)
library(sf)
library(sp)
library(gdistance)
library(geodist)
library(raster)
library(fasterize)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(geosphere)

library("profiling")

# Initialize profiling
.init_profiling(
  script_name = "AlienDetective.R",
  workers = 1,
  data_source = "GBIF",
  species = "Aurelia solida",
  version = "1.0.0"
)

# Source helper functions -------------------------------------------------------
helper_files <- list.files("src", pattern = "^functions_.*\\.R$", full.names = TRUE)
for (f in helper_files) source(f)

# -------------------------------------------------------------------------------
# Workflow
# -------------------------------------------------------------------------------

# 1. Prepare workspace paths and parameters
profile_code("setup workspace", {
    paths <- setup_workspace()
})

# 2. Read species/location tables and metadata
profile_code("get species", {
    species_raw <- get_species(paths)
})

# 3. Augment species object with paths and safe filenames
profile_code("augment species", {
    species <- species_raw
    species$safe_name <- gsub(" ", "_", species$species_vec)
    species$directory <- file.path(paths$output_dir, species$safe_name)
    missing_dirs <- species$directory[!dir.exists(species$directory)]
    if (length(missing_dirs)) {
        for (d in missing_dirs) {
            if (!is.na(d) && nzchar(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
        }
    }
    species$gbif_file <- file.path(species$directory, paste0(species$safe_name, ".csv"))
})

# 4. Spatial data prerequisites (world map raster + cost matrix)
profile_code("get world map", {
    raster_map  <- get_world_map(paths)
})

profile_code("get cost matrix", {
    cost_matrix <- get_cost_matrix(paths, raster_map)
})

# 5. Check input coordinates file
profile_code("check coordinates", {
    species_checked <- species
    species_checked$location_coordinates <- check_coordinates(
        species$location_coordinates, raster_map, cost_matrix
    )
})

# 6. GBIF occurrences download/cache
profile_code("gbif data", {
    gbif_occ <- gbif_data(species_checked)
})

# 7. Coordinate processing & land-to-sea correction
profile_code("process gbif coords", {
    coords <- process_gbif_coords(gbif_occ, raster_map, cost_matrix)
})

# 8. Determine missing locations that need distance computation
profile_code("get location", {
    missing_locs <- get_location(species_checked, gbif_occ)
})

# 9. Prepare distance data table
profile_code("prepare distance data table", {
    if (nrow(missing_locs) == 0) {
        distances_dt <- NULL
    } else {
        distances_dt <- add_missing_dist(species_checked, missing_locs, coords)
    }
})

# 10. Compute seaway & geodesic distances
profile_code("compute distances", {
    if (is.null(distances_dt)) {
    dists <- list(sea_distances = NULL, geodesic_distances = NULL)
    } else {
    dists <- calculate.distances(
        data        = distances_dt,
        raster_map  = raster_map,
        cost_matrix = cost_matrix
    )
    }
})

# 11. Merge distances back into table
profile_code("merge distances", {
    if (is.null(distances_dt)) {
        distances_merged <- NULL
    } else {
        distances_dt[, dist_seaway   := if (!is.null(dists$sea_distances)) dists$sea_distances else NA_real_]
        distances_dt[, dist_geodesic := if (!is.null(dists$geodesic_distances)) dists$geodesic_distances else NA_real_]
        distances_merged <- distances_dt
    }
})

# 12. Write species-specific CSV of occurrences with distances
profile_code("write species-specific CSV", {
    if (is.null(distances_merged)) {
        distances_gbif_list <- NULL
    } else {
        distances_gbif_list <- create_gbif_occurrences_file(
            species_checked, gbif_occ, distances_merged
        )
    }
})

# 13. Final plots (side-effect only)
profile_code("final plots", {
    if (!is.null(distances_gbif_list)) {
        plot_data(distances_gbif_list)
    }
})

# Generate final profiling report without saving
generate_profiling_report(save_report = TRUE,
                          show_report = FALSE)
