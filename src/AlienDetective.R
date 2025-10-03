
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

#  species_selection <-  "Aurelia solida"
species_selection <- c("Fibrocapsa japonica")
# species_selection <- c("Fibrocapsa japonica", "Amphibalanus amphitrite")
# species_selection <- c("Aurelia solida", "Acartia (Acanthacartia) tonsa", "Amphibalanus amphitrite", "Amphibalanus eburneus")
# species_selection <- c("Aurelia solida", "Acartia tonsa", "Amphibalanus amphitrite", "Amphibalanus eburneus")
#  species_selection <- c("Acartia (Acanthacartia) tonsa")
  
# # species_selection <-   c("Aurelia solida", "Acartia (Acanthacartia) tonsa", "Amphibalanus amphitrite", "Amphibalanus eburneus",
# species_selection <-   c("Aurelia solida", "Acartia (Acanthacartia) tonsa", "Amphibalanus amphitrite", "Amphibalanus eburneus",
#                         "Boccardia proboscidea", "Bonnemaisonia hamifera", "Botrylloides violaceus", "Bugula neritina",
#                         "Caprella mutica", "Caprella scaura", "Celleporaria brunnea", "Cephalothrix simula",
#                         "Cordylophora caspia", "Corella eumyota", "Corella sp.", "Crepidula fornicata", "Cutleria multifida",
#                         "Dasysiphonia japonica", "Dreisseninae sp.", "Eucheilota menoni", "Fenestrulina delicia",
#                         "Fibrocapsa japonica", "Ficopomatus enigmaticus", "Gonionemus vertens", "Haloa japonica",
#                         "Halothrix lumbricalis","Hemigrapsus takanoi", "Herdmania momus")

# species_selection <-   c("Aurelia solida", "Acartia (Acanthacartia) tonsa", "Amphibalanus amphitrite", "Amphibalanus eburneus",
#                          "Boccardia proboscidea", "Bonnemaisonia hamifera", "Botrylloides violaceus", "Bugula neritina")

  
# Initialize profiling
.init_profiling(
  script_name = "AlienDetective.R",
  workers = 1,
  data_source = "GBIF",
  species = species_selection,
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
    species_raw <- get_species(paths, species_select = species_selection)
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

####

# 6. GBIF occurrences download/cache (occ_search())
# profile_code("gbif data", {
#    gbif_occ <- gbif_data(species_checked)
# })

profile_code("download data", {
  gbif_occ_down <- download_gbif_data(species_raw$species_vec, user = NULL, pwd = NULL, email = NULL, continent = "europe", has_coords = TRUE)
  
  # gbif_data <- list()
  # for (spec in as.character(unique(gbif_occ_down[, "species"]))) {
  #   #print(spec)
  #   aa <- gbif_occ_down |> filter(species == spec)
  #   bb <- data.table(aa[, "decimalLatitude"],
  #                    aa[, "decimalLongitude"],
  #                    aa[, "year"],
  #                    aa[, "month"],
  #                    aa[, "countryCode"],
  #                    aa[, "basisOfRecord"])
  # 
  #   setnames(bb, old = c("decimalLongitude", "decimalLatitude"), new = c("longitude", "latitude"))
  # 
  #   gbif_data[[spec]] <- bb
  # }
  
  ##############
  setDT(gbif_occ_down)
  
  # rename once globally
  setnames(gbif_occ_down,
           old = c("decimalLongitude", "decimalLatitude"),
           new = c("longitude", "latitude"),
           skip_absent = TRUE)
  
  # keep only relevant columns (species must stay for split)
  gbif_occ_down <- gbif_occ_down[, .(species, latitude, longitude, year, month, countryCode, basisOfRecord)]
  
  # split into list by species
  gbif_data <- split(gbif_occ_down, by = "species", keep.by = FALSE)
  #######
  
  gbif_occ <- list()
  gbif_occ$gbif_occurrences <- gbif_data
  gbif_occ$species <- as.character(unique(gbif_occ_down[, "species"]))
})

# profile_code("load data", {
# 
#   print("loading data...")
#   # gbif_occ <- load_gbif_data(zip_path = "./data", key = "0036430-250920141307145") # 1 species
#   gbif_occ <- load_gbif_data(zip_path = "./data", key = "0036584-250920141307145") # 4 species
#   #gbif_occ <- load_gbif_data(zip_path = "./data", key = "0036452-250920141307145") # over 30 species
# 
#     # gbif_data <- list()
#     # for (spec in as.character(gbif_occ_down$species)) {
#     #   aa <- gbif_occ_down %>% filter(species == spec)
#     #   bb <- data.table(aa[, "decimalLatitude"],
#     #                    aa[, "decimalLongitude"],
#     #                    aa[, "year"],
#     #                    aa[, "month"],
#     #                    aa[, "county"],
#     #                    aa[, "basisOfRecord"])
#     #
#     #   setnames(bb, old = c("decimalLongitude", "decimalLatitude"), new = c("longitude", "latitude"))
#     #
#     #   gbif_data[[spec]] <- bb
#     # }
#     # gbif_occ <- list()
#     # gbif_occ$gbif_occurrences <- gbif_data
#     # gbif_occ$species <- as.character(unique(gbif_occ_down[, "species"]))
# 
#    print("data loaded")
# })

####

# 7. Coordinate processing & land-to-sea correction
profile_code("process gbif coords", {
  print("processing coordinates...")
    coords <- process_gbif_coords(gbif_occ, raster_map, cost_matrix)
    print("coordinates processed")
})

# 8. Determine missing locations that need distance computation
profile_code("get location", {
    print("getting missing locations...")
    missing_locs <- get_location(species_checked, gbif_occ)
    print("missing locations obtained")
})

# 9. Prepare distance data table
profile_code("prepare distance data table", {
    print("preparing distance data table...")
    if (nrow(missing_locs) == 0) {
        distances_dt <- NULL
    } else {
        distances_dt <- add_missing_dist(species_checked, missing_locs, coords)
    }
    print("distance data table prepared")
})

# 10. Compute seaway & geodesic distances
profile_code("compute distances", {
    print("computing distances...")
    if (is.null(distances_dt)) {
    dists <- list(sea_distances = NULL, geodesic_distances = NULL)
    } else {
    dists <- calculate.distances(
        data        = distances_dt,
        raster_map  = raster_map,
        cost_matrix = cost_matrix
    )
    }
    print("distances computed")
})

# 11. Merge distances back into table
profile_code("merge distances", {
    print("merging distances...")
    if (is.null(distances_dt)) {
        distances_merged <- NULL
    } else {
        distances_dt[, dist_seaway   := if (!is.null(dists$sea_distances)) dists$sea_distances else NA_real_]
        distances_dt[, dist_geodesic := if (!is.null(dists$geodesic_distances)) dists$geodesic_distances else NA_real_]
        distances_merged <- distances_dt
    }
    print("distances merged")
})

# 12. Write species-specific CSV of occurrences with distances
profile_code("write species-specific CSV", {
    print("writing species-specific CSV...")
    if (is.null(distances_merged)) {
        distances_gbif_list <- NULL
    } else {
        distances_gbif_list <- create_gbif_occurrences_file(
            species_checked, gbif_occ, distances_merged
        )
    }
    print("species-specific CSV written")
})

# 13. Final plots (side-effect only)
profile_code("final plots", {
    print("creating final plots...")
    if (!is.null(distances_gbif_list)) {
        plot_data(distances_gbif_list)
    }
    print("final plots done")
})

# Generate final profiling report without saving
generate_profiling_report(save_report = TRUE,
                          show_report = FALSE)
