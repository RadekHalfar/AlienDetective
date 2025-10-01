setup_workspace <- function() {
  
  # Load required packages
  if (!requireNamespace("data.table", quietly = TRUE)) {
    install.packages("data.table")
  }
  library(data.table)
  
  setup_start <- Sys.time()
  
  # Reset graphics settings
  graphics.off()
  
  packages <- c(
    "rgbif", "sf", "sp", "gdistance", "geodist", "raster", "fasterize",
    "ggplot2", "rnaturalearth", "rnaturalearthdata", "geosphere"
  )
  
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
  
  paths <- list(species_location_path = args[1],
                location_coordinates_path = args[2],
                rasterized_path = args[3],
                cost_matrix_path = args[4],
                output_dir = args[5])
  return(paths)

}

get_species <- function(paths, species_select = "all") {
  # Read species-location presence/absence matrix using data.table
  species_location      <- data.table::fread(paths$species_location_path,  sep = ";")
  # If there are more than one row per species, keep only the first row for each species
  species_location      <- species_location[!duplicated(species_location, by = names(species_location)[1])]
  location_coordinates  <- data.table::fread(paths$location_coordinates_path, sep = ";")
  # Explicitly ensure both are data.table objects (one-time conversion)
  data.table::setDT(species_location)
  data.table::setDT(location_coordinates)
  
  # Set keys for faster lookups
  data.table::setkeyv(species_location, names(species_location)[1])
  data.table::setkey(location_coordinates, "Observatory.ID")
  
  if(length(species_select) > 1) {
    species_location <- species_location[which(species_location$Specieslist %in% species_select),]
  } else if (species_select != "all") {
    species_location <- species_location[which(species_location$Specieslist %in% species_select),]
  }

  # Create a simple character vector of species names for easy iteration
  species_vec <- as.character(species_location[[1]])
  
  species <- list(species_location = species_location,
                  location_coordinates = location_coordinates,
                  species_vec = species_vec)

  return(species)
}
download_gbif_data <- function(species_vec, user = NULL, pwd = NULL, 
                               email = NULL, continent = NULL, 
                               has_coords = TRUE, output_dir = "./data") {

  if (is.null(user))  user  <- Sys.getenv("GBIF_USER")
  if (is.null(pwd))   pwd   <- Sys.getenv("GBIF_PWD")
  if (is.null(email)) email <- Sys.getenv("GBIF_EMAIL")
  
  if (interactive()) {
    if (user == "")  user  <- readline("Enter GBIF username: ")
    if (pwd == "")   pwd   <- readline("Enter GBIF password: ")
    if (email == "") email <- readline("Enter GBIF email: ")
  }
  
  if (user == "" || pwd == "" || email == "") {
    stop("GBIF credentials (user, pwd, email) must be set as environment variables or provided interactively.")
  }

  # Get the taxon key for Aurelia solida
  key <- name_backbone(name = "Aurelia solida")$usageKey
  # Build predicates for GBIF download
  predicates <- list(
    rgbif::pred_in("taxonKey", key),
    rgbif::pred("continent", continent),
    rgbif::pred("hasCoordinate", has_coords),
    rgbif::pred("basisOfRecord", "HUMAN_OBSERVATION")
  )

  # Use do.call with predicates as first arguments, then credentials as named arguments
  download_key <- do.call(
    rgbif::occ_download,
    c(predicates, list(user = user, pwd = pwd, email = email))
  )
  print(paste("GBIF download key:", download_key))
  rgbif::occ_download_wait(download_key)
  
  dwca_path <- rgbif::occ_download_get(download_key, path = output_dir)
  occ_data <- rgbif::occ_download_import(dwca_path)
  
  print(paste("GBIF data downloaded to:", dwca_path))
  
  return(occ_data)
}

split_gbif_occurrences <- function(occ_data) {
  # Ensure 'species' column exists (rename if needed)
  if (!"species" %in% names(occ_data)) {
    if ("scientificName" %in% names(occ_data)) {
      occ_data[, species := scientificName]
    } else {
      stop("No 'species' or 'scientificName' column found in occ_data.")
    }
  }
  # Get unique species in order of appearance
  species_vec <- unique(occ_data$species)
  # Split data.table by species
  gbif_occurrences <- split(occ_data, by = "species", keep.by = TRUE)
  # Ensure the list is ordered as species_vec
  gbif_occurrences <- gbif_occurrences[species_vec]
  
  # Return as a list with both components
  return(list(
    gbif_occurrences = gbif_occurrences,
    species = species_vec
  ))
}


load_gbif_data <- function(zip_path, key = NULL) {
  # Check if file exists
  if (!file.exists(zip_path)) {
    stop("GBIF zip file not found: ", zip_path)
  }
  # Get the GBIF download object
  gbif_download <- rgbif::occ_download_get(path = zip_path, key = key)
  # Import GBIF Darwin Core Archive using rgbif
  occ_data <- rgbif::occ_download_import(gbif_download)

  occ_data <- data.table::as.data.table(occ_data) # ensure data.table
  occ_data <- occ_data[, .(basisOfRecord,
                           latitude = decimalLatitude,
                           longitude = decimalLongitude,
                           year,
                           month,
                           country = countryCode,
                           species = scientificName)]

  return(split_gbif_occurrences(occ_data))
}

# Load rasterized world map if it exists, otherwise load custom vector shapefile and rasterize it
get_world_map <- function(paths){

  if(file.exists(paths$rasterized_path)) {
    r <- readRDS(paths$rasterized_path)
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
    saveRDS(r, paths$rasterized_path)
    rm(land_polygons)
  }
  return(r)
}

# get cost matrix
get_cost_matrix <- function(paths, raster_map){
  if (file.exists(paths$cost_matrix_path)) {
    cost_matrix <- readRDS(paths$cost_matrix_path)
  } else {
    # Create a transition object for adjacent cells
    cost_matrix <- gdistance::transition(raster_map, transitionFunction = mean, directions = 16)
    # Set infinite costs to NA to prevent travel through these cells
    cost_matrix <- gdistance::geoCorrection(cost_matrix, type = "c", scl = FALSE)
    # Save transition matrix
    saveRDS(cost_matrix, file = paths$cost_matrix_path)
  }
  return(cost_matrix)
}

# Check input coordinates file
check_coordinates <- function(location_coordinates, r, cost_matrix){
  for (i in 1:nrow(location_coordinates)) {
    loc_name <- location_coordinates$Observatory.ID[i]
    longitude <- as.numeric(gsub(",", ".", location_coordinates$Longitude[i]))
    latitude <- as.numeric(gsub(",", ".", location_coordinates$Latitude[i]))  
    if (is_on_land(latitude, longitude, r)) {
      moved <- move_to_sea(latitude, longitude, r, cost_matrix)
      
      if (is.null(moved)) {
        message(loc_name, " is on land, no valid sea coordinates found")
      } else {
        # Update df with coordinates moved point
        location_coordinates$Longitude[i] <- moved$coords[1]
        location_coordinates$Latitude[i] <- moved$coords[2]
        dist <- round((moved$dist/1000), 2)
      }
    } else {
    }
  }
  return(location_coordinates)
}

# load gbif data
gbif_data <- function(species, write_gbif_file = TRUE){
  
  # Columns we want to keep from GBIF
  required_columns <- c("decimalLatitude", "decimalLongitude", "year", "month", "country")
  
  # ------------------------------------------------------------------
  # Helper: given species name & path, load cached file or fetch online
  # ------------------------------------------------------------------
  load_or_fetch <- function(sp_name, file_path){
    if (file.exists(file_path)) {
      occ <- data.table::fread(file_path)
    } else {
      occ <- fetch_gbif_data(sp_name, fields = required_columns)
      if (is.null(occ)) {
        message("[GBIF] No occurrence records for ", sp_name, " – skipping.")
        return(NULL)
      }
      if (write_gbif_file) {
        data.table::fwrite(occ, file = file_path)
      }
    }
    data.table::setDT(occ) # ensure data.table
    return(occ)
  }
  
  occurrences_list <- mapply(load_or_fetch,
                             species$species_vec,
                             species$gbif_file,
                             SIMPLIFY = FALSE)
  
  names(occurrences_list) <- species$species_vec
  return(list(gbif_occurrences = occurrences_list,
              species = species$species_vec))
}

# Determine which locations still need distance computation
get_location <- function(species, gbif_data) {
  
  calc_missing <- function(sp, occ_tbl) {
    # Locations where the species is present (value >= 1)
    row          <- species$species_location[species$species_location$Specieslist == sp]
    pres_values  <- as.numeric(row[, -1])
    detected     <- names(row)[-1][pres_values >= 1 & !is.na(pres_values)]

    # Locations that already have distance columns computed
    processed    <- sub("_(seaway|geodesic)$", "",
                        grep("(_seaway|_geodesic)$", names(occ_tbl), value = TRUE))

    # Return only the locations still missing
    setdiff(detected, processed)
  }

  # Map over species vector and its corresponding GBIF occurrences table
  missing_locs <- Map(calc_missing,
                      species$species_vec,
                      gbif_data$gbif_occurrences)

  # Attach species names for easier downstream access
  names(missing_locs) <- species$species_vec

  # Remove NULLs and empty character vectors
  missing_locs <- Filter(function(x) !is.null(x) && length(x) > 0, missing_locs)

  # convert to data.table
  missing_locs_dt <- rbindlist(
    lapply(names(missing_locs), function(sp) {
      data.table(species = sp, missing_locs = missing_locs[[sp]])
    }),
    use.names = TRUE
  )

  return(missing_locs_dt)
  
}

# Extract unique latitude/longitude pairs for each species
process_gbif_coords <- function(gbif_data, r, cost_matrix) {

  res <- lapply(gbif_data$gbif_occurrences,
                function(tbl) unique(tbl[, c("latitude", "longitude")]))
  # Ensure list is named by species
  names(res) <- names(gbif_data$gbif_occurrences)

  # Remove NULL values
  res <- Filter(Negate(is.null), res)

  # Create data.table with unique latitude/longitude pairs for each species
  result_dt <- rbindlist(
    lapply(names(res), function(nm) {
      dt <- res[[nm]]
      dt[, species := nm]  # Add a new column "name"
      dt[, .(species, latitude, longitude)]  # Keep only desired columns
    }),
    use.names = TRUE
  )

  # Process coordinates, move points on land to sea
  result_dt <- process_coords(result_dt, r, cost_matrix)

  return(result_dt)
}
  
