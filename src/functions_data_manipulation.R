create_gbif_occurrences_file <- function(species, gbif_data, distances_dt, write_csv = TRUE) {
    
    # Step 1: Reshape the data.table
    distances_dt <- data.table::dcast(
        distances_dt,
        species + latitude + longitude + latitude_moved + longitude_moved + dist_moved ~ Observatory.ID,
        value.var = c("dist_seaway", "dist_geodesic")
    )

    # Step 2: Rename new columns to match desired pattern: <Observatory.ID>_<metric>
    old_names <- names(distances_dt)
    new_names <- old_names

    # Rename "dist_seaway_X" to "X_seaway", and "dist_geodesic_X" to "X_geodesic"
    new_names <- gsub("^dist_seaway_(.*)$", "\\1_seaway", new_names)
    new_names <- gsub("^dist_geodesic_(.*)$", "\\1_geodesic", new_names)

    # Apply new names
    data.table::setnames(distances_dt, old = old_names, new = new_names)

    # create list of data.tables per species
    distances_list <- split(distances_dt, by = "species", keep.by = FALSE)

    # merge distances with gbif_data by species
    # distances_gbif_list <- Map(function(dist_dt, gbif_dt) {
    #     merge(gbif_dt, dist_dt, by = c("latitude", "longitude"), all.x = TRUE)
    # }, distances_list, gbif_data$gbif_occurrences)
    # 
    species_in_common <- intersect(names(distances_list), names(gbif_data$gbif_occurrences))
    
    # ensure we only use matching species
    distances_gbif_list <- lapply(species_in_common, function(sp) {
      dist_dt <- distances_list[[sp]]
      gbif_dt <- gbif_data$gbif_occurrences[[sp]]
      
      # convert to data.table just in case
      setDT(dist_dt)
      setDT(gbif_dt)
      
      # fast data.table join instead of merge()
      gbif_dt[dist_dt, on = .(latitude, longitude)]
    })
    
    # name the output list by species
    names(distances_gbif_list) <- species_in_common
    #######################################################
    
    # add gbif file name from species
    file_lookup <- setNames(species$gbif_file, species$species_vec)
    safe_name_lookup <- setNames(species$safe_name, species$species_vec)
    dir_lookup       <- setNames(species$directory, species$species_vec)
    
    # Enrich distances_gbif_list with metadata
    distances_gbif_list <- Map(function(dt, sp) {
      list(
        data      = dt,
        gbif_file = file_lookup[sp],
        safe_name = safe_name_lookup[sp],
        directory = dir_lookup[sp]
      )
    }, distances_gbif_list, names(distances_gbif_list))
    
    # write csv file
    if(write_csv){
      lapply(distances_gbif_list, function(entry) {
        dt <- entry$data
        filename <- entry$gbif_file
        
        data.table::fwrite(dt, file = filename, row.names = FALSE, col.names = TRUE)
      })
    }
    return(distances_gbif_list)
}

# create data.table with missing distance calculation
add_missing_dist <- function(species, missing_locs, unique_coords){
  # add information about lattitude and longitude to missing locations
  setkey(species$location_coordinates, Observatory.ID) # Set keys for fast join
  # Perform the join (left join)
  missing_locs <- species$location_coordinates[
    missing_locs,
    on = .(Observatory.ID = missing_locs)
  ]
  # create data table for row wise calculation of calculate.distances function
  # setnames(missing_locs, c("Longitude", "Latitude"),
  #          c("Longitude_missing_locs", "Latitude_missing_locs")) # Rename Longitude and Latitude in missing_locs
  setnames(missing_locs, c("Longitude", "Latitude"),
           c("Longitude_missing_locs", "Latitude_missing_locs")) # Rename Longitude and Latitude in missing_locs
  
  setkey(unique_coords, species)
  setkey(missing_locs, species)
  distances_dt <- missing_locs[unique_coords, allow.cartesian = TRUE] # Perform the join (many-to-many by species)
  return(distances_dt)
}