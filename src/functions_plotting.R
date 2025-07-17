# Generic helper --------------------------------------------------------------
make_hist_plot <- function(data,
                           x_col,
                           fill_col = NULL,
                           title = "",
                           binwidth = 50,
                           breaks_by = 250,
                           x_label = "Distance (km)",
                           y_label = "Frequency") {
  stopifnot(is.data.frame(data), x_col %in% names(data))

  max_x <- max(data[[x_col]], na.rm = TRUE)

  # Derive maximum bin count for y-axis scaling
  # Use the same binwidth as geom_histogram so the limit is accurate.
  hist_breaks <- seq(0, max_x * 1.1, by = binwidth)
  #max_count <- max(stats::hist(data[[x_col]], breaks = hist_breaks, plot = FALSE)$counts, na.rm = TRUE)
  max_count <- max(hist(data[[x_col]], breaks = hist_breaks, plot = FALSE)$counts, na.rm = TRUE)

  # Dynamically build aesthetic mapping
  mapping <- ggplot2::aes(x = .data[[x_col]])
  if (!is.null(fill_col)) {
    mapping <- ggplot2::aes(x = .data[[x_col]], fill = .data[[fill_col]])
  }

  ggplot2::ggplot(data, mapping) +
    ggplot2::geom_histogram(binwidth = binwidth, boundary = 0, alpha = 0.7,
                            position = "stack", colour = "#e9ecef") +
    ggplot2::labs(title = title, x = x_label, y = y_label, fill = fill_col) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::scale_x_continuous(breaks = seq(0, max_x * 1.1, by = breaks_by),
                                expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(0, max_x * 1.1), ylim = c(0, max_count * 1.1))
}

# -----------------------------------------------------------------------------
# Plot wrappers (all rely on make_hist_plot) ----------------------------------
# -----------------------------------------------------------------------------

country.final <- function(species, data, output_dir) {
  plot <- make_hist_plot(
    data      = data,
    x_col     = "x",
    fill_col  = "location",
    title     = sprintf("Sea-route distances for %s (all ARMS locations)", species)
  )
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir, sprintf("%s_from_all_ARMS_locations.png", gsub(" ", "_", species))),
                  plot, width = 2400, height = 1200, units = "px", dpi = 300)
  invisible(plot)
}

plot.dist.sea <- function(species, location, data, output_dir) {
  plot <- make_hist_plot(
    data      = data,
    x_col     = "x",
    title     = sprintf("Sea-route distances for %s from %s", species, location)
  )
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir, sprintf("%s_from_%s.png", gsub(" ", "_", species), location)),
                  plot, width = 2000, height = 1200, units = "px", dpi = 300)
  invisible(plot)
}

plot.dist.both <- function(species, location, data, output_dir) {
  plot <- make_hist_plot(
    data      = data,
    x_col     = "x",
    fill_col  = "DistanceType",
    title     = sprintf("Sea vs geodesic distances for %s (%s)", species, location)
  ) + ggplot2::labs(fill = "DistanceType")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir, sprintf("%s_from_%s_seadist&geodesic.png", gsub(" ", "_", species), location)),
                  plot, width = 2400, height = 1200, units = "px", dpi = 300)
  invisible(plot)
}

plot.dist.by.country <- function(species, location, data, output_dir) {
  plot <- make_hist_plot(
    data      = data,
    x_col     = "x",
    fill_col  = "country",
    title     = sprintf("Sea-route distances by country for %s (%s)", species, location)
  )
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir, sprintf("%s_from_%s_by_country.png", gsub(" ", "_", species), location)),
                  plot, width = 2400, height = 1200, units = "px", dpi = 300)
  invisible(plot)
}

plot.dist.by.year <- function(species, location, data, output_dir) {
  plot <- make_hist_plot(
    data      = data,
    x_col     = "x",
    fill_col  = "year_category",
    title     = sprintf("Sea-route distances by year for %s (%s)", species, location)
  ) +
    ggplot2::scale_fill_brewer(palette = "YlOrRd", na.value = "black") +
    ggplot2::labs(fill = "Year")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir, sprintf("%s_from_%s_by_year.png", gsub(" ", "_", species), location)),
                  plot, width = 2400, height = 1200, units = "px", dpi = 300)
  invisible(plot)
}

# plot data
plot_data <- function(distances_gbif_list){
  brks <- c(1965, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2025)
  labs <- paste(head(brks, -1), tail(brks, -1), sep = "-")
  
  invisible(lapply(names(distances_gbif_list), function(species) {
    if (is.na(species) || nchar(trimws(species)) == 0) return(NULL)
    
    entry <- distances_gbif_list[[species]]
    safe_name   <- entry$safe_name
    species_dir <- entry$directory
    distance_dt <- copy(entry$data)  # Copy to avoid modifying original
    
    dist_cols <- grep("(_seaway|_geodesic)$", names(distance_dt), value = TRUE)
    if (length(dist_cols) == 0) {
      warning(sprintf("No distance columns found for species '%s'. Skipping.", species))
      return(NULL)
    }
    
    # Ensure numeric and melt
    distance_dt[, (dist_cols) := lapply(.SD, as.numeric), .SDcols = dist_cols]
    long_dt <- melt(distance_dt,
                    measure.vars = dist_cols,
                    variable.name = "loc_type",
                    value.name = "x",
                    variable.factor = FALSE)
    
    long_dt[, `:=`(
      DistanceType = fifelse(grepl("_seaway$", loc_type), "seaway", "geodesic"),
      location     = sub("_(seaway|geodesic)$", "", loc_type)
    )][, loc_type := NULL]
    
    long_dt[, year_category := cut(year, breaks = brks, labels = labs, right = FALSE)]
    
    long_sea <- long_dt[DistanceType == "seaway" & !is.na(x) & is.finite(x)]
    long_geo <- long_dt[DistanceType == "geodesic"]
    
    # Overall plot
    country.final(
      species = species,
      data = long_sea,
      output_dir = species_dir
    )
    
    # Per-location plots
    locs <- unique(long_sea$location)
    lapply(locs, function(loc) {
      sea_loc_data <- long_sea[location == loc]
      geo_loc_data <- long_geo[location == loc]
      combined_distances <- rbindlist(list(sea_loc_data, geo_loc_data), use.names = TRUE)
      
      plot.dist.sea(species, loc, sea_loc_data, species_dir)
      plot.dist.both(species, loc, combined_distances, species_dir)
      plot.dist.by.country(species, loc, sea_loc_data, species_dir)
      plot.dist.by.year(species, loc, sea_loc_data, species_dir)
      NULL
    })
    
    NULL  # Prevents lapply from returning results
  }))
}