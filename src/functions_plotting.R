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
  mapping <- ggplot2::aes_string(x = x_col)
  if (!is.null(fill_col)) {
    mapping <- ggplot2::aes_string(x = x_col, fill = fill_col)
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
