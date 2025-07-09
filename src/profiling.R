#' @importFrom magrittr `%>%`
#' @importFrom utils installed.packages install.packages browseURL
#' @importFrom stats reorder
#' @importFrom methods new

# Profiling functions for AlienDetective scripts
# Enhanced with visualizations similar to modeling.R

# Global variable for profiling data
profiling_data <- NULL

# Check and install required packages
.ensure_packages <- function() {
  required_pkgs <- c("pryr", "tictoc", "knitr", "kableExtra", "ggplot2", "plotly", "DT")
  installed <- rownames(installed.packages())
  to_install <- setdiff(required_pkgs, installed)
  
  if (length(to_install) > 0) {
    message("Installing required packages: ", paste(to_install, collapse = ", "))
    install.packages(to_install, repos = "https://cloud.r-project.org/")
  }
  
  # Load all required packages
  invisible(lapply(required_pkgs, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

# Initialize profiling environment
#
# @param ... Named arguments containing script metadata (e.g., script_name, workers, etc.)
# @return Invisibly returns TRUE if successful
.init_profiling <- function(...) {
  # Ensure all required packages are installed and loaded
  .ensure_packages()
  
  # Create empty profiling data frame
  assign("profiling_data", 
         data.frame(
           Step = character(0),
           Start_Time = as.POSIXct(character(0)),
           End_Time = as.POSIXct(character(0)),
           Duration_sec = numeric(0),
           Memory_Used = character(0),
           Memory_Bytes = numeric(0),
           stringsAsFactors = FALSE
         ),
         envir = .GlobalEnv)
  
  # Store metadata
  metadata <- list(...)
  if (length(metadata) > 0) {
    assign("profiling_metadata", metadata, envir = .GlobalEnv)
  } else {
    assign("profiling_metadata", list(script_name = basename(rstudioapi::getSourceEditorContext()$path)), 
           envir = .GlobalEnv)
  }
  
  # Initialize timing
  assign("profiling_start_time", Sys.time(), envir = .GlobalEnv)
  
  # Initialize memory tracking
  assign("profiling_memory_usage", numeric(), envir = .GlobalEnv)
  
  # Initialize tictoc
  tictoc::tic("Total script execution")
  
  invisible(TRUE)
}

# Start profiling a step
.start_profiling_step <- function(step_name) {
  if (!exists("profiling_data", envir = .GlobalEnv)) {
    .init_profiling()
  }
  
  # Start timer for this step with a unique name
  tictoc::tic(paste0("step_", step_name))
  
  # Record memory before
  mem_before <- pryr::mem_used()
  
  # Store step info in global environment
  assign(paste0("step_", step_name), 
         list(name = step_name, 
              start_time = Sys.time(),
              mem_before = mem_before),
         envir = .GlobalEnv)
}

# End profiling a step
.end_profiling_step <- function(step_name) {
  # Get step info
  step_var <- paste0("step_", step_name)
  if (!exists(step_var, envir = .GlobalEnv)) {
    warning(paste("Profiling step not found:", step_name))
    return()
  }
  
  step_info <- get(step_var, envir = .GlobalEnv)
  
  # Record memory after and calculate usage
  mem_after <- pryr::mem_used()
  mem_used_bytes <- as.numeric(mem_after - step_info$mem_before)
  # Ensure we're working with a single value for formatting
  if (length(mem_used_bytes) > 1) {
    mem_used_bytes <- sum(mem_used_bytes, na.rm = TRUE)
  }
  mem_used_human <- format(structure(mem_used_bytes, class = "object_size"), 
                          units = "auto")
  
  # Stop timer for this step
  tictoc::toc(quiet = TRUE)
  
  # Add to profiling data
  new_row <- data.frame(
    Step = step_name,
    Start_Time = step_info$start_time,
    End_Time = Sys.time(),
    Duration_sec = as.numeric(difftime(Sys.time(), step_info$start_time, units = "secs")),
    Memory_Used = mem_used_human,
    Memory_Bytes = mem_used_bytes,
    stringsAsFactors = FALSE
  )
  
  profiling_data <<- rbind(profiling_data, new_row)
  
  # Clean up
  rm(list = step_var, envir = .GlobalEnv)
  
  return(new_row)
}

# Generate summary statistics for the report
.generate_summary_stats <- function(profiling_data) {
  if (nrow(profiling_data) == 0) return(NULL)
  
  total_time <- sum(profiling_data$Duration_sec, na.rm = TRUE)
  avg_time <- mean(profiling_data$Duration_sec, na.rm = TRUE)
  max_time <- max(profiling_data$Duration_sec, na.rm = TRUE)
  max_mem <- max(profiling_data$Memory_Bytes, na.rm = TRUE)
  
  data.frame(
    Metric = c("Total Execution Time (s)", 
              "Average Step Time (s)",
              "Longest Step (s)",
              "Peak Memory Usage"),
    Value = c(
      round(total_time, 2),
      round(avg_time, 2),
      round(max_time, 2),
      format(structure(max_mem, class = "object_size"), units = "auto")
    ),
    stringsAsFactors = FALSE
  )
}

# Generate interactive plots for the report
.generate_profiling_plots <- function(profiling_data) {
  if (nrow(profiling_data) == 0) return("")
  
  # Create a temporary file for the plot
  plot_file <- tempfile(fileext = ".png")
  
  # Time per step plot
  p1 <- ggplot2::ggplot(profiling_data, 
                       ggplot2::aes(x = stats::reorder(.data$Step, .data$Duration_sec), 
                                   y = .data$Duration_sec, 
                                   fill = .data$Duration_sec)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Execution Time by Step", 
                 x = "Step", 
                 y = "Time (seconds)") +
    ggplot2::theme_minimal() + 
    ggplot2::theme(text = element_text(size = 18),
                   legend.position = "none")
  
  # Memory usage plot
  p2 <- ggplot2::ggplot(profiling_data, 
                       ggplot2::aes(x = stats::reorder(.data$Step, .data$Memory_Bytes), 
                                   y = .data$Memory_Bytes/1024/1024,  # Convert to MB
                                   fill = .data$Memory_Bytes)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Memory Usage by Step", 
                 #x = "Step", 
                 x = "", 
                 y = "Memory (MB)") +
    ggplot2::theme_minimal() + 
    ggplot2::theme(text = element_text(size = 18),
                   legend.position = "none")
  
  # Combine plots
  combined_plot <- gridExtra::grid.arrange(p1, p2, ncol = 2)
  
  # Save to temp file with smaller dimensions
  ggplot2::ggsave(plot_file, combined_plot, width = 25, height = 5, dpi = 300)
  
  # Return base64 encoded image
  knitr::image_uri(plot_file)
}

# Generate HTML report
# Generate a profiling report
#
# @param output_dir Directory to save the report (default: "profiling_reports")
# @param file_name Name of the output HTML file without extension. If NULL, a timestamp-based name will be used.
# @param show_report Logical indicating whether to open the report in browser (default: FALSE)
# @param save_report Logical indicating whether to save the report to a file (default: TRUE)
# @return The HTML content as a character string (invisibly if saved to file)
generate_profiling_report <- function(output_dir = "profiling_reports", file_name = NULL, 
                                    show_report = FALSE, save_report = TRUE) {
  if (!exists("profiling_data", envir = .GlobalEnv) || nrow(profiling_data) == 0) {
    warning("No profiling data available to generate report")
    return(invisible(NULL))
  }
  
  # Stop overall timer and calculate total time
  total_time <- tictoc::toc(quiet = TRUE)
  total_seconds <- if (is.list(total_time)) {
    round(total_time$toc - total_time$tic, 2)
  } else {
    NA_real_
  }
  

  
  # Generate report filename and timestamp
  report_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  display_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  if (is.null(file_name)) {
    file_name <- paste0("profiling_report_", report_timestamp, ".html")
  } else if (!grepl("\\.html$", tolower(file_name))) {
    file_name <- paste0(file_name, ".html")
  }
  
  report_path <- file.path(output_dir, file_name)
  
  # Generate summary statistics
  summary_stats <- .generate_summary_stats(profiling_data)
  
  # Add number of steps to summary
  summary_stats <- rbind(summary_stats, 
                        data.frame(
                          Metric = "Number of Steps",
                          Value = as.character(nrow(profiling_data)),
                          stringsAsFactors = FALSE
                        )
  )
  
  # Generate plots
  plot_data <- tryCatch({
    .generate_profiling_plots(profiling_data)
  }, error = function(e) {
    warning("Failed to generate plots: ", e$message)
    NULL
  })
  
  # Format profiling data for display
  display_data <- profiling_data
  
  # Safely format memory values, handling vector inputs
  format_memory <- function(bytes) {
    if (length(bytes) == 0) return(character(0))
    sapply(bytes, function(x) {
      if (is.na(x) || !is.numeric(x)) return("N/A")
      format(structure(x, class = "object_size"), units = "auto")
    })
  }
  
  display_data$Memory_Used <- format_memory(display_data$Memory_Bytes)
  display_data <- display_data[, c("Step", "Duration_sec", "Memory_Used")]
  
  # Get metadata if it exists
  metadata_html <- ""
  if (exists("profiling_metadata", envir = .GlobalEnv) && length(profiling_metadata) > 0) {
    metadata_html <- '<div class="card mb-4">\n      <div class="card-header">\n        <i class="fas fa-info-circle me-2"></i>Script Information\n      </div>\n      <div class="card-body">\n        <div class="row g-3">\n    '
    for (i in seq_along(profiling_metadata)) {
      name <- names(profiling_metadata)[i]
      if (is.null(name) || name == "") name <- paste0("Parameter ", i)
      value <- profiling_metadata[[i]]
      if (length(value) > 1) value <- paste(value, collapse = ", ")
      
      metadata_html <- paste0(metadata_html, '
        <div class="col-md-4">
          <div class="d-flex align-items-center">
            <strong class="me-2">', tools::toTitleCase(gsub("_", " ", name)), ':</strong>
            <span>', value, '</span>
          </div>
        </div>')
    }
    
    metadata_html <- paste0(metadata_html, '\n        </div>\n      </div>\n    </div>')
  }
  
  # Create HTML report
  html_content <- paste0(
    '<!DOCTYPE html>\n    <html>\n    <head>\n      <title>Profiling Report - ', report_timestamp, '</title>
      <meta charset="UTF-8">
      <meta name="viewport" content="width=device-width, initial-scale=1.0">
      <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
      <link href="https://cdn.datatables.net/1.11.5/css/dataTables.bootstrap5.min.css" rel="stylesheet">
      <style>
        body { font-family: "Segoe UI", Tahoma, Geneva, Verdana, sans-serif; line-height: 1.5; color: #333; font-size: 0.95rem; }
        .header { background-color: #2c3e50; color: white; padding: 1.5rem 0; margin-bottom: 1.5rem; border-radius: 5px; }
        .card { margin-bottom: 1.5rem; box-shadow: 0 2px 4px rgba(0,0,0,0.05); border: 1px solid #eee; }
        .card-header { background-color: #f8f9fa; font-weight: 600; padding: 0.75rem 1.25rem; }
        .card-body { padding: 1rem; }
        .summary-card { background-color: #f8f9fa; border-left: 4px solid #007bff; }
        .plot-container { margin: 1rem 0; text-align: center; }
        .plot-img { max-width: 90%; height: auto; }
        .data-table { margin-top: 1.5rem; font-size: 0.9rem; }
        .footer { margin-top: 2rem; padding: 1rem 0; border-top: 1px solid #eee; color: #6c757d; font-size: 0.85em; }
      </style>
    </head>
    <body>
      <div class="container-fluid">
        <div class="header text-center">
          <h1>Profiling Report</h1>\n          <p class="lead mb-0">', display_timestamp, '</p>
        </div>
        
        <!-- Script Information Section -->
        ', metadata_html, '
        
        <!-- Performance Summary Section -->
        
        <div class="row mb-4">
          <div class="col-12">
            <h2 class="h4 mb-3"><i class="fas fa-chart-line me-2"></i>Performance Summary</h2>
          </div>
          <div class="col-12">
            <div class="row g-3">',
                if (!is.null(summary_stats)) {
                  paste0(sapply(seq_len(nrow(summary_stats)), function(i) {
                    paste0(
                      '<div class="col">
                        <div class="card summary-card h-100 border-0 shadow-sm">
                          <div class="card-body text-center p-3">
                            <h6 class="card-subtitle mb-2 text-muted">', summary_stats$Metric[i], '</h6>
                            <p class="card-title mb-0 fs-4 fw-bold">', summary_stats$Value[i], '</p>
                          </div>
                        </div>
                      </div>'
                    )
                  }), collapse = '\n')
                },
                '</div>
              </div>
            </div>
          </div>
        </div>',
        
        if (!is.null(plot_data)) {
          paste0(
            '<div class="row">
              <div class="col-md-12">
                <div class="card">
                  <div class="card-header">
                    <i class="fas fa-chart-bar me-2"></i>Performance Metrics
                  </div>
                  <div class="card-body">
                    <div class="plot-container">
                      <img src="', plot_data, '" class="img-fluid" alt="Performance Plots">
                    </div>
                  </div>
                </div>
              </div>
            </div>'
          )
        },
        
        '<div class="row">
          <div class="col-md-12">
            <div class="card">
              <div class="card-header">
                <i class="fas fa-table me-2"></i>Detailed Profiling Data
              </div>
              <div class="card-body">
                <div class="table-responsive">',
                {
                  tbl <- knitr::kable(display_data, format = "html", escape = FALSE)
                  tbl <- kableExtra::kable_styling(tbl, 
                                                 bootstrap_options = c("striped", "hover", "responsive"),
                                                 full_width = TRUE,
                                                 position = "left")
                  kableExtra::row_spec(tbl, 0, bold = TRUE, color = "white", background = "#2c3e50")
                },
                '</div>
              </div>
            </div>
          </div>
        </div>
        
        <div class="footer text-center">\n          <p>Report generated on ', display_timestamp, '</p>\n        </div>
      </div>
      
      <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
      <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
      <script src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.min.js"></script>
      <script src="https://cdn.datatables.net/1.11.5/js/dataTables.bootstrap5.min.js"></script>
      <script>
        $(document).ready(function() {
          $("table").DataTable({
            pageLength: 10,
            lengthMenu: [[10, 25, 50, -1], [10, 25, 50, "All"]],
            order: [[1, "desc"]]
          });
        });
      </script>
    </body>
    </html>'
  )
  
  # Add interactive highlighting to the HTML content
  html_content <- sub(
    "</script>",
    "</script>
      <style>
        tr.highlight { background-color: #e8f4f8 !important; }
      </style>
      <script>
        // Add row highlighting on click
        document.addEventListener('DOMContentLoaded', function() {
          const rows = document.querySelectorAll('tbody tr');
          rows.forEach(row => {
            row.addEventListener('click', function() {
              this.classList.toggle('highlight');
            });
          });
        });
      </script>",
    html_content
  )
  
  # Save to file if requested
  if (isTRUE(save_report)) {
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    writeLines(html_content, report_path)
    
    # Open in browser if requested
    if (isTRUE(show_report) && interactive()) {
      utils::browseURL(report_path)
    }
    message("Profiling report saved to: ", normalizePath(report_path))
    return(invisible(html_content))
  } else {
    # Return HTML content directly if not saving to file
    if (isTRUE(show_report) && interactive()) {
      temp_file <- tempfile(fileext = ".html")
      writeLines(html_content, temp_file)
      utils::browseURL(temp_file)
      message("Temporary profiling report opened in browser")
    }
    return(html_content)
  }
}

# Helper function to profile a code block
profile_code <- function(step_name, expr) {
  .start_profiling_step(step_name)
  on.exit(.end_profiling_step(step_name))
  eval.parent(substitute(expr))
}
