# Profiling functions for AlienDetective scripts

# Global variable for profiling data
profiling_data <- NULL

# Initialize profiling environment
.init_profiling <- function() {
  # Load required packages
  if (!requireNamespace("pryr", quietly = TRUE)) {
    install.packages("pryr")
  }
  if (!requireNamespace("tictoc", quietly = TRUE)) {
    install.packages("tictoc")
  }
  if (!requireNamespace("knitr", quietly = TRUE)) {
    install.packages("knitr")
  }
  if (!requireNamespace("kableExtra", quietly = TRUE)) {
    install.packages("kableExtra")
  }
  
  library(pryr)
  library(tictoc)
  
  # Initialize profiling data frame
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
  
  # Start overall timer
  tic("Total script execution")
}

# Start profiling a step
.start_profiling_step <- function(step_name) {
  if (!exists("profiling_data", envir = .GlobalEnv)) {
    .init_profiling()
  }
  
  # Start timer for this step
  tic(step_name)
  
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
  mem_used_human <- format(structure(mem_used_bytes, class = "object_size"), 
                          units = "auto")
  
  # Get timing
  timing <- toc(quiet = TRUE)
  
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

# Generate HTML report
generate_profiling_report <- function(output_dir = "profiling_reports") {
  if (!exists("profiling_data", envir = .GlobalEnv) || nrow(profiling_data) == 0) {
    warning("No profiling data available to generate report")
    return(invisible(NULL))
  }
  
  # Stop overall timer and calculate total time
  if (exists("global_tic", envir = .GlobalEnv)) {
    total_time <- toc(quiet = TRUE)
    total_seconds <- round(total_time$toc - total_time$tic, 2)
  } else {
    # Fallback if global timer wasn't started
    total_seconds <- sum(profiling_data$Duration_sec, na.rm = TRUE)
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Generate timestamp for filename
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  report_path <- file.path(output_dir, paste0("profiling_report_", timestamp, ".html"))
  
  # Create HTML report
  html_content <- paste0(
    '<!DOCTYPE html>
    <html>
    <head>
      <title>Profiling Report - ', timestamp, '</title>
      <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #2c3e50; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        tr:nth-child(even) { background-color: #f9f9f9; }
        .summary { margin: 20px 0; padding: 15px; background-color: #f0f7fb; border-radius: 5px; }
      </style>
    </head>
    <body>
      <h1>Profiling Report</h1>
      <div class="summary">
        <p><strong>Script:</strong> ', basename(sys.frame(1)$ofile), '</p>
        <p><strong>Total execution time:</strong> ', 
        total_seconds, ' seconds</p>
        <p><strong>Number of steps profiled:</strong> ', nrow(profiling_data), '</p>
        <p><strong>Peak memory used:</strong> ', 
        format(structure(max(profiling_data$Memory_Bytes), class = "object_size"), 
               units = "auto"), '</p>
      </div>
      ', knitr::kable(profiling_data, format = "html", escape = FALSE) %>%
         kableExtra::kable_styling("striped", full_width = TRUE) %>%
         kableExtra::column_spec(1, bold = TRUE) %>%
         kableExtra::row_spec(0, bold = TRUE, color = "white", background = "#2c3e50") %>%
         kableExtra::row_spec(which.max(profiling_data$Memory_Bytes), 
                             background = "#fff3cd") %>%
         kableExtra::row_spec(which.max(profiling_data$Duration_sec), 
                             background = "#f8d7da"),
      '
      <script>
        // Add some basic interactivity
        document.addEventListener("DOMContentLoaded", function() {
          const rows = document.querySelectorAll("tr");
          rows.forEach(row => {
            row.addEventListener("click", function() {
              this.classList.toggle("highlight");
            });
          });
        });
      </script>
      <style>
        tr.highlight { background-color: #e8f4f8 !important; }
      </style>
    </body>
    </html>'
  )
  
  # Write to file
  writeLines(html_content, report_path)
  message("Profiling report generated at: ", normalizePath(report_path))
  return(invisible(report_path))
}

# Helper function to profile a code block
profile_code <- function(step_name, expr) {
  .start_profiling_step(step_name)
  on.exit(.end_profiling_step(step_name))
  eval.parent(substitute(expr))
}
