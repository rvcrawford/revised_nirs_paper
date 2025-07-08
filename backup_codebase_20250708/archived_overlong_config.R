# R/config.R - Hemp NIR Analysis Configuration
# ==============================================================================
# 
# Global configuration settings for the hemp NIR spectroscopy analysis pipeline.
# Modify these parameters to customize analysis behavior without changing core code.
# 
# ==============================================================================

# ==============================================================================
# ANALYSIS PARAMETERS
# ==============================================================================

#' Default analysis configuration
#' 
#' @return List containing all analysis parameters
#' @export
get_default_config <- function() {
  list(
    # Iteration parameters for robustness assessment
    preprocessing_iterations = 100,
    final_modeling_iterations = 1000,
    algorithm_comparison_iterations = 100,
    
    # Model optimization parameters
    max_pls_components = 20,
    cv_folds = 10,
    train_test_split = 0.75,
    
    # Data processing parameters
    protein_conversion_factor = 10,  # Convert % to g/kg
    
    # Performance thresholds (RPD-based classification)
    performance_thresholds = list(
      excellent_rpd = 3.0,
      good_rpd_min = 2.0,
      good_rpd_max = 3.0,
      fair_rpd_min = 1.4,
      fair_rpd_max = 2.0,
      
      # Supporting metric thresholds
      excellent_rpiq = 4.1,
      excellent_rsq = 0.9,
      good_rpiq_min = 2.3,
      good_rpiq_max = 4.1,
      good_rsq_min = 0.8,
      good_rsq_max = 0.9,
      fair_rpiq_min = 1.5,
      fair_rpiq_max = 2.3,
      fair_rsq_min = 0.5,
      fair_rsq_max = 0.8
    ),
    
    # Spectral analysis parameters
    vip_threshold = 1.0,
    protein_wavelength_bands = list(
      list(start = 1180, end = 1230, name = "C-H stretch 2nd overtone (amino acids)"),
      list(start = 1480, end = 1530, name = "N-H stretch 1st overtone (peptide bonds)"),
      list(start = 2040, end = 2070, name = "N-H + C-N combination (protein backbone)")
    ),
    
    # Machine learning algorithms for comparison
    comparison_algorithms = c("pls", "svmRadial", "rf"),
    
    # Preprocessing methods to evaluate (1-8)
    preprocessing_methods = 1:8,
    
    # Reproducibility
    random_seeds = list(
      global = 123,
      preprocessing = 456,
      modeling = 789,
      algorithms = 101112
    ),
    
    # Output formatting
    decimal_places = list(
      rmse = 2,
      rsq = 3,
      rpd = 2,
      rpiq = 2,
      percentages = 1
    )
  )
}

# ==============================================================================
# FILE PATHS AND DIRECTORIES
# ==============================================================================

#' Get standardized file paths
#' 
#' @param base_path Base directory path (default: current working directory)
#' @return List of standardized paths
#' @export
get_file_paths <- function(base_path = ".") {
  list(
    # Input data
    input_dir = file.path(base_path, "input_data", "final_data_set"),
    hemp_data = file.path(base_path, "input_data", "final_data_set", "full_hemp_data.csv"),
    preprocessing_key = file.path(base_path, "input_data", "final_data_set", "preprocessing_key.csv"),
    cultivar_table = file.path(base_path, "input_data", "final_data_set", "cultivar_table_clean.csv"),
    
    # Output directories
    output_dir = file.path(base_path, "output"),
    figures_dir = file.path(base_path, "output", "figures"),
    tables_dir = file.path(base_path, "output", "tables"),
    models_dir = file.path(base_path, "output", "models"),
    reports_dir = file.path(base_path, "output", "reports"),
    
    # Manuscript
    manuscript_dir = file.path(base_path, "manuscript"),
    manuscript_file = file.path(base_path, "manuscript", "hemp_nir_paper.qmd"),
    
    # R scripts
    r_dir = file.path(base_path, "R"),
    functions_file = file.path(base_path, "R", "core_functions.R"),
    config_file = file.path(base_path, "R", "config.R")
  )
}

# ==============================================================================
# LOGGING CONFIGURATION
# ==============================================================================

#' Configure logging settings
#' 
#' @param level Minimum log level ("DEBUG", "INFO", "WARN", "ERROR")
#' @param output_file Optional file for log output
#' @param timestamp_format Format string for timestamps
#' @return Logging configuration list
#' @export
get_logging_config <- function(level = "INFO", 
                               output_file = NULL,
                               timestamp_format = "%Y-%m-%d %H:%M:%S") {
  list(
    level = level,
    output_file = output_file,
    timestamp_format = timestamp_format,
    console_output = TRUE,
    file_output = !is.null(output_file)
  )
}

# ==============================================================================
# COMPUTATIONAL RESOURCES
# ==============================================================================

#' Configure computational resources
#' 
#' @param use_parallel Whether to use parallel processing
#' @param n_cores Number of cores to use (NULL for auto-detect)
#' @param memory_limit Memory limit in GB
#' @return Resource configuration list
#' @export
get_resource_config <- function(use_parallel = FALSE,
                                n_cores = NULL,
                                memory_limit = 8) {
  if (use_parallel && is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  list(
    use_parallel = use_parallel,
    n_cores = n_cores,
    memory_limit_gb = memory_limit
  )
}

# ==============================================================================
# VISUALIZATION SETTINGS
# ==============================================================================

#' Configure visualization parameters
#' 
#' @return Visualization configuration list
#' @export
get_visualization_config <- function() {
  list(
    # Figure dimensions (inches)
    figure_width = 8,
    figure_height = 6,
    figure_dpi = 300,
    
    # Color schemes
    colors = list(
      primary = "#2C3E50",
      secondary = "#3498DB",
      accent = "#E74C3C",
      success = "#27AE60",
      warning = "#F39C12",
      neutral = "#95A5A6"
    ),
    
    # Theme settings
    base_theme = "minimal",
    font_family = "Arial",
    base_font_size = 12,
    title_font_size = 14,
    
    # Plot-specific settings
    point_alpha = 0.7,
    line_alpha = 0.8,
    box_alpha = 0.7,
    
    # File formats for saving
    output_formats = c("png", "pdf", "svg")
  )
}

# ==============================================================================
# VALIDATION RULES
# ==============================================================================

#' Define data validation rules
#' 
#' @return List of validation criteria
#' @export
get_validation_rules <- function() {
  list(
    # Minimum sample requirements
    min_samples_total = 50,
    min_samples_per_location = 10,
    min_samples_training = 30,
    
    # Data quality thresholds
    max_missing_proportion = 0.05,
    min_protein_range = 50,  # g/kg
    max_protein_cv = 0.5,    # coefficient of variation
    
    # Spectral data requirements
    required_wavelength_range = c(1100, 2498),
    min_wavelengths = 500,
    
    # Model performance minimums
    min_acceptable_rsq = 0.5,
    min_acceptable_rpd = 1.4,
    max_acceptable_rmse = 50,  # g/kg
    
    # Cross-validation requirements
    min_cv_folds = 5,
    min_iterations = 10
  )
}

# ==============================================================================
# ENVIRONMENT SETUP
# ==============================================================================

#' Initialize analysis environment
#' 
#' Sets up directories, checks dependencies, and validates configuration.
#' 
#' @param config Analysis configuration list
#' @param create_dirs Whether to create missing directories
#' @return Boolean indicating successful setup
#' @export
initialize_environment <- function(config = get_default_config(), create_dirs = TRUE) {
  
  # Get file paths
  paths <- get_file_paths()
  
  # Create output directories if requested
  if (create_dirs) {
    dirs_to_create <- c(
      paths$output_dir, paths$figures_dir, paths$tables_dir,
      paths$models_dir, paths$reports_dir
    )
    
    for (dir in dirs_to_create) {
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
        message("Created directory: ", dir)
      }
    }
  }
  
  # Check for required input files
  required_files <- c(paths$hemp_data, paths$preprocessing_key)
  missing_files <- required_files[!file.exists(required_files)]
  
  if (length(missing_files) > 0) {
    warning("Missing required input files: ", paste(missing_files, collapse = ", "))
    return(FALSE)
  }
  
  # Check R package dependencies
  required_packages <- c(
    "data.table", "caret", "prospectr", "pls", 
    "tidyverse", "kableExtra", "ggplot2", "targets"
  )
  
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    warning("Missing required packages: ", paste(missing_packages, collapse = ", "))
    message("Install with: install.packages(c(", paste0("'", missing_packages, "'", collapse = ", "), "))")
    return(FALSE)
  }
  
  message("Environment initialization completed successfully")
  return(TRUE)
}

# ==============================================================================
# CONFIGURATION LOADING
# ==============================================================================

# Load default configuration on script sourcing
if (!exists("HEMP_NIR_CONFIG")) {
  HEMP_NIR_CONFIG <- get_default_config()
  message("Loaded default hemp NIR analysis configuration")
}