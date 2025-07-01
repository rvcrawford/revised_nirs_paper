# R/config.R - Simple Configuration for Hemp NIR Analysis

# =============================================================================
# ANALYSIS CONFIGURATION
# =============================================================================

# Global configuration settings
.hemp_config <- list(
  # Development vs production mode
  mode = "development",
  
  # Analysis parameters
  n_iterations = 100,           # Preprocessing comparison iterations
  final_iterations = 1000,     # Final modeling iterations
  test_size = 0.25,            # Test set proportion
  max_components = 20,         # Maximum PLS components
  cv_folds = 10,               # Cross-validation folds
  
  # Random seed for reproducibility
  seed = 123,
  
  # File paths
  input_path = "./input_data/final_data_set/",
  output_path = "./output/",
  
  # Preprocessing methods to test
  preprocessing_methods = 1:8,
  
  # Performance thresholds
  good_rpd = 2.0,
  excellent_rpd = 2.5,
  good_rsq = 0.7,
  excellent_rsq = 0.8
)

# =============================================================================
# CONFIGURATION FUNCTIONS
# =============================================================================

get_analysis_config <- function() {
  # Get current analysis configuration
  return(.hemp_config)
}

set_analysis_config <- function(key, value) {
  # Update a configuration value
  .hemp_config[[key]] <<- value
  cat("Config updated:", key, "=", value, "\n")
}

dev_mode <- function() {
  # Switch to development mode (faster, fewer iterations)
  .hemp_config$mode <<- "development"
  .hemp_config$n_iterations <<- 50
  .hemp_config$final_iterations <<- 100
  cat("Development mode: Fast testing with reduced iterations\n")
}

production_mode <- function() {
  # Switch to production mode (full analysis)
  .hemp_config$mode <<- "production"
  .hemp_config$n_iterations <<- 100
  .hemp_config$final_iterations <<- 1000
  cat("Production mode: Full analysis with complete iterations\n")
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

resolve_param <- function(user_value, config_value, param_name) {
  # Resolve parameter: user input overrides config
  if (!is.null(user_value)) {
    return(user_value)
  } else {
    return(config_value)
  }
}

# Initialize in development mode by default
dev_mode()