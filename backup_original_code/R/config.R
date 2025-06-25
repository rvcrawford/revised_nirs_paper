# Step 1: Create R/config.R file
# Copy this entire code block into a new file: R/config.R

# Configuration function for development vs production modes
get_analysis_config <- function(mode = NULL) {
  # Get mode from environment if not specified
  if (is.null(mode)) {
    mode <- Sys.getenv("ANALYSIS_MODE", "development")
  }
  
  configs <- list(
    development = list(
      n_iterations = 10,
      n_components_max = 10,
      parallel = FALSE
    ),
    testing = list(
      n_iterations = 50,
      n_components_max = 15,
      parallel = TRUE
    ),
    production = list(
      n_iterations = 1000,
      n_components_max = 20,
      parallel = TRUE
    )
  )
  
  return(configs[[mode]])
}

# Helper function to resolve parameters with precedence
resolve_param <- function(explicit_value, config_value, param_name) {
  if (!is.null(explicit_value)) {
    # User explicitly set this parameter - use it
    cat("Using explicit", param_name, "=", explicit_value, "\n")
    return(explicit_value)
  } else {
    # Use config default
    cat("Using config default", param_name, "=", config_value, 
        "(mode:", Sys.getenv("ANALYSIS_MODE", "development"), ")\n")
    return(config_value)
  }
}

# Convenience functions for switching modes
dev_mode <- function() {
  Sys.setenv(ANALYSIS_MODE = "development")
  cat("✅ Switched to development mode (fast iterations)\n")
  cat("Current config:\n")
  print(get_analysis_config())
}

prod_mode <- function() {
  Sys.setenv(ANALYSIS_MODE = "production")
  cat("✅ Switched to production mode (full iterations)\n")
  cat("Current config:\n")
  print(get_analysis_config())
}

test_mode <- function() {
  Sys.setenv(ANALYSIS_MODE = "testing")
  cat("✅ Switched to testing mode (medium iterations)\n")
  cat("Current config:\n")
  print(get_analysis_config())
}

# Check current mode
check_mode <- function() {
  current_mode <- Sys.getenv("ANALYSIS_MODE", "development")
  config <- get_analysis_config()
  
  cat("Current analysis mode:", current_mode, "\n")
  cat("Configuration:\n")
  print(config)
  
  return(list(mode = current_mode, config = config))
}