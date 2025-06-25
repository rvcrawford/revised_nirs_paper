# test_streamlined.R - Test the streamlined pipeline

cat("=== TESTING STREAMLINED PIPELINE ===\n")

# Load required libraries
library(targets)

# Check if core files exist
files_to_check <- c("R/config.R", "R/core_functions.R", "_targets.R")
missing_files <- character()

for (file in files_to_check) {
  if (!file.exists(file)) {
    missing_files <- c(missing_files, file)
  }
}

if (length(missing_files) > 0) {
  cat("❌ Missing files:", paste(missing_files, collapse = ", "), "\n")
  cat("Please create these files before testing\n")
  quit()
}

# Source config and core functions
tryCatch({
  source("R/config.R")
  source("R/core_functions.R")
  cat("✅ Successfully loaded core functions\n")
}, error = function(e) {
  cat("❌ Error loading functions:", e$message, "\n")
  quit()
})

# Test configuration
cat("\nTesting configuration system...\n")
dev_mode()
config <- get_analysis_config()
cat("✅ Configuration system working\n")

# Check data files
if (file.exists("./input_data/final_data_set/full_hemp_data.csv")) {
  cat("✅ Hemp data file found\n")
} else {
  cat("⚠️ Hemp data file not found - check path\n")
}

# Test targets pipeline structure
cat("\nTesting targets pipeline...\n")
tryCatch({
  tar_manifest()
  cat("✅ Targets pipeline structure valid\n")
}, error = function(e) {
  cat("❌ Targets pipeline error:", e$message, "\n")
})

cat("\n=== STREAMLINED PIPELINE TEST COMPLETE ===\n")
cat("If all checks passed, you can run: tar_make()\n")

