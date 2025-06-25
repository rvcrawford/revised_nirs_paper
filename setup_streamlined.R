# setup_streamlined.R - Implement Streamlined Hemp NIR Analysis
# Run this script to set up the streamlined version

cat("=== SETTING UP STREAMLINED HEMP NIR ANALYSIS ===\n\n")

# Step 1: Create backup of current code
cat("Step 1: Creating backup of current code...\n")
if (!dir.exists("backup_original_code")) {
  dir.create("backup_original_code")
  
  # Backup R files
  if (dir.exists("R")) {
    file.copy("R", "backup_original_code", recursive = TRUE)
    cat("✅ R/ directory backed up\n")
  }
  
  # Backup _targets.R
  if (file.exists("_targets.R")) {
    file.copy("_targets.R", "backup_original_code/_targets.R")
    cat("✅ _targets.R backed up\n")
  }
  
  # Backup any .qmd files
  qmd_files <- list.files(".", pattern = "\\.qmd$", full.names = TRUE)
  if (length(qmd_files) > 0) {
    file.copy(qmd_files, "backup_original_code/")
    cat("✅", length(qmd_files), ".qmd files backed up\n")
  }
  
} else {
  cat("⚠️ Backup directory already exists - skipping backup\n")
}

# Step 2: Create archive for non-essential code
cat("\nStep 2: Creating archive for non-essential code...\n")
if (!dir.exists("code_archive")) {
  dir.create("code_archive")
  
  # Files to archive (not essential for core paper)
  files_to_archive <- c(
    "R/multi_algorithm_functions.R",
    "R/spectral_analysis_functions.R", 
    "R/weighted_sampling_functions.R",
    "R/weighted_sampling_comparison.R",
    "R/location_weighting_functions.R",
    "R/comprehensive_geneva_findings.R",
    "minimal_test.R",
    "R/test_geneva_analysis.R"
  )
  
  archived_count <- 0
  for (file in files_to_archive) {
    if (file.exists(file)) {
      file.copy(file, file.path("code_archive", basename(file)))
      cat("📦 Archived:", basename(file), "\n")
      archived_count <- archived_count + 1
    }
  }
  
  cat("✅", archived_count, "files archived for future use\n")
  
} else {
  cat("⚠️ Archive directory already exists - skipping archive creation\n")
}

# Step 3: Create manuscript directory
cat("\nStep 3: Setting up manuscript directory...\n")
if (!dir.exists("manuscript")) {
  dir.create("manuscript")
  cat("✅ Created manuscript/ directory\n")
} else {
  cat("ℹ️ manuscript/ directory already exists\n")
}

# Step 4: Check required packages
cat("\nStep 4: Checking required packages...\n")
required_packages <- c(
  "targets", "tarchetypes", "tidyverse", "data.table", 
  "caret", "prospectr", "pls", "kableExtra", "ggplot2"
)

missing_packages <- character()
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) > 0) {
  cat("⚠️ Missing packages detected:", paste(missing_packages, collapse = ", "), "\n")
  cat("Installing missing packages...\n")
  install.packages(missing_packages)
} else {
  cat("✅ All required packages are installed\n")
}

# Step 5: Create streamlined files message
cat("\nStep 5: Ready to create streamlined files\n")
cat("You now need to:\n")
cat("1. Copy the _targets.R content from the artifact into your _targets.R file\n")
cat("2. Create R/core_functions.R with the content from the artifact\n") 
cat("3. Create manuscript/hemp_nir_paper.qmd with the content from the artifact\n")
cat("4. Run the test script below\n\n")

# Create a simple test script
cat("Creating test script...\n")
test_script <- '# test_streamlined.R - Test the streamlined pipeline

cat("=== TESTING STREAMLINED PIPELINE ===\\n")

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
  cat("❌ Missing files:", paste(missing_files, collapse = ", "), "\\n")
  cat("Please create these files before testing\\n")
  quit()
}

# Source config and core functions
tryCatch({
  source("R/config.R")
  source("R/core_functions.R")
  cat("✅ Successfully loaded core functions\\n")
}, error = function(e) {
  cat("❌ Error loading functions:", e$message, "\\n")
  quit()
})

# Test configuration
cat("\\nTesting configuration system...\\n")
dev_mode()
config <- get_analysis_config()
cat("✅ Configuration system working\\n")

# Check data files
if (file.exists("./input_data/final_data_set/full_hemp_data.csv")) {
  cat("✅ Hemp data file found\\n")
} else {
  cat("⚠️ Hemp data file not found - check path\\n")
}

# Test targets pipeline structure
cat("\\nTesting targets pipeline...\\n")
tryCatch({
  tar_manifest()
  cat("✅ Targets pipeline structure valid\\n")
}, error = function(e) {
  cat("❌ Targets pipeline error:", e$message, "\\n")
})

cat("\\n=== STREAMLINED PIPELINE TEST COMPLETE ===\\n")
cat("If all checks passed, you can run: tar_make()\\n")
'

writeLines(test_script, "test_streamlined.R")
cat("✅ Created test_streamlined.R\n")

# Create quick verification script for later
verification_script <- '# verify_streamlined.R - Verify streamlined results

library(targets)

cat("=== VERIFYING STREAMLINED RESULTS ===\\n")

# Check targets status
cat("\\nTargets status:\\n")
print(tar_progress())

# Load and summarize key results
if (tar_exist_objects(c("preprocessing_analysis", "final_model_analysis"))) {
  
  # Preprocessing results
  prep_analysis <- tar_read(preprocessing_analysis)
  cat("\\n📊 PREPROCESSING COMPARISON:\\n")
  cat("- Methods tested:", nrow(prep_analysis$summary), "\\n")
  cat("- Best method:", tar_read(best_method), "\\n")
  
  # Final model results  
  final_analysis <- tar_read(final_model_analysis)
  cat("\\n📊 FINAL MODEL PERFORMANCE:\\n")
  cat("- Mean RMSE:", round(final_analysis$overall_stats$mean_rmse, 2), "g/kg\\n")
  cat("- Mean R²:", round(final_analysis$overall_stats$mean_rsq, 3), "\\n")
  cat("- Mean RPD:", round(final_analysis$overall_stats$mean_rpd, 2), "\\n")
  cat("- Excellent models:", final_analysis$classification$excellent, "\\n")
  
  # Tables and figures
  cat("\\n📊 OUTPUTS GENERATED:\\n")
  outputs <- c("table_sample_summary", "table_preprocessing", 
               "fig_model_calibration", "fig_performance_boxplot", "fig_validation_errors")
  
  for (output in outputs) {
    if (tar_exist_objects(output)) {
      cat("✅", output, "\\n")
    } else {
      cat("❌", output, "\\n")
    }
  }
  
} else {
  cat("⚠️ Key results not found - run tar_make() first\\n")
}

cat("\\n=== VERIFICATION COMPLETE ===\\n")
'

writeLines(verification_script, "verify_streamlined.R")
cat("✅ Created verify_streamlined.R\n")

cat("\n=== SETUP COMPLETE ===\n")
cat("Next steps:\n")
cat("1. Manually create the 3 main files from the artifacts\n")
cat("2. Run: source('test_streamlined.R')\n") 
cat("3. If tests pass, run: tar_make()\n")
cat("4. Run: source('verify_streamlined.R') to check results\n")
cat("\nYour original code is safely backed up in backup_original_code/\n")
cat("Non-essential code is archived in code_archive/ for future papers\n")