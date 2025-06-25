# verify_streamlined.R - Verify streamlined results

library(targets)

cat("=== VERIFYING STREAMLINED RESULTS ===\n")

# Check targets status
cat("\nTargets status:\n")
print(tar_progress())

# Load and summarize key results
if (tar_exist_objects(c("preprocessing_analysis", "final_model_analysis"))) {
  
  # Preprocessing results
  prep_analysis <- tar_read(preprocessing_analysis)
  cat("\n📊 PREPROCESSING COMPARISON:\n")
  cat("- Methods tested:", nrow(prep_analysis$summary), "\n")
  cat("- Best method:", tar_read(best_method), "\n")
  
  # Final model results  
  final_analysis <- tar_read(final_model_analysis)
  cat("\n📊 FINAL MODEL PERFORMANCE:\n")
  cat("- Mean RMSE:", round(final_analysis$overall_stats$mean_rmse, 2), "g/kg\n")
  cat("- Mean R²:", round(final_analysis$overall_stats$mean_rsq, 3), "\n")
  cat("- Mean RPD:", round(final_analysis$overall_stats$mean_rpd, 2), "\n")
  cat("- Excellent models:", final_analysis$classification$excellent, "\n")
  
  # Tables and figures
  cat("\n📊 OUTPUTS GENERATED:\n")
  outputs <- c("table_sample_summary", "table_preprocessing", 
               "fig_model_calibration", "fig_performance_boxplot", "fig_validation_errors")
  
  for (output in outputs) {
    if (tar_exist_objects(output)) {
      cat("✅", output, "\n")
    } else {
      cat("❌", output, "\n")
    }
  }
  
} else {
  cat("⚠️ Key results not found - run tar_make() first\n")
}

cat("\n=== VERIFICATION COMPLETE ===\n")

