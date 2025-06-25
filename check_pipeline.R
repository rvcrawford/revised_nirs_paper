# Hemp NIR Pipeline Status Checker
# Save this as check_pipeline.R in your project root directory

library(targets)

cat("=== HEMP NIR PIPELINE STATUS ===\n\n")

# Check if _targets directory exists
if (!file.exists("_targets")) {
  cat("❌ No _targets directory found\n")
  cat("   Run: tar_make()\n\n")
  stop("Pipeline not initialized")
}

# 1. Check essential targets
essential_targets <- c(
  "hemp_data",
  "best_method", 
  "final_model_analysis",
  "preprocessing_analysis"
)

cat("1. ESSENTIAL TARGETS:\n")
for (target in essential_targets) {
  exists <- tar_exist_objects(target)
  status <- if (exists) "✅" else "❌"
  cat(sprintf("   %s %s\n", status, target))
}

# 2. Check figure targets
figure_targets <- c(
  "fig_model_calibration",
  "fig_performance_boxplot", 
  "fig_validation_errors"
)

cat("\n2. FIGURE TARGETS:\n")
for (target in figure_targets) {
  exists <- tar_exist_objects(target)
  status <- if (exists) "✅" else "❌"
  cat(sprintf("   %s %s\n", status, target))
}

# 3. Check for errors
cat("\n3. TARGET ERRORS:\n")
meta <- tar_meta()
error_targets <- meta[meta$error != "character(0)", c("name", "error")]

if (nrow(error_targets) > 0) {
  cat("❌ Targets with errors:\n")
  for (i in 1:nrow(error_targets)) {
    cat(sprintf("   - %s: %s\n", error_targets$name[i], error_targets$error[i]))
  }
} else {
  cat("✅ No targets have errors\n")
}

# 4. Check recent builds
cat("\n4. BUILD STATUS:\n")
recent_builds <- meta[order(meta$time, decreasing = TRUE), c("name", "time")]
if (nrow(recent_builds) > 0) {
  cat("Most recent builds:\n")
  head(recent_builds, 5) %>% 
    apply(1, function(x) cat(sprintf("   %s: %s\n", x[1], x[2])))
} else {
  cat("❌ No successful builds found\n")
}

# 5. Recommendations
cat("\n5. RECOMMENDATIONS:\n")

missing_essential <- sum(!sapply(essential_targets, tar_exist_objects))
missing_figures <- sum(!sapply(figure_targets, tar_exist_objects))

if (missing_essential > 0) {
  cat("🔧 CRITICAL: Essential targets missing\n")
  cat("   Run: tar_make() to build the pipeline\n")
} else if (missing_figures > 0) {
  cat("🔧 Figures missing but data available\n") 
  cat("   The figures may be failing to generate\n")
  cat("   Check figure creation functions in R/core_functions.R\n")
} else {
  cat("✅ All targets available\n")
  cat("   Figures should display properly\n")
  cat("   Re-render the Quarto document\n")
}

cat("\n=== STATUS CHECK COMPLETE ===\n")