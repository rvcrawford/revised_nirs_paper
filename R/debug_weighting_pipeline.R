# Step 1: Debug the weighting comparison pipeline

# Check what objects exist and their structure
debug_weighting_pipeline <- function() {
  
  cat("=== DEBUGGING WEIGHTING PIPELINE ===\n")
  
  # Check if targets exist
  library(targets)
  
  targets_to_check <- c(
    "balanced_data",
    "weighting_comparison", 
    "weighting_analysis",
    "table_weighting_comparison"
  )
  
  for (target in targets_to_check) {
    tryCatch({
      obj <- tar_read(!!sym(target))
      cat("\n", target, ":\n")
      cat("  Class:", class(obj), "\n")
      cat("  Names:", names(obj), "\n")
      if (is.list(obj)) {
        cat("  Length:", length(obj), "\n")
      }
    }, error = function(e) {
      cat("\n", target, ": ERROR -", e$message, "\n")
    })
  }
  
  # Specifically check weighting_comparison structure
  cat("\n=== DETAILED WEIGHTING_COMPARISON CHECK ===\n")
  tryCatch({
    wc <- tar_read(weighting_comparison)
    
    cat("Weighted results structure:\n")
    if ("weighted" %in% names(wc)) {
      cat("  Names:", names(wc$weighted), "\n")
      if ("metrics" %in% names(wc$weighted)) {
        cat("  Metrics dims:", dim(wc$weighted$metrics), "\n")
        cat("  Metrics cols:", names(wc$weighted$metrics), "\n")
        cat("  First few rows:\n")
        print(head(wc$weighted$metrics, 3))
      }
    }
    
    cat("\nUnweighted results structure:\n")
    if ("unweighted" %in% names(wc)) {
      cat("  Names:", names(wc$unweighted), "\n")
      if ("metrics" %in% names(wc$unweighted)) {
        cat("  Metrics dims:", dim(wc$unweighted$metrics), "\n")
        cat("  Metrics cols:", names(wc$unweighted$metrics), "\n")
      }
    }
    
  }, error = function(e) {
    cat("ERROR reading weighting_comparison:", e$message, "\n")
  })
}

# Run the debug function
debug_weighting_pipeline()