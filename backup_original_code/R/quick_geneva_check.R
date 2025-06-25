# FIXED VERSION OF QUICK GENEVA CHECK
# Replace the print statement in your quick_geneva_check function

quick_geneva_check <- function() {
  library(targets)
  library(dplyr)
  
  cat("=== QUICK GENEVA CHECK ===\n\n")
  
  # Load data
  wc <- tar_read(weighting_comparison)
  bd <- tar_read(balanced_data)
  
  # Calculate location-specific RMSE
  weighted_by_loc <- wc$weighted$predictions %>%
    group_by(location) %>%
    summarise(
      weighted_rmse = sqrt(mean((actual - predicted)^2, na.rm = TRUE)),
      weighted_n = n(),
      .groups = 'drop'
    )
  
  unweighted_by_loc <- wc$unweighted$predictions %>%
    group_by(location) %>%
    summarise(
      unweighted_rmse = sqrt(mean((actual - predicted)^2, na.rm = TRUE)),
      unweighted_n = n(),
      .groups = 'drop'
    )
  
  # Combine results
  comparison <- weighted_by_loc %>%
    full_join(unweighted_by_loc, by = "location") %>%
    mutate(
      rmse_improvement = unweighted_rmse - weighted_rmse,
      pct_improvement = 100 * rmse_improvement / unweighted_rmse
    )
  
  # FIXED: Join with balanced data safely
  tryCatch({
    comparison <- comparison %>%
      left_join(bd$location_summary, by = "location")
  }, error = function(e) {
    cat("Warning: Could not join with location summary:", e$message, "\n")
  })
  
  cat("LOCATION-SPECIFIC RMSE RESULTS:\n")
  cat("================================\n")
  
  # FIXED: Only select columns that actually exist
  available_cols <- intersect(
    c("location", "weight", "weighted_rmse", "unweighted_rmse", 
      "rmse_improvement", "pct_improvement"),
    names(comparison)
  )
  
  if (length(available_cols) > 1) {
    result_table <- comparison %>% 
      select(all_of(available_cols))
    
    if ("weight" %in% available_cols) {
      result_table <- result_table %>% arrange(desc(weight))
    }
    
    print(result_table)
  } else {
    print(comparison)
  }
  
  # Focus on Geneva
  geneva_result <- comparison %>% filter(location == "geneva")
  
  cat("\n\nGENEVA SPECIFIC RESULTS:\n")
  cat("========================\n")
  
  if (nrow(geneva_result) > 0) {
    if ("weight" %in% names(geneva_result)) {
      cat("Sample weight:", geneva_result$weight, "x\n")
    }
    cat("Weighted RMSE:", round(geneva_result$weighted_rmse, 3), "\n")
    cat("Unweighted RMSE:", round(geneva_result$unweighted_rmse, 3), "\n")
    cat("Improvement:", round(geneva_result$rmse_improvement, 3), "\n")
    cat("% Improvement:", round(geneva_result$pct_improvement, 1), "%\n")
    
    geneva_improves <- geneva_result$rmse_improvement > 0
    cat("\nDoes weighting help Geneva?", geneva_improves, "\n")
    
    if (geneva_improves) {
      cat("✓ YES - Geneva RMSE improves by", round(geneva_result$rmse_improvement, 3), "units\n")
      
      # Check if this comes at cost of others
      other_locations <- comparison %>% filter(location != "geneva")
      if (nrow(other_locations) > 0) {
        avg_other_change <- mean(other_locations$rmse_improvement, na.rm = TRUE)
        cat("Average change for other locations:", round(avg_other_change, 3), "\n")
        
        if (avg_other_change < 0) {
          cat("⚠ This improvement comes at the cost of other locations\n")
          cat("  This explains why overall performance doesn't improve\n")
        }
      }
    } else {
      cat("✗ NO - Geneva doesn't benefit from weighting\n")
    }
  } else {
    cat("ERROR: No Geneva results found!\n")
  }
  
  return(comparison)
}