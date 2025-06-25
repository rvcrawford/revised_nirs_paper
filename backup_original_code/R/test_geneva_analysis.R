# Immediate Geneva analysis test using existing targets data

test_geneva_with_existing_data <- function() {
  cat("=== TESTING GENEVA ANALYSIS WITH EXISTING DATA ===\n\n")
  
  library(targets)
  library(dplyr)
  
  # Try to load existing data
  tryCatch({
    wc <- tar_read(weighting_comparison)
    bd <- tar_read(balanced_data)
    
    cat("Loaded existing weighting comparison and balanced data\n")
    
    # Check what's available in predictions
    cat("\nWeighted predictions structure:\n")
    if ("predictions" %in% names(wc$weighted)) {
      cat("Columns:", names(wc$weighted$predictions), "\n")
      cat("First few rows:\n")
      print(head(wc$weighted$predictions, 3))
      
      has_location_weighted <- "location" %in% names(wc$weighted$predictions)
    } else {
      cat("No predictions found in weighted results\n")
      has_location_weighted <- FALSE
    }
    
    cat("\nUnweighted predictions structure:\n")
    if ("predictions" %in% names(wc$unweighted)) {
      cat("Columns:", names(wc$unweighted$predictions), "\n")
      has_location_unweighted <- "location" %in% names(wc$unweighted$predictions)
    } else {
      cat("No predictions found in unweighted results\n")
      has_location_unweighted <- FALSE
    }
    
    cat("\nLocation availability:\n")
    cat("- Weighted has location:", has_location_weighted, "\n")
    cat("- Unweighted has location:", has_location_unweighted, "\n")
    
    if (has_location_weighted && has_location_unweighted) {
      cat("\n✅ BOTH HAVE LOCATION - CAN ANALYZE GENEVA!\n")
      
      # Run Geneva analysis
      geneva_analysis <- analyze_geneva_performance_simple(wc, bd)
      
      return(geneva_analysis)
      
    } else if (has_location_weighted && !has_location_unweighted) {
      cat("\n⚠️ ONLY WEIGHTED HAS LOCATION\n")
      cat("This suggests the unweighted function doesn't preserve location.\n")
      cat("You'll need to re-run with the fixed functions.\n")
      
      # Can still analyze weighted vs expected for Geneva
      geneva_weighted_analysis <- analyze_geneva_weighted_only(wc, bd)
      return(geneva_weighted_analysis)
      
    } else {
      cat("\n❌ NEITHER HAS LOCATION\n")
      cat("Need to re-run weighting comparison with location-preserving functions.\n")
      return(NULL)
    }
    
  }, error = function(e) {
    cat("Error loading data:", e$message, "\n")
    return(NULL)
  })
}

# Simple Geneva analysis if both have location data
analyze_geneva_performance_simple <- function(weighting_comparison, balanced_data) {
  cat("\n=== GENEVA PERFORMANCE ANALYSIS ===\n")
  
  library(dplyr)
  
  # Calculate location-specific RMSE
  weighted_by_loc <- weighting_comparison$weighted$predictions %>%
    group_by(location) %>%
    summarise(
      weighted_rmse = sqrt(mean((actual - predicted)^2, na.rm = TRUE)),
      weighted_n = n(),
      .groups = 'drop'
    )
  
  unweighted_by_loc <- weighting_comparison$unweighted$predictions %>%
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
    ) %>%
    left_join(balanced_data$location_summary, by = "location")
  
  cat("\nLocation-specific results:\n")
  print(comparison %>% 
          select(location, weight, weighted_rmse, unweighted_rmse, 
                 rmse_improvement, pct_improvement))
  
  # Focus on Geneva
  geneva_result <- comparison %>% filter(location == "geneva")
  
  cat("\n=== GENEVA SPECIFIC RESULTS ===\n")
  if (nrow(geneva_result) > 0) {
    cat("Sample weight:", geneva_result$weight, "x\n")
    cat("Weighted RMSE:", round(geneva_result$weighted_rmse, 3), "\n")
    cat("Unweighted RMSE:", round(geneva_result$unweighted_rmse, 3), "\n")
    cat("Improvement:", round(geneva_result$rmse_improvement, 3), "\n")
    cat("% Improvement:", round(geneva_result$pct_improvement, 1), "%\n")
    
    geneva_improves <- geneva_result$rmse_improvement > 0
    cat("\nDoes weighting help Geneva?", geneva_improves, "\n")
    
    if (geneva_improves) {
      cat("✅ YES - Geneva RMSE improves by", round(geneva_result$rmse_improvement, 3), "units\n")
      
      # Statistical test
      geneva_weighted <- weighting_comparison$weighted$predictions %>% 
        filter(location == "geneva") %>%
        mutate(abs_error = abs(actual - predicted))
      
      geneva_unweighted <- weighting_comparison$unweighted$predictions %>% 
        filter(location == "geneva") %>%
        mutate(abs_error = abs(actual - predicted))
      
      if (nrow(geneva_weighted) > 10 && nrow(geneva_unweighted) > 10) {
        wilcox_test <- wilcox.test(geneva_weighted$abs_error, 
                                   geneva_unweighted$abs_error,
                                   alternative = "less")
        
        cat("Statistical significance (p-value):", round(wilcox_test$p.value, 4), "\n")
        cat("Statistically significant:", wilcox_test$p.value < 0.05, "\n")
        
        # Effect size
        mean_diff <- mean(geneva_unweighted$abs_error) - mean(geneva_weighted$abs_error)
        pooled_sd <- sqrt((var(geneva_weighted$abs_error) + var(geneva_unweighted$abs_error)) / 2)
        cohens_d <- mean_diff / pooled_sd
        
        cat("Effect size (Cohen's d):", round(cohens_d, 3), "\n")
        cat("Effect interpretation:", 
            if(abs(cohens_d) < 0.2) "negligible"
            else if(abs(cohens_d) < 0.5) "small"
            else if(abs(cohens_d) < 0.8) "medium"
            else "large", "\n")
        
        # Generate recommendation
        cat("\n=== RECOMMENDATION FOR PAPER ===\n")
        
        if (wilcox_test$p.value < 0.05 && abs(cohens_d) > 0.2) {
          cat("🎯 REPORT THIS AS A POSITIVE FINDING!\n")
          cat("Geneva shows statistically significant improvement with meaningful effect size.\n")
          cat("\nSuggested text:\n")
          cat('"Location-weighted modeling significantly improved predictions for Geneva samples\n')
          cat('(', round(geneva_result$pct_improvement, 1), '% RMSE reduction, p = ', round(wilcox_test$p.value, 3), ', Cohen\'s d = ', round(cohens_d, 2), '),\n')
          cat('demonstrating that sample weighting can effectively address location-specific bias."\n')
          
        } else if (geneva_result$pct_improvement > 2) {
          cat("⚠️ REPORT AS PROMISING BUT NEEDS MORE DATA\n")
          cat("Geneva shows improvement but not statistically significant.\n")
          cat("\nSuggested text:\n")
          cat('"Location-weighted modeling showed promise for Geneva samples\n')
          cat('(', round(geneva_result$pct_improvement, 1), '% RMSE improvement) but larger sample sizes\n')
          cat('may be needed to achieve statistical significance."\n')
          
        } else {
          cat("⚪ REPORT AS MINIMAL EFFECT\n")
          cat("Geneva improvement is small and not significant.\n")
        }
        
        return(list(
          comparison = comparison,
          geneva_result = geneva_result,
          statistical_test = wilcox_test,
          effect_size = cohens_d,
          significant = wilcox_test$p.value < 0.05,
          meaningful = abs(cohens_d) > 0.2
        ))
        
      } else {
        cat("Insufficient data for statistical test\n")
        return(list(comparison = comparison, geneva_result = geneva_result))
      }
      
    } else {
      cat("❌ NO - Geneva doesn't benefit from weighting\n")
      return(list(comparison = comparison, geneva_result = geneva_result))
    }
    
  } else {
    cat("No Geneva results found\n")
    return(NULL)
  }
}

# If only weighted has location, analyze weighted performance
analyze_geneva_weighted_only <- function(weighting_comparison, balanced_data) {
  cat("\n=== GENEVA WEIGHTED-ONLY ANALYSIS ===\n")
  cat("(Can only analyze weighted predictions since unweighted lacks location data)\n")
  
  # Analyze weighted predictions by location
  weighted_by_loc <- weighting_comparison$weighted$predictions %>%
    group_by(location) %>%
    summarise(
      weighted_rmse = sqrt(mean((actual - predicted)^2, na.rm = TRUE)),
      weighted_bias = mean(predicted - actual, na.rm = TRUE),
      n_predictions = n(),
      .groups = 'drop'
    ) %>%
    left_join(balanced_data$location_summary, by = "location")
  
  cat("\nWeighted performance by location:\n")
  print(weighted_by_loc)
  
  geneva_weighted <- weighted_by_loc %>% filter(location == "geneva")
  
  if (nrow(geneva_weighted) > 0) {
    cat("\nGeneva weighted performance:\n")
    cat("RMSE:", round(geneva_weighted$weighted_rmse, 3), "\n")
    cat("Bias:", round(geneva_weighted$weighted_bias, 3), "\n")
    cat("Sample weight:", geneva_weighted$weight, "x\n")
    
    # Compare to other locations
    other_locations <- weighted_by_loc %>% filter(location != "geneva")
    avg_other_rmse <- mean(other_locations$weighted_rmse, na.rm = TRUE)
    
    cat("\nComparison to other locations:\n")
    cat("Geneva RMSE:", round(geneva_weighted$weighted_rmse, 3), "\n")
    cat("Other locations average RMSE:", round(avg_other_rmse, 3), "\n")
    cat("Geneva relative performance:", 
        ifelse(geneva_weighted$weighted_rmse < avg_other_rmse, "Better", "Worse"), "\n")
  }
  
  return(weighted_by_loc)
}

# Run the test
cat("Testing Geneva analysis with existing data...\n")
geneva_test_results <- test_geneva_with_existing_data()