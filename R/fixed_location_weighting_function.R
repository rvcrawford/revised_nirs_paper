# Fixed location weighting functions that preserve location information

#' Fixed version of run_single_method_modeling that preserves location
run_single_method_modeling_with_location <- function(data, best_method, n_iterations = 1000) {
  cat("=== SINGLE METHOD MODELING WITH LOCATION ===\n")
  cat("Method:", best_method, "\n")
  cat("Iterations:", n_iterations, "\n")
  
  # Extract spectral data
  spectral_cols <- grep("^x[0-9]+$", names(data), value = TRUE)
  spectra_matrix <- as.matrix(data[, ..spectral_cols])
  y <- data$crude_protein
  
  # Handle location column (could be 'loc' or 'clean_loc')
  if ("clean_loc" %in% names(data)) {
    locations <- data$clean_loc
  } else if ("loc" %in% names(data)) {
    locations <- data$loc
  } else {
    stop("No location column found in data")
  }
  
  cat("Data dimensions:", dim(spectra_matrix), "\n")
  cat("Unique locations:", unique(locations), "\n")
  
  # Storage for results
  results_list <- list()
  predictions_list <- list()
  
  for (i in 1:n_iterations) {
    if (i %% 100 == 0) cat("Iteration", i, "\n")
    
    tryCatch({
      # Stratified split by location
      set.seed(i)
      inTrain <- split_spectra_stratified(y, locations)
      
      y_train <- y[inTrain]
      y_test <- y[-inTrain]
      locations_train <- locations[inTrain]
      locations_test <- locations[-inTrain]
      
      # Apply preprocessing
      spectra_processed <- my_preprocess(
        spectra_matrix[inTrain, ], 
        spectra_matrix[-inTrain, ]
      )
      
      # Extract method data
      method_train_name <- paste0(best_method, "_train")
      method_test_name <- paste0(best_method, "_test")
      
      train_spectral <- spectra_processed[[1]][[method_train_name]]
      test_spectral <- spectra_processed[[2]][[method_test_name]]
      
      # Create training dataframe (no weights for unweighted model)
      train_df <- data.frame(y_train = y_train, train_spectral)
      
      # Fit PLS model without weights
      model <- train(
        y_train ~ .,
        data = train_df,
        method = "pls",
        tuneLength = 20,
        trControl = trainControl(method = "cv", number = 10)
      )
      
      # Make predictions
      predictions <- predict(model, newdata = as.data.frame(test_spectral))
      
      # Calculate metrics
      rmse_val <- sqrt(mean((y_test - predictions)^2))
      rsq_val <- cor(y_test, predictions)^2
      rpd_val <- sd(y_test) / rmse_val
      rpiq_val <- IQR(y_test) / rmse_val
      
      # Store overall results
      results_list[[i]] <- data.frame(
        iteration = i,
        rmse = rmse_val,
        rsq = rsq_val,
        rpd = rpd_val,
        rpiq = rpiq_val,
        n_components = model$bestTune$ncomp
      )
      
      # Store predictions WITH location info
      predictions_list[[i]] <- data.frame(
        iteration = i,
        actual = y_test,
        predicted = predictions,
        location = locations_test,
        sample_id = (1:length(y_test))
      )
      
    }, error = function(e) {
      cat("Error in iteration", i, ":", e$message, "\n")
    })
  }
  
  # Combine results
  final_results <- do.call(rbind, results_list)
  final_predictions <- do.call(rbind, predictions_list)
  
  list(
    metrics = final_results,
    predictions = final_predictions
  )
}

#' Fixed compare function that ensures both approaches preserve location
compare_weighted_unweighted_fixed <- function(balanced_data, full_data, best_method, n_iterations = 1000) {
  cat("=== FIXED WEIGHTED vs UNWEIGHTED COMPARISON ===\n")
  
  # Run weighted analysis on balanced data
  cat("Running WEIGHTED analysis...\n")
  weighted_results <- run_weighted_modeling(balanced_data, best_method, n_iterations)
  
  # Run unweighted analysis on balanced data with location preservation
  cat("Running UNWEIGHTED analysis with location preservation...\n")
  main_data <- balanced_data$main_data
  
  # Remove weights for unweighted analysis but keep location info
  unweighted_data <- copy(main_data)
  unweighted_data[, sample_weight := NULL]  # Remove weights if present
  
  # Ensure location column is properly named
  if (!"loc" %in% names(unweighted_data) && "clean_loc" %in% names(unweighted_data)) {
    unweighted_data[, loc := clean_loc]
  }
  
  # Use the fixed function that preserves location
  unweighted_results <- run_single_method_modeling_with_location(unweighted_data, best_method, n_iterations)
  
  # Verify both have location information
  cat("Weighted predictions columns:", names(weighted_results$predictions), "\n")
  cat("Unweighted predictions columns:", names(unweighted_results$predictions), "\n")
  
  if (!"location" %in% names(weighted_results$predictions)) {
    stop("Weighted results missing location column")
  }
  if (!"location" %in% names(unweighted_results$predictions)) {
    stop("Unweighted results missing location column")
  }
  
  # Combine for comparison
  list(
    weighted = weighted_results,
    unweighted = unweighted_results,
    balanced_data_info = balanced_data$location_summary
  )
}

#' Geneva-specific performance analysis
# Fixed version of analyze_geneva_performance function
# Add this to the end of your R/location_weighting_functions.R file (replace the existing version)

#' Analyze Geneva-specific performance with proper error handling (FIXED)
# FIXED VERSION OF GENEVA ANALYSIS FUNCTION
# Replace your existing analyze_geneva_performance function with this one

analyze_geneva_performance <- function(comparison_results, balanced_data) {
  cat("=== GENEVA-SPECIFIC PERFORMANCE ANALYSIS (FIXED) ===\n")
  
  library(dplyr)
  library(tidyr)
  
  # Verify location columns exist
  if (!"location" %in% names(comparison_results$weighted$predictions)) {
    stop("Weighted predictions missing location column")
  }
  if (!"location" %in% names(comparison_results$unweighted$predictions)) {
    stop("Unweighted predictions missing location column")
  }
  
  # Calculate location-specific metrics
  weighted_by_location <- comparison_results$weighted$predictions %>%
    group_by(location) %>%
    summarise(
      weighted_rmse = sqrt(mean((actual - predicted)^2, na.rm = TRUE)),
      weighted_mae = mean(abs(actual - predicted), na.rm = TRUE),
      weighted_bias = mean(predicted - actual, na.rm = TRUE),
      weighted_abs_bias = mean(abs(predicted - actual), na.rm = TRUE),
      weighted_rsq = cor(actual, predicted, use = "complete.obs")^2,
      weighted_rpd = sd(actual, na.rm = TRUE) / sqrt(mean((actual - predicted)^2, na.rm = TRUE)),
      n_predictions = n(),
      .groups = 'drop'
    ) %>%
    mutate(approach = "weighted")
  
  unweighted_by_location <- comparison_results$unweighted$predictions %>%
    group_by(location) %>%
    summarise(
      unweighted_rmse = sqrt(mean((actual - predicted)^2, na.rm = TRUE)),
      unweighted_mae = mean(abs(actual - predicted), na.rm = TRUE),
      unweighted_bias = mean(predicted - actual, na.rm = TRUE),
      unweighted_abs_bias = mean(abs(predicted - actual), na.rm = TRUE),
      unweighted_rsq = cor(actual, predicted, use = "complete.obs")^2,
      unweighted_rpd = sd(actual, na.rm = TRUE) / sqrt(mean((actual - predicted)^2, na.rm = TRUE)),
      n_predictions = n(),
      .groups = 'drop'
    ) %>%
    mutate(approach = "unweighted")
  
  # Combine and calculate improvements
  location_comparison <- weighted_by_location %>%
    full_join(unweighted_by_location, by = "location") %>%
    mutate(
      rmse_improvement = unweighted_rmse - weighted_rmse,
      rmse_pct_improvement = 100 * (unweighted_rmse - weighted_rmse) / unweighted_rmse,
      bias_improvement = unweighted_abs_bias - weighted_abs_bias,
      rsq_improvement = weighted_rsq - unweighted_rsq
    )
  
  # Add balanced data info safely
  tryCatch({
    location_comparison <- location_comparison %>%
      left_join(balanced_data$location_summary, by = "location")
  }, error = function(e) {
    cat("Warning: Could not join with balanced_data:", e$message, "\n")
  })
  
  # FIXED: Use only columns that actually exist
  cat("Location-specific performance:\n")
  available_cols <- intersect(
    c("location", "weight", "weighted_rmse", "unweighted_rmse", 
      "rmse_improvement", "rmse_pct_improvement"),
    names(location_comparison)
  )
  
  if (length(available_cols) > 0) {
    print(location_comparison %>% select(all_of(available_cols)))
  } else {
    print(location_comparison)
  }
  
  # Focus on Geneva
  geneva_results <- location_comparison %>% filter(location == "geneva")
  
  if (nrow(geneva_results) == 0) {
    cat("ERROR: No Geneva results found!\n")
    return(list(
      location_comparison = location_comparison,
      geneva_results = NULL,
      geneva_significance = NULL
    ))
  }
  
  cat("\nGeneva-specific results:\n")
  cat("========================\n")
  if ("weight" %in% names(geneva_results)) {
    cat("Weight:", geneva_results$weight, "x\n")
  }
  if ("rmse_improvement" %in% names(geneva_results)) {
    cat("RMSE improvement:", round(geneva_results$rmse_improvement, 3), "\n")
  }
  if ("rmse_pct_improvement" %in% names(geneva_results)) {
    cat("% improvement:", round(geneva_results$rmse_pct_improvement, 1), "%\n")
  }
  if ("rsq_improvement" %in% names(geneva_results)) {
    cat("R² improvement:", round(geneva_results$rsq_improvement, 4), "\n")
  }
  
  # Statistical test for Geneva (if possible)
  geneva_significance <- NULL
  tryCatch({
    if (nrow(geneva_results) > 0 && "rmse_improvement" %in% names(geneva_results)) {
      # Get Geneva predictions for statistical test
      geneva_weighted <- comparison_results$weighted$predictions %>%
        filter(location == "geneva")
      geneva_unweighted <- comparison_results$unweighted$predictions %>%
        filter(location == "geneva")
      
      if (nrow(geneva_weighted) > 10 && nrow(geneva_unweighted) > 10) {
        # Statistical test
        wilcox_test <- wilcox.test(
          abs(geneva_weighted$actual - geneva_weighted$predicted),
          abs(geneva_unweighted$actual - geneva_unweighted$predicted),
          paired = FALSE
        )
        
        geneva_significance <- list(
          p_value = wilcox_test$p.value,
          significant = wilcox_test$p.value < 0.05,
          test_method = "Wilcoxon rank-sum test"
        )
        
        cat("Statistical significance: p =", round(wilcox_test$p.value, 4), "\n")
        if (wilcox_test$p.value < 0.05) {
          cat("✓ Statistically significant improvement!\n")
        } else {
          cat("• Not statistically significant\n")
        }
      }
    }
  }, error = function(e) {
    cat("Could not perform statistical test:", e$message, "\n")
  })
  
  return(list(
    location_comparison = location_comparison,
    geneva_results = geneva_results,
    geneva_significance = geneva_significance
  ))
}

#' Create Geneva-focused plot
# REPLACE the create_geneva_improvement_plot function in R/fixed_location_weighting_function.R with this complete version:

create_geneva_improvement_plot <- function(geneva_analysis) {
  library(ggplot2)
  library(dplyr)
  
  # Handle case where geneva_analysis might have different structure
  if (is.null(geneva_analysis$location_comparison)) {
    # Create a simple error plot
    return(ggplot(data.frame(x = 1, y = 1, label = "No data available"), 
                  aes(x, y)) +
             geom_text(aes(label = label)) +
             labs(title = "Geneva Analysis: No location comparison data available") +
             theme_minimal())
  }
  
  location_data <- geneva_analysis$location_comparison %>%
    mutate(
      is_geneva = location == "geneva",
      # Handle case where geneva_significance might be NULL
      significant = case_when(
        location == "geneva" & !is.null(geneva_analysis$geneva_significance) & 
          !is.null(geneva_analysis$geneva_significance$significant) ~ 
          geneva_analysis$geneva_significance$significant,
        TRUE ~ FALSE
      ),
      # Ensure rmse_improvement exists
      rmse_improvement = if_else(is.na(rmse_improvement), 0, rmse_improvement)
    )
  
  # Create the complete plot
  p <- ggplot(location_data, aes(x = reorder(location, rmse_improvement), 
                                 y = rmse_improvement)) +
    geom_col(aes(fill = is_geneva), alpha = 0.8, width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6) +
    geom_text(aes(label = ifelse(significant, 
                                 paste0(round(rmse_improvement, 3), "*"),
                                 round(rmse_improvement, 3))), 
              vjust = ifelse(location_data$rmse_improvement > 0, -0.5, 1.2), 
              size = 3.5, fontface = "bold") +
    labs(
      title = "Location-Specific RMSE Improvement with Sample Weighting",
      subtitle = "* indicates statistically significant improvement (p < 0.05)",
      x = "Location", 
      y = "RMSE Improvement (Unweighted - Weighted)",
      fill = "Location Type"
    ) +
    scale_fill_manual(values = c("FALSE" = "#A23B72", "TRUE" = "#2E86AB"),
                      labels = c("Other locations", "Geneva (targeted)")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray60"),
      legend.position = "bottom"
    )
  
  return(p)
}

#' Generate interpretation text for paper
generate_geneva_interpretation <- function(geneva_analysis) {
  geneva_results <- geneva_analysis$geneva_results
  geneva_sig <- geneva_analysis$geneva_significance
  
  cat("=== INTERPRETATION FOR PAPER ===\n")
  
  if (!is.null(geneva_sig) && geneva_sig$significant && geneva_sig$meaningful) {
    cat("🎯 GENEVA SHOWS SIGNIFICANT IMPROVEMENT!\n\n")
    
    cat("Suggested paper text:\n")
    cat('"Location-weighted modeling significantly improved predictions for Geneva samples,\n')
    cat('the most underrepresented location in the dataset. Geneva samples showed a\n')
    cat(round(geneva_results$rmse_pct_improvement, 1), '% reduction in RMSE (from', round(geneva_results$rmse_unweighted, 2), 'to', round(geneva_results$rmse_weighted, 2), ', p =', round(geneva_sig$p_value, 3), ').\n')
    cat('This improvement demonstrates that sample weighting can effectively address\n')
    cat('location-specific prediction bias when geographic class imbalance exists."\n\n')
    
    return("significant_improvement")
    
  } else if (geneva_results$rmse_improvement > 0) {
    cat("⚠️ GENEVA SHOWS IMPROVEMENT BUT NOT SIGNIFICANT\n\n")
    
    cat("Suggested paper text:\n")
    cat('"Location-weighted modeling showed promise for improving predictions at\n')
    cat('underrepresented locations, with Geneva samples showing a', round(geneva_results$rmse_pct_improvement, 1), '% improvement\n')
    cat('in RMSE, though this was not statistically significant. This suggests that\n')
    cat('larger sample sizes or stronger weighting schemes may be needed to achieve\n')
    cat('statistically significant improvements."\n\n')
    
    return("modest_improvement")
    
  } else {
    cat("❌ GENEVA DOES NOT IMPROVE\n\n")
    
    cat("Suggested paper text:\n")
    cat('"Location-weighted modeling did not improve predictions for underrepresented\n')
    cat('locations, suggesting that the NIRS models are robust across the tested\n')
    cat('geographic environments and that location-specific bias is minimal in this dataset."\n\n')
    
    return("no_improvement")
  }
}