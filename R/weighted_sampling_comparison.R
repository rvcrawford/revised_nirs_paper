# SIMPLE BASE R VERSION - Replace the entire R/weighted_sampling_comparison.R file with this

library(pls)
library(prospectr)
library(dplyr)

# Simple version using base R pls package
my_comparison_weighted_sampling <- function(data, full_data, seed = 123) {
  library(pls)
  library(prospectr)
  library(dplyr)
  
  # Get or create weights
  if("sample_weight" %in% names(data)) {
    case_weights <- data$sample_weight
  } else {
    # Create location-based weights
    if("clean_loc" %in% names(data)) {
      locations <- data$clean_loc
    } else if("loc" %in% names(data)) {
      locations <- data$loc
    } else {
      stop("No location column found")
    }
    
    location_counts <- table(locations)
    max_count <- max(location_counts)
    location_weights <- max_count / location_counts
    case_weights <- as.numeric(location_weights[locations])
  }
  
  # Prepare spectral data
  spectra_cols <- grepl("^x", names(data))
  spectra_matrix <- as.matrix(data[, spectra_cols, with = FALSE])
  
  # Apply preprocessing
  sg_processed <- prospectr::savitzkyGolay(spectra_matrix, m = 1, p = 3, w = 5)
  final_spectra_matrix <- prospectr::standardNormalVariate(sg_processed)
  
  # Create train/test split
  set.seed(seed)
  n_samples <- nrow(data)
  train_indices <- sample(1:n_samples, size = floor(0.7 * n_samples))
  test_indices <- setdiff(1:n_samples, train_indices)
  
  # Prepare data
  train_spectra <- final_spectra_matrix[train_indices, , drop = FALSE]
  test_spectra <- final_spectra_matrix[test_indices, , drop = FALSE]
  train_y <- data$crude_protein[train_indices]
  test_y <- data$crude_protein[test_indices]
  train_weights <- case_weights[train_indices]
  
  # Create data frames for pls
  train_df <- data.frame(y = train_y, train_spectra)
  test_df <- data.frame(y = test_y, test_spectra)
  
  # Fit unweighted model
  unweighted_model <- pls::plsr(y ~ ., data = train_df, ncomp = 12, validation = "none")
  
  # Fit weighted model using weighted least squares approach
  # Create weights that sum to the number of observations
  normalized_weights <- train_weights / mean(train_weights)
  
  # Use pls with weights argument
  weighted_model <- pls::plsr(y ~ ., data = train_df, ncomp = 12, validation = "none", 
                              weights = normalized_weights)
  
  # Make predictions using 12 components
  pred_unweighted <- as.numeric(predict(unweighted_model, newdata = test_df, ncomp = 12))
  pred_weighted <- as.numeric(predict(weighted_model, newdata = test_df, ncomp = 12))
  
  # Calculate metrics
  actual_clean <- test_y[!is.na(test_y)]
  pred_unweighted_clean <- pred_unweighted[!is.na(test_y)]
  pred_weighted_clean <- pred_weighted[!is.na(test_y)]
  
  # Ensure all vectors are the same length
  min_length <- min(length(actual_clean), length(pred_unweighted_clean), length(pred_weighted_clean))
  actual_clean <- actual_clean[1:min_length]
  pred_unweighted_clean <- pred_unweighted_clean[1:min_length]
  pred_weighted_clean <- pred_weighted_clean[1:min_length]
  
  calc_metrics <- function(actual, predicted) {
    if(length(actual) < 3 || length(predicted) < 3) {
      return(c(rmse = NA, mae = NA, rsq = NA, rpd = NA, rpiq = NA, bias = NA))
    }
    
    rmse <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
    mae <- mean(abs(actual - predicted), na.rm = TRUE)
    rsq <- cor(actual, predicted, use = "complete.obs")^2
    rpd <- sd(actual, na.rm = TRUE) / rmse
    rpiq <- IQR(actual, na.rm = TRUE) / rmse
    bias <- mean(predicted - actual, na.rm = TRUE)
    
    c(rmse = rmse, mae = mae, rsq = rsq, rpd = rpd, rpiq = rpiq, bias = bias)
  }
  
  unweighted_metrics <- calc_metrics(actual_clean, pred_unweighted_clean)
  weighted_metrics <- calc_metrics(actual_clean, pred_weighted_clean)
  
  # Create detailed predictions
  test_data_info <- data[test_indices, ]
  valid_idx <- !is.na(test_y)
  
  detailed_predictions <- data.frame(
    iteration = seed,
    sample_id = if("ith_in_data_set" %in% names(test_data_info)) {
      test_data_info$ith_in_data_set[valid_idx][1:min_length]
    } else {
      test_indices[valid_idx][1:min_length]
    },
    actual = actual_clean,
    pred_unweighted = pred_unweighted_clean,
    pred_weighted = pred_weighted_clean,
    improvement = abs(actual_clean - pred_unweighted_clean) - abs(actual_clean - pred_weighted_clean)
  )
  
  # Add location if available
  if("clean_loc" %in% names(test_data_info)) {
    detailed_predictions$location <- test_data_info$clean_loc[valid_idx][1:min_length]
  } else if("loc" %in% names(test_data_info)) {
    detailed_predictions$location <- test_data_info$loc[valid_idx][1:min_length]
  }
  
  return(list(
    iteration = seed,
    n_samples = length(actual_clean),
    unweighted = unweighted_metrics,
    weighted = weighted_metrics,
    improvements = weighted_metrics - unweighted_metrics,
    detailed_predictions = detailed_predictions,
    models = list(unweighted = unweighted_model, weighted = weighted_model)
  ))
}

# Function to run multiple iterations
run_multiple_comparisons <- function(data, full_data, n_iterations = 1000) {
  library(dplyr)
  
  cat("=== RUNNING", n_iterations, "WEIGHTED vs UNWEIGHTED COMPARISONS ===\n")
  
  # Storage for results
  results_list <- list()
  predictions_list <- list()
  successful_iterations <- 0
  
  # Run iterations
  for(i in 1:n_iterations) {
    if(i %% 5 == 0) cat("Iteration", i, "of", n_iterations, "\n")
    
    tryCatch({
      result <- my_comparison_weighted_sampling(data, full_data, seed = i)
      
      # Check if result is valid
      if(!is.null(result) && !is.null(result$unweighted) && !is.null(result$weighted) &&
         !any(is.na(result$unweighted)) && !any(is.na(result$weighted))) {
        
        results_list[[length(results_list) + 1]] <- data.frame(
          iteration = i,
          n_samples = result$n_samples,
          
          # Unweighted metrics
          unw_rmse = result$unweighted["rmse"],
          unw_mae = result$unweighted["mae"], 
          unw_rsq = result$unweighted["rsq"],
          unw_rpd = result$unweighted["rpd"],
          unw_rpiq = result$unweighted["rpiq"],
          unw_bias = result$unweighted["bias"],
          
          # Weighted metrics
          w_rmse = result$weighted["rmse"],
          w_mae = result$weighted["mae"],
          w_rsq = result$weighted["rsq"], 
          w_rpd = result$weighted["rpd"],
          w_rpiq = result$weighted["rpiq"],
          w_bias = result$weighted["bias"],
          
          # Improvements
          rmse_improvement = result$unweighted["rmse"] - result$weighted["rmse"],
          mae_improvement = result$unweighted["mae"] - result$weighted["mae"],
          rsq_improvement = result$improvements["rsq"],
          rpd_improvement = result$improvements["rpd"],
          rpiq_improvement = result$improvements["rpiq"],
          bias_improvement = abs(result$unweighted["bias"]) - abs(result$weighted["bias"])
        )
        
        # Store predictions
        if(!is.null(result$detailed_predictions) && nrow(result$detailed_predictions) > 0) {
          predictions_list[[length(predictions_list) + 1]] <- result$detailed_predictions
        }
        
        successful_iterations <- successful_iterations + 1
      }
      
    }, error = function(e) {
      cat("Error in iteration", i, ":", e$message, "\n")
    })
  }
  
  cat("Completed", successful_iterations, "successful iterations out of", n_iterations, "\n")
  
  # Check if we have any results
  if(length(results_list) == 0) {
    cat("ERROR: No successful iterations!\n")
    return(NULL)
  }
  
  # Combine results
  all_metrics <- do.call(rbind, results_list)
  all_predictions <- if(length(predictions_list) > 0) {
    do.call(rbind, predictions_list)
  } else {
    data.frame()
  }
  
  cat("Final metrics dimensions:", dim(all_metrics), "\n")
  cat("Final predictions dimensions:", dim(all_predictions), "\n")
  
  return(list(
    metrics = all_metrics,
    predictions = all_predictions
  ))
}

# Function to summarize multiple comparison results
summarize_multiple_comparisons <- function(results) {
  library(dplyr)
  
  # Check if results is NULL or missing metrics
  if(is.null(results)) {
    stop("Results object is NULL")
  }
  
  if(is.null(results$metrics)) {
    stop("Results$metrics is NULL")
  }
  
  metrics <- results$metrics
  
  if(nrow(metrics) == 0) {
    stop("No data in metrics")
  }
  
  cat("Processing", nrow(metrics), "iterations of results\n")
  
  # Calculate summary statistics
  summary_stats <- metrics %>%
    summarise(
      n_iterations = n(),
      
      # Mean performance
      mean_unw_rmse = round(mean(unw_rmse, na.rm = TRUE), 3),
      mean_w_rmse = round(mean(w_rmse, na.rm = TRUE), 3),
      mean_unw_rsq = round(mean(unw_rsq, na.rm = TRUE), 3),
      mean_w_rsq = round(mean(w_rsq, na.rm = TRUE), 3),
      mean_unw_rpd = round(mean(unw_rpd, na.rm = TRUE), 3),
      mean_w_rpd = round(mean(w_rpd, na.rm = TRUE), 3),
      
      # Mean improvements
      mean_rmse_improvement = round(mean(rmse_improvement, na.rm = TRUE), 3),
      mean_rsq_improvement = round(mean(rsq_improvement, na.rm = TRUE), 3),
      mean_rpd_improvement = round(mean(rpd_improvement, na.rm = TRUE), 3),
      
      # Standard deviations
      sd_rmse_improvement = round(sd(rmse_improvement, na.rm = TRUE), 3),
      sd_rsq_improvement = round(sd(rsq_improvement, na.rm = TRUE), 3),
      sd_rpd_improvement = round(sd(rpd_improvement, na.rm = TRUE), 3),
      
      # Percentage of iterations where weighted performed better
      pct_rmse_better = round(100 * mean(rmse_improvement > 0, na.rm = TRUE), 1),
      pct_rsq_better = round(100 * mean(rsq_improvement > 0, na.rm = TRUE), 1),
      pct_rpd_better = round(100 * mean(rpd_improvement > 0, na.rm = TRUE), 1),
      
      # Confidence intervals (95%)
      rmse_improvement_ci_lower = round(quantile(rmse_improvement, 0.025, na.rm = TRUE), 3),
      rmse_improvement_ci_upper = round(quantile(rmse_improvement, 0.975, na.rm = TRUE), 3),
      rsq_improvement_ci_lower = round(quantile(rsq_improvement, 0.025, na.rm = TRUE), 3),
      rsq_improvement_ci_upper = round(quantile(rsq_improvement, 0.975, na.rm = TRUE), 3)
    )
  
  # Create performance comparison table
  performance_table <- data.frame(
    Metric = c("RMSE", "R²", "RPD", "RPIQ"),
    Unweighted_Mean = c(summary_stats$mean_unw_rmse, summary_stats$mean_unw_rsq, 
                        summary_stats$mean_unw_rpd, round(mean(metrics$unw_rpiq, na.rm = TRUE), 3)),
    Weighted_Mean = c(summary_stats$mean_w_rmse, summary_stats$mean_w_rsq,
                      summary_stats$mean_w_rpd, round(mean(metrics$w_rpiq, na.rm = TRUE), 3)),
    Mean_Improvement = c(summary_stats$mean_rmse_improvement, summary_stats$mean_rsq_improvement,
                         summary_stats$mean_rpd_improvement, round(mean(metrics$rpiq_improvement, na.rm = TRUE), 3)),
    Pct_Better = c(summary_stats$pct_rmse_better, summary_stats$pct_rsq_better, 
                   summary_stats$pct_rpd_better, round(100 * mean(metrics$rpiq_improvement > 0, na.rm = TRUE), 1))
  )
  
  # Print results
  cat("=== SUMMARY OF", summary_stats$n_iterations, "ITERATIONS ===\n\n")
  
  cat("PERFORMANCE COMPARISON:\n")
  print(performance_table)
  
  cat("\nKEY FINDINGS:\n")
  cat("- Mean RMSE improvement:", summary_stats$mean_rmse_improvement, "±", summary_stats$sd_rmse_improvement, "\n")
  cat("- Mean R² improvement:", summary_stats$mean_rsq_improvement, "±", summary_stats$sd_rsq_improvement, "\n")
  cat("- Mean RPD improvement:", summary_stats$mean_rpd_improvement, "±", summary_stats$sd_rpd_improvement, "\n")
  
  cat("\nCONFIDENCE INTERVALS (95%):\n")
  cat("- RMSE improvement: [", summary_stats$rmse_improvement_ci_lower, ",", summary_stats$rmse_improvement_ci_upper, "]\n")
  cat("- R² improvement: [", summary_stats$rsq_improvement_ci_lower, ",", summary_stats$rsq_improvement_ci_upper, "]\n")
  
  cat("\nCONSISTENCY:\n")
  cat("- Weighted better RMSE in", summary_stats$pct_rmse_better, "% of iterations\n")
  cat("- Weighted better R² in", summary_stats$pct_rsq_better, "% of iterations\n")
  cat("- Weighted better RPD in", summary_stats$pct_rpd_better, "% of iterations\n")
  
  return(list(
    summary_stats = summary_stats,
    performance_table = performance_table,
    full_metrics = metrics,
    full_predictions = results$predictions
  ))
}

# Test function
test_simple_version <- function() {
  cat("=== TESTING SIMPLE BASE R VERSION ===\n")
  
  # Load data
  balanced_data <- tar_read(balanced_data)
  full_data <- tar_read(full_data)
  main_data <- balanced_data$main_data
  
  cat("Testing single iteration...\n")
  tryCatch({
    result <- my_comparison_weighted_sampling(main_data, full_data, seed = 1)
    cat("✓ Simple version SUCCESS!\n")
    cat("- n_samples:", result$n_samples, "\n")
    cat("- Unweighted RMSE:", round(result$unweighted["rmse"], 3), "\n")
    cat("- Weighted RMSE:", round(result$weighted["rmse"], 3), "\n")
    cat("- RMSE improvement:", round(result$improvements["rmse"], 3), "\n")
    
    # Test small comparison
    cat("\nTesting small comparison (2 iterations)...\n")
    small_result <- run_multiple_comparisons(main_data, full_data, n_iterations = 2)
    if(!is.null(small_result)) {
      cat("✓ Small comparison SUCCESS!\n")
      cat("- Got", nrow(small_result$metrics), "successful iterations\n")
      return(TRUE)
    } else {
      cat("✗ Small comparison returned NULL\n")
      return(FALSE)
    }
    
  }, error = function(e) {
    cat("✗ ERROR:", e$message, "\n")
    print(e)
    return(FALSE)
  })
}