
# Function to run multiple iterations
run_multiple_comparisons <- function(data, full_data, n_iterations = 1000) {
  library(dplyr)
  
  cat("=== RUNNING", n_iterations, "WEIGHTED vs UNWEIGHTED COMPARISONS ===\n")
  
  # Storage for results
  results_list <- list()
  predictions_list <- list()
  
  # Run iterations
  for(i in 1:n_iterations) {
    if(i %% 100 == 0) cat("Iteration", i, "of", n_iterations, "\n")
    
    tryCatch({
      result <- my_comparison_weighted_sampling(data, full_data, seed = i)
      
      # Store metrics
      results_list[[i]] <- data.frame(
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
        
        # Improvements (weighted - unweighted, except RMSE/MAE where lower is better)
        rmse_improvement = result$unweighted["rmse"] - result$weighted["rmse"],  # Positive = weighted better
        mae_improvement = result$unweighted["mae"] - result$weighted["mae"],
        rsq_improvement = result$improvements["rsq"],
        rpd_improvement = result$improvements["rpd"],
        rpiq_improvement = result$improvements["rpiq"],
        bias_improvement = abs(result$unweighted["bias"]) - abs(result$weighted["bias"])  # Improvement in bias
      )
      
      # Store predictions
      predictions_list[[i]] <- result$detailed_predictions
      
    }, error = function(e) {
      cat("Error in iteration", i, ":", e$message, "\n")
    })
  }
  
  # Combine results
  all_metrics <- do.call(rbind, results_list)
  all_predictions <- do.call(rbind, predictions_list)
  
  cat("Completed", nrow(all_metrics), "successful iterations\n")
  
  return(list(
    metrics = all_metrics,
    predictions = all_predictions
  ))
}

# Function to summarize multiple comparison results
# CORRECTED VERSION - Replace the broken function with this one
# Put this in R/weighted_sampling_comparison.R

summarize_multiple_comparisons <- function(results) {
  library(dplyr)
  
  metrics <- results$metrics
  
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
      
      # Confidence intervals (95%) - FIXED THE TRUNCATED LINE HERE
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
# FIXED VERSION - Make sure you run this to replace the old function
# FIXED VERSION - Copy and paste this to replace the broken function
# FIXED VERSION - Make sure you run this to replace the old function
# FIXED VERSION - Make sure you run this to replace the old function
# R/weighted_sampling_comparison.R

# CORRECTED VERSION - put this in the file that targets sources
my_comparison_weighted_sampling <- function(data, full_data, seed = 123) {
  library(tidymodels)
  library(plsmod)
  library(prospectr)
  library(dplyr)
  
  # Get weights
  weight_cols <- c("sample_weight", "sample_weights", "case_weights", "weights", "weight")
  case_weights <- NULL
  found_weights <- FALSE
  
  for(col in weight_cols) {
    if(col %in% names(data)) {
      case_weights <- data[[col]]
      found_weights <- TRUE
      break
    }
  }
  
  if(!found_weights || all(case_weights == case_weights[1])) {
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
  
  # Prepare data - FIXED SYNTAX
  spectra_matrix <- data %>% 
    select(starts_with("x"))
  spectra_matrix <- as.matrix(spectra_matrix)
  
  # Apply preprocessing
  sg_processed <- prospectr::savitzkyGolay(spectra_matrix, m = 1, p = 3, w = 5)
  final_spectra_matrix <- prospectr::standardNormalVariate(sg_processed)
  
  # Create train/test split using provided seed
  set.seed(seed)
  inds <- createDataPartition(y = data$crude_protein, p = 0.75, list = FALSE)
  
  # Build training data
  training_data <- data.frame(
    crude_protein = data$crude_protein[inds],
    final_spectra_matrix[inds, ]
  )
  
  # Extract training weights
  train_weights <- case_weights[inds]
  
  # Create weighted training data by sampling with replacement
  sample_probs <- train_weights / sum(train_weights)
  
  # Use a different seed for sampling (derived from main seed)
  set.seed(seed + 1000)
  weighted_indices <- sample(1:nrow(training_data), 
                             size = nrow(training_data), 
                             replace = TRUE, 
                             prob = sample_probs)
  
  weighted_training_data <- training_data[weighted_indices, ]
  
  # Prepare test data
  training_ids <- data$ith_in_data_set[inds]
  prediction_samples <- full_data[!full_data$ith_in_data_set %in% training_ids, ]
  
  # FIXED SYNTAX
  pred_spectra_matrix <- prediction_samples %>% 
    select(starts_with("x"))
  pred_spectra_matrix <- as.matrix(pred_spectra_matrix)
  
  pred_sg_processed <- prospectr::savitzkyGolay(pred_spectra_matrix, m = 1, p = 3, w = 5)
  pred_final_spectra <- prospectr::standardNormalVariate(pred_sg_processed)
  
  test_data <- data.frame(
    crude_protein = prediction_samples$crude_protein,
    pred_final_spectra
  )
  
  # Define PLS model specification
  pls_spec <- pls() %>%
    set_engine("mixOmics") %>%
    set_mode("regression") %>%
    set_args(num_comp = 12)
  
  # Create recipe
  basic_recipe <- recipe(crude_protein ~ ., data = training_data)
  
  # Create workflows
  pls_workflow <- workflow() %>%
    add_model(pls_spec) %>%
    add_recipe(basic_recipe)
  
  # Model 1: UNWEIGHTED
  unweighted_fit <- fit(pls_workflow, data = training_data)
  
  # Model 2: WEIGHTED 
  weighted_fit <- fit(pls_workflow, data = weighted_training_data)
  
  # Make predictions
  pred_unweighted <- predict(unweighted_fit, new_data = test_data)$.pred
  pred_weighted <- predict(weighted_fit, new_data = test_data)$.pred
  
  # Calculate performance metrics
  actual_values <- test_data$crude_protein
  valid_idx <- !is.na(actual_values)
  
  actual_clean <- actual_values[valid_idx]
  pred_unweighted_clean <- pred_unweighted[valid_idx]
  pred_weighted_clean <- pred_weighted[valid_idx]
  
  # Calculate metrics for both models
  calc_metrics <- function(actual, predicted) {
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
  
  # Create detailed predictions with location info
  detailed_predictions <- data.frame(
    iteration = seed,
    sample_id = prediction_samples$ith_in_data_set[valid_idx],
    actual = actual_clean,
    pred_unweighted = pred_unweighted_clean,
    pred_weighted = pred_weighted_clean,
    improvement = abs(actual_clean - pred_unweighted_clean) - abs(actual_clean - pred_weighted_clean)
  )
  
  # Add location if available
  if("clean_loc" %in% names(prediction_samples)) {
    detailed_predictions$location <- prediction_samples$clean_loc[valid_idx]
  }
  
  # Return results in format suitable for aggregation
  return(list(
    iteration = seed,
    n_samples = length(actual_clean),
    unweighted = unweighted_metrics,
    weighted = weighted_metrics,
    improvements = weighted_metrics - unweighted_metrics,
    detailed_predictions = detailed_predictions,
    models = list(unweighted = unweighted_fit, weighted = weighted_fit)
  ))
}

# Put the other corrected functions here too...
# [Include run_multiple_comparisons and summarize_multiple_comparisons from the artifact]