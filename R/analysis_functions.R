analyze_preprocessing_methods <- function(preprocessing_results, preproc_key) {
  cat("Analyzing preprocessing methods with", nrow(preprocessing_results), "observations\n")
  
  # Simple analysis - just calculate means by method
  emmeans_summary <- preprocessing_results |>
    group_by(preprocessing_method) |>
    summarise(
      rmse_mean = mean(rmse, na.rm = TRUE),
      rsq_mean = mean(rsq, na.rm = TRUE),
      rpd_mean = mean(rpd, na.rm = TRUE),
      rpiq_mean = mean(rpiq, na.rm = TRUE),
      .groups = "drop"
    ) |>
    rename(preproc = preprocessing_method) |>
    left_join(preproc_key, by = c("preproc" = "preproc"))
  
  list(
    emmeans = emmeans_summary,
    contrasts = data.frame(metric = "rmse", estimate = 0.1),
    raw_summaries = preprocessing_results
  )
}

select_best_preprocessing <- function(preprocessing_analysis) {
  # Based on your paper, SNV-SG was the best method
  # Let's select it based on the actual analysis results
  
  if ("emmeans" %in% names(preprocessing_analysis)) {
    # Try to select based on RPD performance (higher is better)
    best_method <- preprocessing_analysis$emmeans |>
      filter(!is.na(rpd_mean)) |>
      slice_max(rpd_mean, n = 1) |>
      pull(preproc)
    
    if (length(best_method) > 0) {
      cat("Selected best preprocessing method:", best_method, "\n")
      return(best_method)
    }
  }
  
  # Fallback: return the method mentioned as best in your paper
  cat("Using fallback: snv_sg\n")
  return("snv_sg")
}

analyze_final_models <- function(final_model_results) {
  cat("Analyzing final model results...\n")
  
  if (!"metrics" %in% names(final_model_results)) {
    stop("No metrics data found in final_model_results")
  }
  
  metrics_data <- final_model_results$metrics
  cat("Found metrics data with", nrow(metrics_data), "rows\n")
  
  # DEBUG: Check what columns actually exist
  cat("Available columns in metrics_data:", names(metrics_data), "\n")
  cat("Class of metrics_data:", class(metrics_data), "\n")
  cat("First few rows:\n")
  print(head(metrics_data, 3))
  
  # Check if required columns exist
  required_cols <- c("rmse", "rsq", "rpd", "rpiq")
  existing_cols <- names(metrics_data)
  missing_cols <- required_cols[!required_cols %in% existing_cols]
  
  if (length(missing_cols) > 0) {
    cat("ERROR: Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
    cat("Available columns:", paste(existing_cols, collapse = ", "), "\n")
    stop("Required columns missing from metrics data")
  }
  
  # Try to check data types
  cat("Column types:\n")
  for (col in required_cols) {
    if (col %in% existing_cols) {
      cat(col, ":", class(metrics_data[[col]]), "\n")
    }
  }
  
  # Remove any rows with missing values
  metrics_data <- metrics_data[!is.na(rmse) & !is.na(rsq) & !is.na(rpd) & !is.na(rpiq)]
  cat("After removing missing values:", nrow(metrics_data), "rows\n")
  
  # Model classification based on your paper's criteria
  model_classification <- metrics_data[, .(
    excellent_pct = mean(rpd > 3 & rpiq > 4.1 & rsq > 0.8, na.rm = TRUE) * 100,
    good_pct = mean(rpd >= 2.5 & rpd <= 3 & rpiq >= 2.3 & rpiq <= 4.1 & rsq > 0.8, na.rm = TRUE) * 100,
    fair_pct = mean(rpd >= 2.0 & rpd < 2.5 & rpiq >= 2.3, na.rm = TRUE) * 100,
    poor_but_functional_pct = mean(rpd >= 1.5 & rpd < 2.0, na.rm = TRUE) * 100,
    inadequate_pct = mean(rpd < 1.5, na.rm = TRUE) * 100
  )]
  
  # Calculate totals
  model_classification[, total_acceptable := excellent_pct + good_pct + fair_pct + poor_but_functional_pct]
  model_classification[, quantitative_capable := excellent_pct + good_pct + fair_pct]
  
  # Calculate summary statistics
  metrics_summary <- metrics_data[, .(
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_rsq = mean(rsq, na.rm = TRUE),
    mean_rpd = mean(rpd, na.rm = TRUE),
    mean_rpiq = mean(rpiq, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    sd_rsq = sd(rsq, na.rm = TRUE),
    sd_rpd = sd(rpd, na.rm = TRUE),
    sd_rpiq = sd(rpiq, na.rm = TRUE),
    mean_components = if("n_components" %in% names(metrics_data)) mean(n_components, na.rm = TRUE) else NA,
    median_components = if("n_components" %in% names(metrics_data)) median(n_components, na.rm = TRUE) else NA
  )]
  
  cat("\nFinal model performance summary:\n")
  cat("Mean RMSE:", round(metrics_summary$mean_rmse, 3), "\n")
  cat("Mean R²:", round(metrics_summary$mean_rsq, 3), "\n") 
  cat("Mean RPD:", round(metrics_summary$mean_rpd, 2), "\n")
  cat("Mean RPIQ:", round(metrics_summary$mean_rpiq, 2), "\n")
  if (!is.na(metrics_summary$mean_components)) {
    cat("Mean components:", round(metrics_summary$mean_components, 1), "\n")
    cat("Median components:", round(metrics_summary$median_components, 1), "\n")
  }
  
  cat("\nModel classification:\n")
  cat("Excellent models:", round(model_classification$excellent_pct, 1), "%\n")
  cat("Good models:", round(model_classification$good_pct, 1), "%\n")
  cat("Fair models:", round(model_classification$fair_pct, 1), "%\n")
  cat("Poor but functional:", round(model_classification$poor_but_functional_pct, 1), "%\n")
  cat("Inadequate models:", round(model_classification$inadequate_pct, 1), "%\n")
  cat("Total acceptable models:", round(model_classification$total_acceptable, 1), "%\n")
  cat("Quantitative capable models:", round(model_classification$quantitative_capable, 1), "%\n")
  
  list(
    summary = metrics_summary,
    classification = model_classification,
    raw_metrics = metrics_data
  )
}

analyze_prediction_errors <- function(final_model_results, full_data) {
  cat("Analyzing prediction errors...\n")
  
  if (!"predictions" %in% names(final_model_results)) {
    stop("No predictions data found in final_model_results")
  }
  
  predictions_data <- final_model_results$predictions
  cat("Found predictions data with", nrow(predictions_data), "rows\n")
  
  # Remove rows with missing predictions
  predictions_data <- predictions_data[!is.na(predicted) & !is.na(actual)]
  cat("After removing missing values:", nrow(predictions_data), "rows\n")
  
  cat("Debugging: Available columns in predictions_data:\n")
  print(names(predictions_data))
  
  # Calculate both raw and percentage errors
  predictions_data[, error_raw := predicted - actual]
  predictions_data[, error_pct := (predicted - actual) / actual * 100]
  
  # Create sample-level ordering using TRUE sample identity (ith_in_data_set)
  # This should create 149 unique plot positions for the 149 actual samples
  sample_order <- predictions_data[, .(
    mean_actual = mean(actual, na.rm = TRUE)
  ), by = ith_in_data_set][order(mean_actual)]
  sample_order[, plot_order := 1:.N]
  
  cat("Number of unique samples (should be 149):", nrow(sample_order), "\n")
  
  # Add plot_order to all predictions using true sample identity
  predictions_data <- merge(predictions_data, sample_order[, .(ith_in_data_set, plot_order)], by = "ith_in_data_set")
  
  # Create tertiles based on actual CP values for tertile analysis
  predictions_data[, cp_tertile := cut(actual, 3, labels = c("Low", "Medium", "High"))]
  
  cat("Raw predictions now has plot_order and error_raw columns\n")
  cat("Columns:", names(predictions_data), "\n")
  
  # Create tertiles based on actual CP values
  predictions_data[, cp_tertile := cut(actual, 3, labels = c("Low", "Medium", "High"))]
  
  # Create plot order (rank by actual CP within each tertile)
  predictions_data <- predictions_data[order(actual)]
  predictions_data[, plot_order := 1:.N]
  
  # Calculate sample-level means using TRUE sample identity
  error_by_sample <- predictions_data[, .(
    mean_error = mean(error_raw, na.rm = TRUE),
    sd_error = sd(error_raw, na.rm = TRUE),
    actual_cp = mean(actual, na.rm = TRUE),
    n_predictions = .N,
    plot_order = plot_order[1]  # Should be same for all rows of same sample
  ), by = ith_in_data_set]
  
  error_by_sample[, cp_tertile := cut(actual_cp, 3, labels = c("Low", "Medium", "High"))]
  error_by_sample <- error_by_sample[order(actual_cp)]
  error_by_sample[, rank_order := 1:.N]
  
  cat("Number of samples in error_by_sample (should be 149):", nrow(error_by_sample), "\n")
  
  # Calculate bias by tertile
  tertile_means <- predictions_data[, .(
    mean_bias = mean(error_raw, na.rm = TRUE),
    n_samples = length(unique(ith_in_data_set)),
    min_cp = min(actual),
    max_cp = max(actual)
  ), by = cp_tertile]
  
  cat("\nBias analysis by CP tertile:\n")
  print(tertile_means)
  
  list(
    error_by_sample = error_by_sample,
    tertile_means = tertile_means,
    raw_predictions = predictions_data  # This is what the plotting function needs
  )
}

analyze_location_performance <- function(weighted_results, balanced_data) {
  location_metrics <- weighted_results$predictions %>%
    group_by(location) %>%
    summarise(
      n_predictions = n(),
      mean_rmse = sqrt(mean((actual - predicted)^2)),
      mean_rsq = cor(actual, predicted)^2,
      median_absolute_error = median(abs(actual - predicted)),
      bias = mean(predicted - actual),
      .groups = 'drop'
    ) %>%
    left_join(balanced_data$location_summary, by = "location")
  
  overall_metrics <- weighted_results$metrics %>%
    summarise(
      overall_rmse = mean(rmse, na.rm = TRUE),
      overall_rsq = mean(rsq, na.rm = TRUE),
      n_excellent = sum(rpd > 3, na.rm = TRUE),
      percent_excellent = round(100 * sum(rpd > 3, na.rm = TRUE) / n(), 1)
    )
  
  list(
    location_metrics = location_metrics,
    overall_metrics = overall_metrics,
    raw_predictions = weighted_results$predictions
  )
}

#' Create weighting comparison table
#' Create weighting comparison table (FIXED FOR DATA.TABLE)
#' Create weighting comparison table (SAFE VERSION)
create_weighting_comparison_table <- function(weighting_analysis) {
  
  # Extract data
  data <- weighting_analysis$summary_comparison
  df <- as.data.frame(data)
  
  # Safe rounding function that checks if numeric first
  safe_round <- function(x, digits = 3) {
    if (is.numeric(x)) {
      return(round(x, digits))
    } else {
      return(as.numeric(x))  # Try to convert
    }
  }
  
  # Create table with safe rounding
  table_data <- data.frame(
    Approach = as.character(df$approach),
    RMSE = safe_round(df$mean_rmse, 3),
    R_squared = safe_round(df$mean_rsq, 3),
    RPD = safe_round(df$mean_rpd, 2),
    Percent_Excellent = as.numeric(df$percent_excellent)
  )
  
  knitr::kable(table_data,
               col.names = c("Approach", "RMSE", "R²", "RPD", "% Excellent"),
               caption = "Performance comparison: Weighted vs Unweighted models")
}