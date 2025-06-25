# R/core_functions.R - Hemp NIR Analysis Core Functions
# All essential functions for the paper, consolidated and streamlined

library(data.table)
library(caret)
library(prospectr)
library(pls)
library(tidyverse)
library(kableExtra)

# =============================================================================
# DATA PREPARATION FUNCTIONS
# =============================================================================

prepare_hemp_data <- function(raw_data) {
  # Prepare hemp data for analysis
  
  data <- copy(raw_data)
  
  # Convert protein from % to g/kg (multiply by 10 for paper units)
  data[, crude_protein := crude_protein * 10]
  
  # Create clean location names
  data[, location := fcase(
    loc %in% c("ithaca", "eith"), "Ithaca",
    loc %in% c("freev", "free"), "Freeville", 
    loc == "chazy", "Chazy",
    loc == "willsboro", "Willsboro",
    loc %in% c("rsgeneva", "rn041_gen", "rn028_gen"), "Geneva",
    default = tools::toTitleCase(tolower(loc))
  )]
  
  cat("Hemp data prepared:\n")
  cat("- Samples:", nrow(data), "\n")
  cat("- Protein range:", round(range(data$crude_protein)), "g/kg\n")
  cat("- Locations:", paste(unique(data$location), collapse = ", "), "\n")
  
  return(data)
}

prepare_preprocessing_key <- function(raw_key) {
  # Create preprocessing method lookup table
  
  key <- copy(raw_key)
  
  # Add human-readable method names
  key[, method_name := c(
    "Raw spectra",
    "First derivative", 
    "Savitzky-Golay smoothing",
    "Gap-segment derivative",
    "Standard normal variate (SNV)",
    "SNV + Savitzky-Golay", 
    "SNV + detrending",
    "Multiplicative scatter correction"
  )]
  
  return(key)
}

# =============================================================================
# CORE MODELING FUNCTIONS
# =============================================================================

run_preprocessing_comparison <- function(data, n_iterations = NULL) {
  # Compare 8 preprocessing methods - core analysis for paper
  
  # Use configuration for iteration count
  config <- get_analysis_config()
  n_iterations <- resolve_param(n_iterations, config$n_iterations, "n_iterations")
  
  cat("=== PREPROCESSING COMPARISON ===\n")
  cat("Methods: 8 preprocessing approaches\n")
  cat("Iterations:", n_iterations, "\n")
  
  # Extract spectral data and response
  spectral_cols <- grep("^x[0-9]+$", names(data), value = TRUE)
  spectra_matrix <- as.matrix(data[, ..spectral_cols])
  y <- data$crude_protein
  
  cat("Data:", nrow(data), "samples,", length(spectral_cols), "wavelengths\n")
  cat("Protein range:", round(range(y)), "g/kg\n")
  
  # Run iterations
  all_results <- list()
  successful_iterations <- 0
  
  for (i in 1:n_iterations) {
    if (i %% max(1, n_iterations %/% 10) == 0) cat("Iteration", i, "\n")
    
    tryCatch({
      # Single iteration comparison
      iteration_result <- compare_preprocessing_single_iteration(spectra_matrix, y)
      iteration_result[, iteration := i]
      all_results[[length(all_results) + 1]] <- iteration_result
      successful_iterations <- successful_iterations + 1
    }, error = function(e) {
      cat("Error in iteration", i, ":", e$message, "\n")
    })
  }
  
  cat("Completed", successful_iterations, "successful iterations\n")
  
  if (length(all_results) == 0) {
    stop("No successful iterations completed")
  }
  
  # Combine all results
  combined_results <- rbindlist(all_results)
  
  return(combined_results)
}

compare_preprocessing_single_iteration <- function(spectra, y) {
  # Single iteration of preprocessing method comparison
  
  # Create train/test split
  set.seed(sample.int(10000, 1))  # Random seed for each iteration
  inTrain <- createDataPartition(y, p = 0.75, list = FALSE, groups = 3)
  
  y_train <- y[inTrain]
  y_test <- y[-inTrain]
  
  # Apply all preprocessing methods
  processed_data <- apply_all_preprocessing_methods(
    spectra[inTrain, , drop = FALSE], 
    spectra[-inTrain, , drop = FALSE]
  )
  
  # Test each preprocessing method
  method_results <- list()
  
  for (method_name in names(processed_data$train)) {
    tryCatch({
      # Get processed data for this method
      X_train <- processed_data$train[[method_name]]
      X_test <- processed_data$test[[method_name]]
      
      # Skip if preprocessing failed (dimensions don't match)
      if (ncol(X_train) != ncol(X_test) || ncol(X_train) < 10) {
        next
      }
      
      # Fit PLS model
      model <- train(
        x = X_train,
        y = y_train,
        method = "pls",
        tuneLength = 15,  # Reasonable component range
        trControl = trainControl(
          method = "cv", 
          number = 5,  # 5-fold CV for speed
          verboseIter = FALSE
        )
      )
      
      # Make predictions
      predictions <- predict(model, newdata = X_test)
      
      # Calculate metrics
      metrics <- calculate_performance_metrics(y_test, predictions)
      metrics[, preprocessing_method := method_name]
      metrics[, n_components := model$bestTune$ncomp]
      
      method_results[[method_name]] <- metrics
      
    }, error = function(e) {
      cat("Error with method", method_name, ":", e$message, "\n")
    })
  }
  
  if (length(method_results) == 0) {
    stop("No preprocessing methods completed successfully")
  }
  
  return(rbindlist(method_results))
}

apply_all_preprocessing_methods <- function(train_spectra, test_spectra) {
  # Apply all 8 preprocessing methods
  
  train_methods <- list()
  test_methods <- list()
  
  # Method 1: Raw spectra
  train_methods$raw <- train_spectra
  test_methods$raw <- test_spectra
  
  # Method 2: First derivative  
  tryCatch({
    train_methods$first_derivative <- t(diff(t(train_spectra), differences = 1))
    test_methods$first_derivative <- t(diff(t(test_spectra), differences = 1))
  }, error = function(e) {
    train_methods$first_derivative <- train_spectra[, -ncol(train_spectra)]
    test_methods$first_derivative <- test_spectra[, -ncol(test_spectra)]
  })
  
  # Method 3: Savitzky-Golay
  tryCatch({
    train_methods$sav_gol <- prospectr::savitzkyGolay(train_spectra, m = 1, p = 3, w = 5)
    test_methods$sav_gol <- prospectr::savitzkyGolay(test_spectra, m = 1, p = 3, w = 5)
  }, error = function(e) {
    train_methods$sav_gol <- train_spectra
    test_methods$sav_gol <- test_spectra
  })
  
  # Method 4: Gap-segment derivative
  tryCatch({
    train_methods$gap_der <- prospectr::gapDer(train_spectra, m = 1, w = 11, s = 5)
    test_methods$gap_der <- prospectr::gapDer(test_spectra, m = 1, w = 11, s = 5)
  }, error = function(e) {
    train_methods$gap_der <- train_spectra
    test_methods$gap_der <- test_spectra
  })
  
  # Method 5: Standard Normal Variate (SNV)
  tryCatch({
    train_methods$snv <- prospectr::standardNormalVariate(train_spectra)
    test_methods$snv <- prospectr::standardNormalVariate(test_spectra)
  }, error = function(e) {
    train_methods$snv <- train_spectra
    test_methods$snv <- test_spectra
  })
  
  # Method 6: SNV + Savitzky-Golay (likely best method)
  tryCatch({
    train_sg <- prospectr::savitzkyGolay(train_spectra, m = 1, p = 3, w = 5)
    test_sg <- prospectr::savitzkyGolay(test_spectra, m = 1, p = 3, w = 5)
    train_methods$snv_sg <- prospectr::standardNormalVariate(train_sg)
    test_methods$snv_sg <- prospectr::standardNormalVariate(test_sg)
  }, error = function(e) {
    train_methods$snv_sg <- train_spectra
    test_methods$snv_sg <- test_spectra
  })
  
  # Method 7: SNV + Detrending  
  tryCatch({
    train_snv <- prospectr::standardNormalVariate(train_spectra)
    test_snv <- prospectr::standardNormalVariate(test_spectra)
    train_methods$snv_detrend <- prospectr::detrend(train_snv, wav = 1:ncol(train_snv))
    test_methods$snv_detrend <- prospectr::detrend(test_snv, wav = 1:ncol(test_snv))
  }, error = function(e) {
    train_methods$snv_detrend <- train_spectra
    test_methods$snv_detrend <- test_spectra
  })
  
  # Method 8: Multiplicative Scatter Correction (MSC)
  tryCatch({
    train_methods$msc <- prospectr::msc(train_spectra)
    # MSC on test set using training reference
    test_methods$msc <- prospectr::msc(test_spectra, reference = attr(train_methods$msc, "Reference spectrum"))
  }, error = function(e) {
    train_methods$msc <- train_spectra
    test_methods$msc <- test_spectra
  })
  
  return(list(train = train_methods, test = test_methods))
}

calculate_performance_metrics <- function(actual, predicted) {
  # Calculate NIR performance metrics
  
  # Remove any NA values
  valid_idx <- !is.na(actual) & !is.na(predicted)
  actual <- actual[valid_idx]
  predicted <- predicted[valid_idx]
  
  if (length(actual) < 3) {
    return(data.table(rmse = NA_real_, rsq = NA_real_, rpd = NA_real_, rpiq = NA_real_))
  }
  
  # Core metrics
  rmse <- sqrt(mean((actual - predicted)^2))
  rsq <- cor(actual, predicted)^2
  
  # NIR-specific metrics
  rpd <- sd(actual) / rmse  # Relative Prediction Deviation
  rpiq <- IQR(actual) / rmse  # Ratio of Performance to InterQuartile distance
  
  data.table(
    rmse = rmse,
    rsq = rsq,
    rpd = rpd,
    rpiq = rpiq
  )
}

run_final_modeling <- function(data, best_method, n_iterations = NULL) {
  # Run final modeling with best preprocessing method
  
  config <- get_analysis_config()
  n_iterations <- resolve_param(n_iterations, config$n_iterations, "n_iterations")
  
  cat("=== FINAL MODELING ===\n")
  cat("Best method:", best_method, "\n")
  cat("Iterations:", n_iterations, "\n")
  
  # Extract data
  spectral_cols <- grep("^x[0-9]+$", names(data), value = TRUE)
  spectra_matrix <- as.matrix(data[, ..spectral_cols])
  y <- data$crude_protein
  
  # Storage for results
  metrics_list <- list()
  predictions_list <- list()
  
  for (i in 1:n_iterations) {
    if (i %% max(1, n_iterations %/% 10) == 0) cat("Iteration", i, "\n")
    
    tryCatch({
      result <- run_single_final_modeling_iteration(spectra_matrix, y, best_method, i)
      metrics_list[[i]] <- result$metrics
      predictions_list[[i]] <- result$predictions
    }, error = function(e) {
      cat("Error in iteration", i, ":", e$message, "\n")
    })
  }
  
  # Combine results
  metrics <- rbindlist(metrics_list, idcol = "iteration")
  predictions <- rbindlist(predictions_list, idcol = "iteration")
  
  cat("Final modeling completed:", nrow(metrics), "successful iterations\n")
  
  return(list(metrics = metrics, predictions = predictions))
}

run_single_final_modeling_iteration <- function(spectra, y, method, seed) {
  # Single iteration of final modeling
  
  set.seed(seed)
  inTrain <- createDataPartition(y, p = 0.75, list = FALSE, groups = 3)
  
  # Apply preprocessing (just the best method)
  processed <- apply_single_preprocessing_method(
    spectra[inTrain, , drop = FALSE], 
    spectra[-inTrain, , drop = FALSE], 
    method
  )
  
  # Fit model with more comprehensive tuning for final analysis
  model <- train(
    x = processed$train,
    y = y[inTrain],
    method = "pls",
    tuneLength = 20,  # More components for final model
    trControl = trainControl(
      method = "cv", 
      number = 10,  # 10-fold CV for final model
      verboseIter = FALSE
    )
  )
  
  # Predictions
  predictions <- predict(model, processed$test)
  actual <- y[-inTrain]
  
  # Metrics
  metrics <- calculate_performance_metrics(actual, predictions)
  metrics[, n_components := model$bestTune$ncomp]
  
  # Detailed predictions for error analysis
  pred_details <- data.table(
    actual = actual,
    predicted = predictions,
    sample_id = which(!seq_along(y) %in% inTrain),
    residual = actual - predictions
  )
  
  return(list(metrics = metrics, predictions = pred_details))
}

apply_single_preprocessing_method <- function(train_spectra, test_spectra, method) {
  # Apply single preprocessing method
  
  switch(method,
         "raw" = list(train = train_spectra, test = test_spectra),
         "snv_sg" = {
           # Most likely best method based on your paper
           train_sg <- prospectr::savitzkyGolay(train_spectra, m = 1, p = 3, w = 5)
           test_sg <- prospectr::savitzkyGolay(test_spectra, m = 1, p = 3, w = 5)
           list(
             train = prospectr::standardNormalVariate(train_sg),
             test = prospectr::standardNormalVariate(test_sg)
           )
         },
         "snv" = {
           list(
             train = prospectr::standardNormalVariate(train_spectra),
             test = prospectr::standardNormalVariate(test_spectra)
           )
         },
         stop("Preprocessing method not implemented: ", method)
  )
}

# =============================================================================
# ANALYSIS FUNCTIONS
# =============================================================================

analyze_preprocessing_methods <- function(results, preproc_key) {
  # Analyze preprocessing comparison results
  
  cat("Analyzing preprocessing results with", nrow(results), "observations\n")
  
  # Calculate summary statistics by method
  summary_stats <- results[, .(
    rmse_mean = mean(rmse, na.rm = TRUE),
    rmse_sd = sd(rmse, na.rm = TRUE),
    rsq_mean = mean(rsq, na.rm = TRUE), 
    rsq_sd = sd(rsq, na.rm = TRUE),
    rpd_mean = mean(rpd, na.rm = TRUE),
    rpd_sd = sd(rpd, na.rm = TRUE),
    rpiq_mean = mean(rpiq, na.rm = TRUE),
    rpiq_sd = sd(rpiq, na.rm = TRUE),
    n_iterations = .N
  ), by = preprocessing_method]
  
  # Merge with method names
  analysis <- merge(summary_stats, preproc_key, 
                    by.x = "preprocessing_method", by.y = "preproc", all.x = TRUE)
  
  # Sort by RPD performance (higher is better)
  setorder(analysis, -rpd_mean)
  
  cat("Top 3 preprocessing methods by RPD:\n")
  print(analysis[1:min(3, nrow(analysis)), .(preprocessing_method, rmse_mean, rsq_mean, rpd_mean)])
  
  return(list(summary = analysis, raw_data = results))
}

select_best_preprocessing_method <- function(analysis) {
  # Select best preprocessing method based on RPD performance
  
  best_method <- analysis$summary[which.max(rpd_mean), preprocessing_method]
  best_stats <- analysis$summary[preprocessing_method == best_method]
  
  cat("=== BEST PREPROCESSING METHOD ===\n")
  cat("Method:", best_method, "\n")
  cat("RPD:", round(best_stats$rpd_mean, 2), "±", round(best_stats$rpd_sd, 2), "\n")
  cat("RMSE:", round(best_stats$rmse_mean, 1), "±", round(best_stats$rmse_sd, 1), "g/kg\n")
  cat("R²:", round(best_stats$rsq_mean, 3), "±", round(best_stats$rsq_sd, 3), "\n")
  
  return(best_method)
}

analyze_final_model_performance <- function(results) {
  # Analyze final model performance and classify models
  
  metrics <- results$metrics
  
  # Overall statistics
  overall_stats <- list(
    mean_rmse = mean(metrics$rmse, na.rm = TRUE),
    sd_rmse = sd(metrics$rmse, na.rm = TRUE),
    mean_rsq = mean(metrics$rsq, na.rm = TRUE),
    sd_rsq = sd(metrics$rsq, na.rm = TRUE),
    mean_rpd = mean(metrics$rpd, na.rm = TRUE),
    sd_rpd = sd(metrics$rpd, na.rm = TRUE),
    mean_rpiq = mean(metrics$rpiq, na.rm = TRUE),
    sd_rpiq = sd(metrics$rpiq, na.rm = TRUE),
    mean_components = mean(metrics$n_components, na.rm = TRUE)
  )
  
  # Model classification based on RPD and RPIQ thresholds
  classification <- metrics[, .(
    excellent = sum(rpd > 3 & rpiq > 4.1, na.rm = TRUE),
    good = sum(rpd >= 2.5 & rpd <= 3 & rpiq >= 2.3 & rpiq <= 4.1, na.rm = TRUE),
    approximate = sum(rpd >= 2.0 & rpd < 2.5, na.rm = TRUE),
    screening = sum(rpd >= 1.5 & rpd < 2.0, na.rm = TRUE),
    poor = sum(rpd < 1.5, na.rm = TRUE),
    total = .N
  )]
  
  # Add percentages
  classification[, excellent_pct := round(100 * excellent / total, 1)]
  classification[, good_pct := round(100 * good / total, 1)]
  classification[, approximate_pct := round(100 * approximate / total, 1)]
  classification[, screening_pct := round(100 * screening / total, 1)]
  
  cat("=== FINAL MODEL PERFORMANCE ===\n")
  cat("RMSE:", round(overall_stats$mean_rmse, 2), "±", round(overall_stats$sd_rmse, 2), "g/kg\n")
  cat("R²:", round(overall_stats$mean_rsq, 3), "±", round(overall_stats$sd_rsq, 3), "\n")
  cat("RPD:", round(overall_stats$mean_rpd, 2), "±", round(overall_stats$sd_rpd, 2), "\n")
  cat("RPIQ:", round(overall_stats$mean_rpiq, 2), "±", round(overall_stats$sd_rpiq, 2), "\n")
  cat("Mean components:", round(overall_stats$mean_components, 1), "\n")
  cat("\nModel Classification:\n")
  cat("- Excellent:", classification$excellent, "(", classification$excellent_pct, "%)\n")
  cat("- Good:", classification$good, "(", classification$good_pct, "%)\n")
  cat("- Approximate:", classification$approximate, "(", classification$approximate_pct, "%)\n")
  
  return(list(
    overall_stats = overall_stats,
    classification = classification,
    raw_metrics = metrics
  ))
}

analyze_prediction_errors <- function(results, hemp_data) {
  # Analyze systematic bias in predictions
  
  predictions <- results$predictions
  
  # Calculate mean error per sample across all iterations
  sample_errors <- predictions[, .(
    mean_error = mean(residual, na.rm = TRUE),
    actual_value = mean(actual, na.rm = TRUE),
    n_predictions = .N
  ), by = sample_id]
  
  # Only keep samples that appeared in multiple iterations
  sample_errors <- sample_errors[n_predictions >= 5]
  
  # Divide into tertiles by protein content for analysis
  sample_errors[, tertile := cut(actual_value, breaks = 3, labels = c("Low", "Medium", "High"))]
  
  # Calculate bias by tertile
  tertile_bias <- sample_errors[, .(
    mean_bias = mean(mean_error, na.rm = TRUE),
    sd_bias = sd(mean_error, na.rm = TRUE),
    n_samples = .N,
    protein_range = paste(round(range(actual_value)), collapse = "-")
  ), by = tertile]
  
  cat("=== SYSTEMATIC BIAS ANALYSIS ===\n")
  print(tertile_bias)
  
  return(list(
    sample_errors = sample_errors,
    tertile_bias = tertile_bias
  ))
}

# =============================================================================
# TABLE AND FIGURE FUNCTIONS
# =============================================================================

create_sample_summary_table <- function(data) {
  # Table 1: Sample characteristics by location
  
  summary_by_location <- data[, .(
    n_samples = .N,
    mean_protein = round(mean(crude_protein), 0),
    sd_protein = round(sd(crude_protein), 0),  
    min_protein = round(min(crude_protein), 0),
    max_protein = round(max(crude_protein), 0),
    n_cultivars = length(unique(cultivar))
  ), by = location]
  
  # Add total row
  total_row <- data[, .(
    location = "Total",
    n_samples = .N,
    mean_protein = round(mean(crude_protein), 0),
    sd_protein = round(sd(crude_protein), 0),
    min_protein = round(min(crude_protein), 0),
    max_protein = round(max(crude_protein), 0),
    n_cultivars = length(unique(cultivar))
  )]
  
  final_table <- rbind(summary_by_location, total_row)
  setnames(final_table, c("Location", "Samples", "Mean CP (g/kg)", "SD CP (g/kg)", 
                          "Min CP (g/kg)", "Max CP (g/kg)", "Cultivars"))
  
  kable(final_table, caption = "Sample characteristics by location") %>%
    kable_styling(bootstrap_options = c("striped", "hover"))
}

create_preprocessing_comparison_table <- function(analysis) {
  # Table 2: Preprocessing method comparison
  
  table_data <- analysis$summary[, .(
    Method = ifelse(is.na(method_name), preprocessing_method, method_name),
    RMSE = paste0(round(rmse_mean, 1), " ± ", round(rmse_sd, 1)),
    R2 = paste0(round(rsq_mean, 3), " ± ", round(rsq_sd, 3)),
    RPD = paste0(round(rpd_mean, 2), " ± ", round(rpd_sd, 2)),
    RPIQ = paste0(round(rpiq_mean, 2), " ± ", round(rpiq_sd, 2))
  )]
  
  setnames(table_data, c("Method", "RMSE (g/kg)", "R²", "RPD", "RPIQ"))
  
  kable(table_data, caption = "Preprocessing method comparison (mean ± SD)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"))
}

create_calibration_plot <- function(results) {
  # Figure 1: Model calibration showing RMSE vs components
  
  # Calculate mean RMSE by number of components
  component_summary <- results$metrics[, .(
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    n_models = .N
  ), by = n_components]
  
  setorder(component_summary, n_components)
  
  ggplot(component_summary, aes(x = n_components, y = mean_rmse)) +
    geom_line(size = 1.2, color = "blue") +
    geom_point(size = 2.5, color = "blue") +
    geom_errorbar(aes(ymin = mean_rmse - sd_rmse, ymax = mean_rmse + sd_rmse), 
                  width = 0.2, color = "blue", alpha = 0.7) +
    labs(
      title = "Model Calibration: RMSE vs Number of PLS Components",
      x = "Number of PLS Components", 
      y = "Root Mean Square Error (g/kg)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10)
    ) +
    scale_x_continuous(breaks = seq(0, max(component_summary$n_components), by = 2))
}

create_performance_boxplot <- function(analysis) {
  # Figure 2: Final model performance distribution
  
  metrics_long <- analysis$raw_metrics %>%
    select(rmse, rsq, rpd, rpiq) %>%
    pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
    mutate(
      metric_label = case_when(
        metric == "rmse" ~ "RMSE (g/kg)",
        metric == "rsq" ~ "R²",
        metric == "rpd" ~ "RPD", 
        metric == "rpiq" ~ "RPIQ"
      )
    )
  
  ggplot(metrics_long, aes(x = metric_label, y = value)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.alpha = 0.5) +
    facet_wrap(~metric_label, scales = "free", ncol = 2) +
    labs(
      title = "Final Model Performance Distribution",
      y = "Metric Value"
    ) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12),
      strip.text = element_text(size = 10)
    )
}

create_validation_error_plot <- function(error_analysis) {
  # Figure 3: Systematic bias validation plot
  
  ggplot(error_analysis$sample_errors, aes(x = actual_value, y = mean_error)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_smooth(method = "lm", color = "red", se = TRUE, size = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.7) +
    facet_wrap(~tertile, scales = "free_x", labeller = label_both) +
    labs(
      title = "Systematic Bias in Predictions by Protein Tertile",
      x = "Actual Protein Content (g/kg)",
      y = "Mean Prediction Error (g/kg)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 10)
    )
}