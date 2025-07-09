# R/core_functions.R - SIMPLIFIED Hemp NIR Analysis Functions
# Focus: Essential functions only for paper pipeline

library(data.table)
library(caret)
library(prospectr)
library(pls)
library(tidyverse)
library(kableExtra)

# =============================================================================
# HELPER FUNCTIONS (for legacy compatibility)
# =============================================================================

# Utility function for skewness calculation
calculate_skewness <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  sd_x <- sqrt(sum((x - mean_x)^2) / n)
  z <- (x - mean_x) / sd_x
  sum(z^3) / n
}

my_preprocess <- function(spectra_train, spectra_test) {
  # Legacy preprocessing function for multi-algorithm comparison
  # Returns list with train and test data for each method
  
  train_methods <- list(
    raw = spectra_train,
    first_derivative = apply_preprocessing(spectra_train, 2),
    sav_gol = apply_preprocessing(spectra_train, 3),
    gap_der = apply_preprocessing(spectra_train, 4),
    snv = apply_preprocessing(spectra_train, 5),
    snv_sg = apply_preprocessing(spectra_train, 6),
    snv_detrend = apply_preprocessing(spectra_train, 7),
    msc = apply_preprocessing(spectra_train, 8)
  )
  
  test_methods <- list(
    raw = spectra_test,
    first_derivative = apply_preprocessing(spectra_test, 2),
    sav_gol = apply_preprocessing(spectra_test, 3),
    gap_der = apply_preprocessing(spectra_test, 4),
    snv = apply_preprocessing(spectra_test, 5),
    snv_sg = apply_preprocessing(spectra_test, 6),
    snv_detrend = apply_preprocessing(spectra_test, 7),
    msc = apply_preprocessing(spectra_test, 8)
  )
  
  return(list(train_methods, test_methods))
}

split_spectra <- function(y, p = 0.75) {
  # Legacy function for creating train/test splits
  # Returns logical indices for training set
  train_indices <- createDataPartition(y, p = p, list = FALSE, groups = 3)[,1]
  return(train_indices)
}

# =============================================================================
# CONFIGURATION AND SETUP
# =============================================================================

get_analysis_config <- function() {
  # Get analysis configuration
  list(
    n_iterations = 100,
    test_size = 0.25,
    max_components = 20,
    cv_folds = 10
  )
}

resolve_param <- function(user_value, config_value, param_name) {
  # Resolve parameter value (user > config > default)
  if (!is.null(user_value)) {
    return(user_value)
  } else {
    return(config_value)
  }
}

# =============================================================================
# DATA PREPARATION FUNCTIONS
# =============================================================================

prepare_hemp_data <- function(raw_data) {
  # Clean and prepare hemp data
  data <- copy(raw_data)
  
  # Convert protein from % to g/kg
  data[, crude_protein := crude_protein * 10]
  
  # Standardize location names
  data[, location := fcase(
    loc %in% c("ithaca", "eith"), "Ithaca",
    loc %in% c("freev", "free"), "Freeville", 
    loc == "chazy", "Chazy",
    loc == "willsboro", "Willsboro",
    loc %in% c("rsgeneva", "rn041_gen", "rn028_gen"), "Geneva",
    default = tools::toTitleCase(tolower(loc))
  )]
  
  cat("Hemp data prepared: ", nrow(data), " samples\n")
  cat("Protein range: ", round(range(data$crude_protein)), " g/kg\n")
  
  return(data)
}

prepare_preprocessing_key <- function(raw_key) {
  # Create preprocessing method lookup
  key <- copy(raw_key)
  
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
# PREPROCESSING FUNCTIONS
# =============================================================================

apply_preprocessing <- function(spectra_matrix, method = 1) {
  # Apply single preprocessing method to spectral data
  
  # Convert method name to numeric ID if needed
  if (is.character(method)) {
    method_map <- c(
      "Raw spectra" = 1,
      "First derivative" = 2, 
      "Savitzky-Golay smoothing" = 3,
      "Gap-segment derivative" = 4,
      "Standard normal variate (SNV)" = 5,
      "SNV + Savitzky-Golay" = 6, 
      "SNV + detrending" = 7,
      "Multiplicative scatter correction" = 8
    )
    method <- method_map[method]
  }
  
  # Apply preprocessing
  processed <- switch(method,
                      "1" = spectra_matrix,  # Raw
                      "2" = t(diff(t(spectra_matrix), diff = 1)),  # First derivative
                      "3" = prospectr::savitzkyGolay(spectra_matrix, m = 1, p = 3, w = 5),  # SG smoothing
                      "4" = prospectr::gapDer(spectra_matrix, m = 1, w = 11, s = 5),  # Gap-segment
                      "5" = prospectr::standardNormalVariate(spectra_matrix),  # SNV
                      "6" = {  # SNV + SG
                        sg <- prospectr::savitzkyGolay(spectra_matrix, m = 1, p = 3, w = 5)
                        prospectr::standardNormalVariate(sg)
                      },
                      "7" = {  # SNV + detrend
                        snv <- prospectr::standardNormalVariate(spectra_matrix)
                        prospectr::detrend(snv, wav = 1:ncol(spectra_matrix))
                      },
                      "8" = prospectr::msc(spectra_matrix),  # MSC
                      spectra_matrix  # Default to raw
  )
  
  return(processed)
}

apply_all_preprocessing <- function(spectra_matrix) {
  # Apply all 8 preprocessing methods at once
  
  methods <- list(
    method_1 = apply_preprocessing(spectra_matrix, 1),
    method_2 = apply_preprocessing(spectra_matrix, 2),
    method_3 = apply_preprocessing(spectra_matrix, 3),
    method_4 = apply_preprocessing(spectra_matrix, 4),
    method_5 = apply_preprocessing(spectra_matrix, 5),
    method_6 = apply_preprocessing(spectra_matrix, 6),
    method_7 = apply_preprocessing(spectra_matrix, 7),
    method_8 = apply_preprocessing(spectra_matrix, 8)
  )
  
  return(methods)
}

# =============================================================================
# MODELING FUNCTIONS
# =============================================================================

fit_pls_model <- function(x_train, y_train, x_test, y_test, max_components = 20) {
  # Fit single PLS model and make predictions
  
  # Ensure matrix format
  if (!is.matrix(x_train)) x_train <- as.matrix(x_train)
  if (!is.matrix(x_test)) x_test <- as.matrix(x_test)
  
  # Fit PLS with cross-validation
  train_data <- data.frame(y = y_train, x_train)
  pls_model <- plsr(y ~ ., data = train_data, ncomp = max_components, validation = "CV")
  
  # Find optimal components
  cv_results <- RMSEP(pls_model, estimate = "CV")
  optimal_comps <- which.min(cv_results$val[1, 1, -1])
  
  # Make predictions
  predictions <- predict(pls_model, newdata = data.frame(x_test), ncomp = optimal_comps)
  predictions <- as.numeric(predictions)
  
  # Calculate metrics
  rmse <- sqrt(mean((y_test - predictions)^2))
  rsq <- cor(y_test, predictions)^2
  rpd <- sd(y_test) / rmse
  rpiq <- (quantile(y_test, 0.75) - quantile(y_test, 0.25)) / rmse
  
  return(list(
    predictions = predictions,
    optimal_components = optimal_comps,
    rmse = rmse,
    rsq = rsq,
    rpd = rpd,
    rpiq = rpiq,
    observed = y_test
  ))
}

run_single_iteration <- function(hemp_data, preprocessing_method_id, iteration) {
  # Run single modeling iteration
  
  # Extract data
  spectral_cols <- grep("^x[0-9]+$", names(hemp_data), value = TRUE)
  spectra_matrix <- as.matrix(hemp_data[, ..spectral_cols])
  y <- hemp_data$crude_protein
  
  # Train/test split
  set.seed(iteration)
  train_indices <- createDataPartition(y, p = 0.8, list = FALSE, groups = 3)[,1]
  
  x_train <- spectra_matrix[train_indices, ]
  y_train <- y[train_indices]
  x_test <- spectra_matrix[-train_indices, ]
  y_test <- y[-train_indices]
  
  # Apply preprocessing
  x_train_proc <- apply_preprocessing(x_train, preprocessing_method_id)
  x_test_proc <- apply_preprocessing(x_test, preprocessing_method_id)
  
  # Fit model and get results
  results <- fit_pls_model(x_train_proc, y_train, x_test_proc, y_test)
  
  # Return standardized results
  data.table(
    iteration = iteration,
    preprocessing_method = preprocessing_method_id,
    optimal_components = results$optimal_components,
    rmse = results$rmse,
    rsq = results$rsq,
    rpd = results$rpd,
    rpiq = results$rpiq
  )
}

# =============================================================================
# ANALYSIS WORKFLOW FUNCTIONS
# =============================================================================

run_preprocessing_comparison <- function(hemp_data, n_iterations = 100) {
  # Compare all preprocessing methods
  
  cat("Starting preprocessing comparison with", n_iterations, "iterations\n")
  
  # Run all combinations
  all_results <- list()
  
  for (iteration in 1:n_iterations) {
    if (iteration %% 20 == 0) cat("Iteration", iteration, "\n")
    
    for (method_id in 1:8) {
      tryCatch({
        result <- run_single_iteration(hemp_data, method_id, iteration)
        all_results[[length(all_results) + 1]] <- result
      }, error = function(e) {
        cat("Error in iteration", iteration, "method", method_id, ":", e$message, "\n")
      })
    }
  }
  
  # Combine results
  metrics <- rbindlist(all_results)
  
  cat("Preprocessing comparison complete:", nrow(metrics), "results\n")
  return(metrics)
}

analyze_preprocessing_methods <- function(preprocessing_comparison, preproc_key) {
  # Analyze preprocessing comparison results
  
  metrics <- preprocessing_comparison
  
  # Calculate summary statistics
  summary_stats <- metrics[, .(
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    mean_rsq = mean(rsq, na.rm = TRUE),
    sd_rsq = sd(rsq, na.rm = TRUE),
    mean_rpd = mean(rpd, na.rm = TRUE),
    sd_rpd = sd(rpd, na.rm = TRUE),
    mean_rpiq = mean(rpiq, na.rm = TRUE),
    sd_rpiq = sd(rpiq, na.rm = TRUE),
    mean_components = mean(optimal_components, na.rm = TRUE),
    n_iterations = .N
  ), by = preprocessing_method]
  
  # Add method names
  summary_with_names <- merge(summary_stats, preproc_key, 
                              by.x = "preprocessing_method", by.y = "preproc_id")
  
  cat("Preprocessing analysis complete for", nrow(summary_with_names), "methods\n")
  
  return(list(
    summary = summary_with_names,
    raw_metrics = metrics
  ))
}

select_best_preprocessing_method <- function(analysis) {
  # Select best method based on RMSE
  
  summary_data <- analysis$summary
  best_row <- summary_data[which.min(mean_rmse)]
  best_method <- best_row$method_name
  
  cat("Best preprocessing method:", best_method, "\n")
  cat("RMSE:", round(best_row$mean_rmse, 2), "±", round(best_row$sd_rmse, 2), "g/kg\n")
  
  return(best_method)
}

run_final_modeling <- function(hemp_data, best_method, n_iterations = 10) {
  
  cat("Running final modeling with", best_method, "for", n_iterations, "iterations\n")
  
  # Extract spectral data and preprocessing
  spec_cols <- grep("^x[0-9]+", names(hemp_data), value = TRUE)
  spectral_data <- as.matrix(hemp_data[, ..spec_cols])
  protein_data <- hemp_data$crude_protein
  
  # Apply preprocessing method
  # (Your existing preprocessing code here)
  processed_spectra <- apply_preprocessing(spectral_data, best_method)
  
  # Initialize storage
  metrics <- data.table()
  predictions <- data.table()
  model_n_comp_statistics <- data.table()  # ADD THIS LINE
  
  # Configuration
  max_components <- 20  # ADD THIS LINE
  
  for (i in 1:n_iterations) {
    
    tryCatch({
      # Your existing train/test split code
      train_idx <- sample(nrow(hemp_data), size = floor(0.75 * nrow(hemp_data)))
      test_idx <- setdiff(1:nrow(hemp_data), train_idx)
      
      X_train <- processed_spectra[train_idx, ]
      y_train <- protein_data[train_idx]
      X_test <- processed_spectra[test_idx, ]
      y_test <- protein_data[test_idx]
      
      # ADD THIS ENTIRE SECTION: Test each number of components
      iteration_rmse <- numeric(max_components)
      
      for (ncomp in 1:max_components) {
        # Fit PLS model with ncomp components
        pls_model <- pls::plsr(y_train ~ X_train, ncomp = ncomp, validation = "none")
        
        # Predict on test set
        test_preds <- predict(pls_model, newdata = X_test, ncomp = ncomp)[,,1]
        
        # Calculate RMSE for this component number
        iteration_rmse[ncomp] <- sqrt(mean((y_test - test_preds)^2))
        
        # Store in component progression data (THIS IS THE KEY FOR FIGURE 2)
        model_n_comp_statistics <- rbind(model_n_comp_statistics,
                                         data.table(id = i, ncomp = ncomp, RMSE = iteration_rmse[ncomp]))
      }
      
      # Find optimal number of components (minimum RMSE)
      optimal_ncomp <- which.min(iteration_rmse)
      # END OF ADDED SECTION
      
      # Your existing final model fitting (modify to use optimal_ncomp)
      final_model <- pls::plsr(y_train ~ X_train, ncomp = optimal_ncomp, validation = "none")
      final_preds <- predict(final_model, newdata = X_test, ncomp = optimal_ncomp)[,,1]
      
      # Your existing metrics calculation
      final_rmse <- sqrt(mean((y_test - final_preds)^2))
      final_rsq <- cor(y_test, final_preds)^2
      final_rpd <- sd(y_test) / final_rmse
      final_rpiq <- (quantile(y_test, 0.75) - quantile(y_test, 0.25)) / final_rmse
      
      # Store metrics (modify to include optimal_components)
      iteration_metrics <- data.table(
        iteration = i,
        optimal_components = optimal_ncomp,  # ADD THIS LINE
        rmse = final_rmse,
        rsq = final_rsq,
        rpd = final_rpd,
        rpiq = final_rpiq
      )
      metrics <- rbind(metrics, iteration_metrics)
      
      # Your existing predictions storage
      iteration_predictions <- data.table(
        iteration = i,
        sample_id = test_idx,
        actual = y_test,
        predicted = final_preds
      )
      predictions <- rbind(predictions, iteration_predictions)
      
    }, error = function(e) {
      cat("Error in iteration", i, ":", e$message, "\n")
    })
  }
  
  cat("Final modeling complete:", nrow(metrics), "iterations\n")
  
  return(list(
    metrics = metrics,
    predictions = predictions,
    model_n_comp_statistics = model_n_comp_statistics,  # ADD THIS LINE
    method = best_method
  ))
}

analyze_final_model <- function(final_model_results) {
  cat("=== USING FRESH RPD-PRIMARY ANALYSIS ===\n")
  
  metrics <- final_model_results$metrics
  predictions <- final_model_results$predictions
  
  # DEBUG: Check data
  cat("Total models:", nrow(metrics), "\n")
  cat("RPD range:", paste(round(range(metrics$rpd, na.rm = TRUE), 2), collapse = " to "), "\n")
  cat("Mean RPD:", round(mean(metrics$rpd, na.rm = TRUE), 2), "\n")
  
  # Overall statistics
  overall_stats <- metrics[, .(
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    mean_rsq = mean(rsq, na.rm = TRUE),
    sd_rsq = sd(rsq, na.rm = TRUE),
    mean_rpd = mean(rpd, na.rm = TRUE),
    sd_rpd = sd(rpd, na.rm = TRUE),
    mean_rpiq = mean(rpiq, na.rm = TRUE),
    sd_rpiq = sd(rpiq, na.rm = TRUE),
    n_successful = .N
  )]
  
  # Classification using RPD ONLY (primary approach)
  # Create classification vectors
  excellent_mask <- metrics$rpd > 3.0
  good_mask <- metrics$rpd >= 2.0 & metrics$rpd <= 3.0
  fair_mask <- metrics$rpd >= 1.5 & metrics$rpd < 2.0
  poor_mask <- metrics$rpd < 1.5
  
  # Count each category
  excellent_count <- sum(excellent_mask, na.rm = TRUE)
  good_count <- sum(good_mask, na.rm = TRUE)
  fair_count <- sum(fair_mask, na.rm = TRUE)
  poor_count <- sum(poor_mask, na.rm = TRUE)
  
  # DEBUG: Show classification counts
  cat("Classification counts:\n")
  cat("- Excellent (RPD > 3.0):", excellent_count, "\n")
  cat("- Good (RPD 2.0-3.0):", good_count, "\n")
  cat("- Fair (RPD 1.4-2.0):", fair_count, "\n")
  cat("- Poor (RPD < 1.4):", poor_count, "\n")
  cat("- Total:", excellent_count + good_count + fair_count + poor_count, "\n")
  
  # Supporting metric analysis (for those in each RPD category)
  excellent_with_support <- sum(excellent_mask & metrics$rpiq > 4.1 & metrics$rsq > 0.9, na.rm = TRUE)
  good_with_support <- sum(good_mask & metrics$rpiq >= 2.3 & metrics$rsq >= 0.8, na.rm = TRUE)
  fair_with_support <- sum(fair_mask & metrics$rpiq >= 1.5 & metrics$rsq >= 0.5, na.rm = TRUE)
  
  # Calculate percentages
  total_models <- nrow(metrics)
  excellent_pct <- round(100 * excellent_count / total_models, 1)
  good_pct <- round(100 * good_count / total_models, 1)
  fair_pct <- round(100 * fair_count / total_models, 1)
  poor_pct <- round(100 * poor_count / total_models, 1)
  
  # Supporting metric concordance rates
  excellent_support_pct <- ifelse(excellent_count > 0, round(100 * excellent_with_support / excellent_count, 1), 0)
  good_support_pct <- ifelse(good_count > 0, round(100 * good_with_support / good_count, 1), 0)
  fair_support_pct <- ifelse(fair_count > 0, round(100 * fair_with_support / fair_count, 1), 0)
  
  # Summary metrics
  quantitative_capable <- excellent_pct + good_pct
  total_acceptable <- excellent_pct + good_pct + fair_pct
  
  # Create results structure
  performance_summary <- list(
    total_models = total_models,
    excellent_models = excellent_count,
    good_models = good_count,
    fair_models = fair_count,
    poor_models = poor_count,
    excellent_percent = excellent_pct,
    good_percent = good_pct,
    fair_percent = fair_pct,
    poor_percent = poor_pct,
    excellent_support_rate = excellent_support_pct,
    good_support_rate = good_support_pct,
    fair_support_rate = fair_support_pct,
    quantitative_capable = quantitative_capable,
    total_acceptable = total_acceptable
  )
  
  # Print results
  cat("\n=== FINAL RESULTS ===\n")
  cat("Performance metrics (mean ± SD):\n")
  cat("- RMSE:", round(overall_stats$mean_rmse, 2), "±", round(overall_stats$sd_rmse, 2), "g/kg\n")
  cat("- R²:", round(overall_stats$mean_rsq, 3), "±", round(overall_stats$sd_rsq, 3), "\n")
  cat("- RPD:", round(overall_stats$mean_rpd, 2), "±", round(overall_stats$sd_rpd, 2), "\n")
  cat("- RPIQ:", round(overall_stats$mean_rpiq, 2), "±", round(overall_stats$sd_rpiq, 2), "\n\n")
  
  cat("Model classification (RPD-primary):\n")
  cat("- Excellent (RPD > 3.0):", excellent_count, "(", excellent_pct, "%)\n")
  cat("- Good (RPD 2.0-3.0):", good_count, "(", good_pct, "%)\n")
  cat("- Fair (RPD 1.5-2.0):", fair_count, "(", fair_pct, "%)\n")
  cat("- Poor (RPD < 1.5):", poor_count, "(", poor_pct, "%)\n\n")
  
  cat("Supporting metric concordance:\n")
  cat("- Excellent with full support:", excellent_support_pct, "%\n")
  cat("- Good with full support:", good_support_pct, "%\n")
  cat("- Fair with full support:", fair_support_pct, "%\n\n")
  
  cat("Summary:\n")
  cat("- Quantitative capable:", quantitative_capable, "%\n")
  cat("- Total acceptable:", total_acceptable, "%\n")
  
  return(list(
    overall_stats = overall_stats,
    performance_summary = performance_summary,
    metrics = metrics,
    predictions = predictions,
    classification_method = "rpd_primary_fresh"
  ))
}

analyze_prediction_errors <- function(final_model_results, hemp_data) {
  cat("=== ANALYZING REAL PREDICTION ERRORS ===\n")
  
  # Check that we have the predictions component
  if (!"predictions" %in% names(final_model_results)) {
    cat("❌ No predictions component found in final_model_results\n")
    cat("Available components:", names(final_model_results), "\n")
    return(list(
      sample_errors = data.table::data.table(),
      systematic_bias = "No predictions data available",
      raw_predictions = data.table::data.table()
    ))
  }
  
  predictions_data <- data.table::copy(final_model_results$predictions)
  cat("Found predictions data with", nrow(predictions_data), "observations\n")
  
  # DEBUG: Show actual column names
  cat("Actual columns in predictions data:", paste(names(predictions_data), collapse = ", "), "\n")
  
  # FLEXIBLE COLUMN DETECTION: Handle both sample_id and ith_in_data_set
  sample_id_col <- NULL
  if ("ith_in_data_set" %in% names(predictions_data)) {
    sample_id_col <- "ith_in_data_set"
  } else if ("sample_id" %in% names(predictions_data)) {
    sample_id_col <- "sample_id"
  } else {
    # Look for any column that might be sample ID
    possible_cols <- names(predictions_data)[grepl("sample|id", names(predictions_data), ignore.case = TRUE)]
    if (length(possible_cols) > 0) {
      sample_id_col <- possible_cols[1]
      cat("Using", sample_id_col, "as sample identifier\n")
    }
  }
  
  # Check for required columns with flexible sample ID
  required_cols <- c(sample_id_col, "actual", "predicted", "iteration")
  missing_cols <- setdiff(required_cols, names(predictions_data))
  
  if (length(missing_cols) > 0 || is.null(sample_id_col)) {
    cat("❌ Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
    if (is.null(sample_id_col)) {
      cat("❌ No sample ID column found. Looking for: ith_in_data_set, sample_id\n")
    }
    return(list(
      sample_errors = data.table::data.table(),
      systematic_bias = "Required columns missing from predictions data",
      raw_predictions = data.table::data.table()
    ))
  }
  
  # Standardize column name to ith_in_data_set for consistency
  if (sample_id_col != "ith_in_data_set") {
    data.table::setnames(predictions_data, sample_id_col, "ith_in_data_set")
  }
  
  cat("Covering", length(unique(predictions_data$ith_in_data_set)), "unique samples\n")
  cat("Across", length(unique(predictions_data$iteration)), "iterations\n")
  
  # Calculate additional error metrics if not present
  if (!"error_raw" %in% names(predictions_data)) {
    predictions_data[, error_raw := predicted - actual]
  }
  if (!"error_pct" %in% names(predictions_data)) {
    predictions_data[, error_pct := (predicted - actual) / actual * 100]
  }
  if (!"abs_residual" %in% names(predictions_data)) {
    predictions_data[, abs_residual := abs(error_raw)]
  }
  
  # Remove any missing values
  clean_data <- predictions_data[!is.na(actual) & !is.na(predicted)]
  cat("Clean observations after removing missing values:", nrow(clean_data), "\n")
  
  if (nrow(clean_data) == 0) {
    cat("❌ No valid observations after cleaning\n")
    return(list(
      sample_errors = data.table::data.table(),
      systematic_bias = "No valid observations after data cleaning",
      raw_predictions = data.table::data.table()
    ))
  }
  
  # Create sample-level ordering using TRUE sample identity
  sample_order <- clean_data[, .(
    mean_actual = mean(actual, na.rm = TRUE)
  ), by = ith_in_data_set][order(mean_actual)]
  sample_order[, plot_order := 1:.N]
  
  cat("Number of unique samples:", nrow(sample_order), "\n")
  
  # Add plot_order to all predictions using true sample identity
  clean_data <- merge(clean_data, sample_order[, .(ith_in_data_set, plot_order)], by = "ith_in_data_set")
  
  # Create tertiles based on actual CP values for tertile analysis
  clean_data[, cp_tertile := cut(actual, 3, labels = c("Low", "Medium", "High"))]
  
  cat("Raw predictions now has plot_order and error_raw columns\n")
  cat("Columns:", names(clean_data), "\n")
  
  # Calculate sample-level means using TRUE sample identity
  error_by_sample <- clean_data[, .(
    mean_error = mean(error_raw, na.rm = TRUE),
    sd_error = sd(error_raw, na.rm = TRUE),
    actual_cp = mean(actual, na.rm = TRUE),
    n_predictions = .N,
    plot_order = plot_order[1]  # Should be same for all rows of same sample
  ), by = ith_in_data_set]
  
  error_by_sample[, cp_tertile := cut(actual_cp, 3, labels = c("Low", "Medium", "High"))]
  error_by_sample <- error_by_sample[order(actual_cp)]
  error_by_sample[, rank_order := 1:.N]
  
  cat("Number of samples in error_by_sample:", nrow(error_by_sample), "\n")
  
  # Calculate bias by tertile
  tertile_means <- clean_data[, .(
    mean_error = mean(error_raw, na.rm = TRUE),
    mean_actual = mean(actual, na.rm = TRUE),
    n_obs = .N
  ), by = cp_tertile]
  
  cat("Tertile bias analysis:\n")
  print(tertile_means)
  
  # Calculate systematic bias using linear model
  lm_bias <- lm(mean_error ~ actual_cp, data = error_by_sample)
  bias_summary <- list(
    slope = coef(lm_bias)[2],
    intercept = coef(lm_bias)[1],
    r_squared = summary(lm_bias)$r.squared,
    tertile_biases = tertile_means
  )
  
  cat("Systematic bias analysis complete\n")
  cat("- Linear slope:", round(bias_summary$slope, 4), "\n")
  cat("- R-squared:", round(bias_summary$r_squared, 3), "\n")
  
  # Return comprehensive error analysis
  return(list(
    raw_predictions = clean_data,
    sample_errors = error_by_sample,
    systematic_bias = bias_summary,
    tertile_analysis = tertile_means,
    summary = list(
      total_predictions = nrow(clean_data),
      unique_samples = length(unique(clean_data$ith_in_data_set)),
      iterations = length(unique(clean_data$iteration)),
      mean_absolute_error = mean(abs(clean_data$error_raw)),
      rmse = sqrt(mean(clean_data$error_raw^2))
    )
  ))
}

# =============================================================================
# MULTI-ALGORITHM COMPARISON FUNCTIONS
# =============================================================================

run_multi_algorithm_comparison <- function(data, best_method, n_iterations = 100, 
                                           algorithms = c("pls", "svmRadial", "rf")) {
  
  cat("=== MULTI-ALGORITHM COMPARISON ===\n")
  cat("Algorithms:", paste(algorithms, collapse = ", "), "\n")
  cat("Preprocessing method:", best_method, "\n")
  cat("Iterations:", n_iterations, "\n\n")
  
  # Extract spectral data
  spectral_cols <- grep("^x[0-9]+$", names(data), value = TRUE)
  spectra_matrix <- as.matrix(data[, ..spectral_cols])
  y <- data$crude_protein
  
  # Convert method name to numeric if needed
  if (is.character(best_method)) {
    method_map <- c(
      "Raw spectra" = 1,
      "First derivative" = 2, 
      "Savitzky-Golay smoothing" = 3,
      "Gap-segment derivative" = 4,
      "Standard normal variate (SNV)" = 5,
      "SNV + Savitzky-Golay" = 6, 
      "SNV + detrending" = 7,
      "Multiplicative scatter correction" = 8
    )
    method_id <- method_map[best_method]
    if (is.na(method_id)) method_id <- 1  # Default to raw
  } else {
    method_id <- best_method
  }
  
  # Storage for results
  all_results <- list()
  
  for (i in 1:n_iterations) {
    if (i %% 25 == 0) cat("Iteration", i, "\n")
    
    tryCatch({
      # Same train/test split for fair comparison
      set.seed(i)
      train_indices <- createDataPartition(y, p = 0.75, list = FALSE, groups = 3)[,1]
      
      x_train <- spectra_matrix[train_indices, ]
      y_train <- y[train_indices]
      x_test <- spectra_matrix[-train_indices, ]
      y_test <- y[-train_indices]
      
      # Apply preprocessing
      x_train_proc <- apply_preprocessing(x_train, method_id)
      x_test_proc <- apply_preprocessing(x_test, method_id)
      
      # Convert to data frames
      train_data <- data.frame(y = y_train, x_train_proc)
      test_data <- data.frame(x_test_proc)
      
      # Train control
      train_control <- trainControl(method = "cv", number = 5, savePredictions = "final")
      
      # Test each algorithm
      for (algo in algorithms) {
        algo_result <- tryCatch({
          
          # Fit model
          if (algo == "pls") {
            model <- train(y ~ ., data = train_data, method = "pls", 
                           trControl = train_control, tuneLength = 15)
          } else if (algo == "svmRadial") {
            model <- train(y ~ ., data = train_data, method = "svmRadial", 
                           trControl = train_control, tuneLength = 8)
          } else if (algo == "rf") {
            model <- train(y ~ ., data = train_data, method = "rf", 
                           trControl = train_control, tuneLength = 5, ntree = 200)
          }
          
          # Make predictions
          predictions <- predict(model, test_data)
          
          # Calculate metrics
          rmse <- sqrt(mean((y_test - predictions)^2))
          rsq <- cor(y_test, predictions)^2
          rpd <- sd(y_test) / rmse
          rpiq <- (quantile(y_test, 0.75) - quantile(y_test, 0.25)) / rmse
          
          data.table(
            iteration = i,
            algorithm = algo,
            rmse = rmse,
            rsq = rsq,
            rpd = rpd,
            rpiq = rpiq
          )
          
        }, error = function(e) {
          cat("Error with", algo, "in iteration", i, ":", e$message, "\n")
          NULL
        })
        
        if (!is.null(algo_result)) {
          all_results[[length(all_results) + 1]] <- algo_result
        }
      }
      
    }, error = function(e) {
      cat("Error in iteration", i, ":", e$message, "\n")
    })
  }
  
  # Combine results
  if (length(all_results) > 0) {
    combined_results <- rbindlist(all_results)
    cat("Multi-algorithm comparison complete:", nrow(combined_results), "results\n")
    return(combined_results)
  } else {
    stop("No successful results from multi-algorithm comparison")
  }
}

analyze_multi_algorithm_results <- function(multi_algorithm_results) {
  # Analyze multi-algorithm comparison results using RPD-primary unified evaluation criteria
  # RPD determines the category, with optional reporting of RPIQ and R² support
  
  # Summary statistics by algorithm
  summary_stats <- multi_algorithm_results[, .(
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    mean_rsq = mean(rsq, na.rm = TRUE),
    sd_rsq = sd(rsq, na.rm = TRUE),
    mean_rpd = mean(rpd, na.rm = TRUE),
    sd_rpd = sd(rpd, na.rm = TRUE),
    mean_rpiq = mean(rpiq, na.rm = TRUE),
    sd_rpiq = sd(rpiq, na.rm = TRUE),
    n_iterations = .N
  ), by = algorithm]
  
  # Model quality assessment using RPD as primary classifier
  model_quality <- multi_algorithm_results[, {
    
    # Primary classification by RPD thresholds
    excellent_rpd <- rpd > 3.0
    good_rpd <- rpd >= 2.0 & rpd <= 3.0
    fair_rpd <- rpd >= 1.5 & rpd < 2.0
    poor_rpd <- rpd < 1.5
    
    # Supporting metrics analysis (for reporting, not classification)
    excellent_rpiq_support <- rpd > 3.0 & rpiq > 4.1
    excellent_rsq_support <- rpd > 3.0 & rsq > 0.9
    excellent_full_support <- rpd > 3.0 & rpiq > 4.1 & rsq > 0.9
    
    good_rpiq_support <- rpd >= 2.0 & rpd <= 3.0 & rpiq >= 2.3 & rpiq <= 4.1
    good_rsq_support <- rpd >= 2.0 & rpd <= 3.0 & rsq >= 0.8 & rsq <= 0.9
    good_full_support <- rpd >= 2.0 & rpd <= 3.0 & rpiq >= 2.3 & rpiq <= 4.1 & rsq >= 0.8 & rsq <= 0.9
    
    fair_rpiq_support <- rpd >= 1.5 & rpd < 2.0 & rpiq >= 1.5 & rpiq < 2.3
    fair_rsq_support <- rpd >= 1.5 & rpd < 2.0 & rsq >= 0.5 & rsq < 0.8
    fair_full_support <- rpd >= 1.5 & rpd < 2.0 & rpiq >= 1.5 & rpiq < 2.3 & rsq >= 0.5 & rsq < 0.8
    
    .(
      # Primary counts (based on RPD only)
      excellent_models = sum(excellent_rpd, na.rm = TRUE),
      good_models = sum(good_rpd, na.rm = TRUE),
      fair_models = sum(fair_rpd, na.rm = TRUE),
      poor_models = sum(poor_rpd, na.rm = TRUE),
      total_models = .N,
      
      # Supporting metric concordance (for quality assessment)
      excellent_with_rpiq = sum(excellent_rpiq_support, na.rm = TRUE),
      excellent_with_rsq = sum(excellent_rsq_support, na.rm = TRUE),
      excellent_with_all = sum(excellent_full_support, na.rm = TRUE),
      
      good_with_rpiq = sum(good_rpiq_support, na.rm = TRUE),
      good_with_rsq = sum(good_rsq_support, na.rm = TRUE),
      good_with_all = sum(good_full_support, na.rm = TRUE),
      
      fair_with_rpiq = sum(fair_rpiq_support, na.rm = TRUE),
      fair_with_rsq = sum(fair_rsq_support, na.rm = TRUE),
      fair_with_all = sum(fair_full_support, na.rm = TRUE),
      
      # Overall quality indicators
      high_quality_models = sum(rpd >= 2.0 & rpiq >= 2.3 & rsq >= 0.8, na.rm = TRUE),
      na_models = sum(is.na(rpd) | is.na(rpiq) | is.na(rsq))
    )
  }, by = algorithm]
  
  # Calculate percentages
  model_quality[, `:=`(
    excellent_pct = round(100 * excellent_models / total_models, 1),
    good_pct = round(100 * good_models / total_models, 1),
    fair_pct = round(100 * fair_models / total_models, 1),
    poor_pct = round(100 * poor_models / total_models, 1),
    
    # Supporting metric concordance rates
    excellent_rpiq_concordance = ifelse(excellent_models > 0, round(100 * excellent_with_rpiq / excellent_models, 1), 0),
    excellent_rsq_concordance = ifelse(excellent_models > 0, round(100 * excellent_with_rsq / excellent_models, 1), 0),
    excellent_full_concordance = ifelse(excellent_models > 0, round(100 * excellent_with_all / excellent_models, 1), 0),
    
    good_rpiq_concordance = ifelse(good_models > 0, round(100 * good_with_rpiq / good_models, 1), 0),
    good_rsq_concordance = ifelse(good_models > 0, round(100 * good_with_rsq / good_models, 1), 0),
    good_full_concordance = ifelse(good_models > 0, round(100 * good_with_all / good_models, 1), 0),
    
    # Summary categories
    quantitative_capable = round(100 * (excellent_models + good_models) / total_models, 1),  # Excellent + Good
    qualitative_capable = round(100 * fair_models / total_models, 1),  # Fair only
    total_acceptable = round(100 * (excellent_models + good_models + fair_models) / total_models, 1),  # All except Poor
    high_quality_pct = round(100 * high_quality_models / total_models, 1)  # Meet all criteria regardless of category
  )]
  
  cat("Multi-algorithm analysis complete using RPD-primary unified criteria for", length(unique(summary_stats$algorithm)), "algorithms\n")
  cat("Primary classification by RPD thresholds:\n")
  cat("- Excellent: RPD > 3.0 (excellent quantitative prediction)\n")
  cat("- Good: RPD 2.0-3.0 (good quantitative prediction)\n") 
  cat("- Fair: RPD 1.5-2.0 (qualitative/screening prediction)\n")
  cat("- Poor: RPD < 1.5 (unreliable)\n\n")
  
  # Print classification summary
  for (i in 1:nrow(model_quality)) {
    algo <- model_quality$algorithm[i]
    total <- model_quality$total_models[i]
    excellent <- model_quality$excellent_models[i]
    good <- model_quality$good_models[i]
    fair <- model_quality$fair_models[i]
    poor <- model_quality$poor_models[i]
    
    cat(sprintf("%s: %d total models\n", algo, total))
    cat(sprintf("  Excellent: %d (%.1f%%), Good: %d (%.1f%%), Fair: %d (%.1f%%), Poor: %d (%.1f%%)\n", 
                excellent, model_quality$excellent_pct[i],
                good, model_quality$good_pct[i], 
                fair, model_quality$fair_pct[i],
                poor, model_quality$poor_pct[i]))
    
    if (excellent > 0) {
      cat(sprintf("  Excellent models with RPIQ support: %.1f%%, R² support: %.1f%%, both: %.1f%%\n",
                  model_quality$excellent_rpiq_concordance[i],
                  model_quality$excellent_rsq_concordance[i],
                  model_quality$excellent_full_concordance[i]))
    }
    cat("\n")
  }
  
  return(list(
    summary_stats = summary_stats,
    model_quality = model_quality,
    raw_results = multi_algorithm_results,
    classification_method = "unified_criteria_rpd_primary"
  ))
}

# analyze_multi_algorithm_results <- function(multi_algorithm_results) {
#   # Analyze multi-algorithm comparison results using STRICT unified evaluation criteria
#   # A model must meet ALL three criteria (RPD AND RPIQ AND R²) to be classified in a category
#   
#   # Summary statistics by algorithm
#   summary_stats <- multi_algorithm_results[, .(
#     mean_rmse = mean(rmse, na.rm = TRUE),
#     sd_rmse = sd(rmse, na.rm = TRUE),
#     mean_rsq = mean(rsq, na.rm = TRUE),
#     sd_rsq = sd(rsq, na.rm = TRUE),
#     mean_rpd = mean(rpd, na.rm = TRUE),
#     sd_rpd = sd(rpd, na.rm = TRUE),
#     mean_rpiq = mean(rpiq, na.rm = TRUE),
#     sd_rpiq = sd(rpiq, na.rm = TRUE),
#     n_iterations = .N
#   ), by = algorithm]
#   
#   # Strict model quality assessment - ALL criteria must be met
#   # Categories are mutually exclusive and hierarchical
#   model_quality <- multi_algorithm_results[, {
#     
#     # Apply strict hierarchical classification
#     # Excellent: RPD > 3.0 AND RPIQ > 4.1 AND R² > 0.9
#     excellent_mask <- rpd > 3.0 & rpiq > 4.1 & rsq > 0.9
#     
#     # Good: RPD 2.0-3.0 AND RPIQ 2.3-4.1 AND R² 0.8-0.9 (and not excellent)
#     good_mask <- !excellent_mask & 
#       rpd >= 2.0 & rpd <= 3.0 & 
#       rpiq >= 2.3 & rpiq <= 4.1 & 
#       rsq >= 0.8 & rsq <= 0.9
#     
#     # Fair: RPD 1.4-2.0 AND RPIQ 1.5-2.3 AND R² 0.5-0.8 (and not excellent or good)
#     fair_mask <- !excellent_mask & !good_mask &
#       rpd >= 1.4 & rpd < 2.0 & 
#       rpiq >= 1.5 & rpiq < 2.3 & 
#       rsq >= 0.5 & rsq < 0.8
#     
#     # Poor: Everything else (fails to meet any complete category)
#     poor_mask <- !excellent_mask & !good_mask & !fair_mask
#     
#     # Count models in each category
#     .(
#       excellent_models = sum(excellent_mask, na.rm = TRUE),
#       good_models = sum(good_mask, na.rm = TRUE),
#       fair_models = sum(fair_mask, na.rm = TRUE),
#       poor_models = sum(poor_mask, na.rm = TRUE),
#       total_models = .N,
#       
#       # Additional diagnostics
#       na_models = sum(is.na(rpd) | is.na(rpiq) | is.na(rsq))
#     )
#   }, by = algorithm]
#   
#   # Calculate percentages
#   model_quality[, `:=`(
#     excellent_pct = round(100 * excellent_models / total_models, 1),
#     good_pct = round(100 * good_models / total_models, 1),
#     fair_pct = round(100 * fair_models / total_models, 1),
#     poor_pct = round(100 * poor_models / total_models, 1),
#     
#     # Summary categories
#     quantitative_capable = round(100 * (excellent_models + good_models) / total_models, 1),  # Excellent + Good
#     qualitative_capable = round(100 * fair_models / total_models, 1),  # Fair only  
#     total_acceptable = round(100 * (excellent_models + good_models + fair_models) / total_models, 1)  # All except Poor
#   )]
#   
#   # Verification: check that all models are classified
#   model_quality[, total_classified := excellent_models + good_models + fair_models + poor_models]
#   
#   cat("Multi-algorithm analysis complete using STRICT unified criteria for", length(unique(summary_stats$algorithm)), "algorithms\n")
#   cat("Classification requires ALL three metrics to meet thresholds:\n")
#   cat("- Excellent: RPD > 3.0 AND RPIQ > 4.1 AND R² > 0.9\n")
#   cat("- Good: RPD 2.0-3.0 AND RPIQ 2.3-4.1 AND R² 0.8-0.9\n") 
#   cat("- Fair: RPD 1.4-2.0 AND RPIQ 1.5-2.3 AND R² 0.5-0.8\n")
#   cat("- Poor: Fails to meet any complete category criteria\n\n")
#   
#   # Print verification
#   for (i in 1:nrow(model_quality)) {
#     algo <- model_quality$algorithm[i]
#     total <- model_quality$total_models[i]
#     classified <- model_quality$total_classified[i]
#     cat(sprintf("%s: %d/%d models classified", algo, classified, total), "\n")
#   }
#   
#   return(list(
#     summary_stats = summary_stats,
#     model_quality = model_quality,
#     raw_results = multi_algorithm_results,
#     classification_method = "unified_criteria_strict"
#   ))
# }

# =============================================================================
# SPECTRAL ANALYSIS FUNCTIONS
# =============================================================================

extract_wavelengths <- function(data) {
  # Extract wavelength information from column names
  spectral_cols <- grep("^x[0-9]+$", names(data), value = TRUE)
  wavelengths <- as.numeric(gsub("^x", "", spectral_cols))
  
  list(
    column_names = spectral_cols,
    wavelengths = wavelengths,
    n_wavelengths = length(wavelengths),
    range = range(wavelengths)
  )
}

select_protein_wavelengths <- function(data) {
  # Select protein-specific wavelengths
  
  protein_ranges <- list(
    list(start = 1200, end = 1250, name = "C-H stretch 2nd overtone (amino acids)"),
    list(start = 1500, end = 1550, name = "N-H stretch 1st overtone (peptide bonds)"),
    # list(start = 1660, end = 1700, name = "C-H stretch 1st overtone (aliphatic amino acids)"),
    list(start = 2040, end = 2090, name = "N-H + C-N combination (protein backbone)")
    # list(start = 2270, end = 2310, name = "C-H + C-H combination (amino acid structure)")
  )
  
  wavelength_info <- extract_wavelengths(data)
  all_wavelengths <- wavelength_info$wavelengths
  
  # Select wavelengths in protein bands
  protein_wavelengths <- c()
  protein_assignments <- c()
  
  for (band in protein_ranges) {
    in_band <- all_wavelengths >= band$start & all_wavelengths <= band$end
    band_wavelengths <- all_wavelengths[in_band]
    
    if (length(band_wavelengths) > 0) {
      protein_wavelengths <- c(protein_wavelengths, band_wavelengths)
      protein_assignments <- c(protein_assignments, rep(band$name, length(band_wavelengths)))
    }
  }
  
  selected_columns <- paste0("x", protein_wavelengths)
  available_columns <- selected_columns[selected_columns %in% names(data)]
  
  cat("Selected", length(protein_wavelengths), "protein-specific wavelengths\n")
  cat("Available in data:", length(available_columns), "wavelengths\n")
  
  list(
    wavelengths = protein_wavelengths[selected_columns %in% names(data)],
    column_names = available_columns,
    assignments = protein_assignments[selected_columns %in% names(data)],
    n_selected = length(available_columns)
  )
}

fit_pls_for_spectral_analysis <- function(data, best_method) {
  # Fit PLS model for spectral analysis
  
  # Extract spectral data
  spectral_cols <- grep("^x[0-9]+$", names(data), value = TRUE)
  spectra_matrix <- as.matrix(data[, ..spectral_cols])
  y <- data$crude_protein
  
  # Create train/test split
  set.seed(123)
  train_indices <- createDataPartition(y, p = 0.75, list = FALSE, groups = 3)[,1]
  
  x_train <- spectra_matrix[train_indices, ]
  y_train <- y[train_indices]
  
  # Apply preprocessing
  x_train_proc <- apply_preprocessing(x_train, best_method)
  
  # Convert to data frame
  train_data <- data.frame(y = y_train, x_train_proc)
  
  # Train PLS model
  train_control <- trainControl(method = "cv", number = 10)
  model <- train(y ~ ., data = train_data, method = "pls", 
                 trControl = train_control, tuneLength = 20)
  
  cat("PLS model fitted with", model$bestTune$ncomp, "components\n")
  
  # Get wavelength info
  wavelength_info <- extract_wavelengths(data)
  
  list(
    model = model,
    wavelengths = wavelength_info$wavelengths,
    X_train = train_data[, -1],  # Remove y column
    y_train = y_train,
    train_indices = train_indices,
    preprocessing_method = best_method
  )
}

extract_pls_coefficients <- function(spectral_fit) {
  # Extract PLS coefficients
  
  model <- spectral_fit$model
  coefficients <- coef(model$finalModel, ncomp = model$bestTune$ncomp)
  coefficients <- as.numeric(coefficients)
  
  # Get wavelengths
  wavelengths <- spectral_fit$wavelengths
  n_coeff <- length(coefficients)
  n_wavelengths <- length(wavelengths)
  
  # Handle length mismatch (common with derivatives)
  if (n_coeff == n_wavelengths) {
    used_wavelengths <- wavelengths
  } else {
    used_wavelengths <- wavelengths[1:n_coeff]
  }
  
  data.frame(
    wavelength = used_wavelengths,
    coefficient = coefficients,
    abs_coefficient = abs(coefficients)
  )
}

calculate_vip_scores <- function(spectral_fit) {
  # Calculate Variable Importance in Projection scores
  
  model <- spectral_fit$model
  pls_object <- model$finalModel
  optimal_ncomp <- model$bestTune$ncomp
  
  # Get PLS components
  W <- pls_object$loading.weights[, 1:optimal_ncomp, drop = FALSE]
  T <- pls_object$scores[, 1:optimal_ncomp, drop = FALSE]
  
  # Calculate VIP scores
  SS <- colSums(T^2)
  SS_total <- sum(SS)
  p <- nrow(W)
  vip_scores <- sqrt(p * rowSums((W^2) %*% diag(SS / SS_total)) / optimal_ncomp)
  
  # Align with wavelengths
  wavelengths <- spectral_fit$wavelengths
  n_vip <- length(vip_scores)
  n_wavelengths <- length(wavelengths)
  
  if (n_vip == n_wavelengths) {
    used_wavelengths <- wavelengths
  } else {
    used_wavelengths <- wavelengths[1:n_vip]
  }
  
  data.frame(
    wavelength = used_wavelengths,
    vip_score = vip_scores,
    important = vip_scores > 1.0
  )
}

run_spectral_analysis <- function(hemp_data, best_method, n_iterations = 1000) {
  # Run full spectrum validation with VIP calculation for each iteration
  
  cat("Starting full spectrum validation with VIP calculation for", n_iterations, "iterations\n")
  
  # Extract spectral data
  spectral_cols <- grep("^x[0-9]+$", names(hemp_data), value = TRUE)
  spectra_matrix <- as.matrix(hemp_data[, ..spectral_cols])
  y <- hemp_data$crude_protein
  
  # Convert method to numeric ID
  method_id <- ifelse(is.character(best_method), 
                      c("Raw spectra"=1, "First derivative"=2, "Savitzky-Golay smoothing"=3, 
                        "Gap-segment derivative"=4, "Standard normal variate (SNV)"=5, 
                        "SNV + Savitzky-Golay"=6, "SNV + detrending"=7, 
                        "Multiplicative scatter correction"=8)[best_method], 
                      best_method)
  
  # Get wavelengths
  wavelengths <- as.numeric(gsub("^x", "", spectral_cols))
  
  metrics <- data.table()
  vip_scores_all <- data.table()
  
  for (i in 1:n_iterations) {
    if (i %% 100 == 0) cat("Iteration", i, "\n")
    
    tryCatch({
      # SAME split as protein-focused analysis
      set.seed(i)
      train_idx <- sample(nrow(hemp_data), size = floor(0.75 * nrow(hemp_data)))
      test_idx <- setdiff(1:nrow(hemp_data), train_idx)
      
      x_train <- spectra_matrix[train_idx, ]
      y_train <- y[train_idx]
      x_test <- spectra_matrix[test_idx, ]
      y_test <- y[test_idx]
      
      # Apply preprocessing and fit model
      x_train_proc <- apply_preprocessing(x_train, method_id)
      x_test_proc <- apply_preprocessing(x_test, method_id)
      
      train_data <- data.frame(y = y_train, x_train_proc)
      model <- train(y ~ ., data = train_data, method = "pls", 
                     trControl = trainControl(method = "cv", number = 10), 
                     tuneLength = 15)
      
      final_preds <- predict(model, data.frame(x_test_proc))
      
      # Calculate metrics
      final_rmse <- sqrt(mean((y_test - final_preds)^2))
      final_rsq <- cor(y_test, final_preds)^2
      final_rpd <- sd(y_test) / final_rmse
      final_rpiq <- (quantile(y_test, 0.75) - quantile(y_test, 0.25)) / final_rmse
      
      metrics <- rbind(metrics, data.table(
        iteration = i, optimal_components = model$bestTune$ncomp,
        rmse = final_rmse, rsq = final_rsq, rpd = final_rpd, rpiq = final_rpiq
      ))
      
      # Calculate VIP scores
      pls_object <- model$finalModel
      optimal_ncomp <- model$bestTune$ncomp
      
      W <- pls_object$loading.weights[, 1:optimal_ncomp, drop = FALSE]
      T <- pls_object$scores[, 1:optimal_ncomp, drop = FALSE]
      
      SS <- colSums(T^2)
      SS_total <- sum(SS)
      p <- nrow(W)
      vip_scores <- sqrt(p * rowSums((W^2) %*% diag(SS / SS_total)) / optimal_ncomp)
      
      # Handle length mismatch (common with derivatives)
      n_vars <- length(vip_scores)
      used_wavelengths <- wavelengths[1:n_vars]
      
      # Store VIP scores
      iteration_vip <- data.table(
        iteration = i,
        wavelength = used_wavelengths,
        vip_score = vip_scores,
        important = vip_scores > 1.0
      )
      vip_scores_all <- rbind(vip_scores_all, iteration_vip)
      
    }, error = function(e) cat("Error in iteration", i, ":", e$message, "\n"))
  }
  
  return(list(
    metrics = metrics,
    vip_scores = vip_scores_all,
    method = best_method,
    n_wavelengths = length(wavelengths)
  ))
}

run_protein_focused_analysis <- function(data, best_method) {
  # Run protein-focused spectral analysis
  
  cat("Running protein-focused spectral analysis...\n")
  
  # Select protein wavelengths
  protein_selection <- select_protein_wavelengths(data)
  
  # Create subset data with only protein wavelengths
  protein_data <- data[, c("crude_protein", protein_selection$column_names), with = FALSE]
  
  # Fit PLS model on protein wavelengths only
  spectral_cols <- protein_selection$column_names
  spectra_matrix <- as.matrix(protein_data[, ..spectral_cols])
  y <- protein_data$crude_protein
  
  # Train/test split
  set.seed(123)
  train_indices <- createDataPartition(y, p = 0.75, list = FALSE, groups = 3)[,1]
  
  x_train <- spectra_matrix[train_indices, ]
  y_train <- y[train_indices]
  
  # Apply preprocessing
  x_train_proc <- apply_preprocessing(x_train, best_method)
  
  # Train model
  train_data <- data.frame(y = y_train, x_train_proc)
  train_control <- trainControl(method = "cv", number = 10)
  model <- train(y ~ ., data = train_data, method = "pls", 
                 trControl = train_control, tuneLength = 15)
  
  # Create spectral fit object
  spectral_fit <- list(
    model = model,
    wavelengths = protein_selection$wavelengths,
    X_train = train_data[, -1],
    y_train = y_train,
    train_indices = train_indices,
    preprocessing_method = best_method
  )
  
  # Extract analysis data
  coeff_data <- extract_pls_coefficients(spectral_fit)
  vip_data <- calculate_vip_scores(spectral_fit)
  
  n_important_vip <- sum(vip_data$vip_score > 1.0)
  
  cat("Protein-focused analysis complete:\n")
  cat("- Protein wavelengths used:", length(protein_selection$wavelengths), "\n")
  cat("- Important wavelengths (VIP > 1.0):", n_important_vip, "\n")
  cat("- Model components:", model$bestTune$ncomp, "\n")
  
  list(
    spectral_fit = spectral_fit,
    coefficients = coeff_data,
    vip_scores = vip_data,
    protein_selection = protein_selection,
    analysis_stats = list(
      n_important_vip = n_important_vip,
      optimal_components = model$bestTune$ncomp,
      n_protein_wavelengths = length(protein_selection$wavelengths)
    )
  )
}

compare_full_vs_protein_models <- function(spectral_analysis, protein_focused_analysis) {
  # Compare full spectrum vs protein-focused models
  
  cat("Comparing full spectrum vs protein-focused models...\n")
  
  # Extract model performance
  full_model <- spectral_analysis$spectral_fit$model
  protein_model <- protein_focused_analysis$spectral_fit$model
  
  # Get best performance metrics
  full_rmse <- min(full_model$results$RMSE)
  full_rsq <- max(full_model$results$Rsquared)
  protein_rmse <- min(protein_model$results$RMSE)
  protein_rsq <- max(protein_model$results$Rsquared)
  
  # Create comparison table
  comparison_data <- data.frame(
    Model = c("Full Spectrum", "Protein-Focused"),
    Wavelengths = c(length(spectral_analysis$spectral_fit$wavelengths),
                    length(protein_focused_analysis$protein_selection$wavelengths)),
    Components = c(full_model$bestTune$ncomp, protein_model$bestTune$ncomp),
    RMSE = c(full_rmse, protein_rmse),
    R_squared = c(full_rsq, protein_rsq),
    Efficiency = c("Baseline", 
                   paste0(round(100 * length(protein_focused_analysis$protein_selection$wavelengths) / 
                                  length(spectral_analysis$spectral_fit$wavelengths), 1), "% of variables"))
  )
  
  cat("Model comparison complete:\n")
  cat("- Full spectrum RMSE:", round(full_rmse, 3), "\n")
  cat("- Protein-focused RMSE:", round(protein_rmse, 3), "\n")
  cat("- Variable reduction:", round(100 * (1 - nrow(comparison_data[2,]) / nrow(comparison_data[1,]))), "%\n")
  
  list(
    comparison_table = comparison_data,
    full_analysis = spectral_analysis,
    protein_analysis = protein_focused_analysis
  )
}

# =============================================================================
# TABLE AND VISUALIZATION FUNCTIONS
# =============================================================================

create_sample_summary_table <- function(hemp_data) {
  # Create sample summary table
  
  protein_summary <- hemp_data[, .(
    Mean = round(mean(crude_protein, na.rm = TRUE), 0),
    SD = round(sd(crude_protein, na.rm = TRUE), 0),
    Minimum = round(min(crude_protein, na.rm = TRUE), 0),
    Q1 = round(quantile(crude_protein, 0.25, na.rm = TRUE), 0),
    Median = round(median(crude_protein, na.rm = TRUE), 0),
    Q3 = round(quantile(crude_protein, 0.75, na.rm = TRUE), 0),
    Maximum = round(max(crude_protein, na.rm = TRUE), 0)
  )]
  
  knitr::kable(
    protein_summary,
    caption = "Summary of Laboratory Assayed CP Values (g kg^-1^)",
    col.names = c("Mean", "SD", "Minimum", "First Quartile", "Median", "Third Quartile", "Maximum")
  )
}

create_preprocessing_comparison_table <- function(analysis) {
  # Create preprocessing comparison table
  
  table_data <- analysis$summary[, .(
    Method = method_name,
    RMSE = paste0(sprintf("%.1f", mean_rmse), " ± ", sprintf("%.1f", sd_rmse)),
    R_squared = paste0(sprintf("%.3f", mean_rsq), " ± ", sprintf("%.3f", sd_rsq)),
    RPD = paste0(sprintf("%.2f", mean_rpd), " ± ", sprintf("%.2f", sd_rpd)),
    RPIQ = paste0(sprintf("%.2f", mean_rpiq), " ± ", sprintf("%.2f", sd_rpiq))
  )]
  
  # Order by RMSE (best first)
  table_data <- table_data[order(as.numeric(gsub(" .*", "", RMSE)))]
  
  knitr::kable(
    table_data,
    caption = "Preprocessing Method Comparison (Mean ± SD)",
    col.names = c("Preprocessing Method", "RMSE (g/kg)", "R²", "RPD", "RPIQ")
  )
}

create_algorithm_comparison_table <- function(multi_algo_analysis) {
  # Create algorithm comparison table - DOCX compatible version
  
  summary_stats <- multi_algo_analysis$summary_stats
  model_quality <- multi_algo_analysis$model_quality
  
  # Combine data
  combined_table <- merge(summary_stats, model_quality, by = "algorithm")
  
  # Format for manuscript
  formatted_table <- combined_table %>%
    mutate(
      Algorithm = case_when(
        algorithm == "pls" ~ "Partial Least Squares",
        algorithm == "svmRadial" ~ "Support Vector Machine", 
        algorithm == "rf" ~ "Random Forest",
        TRUE ~ str_to_title(algorithm)
      ),
      RMSE_formatted = paste0(round(mean_rmse, 3), " (±", round(sd_rmse, 3), ")"),
      R2_formatted = paste0(round(mean_rsq, 3), " (±", round(sd_rsq, 3), ")"),
      RPD_formatted = paste0(round(mean_rpd, 2), " (±", round(sd_rpd, 2), ")"),
      RPIQ_formatted = paste0(round(mean_rpiq, 2), " (±", round(sd_rpiq, 2), ")"),
      
      # RPD-based categories
      Excellent_formatted = paste0(excellent_pct, "%"),
      Good_formatted = paste0(good_pct, "%"),
      Fair_formatted = paste0(fair_pct, "%"),
      Poor_formatted = paste0(poor_pct, "%"),
      
      # Summary categories  
      Quantitative_formatted = paste0(quantitative_capable, "%"),
      Qualitative_formatted = paste0(qualitative_capable, "%"),
      Total_formatted = paste0(total_acceptable, "%"),
      HighQuality_formatted = paste0(high_quality_pct, "%")
    )
  
  # Create final table
  final_table <- data.frame(
    "Algorithm" = formatted_table$Algorithm,
    "RMSE (±SD)" = formatted_table$RMSE_formatted,
    "R² (±SD)" = formatted_table$R2_formatted,
    "RPD (±SD)" = formatted_table$RPD_formatted,
    "RPIQ (±SD)" = formatted_table$RPIQ_formatted,
    "Excellent (%)" = formatted_table$Excellent_formatted,
    "Good (%)" = formatted_table$Good_formatted,
    "Fair (%)" = formatted_table$Fair_formatted,
    "Poor (%)" = formatted_table$Poor_formatted,
    "Quantitative (%)" = formatted_table$Quantitative_formatted,
    "Qualitative (%)" = formatted_table$Qualitative_formatted,
    "Total Usable (%)" = formatted_table$Total_formatted,
    check.names = FALSE
  )
  
  # DOCX-compatible table - no kableExtra styling
  knitr::kable(
    final_table,
    caption = "Algorithm performance using RPD-primary unified evaluation criteria. Primary classification by RPD: Excellent (>3.0), Good (2.0-3.0), Fair (1.4-2.0), Poor (<1.4). RPIQ and R² provide supporting evidence.",
    row.names = FALSE,
    align = c("l", rep("c", 11))
  )
}


# =============================================================================
# OPTIONAL: FUNCTION TO DIAGNOSE CLASSIFICATION MISMATCHES  
# =============================================================================

diagnose_model_classification <- function(multi_algorithm_results, algorithm_name = "pls") {
  # Diagnostic function to understand why models fall into different categories
  
  data_subset <- multi_algorithm_results[algorithm == algorithm_name]
  
  # Apply the same classification logic
  data_subset[, category := {
    excellent_mask <- rpd > 3.0 & rpiq > 4.1 & rsq > 0.9
    good_mask <- !excellent_mask & 
      rpd >= 2.0 & rpd <= 3.0 & 
      rpiq >= 2.3 & rpiq <= 4.1 & 
      rsq >= 0.8 & rsq <= 0.9
    fair_mask <- !excellent_mask & !good_mask &
      rpd >= 1.5 & rpd < 2.0 & 
      rpiq >= 1.5 & rpiq < 2.3 & 
      rsq >= 0.5 & rsq < 0.8
    
    case_when(
      excellent_mask ~ "Excellent",
      good_mask ~ "Good", 
      fair_mask ~ "Fair",
      TRUE ~ "Poor"
    )
  }]
  
  # Show breakdown
  cat("Classification breakdown for", algorithm_name, ":\n")
  print(table(data_subset$category))
  
  # Show examples of each category
  cat("\nSample models from each category:\n")
  for (cat_name in c("Excellent", "Good", "Fair", "Poor")) {
    sample_data <- data_subset[category == cat_name][1:min(3, .N)]
    if (nrow(sample_data) > 0) {
      cat("\n", cat_name, "examples:\n")
      print(sample_data[, .(iteration, rpd, rpiq, rsq)])
    }
  }
  
  return(data_subset)
}

# =============================================================================
# BACKWARDS COMPATIBILITY: UPDATED analyze_final_model FUNCTION
# =============================================================================

analyze_final_model <- function(final_model_results) {
  # Analyze final model performance using RPD-primary unified criteria
  # This function uses a completely different name to avoid conflicts
  
  cat("=== USING FRESH RPD-PRIMARY ANALYSIS ===\n")
  
  metrics <- final_model_results$metrics
  predictions <- final_model_results$predictions
  
  # DEBUG: Check data
  cat("Total models:", nrow(metrics), "\n")
  cat("RPD range:", paste(round(range(metrics$rpd, na.rm = TRUE), 2), collapse = " to "), "\n")
  cat("Mean RPD:", round(mean(metrics$rpd, na.rm = TRUE), 2), "\n")
  
  # Overall statistics
  overall_stats <- metrics[, .(
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    mean_rsq = mean(rsq, na.rm = TRUE),
    sd_rsq = sd(rsq, na.rm = TRUE),
    mean_rpd = mean(rpd, na.rm = TRUE),
    sd_rpd = sd(rpd, na.rm = TRUE),
    mean_rpiq = mean(rpiq, na.rm = TRUE),
    sd_rpiq = sd(rpiq, na.rm = TRUE),
    n_successful = .N
  )]
  
  # Classification using RPD ONLY (primary approach)
  # Create classification vectors
  excellent_mask <- metrics$rpd > 3.0
  good_mask <- metrics$rpd >= 2.0 & metrics$rpd <= 3.0
  fair_mask <- metrics$rpd >= 1.5 & metrics$rpd < 2.0
  poor_mask <- metrics$rpd < 1.5
  
  # Count each category
  excellent_count <- sum(excellent_mask, na.rm = TRUE)
  good_count <- sum(good_mask, na.rm = TRUE)
  fair_count <- sum(fair_mask, na.rm = TRUE)
  poor_count <- sum(poor_mask, na.rm = TRUE)
  
  # DEBUG: Show classification counts
  cat("Classification counts:\n")
  cat("- Excellent (RPD > 3.0):", excellent_count, "\n")
  cat("- Good (RPD 2.0-3.0):", good_count, "\n")
  cat("- Fair (RPD 1.5-2.0):", fair_count, "\n")
  cat("- Poor (RPD < 1.5):", poor_count, "\n")
  cat("- Total:", excellent_count + good_count + fair_count + poor_count, "\n")
  
  # Supporting metric analysis (for those in each RPD category)
  excellent_with_support <- sum(excellent_mask & metrics$rpiq > 4.1 & metrics$rsq > 0.9, na.rm = TRUE)
  good_with_support <- sum(good_mask & metrics$rpiq >= 2.3 & metrics$rsq >= 0.8, na.rm = TRUE)
  fair_with_support <- sum(fair_mask & metrics$rpiq >= 1.5 & metrics$rsq >= 0.5, na.rm = TRUE)
  
  # Calculate percentages
  total_models <- nrow(metrics)
  excellent_pct <- round(100 * excellent_count / total_models, 1)
  good_pct <- round(100 * good_count / total_models, 1)
  fair_pct <- round(100 * fair_count / total_models, 1)
  poor_pct <- round(100 * poor_count / total_models, 1)
  
  # Supporting metric concordance rates
  excellent_support_pct <- ifelse(excellent_count > 0, round(100 * excellent_with_support / excellent_count, 1), 0)
  good_support_pct <- ifelse(good_count > 0, round(100 * good_with_support / good_count, 1), 0)
  fair_support_pct <- ifelse(fair_count > 0, round(100 * fair_with_support / fair_count, 1), 0)
  
  # Summary metrics
  quantitative_capable <- excellent_pct + good_pct
  total_acceptable <- excellent_pct + good_pct + fair_pct
  
  # Create results structure
  performance_summary <- list(
    total_models = total_models,
    excellent_models = excellent_count,
    good_models = good_count,
    fair_models = fair_count,
    poor_models = poor_count,
    excellent_percent = excellent_pct,
    good_percent = good_pct,
    fair_percent = fair_pct,
    poor_percent = poor_pct,
    excellent_support_rate = excellent_support_pct,
    good_support_rate = good_support_pct,
    fair_support_rate = fair_support_pct,
    quantitative_capable = quantitative_capable,
    total_acceptable = total_acceptable
  )
  
  # Print results
  cat("\n=== FINAL RESULTS ===\n")
  cat("Performance metrics (mean ± SD):\n")
  cat("- RMSE:", round(overall_stats$mean_rmse, 2), "±", round(overall_stats$sd_rmse, 2), "g/kg\n")
  cat("- R²:", round(overall_stats$mean_rsq, 3), "±", round(overall_stats$sd_rsq, 3), "\n")
  cat("- RPD:", round(overall_stats$mean_rpd, 2), "±", round(overall_stats$sd_rpd, 2), "\n")
  cat("- RPIQ:", round(overall_stats$mean_rpiq, 2), "±", round(overall_stats$sd_rpiq, 2), "\n\n")
  
  cat("Model classification (RPD-primary):\n")
  cat("- Excellent (RPD > 3.0):", excellent_count, "(", excellent_pct, "%)\n")
  cat("- Good (RPD 2.0-3.0):", good_count, "(", good_pct, "%)\n")
  cat("- Fair (RPD 1.5-2.0):", fair_count, "(", fair_pct, "%)\n")
  cat("- Poor (RPD < 1.5):", poor_count, "(", poor_pct, "%)\n\n")
  
  cat("Supporting metric concordance:\n")
  cat("- Excellent with full support:", excellent_support_pct, "%\n")
  cat("- Good with full support:", good_support_pct, "%\n")
  cat("- Fair with full support:", fair_support_pct, "%\n\n")
  
  cat("Summary:\n")
  cat("- Quantitative capable:", quantitative_capable, "%\n")
  cat("- Total acceptable:", total_acceptable, "%\n")
  
  return(list(
    overall_stats = overall_stats,
    performance_summary = performance_summary,
    metrics = metrics,
    predictions = predictions,
    classification_method = "rpd_primary_fresh"
  ))
}


create_model_comparison_table <- function(model_comparison) {
  # Create model comparison table
  
  knitr::kable(
    model_comparison$comparison_table,
    caption = "Comparison of Full Spectrum vs Protein-Focused Models",
    row.names = FALSE
  )
}

create_calibration_plot <- function(results) {
  if ("model_n_comp_statistics" %in% names(results)) {
    results$model_n_comp_statistics |> 
      ggplot(aes(as.factor(ncomp), RMSE)) + 
      geom_line(aes(group = id), alpha = 0.03) + 
      theme_classic() + 
      xlab("Crude Protein Model Number of Components") + 
      ylab("Crude Protein Model Root Mean Squared Error")
  } else {
    ggplot() + 
      geom_text(aes(x = 0.5, y = 0.5, 
                    label = "model_n_comp_statistics missing\nNeed component progression in modeling"), 
                size = 4) +
      theme_void()
  }
}
create_performance_boxplot <- function(final_model_analysis) {
  # Create performance boxplot
  
  metrics_long <- melt(final_model_analysis$metrics, 
                       id.vars = "iteration",
                       measure.vars = c("rmse", "rsq", "rpd", "rpiq"),
                       variable.name = "metric", value.name = "value")
  
  ggplot(metrics_long, aes(x = metric, y = value)) +
    geom_boxplot(alpha = 0.7) +
    facet_wrap(~metric, scales = "free_y") +
    labs(
      title = "Model Performance Distribution",
      x = "Performance Metric",
      y = "Value"
    ) +
    theme_minimal()
}

create_algorithm_comparison_plot <- function(multi_algo_analysis) {
  # Create algorithm comparison plot
  
  plot_data <- multi_algo_analysis$raw_results
  
  ggplot(plot_data, aes(x = algorithm, y = rmse)) +
    geom_boxplot(alpha = 0.7) +
    labs(
      title = "Algorithm Performance Comparison",
      x = "Algorithm",
      y = "RMSE (g/kg)",
      fill = "Algorithm"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}

create_vip_plot <- function(spectral_analysis, threshold = 1.0, use_points = FALSE) {
  # Create VIP scores plot with option for points or lines
  
  vip_data <- spectral_analysis$vip_scores
  
  p <- ggplot(vip_data, aes(x = wavelength, y = vip_score))
  
  # Use points or lines based on parameter
  if (use_points) {
    p <- p + geom_point(color = "steelblue", alpha = 0.7, size = 1.5)
  } else {
    p <- p + geom_line(color = "steelblue", alpha = 0.7)
  }
  
  p <- p +
    geom_point(data = subset(vip_data, important), 
               aes(color = "VIP ≥ 1.0"), size = 2) +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red", alpha = 0.8) +
    scale_color_manual(values = c("VIP ≥ 1.0" = "red")) +
    labs(
      title = "Variable Importance in Projection (VIP) Scores",
      x = "Wavelength (nm)",
      y = "VIP Score",
      color = ""
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}
create_model_comparison_plot <- function(spectral_analysis, protein_focused_analysis) {
  # Create side-by-side model comparison plot with lines for full spectrum, points for protein-focused
  
  # Combine VIP data from both models
  full_vip <- spectral_analysis$vip_scores
  full_vip$model <- "Full Spectrum"
  
  protein_vip <- protein_focused_analysis$vip_scores
  protein_vip$model <- "Protein-Focused"
  
  combined_vip <- rbind(full_vip, protein_vip)
  
  ggplot(combined_vip, aes(x = wavelength, y = vip_score)) +
    # Lines for full spectrum
    geom_line(data = subset(combined_vip, model == "Full Spectrum"), 
              color = "#FF6B6B", alpha = 0.7) +
    # Points for protein-focused
    geom_point(data = subset(combined_vip, model == "Protein-Focused"), 
               color = "#4ECDC4", alpha = 0.8, size = 2) +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", alpha = 0.8) +
    facet_wrap(~model, scales = "free_x") +
    labs(
      title = "VIP Scores: Full Spectrum vs Protein-Focused Models",
      x = "Wavelength (nm)",
      y = "VIP Score"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
}
# =============================================================================
# FIXED PLOTTING FUNCTIONS FOR FIGURES 2 AND 3
# Add these functions to your R/core_functions.R file
# =============================================================================

# FIGURE 2: Performance Boxplot (Fixed version from backup)
create_performance_boxplot <- function(analysis) {
  analysis$metrics |>
    pivot_longer(
      cols = c(rmse, rsq, rpd, rpiq),
      names_to = "metric",
      values_to = "value"
    ) |>
    mutate(metric_label = toupper(metric)) |>
    ggplot(aes(x = metric, y = value)) +
    theme_classic() +
    geom_boxplot() +
    facet_wrap(
      ~ factor(metric_label, levels = c("RMSE", "RSQ", "RPD", "RPIQ")),
      scales = "free",
      nrow = 1,
      labeller = as_labeller(c(
        "RMSE" = "RMSE", "RSQ" = "R^2", "RPIQ" = "RPIQ", "RPD" = "RPD"
      ))
    ) +
    labs(
      title = "Model Performance Distribution",
      y = "Estimate"
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
}

# FIGURE 3: Validation Error Plot (Fixed version from backup)
create_validation_error_plot <- function(error_analysis) {
  
  # Check for required data
  if (!"raw_predictions" %in% names(error_analysis)) {
    cat("❌ raw_predictions not found in error_analysis\n")
    return(create_error_placeholder())
  }
  
  # Start with clean copy
  raw_data <- data.table::copy(error_analysis$raw_predictions)
  
  cat("DEBUG: Starting with", nrow(raw_data), "rows\n")
  cat("DEBUG: Columns:", paste(names(raw_data), collapse = ", "), "\n")
  
  # Ensure basic error calculation
  if (!"error_raw" %in% names(raw_data)) {
    raw_data[, error_raw := predicted - actual]
  }
  
  # STEP 1: Remove existing plot_order to avoid merge conflicts
  if ("plot_order" %in% names(raw_data)) {
    raw_data[, plot_order := NULL]
    cat("DEBUG: Removed existing plot_order column\n")
  }
  
  # STEP 2: Create sample-level summaries with fresh plot_order
  sample_summaries <- raw_data[, .(
    mean_actual = mean(actual, na.rm = TRUE),
    mean_error = mean(error_raw, na.rm = TRUE),
    n_obs = .N
  ), by = ith_in_data_set]
  
  # Add plot order to sample summaries
  sample_summaries <- sample_summaries[order(mean_actual)]
  sample_summaries[, plot_order := 1:.N]
  
  # Add tertiles to sample summaries
  cutpoints <- quantile(sample_summaries$mean_actual, probs = c(0, 1/3, 2/3, 1))
  sample_summaries[, cutpoints := cut(mean_actual, breaks = cutpoints, include.lowest = TRUE)]
  levels(sample_summaries$cutpoints) <- c("Lowest~Tertile", "Middle~Tertile", "Highest~Tertile")
  
  # Add systematic bias to sample summaries
  lm_mod <- lm(mean_error ~ mean_actual, data = sample_summaries)
  sample_summaries[, systematic_bias := predict(lm_mod)]
  
  cat("DEBUG: Sample summaries complete with", nrow(sample_summaries), "rows\n")
  
  # STEP 3: Clean merge - now no plot_order conflicts
  plot_data <- merge(
    raw_data, 
    sample_summaries[, .(ith_in_data_set, plot_order, cutpoints, systematic_bias)], 
    by = "ith_in_data_set",
    all.x = TRUE
  )
  
  cat("DEBUG: After merge, plot_data has", nrow(plot_data), "rows\n")
  cat("DEBUG: Plot data columns:", paste(names(plot_data), collapse = ", "), "\n")
  
  # Verify plot_order exists (should now work!)
  if (!"plot_order" %in% names(plot_data)) {
    cat("❌ plot_order STILL missing! Using fallback\n")
    return(create_error_placeholder())
  }
  
  cat("✅ plot_order successfully created! Range:", range(plot_data$plot_order, na.rm = TRUE), "\n")
  
  # STEP 4: Create the plot
  plot_data <- plot_data[order(plot_order)]
  
  p <- ggplot(plot_data) +
    # Background points (very light triangles)
    geom_point(aes(x = plot_order, y = error_raw), alpha = 0.05, shape = 2) +
    # Heavy dashed horizontal line at zero
    geom_hline(yintercept = 0, linewidth = 2, lty = 2) +
    # CROSSES for systematic bias
    geom_point(data = sample_summaries, 
               aes(x = plot_order, y = systematic_bias), 
               shape = 3, size = 1.2, color = "black") +
    # Faceting
    facet_wrap(~cutpoints, 
               labeller = label_parsed, 
               scales = "free_x") +
    # Y-axis label matching backup
    ylab("Crude Protein Predicted Percent Difference\nfrom Assayed Value") +
    theme_classic() +
    # Remove x-axis labeling
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  return(p)
}

# Helper function for placeholder plot
create_error_placeholder <- function() {
  placeholder_data <- data.frame(
    x = seq(200, 350, length.out = 50),
    y = rnorm(50, 0, 3) + 0.01 * (seq(200, 350, length.out = 50) - 275)
  )
  
  ggplot(placeholder_data, aes(x = x, y = y)) +
    geom_point(alpha = 0.6, color = "gray50") +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.3, color = "orange") +
    theme_classic() +
    labs(
      title = "Validation Error Analysis (Placeholder)",
      subtitle = "Run diagnostic to understand your error_analysis structure",
      x = "Actual Protein Concentration (g/kg)",
      y = "Prediction Error (g/kg)",
      caption = "Note: This is placeholder data - check error_analysis structure"
    ) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      plot.caption = element_text(face = "italic", color = "red")
    ) +
    annotate("text", x = 275, y = 0, 
             label = "Check error_analysis structure\nwith debug script", 
             hjust = 0.5, vjust = -1, size = 4, color = "red")
}

# =============================================================================
# FIGURE GENERATION FUNCTION
# =============================================================================

generate_all_figures_safely <- function() {
  # Safely generate all figures with error handling
  
  cat("=== GENERATING ALL FIGURES ===\n")
  
  figures <- list()
  
  # Figure 1: Model Calibration
  tryCatch({
    if (tar_exist_objects("final_model_results")) {
      final_model_results <- tar_read(final_model_results)
      figures$calibration <- create_calibration_plot(final_model_results)
      cat("✅ Calibration plot generated\n")
    } else {
      cat("❌ final_model_results not available\n")
    }
  }, error = function(e) {
    cat("❌ Error creating calibration plot:", e$message, "\n")
  })
  
  # Figure 2: Performance Distribution
  tryCatch({
    if (tar_exist_objects("final_model_analysis")) {
      final_model_analysis <- tar_read(final_model_analysis)
      figures$performance <- create_performance_boxplot(final_model_analysis)
      cat("✅ Performance boxplot generated\n")
    } else {
      cat("❌ final_model_analysis not available\n")
    }
  }, error = function(e) {
    cat("❌ Error creating performance plot:", e$message, "\n")
  })
  
  # Figure 3: Validation Errors
  tryCatch({
    if (tar_exist_objects("error_analysis")) {
      error_analysis <- tar_read(error_analysis)
      figures$validation <- create_validation_error_plot(error_analysis)
      cat("✅ Validation error plot generated\n")
    } else {
      cat("❌ error_analysis not available\n")
    }
  }, error = function(e) {
    cat("❌ Error creating validation plot:", e$message, "\n")
  })
  
  cat("=== FIGURE GENERATION COMPLETE ===\n")
  return(figures)
}

# =============================================================================
# UTILITY FUNCTIONS FOR DEBUGGING
# =============================================================================

check_figure_dependencies <- function() {
  # Check if all dependencies for figures are available
  
  cat("Checking figure dependencies...\n")
  
  required_targets <- c(
    "final_model_results",
    "final_model_analysis", 
    "error_analysis"
  )
  
  results <- list()
  for (target in required_targets) {
    exists <- tar_exist_objects(target)
    results[[target]] <- exists
    status <- if (exists) "✅" else "❌"
    cat(sprintf("  %s %s\n", status, target))
  }
  
  return(results)
}

run_protein_focused_validation <- function(hemp_data, best_method, n_iterations = 1000) {
  # Run protein-focused validation with VIP calculations
  
  cat("Starting protein-focused validation with", n_iterations, "iterations\n")
  
  # Select protein wavelengths  
  protein_selection <- select_protein_wavelengths(hemp_data)
  protein_data <- hemp_data[, c("crude_protein", protein_selection$column_names), with = FALSE]
  
  spectral_cols <- protein_selection$column_names
  spectra_matrix <- as.matrix(protein_data[, ..spectral_cols])
  y <- protein_data$crude_protein
  
  # Convert method to numeric ID
  method_id <- ifelse(is.character(best_method), 
                      c("Raw spectra"=1, "First derivative"=2, "Savitzky-Golay smoothing"=3, 
                        "Gap-segment derivative"=4, "Standard normal variate (SNV)"=5, 
                        "SNV + Savitzky-Golay"=6, "SNV + detrending"=7, 
                        "Multiplicative scatter correction"=8)[best_method], 
                      best_method)
  
  metrics <- data.table()
  vip_scores_all <- data.table()
  coefficients_all <- data.table()
  
  for (i in 1:n_iterations) {
    if (i %% 100 == 0) cat("Iteration", i, "\n")
    
    tryCatch({
      # SAME split as full spectrum analysis
      set.seed(i)
      train_idx <- sample(nrow(hemp_data), size = floor(0.75 * nrow(hemp_data)))
      test_idx <- setdiff(1:nrow(hemp_data), train_idx)
      
      x_train <- spectra_matrix[train_idx, ]
      y_train <- y[train_idx]
      x_test <- spectra_matrix[test_idx, ]
      y_test <- y[test_idx]
      
      # Apply preprocessing and fit model
      x_train_proc <- apply_preprocessing(x_train, method_id)
      x_test_proc <- apply_preprocessing(x_test, method_id)
      
      train_data <- data.frame(y = y_train, x_train_proc)
      model <- train(y ~ ., data = train_data, method = "pls", 
                     trControl = trainControl(method = "cv", number = 10), 
                     tuneLength = 15)
      
      final_preds <- predict(model, data.frame(x_test_proc))
      
      # Calculate metrics
      final_rmse <- sqrt(mean((y_test - final_preds)^2))
      final_rsq <- cor(y_test, final_preds)^2
      final_rpd <- sd(y_test) / final_rmse
      final_rpiq <- (quantile(y_test, 0.75) - quantile(y_test, 0.25)) / final_rmse
      
      metrics <- rbind(metrics, data.table(
        iteration = i, optimal_components = model$bestTune$ncomp,
        rmse = final_rmse, rsq = final_rsq, rpd = final_rpd, rpiq = final_rpiq
      ))
      
      # Calculate VIP scores
      pls_object <- model$finalModel
      optimal_ncomp <- model$bestTune$ncomp
      
      # Get PLS components
      W <- pls_object$loading.weights[, 1:optimal_ncomp, drop = FALSE]
      T <- pls_object$scores[, 1:optimal_ncomp, drop = FALSE]
      
      # Calculate VIP scores
      SS <- colSums(T^2)
      SS_total <- sum(SS)
      p <- nrow(W)
      vip_scores <- sqrt(p * rowSums((W^2) %*% diag(SS / SS_total)) / optimal_ncomp)
      
      # Get coefficients
      coefficients <- coef(pls_object, ncomp = optimal_ncomp)
      coefficients <- as.numeric(coefficients)
      
      # Handle length mismatch (common with derivatives)
      n_vars <- length(vip_scores)
      used_wavelengths <- protein_selection$wavelengths[1:n_vars]
      used_assignments <- protein_selection$assignments[1:n_vars]
      
      # Store VIP scores
      iteration_vip <- data.table(
        iteration = i,
        wavelength = used_wavelengths,
        assignment = used_assignments,
        vip_score = vip_scores,
        important = vip_scores > 1.0
      )
      vip_scores_all <- rbind(vip_scores_all, iteration_vip)
      
      # Store coefficients
      iteration_coeff <- data.table(
        iteration = i,
        wavelength = used_wavelengths,
        assignment = used_assignments,
        coefficient = coefficients[1:n_vars],
        abs_coefficient = abs(coefficients[1:n_vars])
      )
      coefficients_all <- rbind(coefficients_all, iteration_coeff)
      
    }, error = function(e) cat("Error in iteration", i, ":", e$message, "\n"))
  }
  
  return(list(
    metrics = metrics, 
    vip_scores = vip_scores_all,
    coefficients = coefficients_all,
    protein_selection = protein_selection, 
    method = best_method, 
    n_wavelengths = length(protein_selection$wavelengths)
  ))
}

analyze_protein_focused_results <- function(protein_results) {
  # Analyze protein-focused results with variable importance
  
  metrics <- protein_results$metrics
  vip_scores <- protein_results$vip_scores
  coefficients <- protein_results$coefficients
  
  # Overall statistics
  overall_stats <- metrics[, .(
    mean_rmse = mean(rmse, na.rm = TRUE), sd_rmse = sd(rmse, na.rm = TRUE),
    mean_rsq = mean(rsq, na.rm = TRUE), sd_rsq = sd(rsq, na.rm = TRUE),
    mean_rpd = mean(rpd, na.rm = TRUE), sd_rpd = sd(rpd, na.rm = TRUE),
    mean_rpiq = mean(rpiq, na.rm = TRUE), sd_rpiq = sd(rpiq, na.rm = TRUE),
    mean_components = mean(optimal_components, na.rm = TRUE), n_successful = .N
  )]
  
  # RPD-based classification
  excellent_count <- sum(metrics$rpd > 3.0, na.rm = TRUE)
  good_count <- sum(metrics$rpd >= 2.0 & metrics$rpd <= 3.0, na.rm = TRUE)
  fair_count <- sum(metrics$rpd >= 1.4 & metrics$rpd < 2.0, na.rm = TRUE)
  poor_count <- sum(metrics$rpd < 1.4, na.rm = TRUE)
  
  total_models <- nrow(metrics)
  performance_summary <- list(
    excellent_percent = round(100 * excellent_count / total_models, 1),
    good_percent = round(100 * good_count / total_models, 1),
    fair_percent = round(100 * fair_count / total_models, 1),
    poor_percent = round(100 * poor_count / total_models, 1),
    quantitative_capable = round(100 * (excellent_count + good_count) / total_models, 1),
    total_acceptable = round(100 * (excellent_count + good_count + fair_count) / total_models, 1),
    n_wavelengths = protein_results$n_wavelengths
  )
  
  # Variable importance analysis
  vip_summary <- vip_scores[, .(
    mean_vip = mean(vip_score, na.rm = TRUE),
    sd_vip = sd(vip_score, na.rm = TRUE),
    median_vip = median(vip_score, na.rm = TRUE),
    importance_frequency = mean(important, na.rm = TRUE) * 100,  # % of times VIP > 1.0
    n_iterations = .N
  ), by = .(wavelength, assignment)]
  
  # Sort by importance frequency
  vip_summary <- vip_summary[order(-importance_frequency)]
  
  # Coefficient analysis
  coeff_summary <- coefficients[, .(
    mean_coeff = mean(coefficient, na.rm = TRUE),
    sd_coeff = sd(coefficient, na.rm = TRUE),
    mean_abs_coeff = mean(abs_coefficient, na.rm = TRUE),
    sd_abs_coeff = sd(abs_coefficient, na.rm = TRUE)
  ), by = .(wavelength, assignment)]
  
  # Combine VIP and coefficient summaries
  variable_importance <- merge(vip_summary, coeff_summary, by = c("wavelength", "assignment"))
  variable_importance <- variable_importance[order(-importance_frequency)]
  
  cat("Protein-focused analysis complete:", total_models, "models\n")
  cat("Wavelengths used:", protein_results$n_wavelengths, "\n")
  cat("Quantitative capable:", performance_summary$quantitative_capable, "%\n")
  cat("Top important wavelengths (>80% frequency):", sum(variable_importance$importance_frequency > 80), "\n")
  
  return(list(
    overall_stats = overall_stats, 
    performance_summary = performance_summary, 
    metrics = metrics,
    variable_importance = variable_importance,
    vip_scores = vip_scores,
    coefficients = coefficients
  ))
}

create_protein_focused_table <- function(protein_analysis) {
  # Create summary table for protein-focused results
  
  stats <- protein_analysis$overall_stats
  summary <- protein_analysis$performance_summary
  
  table_data <- data.frame(
    "Metric" = c("RMSE (g/kg)", "R²", "RPD", "RPIQ", "Excellent (%)", "Good (%)", 
                 "Fair (%)", "Poor (%)", "Quantitative (%)", "Total Usable (%)"),
    "Value" = c(
      paste0(round(stats$mean_rmse, 2), " ± ", round(stats$sd_rmse, 2)),
      paste0(round(stats$mean_rsq, 3), " ± ", round(stats$sd_rsq, 3)),
      paste0(round(stats$mean_rpd, 2), " ± ", round(stats$sd_rpd, 2)),
      paste0(round(stats$mean_rpiq, 2), " ± ", round(stats$sd_rpiq, 2)),
      paste0(summary$excellent_percent, "%"),
      paste0(summary$good_percent, "%"),
      paste0(summary$fair_percent, "%"),
      paste0(summary$poor_percent, "%"),
      paste0(summary$quantitative_capable, "%"),
      paste0(summary$total_acceptable, "%")
    ),
    check.names = FALSE
  )
  
  return(knitr::kable(table_data, caption = "Protein-Focused Model Performance Summary", 
                      row.names = FALSE))
}

create_variable_importance_table <- function(protein_analysis, top_n = 10) {
  # Create table showing most important variables
  
  var_importance <- protein_analysis$variable_importance[1:top_n]
  
  table_data <- data.frame(
    "Rank" = 1:nrow(var_importance),
    "Wavelength (nm)" = var_importance$wavelength,
    "Assignment" = var_importance$assignment,
    "VIP Frequency (%)" = round(var_importance$importance_frequency, 1),
    "Mean VIP" = round(var_importance$mean_vip, 2),
    "Mean |Coefficient|" = round(var_importance$mean_abs_coeff, 4),
    check.names = FALSE
  )
  
  return(knitr::kable(table_data, 
                      caption = paste0("Top ", top_n, " Most Important Protein Wavelengths"), 
                      row.names = FALSE))
}
# Wrapper function that combines everything
run_complete_protein_focused_analysis <- function(hemp_data, best_method, n_iterations = 1000) {
  # Complete protein-focused analysis workflow
  
  cat("=== COMPLETE PROTEIN-FOCUSED ANALYSIS ===\n")
  
  # Run validation
  validation_results <- run_protein_focused_validation(hemp_data, best_method, n_iterations)
  
  # Analyze results
  analysis_results <- analyze_protein_focused_results(validation_results)
  
  # Create summary table
  summary_table <- create_protein_focused_summary_table(analysis_results)
  
  return(list(
    validation_results = validation_results,
    analysis_results = analysis_results,
    summary_table = summary_table
  ))
}

create_combined_vip_distribution_plot <- function(spectral_analysis, protein_analysis, alpha_level = 0.1) {
  # Plot VIP scores from both full spectrum and protein-focused analyses
  
  # Prepare full spectrum data
  full_vip <- copy(spectral_analysis$vip_scores)
  full_vip$model <- "Full Spectrum"
  
  # Prepare protein-focused data
  protein_vip <- copy(protein_analysis$vip_scores)
  protein_vip$model <- "Protein-Focused"
  
  # Select only common columns (no assignment column)
  common_cols <- c("iteration", "wavelength", "vip_score", "important", "model")
  
  full_vip_clean <- full_vip[, ..common_cols]
  protein_vip_clean <- protein_vip[, ..common_cols]
  
  # Combine data
  combined_vip <- rbind(full_vip_clean, protein_vip_clean)
  
  # Create the plot
  p <- ggplot(combined_vip, aes(x = wavelength, y = vip_score)) +
    geom_point(alpha = alpha_level, size = 0.6, color = "black") +
    geom_hline(yintercept = 1.0,  linewidth = 2, lty = 2) +
    facet_wrap(~model, scales = "free_x") +
    labs(
      title = "VIP Score Distribution: Full Spectrum vs Protein-Focused Models",
      x = "Wavelength (nm)",
      y = "VIP Score"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      strip.text = element_text(size = 12, face = "bold")
    )
  
  return(p)
}

library(data.table)
library(knitr)

# Function to categorize model performance based on RPD values only
categorize_performance <- function(rpd) {
  # Based on the paper's established RPD thresholds
  ifelse(rpd > 3.0, "Excellent",
         ifelse(rpd >= 2.5, "Good", 
                ifelse(rpd >= 2.0, "Approximate Quantitative",
                       ifelse(rpd >= 1.5, "Qualitative Screening", "Poor"))))
}

# Process the analysis results
process_analysis_results <- function(protein_analysis, spectral_analysis) {
  
  # Convert to data.table for easier processing
  protein_metrics <- as.data.table(protein_analysis$metrics)
  spectral_metrics <- as.data.table(spectral_analysis$metrics)
  
  # Add performance categories based on RPD only
  protein_metrics[, performance_category := categorize_performance(rpd)]
  spectral_metrics[, performance_category := categorize_performance(rpd)]
  
  # Calculate summary statistics for each analysis
  protein_summary <- protein_metrics[, .(
    Analysis = "Protein-Focused",
    Mean_RMSE = round(mean(rmse), 2),
    SD_RMSE = round(sd(rmse), 2),
    Mean_R2 = round(mean(rsq), 3),
    SD_R2 = round(sd(rsq), 3),
    Mean_RPD = round(mean(rpd), 2),
    SD_RPD = round(sd(rpd), 2),
    Mean_RPIQ = round(mean(rpiq), 2),
    SD_RPIQ = round(sd(rpiq), 2),
    Mean_Components = round(mean(optimal_components), 1),
    SD_Components = round(sd(optimal_components), 1)
  )]
  
  spectral_summary <- spectral_metrics[, .(
    Analysis = "Full Spectrum",
    Mean_RMSE = round(mean(rmse), 2),
    SD_RMSE = round(sd(rmse), 2),
    Mean_R2 = round(mean(rsq), 3),
    SD_R2 = round(sd(rsq), 3),
    Mean_RPD = round(mean(rpd), 2),
    SD_RPD = round(sd(rpd), 2),
    Mean_RPIQ = round(mean(rpiq), 2),
    SD_RPIQ = round(sd(rpiq), 2),
    Mean_Components = round(mean(optimal_components), 1),
    SD_Components = round(sd(optimal_components), 1)
  )]
  
  # Combine summaries
  combined_summary <- rbind(protein_summary, spectral_summary)
  
  # Create performance distribution table
  protein_performance <- protein_metrics[, .N, by = performance_category][
    , .(Category = performance_category, 
        Protein_Count = N, 
        Protein_Percent = round(N/1000 * 100, 1))]
  
  spectral_performance <- spectral_metrics[, .N, by = performance_category][
    , .(Category = performance_category, 
        Spectral_Count = N, 
        Spectral_Percent = round(N/1000 * 100, 1))]
  
  # Merge performance distributions
  performance_dist <- merge(protein_performance, spectral_performance, 
                            by = "Category", all = TRUE)
  performance_dist[is.na(performance_dist)] <- 0
  
  # Order by performance level
  performance_order <- c("Excellent", "Good", "Approximate Quantitative", 
                         "Qualitative Screening", "Poor")
  performance_dist <- performance_dist[match(performance_order, performance_dist$Category)]
  performance_dist <- performance_dist[!is.na(performance_dist$Category)]
  
  return(list(
    summary_table = combined_summary,
    performance_distribution = performance_dist,
    protein_metrics = protein_metrics,
    spectral_metrics = spectral_metrics
  ))
}

# Create formatted tables
create_summary_table <- function(results) {
  kable(results$summary_table, 
        caption = "Model Performance Summary Statistics (n=1000 iterations)",
        col.names = c("Analysis Type", "RMSE (g/kg)", "SD", "R²", "SD", 
                      "RPD", "SD", "RPIQ", "SD", "Components", "SD"),
        align = c("l", rep("c", 10)))
}

create_performance_distribution_table <- function(results) {
  # Add RPD ranges for clarity - based on established thresholds
  rpd_ranges <- c(
    "RPD > 3.0",
    "RPD 2.5-3.0", 
    "RPD 2.0-2.5",
    "RPD 1.5-2.0",
    "RPD < 1.5"
  )
  
  table_data <- results$performance_distribution
  table_data$RPD_Range <- rpd_ranges[1:nrow(table_data)]
  
  # Select only the 4 desired columns
  table_data <- table_data[, c("Category", "RPD_Range", "Protein_Percent", "Spectral_Percent")]
  
  kable(table_data,
        caption = "Distribution of Model Performance Categories (Based on RPD)",
        col.names = c("Performance Category", "RPD Criteria", 
                      "Protein-Focused (%)", "Full Spectrum (%)"),
        align = c("l", "l", "c", "c"))
}
