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
                      "2" = prospectr::savitzkyGolay(spectra_matrix, m = 1, p = 3, w = 5),  # First derivative
                      "3" = prospectr::savitzkyGolay(spectra_matrix, m = 0, p = 3, w = 5),  # SG smoothing
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

run_final_modeling <- function(hemp_data, best_method, n_iterations = 1000) {
  # Run final modeling with best method
  
  cat("Running final modeling with", best_method, "for", n_iterations, "iterations\n")
  
  # Extract data
  spectral_cols <- grep("^x[0-9]+", names(hemp_data), value = TRUE)
  spectra_matrix <- as.matrix(hemp_data[, ..spectral_cols])
  protein_data <- hemp_data$crude_protein
  
  # Storage for results
  metrics <- data.table()
  predictions <- data.table()
  
  set.seed(123)
  
  for(i in 1:n_iterations) {
    if (i %% 100 == 0) cat("Iteration", i, "\n")
    
    tryCatch({
      # Train/test split
      train_idx <- sample(nrow(hemp_data), size = floor(0.75 * nrow(hemp_data)))
      test_idx <- setdiff(1:nrow(hemp_data), train_idx)
      
      x_train <- spectra_matrix[train_idx, ]
      y_train <- protein_data[train_idx]
      x_test <- spectra_matrix[test_idx, ]
      y_test <- protein_data[test_idx]
      
      # Apply preprocessing
      x_train_proc <- apply_preprocessing(x_train, best_method)
      x_test_proc <- apply_preprocessing(x_test, best_method)
      
      # Fit model
      results <- fit_pls_model(x_train_proc, y_train, x_test_proc, y_test)
      
      # Store metrics
      iteration_metrics <- data.table(
        iteration = i,
        n_components = results$optimal_components,
        rmse = results$rmse,
        rsq = results$rsq,
        rpd = results$rpd,
        rpiq = results$rpiq
      )
      metrics <- rbind(metrics, iteration_metrics)
      
      # Store predictions
      iteration_predictions <- data.table(
        iteration = i,
        sample_id = test_idx,
        actual = results$observed,
        predicted = results$predictions
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
    method = best_method
  ))
}

analyze_final_model <- function(final_model_results) {
  # Analyze final model performance
  
  metrics <- final_model_results$metrics
  predictions <- final_model_results$predictions
  
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
  
  # Performance categories
  good_models <- metrics[rpd >= 2.0 & rsq >= 0.7, .N]
  excellent_models <- metrics[rpd >= 2.5 & rsq >= 0.8, .N]
  total_models <- nrow(metrics)
  
  performance_summary <- list(
    total_models = total_models,
    good_models = good_models,
    excellent_models = excellent_models,
    good_percent = round(100 * good_models / total_models, 1),
    excellent_percent = round(100 * excellent_models / total_models, 1)
  )
  
  cat("Final model analysis:\n")
  cat("- Total successful models:", total_models, "\n")
  cat("- Good models (RPD≥2.0, R²≥0.7):", good_models, "(", performance_summary$good_percent, "%)\n")
  cat("- Excellent models (RPD≥2.5, R²≥0.8):", excellent_models, "(", performance_summary$excellent_percent, "%)\n")
  
  return(list(
    overall_stats = overall_stats,
    performance_summary = performance_summary,
    metrics = metrics,
    predictions = predictions
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
  # Analyze multi-algorithm comparison results
  
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
  
  # Model quality assessment
  model_quality <- multi_algorithm_results[, .(
    excellent_models = sum(rpd >= 2.5 & rsq >= 0.8, na.rm = TRUE),
    good_models = sum(rpd >= 2.0 & rsq >= 0.7, na.rm = TRUE),
    quantitative_capable = round(100 * sum(rpd >= 2.0, na.rm = TRUE) / .N, 1),
    total_acceptable = round(100 * sum(rpd >= 1.5, na.rm = TRUE) / .N, 1)
  ), by = algorithm]
  
  cat("Multi-algorithm analysis complete for", length(unique(summary_stats$algorithm)), "algorithms\n")
  
  return(list(
    summary_stats = summary_stats,
    model_quality = model_quality,
    raw_results = multi_algorithm_results
  ))
}

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
    list(start = 1180, end = 1230, name = "C-H stretch 2nd overtone (amino acids)"),
    list(start = 1480, end = 1530, name = "N-H stretch 1st overtone (peptide bonds)"),
    list(start = 1660, end = 1700, name = "C-H stretch 1st overtone (aliphatic amino acids)"),
    list(start = 2040, end = 2070, name = "N-H + C-N combination (protein backbone)"),
    list(start = 2270, end = 2310, name = "C-H + C-H combination (amino acid structure)")
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

run_spectral_analysis <- function(data, best_method) {
  # Run full spectrum spectral analysis
  
  cat("Running full spectrum spectral analysis...\n")
  
  spectral_fit <- fit_pls_for_spectral_analysis(data, best_method)
  coeff_data <- extract_pls_coefficients(spectral_fit)
  vip_data <- calculate_vip_scores(spectral_fit)
  
  n_important_vip <- sum(vip_data$vip_score > 1.0)
  
  cat("Full spectrum analysis complete:\n")
  cat("- Important wavelengths (VIP > 1.0):", n_important_vip, "\n")
  cat("- Model components:", spectral_fit$model$bestTune$ncomp, "\n")
  
  list(
    spectral_fit = spectral_fit,
    coefficients = coeff_data,
    vip_scores = vip_data,
    analysis_stats = list(
      n_important_vip = n_important_vip,
      optimal_components = spectral_fit$model$bestTune$ncomp
    )
  )
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
  # Create algorithm comparison table
  
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
      Quantitative_formatted = paste0(round(quantitative_capable, 1), "%"),
      Total_formatted = paste0(round(total_acceptable, 1), "%")
    )
  
  final_table <- data.frame(
    "Algorithm" = formatted_table$Algorithm,
    "RMSE (±SD)" = formatted_table$RMSE_formatted,
    "R² (±SD)" = formatted_table$R2_formatted,
    "RPD (±SD)" = formatted_table$RPD_formatted,
    "RPIQ (±SD)" = formatted_table$RPIQ_formatted,
    "Quantitative Capable (%)" = formatted_table$Quantitative_formatted,
    "Total Acceptable (%)" = formatted_table$Total_formatted,
    check.names = FALSE
  )
  
  knitr::kable(
    final_table,
    caption = "Performance comparison of machine learning algorithms for hemp grain protein prediction",
    row.names = FALSE,
    align = c("l", rep("c", 6))
  )
}

create_model_comparison_table <- function(model_comparison) {
  # Create model comparison table
  
  knitr::kable(
    model_comparison$comparison_table,
    caption = "Comparison of Full Spectrum vs Protein-Focused Models",
    row.names = FALSE
  )
}

create_calibration_plot <- function(final_model_analysis) {
  # Create calibration plot
  
  # Sample predictions for plotting (to avoid overplotting)
  plot_data <- final_model_analysis$predictions[iteration <= 10]
  
  ggplot(plot_data, aes(x = actual, y = predicted)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    labs(
      title = "Model Calibration: Predicted vs Actual Protein Content",
      x = "Laboratory Measured Protein (g/kg)",
      y = "NIR Predicted Protein (g/kg)"
    ) +
    theme_minimal() +
    coord_equal()
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
  
  ggplot(plot_data, aes(x = algorithm, y = rmse, fill = algorithm)) +
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

create_vip_plot <- function(spectral_analysis, threshold = 1.0) {
  # Create VIP scores plot
  
  vip_data <- spectral_analysis$vip_scores
  
  ggplot(vip_data, aes(x = wavelength, y = vip_score)) +
    geom_line(color = "steelblue", alpha = 0.7) +
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
}

create_model_comparison_plot <- function(spectral_analysis, protein_focused_analysis) {
  # Create side-by-side model comparison plot
  
  # Combine VIP data from both models
  full_vip <- spectral_analysis$vip_scores
  full_vip$model <- "Full Spectrum"
  
  protein_vip <- protein_focused_analysis$vip_scores
  protein_vip$model <- "Protein-Focused"
  
  combined_vip <- rbind(full_vip, protein_vip)
  
  ggplot(combined_vip, aes(x = wavelength, y = vip_score, color = model)) +
    geom_line(alpha = 0.7) +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", alpha = 0.8) +
    facet_wrap(~model, scales = "free_x") +
    labs(
      title = "VIP Scores: Full Spectrum vs Protein-Focused Models",
      x = "Wavelength (nm)",
      y = "VIP Score",
      color = "Model"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
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
  if ("raw_predictions" %in% names(error_analysis)) {
    
    raw_data <- data.table::copy(error_analysis$raw_predictions)
    
    # Check if error_raw already exists, if not create it
    if (!"error_raw" %in% names(raw_data)) {
      raw_data[, error_raw := predicted - actual]
    }
    
    # Create sample-level data first, then assign plot_order
    # Aggregate to one row per sample
    sample_data <- raw_data[, .(
      mean_actual = mean(actual, na.rm = TRUE),
      mean_error = mean(error_raw, na.rm = TRUE)
    ), by = ith_in_data_set]
    
    # Create plot order based on sample ranking
    sample_data[, plot_order := rank(mean_actual)]
    
    # Create tertiles for faceting
    cutpoints <- quantile(sample_data$mean_actual, probs = c(0, 1/3, 2/3, 1))
    sample_data[, cutpoints := cut(mean_actual, breaks = cutpoints, include.lowest = TRUE)]
    levels(sample_data$cutpoints) <- c("Lowest~Tertile", "Middle~Tertile", "Highest~Tertile")
    
    # Add plot_order back to raw_data for the scatter points
    raw_data <- merge(raw_data, sample_data[, .(ith_in_data_set, plot_order, cutpoints)], 
                      by = "ith_in_data_set")
    
    # Fit linear model for systematic bias using sample-level data
    lm_mod <- lm(mean_error ~ mean_actual, data = sample_data)
    sample_data[, systematic_bias := predict(lm_mod, newdata = sample_data)]
    
    # Merge systematic bias back using sample data
    plot_data <- merge(raw_data, sample_data[, .(ith_in_data_set, systematic_bias)], 
                       by = "ith_in_data_set")
    
    # Create the plot
    p <- ggplot(plot_data, aes(x = plot_order, y = error_raw)) +
      geom_point(alpha = 0.3, size = 0.8, color = "gray60") +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.7) +
      geom_smooth(aes(y = systematic_bias), method = "lm", se = TRUE, 
                  color = "orange", linewidth = 1.2, alpha = 0.8) +
      facet_grid(
        . ~ cutpoints,
        scales = "free_x",
        space = "free_x",
        labeller = label_parsed
      ) +
      theme_classic() +
      labs(
        title = "Testing Set Prediction Errors by Sample",
        subtitle = "Samples ranked from lowest to highest actual CP concentration",
        x = "Sample Rank Order",
        y = "Prediction Error (g/kg)",
        caption = "Orange line shows systematic bias trend"
      ) +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        plot.caption = element_text(face = "italic", color = "gray50"),
        strip.text = element_text(size = 10)
      )
    
    return(p)
    
  } else {
    # Fallback to placeholder if data structure doesn't match
    return(create_error_placeholder())
  }
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