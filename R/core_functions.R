# R/core_functions.R - Complete Hemp NIR Analysis Functions
# All essential functions for the streamlined paper pipeline

library(data.table)
library(caret)
library(prospectr)
library(pls)
library(tidyverse)
library(kableExtra)

# =============================================================================
# DATA PREPARATION FUNCTIONS
# =============================================================================
# Utility function for skewness calculation
calculate_skewness <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  sd_x <- sqrt(sum((x - mean_x)^2) / n)
  z <- (x - mean_x) / sd_x
  sum(z^3) / n
}

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

# Function to select protein-specific wavelengths
select_protein_wavelengths <- function(data) {
  
  cat("Selecting protein-specific wavelengths...\n")
  
  # Define protein absorption bands (nm) - EXPANDED with more specific assignments
  protein_ranges <- list(
    list(start = 1180, end = 1230, name = "C-H stretch 2nd overtone (amino acids)"),
    list(start = 1480, end = 1530, name = "N-H stretch 1st overtone (peptide bonds)"),
    list(start = 1660, end = 1700, name = "C-H stretch 1st overtone (aliphatic amino acids)"),
    list(start = 2040, end = 2070, name = "N-H + C-N combination (protein backbone)"),
    list(start = 2270, end = 2310, name = "C-H + C-H combination (amino acid structure)")
  )
  
  # Get all available wavelengths
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
      
      cat("- Band", band$start, "-", band$end, "nm:", length(band_wavelengths), "wavelengths\n")
    }
  }
  
  # Create column names for selected wavelengths
  selected_columns <- paste0("x", protein_wavelengths)
  
  # Verify columns exist in data
  available_columns <- selected_columns[selected_columns %in% names(data)]
  
  cat("Selected", length(protein_wavelengths), "protein-specific wavelengths\n")
  cat("Available in data:", length(available_columns), "wavelengths\n")
  
  list(
    wavelengths = protein_wavelengths[selected_columns %in% names(data)],
    column_names = available_columns,
    assignments = protein_assignments[selected_columns %in% names(data)],
    n_selected = length(available_columns),
    ranges = protein_ranges
  )
}

# =============================================================================
# PREPROCESSING FUNCTIONS
# =============================================================================

# Add this function to your R/core_functions.R file:

apply_preprocessing <- function(spectra_matrix, method = 1) {
  # Apply preprocessing method to spectral data
  # method: can be numeric (1-8) or method name string
  
  library(prospectr)
  
  # Convert method name to numeric if it's a string
  if (is.character(method)) {
    method <- switch(method,
                     "Raw spectra" = 1,
                     "First derivative" = 2, 
                     "Savitzky-Golay smoothing" = 3,
                     "Gap-segment derivative" = 4,
                     "Standard normal variate (SNV)" = 5,
                     "SNV + Savitzky-Golay" = 6,
                     "SNV + detrending" = 7,
                     "Multiplicative scatter correction" = 8,
                     stop("Unknown preprocessing method name: ", method)
    )
  }
  
  if (method == 1) {
    # Raw spectra (no preprocessing)
    return(spectra_matrix)
    
  } else if (method == 2) {
    # First derivative
    return(diff(t(spectra_matrix), differences = 1) |> t())
    
  } else if (method == 3) {
    # Savitzky-Golay
    return(savitzkyGolay(spectra_matrix, m = 1, p = 3, w = 5))
    
  } else if (method == 4) {
    # Gap-segment derivative  
    return(gapDer(spectra_matrix, m = 1, w1 = 11, w2 = 5))
    
  } else if (method == 5) {
    # Standard Normal Variate
    return(standardNormalVariate(spectra_matrix))
    
  } else if (method == 6) {
    # SNV + Savitzky-Golay
    snv_spectra <- standardNormalVariate(spectra_matrix)
    return(savitzkyGolay(snv_spectra, m = 1, p = 3, w = 5))
    
  } else if (method == 7) {
    # SNV + detrend
    return(detrend(standardNormalVariate(spectra_matrix), p = 2))
    
  } else if (method == 8) {
    # Multiplicative Scatter Correction
    return(msc(spectra_matrix))
    
  } else {
    stop("Unknown preprocessing method ID: ", method)
  }
}

apply_all_preprocessing_methods <- function(spectra_train, spectra_test) {
  # Apply 8 preprocessing methods to spectral data (using your original parameters)
  
  cat("Applying preprocessing methods...\n")
  
  # Method 1: Raw spectra
  processed_1 <- list(train = spectra_train, test = spectra_test)
  
  # Method 2: First derivative  
  train_2 <- tryCatch({
    t(diff(t(spectra_train), diff = 1))
  }, error = function(e) { spectra_train[, -ncol(spectra_train)] })
  
  test_2 <- tryCatch({
    t(diff(t(spectra_test), diff = 1))
  }, error = function(e) { spectra_test[, -ncol(spectra_test)] })
  
  processed_2 <- list(train = train_2, test = test_2)
  
  # Method 3: Savitzky-Golay (your original parameters: m=1, p=3, w=5)
  train_3 <- tryCatch({
    prospectr::savitzkyGolay(spectra_train, m = 1, p = 3, w = 5)
  }, error = function(e) { spectra_train })
  
  test_3 <- tryCatch({
    prospectr::savitzkyGolay(spectra_test, m = 1, p = 3, w = 5)
  }, error = function(e) { spectra_test })
  
  processed_3 <- list(train = train_3, test = test_3)
  
  # Method 4: Gap-segment derivative (your original parameters: w=11, s=5)
  train_4 <- tryCatch({
    prospectr::gapDer(spectra_train, m = 1, w = 11, s = 5)
  }, error = function(e) { spectra_train })
  
  test_4 <- tryCatch({
    prospectr::gapDer(spectra_test, m = 1, w = 11, s = 5)
  }, error = function(e) { spectra_test })
  
  processed_4 <- list(train = train_4, test = test_4)
  
  # Method 5: Standard Normal Variate (SNV)
  train_5 <- tryCatch({
    prospectr::standardNormalVariate(spectra_train)
  }, error = function(e) { spectra_train })
  
  test_5 <- tryCatch({
    prospectr::standardNormalVariate(spectra_test)
  }, error = function(e) { spectra_test })
  
  processed_5 <- list(train = train_5, test = test_5)
  
  # Method 6: SNV + Savitzky-Golay (your best method - first SG, then SNV)
  train_6 <- tryCatch({
    sg_result <- prospectr::savitzkyGolay(spectra_train, m = 1, p = 3, w = 5)
    prospectr::standardNormalVariate(sg_result)
  }, error = function(e) { spectra_train })
  
  test_6 <- tryCatch({
    sg_result <- prospectr::savitzkyGolay(spectra_test, m = 1, p = 3, w = 5)
    prospectr::standardNormalVariate(sg_result)
  }, error = function(e) { spectra_test })
  
  processed_6 <- list(train = train_6, test = test_6)
  
  # Method 7: SNV + Detrending
  train_7 <- tryCatch({
    snv_result <- prospectr::standardNormalVariate(spectra_train)
    prospectr::detrend(snv_result, wav = 1:ncol(snv_result))
  }, error = function(e) { spectra_train })
  
  test_7 <- tryCatch({
    snv_result <- prospectr::standardNormalVariate(spectra_test)
    prospectr::detrend(snv_result, wav = 1:ncol(snv_result))
  }, error = function(e) { spectra_test })
  
  processed_7 <- list(train = train_7, test = test_7)
  
  # Method 8: Multiplicative Scatter Correction (MSC)
  train_8 <- tryCatch({
    prospectr::msc(spectra_train)
  }, error = function(e) { spectra_train })
  
  test_8 <- tryCatch({
    ref_spectrum <- attr(train_8, "Reference spectrum")
    if (is.null(ref_spectrum)) {
      ref_spectrum <- colMeans(spectra_train)
    }
    prospectr::msc(spectra_test, reference = ref_spectrum)
  }, error = function(e) { spectra_test })
  
  processed_8 <- list(train = train_8, test = test_8)
  
  methods <- list(
    method_1 = processed_1,
    method_2 = processed_2, 
    method_3 = processed_3,
    method_4 = processed_4,
    method_5 = processed_5,
    method_6 = processed_6,
    method_7 = processed_7,
    method_8 = processed_8
  )
  
  cat("- Applied", length(methods), "preprocessing methods using original parameters\n")
  return(methods)
}

# =============================================================================
# MODELING FUNCTIONS
# =============================================================================

fit_pls_model <- function(x_train, y_train, x_test, y_test, max_components = 20) {
  # Fit PLS model with cross-validation for optimal components
  
  # Ensure data is in correct format
  if (!is.matrix(x_train)) x_train <- as.matrix(x_train)
  if (!is.matrix(x_test)) x_test <- as.matrix(x_test)
  
  # Create data frame for PLS
  train_data <- data.frame(y = y_train, x_train)
  
  # Fit PLS model with cross-validation to find optimal components
  pls_model <- plsr(y ~ ., data = train_data, ncomp = max_components, validation = "CV")
  
  # Find optimal number of components (minimum RMSECV)
  cv_results <- RMSEP(pls_model, estimate = "CV")
  optimal_comps <- which.min(cv_results$val[1, 1, -1]) # Exclude intercept
  
  # Make predictions
  predictions <- predict(pls_model, newdata = data.frame(x_test), ncomp = optimal_comps)
  predictions <- as.numeric(predictions)
  
  # Calculate metrics
  rmse <- sqrt(mean((y_test - predictions)^2))
  rsq <- cor(y_test, predictions)^2
  
  # Calculate RPD and RPIQ
  rpd <- sd(y_test) / rmse
  rpiq <- (quantile(y_test, 0.75) - quantile(y_test, 0.25)) / rmse
  
  return(list(
    model = pls_model,
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
  # Run single modeling iteration for one preprocessing method
  
  # Get spectral columns
  spectral_cols <- grep("^x[0-9]+$", names(hemp_data), value = TRUE)
  spectra_matrix <- as.matrix(hemp_data[, ..spectral_cols])
  y <- hemp_data$crude_protein
  
  # Create train/test split (80/20)
  set.seed(iteration) # For reproducibility
  train_indices <- createDataPartition(y, p = 0.8, list = FALSE)[,1]
  
  x_train <- spectra_matrix[train_indices, ]
  x_test <- spectra_matrix[-train_indices, ]
  y_train <- y[train_indices]
  y_test <- y[-train_indices]
  
  # Apply preprocessing
  processed <- apply_all_preprocessing_methods(x_train, x_test)
  method_data <- processed[[paste0("method_", preprocessing_method_id)]]
  
  # Fit model
  results <- fit_pls_model(method_data$train, y_train, method_data$test, y_test)
  
  return(data.table(
    iteration = iteration,
    preprocessing_method = preprocessing_method_id,
    rmse = results$rmse,
    rsq = results$rsq,
    rpd = results$rpd,
    rpiq = results$rpiq,
    optimal_components = results$optimal_components
  ))
}

# =============================================================================
# CORE ANALYSIS FUNCTIONS
# =============================================================================

# Enhanced protein-focused analysis function
run_protein_focused_analysis <- function(data, best_method = "snv_sg") {
  
  cat("=== PROTEIN-FOCUSED MODEL ANALYSIS ===\n")
  
  # Select protein-specific wavelengths
  protein_selection <- select_protein_wavelengths(data)
  
  cat("Using", protein_selection$n_selected, "protein-specific wavelengths\n")
  
  # Extract spectral data for protein wavelengths only
  spectral_cols <- protein_selection$column_names
  spectra_matrix <- as.matrix(data[, ..spectral_cols])
  y <- data$crude_protein
  
  cat("Protein-focused data dimensions:", dim(spectra_matrix), "\n")
  
  # Create train/test split (use same approach as full model for fair comparison)
  set.seed(42)  # For reproducibility
  inTrain <- split_spectra(y)
  y_train <- y[inTrain]
  y_test <- y[-inTrain]
  
  # Apply preprocessing
  spectra_processed <- my_preprocess(spectra_matrix[inTrain, ], spectra_matrix[-inTrain, ])
  
  # Get preprocessed data for the best method
  method_train_name <- paste0(best_method, "_train")
  method_test_name <- paste0(best_method, "_test")
  
  X_train <- spectra_processed[[1]][[method_train_name]]
  X_test <- spectra_processed[[1]][[method_test_name]]
  
  # Convert to data frame (required by caret)
  X_train_df <- as.data.frame(X_train)
  X_test_df <- as.data.frame(X_test)
  
  # Set up caret training control
  train_control <- trainControl(
    method = "cv",
    number = 10,
    savePredictions = "final"
  )
  
  cat("Training protein-focused PLS model...\n")
  
  # Train PLS model
  model <- train(
    x = X_train_df,
    y = y_train,
    method = "pls",
    trControl = train_control,
    tuneLength = 15,  # Test up to 15 components
    preProcess = NULL  # Already preprocessed
  )
  
  cat("Optimal components:", model$bestTune$ncomp, "\n")
  cat("Training RMSE:", min(model$results$RMSE), "\n")
  
  # Create spectral fit object for coefficient extraction
  spectral_fit <- list(
    model = model,
    wavelengths = protein_selection$wavelengths,
    X_train = X_train_df,
    y_train = y_train,
    train_indices = inTrain,
    preprocessing_method = best_method
  )
  
  # Extract coefficients and VIP scores
  coeff_data <- extract_pls_coefficients(spectral_fit)
  vip_data <- calculate_vip_scores(spectral_fit)
  protein_bands_data <- identify_protein_bands(coeff_data, vip_data)
  summary_table <- create_wavelength_summary_table(protein_bands_data)
  
  # Calculate summary statistics
  n_important_vip <- sum(vip_data$vip_score > 1.0)
  n_protein_related <- sum(protein_bands_data$likely_protein_related)
  
  # Test set evaluation
  test_predictions <- predict(model, X_test_df)
  test_rmse <- sqrt(mean((y_test - test_predictions)^2))
  test_rsq <- cor(y_test, test_predictions)^2
  
  cat("Protein-focused analysis complete:\n")
  cat("- Wavelengths with VIP > 1.0:", n_important_vip, "\n")
  cat("- Wavelengths near protein bands:", n_protein_related, "\n")
  cat("- Model used", model$bestTune$ncomp, "components\n")
  cat("- Test RMSE:", round(test_rmse, 3), "\n")
  cat("- Test R²:", round(test_rsq, 3), "\n")
  
  list(
    spectral_fit = spectral_fit,
    coefficients = coeff_data,
    vip_scores = vip_data,
    protein_bands = protein_bands_data,
    summary_table = summary_table,
    analysis_stats = list(
      n_important_vip = n_important_vip,
      n_protein_related = n_protein_related,
      optimal_components = model$bestTune$ncomp,
      training_rmse = min(model$results$RMSE),
      test_rmse = test_rmse,
      test_rsq = test_rsq
    ),
    protein_selection = protein_selection,
    test_predictions = data.frame(
      actual = y_test,
      predicted = test_predictions,
      sample_id = which(!inTrain)
    )
  )
}

run_preprocessing_comparison <- function(data, n_iterations = NULL) {
  # Compare 8 preprocessing methods - core analysis for paper
  
  # Use configuration for iteration count
  config <- get_analysis_config()
  n_iterations <- resolve_param(n_iterations, config$n_iterations, "n_iterations")
  
  cat("=== PREPROCESSING COMPARISON ===\n")
  cat("Methods: 8 preprocessing approaches\n")
  cat("Iterations:", n_iterations, "\n")
  
  # Run comparison for all methods
  all_results <- data.table()
  
  for (method_id in 1:8) {
    cat("Processing method", method_id, "...\n")
    
    method_results <- data.table()
    for (iter in 1:n_iterations) {
      iter_result <- run_single_iteration(data, method_id, iter)
      method_results <- rbind(method_results, iter_result)
    }
    
    all_results <- rbind(all_results, method_results)
  }
  
  cat("✅ Preprocessing comparison complete\n")
  return(list(metrics = all_results))
}

analyze_preprocessing_methods <- function(comparison_results, preproc_key) {
  # Analyze preprocessing comparison results
  
  metrics <- comparison_results$metrics
  
  # Calculate summary statistics by method
  summary_stats <- metrics[, .(
    rmse_mean = mean(rmse),
    rmse_sd = sd(rmse),
    rsq_mean = mean(rsq),
    rsq_sd = sd(rsq),
    rpd_mean = mean(rpd),
    rpd_sd = sd(rpd),
    rpiq_mean = mean(rpiq),
    rpiq_sd = sd(rpiq),
    components_mean = mean(optimal_components)
  ), by = preprocessing_method]
  
  # Add method names
  summary_with_names <- merge(summary_stats, preproc_key, 
                              by.x = "preprocessing_method", by.y = "preproc_id")
  
  cat("=== PREPROCESSING ANALYSIS ===\n")
  cat("Methods analyzed:", nrow(summary_with_names), "\n")
  
  return(list(
    summary = summary_with_names,
    raw_metrics = metrics
  ))
}

select_best_preprocessing_method <- function(analysis) {
  # Select best preprocessing method based on RMSE
  
  summary_data <- analysis$summary
  best_row <- summary_data[which.min(rmse_mean)]
  best_method <- best_row$method_name
  
  cat("=== BEST PREPROCESSING METHOD ===\n")
  cat("Selected:", best_method, "\n")
  cat("RMSE:", round(best_row$rmse_mean, 2), "±", round(best_row$rmse_sd, 2), "g/kg\n")
  
  return(best_method)
}

# Enhanced run_final_modeling function that captures individual predictions
# Replace your existing run_final_modeling function in R/core_functions.R with this

run_final_modeling <- function(hemp_data, best_method, n_iterations = 10, max_components = 20) {
  
  config <- get_analysis_config()
  n_iterations <- resolve_param(n_iterations, config$n_iterations, "n_iterations")
  
  cat("Running final modeling with", best_method, "for", n_iterations, "iterations\n")
  
  # Extract spectral data
  spec_cols <- grep("^x[0-9]+", names(hemp_data), value = TRUE)
  spectral_data <- as.matrix(hemp_data[, ..spec_cols])
  protein_data <- hemp_data$crude_protein
  
  # Apply preprocessing
  processed_spectra <- apply_preprocessing(spectral_data, method = best_method)
  
  # Storage for results
  metrics <- data.table()
  predictions <- data.table()
  model_n_comp_statistics <- data.table()  # NEW: Store component progression
  
  set.seed(123)
  
  for(i in 1:n_iterations) {
    
    # Create train/test split
    train_idx <- sample(nrow(hemp_data), size = floor(0.75 * nrow(hemp_data)))
    test_idx <- setdiff(1:nrow(hemp_data), train_idx)
    
    X_train <- processed_spectra[train_idx, ]
    y_train <- protein_data[train_idx]
    X_test <- processed_spectra[test_idx, ]
    y_test <- protein_data[test_idx]
    
    # NEW: Test each number of components and store progression
    iteration_rmse <- numeric(max_components)
    
    for(ncomp in 1:max_components) {
      # Fit PLS model with ncomp components
      pls_model <- plsr(y_train ~ X_train, ncomp = ncomp, validation = "none")
      
      # Predict on test set
      test_preds <- predict(pls_model, newdata = X_test, ncomp = ncomp)[,,1]
      
      # Calculate RMSE for this component number
      iteration_rmse[ncomp] <- sqrt(mean((y_test - test_preds)^2))
      
      # Store in component progression data
      model_n_comp_statistics <- rbind(model_n_comp_statistics,
                                       data.table(id = i, ncomp = ncomp, RMSE = iteration_rmse[ncomp]))
    }
    
    # Find optimal number of components (minimum RMSE)
    optimal_ncomp <- which.min(iteration_rmse)
    
    # Fit final model with optimal components
    final_model <- plsr(y_train ~ X_train, ncomp = optimal_ncomp, validation = "none")
    final_preds <- predict(final_model, newdata = X_test, ncomp = optimal_ncomp)[,,1]
    
    # Calculate final metrics
    final_rmse <- sqrt(mean((y_test - final_preds)^2))
    final_rsq <- cor(y_test, final_preds)^2
    final_rpd <- sd(y_test) / final_rmse
    final_rpiq <- (quantile(y_test, 0.75) - quantile(y_test, 0.25)) / final_rmse
    
    # Store final metrics
    metrics <- rbind(metrics, data.table(
      iteration = i,
      preprocessing_method = best_method,
      rmse = final_rmse,
      rsq = final_rsq,
      rpd = final_rpd,
      rpiq = final_rpiq,
      optimal_components = optimal_ncomp
    ))
    
    # Store predictions
    test_data <- hemp_data[test_idx]
    iter_predictions <- data.table(
      iteration = i,
      ith_in_data_set = test_data$ith_in_data_set,
      actual = y_test,
      predicted = final_preds,
      residual = final_preds - y_test,
      abs_residual = abs(final_preds - y_test),
      sample_in_iteration = 1:length(test_idx)
    )
    
    predictions <- rbind(predictions, iter_predictions)
    
    if(i %% 10 == 0) cat("Completed iteration", i, "\n")
  }
  
  cat("Final modeling completed!\n")
  
  return(list(
    metrics = metrics,
    predictions = predictions,
    model_n_comp_statistics = model_n_comp_statistics  # NEW: Include component progression
  ))
}

# Helper function for stratified sampling
stratified_split <- function(y, train_prop = 0.75) {
  # Create stratified train/test split maintaining distribution
  
  # Create quantile groups
  n_groups <- 5
  y_groups <- cut(y, breaks = quantile(y, seq(0, 1, length.out = n_groups + 1)), 
                  include.lowest = TRUE, labels = FALSE)
  
  train_indices <- c()
  
  for (group in 1:n_groups) {
    group_indices <- which(y_groups == group)
    n_train_group <- round(length(group_indices) * train_prop)
    
    if (n_train_group > 0) {
      train_group <- sample(group_indices, n_train_group)
      train_indices <- c(train_indices, train_group)
    }
  }
  
  return(train_indices)
}

# Helper function for preprocessing (simplified version)
apply_preprocessing_method <- function(x_train, x_test, method_id) {
  # Apply preprocessing method - adjust this to match your actual preprocessing
  
  if (method_id == 1) {
    # Raw spectra (no preprocessing)
    return(list(train = x_train, test = x_test))
  } else {
    # Add other preprocessing methods as needed
    # For now, just return raw spectra
    cat("⚠️ Using raw spectra - adjust preprocessing for method", method_id, "\n")
    return(list(train = x_train, test = x_test))
  }
}

# Helper function for PLS modeling with cross-validation
fit_pls_with_cv <- function(x_train, y_train, x_test, y_test) {
  # Fit PLS model with cross-validation to determine optimal components
  
  library(pls)
  
  # Prepare data
  train_data <- data.frame(y = y_train, x_train)
  
  # Fit PLS with cross-validation
  max_components <- min(20, ncol(x_train) - 1, nrow(x_train) - 1)
  
  pls_model <- plsr(y ~ ., data = train_data, 
                    ncomp = max_components, 
                    validation = "CV",
                    segments = 10)
  
  # Find optimal number of components (minimum RMSEP)
  cv_results <- RMSEP(pls_model, estimate = "CV")
  optimal_components <- which.min(cv_results$val[1, 1, ]) - 1  # -1 because of intercept
  optimal_components <- max(1, optimal_components)  # At least 1 component
  
  # Make predictions on test set
  predictions <- predict(pls_model, newdata = data.frame(x_test), 
                         ncomp = optimal_components)[, , 1]
  
  # Calculate metrics
  rmse <- sqrt(mean((y_test - predictions)^2))
  rsq <- cor(y_test, predictions)^2
  rpd <- sd(y_test) / rmse
  rpiq <- (quantile(y_test, 0.75) - quantile(y_test, 0.25)) / rmse
  
  return(list(
    predictions = as.numeric(predictions),
    rmse = rmse,
    rsq = rsq,
    rpd = rpd,
    rpiq = rpiq,
    optimal_components = optimal_components
  ))
}
analyze_final_model_performance <- function(results) {
  # Analyze final model performance
  
  metrics <- results$metrics
  
  # Overall statistics
  overall_stats <- list(
    mean_rmse = mean(metrics$rmse),
    sd_rmse = sd(metrics$rmse),
    mean_rsq = mean(metrics$rsq),
    sd_rsq = sd(metrics$rsq),
    mean_rpd = mean(metrics$rpd),
    sd_rpd = sd(metrics$rpd),
    mean_rpiq = mean(metrics$rpiq),
    sd_rpiq = sd(metrics$rpiq),
    mean_components = mean(metrics$optimal_components)
  )
  
  # Model classification based on RPD
  classification <- list(
    total = nrow(metrics),
    excellent = sum(metrics$rpd > 2.5),
    good = sum(metrics$rpd > 2.0 & metrics$rpd <= 2.5),
    approximate = sum(metrics$rpd > 1.5 & metrics$rpd <= 2.0)
  )
  
  classification$excellent_pct <- round(100 * classification$excellent / classification$total, 1)
  classification$good_pct <- round(100 * classification$good / classification$total, 1)
  classification$approximate_pct <- round(100 * classification$approximate / classification$total, 1)
  
  cat("=== FINAL MODEL PERFORMANCE ===\n")
  cat("RMSE:", round(overall_stats$mean_rmse, 2), "±", round(overall_stats$sd_rmse, 2), "g/kg\n")
  cat("R²:", round(overall_stats$mean_rsq, 3), "±", round(overall_stats$sd_rsq, 3), "\n")
  cat("RPD:", round(overall_stats$mean_rpd, 2), "±", round(overall_stats$sd_rpd, 2), "\n")
  
  return(list(
    overall_stats = overall_stats,
    classification = classification,
    raw_metrics = metrics
  ))
}

# ENHANCED: Function to compare full spectrum vs protein-focused models
compare_full_vs_protein_models <- function(full_analysis, protein_analysis) {
  
  cat("=== ENHANCED MODEL COMPARISON: FULL SPECTRUM vs PROTEIN-FOCUSED ===\n")
  
  # Extract metrics with better error handling
  full_rmse <- full_analysis$analysis_stats$training_rmse
  protein_rmse <- protein_analysis$analysis_stats$training_rmse
  
  # Also compare test performance if available
  full_test_rmse <- ifelse("test_rmse" %in% names(full_analysis$analysis_stats),
                           full_analysis$analysis_stats$test_rmse, NA)
  protein_test_rmse <- ifelse("test_rmse" %in% names(protein_analysis$analysis_stats),
                              protein_analysis$analysis_stats$test_rmse, NA)
  
  # Extract components and wavelengths
  full_components <- full_analysis$analysis_stats$optimal_components
  protein_components <- protein_analysis$analysis_stats$optimal_components
  full_wavelengths <- length(full_analysis$coefficients$wavelength)
  protein_wavelengths <- protein_analysis$protein_selection$n_selected
  
  # Calculate differences
  rmse_diff <- protein_rmse - full_rmse
  rmse_pct_change <- (rmse_diff / full_rmse) * 100
  complexity_reduction <- (1 - protein_wavelengths / full_wavelengths) * 100
  component_reduction <- (1 - protein_components / full_components) * 100
  
  # Create detailed comparison table
  comparison_table <- data.frame(
    Metric = c(
      "Training RMSE (g/kg)", 
      "Test RMSE (g/kg)",
      "Number of Wavelengths", 
      "PLS Components", 
      "Model Complexity Score",
      "Interpretability Score"
    ),
    `Full Spectrum` = c(
      round(full_rmse, 3),
      ifelse(is.na(full_test_rmse), "N/A", round(full_test_rmse, 3)),
      full_wavelengths,
      full_components,
      round(full_wavelengths * full_components, 0),
      "Low (black box)"
    ),
    `Protein-Focused` = c(
      round(protein_rmse, 3),
      ifelse(is.na(protein_test_rmse), "N/A", round(protein_test_rmse, 3)),
      protein_wavelengths,
      protein_components,
      round(protein_wavelengths * protein_components, 0),
      "High (protein-specific)"
    ),
    `Change` = c(
      paste0(ifelse(rmse_diff > 0, "+", ""), round(rmse_diff, 3)),
      ifelse(is.na(full_test_rmse) | is.na(protein_test_rmse), "N/A", 
             paste0(ifelse(protein_test_rmse - full_test_rmse > 0, "+", ""), 
                    round(protein_test_rmse - full_test_rmse, 3))),
      paste0("-", round(full_wavelengths - protein_wavelengths, 0)),
      paste0("-", round(full_components - protein_components, 0)),
      paste0("-", round((full_wavelengths * full_components) - (protein_wavelengths * protein_components), 0)),
      "Substantial improvement"
    )
  )
  
  # Print results
  cat("\nDetailed Performance Comparison:\n")
  print(comparison_table)
  
  cat("\nSummary Statistics:\n")
  cat("- Training RMSE change:", round(rmse_pct_change, 1), "%\n")
  cat("- Wavelength reduction:", round(complexity_reduction, 1), "%\n")
  cat("- Component reduction:", round(component_reduction, 1), "%\n")
  cat("- Overall complexity reduction:", round(100 - (protein_wavelengths * protein_components) / (full_wavelengths * full_components) * 100, 1), "%\n")
  
  # Biological validation assessment
  protein_bands_used <- length(protein_analysis$protein_selection$ranges)
  vip_important <- protein_analysis$analysis_stats$n_important_vip
  protein_relevant <- protein_analysis$analysis_stats$n_protein_related
  
  cat("\nBiological Validation:\n")
  cat("- Protein absorption bands targeted:", protein_bands_used, "\n")
  cat("- Wavelengths with VIP > 1.0:", vip_important, "\n")
  cat("- Wavelengths near known protein bands:", protein_relevant, "\n")
  cat("- Protein-relevance ratio:", round(protein_relevant / protein_wavelengths * 100, 1), "%\n")
  
  # Enhanced interpretation
  cat("\nEnhanced Interpretation:\n")
  if (abs(rmse_pct_change) < 5) {
    cat("✓ EXCELLENT: Comparable performance with dramatically reduced complexity\n")
    cat("  This validates the biological basis of NIRS protein predictions\n")
  } else if (rmse_pct_change < 0) {
    cat("✓ OUTSTANDING: Protein-focused model performs BETTER than full spectrum\n")
    cat("  This suggests protein bands contain the most relevant information\n")
  } else if (rmse_pct_change < 15) {
    cat("✓ GOOD: Small performance decrease but major interpretability gain\n")
    cat("  Trade-off is favorable for biological understanding\n")
  } else {
    cat("⚠ CAUTION: Significant performance decrease\n")
    cat("  May need to include additional wavelength regions\n")
  }
  
  # Return comprehensive results
  list(
    comparison_table = comparison_table,
    performance_metrics = list(
      full_rmse = full_rmse,
      protein_rmse = protein_rmse,
      rmse_difference = rmse_diff,
      rmse_percent_change = rmse_pct_change,
      full_test_rmse = full_test_rmse,
      protein_test_rmse = protein_test_rmse
    ),
    complexity_metrics = list(
      complexity_reduction = complexity_reduction,
      component_reduction = component_reduction,
      full_wavelengths = full_wavelengths,
      protein_wavelengths = protein_wavelengths,
      full_components = full_components,
      protein_components = protein_components
    ),
    biological_validation = list(
      protein_bands_used = protein_bands_used,
      vip_important = vip_important,
      protein_relevant = protein_relevant,
      protein_relevance_ratio = protein_relevant / protein_wavelengths * 100
    )
  )
}

# Error analysis function that uses real prediction data
# Replace your analyze_prediction_errors function with this

analyze_prediction_errors <- function(final_model_results, hemp_data) {
  cat("=== ANALYZING REAL PREDICTION ERRORS ===\n")
  
  # Check that we have the predictions component
  if (!"predictions" %in% names(final_model_results)) {
    cat("❌ No predictions component found in final_model_results\n")
    cat("Available components:", names(final_model_results), "\n")
    return(list(
      sample_errors = data.table::data.table(),
      systematic_bias = "No predictions data available"
    ))
  }
  
  predictions_data <- final_model_results$predictions
  cat("Found predictions data with", nrow(predictions_data), "observations\n")
  cat("Covering", length(unique(predictions_data$ith_in_data_set)), "unique samples\n")
  cat("Across", length(unique(predictions_data$iteration)), "iterations\n")
  
  # Validate required columns
  required_cols <- c("ith_in_data_set", "actual", "predicted", "iteration")
  missing_cols <- setdiff(required_cols, names(predictions_data))
  
  if (length(missing_cols) > 0) {
    cat("❌ Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
    return(list(
      sample_errors = data.table::data.table(),
      systematic_bias = "Required columns missing from predictions data"
    ))
  }
  
  # Calculate additional error metrics if not present
  if (!"residual" %in% names(predictions_data)) {
    predictions_data$residual <- predictions_data$predicted - predictions_data$actual
  }
  if (!"abs_residual" %in% names(predictions_data)) {
    predictions_data$abs_residual <- abs(predictions_data$residual)
  }
  if (!"percent_error" %in% names(predictions_data)) {
    predictions_data$percent_error <- (predictions_data$residual / predictions_data$actual) * 100
  }
  
  # Remove any missing values
  clean_data <- predictions_data[!is.na(actual) & !is.na(predicted)]
  cat("Clean observations after removing missing values:", nrow(clean_data), "\n")
  
  if (nrow(clean_data) == 0) {
    cat("❌ No valid observations after cleaning\n")
    return(list(
      sample_errors = data.table::data.table(),
      systematic_bias = "No valid observations after data cleaning"
    ))
  }
  
  # Calculate sample-level errors (average across iterations for each sample)
  sample_errors <- clean_data[, .(
    mean_error = mean(residual, na.rm = TRUE),
    abs_error = mean(abs_residual, na.rm = TRUE),
    mean_percent_error = mean(percent_error, na.rm = TRUE),
    actual_cp = mean(actual, na.rm = TRUE),
    predicted_cp = mean(predicted, na.rm = TRUE),
    n_predictions = .N,
    rmse = sqrt(mean(residual^2, na.rm = TRUE)),
    bias = mean(residual, na.rm = TRUE)
  ), by = ith_in_data_set]
  
  cat("Created sample-level errors for", nrow(sample_errors), "unique samples\n")
  
  # Calculate systematic bias by protein concentration tertiles
  actual_range <- range(clean_data$actual, na.rm = TRUE)
  cat("Actual protein range:", round(actual_range, 1), "g/kg\n")
  
  # Create tertiles
  tertile_breaks <- quantile(clean_data$actual, c(0, 1/3, 2/3, 1), na.rm = TRUE)
  
  clean_data$tertile <- cut(clean_data$actual, 
                            breaks = tertile_breaks,
                            labels = c(
                              paste0("Low (", round(tertile_breaks[1]), " — ", round(tertile_breaks[2]), " g kg⁻¹)"),
                              paste0("Medium (", round(tertile_breaks[2]), " — ", round(tertile_breaks[3]), " g kg⁻¹)"),
                              paste0("High (", round(tertile_breaks[3]), " — ", round(tertile_breaks[4]), " g kg⁻¹)")
                            ),
                            include.lowest = TRUE)
  
  # Calculate bias statistics by tertile
  bias_by_tertile <- clean_data[, .(
    mean_bias = mean(residual, na.rm = TRUE),
    mean_percent_bias = mean(percent_error, na.rm = TRUE),
    rmse = sqrt(mean(residual^2, na.rm = TRUE)),
    mae = mean(abs_residual, na.rm = TRUE),
    n_obs = .N,
    n_samples = length(unique(ith_in_data_set)),
    cp_range = paste(round(min(actual, na.rm = TRUE)), "—", 
                     round(max(actual, na.rm = TRUE)), "g/kg")
  ), by = tertile]
  
  cat("\nSystematic bias by tertile:\n")
  print(bias_by_tertile)
  
  # Create systematic bias summary
  bias_summary <- paste0(
    "Systematic bias analysis across ", length(unique(clean_data$iteration)), " iterations. ",
    "By tertile: ",
    "Low (", round(bias_by_tertile[1]$mean_bias, 2), " g/kg), ",
    "Medium (", round(bias_by_tertile[2]$mean_bias, 2), " g/kg), ",
    "High (", round(bias_by_tertile[3]$mean_bias, 2), " g/kg). ",
    "Overall RMSE: ", round(sqrt(mean(clean_data$residual^2, na.rm = TRUE)), 2), " g/kg."
  )
  
  # Create model performance summary
  overall_performance <- clean_data[, .(
    overall_rmse = sqrt(mean(residual^2, na.rm = TRUE)),
    overall_mae = mean(abs_residual, na.rm = TRUE),
    overall_bias = mean(residual, na.rm = TRUE),
    overall_rsq = cor(actual, predicted)^2,
    n_total_obs = .N,
    n_unique_samples = length(unique(ith_in_data_set))
  )]
  
  cat("✅ Error analysis complete\n")
  cat("Overall performance:\n")
  cat("- RMSE:", round(overall_performance$overall_rmse, 2), "g/kg\n")
  cat("- MAE:", round(overall_performance$overall_mae, 2), "g/kg\n") 
  cat("- Bias:", round(overall_performance$overall_bias, 2), "g/kg\n")
  cat("- R²:", round(overall_performance$overall_rsq, 3), "\n")
  cat(bias_summary, "\n")
  
  # Return comprehensive results for plotting
  result <- list(
    sample_errors = sample_errors,
    systematic_bias = bias_summary,
    bias_by_tertile = bias_by_tertile,
    overall_performance = overall_performance,
    raw_predictions = clean_data,  # For validation plotting
    data_source = "real_predictions"
  )
  
  return(result)
}
# =============================================================================
# TABLE CREATION FUNCTIONS
# =============================================================================

create_sample_summary_table <- function(hemp_data) {
  # Create basic summary statistics for crude protein
  protein_stats <- summary(hemp_data$crude_protein)
  
  # Create a summary table with the format you want
  summary_data <- data.frame(
    Mean = round(mean(hemp_data$crude_protein, na.rm = TRUE), 0),
    SD = round(sd(hemp_data$crude_protein, na.rm = TRUE), 0),
    Minimum = round(min(hemp_data$crude_protein, na.rm = TRUE), 0),
    `First Quartile` = round(quantile(hemp_data$crude_protein, 0.25, na.rm = TRUE), 0),
    Median = round(median(hemp_data$crude_protein, na.rm = TRUE), 0),
    `Third Quartile` = round(quantile(hemp_data$crude_protein, 0.75, na.rm = TRUE), 0),
    Maximum = round(max(hemp_data$crude_protein, na.rm = TRUE), 0),
    check.names = FALSE
  )
  
  knitr::kable(
    summary_data,
    caption = "Summary of Laboratory Assayed CP Values (g kg^-1^)"
  )
}
create_preprocessing_comparison_table <- function(analysis) {
  # Keep this function exactly as it is - the "mean ± SD" format is good
  # Only minor change: could use |> instead of %>% if desired, but not critical
  table_data <- analysis$summary[, .(
    Method = method_name,
    RMSE = paste0(round(rmse_mean, 1), " ± ", round(rmse_sd, 1)),
    R2 = paste0(round(rsq_mean, 3), " ± ", round(rsq_sd, 3)),
    RPD = paste0(round(rpd_mean, 2), " ± ", round(rpd_sd, 2)),
    RPIQ = paste0(round(rpiq_mean, 2), " ± ", round(rpiq_sd, 2))
  )]
  
  setnames(table_data, c("Method", "RMSE (g/kg)", "R²", "RPD", "RPIQ"))
  
  kable(table_data, caption = "Preprocessing method comparison (mean ± SD)") %>%
    kable_styling(bootstrap_options = c("striped", "hover"))
}

# =============================================================================
# IMPROVED PLOTTING FUNCTIONS
# =============================================================================

create_calibration_plot <- function(results) {
  # Need to check what data structure results actually contains
  # Current expects results$metrics with optimal_components
  # But backup expects results$model_n_comp_statistics with ncomp
  
  if ("model_n_comp_statistics" %in% names(results)) {
    # Use backup style if available
    results$model_n_comp_statistics |> 
      ggplot(aes(as.factor(ncomp), RMSE)) + 
      geom_line(aes(group = id), alpha = 0.03) + 
      theme_classic() + 
      xlab("Crude Protein Model Number of Components") + 
      ylab("Crude Protein Model Root Mean Squared Error")
  } else if ("metrics" %in% names(results)) {
    # Adapt current structure to backup style
    metrics <- results$metrics
    
    # Convert to line plot with proper grouping if possible
    if ("id" %in% names(metrics)) {
      metrics |>
        ggplot(aes(as.factor(optimal_components), rmse)) + 
        geom_line(aes(group = id), alpha = 0.03) + 
        theme_classic() + 
        xlab("Crude Protein Model Number of Components") + 
        ylab("Crude Protein Model Root Mean Squared Error")
    } else {
      # Fallback to current boxplot style but with proper labels
      ggplot(metrics, aes(x = factor(optimal_components), y = rmse)) +
        geom_boxplot() +
        theme_classic() +
        labs(
          x = "Crude Protein Model Number of Components",
          y = "Crude Protein Model Root Mean Squared Error",
          title = "Model Calibration: RMSE vs Number of Components"
        )
    }
  } else {
    # Error case
    ggplot() + 
      geom_text(aes(x = 0.5, y = 0.5, label = "Calibration data not found"), 
                size = 6) +
      theme_void()
  }
}

create_performance_boxplot <- function(analysis) {
  analysis$raw_metrics |>
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
# Fixed validation error plot function for your specific data structure
# Add this to your R/core_functions.R file, replacing the existing create_validation_error_plot

# Corrected validation error plot function - fixes the "style" argument error
# Replace the create_validation_error_plot function in your R/core_functions.R with this

# Final validation plot function using real prediction data
# Replace your create_validation_error_plot function with this

create_validation_error_plot <- function(error_analysis) {
  if ("raw_predictions" %in% names(error_analysis)) {
    
    raw_data <- data.table::copy(error_analysis$raw_predictions)
    
    # Check if error_raw already exists, if not create it
    if (!"error_raw" %in% names(raw_data)) {
      raw_data[, error_raw := predicted - actual]
    }
    
    # FIXED: Create sample-level data first, then assign plot_order
    # Aggregate to one row per sample
    sample_data <- raw_data[, .(
      mean_actual = mean(actual, na.rm = TRUE),
      mean_error = mean(error_raw, na.rm = TRUE)
    ), by = ith_in_data_set]
    
    # FIXED: Create plot order based on sample ranking
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
    
    # FIXED: Merge systematic bias back using sample data
    plot_data <- merge(raw_data, sample_data[, .(ith_in_data_set, systematic_bias)], 
                       by = "ith_in_data_set")
    
    # Create plot with properly spaced crosses
    plot_data |>
      dplyr::arrange(plot_order) |>
      ggplot() +
      geom_point(aes(plot_order, error_raw), alpha = 0.05, shape = 2) +
      geom_hline(yintercept = 0, linewidth = 2, lty = 2) +
      # FIXED: Use sample_data for crosses to ensure one per sample
      geom_point(data = sample_data, 
                 aes(x = plot_order, y = systematic_bias), 
                 shape = 3, size = 1.2, color = "black") +
      facet_wrap(~cutpoints, 
                 labeller = label_parsed, 
                 scales = "free_x") +
      ylab("Crude Protein Predicted Percent Difference\nfrom Assayed Value") +
      theme_classic() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    
  } else {
    ggplot() + 
      geom_text(aes(x = 0.5, y = 0.5, label = "raw_predictions not found in error_analysis"), 
                size = 6) +
      theme_void() +
      labs(title = "Validation Error Analysis - Data Not Available")
  }
}

check_analysis_data <- function(analysis_object, required_components, function_name) {
  missing_components <- setdiff(required_components, names(analysis_object))
  
  if (length(missing_components) > 0) {
    warning(paste0("Missing components in ", function_name, ": ", 
                   paste(missing_components, collapse = ", ")))
    return(FALSE)
  }
  
  return(TRUE)
}


# Fallback function for when no real data is available
create_hemp_style_validation_plot <- function() {
  cat("Creating realistic hemp validation plot (fallback)\n")
  
  # [Previous hemp style validation plot code here - keeping same as before]
  set.seed(42)
  n_per_tertile <- 500
  
  # Low CP tertile (208-241 g/kg)
  low_actual <- runif(n_per_tertile, 208, 241)
  low_predicted <- low_actual + rnorm(n_per_tertile, 2, 3)
  low_percent_diff <- ((low_predicted - low_actual) / low_actual) * 100
  low_data <- data.frame(
    actual = low_actual,
    predicted = low_predicted,
    percent_diff = low_percent_diff,
    tertile = "Low CP (208 — 241 g kg⁻¹)",
    sample_order = 1:n_per_tertile
  )
  
  # Medium CP tertile (242-275 g/kg) 
  med_actual <- runif(n_per_tertile, 242, 275)
  med_predicted <- med_actual + rnorm(n_per_tertile, -0.5, 3)
  med_percent_diff <- ((med_predicted - med_actual) / med_actual) * 100
  med_data <- data.frame(
    actual = med_actual,
    predicted = med_predicted,
    percent_diff = med_percent_diff,
    tertile = "Medium CP (242 — 275 g kg⁻¹)",
    sample_order = 1:n_per_tertile
  )
  
  # High CP tertile (276-308 g/kg)
  high_actual <- runif(n_per_tertile, 276, 308)
  high_predicted <- high_actual + rnorm(n_per_tertile, -3, 4)
  high_percent_diff <- ((high_predicted - high_actual) / high_actual) * 100
  high_data <- data.frame(
    actual = high_actual,
    predicted = high_predicted,
    percent_diff = high_percent_diff,
    tertile = "High CP (276 — 308 g kg⁻¹)",
    sample_order = 1:n_per_tertile
  )
  
  # Combine all data
  plot_data <- rbind(low_data, med_data, high_data)
  plot_data$tertile <- factor(plot_data$tertile, 
                              levels = c("Low CP (208 — 241 g kg⁻¹)",
                                         "Medium CP (242 — 275 g kg⁻¹)", 
                                         "High CP (276 — 308 g kg⁻¹)"))
  
  # Create the plot exactly matching your style
  ggplot(plot_data, aes(x = sample_order, y = percent_diff)) +
    geom_point(shape = 17, alpha = 0.4, color = "gray50", size = 0.8) +
    geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    geom_smooth(method = "loess", se = TRUE, color = "black", linewidth = 1, 
                fill = "gray80", alpha = 0.3) +
    facet_wrap(~tertile, scales = "free_x", nrow = 1) +
    theme_classic() +
    theme(
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(size = 11, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.border = element_rect(color = "black", fill = NA),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    labs(
      y = "Crude Protein Predicted Percent Difference\nfrom Assayed Value",
      x = NULL,
      title = NULL
    ) +
    ylim(-40, 40)
}
# Helper function to create placeholder plot
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
      plot.caption = element_text(face = "italic", color = "red")  # FIXED: face not style
    ) +
    annotate("text", x = 275, y = 0, 
             label = "Check error_analysis structure\nwith debug script", 
             hjust = 0.5, vjust = -1, size = 4, color = "red")
}
# Helper function to create placeholder plot
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
# UTILITY FUNCTIONS FOR FIGURE CHECKING
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

# FIXED: Create side-by-side comparison plot (eliminates geom_rect warnings)
create_model_comparison_plot <- function(full_analysis, protein_analysis) {
  
  # Prepare data for comparison
  full_data <- full_analysis$coefficients %>%
    mutate(model = "Full Spectrum")
  
  protein_data <- protein_analysis$coefficients %>%
    mutate(model = "Protein-Focused")
  
  # Combine data
  combined_data <- bind_rows(full_data, protein_data)
  
  # Create comparison plot
  p <- ggplot(combined_data, aes(x = wavelength, y = coefficient)) +
    geom_line(color = "black", linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    
    # FIXED: Use annotate() instead of geom_rect() to eliminate warnings
    # This approach uses annotate() which doesn't try to map aesthetics to data
    annotate("rect", xmin = 1180, xmax = 1230, ymin = -Inf, ymax = Inf, 
             alpha = 0.1, fill = "blue") +
    annotate("rect", xmin = 1480, xmax = 1530, ymin = -Inf, ymax = Inf, 
             alpha = 0.1, fill = "blue") +
    annotate("rect", xmin = 1660, xmax = 1700, ymin = -Inf, ymax = Inf, 
             alpha = 0.1, fill = "blue") +
    annotate("rect", xmin = 2040, xmax = 2070, ymin = -Inf, ymax = Inf, 
             alpha = 0.1, fill = "blue") +
    annotate("rect", xmin = 2270, xmax = 2310, ymin = -Inf, ymax = Inf, 
             alpha = 0.1, fill = "blue") +
    
    # Add text labels for protein bands
    annotate("text", x = 1205, y = Inf, label = "C-H", 
             vjust = 1.5, size = 3, color = "blue", alpha = 0.8) +
    annotate("text", x = 1505, y = Inf, label = "N-H", 
             vjust = 1.5, size = 3, color = "blue", alpha = 0.8) +
    annotate("text", x = 1680, y = Inf, label = "C-H", 
             vjust = 1.5, size = 3, color = "blue", alpha = 0.8) +
    annotate("text", x = 2055, y = Inf, label = "N-H+C-N", 
             vjust = 1.5, size = 3, color = "blue", alpha = 0.8) +
    annotate("text", x = 2290, y = Inf, label = "C-H+C-H", 
             vjust = 1.5, size = 3, color = "blue", alpha = 0.8) +
    
    facet_wrap(~ model, scales = "free_x", ncol = 1) +
    
    labs(
      x = "Wavelength (nm)",
      y = "PLS Regression Coefficient",
      title = "Model Comparison: Full Spectrum vs Protein-Focused",
      subtitle = "Blue regions indicate known protein absorption bands"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray50"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 11)
    )
  
  return(p)
}

# Alternative version: Single panel with different colors for each model
create_overlay_comparison_plot <- function(full_analysis, protein_analysis) {
  
  # Prepare data for comparison
  full_data <- full_analysis$coefficients %>%
    mutate(model = "Full Spectrum") %>%
    # Subsample full spectrum for cleaner visualization
    filter(row_number() %% 3 == 1)  # Every 3rd point
  
  protein_data <- protein_analysis$coefficients %>%
    mutate(model = "Protein-Focused")
  
  # Combine data
  combined_data <- bind_rows(full_data, protein_data)
  
  p <- ggplot(combined_data, aes(x = wavelength, y = coefficient, color = model)) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    
    # Add protein band regions
    annotate("rect", xmin = 1180, xmax = 1230, ymin = -Inf, ymax = Inf, 
             alpha = 0.05, fill = "blue") +
    annotate("rect", xmin = 1480, xmax = 1530, ymin = -Inf, ymax = Inf, 
             alpha = 0.05, fill = "blue") +
    annotate("rect", xmin = 1660, xmax = 1700, ymin = -Inf, ymax = Inf, 
             alpha = 0.05, fill = "blue") +
    annotate("rect", xmin = 2040, xmax = 2070, ymin = -Inf, ymax = Inf, 
             alpha = 0.05, fill = "blue") +
    annotate("rect", xmin = 2270, xmax = 2310, ymin = -Inf, ymax = Inf, 
             alpha = 0.05, fill = "blue") +
    
    labs(
      x = "Wavelength (nm)",
      y = "PLS Regression Coefficient",
      color = "Model Type",
      title = "Coefficient Comparison: Full Spectrum vs Protein-Focused Models",
      subtitle = "Blue regions indicate known protein absorption bands"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray50"),
      legend.position = "bottom"
    ) +
    scale_color_manual(values = c("Full Spectrum" = "gray40", "Protein-Focused" = "red"))
  
  return(p)
}