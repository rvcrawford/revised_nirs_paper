map_location_codes <- function(loc_code) {
  # Based on user clarification - should give us 76 Ithaca, 41 Geneva, 24 Freeville
  location_mapping <- c(
    # Ithaca variants (76 total)
    "ithaca" = "ithaca",      # 28 samples
    "eith" = "ithaca",        # 13 samples  
    "mg4e" = "ithaca",        # 13 samples
    "mg4l" = "ithaca",        # 13 samples
    "mcg" = "ithaca",         # 9 samples
    # Total: 76 samples ✓
    
    # Geneva variants (41 total)
    "rsgeneva" = "geneva",    # 2 samples
    "rn041_gen" = "geneva",   # 8 samples
    "rn028_gen" = "geneva",   # 2 samples  
    "cnoll" = "geneva",       # 13 samples
    "rnooa" = "geneva",       # 13 samples
    "cn north" = "geneva",    # 3 samples
    # Total: 41 samples ✓
    
    # Freeville variants (24 total)
    "freev" = "freeville",    # 13 samples
    "free" = "freeville",     # 11 samples
    # Total: 24 samples ✓
    
    # Tiny locations
    "chazy" = "chazy",        # 5 samples
    "willsboro" = "willsboro" # 3 samples
  )
  
  mapped <- location_mapping[loc_code]
  mapped[is.na(mapped)] <- "unknown"
  return(mapped)
}

#' Filter data to main locations and add class weights
#' Filter data to main locations and add class weights (FIXED)
prepare_location_balanced_data <- function(data) {
  cat("=== LOCATION-BALANCED DATA PREPARATION ===\n")
  
  # Show original distribution
  location_counts <- table(data$loc)
  cat("Original location distribution:\n")
  print(location_counts)
  
  # Map location codes to clean names
  data[, clean_loc := map_location_codes(loc)]
  
  # Show mapped distribution
  clean_location_counts <- table(data$clean_loc)
  cat("\nMapped location distribution:\n")
  print(clean_location_counts)
  
  # Define main locations (sufficient samples for training)
  main_locations <- c("ithaca", "geneva", "freeville")
  tiny_locations <- c("chazy", "willsboro")
  
  # Filter to main locations only
  main_data <- data[clean_loc %in% main_locations]
  tiny_data <- data[clean_loc %in% tiny_locations]
  
  cat("\nFiltered to main locations:\n")
  main_counts <- table(main_data$clean_loc)
  print(main_counts)
  
  # FIXED: Calculate class weights using simple ratio method
  # Use ithaca as baseline (weight = 1.0)
  ithaca_count <- sum(main_data$clean_loc == "ithaca")
  
  weights_list <- list()
  for (loc in main_locations) {
    count <- sum(main_data$clean_loc == loc)
    if (count > 0) {
      weight <- ithaca_count / count  # Simple ratio: larger weights for smaller classes
      weights_list[[loc]] <- weight
      cat("Location", loc, ":", count, "samples, weight =", round(weight, 3), "\n")
    } else {
      cat("WARNING: Location", loc, "has 0 samples after filtering\n")
    }
  }
  
  # Add weight column to main data using data.table syntax
  main_data[, sample_weight := weights_list[[clean_loc]], by = clean_loc]
  
  # Create clean location summary
  location_summary <- data.frame(
    location = names(main_counts),
    count = as.numeric(main_counts),
    weight = sapply(names(main_counts), function(x) weights_list[[x]]),
    stringsAsFactors = FALSE
  )
  
  # Return structured results
  list(
    main_data = main_data,
    tiny_data = tiny_data,
    weights = weights_list,
    location_summary = location_summary
  )
}
#' Stratified train/test split for locations
split_spectra_stratified <- function(y, locations = NULL, p = 0.75) {
  if (is.null(locations)) {
    # Original behavior
    inTrain <- createDataPartition(y = y, p = p, list = FALSE, groups = 3)
  } else {
    # Stratified by location
    cat("Using location-stratified split\n")
    inTrain <- createDataPartition(y = locations, p = p, list = FALSE)
  }
  return(inTrain)
}

#' Run weighted modeling analysis (fast version for testing)
run_weighted_modeling_fast <- function(balanced_data, best_method, n_iterations = 100) {
  cat("=== FAST WEIGHTED MODELING ANALYSIS ===\n")
  cat("Method:", best_method, "\n")
  cat("Iterations:", n_iterations, "\n")
  
  main_data <- balanced_data$main_data
  
  # Extract spectral data
  spectral_cols <- grep("^x[0-9]+$", names(main_data), value = TRUE)
  spectra_matrix <- as.matrix(main_data[, ..spectral_cols])
  y <- main_data$crude_protein
  locations <- main_data$clean_loc
  weights_lookup <- balanced_data$weights
  
  cat("Data dimensions:", dim(spectra_matrix), "\n")
  
  # Storage for results
  results_list <- list()
  predictions_list <- list()
  
  for (i in 1:n_iterations) {
    if (i %% 25 == 0) cat("Iteration", i, "\n")
    
    tryCatch({
      # Stratified split
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
      
      # Apply weights to training data
      weights_train <- sapply(locations_train, function(loc) weights_lookup[[loc]])
      
      # Create training dataframe
      train_df <- data.frame(y_train = y_train, train_spectral)
      
      # Fit FAST weighted PLS model (reduced CV and components for speed)
      model <- train(
        y_train ~ .,
        data = train_df,
        method = "pls",
        tuneLength = 10,        # Reduced from 20
        weights = weights_train,
        trControl = trainControl(
          method = "cv", 
          number = 5,           # Reduced from 10 folds
          verboseIter = FALSE   # Reduce output
        )
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
        ncomp = model$bestTune$ncomp
      )
      
      # Store predictions with location info
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

#' Run weighted modeling analysis
run_weighted_modeling <- function(balanced_data, best_method, n_iterations = 1000, n_components = NULL) {
  cat("=== WEIGHTED MODELING ANALYSIS ===\n")
  cat("Method:", best_method, "\n")
  cat("Iterations:", n_iterations, "\n")
  
  main_data <- balanced_data$main_data
  
  # Extract spectral data
  spectral_cols <- grep("^x[0-9]+$", names(main_data), value = TRUE)
  spectra_matrix <- as.matrix(main_data[, ..spectral_cols])
  y <- main_data$crude_protein
  locations <- main_data$clean_loc
  weights_lookup <- balanced_data$weights
  
  cat("Data dimensions:", dim(spectra_matrix), "\n")
  
  # Storage for results
  results_list <- list()
  predictions_list <- list()
  
  for (i in 1:n_iterations) {
    if (i %% 100 == 0) cat("Iteration", i, "\n")
    
    tryCatch({
      # Stratified split
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
      
      # Apply weights to training data
      weights_train <- sapply(locations_train, function(loc) weights_lookup[[loc]])
      
      # Create training dataframe
      train_df <- data.frame(y_train = y_train, train_spectral)
      
      # Fit weighted PLS model with optional fixed components
      if (!is.null(n_components)) {
        # Use fixed number of components (faster)
        model <- train(
          y_train ~ .,
          data = train_df,
          method = "pls",
          tuneGrid = data.frame(ncomp = n_components),
          weights = weights_train,
          trControl = trainControl(method = "none")  # Skip CV for speed
        )
      } else {
        # Use cross-validation to find optimal components (slower but thorough)
        model <- train(
          y_train ~ .,
          data = train_df,
          method = "pls",
          tuneLength = 20,
          weights = weights_train,
          trControl = trainControl(method = "cv", number = 10)
        )
      }
      
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
        ncomp = model$bestTune$ncomp
      )
      
      # Store predictions with location info
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

#' Compare weighted vs unweighted approaches
compare_weighted_unweighted <- function(balanced_data, full_data, best_method, n_iterations = 500) {
  cat("=== WEIGHTED vs UNWEIGHTED COMPARISON ===\n")
  
  # Run weighted analysis on balanced data
  cat("Running WEIGHTED analysis...\n")
  weighted_results <- run_weighted_modeling(balanced_data, best_method, n_iterations)
  
  # Run unweighted analysis on balanced data (same locations, no weights)
  cat("Running UNWEIGHTED analysis...\n")
  main_data <- balanced_data$main_data
  main_data[, sample_weight := NULL]  # Remove weights
  
  # Create a temporary version for the existing function (which expects 'loc' column)
  main_data[, loc := clean_loc]  # Use clean location names
  
  # Use your existing function but on filtered data
  unweighted_results <- run_single_method_modeling(main_data, best_method, n_iterations)
  
  # Combine for comparison
  list(
    weighted = weighted_results,
    unweighted = unweighted_results,
    balanced_data_info = balanced_data$location_summary
  )
}

#' Analyze weighted vs unweighted results
#' Fixed analyze_weighting_comparison function
#' Analyze weighted vs unweighted results (FIXED VERSION)
#' Analyze weighted vs unweighted results (ROBUST DATA.TABLE VERSION)
analyze_weighting_comparison <- function(comparison_results) {
  cat("=== ANALYZING WEIGHTING COMPARISON ===\n")
  
  # Load required libraries
  library(data.table)
  library(dplyr)
  
  # Extract metrics and ensure they're data.tables
  weighted_metrics <- as.data.table(comparison_results$weighted$metrics)
  unweighted_metrics <- as.data.table(comparison_results$unweighted$metrics)
  
  cat("Weighted metrics columns:", names(weighted_metrics), "\n")
  cat("Unweighted metrics columns:", names(unweighted_metrics), "\n")
  cat("Weighted dimensions:", dim(weighted_metrics), "\n")
  cat("Unweighted dimensions:", dim(unweighted_metrics), "\n")
  
  # FIX: Handle column name differences using data.table syntax
  if ("ncomp" %in% names(weighted_metrics)) {
    setnames(weighted_metrics, "ncomp", "n_components")
    cat("Renamed ncomp to n_components in weighted data\n")
  }
  
  # Verify required columns exist
  required_cols <- c("iteration", "rmse", "rsq", "rpd", "rpiq")
  missing_weighted <- setdiff(required_cols, names(weighted_metrics))
  missing_unweighted <- setdiff(required_cols, names(unweighted_metrics))
  
  if (length(missing_weighted) > 0) {
    stop("Missing columns in weighted_metrics: ", paste(missing_weighted, collapse = ", "))
  }
  if (length(missing_unweighted) > 0) {
    stop("Missing columns in unweighted_metrics: ", paste(missing_unweighted, collapse = ", "))
  }
  
  cat("All required columns present\n")
  
  # Calculate summary statistics using data.table syntax
  cat("Calculating weighted summary...\n")
  weighted_summary <- weighted_metrics[, .(
    approach = "weighted",
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_rsq = mean(rsq, na.rm = TRUE),
    mean_rpd = mean(rpd, na.rm = TRUE),
    median_rmse = median(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    q75_rmse = quantile(rmse, 0.75, na.rm = TRUE),
    n_excellent = sum(rpd > 3, na.rm = TRUE),
    percent_excellent = round(100 * sum(rpd > 3, na.rm = TRUE) / .N, 1)
  )]
  
  cat("Calculating unweighted summary...\n")
  unweighted_summary <- unweighted_metrics[, .(
    approach = "unweighted",
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_rsq = mean(rsq, na.rm = TRUE),
    mean_rpd = mean(rpd, na.rm = TRUE),
    median_rmse = median(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    q75_rmse = quantile(rmse, 0.75, na.rm = TRUE),
    n_excellent = sum(rpd > 3, na.rm = TRUE),
    percent_excellent = round(100 * sum(rpd > 3, na.rm = TRUE) / .N, 1)
  )]
  
  cat("Weighted summary:\n")
  print(weighted_summary)
  cat("Unweighted summary:\n") 
  print(unweighted_summary)
  
  # Combine summaries using data.table rbind
  cat("Combining summaries...\n")
  summary_comparison <- rbind(weighted_summary, unweighted_summary)
  
  if (nrow(summary_comparison) == 0) {
    stop("Failed to create summary_comparison - rbind resulted in 0 rows")
  }
  
  cat("Summary comparison created with", nrow(summary_comparison), "rows\n")
  cat("Summary comparison columns:", names(summary_comparison), "\n")
  print(summary_comparison)
  
  # Calculate location-specific performance for weighted approach
  cat("Calculating location performance...\n")
  location_performance <- comparison_results$weighted$predictions %>%
    group_by(location) %>%
    summarise(
      n_predictions = n(),
      mean_rmse = sqrt(mean((actual - predicted)^2)),
      mean_rsq = cor(actual, predicted)^2,
      mean_bias = mean(predicted - actual),
      .groups = 'drop'
    )
  
  # Calculate improvements
  rmse_improvement <- weighted_summary$mean_rmse - unweighted_summary$mean_rmse
  rsq_improvement <- weighted_summary$mean_rsq - unweighted_summary$mean_rsq
  rpd_improvement <- weighted_summary$mean_rpd - unweighted_summary$mean_rpd
  
  cat("RMSE improvement:", rmse_improvement, "\n")
  cat("R² improvement:", rsq_improvement, "\n")
  cat("RPD improvement:", rpd_improvement, "\n")
  
  result <- list(
    summary_comparison = summary_comparison,
    location_performance = location_performance,
    improvement = list(
      rmse_improvement = rmse_improvement,
      rsq_improvement = rsq_improvement,
      rpd_improvement = rpd_improvement
    )
  )
  
  cat("Analysis complete - returning result\n")
  return(result)
}

# Append this to the end of your R/location_weighting_functions.R file

#' CORRECTED: Compare weighted vs unweighted with forced differences
compare_weighted_unweighted_corrected <- function(balanced_data, full_data, best_method, n_iterations = 1000) {
  cat("=== CORRECTED WEIGHTED vs UNWEIGHTED COMPARISON ===\n")
  
  library(caret)
  library(data.table)
  
  # 1. Run weighted analysis (use existing function)
  cat("Running WEIGHTED analysis...\n")
  weighted_results <- run_weighted_modeling(balanced_data, best_method, n_iterations)
  
  # 2. Create TRULY unweighted data for comparison
  cat("Preparing UNWEIGHTED data...\n")
  
  # Start with the same balanced data but FORCE no weights
  unweighted_data <- copy(balanced_data$main_data)
  
  # CRITICAL: Remove ALL weight-related columns
  weight_cols <- c("sample_weight", "weight")
  for (col in weight_cols) {
    if (col %in% names(unweighted_data)) {
      unweighted_data[, (col) := NULL]
    }
  }
  
  # Add location column for predictions
  if (!"loc" %in% names(unweighted_data)) {
    unweighted_data[, loc := clean_loc]
  }
  
  cat("Unweighted data prepared - columns:", names(unweighted_data), "\n")
  cat("No weight columns present:", !any(c("sample_weight", "weight") %in% names(unweighted_data)), "\n")
  
  # 3. Run MANUAL unweighted analysis to ensure no weights
  cat("Running MANUAL UNWEIGHTED analysis...\n")
  
  # Extract components manually
  spectral_cols <- grep("^x[0-9]+$", names(unweighted_data), value = TRUE)
  spectra_matrix <- as.matrix(unweighted_data[, ..spectral_cols])
  y <- unweighted_data$crude_protein
  locations <- unweighted_data$loc
  
  cat("Manual unweighted setup:\n")
  cat("- Spectral dimensions:", dim(spectra_matrix), "\n")
  cat("- Response length:", length(y), "\n")
  cat("- Locations:", unique(locations), "\n")
  
  # Storage for unweighted results
  unweighted_results_list <- list()
  unweighted_predictions_list <- list()
  
  for (i in 1:n_iterations) {
    if (i %% 100 == 0) cat("Unweighted iteration", i, "\n")
    
    tryCatch({
      # Use same random seed as weighted for fair comparison
      set.seed(i)
      inTrain <- split_spectra_stratified(y, locations)
      
      y_train <- y[inTrain]
      y_test <- y[-inTrain]
      locations_train <- locations[inTrain]
      locations_test <- locations[-inTrain]
      
      # Apply preprocessing (same as weighted)
      spectra_processed <- my_preprocess(
        spectra_matrix[inTrain, ], 
        spectra_matrix[-inTrain, ]
      )
      
      # Extract method data
      method_train_name <- paste0(best_method, "_train")
      method_test_name <- paste0(best_method, "_test")
      
      train_spectral <- spectra_processed[[1]][[method_train_name]]
      test_spectral <- spectra_processed[[2]][[method_test_name]]
      
      # Create training dataframe - CRITICAL: NO WEIGHTS
      train_df <- data.frame(y_train = y_train, train_spectral)
      
      # Fit PLS model WITHOUT weights - FORCE this to be different
      model <- train(
        y_train ~ .,
        data = train_df,
        method = "pls",
        tuneLength = 20,
        # CRITICAL: NO weights parameter at all
        trControl = trainControl(method = "cv", number = 10)
      )
      
      # Make predictions
      predictions <- predict(model, newdata = as.data.frame(test_spectral))
      
      # Calculate metrics
      rmse_val <- sqrt(mean((y_test - predictions)^2))
      rsq_val <- cor(y_test, predictions)^2
      rpd_val <- sd(y_test) / rmse_val
      rpiq_val <- IQR(y_test) / rmse_val
      
      # Store results
      unweighted_results_list[[i]] <- data.frame(
        iteration = i,
        rmse = rmse_val,
        rsq = rsq_val,
        rpd = rpd_val,
        rpiq = rpiq_val,
        n_components = model$bestTune$ncomp
      )
      
      # Store predictions WITH location info
      unweighted_predictions_list[[i]] <- data.frame(
        iteration = i,
        actual = y_test,
        predicted = predictions,
        location = locations_test,
        sample_id = (1:length(y_test))
      )
      
    }, error = function(e) {
      cat("Error in unweighted iteration", i, ":", e$message, "\n")
    })
  }
  
  # Combine unweighted results
  unweighted_metrics <- do.call(rbind, unweighted_results_list)
  unweighted_predictions <- do.call(rbind, unweighted_predictions_list)
  
  unweighted_results <- list(
    metrics = unweighted_metrics,
    predictions = unweighted_predictions
  )
  
  # 4. Verify they're actually different
  cat("\n=== VERIFICATION ===\n")
  cat("Weighted predictions (first 5):", head(weighted_results$predictions$predicted, 5), "\n")
  cat("Unweighted predictions (first 5):", head(unweighted_results$predictions$predicted, 5), "\n")
  
  identical_check <- identical(
    head(weighted_results$predictions$predicted, 100),
    head(unweighted_results$predictions$predicted, 100)
  )
  cat("Are first 100 predictions identical?", identical_check, "\n")
  
  if (identical_check) {
    stop("ERROR: Predictions are still identical! Bug not fixed.")
  } else {
    cat("SUCCESS: Predictions are different!\n")
  }
  
  # Return comparison
  list(
    weighted = weighted_results,
    unweighted = unweighted_results,
    balanced_data_info = balanced_data$location_summary
  )
}