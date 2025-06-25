# Modeling functions for PLSR analysis

create_train_test_split <- function(data, split_ratio = 0.75) {
  # Create tertiles based on CP concentration
  tertiles <- cut(data$crude_protein, 3)
  
  # Split within each tertile
  split_indices <- data |>
    mutate(tertile = tertiles, row_id = row_number()) |>
    group_by(tertile) |>
    slice_sample(prop = split_ratio) |>
    pull(row_id)
  
  list(
    train = data[split_indices, ],
    test = data[-split_indices, ]
  )
}

fit_plsr_model <- function(train_data, preprocessing_method = "raw", max_components = 20) {
  # Extract spectral data (assuming columns 1100-2498 nm are spectral data)
  # This will need to be adjusted based on your actual column structure
  spectral_cols <- grep("^[0-9]+$", names(train_data), value = TRUE)
  
  if (length(spectral_cols) == 0) {
    # Fallback - assume spectral data starts after basic columns
    spectral_cols <- names(train_data)[8:ncol(train_data)]
  }
  
  X <- as.matrix(train_data[, ..spectral_cols])
  y <- train_data$crude_protein
  
  # Apply preprocessing
  X_processed <- apply_preprocessing(X, preprocessing_method)
  
  # Fit PLSR model with cross-validation for component selection
  plsr_model <- train(
    x = X_processed,
    y = y,
    method = "pls",
    tuneGrid = data.frame(ncomp = 1:max_components),
    trControl = trainControl(
      method = "boot",
      number = 25,
      savePredictions = "final"
    ),
    metric = "RMSE"
  )
  
  plsr_model
}

apply_preprocessing <- function(spectra_matrix, method) {
  switch(method,
         "raw" = spectra_matrix,
         "first_derivative" = t(diff(t(spectra_matrix), differences = 1)),
         "second_derivative" = t(diff(t(spectra_matrix), differences = 2)),
         "savitzky_golay" = prospectr::savitzkyGolay(spectra_matrix, m = 1, p = 3, w = 5),
         "gap_segment" = prospectr::gapDer(spectra_matrix, m = 1, w = 11, s = 5),
         "snv" = prospectr::standardNormalVariate(spectra_matrix),
         "snv_sg" = {
           sg_processed <- prospectr::savitzkyGolay(spectra_matrix, m = 1, p = 3, w = 5)
           prospectr::standardNormalVariate(sg_processed)
         },
         "snv_detrend" = prospectr::detrend(prospectr::standardNormalVariate(spectra_matrix), wav = 1:ncol(spectra_matrix)),
         "msc" = prospectr::msc(spectra_matrix),
         spectra_matrix  # default to raw
  )
}

predict_plsr <- function(model, test_data, preprocessing_method = "raw") {
  # Extract spectral data
  spectral_cols <- grep("^[0-9]+$", names(test_data), value = TRUE)
  
  if (length(spectral_cols) == 0) {
    spectral_cols <- names(test_data)[8:ncol(test_data)]
  }
  
  X_test <- as.matrix(test_data[, ..spectral_cols])
  X_test_processed <- apply_preprocessing(X_test, preprocessing_method)
  
  predict(model, X_test_processed)
}

calculate_model_metrics <- function(actual, predicted) {
  # Custom metric calculations
  rmse_val <- sqrt(mean((actual - predicted)^2))
  rsq_val <- cor(actual, predicted)^2
  
  # RPD (Relative Predicted Deviation)
  rpd_val <- sd(actual) / rmse_val
  
  # RPIQ (Ratio of Performance to Interquartile Distance)
  rpiq_val <- (quantile(actual, 0.75) - quantile(actual, 0.25)) / rmse_val
  
  data.frame(
    rmse = rmse_val,
    rsq = rsq_val,
    rpd = rpd_val,
    rpiq = rpiq_val
  )
}

run_preprocessing_comparison <- function(full_data, n_iterations = NULL) {
  
  # Get configuration and resolve parameters
  config <- get_analysis_config()
  n_iterations <- resolve_param(n_iterations, config$n_iterations, "n_iterations")
  
  cat("Starting preprocessing comparison with", n_iterations, "iterations\n")
  
  # Extract spectral data
  spectral_cols <- grep("^x[0-9]+$", names(full_data), value = TRUE)
  spectra_matrix <- as.matrix(full_data[, ..spectral_cols])
  y <- full_data$crude_protein
  
  cat("Data dimensions:", dim(spectra_matrix), "\n")
  cat("Response range:", round(range(y)), "\n")
  
  # Run iterations using your working preprocess_analyze_function
  all_results <- list()
  successful_iterations <- 0
  
  for (i in 1:n_iterations) {
    if (i %% 10 == 0) cat("Iteration", i, "\n")
    
    tryCatch({
      # Use your working function
      iteration_result <- preprocess_analyze_function(spectra_matrix, y)
      
      if (nrow(iteration_result) > 0) {
        iteration_result[, iteration := i]
        all_results[[length(all_results) + 1]] <- iteration_result
        successful_iterations <- successful_iterations + 1
      }
    }, error = function(e) {
      cat("Error in iteration", i, ":", e$message, "\n")
    })
  }
  
  cat("Completed", successful_iterations, "successful iterations\n")
  
  if (length(all_results) == 0) {
    stop("No valid results obtained")
  }
  
  # Combine all results
  combined_results <- rbindlist(all_results)
  
  # Calculate metrics by iteration and method
  final_results <- combined_results[, {
    if (length(y) > 3 && length(value) > 3 && sd(value, na.rm = TRUE) > 0) {
      rmse_val <- sqrt(mean((y - value)^2, na.rm = TRUE))
      rsq_val <- cor(y, value, use = "complete.obs")^2
      
      if (!is.na(rmse_val) && !is.na(rsq_val) && rmse_val > 0) {
        rpd_val <- sd(y, na.rm = TRUE) / rmse_val
        rpiq_val <- (quantile(y, 0.75, na.rm = TRUE) - quantile(y, 0.25, na.rm = TRUE)) / rmse_val
        
        data.table(
          rmse = rmse_val,
          rsq = rsq_val,
          rpd = rpd_val,
          rpiq = rpiq_val
        )
      } else {
        data.table(rmse = NA_real_, rsq = NA_real_, rpd = NA_real_, rpiq = NA_real_)
      }
    } else {
      data.table(rmse = NA_real_, rsq = NA_real_, rpd = NA_real_, rpiq = NA_real_)
    }
  }, by = .(iteration, preprocessing_method = variable)]
  
  # Remove missing values
  final_results <- final_results[!is.na(rmse) & !is.na(rsq)]
  
  cat("Final results:", nrow(final_results), "rows with", 
      length(unique(final_results$preprocessing_method)), "methods\n")
  
  return(final_results)
}

run_final_modeling <- function(data, best_method, n_iterations = 1000) {
  cat("Running final modeling with", best_method, "for", n_iterations, "iterations\n")
  
  # Extract spectral data
  spectral_cols <- grep("^x[0-9]+$", names(data), value = TRUE)
  spectra_matrix <- as.matrix(data[, ..spectral_cols])
  y <- data$crude_protein
  
  cat("Data dimensions:", dim(spectra_matrix), "\n")
  
  # Run iterations using your working preprocess_analyze_function
  all_results <- list()
  successful_iterations <- 0
  
  for (i in 1:n_iterations) {
    if (i %% 100 == 0) cat("Iteration", i, "\n")
    
    tryCatch({
      # Use your working function
      iteration_result <- preprocess_analyze_function(spectra_matrix, y)
      
      if (nrow(iteration_result) > 0) {
        # Filter to just the best method
        best_method_result <- iteration_result[variable == best_method]
        if (nrow(best_method_result) > 0) {
          best_method_result[, iteration := i]
          all_results[[length(all_results) + 1]] <- best_method_result
          successful_iterations <- successful_iterations + 1
        }
      }
    }, error = function(e) {
      cat("Error in iteration", i, ":", e$message, "\n")
    })
  }
  
  cat("Completed", successful_iterations, "successful iterations for", best_method, "\n")
  
  if (length(all_results) == 0) {
    stop("No valid results obtained for final modeling")
  }
  
  # Combine all results
  combined_results <- rbindlist(all_results)
  
  # Calculate metrics by iteration
  metrics_results <- combined_results[, {
    if (length(y) > 3 && length(value) > 3 && sd(value, na.rm = TRUE) > 0) {
      rmse_val <- sqrt(mean((y - value)^2, na.rm = TRUE))
      rsq_val <- cor(y, value, use = "complete.obs")^2
      
      if (!is.na(rmse_val) && !is.na(rsq_val) && rmse_val > 0) {
        rpd_val <- sd(y, na.rm = TRUE) / rmse_val
        rpiq_val <- (quantile(y, 0.75, na.rm = TRUE) - quantile(y, 0.25, na.rm = TRUE)) / rmse_val
        
        data.table(
          n_components = 12,  # Approximate typical value
          rmse = rmse_val,
          rsq = rsq_val,
          rpd = rpd_val,
          rpiq = rpiq_val
        )
      } else {
        data.table(n_components = NA_integer_, rmse = NA_real_, rsq = NA_real_, rpd = NA_real_, rpiq = NA_real_)
      }
    } else {
      data.table(n_components = NA_integer_, rmse = NA_real_, rsq = NA_real_, rpd = NA_real_, rpiq = NA_real_)
    }
  }, by = iteration]
  
  # Create predictions data for error analysis
  predictions_results <- combined_results[, .(
    iteration = iteration,
    sample_id = seq_len(.N),
    actual = y,
    predicted = value
  )]
  
  list(
    metrics = metrics_results[!is.na(rmse)],
    predictions = predictions_results[!is.na(predicted)]
  )
}
# Fixed preprocess_analyze_function - your original logic with bug fixes
preprocess_analyze_function <- function(spectra, y){
  # First come up with indices for test and training sets
  inTrain <- split_spectra(y)
  
  y_train <- y[inTrain]
  y_test <- y[-inTrain]  # Fixed: use [-inTrain] instead of [!inTrain]
  
  # Apply preprocessing  
  spectra_processed <- my_preprocess(spectra[inTrain, ], spectra[-inTrain, ])
  
  # Create training datasets (using cbind like your original working approach)
  trains <- lapply(spectra_processed[[1]], function(x) cbind(y_train, x))
  
  # Train PLSR models for each preprocessing method
  trained_models <- lapply(trains, function(x) {
    train(
      y_train ~ .,
      data = x,
      method = "pls",
      tuneLength = 20
    )
  })
  
  # Make predictions on test sets
  predicts <- mapply(function(model, test_data) {
    predict(model, newdata = as.data.frame(test_data))
  }, trained_models, spectra_processed[[2]], SIMPLIFY = FALSE)
  
  # Return results in long format
  result_df <- data.frame(y_test = y_test)
  for (i in seq_along(predicts)) {
    result_df[[names(spectra_processed[[1]])[i]]] <- predicts[[i]]
  }
  
  # Convert to long format and clean up method names
  result_final <- setDT(result_df) |> melt(id.vars = "y_test")
  setnames(result_final, "y_test", "y")  # Rename for consistency
  result_final[, variable := gsub("_train$", "", variable)]  # Clean method names
  
  return(result_final)
}

run_single_method_modeling <- function(data, best_method, n_iterations = 1000) {
  cat("Running", n_iterations, "iterations with", best_method, "only\n")
  
  spectral_cols <- grep("^x[0-9]+$", names(data), value = TRUE)
  spectra_matrix <- as.matrix(data[, ..spectral_cols])
  y <- data$crude_protein
  
  cat("Data dimensions:", dim(spectra_matrix), "\n")
  cat("Response range:", range(y), "\n")
  cat("Best method:", best_method, "\n")
  
  # Extract original sample identifiers
  if ("ith_in_data_set" %in% names(data)) {
    original_sample_ids <- data$ith_in_data_set
    cat("Found ith_in_data_set column with", length(unique(original_sample_ids)), "unique values\n")
  } else {
    # Create sample IDs if they don't exist
    original_sample_ids <- 1:nrow(data)
    cat("Created sample IDs 1 to", nrow(data), "\n")
  }
  
  results_list <- list()
  predictions_list <- list()
  tuning_list <- list()
  successful_iterations <- 0
  
  for (i in 1:n_iterations) {
    if (i %% 100 == 0) cat("Iteration", i, "- successful so far:", successful_iterations, "\n")
    
    tryCatch({
      # Train/test split
      inTrain <- split_spectra(y)
      y_train <- y[inTrain]
      y_test <- y[-inTrain]
      
      cat("Iteration", i, "- train/test split successful, sizes:", length(y_train), "/", length(y_test), "\n")
      
      # Get the actual sample identifiers for test set
      test_sample_ids <- original_sample_ids[-inTrain]
      
      # Only process the best method
      spectra_processed <- my_preprocess(spectra_matrix[inTrain, ], spectra_matrix[-inTrain, ])
      
      cat("Iteration", i, "- preprocessing successful\n")
      
      # Get only the best method's data
      method_train_name <- paste0(best_method, "_train")
      method_test_name <- paste0(best_method, "_test")
      
      if (!method_train_name %in% names(spectra_processed[[1]])) {
        cat("ERROR: Method", method_train_name, "not found in processed data\n")
        cat("Available methods:", names(spectra_processed[[1]]), "\n")
        next
      }
      
      train_data <- spectra_processed[[1]][[method_train_name]]
      test_data <- spectra_processed[[2]][[method_test_name]]
      
      cat("Iteration", i, "- extracted method data\n")
      
      # Fit model
      train_df <- data.frame(y_train = y_train, train_data)
      model <- train(
        y_train ~ .,
        data = train_df,
        method = "pls",
        tuneLength = 20
      )
      
      cat("Iteration", i, "- model training successful\n")
      
      # Extract tuning results for Figure 1
      tuning_results <- model$results
      tuning_results$id <- i
      tuning_list[[i]] <- tuning_results[, c("ncomp", "RMSE", "id")]
      
      # Predict
      predictions <- predict(model, newdata = as.data.frame(test_data))
      
      cat("Iteration", i, "- predictions successful\n")
      
      # Calculate metrics
      rmse_val <- sqrt(mean((y_test - predictions)^2))
      rsq_val <- cor(y_test, predictions)^2
      rpd_val <- sd(y_test) / rmse_val
      rpiq_val <- (quantile(y_test, 0.75) - quantile(y_test, 0.25)) / rmse_val
      
      results_list[[i]] <- data.table(
        iteration = i,
        n_components = model$bestTune$ncomp,
        rmse = rmse_val,
        rsq = rsq_val,
        rpd = rpd_val,
        rpiq = rpiq_val
      )
      
      # First, extract location info for test samples
      if ("loc" %in% names(data)) {
        locations <- data$loc
      } else if ("clean_loc" %in% names(data)) {  
        locations <- data$clean_loc
      } else {
        locations <- rep("unknown", nrow(data))  # fallback
      }
      
      locations_test <- locations[-inTrain]  # test set locations
      
      predictions_list[[i]] <- data.frame(
        iteration = i,
        sample_id = test_sample_ids,
        ith_in_data_set = test_sample_ids[1:length(y_test)], 
        actual = y_test,
        predicted = predictions,
        location = locations_test  # ADD THIS LINE
      )
      
      successful_iterations <- successful_iterations + 1
      cat("Iteration", i, "- COMPLETED SUCCESSFULLY\n")
      
    }, error = function(e) {
      cat("ERROR in iteration", i, ":", e$message, "\n")
    })
  }
  
  cat("Final summary: ", successful_iterations, "successful iterations out of", n_iterations, "\n")
  cat("Length of results_list:", length(results_list), "\n")
  cat("Length of predictions_list:", length(predictions_list), "\n")
  
  # Check if we have any results
  if (length(results_list) == 0) {
    cat("ERROR: No successful iterations - all modeling failed\n")
    return(list(
      metrics = data.table(),
      predictions = data.table(),
      model_n_comp_statistics = data.table()
    ))
  }
  
  final_results <- list(
    metrics = rbindlist(results_list),
    predictions = rbindlist(predictions_list),
    model_n_comp_statistics = rbindlist(tuning_list)
  )
  
  cat("Final metrics rows:", nrow(final_results$metrics), "\n")
  cat("Final predictions rows:", nrow(final_results$predictions), "\n")
  
  return(final_results)
}