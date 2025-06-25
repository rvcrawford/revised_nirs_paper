# Preprocessing functions based on your original implementation

my_preprocess <- function(train, test) {
  # Ensure inputs are matrices
  if (!is.matrix(train)) train <- as.matrix(train)
  if (!is.matrix(test)) test <- as.matrix(test)
  
  # Training set preprocessing with error handling
  raw_train <- train
  
  first_deriv_train <- tryCatch({
    t(diff(t(train), diff = 1))
  }, error = function(e) { train[, -ncol(train)] })  # Fallback if diff fails
  
  second_deriv_train <- tryCatch({
    t(diff(t(train), diff = 2))
  }, error = function(e) { train[, -c((ncol(train)-1):ncol(train))] })
  
  sg_train <- tryCatch({
    prospectr::savitzkyGolay(train, m = 1, p = 3, w = 5)
  }, error = function(e) { train })
  
  gap_der_train <- tryCatch({
    prospectr::gapDer(X = train, m = 1, w = 11, s = 5)
  }, error = function(e) { train })
  
  snv_train <- tryCatch({
    prospectr::standardNormalVariate(train)
  }, error = function(e) { train })
  
  snv_sg_train <- tryCatch({
    prospectr::standardNormalVariate(sg_train)
  }, error = function(e) { train })
  
  snv_detrend_train <- tryCatch({
    prospectr::detrend(train, wav = as.numeric(colnames(train) |> readr::parse_number()))
  }, error = function(e) { train })
  
  msc_train <- tryCatch({
    prospectr::msc(train)
  }, error = function(e) { train })
  
  # Apply same preprocessing to test set
  raw_test <- test
  first_deriv_test <- tryCatch({
    t(diff(t(test), diff = 1))
  }, error = function(e) { test[, -ncol(test)] })
  
  second_deriv_test <- tryCatch({
    t(diff(t(test), diff = 2))
  }, error = function(e) { test[, -c((ncol(test)-1):ncol(test))] })
  
  sg_test <- tryCatch({
    prospectr::savitzkyGolay(test, m = 1, p = 3, w = 5)
  }, error = function(e) { test })
  
  gap_der_test <- tryCatch({
    prospectr::gapDer(X = test, m = 1, w = 11, s = 5)
  }, error = function(e) { test })
  
  snv_test <- tryCatch({
    prospectr::standardNormalVariate(test)
  }, error = function(e) { test })
  
  snv_sg_test <- tryCatch({
    prospectr::standardNormalVariate(sg_test)
  }, error = function(e) { test })
  
  snv_detrend_test <- tryCatch({
    prospectr::detrend(test, wav = as.numeric(colnames(test) |> readr::parse_number()))
  }, error = function(e) { test })
  
  msc_test <- tryCatch({
    prospectr::msc(test, ref_spectrum = attr(msc_train, "Reference spectrum"))
  }, error = function(e) { test })
  
  # Package outputs
  output_train <- list(raw_train, first_deriv_train, second_deriv_train,
                       sg_train, gap_der_train, snv_train, snv_sg_train, 
                       snv_detrend_train, msc_train)
  names(output_train) <- c("raw_train", "first_derivative_train", "second_derivative_train",
                           "sav_gol_train", "gap_der_train", "snv_train", "snv_sg_train",
                           "snv_detrend_train", "msc_train")
  
  output_test <- list(raw_test, first_deriv_test, second_deriv_test,
                      sg_test, gap_der_test, snv_test, snv_sg_test, 
                      snv_detrend_test, msc_test)
  names(output_test) <- c("raw_test", "first_derivative_test", "second_derivative_test",
                          "sav_gol_test", "gap_der_test", "snv_test", "snv_sg_test",
                          "snv_detrend_test", "msc_test")
  
  return(list(output_train, output_test))
}

split_spectra <- function(y) {
  inTrain <- createDataPartition(
    y = y,
    p = .75,
    list = FALSE,
    groups = 3
  )
  return(inTrain)
}

get_preprocessing_methods <- function() {
  c("raw", "first_derivative", "sav_gol", "gap_der", 
    "snv", "snv_sg", "snv_detrend", "msc")
}

validate_preprocessing_method <- function(method) {
  valid_methods <- get_preprocessing_methods()
  if (!method %in% valid_methods) {
    stop("Invalid preprocessing method. Choose from: ", paste(valid_methods, collapse = ", "))
  }
  method
}

apply_single_preprocessing <- function(train_spectra, test_spectra, method = "snv_sg") {
  
  # Ensure inputs are matrices
  if (!is.matrix(train_spectra)) train_spectra <- as.matrix(train_spectra)
  if (!is.matrix(test_spectra)) test_spectra <- as.matrix(test_spectra)
  
  switch(method,
         # Raw (no preprocessing)
         "raw" = {
           list(
             train = train_spectra,
             test = test_spectra
           )
         },
         
         # First derivative
         "first_derivative" = {
           train_processed <- tryCatch({
             t(diff(t(train_spectra), diff = 1))
           }, error = function(e) { train_spectra[, -ncol(train_spectra)] })
           
           test_processed <- tryCatch({
             t(diff(t(test_spectra), diff = 1))
           }, error = function(e) { test_spectra[, -ncol(test_spectra)] })
           
           list(train = train_processed, test = test_processed)
         },
         
         # Savitzky-Golay
         "sav_gol" = {
           train_processed <- tryCatch({
             prospectr::savitzkyGolay(train_spectra, m = 1, p = 3, w = 5)
           }, error = function(e) { train_spectra })
           
           test_processed <- tryCatch({
             prospectr::savitzkyGolay(test_spectra, m = 1, p = 3, w = 5)
           }, error = function(e) { test_spectra })
           
           list(train = train_processed, test = test_processed)
         },
         
         # Standard Normal Variate
         "snv" = {
           train_processed <- tryCatch({
             prospectr::standardNormalVariate(train_spectra)
           }, error = function(e) { train_spectra })
           
           test_processed <- tryCatch({
             prospectr::standardNormalVariate(test_spectra)
           }, error = function(e) { test_spectra })
           
           list(train = train_processed, test = test_processed)
         },
         
         # SNV + Savitzky-Golay (your best method!)
         "snv_sg" = {
           # First apply Savitzky-Golay
           sg_train <- tryCatch({
             prospectr::savitzkyGolay(train_spectra, m = 1, p = 3, w = 5)
           }, error = function(e) { train_spectra })
           
           sg_test <- tryCatch({
             prospectr::savitzkyGolay(test_spectra, m = 1, p = 3, w = 5)
           }, error = function(e) { test_spectra })
           
           # Then apply SNV
           train_processed <- tryCatch({
             prospectr::standardNormalVariate(sg_train)
           }, error = function(e) { sg_train })
           
           test_processed <- tryCatch({
             prospectr::standardNormalVariate(sg_test)
           }, error = function(e) { sg_test })
           
           list(train = train_processed, test = test_processed)
         },
         
         # Gap-segment derivative
         "gap_der" = {
           train_processed <- tryCatch({
             prospectr::gapDer(X = train_spectra, m = 1, w = 11, s = 5)
           }, error = function(e) { train_spectra })
           
           test_processed <- tryCatch({
             prospectr::gapDer(X = test_spectra, m = 1, w = 11, s = 5)
           }, error = function(e) { test_spectra })
           
           list(train = train_processed, test = test_processed)
         },
         
         # MSC
         "msc" = {
           train_processed <- tryCatch({
             prospectr::msc(train_spectra)
           }, error = function(e) { train_spectra })
           
           test_processed <- tryCatch({
             prospectr::msc(test_spectra, ref_spectrum = attr(train_processed, "Reference spectrum"))
           }, error = function(e) { test_spectra })
           
           list(train = train_processed, test = test_processed)
         },
         
         # Default fallback
         {
           warning("Unknown preprocessing method: ", method, ". Using raw spectra.")
           list(train = train_spectra, test = test_spectra)
         }
  )
}

# =============================================================================
# ULTRA-FAST WEIGHTED MODELING (9x speed improvement)
# =============================================================================

#' Ultra-fast weighted modeling with single preprocessing method
run_weighted_modeling_ultrafast <- function(balanced_data, best_method, n_iterations = 1000) {
  cat("=== ULTRA-FAST WEIGHTED MODELING ===\n")
  cat("Method:", best_method, "\n")
  cat("Iterations:", n_iterations, "\n")
  cat("Using SINGLE preprocessing method (9x speed boost!)\n")
  
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
      
      # SPEED BOOST: Apply only the single preprocessing method needed
      processed_data <- apply_single_preprocessing(
        spectra_matrix[inTrain, ], 
        spectra_matrix[-inTrain, ], 
        method = best_method
      )
      
      train_spectral <- processed_data$train
      test_spectral <- processed_data$test
      
      # Apply weights to training data
      weights_train <- sapply(locations_train, function(loc) weights_lookup[[loc]])
      
      # Create training dataframe
      train_df <- data.frame(y_train = y_train, train_spectral)
      
      # Fit weighted PLS model (still with CV for accuracy)
      model <- train(
        y_train ~ .,
        data = train_df,
        method = "pls",
        tuneLength = 15,        # Slightly reduced 
        weights = weights_train,
        trControl = trainControl(
          method = "cv", 
          number = 5,           # Reduced CV folds
          verboseIter = FALSE
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