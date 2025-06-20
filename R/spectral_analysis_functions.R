# ============================================================================
# SPECTRAL FEATURE ANALYSIS FUNCTIONS
# Save this as: spectral_analysis_functions.R
# ============================================================================

library(pls)
library(tidyverse)
library(caret)

# 1. Extract wavelength information from your spectral data
# ============================================================================
# REPLACE THESE TWO FUNCTIONS in your R/spectral_analysis_functions.R
# ============================================================================

# FIXED: Extract wavelength information 
extract_wavelengths <- function(spectral_data) {
  cat("Extracting wavelength information...\n")
  
  # Your spectral columns are named like "x1100", "x1102", etc.
  spectral_cols <- names(spectral_data)[grepl("^x[0-9]+$", names(spectral_data))]
  
  if (length(spectral_cols) == 0) {
    cat("ERROR: No spectral columns found with pattern '^x[0-9]+$'\n")
    cat("Available columns (first 10):", paste(head(names(spectral_data), 10), collapse = ", "), "\n")
    stop("No spectral columns detected")
  }
  
  wavelengths <- as.numeric(gsub("x", "", spectral_cols))
  
  cat("- Found", length(spectral_cols), "spectral columns\n")
  cat("- Wavelength range:", range(wavelengths), "nm\n")
  
  list(
    wavelengths = wavelengths,
    spectral_cols = spectral_cols,
    n_wavelengths = length(wavelengths),
    range = range(wavelengths)
  )
}
# 2. Fit a PLS model specifically for spectral analysis
fit_pls_for_spectral_analysis <- function(data, best_method = "snv_sg") {
  cat("Fitting PLS model for spectral analysis using", best_method, "\n")
  
  # Extract spectral data with better error checking
  wavelength_info <- extract_wavelengths(data)
  
  # FIXED: Use proper data.table column selection
  if (inherits(data, "data.table")) {
    # For data.table, use this syntax
    spectral_matrix <- as.matrix(data[, wavelength_info$spectral_cols, with = FALSE])
  } else {
    # For regular data.frame
    spectral_matrix <- as.matrix(data[, wavelength_info$spectral_cols])
  }
  
  y <- data$crude_protein
  
  cat("Spectral matrix dimensions:", dim(spectral_matrix), "\n")
  cat("Response variable range:", range(y, na.rm = TRUE), "\n")
  
  # Verify dimensions are correct (should be samples x wavelengths)
  expected_samples <- length(y)
  expected_wavelengths <- length(wavelength_info$wavelengths)
  
  if (nrow(spectral_matrix) != expected_samples) {
    stop("Dimension ERROR: spectral matrix has ", nrow(spectral_matrix), 
         " rows but expected ", expected_samples, " samples")
  }
  
  if (ncol(spectral_matrix) != expected_wavelengths) {
    stop("Dimension ERROR: spectral matrix has ", ncol(spectral_matrix), 
         " columns but expected ", expected_wavelengths, " wavelengths")
  }
  
  cat("✓ Dimensions correct:", nrow(spectral_matrix), "samples x", ncol(spectral_matrix), "wavelengths\n")
  
  # Create train/test split (consistent with your existing approach)
  set.seed(123)  # For reproducibility of wavelength analysis
  train_indices <- split_spectra(y)
  
  cat("Train/test split: ", length(train_indices), "/", nrow(spectral_matrix) - length(train_indices), "\n")
  
  # Apply preprocessing using your existing function
  processed_data <- my_preprocess(
    spectral_matrix[train_indices, ], 
    spectral_matrix[-train_indices, ]
  )
  
  # Get the best method's training data
  method_name <- paste0(best_method, "_train")
  if (!method_name %in% names(processed_data[[1]])) {
    cat("Available methods:", names(processed_data[[1]]), "\n")
    stop("Method ", method_name, " not found. Available: ", 
         paste(names(processed_data[[1]]), collapse = ", "))
  }
  
  X_train <- processed_data[[1]][[method_name]]
  y_train <- y[train_indices]
  
  cat("Processed training data dimensions:", dim(X_train), "\n")
  
  # Verify processed dimensions make sense
  if (nrow(X_train) != length(y_train)) {
    stop("After preprocessing: X_train has ", nrow(X_train), 
         " rows but y_train has ", length(y_train), " values")
  }
  
  if (ncol(X_train) == 0) {
    stop("After preprocessing: X_train has 0 columns - preprocessing failed")
  }
  
  # CRITICAL FIX: Ensure column names are preserved
  if (is.null(colnames(X_train))) {
    cat("Adding missing column names...\n")
    n_cols <- ncol(X_train)
    if (n_cols == length(wavelength_info$wavelengths)) {
      colnames(X_train) <- paste0("w", wavelength_info$wavelengths)
    } else {
      # If preprocessing changed the number of columns (e.g., derivatives)
      colnames(X_train) <- paste0("w", 1:n_cols)
    }
  }
  
  # Convert to data frame for caret (ensures proper format)
  X_train_df <- as.data.frame(X_train)
  
  cat("Final training data for PLS:\n")
  cat("- Dimensions:", dim(X_train_df), "\n")
  cat("- Column names:", head(names(X_train_df), 3), "...", tail(names(X_train_df), 3), "\n")
  cat("- Any missing values:", any(is.na(X_train_df)), "\n")
  
  # Fit PLS model using caret (consistent with your approach)
  cat("Fitting PLS model...\n")
  model <- train(
    x = X_train_df,  # Use data frame with proper column names
    y = y_train,
    method = "pls",
    tuneLength = 20,
    trControl = trainControl(method = "cv", number = 10),
    metric = "RMSE"
  )
  
  cat("PLS model fitted successfully with", model$bestTune$ncomp, "components\n")
  
  list(
    model = model,
    wavelengths = wavelength_info$wavelengths,
    X_train = X_train_df,
    y_train = y_train,
    train_indices = train_indices,
    preprocessing_method = best_method
  )
}
# 3. Extract PLS regression coefficients

# ALSO UPDATE the extract_pls_coefficients function to handle column name issues
extract_pls_coefficients <- function(spectral_fit) {
  model <- spectral_fit$model
  wavelengths <- spectral_fit$wavelengths
  
  # Get optimal number of components
  optimal_ncomp <- model$bestTune$ncomp
  
  # Extract coefficients from the final PLS model
  pls_object <- model$finalModel
  coeffs <- coef(pls_object, ncomp = optimal_ncomp, intercept = FALSE)
  
  cat("Coefficient extraction:\n")
  cat("- Optimal components:", optimal_ncomp, "\n")
  cat("- Coefficient matrix dimensions:", dim(coeffs), "\n")
  
  # Handle case where preprocessing changed the number of wavelengths
  n_coeffs <- length(as.vector(coeffs))
  n_wavelengths <- length(wavelengths)
  
  if (n_coeffs == n_wavelengths) {
    # Perfect match - use original wavelengths
    used_wavelengths <- wavelengths
  } else {
    cat("Warning: Number of coefficients (", n_coeffs, ") != number of wavelengths (", n_wavelengths, ")\n")
    cat("This is normal for derivative preprocessing. Creating sequential wavelength labels.\n")
    # Create sequential wavelengths for the available coefficients
    used_wavelengths <- wavelengths[1:n_coeffs]
  }
  
  # Create coefficient data frame
  coeff_data <- data.frame(
    wavelength = used_wavelengths,
    coefficient = as.vector(coeffs),
    abs_coefficient = abs(as.vector(coeffs))
  )
  
  # Add ranking by importance
  coeff_data$rank <- rank(-coeff_data$abs_coefficient)
  
  cat("Extracted coefficients for", optimal_ncomp, "components\n")
  cat("Coefficient range:", round(range(coeff_data$coefficient), 4), "\n")
  cat("Using", nrow(coeff_data), "wavelength points\n")
  
  return(coeff_data)
}

# 4. Calculate Variable Importance in Projection (VIP) scores
  # ============================================================================
  # REPLACE the calculate_vip_scores function in your R/spectral_analysis_functions.R
  # ============================================================================
  
  # FIXED: Calculate Variable Importance in Projection (VIP) scores
  calculate_vip_scores <- function(spectral_fit) {
    model <- spectral_fit$model
    
    # Extract PLS components
    pls_object <- model$finalModel
    optimal_ncomp <- model$bestTune$ncomp
    
    cat("VIP calculation:\n")
    cat("- Optimal components:", optimal_ncomp, "\n")
    
    # Calculate VIP scores manually (since caret doesn't export VIP function)
    # VIP formula: sqrt(p * sum((w_a * (SS_a / SS_total))^2) / A)
    
    W <- pls_object$loading.weights[, 1:optimal_ncomp, drop = FALSE]
    T <- pls_object$scores[, 1:optimal_ncomp, drop = FALSE]
    Y <- spectral_fit$y_train
    
    cat("- Loading weights dimensions:", dim(W), "\n")
    cat("- Scores dimensions:", dim(T), "\n")
    
    # Calculate sum of squares for each component
    SS <- colSums(T^2)
    SS_total <- sum(SS)
    
    # Calculate VIP scores
    p <- nrow(W)
    vip_scores <- sqrt(p * rowSums((W^2) %*% diag(SS / SS_total)) / optimal_ncomp)
    
    cat("- VIP scores calculated for", length(vip_scores), "variables\n")
    
    # CRITICAL FIX: Handle wavelength alignment (same as coefficients)
    wavelengths <- spectral_fit$wavelengths
    n_vip <- length(vip_scores)
    n_wavelengths <- length(wavelengths)
    
    cat("- Original wavelengths:", n_wavelengths, "\n")
    cat("- VIP scores calculated:", n_vip, "\n")
    
    if (n_vip == n_wavelengths) {
      # Perfect match - use original wavelengths
      used_wavelengths <- wavelengths
      cat("- Using original wavelengths (perfect match)\n")
    } else {
      cat("- Wavelength mismatch detected (normal for derivative preprocessing)\n")
      # Create sequential wavelengths for the available VIP scores
      # Use the same logic as coefficients to ensure consistency
      used_wavelengths <- wavelengths[1:n_vip]
      cat("- Using first", n_vip, "wavelengths:", range(used_wavelengths), "\n")
    }
    
    # Ensure we have the same number of wavelengths and VIP scores
    if (length(used_wavelengths) != length(vip_scores)) {
      stop("INTERNAL ERROR: wavelengths (", length(used_wavelengths), 
           ") and VIP scores (", length(vip_scores), ") still don't match")
    }
    
    vip_data <- data.frame(
      wavelength = used_wavelengths,
      vip_score = vip_scores,
      important = vip_scores > 1.0
    )
    
    cat("- Final VIP data dimensions:", dim(vip_data), "\n")
    cat("- Important wavelengths (VIP > 1.0):", sum(vip_data$important), "\n")
    
    return(vip_data)
  }

  # 5. Identify chemically-relevant wavelengths
identify_protein_bands <- function(coeff_data, vip_data) {
  # Known protein absorption bands in NIR
  protein_bands <- data.frame(
    band_center = c(1190, 1215, 1395, 1440, 1500, 1515, 1680, 1725, 1765, 1940, 2055, 2180, 2290, 2350),
    band_range_low = c(1180, 1200, 1380, 1420, 1480, 1500, 1660, 1710, 1750, 1920, 2040, 2160, 2270, 2330),
    band_range_high = c(1200, 1230, 1410, 1460, 1520, 1530, 1700, 1740, 1780, 1960, 2070, 2200, 2310, 2370),
    assignment = c(
      "C-H stretch 2nd overtone",
      "C-H stretch 2nd overtone", 
      "O-H stretch 1st overtone",
      "O-H stretch + C-H deformation",
      "N-H stretch 1st overtone",
      "N-H stretch 1st overtone",
      "C-H stretch 1st overtone",
      "C-H stretch 1st overtone",
      "C-H stretch 1st overtone",
      "O-H stretch + C-H deformation combination",
      "N-H stretch + C-N stretch combination",
      "C-H stretch + C-C stretch combination",
      "C-H stretch + C-H deformation combination",
      "C-H stretch + C-H deformation combination"
    ),
    protein_relevance = c(
      "Amino acid side chains",
      "Amino acid side chains",
      "Hydroxyl groups in amino acids",
      "Mixed amino acid vibrations",
      "Peptide bonds, amino groups",
      "Peptide bonds, amino groups", 
      "Aliphatic amino acids",
      "Aliphatic amino acids",
      "Aliphatic amino acids",
      "Protein-water interactions",
      "Protein backbone structure",
      "Amino acid structural vibrations",
      "Complex amino acid interactions",
      "Complex amino acid interactions"
    )
  )
  
  # Combine coefficient and VIP data
  combined_data <- merge(coeff_data, vip_data, by = "wavelength")
  
  # Map to protein bands
  combined_data$closest_band <- sapply(combined_data$wavelength, function(w) {
    distances <- abs(protein_bands$band_center - w)
    protein_bands$band_center[which.min(distances)]
  })
  
  combined_data$assignment <- sapply(combined_data$wavelength, function(w) {
    distances <- abs(protein_bands$band_center - w)
    protein_bands$assignment[which.min(distances)]
  })
  
  combined_data$protein_relevance <- sapply(combined_data$wavelength, function(w) {
    distances <- abs(protein_bands$band_center - w)
    protein_bands$protein_relevance[which.min(distances)]
  })
  
  combined_data$distance_to_band <- sapply(combined_data$wavelength, function(w) {
    distances <- abs(protein_bands$band_center - w)
    min(distances)
  })
  
  # Flag wavelengths as protein-related if within 20 nm of known bands
  combined_data$likely_protein_related <- combined_data$distance_to_band <= 20
  
  return(combined_data)
}

# 6. Create summary table for manuscript
# ============================================================================
# REPLACE the create_wavelength_summary_table function in your R/spectral_analysis_functions.R
# ============================================================================

# FIXED: Create summary table for manuscript
create_wavelength_summary_table <- function(protein_bands_data, n_top = 15) {
  
  cat("Creating wavelength summary table...\n")
  cat("- Input data dimensions:", dim(protein_bands_data), "\n")
  cat("- Protein-related wavelengths:", sum(protein_bands_data$likely_protein_related, na.rm = TRUE), "\n")
  
  # Get top wavelengths by coefficient magnitude or VIP score
  # Use a more robust approach for the filtering and selection
  
  # First, create the filtered dataset
  candidate_wavelengths <- protein_bands_data %>%
    filter(likely_protein_related == TRUE | rank <= 10) %>%
    arrange(rank)
  
  cat("- Candidate wavelengths after filtering:", nrow(candidate_wavelengths), "\n")
  
  # Take the top n_top
  if (nrow(candidate_wavelengths) == 0) {
    cat("WARNING: No wavelengths found meeting criteria. Using top 10 by rank.\n")
    candidate_wavelengths <- protein_bands_data %>%
      arrange(rank) %>%
      slice_head(n = min(10, nrow(protein_bands_data)))
  }
  
  top_wavelengths <- candidate_wavelengths %>%
    slice_head(n = n_top)
  
  cat("- Final wavelengths selected:", nrow(top_wavelengths), "\n")
  
  # Create the formatted table with simpler column creation
  # Avoid the select() issue by creating columns step by step
  formatted_table <- top_wavelengths %>%
    mutate(
      wavelength_nm = wavelength,
      pls_coefficient = sprintf("%.4f", coefficient),
      vip_score_formatted = sprintf("%.2f", vip_score),
      chemical_assignment = assignment,
      protein_relevance_text = protein_relevance,
      distance_to_band = round(distance_to_band, 0)
    )
  
  # Create final data frame with desired column names
  final_table <- data.frame(
    "Wavelength (nm)" = formatted_table$wavelength_nm,
    "PLS Coefficient" = formatted_table$pls_coefficient,
    "VIP Score" = formatted_table$vip_score_formatted,
    "Chemical Assignment" = formatted_table$chemical_assignment,
    "Protein Relevance" = formatted_table$protein_relevance_text,
    "Distance to Known Band (nm)" = formatted_table$distance_to_band,
    check.names = FALSE  # Keep the column names with spaces
  )
  
  cat("- Final table dimensions:", dim(final_table), "\n")
  cat("- Column names:", names(final_table), "\n")
  
  return(final_table)
}


# 7. Main spectral analysis function (integrates with your workflow)
run_spectral_analysis <- function(data, best_method = "snv_sg") {
  cat("Running spectral feature analysis...\n")
  
  # Fit PLS model for spectral analysis
  spectral_fit <- fit_pls_for_spectral_analysis(data, best_method)
  
  # Extract coefficients and VIP scores
  coeff_data <- extract_pls_coefficients(spectral_fit)
  vip_data <- calculate_vip_scores(spectral_fit)
  
  # Identify protein-relevant bands
  protein_bands_data <- identify_protein_bands(coeff_data, vip_data)
  
  # Create summary table
  summary_table <- create_wavelength_summary_table(protein_bands_data)
  
  # Calculate some summary statistics
  n_important_vip <- sum(vip_data$vip_score > 1.0)
  n_protein_related <- sum(protein_bands_data$likely_protein_related)
  
  cat("Analysis complete:\n")
  cat("- Wavelengths with VIP > 1.0:", n_important_vip, "\n")
  cat("- Wavelengths near protein bands:", n_protein_related, "\n")
  cat("- Model used", spectral_fit$model$bestTune$ncomp, "components\n")
  
  list(
    spectral_fit = spectral_fit,
    coefficients = coeff_data,
    vip_scores = vip_data,
    protein_bands = protein_bands_data,
    summary_table = summary_table,
    analysis_stats = list(
      n_important_vip = n_important_vip,
      n_protein_related = n_protein_related,
      optimal_components = spectral_fit$model$bestTune$ncomp
    )
  )
}

# ============================================================================
# COPY AND PASTE THESE 3 FUNCTIONS TO THE END OF YOUR R/spectral_analysis_functions.R FILE
# ============================================================================

# Function to select protein-specific wavelengths
select_protein_wavelengths <- function(data) {
  
  cat("Selecting protein-specific wavelengths...\n")
  
  # Define protein absorption bands (nm)
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
  
  # Create corresponding column names
  protein_cols <- paste0("x", protein_wavelengths)
  protein_cols <- protein_cols[protein_cols %in% wavelength_info$spectral_cols]
  
  cat("Summary:\n")
  cat("- Total available wavelengths:", length(all_wavelengths), "\n")
  cat("- Protein-specific wavelengths:", length(protein_wavelengths), "\n")
  cat("- Reduction:", round((1 - length(protein_wavelengths)/length(all_wavelengths)) * 100, 1), "%\n")
  
  list(
    wavelengths = protein_wavelengths,
    column_names = protein_cols,
    assignments = protein_assignments,
    n_selected = length(protein_wavelengths),
    n_total = length(all_wavelengths)
  )
}

# Modified spectral analysis function for protein-only wavelengths
run_protein_focused_analysis <- function(data, best_method = "snv_sg") {
  
  cat("Running PROTEIN-FOCUSED spectral analysis...\n")
  
  # Select protein wavelengths
  protein_selection <- select_protein_wavelengths(data)
  
  if (length(protein_selection$column_names) < 10) {
    stop("Not enough protein wavelengths found (", length(protein_selection$column_names), 
         "). Need at least 10 for modeling.")
  }
  
  # Extract only protein-specific spectral data
  if (inherits(data, "data.table")) {
    spectral_matrix <- as.matrix(data[, protein_selection$column_names, with = FALSE])
  } else {
    spectral_matrix <- as.matrix(data[, protein_selection$column_names])
  }
  
  y <- data$crude_protein
  
  cat("Protein-focused matrix dimensions:", dim(spectral_matrix), "\n")
  cat("Response variable range:", range(y, na.rm = TRUE), "\n")
  
  # Apply same preprocessing as full model
  set.seed(123)  # Same seed for fair comparison
  train_indices <- split_spectra(y)
  
  processed_data <- my_preprocess(
    spectral_matrix[train_indices, ], 
    spectral_matrix[-train_indices, ]
  )
  
  # Get processed training data
  method_name <- paste0(best_method, "_train")
  X_train <- processed_data[[1]][[method_name]]
  y_train <- y[train_indices]
  
  cat("Processed protein-focused data dimensions:", dim(X_train), "\n")
  
  # Add column names if missing
  if (is.null(colnames(X_train))) {
    n_cols <- ncol(X_train)
    if (n_cols == length(protein_selection$wavelengths)) {
      colnames(X_train) <- paste0("w", protein_selection$wavelengths)
    } else {
      colnames(X_train) <- paste0("w", 1:n_cols)
    }
  }
  
  X_train_df <- as.data.frame(X_train)
  
  # Fit PLS model
  cat("Fitting protein-focused PLS model...\n")
  model <- train(
    x = X_train_df,
    y = y_train,
    method = "pls",
    tuneLength = 20,
    trControl = trainControl(method = "cv", number = 10),
    metric = "RMSE"
  )
  
  cat("Protein-focused model fitted with", model$bestTune$ncomp, "components\n")
  cat("Training RMSE:", round(min(model$results$RMSE), 3), "\n")
  
  # Extract coefficients and VIP scores (same functions work)
  spectral_fit <- list(
    model = model,
    wavelengths = protein_selection$wavelengths,
    X_train = X_train_df,
    y_train = y_train,
    train_indices = train_indices,
    preprocessing_method = best_method
  )
  
  coeff_data <- extract_pls_coefficients(spectral_fit)
  vip_data <- calculate_vip_scores(spectral_fit)
  protein_bands_data <- identify_protein_bands(coeff_data, vip_data)
  summary_table <- create_wavelength_summary_table(protein_bands_data)
  
  # Calculate summary statistics
  n_important_vip <- sum(vip_data$vip_score > 1.0)
  n_protein_related <- sum(protein_bands_data$likely_protein_related)
  
  cat("Protein-focused analysis complete:\n")
  cat("- Wavelengths with VIP > 1.0:", n_important_vip, "\n")
  cat("- Wavelengths near protein bands:", n_protein_related, "\n")
  cat("- Model used", model$bestTune$ncomp, "components\n")
  
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
      training_rmse = min(model$results$RMSE)
    ),
    protein_selection = protein_selection
  )
}

# Function to compare full spectrum vs protein-focused models
# ============================================================================
# REPLACE the compare_full_vs_protein_models function in your R/spectral_analysis_functions.R
# ============================================================================

# FIXED: Function to compare full spectrum vs protein-focused models
compare_full_vs_protein_models <- function(full_analysis, protein_analysis) {
  
  cat("=== MODEL COMPARISON: FULL SPECTRUM vs PROTEIN-FOCUSED ===\n")
  
  # Extract key metrics with error handling
  # For full analysis, extract RMSE from the model results
  if ("training_rmse" %in% names(full_analysis$analysis_stats)) {
    full_rmse <- full_analysis$analysis_stats$training_rmse
  } else {
    # Extract from the spectral_fit model
    full_rmse <- min(full_analysis$spectral_fit$model$results$RMSE, na.rm = TRUE)
  }
  
  # For protein analysis, should have training_rmse
  if ("training_rmse" %in% names(protein_analysis$analysis_stats)) {
    protein_rmse <- protein_analysis$analysis_stats$training_rmse
  } else {
    protein_rmse <- min(protein_analysis$spectral_fit$model$results$RMSE, na.rm = TRUE)
  }
  
  # Check if RMSE values are numeric
  if (!is.numeric(full_rmse) || !is.numeric(protein_rmse) || 
      is.na(full_rmse) || is.na(protein_rmse)) {
    cat("ERROR: Could not extract valid RMSE values\n")
    cat("Full RMSE:", full_rmse, "class:", class(full_rmse), "\n")
    cat("Protein RMSE:", protein_rmse, "class:", class(protein_rmse), "\n")
    stop("Invalid RMSE values for comparison")
  }
  
  # Extract components
  full_components <- full_analysis$analysis_stats$optimal_components
  protein_components <- protein_analysis$analysis_stats$optimal_components
  
  # Extract wavelength counts
  full_wavelengths <- length(full_analysis$coefficients$wavelength)
  protein_wavelengths <- protein_analysis$protein_selection$n_selected
  
  # Calculate differences
  rmse_diff <- protein_rmse - full_rmse
  rmse_pct_change <- (rmse_diff / full_rmse) * 100
  complexity_reduction <- (1 - protein_wavelengths / full_wavelengths) * 100
  
  cat("\nPerformance Comparison:\n")
  cat("- Full spectrum RMSE:", round(full_rmse, 3), "\n")
  cat("- Protein-focused RMSE:", round(protein_rmse, 3), "\n")
  cat("- Difference:", round(rmse_diff, 3), "(", round(rmse_pct_change, 1), "%)\n")
  
  cat("\nModel Complexity:\n")
  cat("- Full spectrum wavelengths:", full_wavelengths, "\n")
  cat("- Protein-focused wavelengths:", protein_wavelengths, "\n")
  cat("- Complexity reduction:", round(complexity_reduction, 1), "%\n")
  
  cat("\nComponents:\n")
  cat("- Full spectrum components:", full_components, "\n")
  cat("- Protein-focused components:", protein_components, "\n")
  
  # Interpretation
  cat("\nInterpretation:\n")
  if (abs(rmse_pct_change) < 5) {
    cat("✓ Similar performance with much simpler model\n")
  } else if (rmse_pct_change < 0) {
    cat("✓ Protein-focused model performs BETTER!\n")
  } else if (rmse_pct_change < 10) {
    cat("○ Slight performance decrease but much better interpretability\n")
  } else {
    cat("⚠ Significant performance decrease - may need different approach\n")
  }
  
  list(
    full_rmse = full_rmse,
    protein_rmse = protein_rmse,
    rmse_difference = rmse_diff,
    rmse_percent_change = rmse_pct_change,
    complexity_reduction = complexity_reduction,
    full_wavelengths = full_wavelengths,
    protein_wavelengths = protein_wavelengths,
    full_components = full_components,
    protein_components = protein_components
  )
}