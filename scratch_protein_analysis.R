# =============================================================================
# PROTEIN-FOCUSED MODEL VALIDATION FUNCTIONS
# =============================================================================

# set library and sources

library(targets)
tar_load(c(hemp_data, best_method))
source("./R/core_functions.R")

# used
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

# used
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

# tar_load
mods <- run_complete_protein_focused_analysis(hemp_data, best_method, n_iterations = 100)

mods$summary_table

protein_analysis <- analyze_protein_focused_results(run_protein_focused_validation(hemp_data, best_method, n_iterations = 100))
create_variable_importance_table(protein_analysis = protein_analysis)

tar_load(full_spectral_analysis)


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
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", size = 0.8) +
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
tar_load(spectral_analysis)

create_combined_vip_distribution_plot(spectral_analysis, protein_analysis, alpha_level = 0.03)
