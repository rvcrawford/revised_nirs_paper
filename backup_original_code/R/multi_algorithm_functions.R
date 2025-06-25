# Multi-algorithm modeling functions for NIRS hemp protein prediction
# Add this as: R/multi_algorithm_functions.R

library(caret)
library(data.table)
library(tidyverse)

# 1. Main function to compare multiple algorithms
run_multi_algorithm_comparison <- function(data, best_method, n_iterations = 100, algorithms = c("pls", "svmRadial", "rf")) {
  
  cat("=== MULTI-ALGORITHM COMPARISON ===\n")
  cat("Algorithms:", paste(algorithms, collapse = ", "), "\n")
  cat("Preprocessing method:", best_method, "\n")
  cat("Iterations:", n_iterations, "\n\n")
  
  # Extract spectral data
  spectral_cols <- grep("^x[0-9]+$", names(data), value = TRUE)
  spectra_matrix <- as.matrix(data[, ..spectral_cols])
  y <- data$crude_protein
  
  # Extract original sample identifiers  
  if ("ith_in_data_set" %in% names(data)) {
    original_sample_ids <- data$ith_in_data_set
  } else {
    original_sample_ids <- 1:nrow(data)
  }
  
  cat("Data dimensions:", dim(spectra_matrix), "\n")
  cat("Response range:", range(y), "\n\n")
  
  # Storage for results
  results_list <- list()
  predictions_list <- list()
  successful_iterations <- 0
  
  for (i in 1:n_iterations) {
    if (i %% 25 == 0) cat("Iteration", i, "- successful so far:", successful_iterations, "\n")
    
    tryCatch({
      # Use SAME train/test split for all algorithms (critical for fair comparison)
      set.seed(i)  # Ensure reproducible splits
      inTrain <- split_spectra(y)
      y_train <- y[inTrain]
      y_test <- y[-inTrain]
      test_sample_ids <- original_sample_ids[-inTrain]
      
      # Apply preprocessing (same for all algorithms)
      spectra_processed <- my_preprocess(spectra_matrix[inTrain, ], spectra_matrix[-inTrain, ])
      
      # Get preprocessed data for the best method
      method_train_name <- paste0(best_method, "_train")
      method_test_name <- paste0(best_method, "_test")
      
      if (!method_train_name %in% names(spectra_processed[[1]])) {
        stop("Method ", method_train_name, " not found")
      }
      
      train_data <- spectra_processed[[1]][[method_train_name]]
      test_data <- spectra_processed[[2]][[method_test_name]]
      
      # Prepare training data frame
      train_df <- data.frame(y_train = y_train, train_data)
      test_df <- data.frame(test_data)
      
      # Fit each algorithm on the SAME data
      iteration_results <- list()
      iteration_predictions <- list()
      
      for (alg in algorithms) {
        
        # Algorithm-specific parameters
        if (alg == "pls") {
          model <- train(
            y_train ~ .,
            data = train_df,
            method = "pls",
            tuneLength = 20,
            trControl = trainControl(method = "cv", number = 5, verboseIter = FALSE)
          )
        } else if (alg == "svmRadial") {
          model <- train(
            y_train ~ .,
            data = train_df,
            method = "svmRadial",
            tuneLength = 10,
            trControl = trainControl(method = "cv", number = 5, verboseIter = FALSE),
            preProcess = c("center", "scale")  # SVM benefits from scaling
          )
        } else if (alg == "rf") {
          model <- train(
            y_train ~ .,
            data = train_df,
            method = "rf",
            tuneLength = 10,
            trControl = trainControl(method = "cv", number = 5, verboseIter = FALSE),
            importance = TRUE  # For feature importance
          )
        } else {
          # Generic approach for other algorithms
          model <- train(
            y_train ~ .,
            data = train_df,
            method = alg,
            tuneLength = 10,
            trControl = trainControl(method = "cv", number = 5, verboseIter = FALSE)
          )
        }
        
        # Make predictions
        predictions <- predict(model, newdata = test_df)
        
        # Calculate metrics
        rmse_val <- sqrt(mean((y_test - predictions)^2))
        rsq_val <- cor(y_test, predictions)^2
        rpd_val <- sd(y_test) / rmse_val
        rpiq_val <- (quantile(y_test, 0.75) - quantile(y_test, 0.25)) / rmse_val
        
        # Store results
        iteration_results[[alg]] <- data.table(
          iteration = i,
          algorithm = alg,
          rmse = rmse_val,
          rsq = rsq_val,
          rpd = rpd_val,
          rpiq = rpiq_val,
          n_predictors = ncol(train_data),
          best_tune = if("bestTune" %in% names(model)) paste(names(model$bestTune), model$bestTune, sep="=", collapse="; ") else "default"
        )
        
        # Store predictions
        iteration_predictions[[alg]] <- data.table(
          iteration = i,
          algorithm = alg,
          sample_id = 1:length(y_test),
          ith_in_data_set = test_sample_ids,
          actual = y_test,
          predicted = predictions
        )
      }
      
      # Combine results for this iteration
      results_list[[length(results_list) + 1]] <- rbindlist(iteration_results)
      predictions_list[[length(predictions_list) + 1]] <- rbindlist(iteration_predictions)
      
      successful_iterations <- successful_iterations + 1
      
    }, error = function(e) {
      cat("ERROR in iteration", i, ":", e$message, "\n")
    })
  }
  
  cat("\nCompleted", successful_iterations, "successful iterations\n")
  
  if (length(results_list) == 0) {
    stop("No successful iterations")
  }
  
  # Combine all results
  final_results <- list(
    metrics = rbindlist(results_list),
    predictions = rbindlist(predictions_list),
    summary_stats = NULL  # Will be calculated by analysis function
  )
  
  cat("Final metrics rows:", nrow(final_results$metrics), "\n")
  cat("Final predictions rows:", nrow(final_results$predictions), "\n")
  
  return(final_results)
}

# 2. Analysis function for multi-algorithm results
analyze_multi_algorithm_results <- function(multi_algo_results) {
  
  cat("=== ANALYZING MULTI-ALGORITHM RESULTS ===\n")
  
  metrics_data <- multi_algo_results$metrics
  predictions_data <- multi_algo_results$predictions
  
  # Calculate summary statistics by algorithm
  summary_stats <- metrics_data[, .(
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
  
  # Calculate model quality classifications by algorithm
  model_quality <- metrics_data[, .(
    excellent_pct = mean(rpd > 3 & rpiq > 4.1 & rsq > 0.8, na.rm = TRUE) * 100,
    good_pct = mean(rpd >= 2.5 & rpd <= 3 & rpiq >= 2.3 & rpiq <= 4.1 & rsq > 0.8, na.rm = TRUE) * 100,
    fair_pct = mean(rpd >= 2.0 & rpd < 2.5 & rpiq >= 2.3, na.rm = TRUE) * 100,
    poor_but_functional_pct = mean(rpd >= 1.5 & rpd < 2.0, na.rm = TRUE) * 100,
    inadequate_pct = mean(rpd < 1.5, na.rm = TRUE) * 100
  ), by = algorithm]
  
  model_quality[, total_acceptable := excellent_pct + good_pct + fair_pct + poor_but_functional_pct]
  model_quality[, quantitative_capable := excellent_pct + good_pct + fair_pct]
  
  # Statistical comparisons (pairwise t-tests)
  algorithms <- unique(metrics_data$algorithm)
  pairwise_comparisons <- data.table()
  
  if (length(algorithms) > 1) {
    for (metric in c("rmse", "rsq", "rpd", "rpiq")) {
      for (i in 1:(length(algorithms)-1)) {
        for (j in (i+1):length(algorithms)) {
          alg1 <- algorithms[i]
          alg2 <- algorithms[j]
          
          values1 <- metrics_data[algorithm == alg1][[metric]]
          values2 <- metrics_data[algorithm == alg2][[metric]]
          
          # Paired t-test (same iterations)
          test_result <- t.test(values1, values2, paired = TRUE)
          
          pairwise_comparisons <- rbind(pairwise_comparisons, data.table(
            metric = metric,
            algorithm1 = alg1,
            algorithm2 = alg2,
            mean_diff = mean(values1 - values2),
            p_value = test_result$p.value,
            significant = test_result$p.value < 0.05
          ))
        }
      }
    }
  }
  
  # Print summary
  cat("\nSUMMARY STATISTICS BY ALGORITHM:\n")
  print(summary_stats)
  
  cat("\nMODEL QUALITY CLASSIFICATIONS:\n")
  print(model_quality)
  
  if (nrow(pairwise_comparisons) > 0) {
    cat("\nSIGNIFICANT DIFFERENCES (p < 0.05):\n")
    significant_diffs <- pairwise_comparisons[significant == TRUE]
    if (nrow(significant_diffs) > 0) {
      print(significant_diffs[, .(metric, algorithm1, algorithm2, mean_diff, p_value)])
    } else {
      cat("No significant differences found\n")
    }
  }
  
  return(list(
    summary_stats = summary_stats,
    model_quality = model_quality,
    pairwise_comparisons = pairwise_comparisons,
    raw_metrics = metrics_data,
    raw_predictions = predictions_data
  ))
}

# 3. Create comparison plots
create_algorithm_comparison_plot <- function(multi_algo_analysis) {
  
  metrics_data <- multi_algo_analysis$raw_metrics
  
  # Create comparison boxplot
  metrics_long <- metrics_data %>%
    pivot_longer(
      cols = c(rmse, rsq, rpd, rpiq),
      names_to = "metric",
      values_to = "value"
    ) %>%
    mutate(
      metric_label = case_when(
        metric == "rmse" ~ "RMSE",
        metric == "rsq" ~ "R²",
        metric == "rpd" ~ "RPD", 
        metric == "rpiq" ~ "RPIQ"
      ),
      algorithm_label = case_when(
        algorithm == "pls" ~ "PLS",
        algorithm == "svmRadial" ~ "SVM",
        algorithm == "rf" ~ "Random Forest",
        TRUE ~ toupper(algorithm)
      )
    )
  
  p <- ggplot(metrics_long, aes(x = algorithm_label, y = value, fill = algorithm_label)) +
    geom_boxplot(alpha = 0.7) +
    facet_wrap(~ factor(metric_label, levels = c("RMSE", "R²", "RPD", "RPIQ")), 
               scales = "free_y", 
               nrow = 2) +
    labs(
      title = "Algorithm Comparison: PLS vs SVM vs Random Forest",
      subtitle = "Performance on hemp grain protein prediction (100 iterations)",
      x = "Algorithm",
      y = "Metric Value",
      fill = "Algorithm"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12)
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")
  
  return(p)
}

# 4. Create summary table for manuscript
# REPLACE the create_algorithm_comparison_table function in your R/multi_algorithm_functions.R with this:

create_algorithm_comparison_table <- function(multi_algo_analysis) {
  
  summary_stats <- multi_algo_analysis$summary_stats
  model_quality <- multi_algo_analysis$model_quality
  
  # Combine summary stats with quality metrics
  combined_table <- merge(summary_stats, model_quality, by = "algorithm")
  
  # Format for manuscript - FIXED VERSION
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
  
  # Create final table with proper column names
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
  
  # Create kable table
  table_output <- final_table %>%
    knitr::kable(
      caption = "Performance comparison of machine learning algorithms for hemp grain protein prediction",
      row.names = FALSE,
      align = c("l", rep("c", 6))
    )
  
  return(table_output)
}
# 5. Generate interpretation text
generate_algorithm_interpretation <- function(multi_algo_analysis) {
  
  summary_stats <- multi_algo_analysis$summary_stats
  pairwise_comparisons <- multi_algo_analysis$pairwise_comparisons
  
  # Find best performing algorithm by RMSE
  best_algorithm <- summary_stats[which.min(mean_rmse)]$algorithm
  best_rmse <- round(summary_stats[which.min(mean_rmse)]$mean_rmse, 3)
  
  # Find best performing algorithm by RPD
  best_rpd_algorithm <- summary_stats[which.max(mean_rpd)]$algorithm
  best_rpd <- round(summary_stats[which.max(mean_rpd)]$mean_rpd, 2)
  
  # Check for significant differences
  significant_diffs <- pairwise_comparisons[significant == TRUE]
  
  # Algorithm name mapping
  alg_names <- list(
    "pls" = "PLS",
    "svmRadial" = "SVM", 
    "rf" = "Random Forest"
  )
  
  interpretation <- paste0(
    "Comparison of three machine learning algorithms revealed that ",
    alg_names[[best_algorithm]], " achieved the lowest RMSE (", best_rmse, "), ",
    "while ", alg_names[[best_rpd_algorithm]], " achieved the highest RPD (", best_rpd, "). "
  )
  
  if (nrow(significant_diffs) > 0) {
    interpretation <- paste0(interpretation,
                             "Statistical analysis revealed significant differences between algorithms for ",
                             length(unique(significant_diffs$metric)), " out of 4 performance metrics. "
    )
  } else {
    interpretation <- paste0(interpretation,
                             "Statistical analysis revealed no significant differences between algorithms, ",
                             "suggesting that the choice of algorithm may be less important than proper preprocessing ",
                             "and feature selection for this application. "
    )
  }
  
  interpretation <- paste0(interpretation,
                           "These results demonstrate that multiple machine learning approaches can successfully ",
                           "predict hemp grain protein content from NIR spectra, providing researchers with ",
                           "flexibility in model selection based on specific requirements for interpretability, ",
                           "computational efficiency, or prediction accuracy."
  )
  
  return(interpretation)
}