
# Function to identify and analyze the worst predicted samples

identify_worst_predicted_samples <- function(weighted_sampling_results, full_data = NULL, n_worst = 15, metric = "mae") {
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(gridExtra)
  
  cat("=== IDENTIFYING WORST PREDICTED SAMPLES ===\n")
  cat("Metric:", toupper(metric), "\n")
  cat("Number of worst samples to identify:", n_worst, "\n")
  
  # Extract predictions data
  if (is.null(weighted_sampling_results$predictions) || nrow(weighted_sampling_results$predictions) == 0) {
    stop("No predictions data found in weighted_sampling_results")
  }
  
  predictions_data <- as.data.table(weighted_sampling_results$predictions)
  
  # Remove rows with missing values
  predictions_data <- predictions_data[!is.na(actual) & !is.na(pred_unweighted) & !is.na(pred_weighted)]
  
  # Calculate prediction errors
  predictions_data[, error_unweighted := pred_unweighted - actual]
  predictions_data[, error_weighted := pred_weighted - actual]
  predictions_data[, abs_error_unweighted := abs(error_unweighted)]
  predictions_data[, abs_error_weighted := abs(error_weighted)]
  predictions_data[, sq_error_unweighted := error_unweighted^2]
  predictions_data[, sq_error_weighted := error_weighted^2]
  
  # Map locations (using the corrected mapping from previous function)
  predictions_data[, clean_location := fcase(
    # Ithaca variants
    location %in% c("Ithaca", "Eith", "Mg4e", "Mg4l", "Mcg"), "Ithaca",
    # Geneva variants  
    location %in% c("Geneva", "Cornell", "Rnooa", "Cn North"), "Geneva",
    # Freeville variants
    location %in% c("Freeville", "Free"), "Freeville",
    # Small locations -> Other
    location %in% c("Chazy", "Willsboro"), "Other",
    default = "Other"
  )]
  
  # Calculate sample-level statistics (averaging across all iterations for each sample)
  if (metric == "mae") {
    sample_stats <- predictions_data[, .(
      mean_actual = mean(actual, na.rm = TRUE),
      mean_pred_unweighted = mean(pred_unweighted, na.rm = TRUE),
      mean_pred_weighted = mean(pred_weighted, na.rm = TRUE),
      mae_unweighted = mean(abs_error_unweighted, na.rm = TRUE),
      mae_weighted = mean(abs_error_weighted, na.rm = TRUE),
      rmse_unweighted = sqrt(mean(sq_error_unweighted, na.rm = TRUE)),
      rmse_weighted = sqrt(mean(sq_error_weighted, na.rm = TRUE)),
      bias_unweighted = mean(error_unweighted, na.rm = TRUE),
      bias_weighted = mean(error_weighted, na.rm = TRUE),
      n_predictions = .N,
      clean_location = clean_location[1]
    ), by = sample_id]
    
    # Sort by MAE (unweighted as primary, then weighted)
    sample_stats[, primary_error := mae_unweighted]
    sample_stats[, secondary_error := mae_weighted]
    
  } else if (metric == "rmse") {
    sample_stats <- predictions_data[, .(
      mean_actual = mean(actual, na.rm = TRUE),
      mean_pred_unweighted = mean(pred_unweighted, na.rm = TRUE),
      mean_pred_weighted = mean(pred_weighted, na.rm = TRUE),
      mae_unweighted = mean(abs_error_unweighted, na.rm = TRUE),
      mae_weighted = mean(abs_error_weighted, na.rm = TRUE),
      rmse_unweighted = sqrt(mean(sq_error_unweighted, na.rm = TRUE)),
      rmse_weighted = sqrt(mean(sq_error_weighted, na.rm = TRUE)),
      bias_unweighted = mean(error_unweighted, na.rm = TRUE),
      bias_weighted = mean(error_weighted, na.rm = TRUE),
      n_predictions = .N,
      clean_location = clean_location[1]
    ), by = sample_id]
    
    # Sort by RMSE
    sample_stats[, primary_error := rmse_unweighted]
    sample_stats[, secondary_error := rmse_weighted]
  }
  
  # Add improvement metrics
  sample_stats[, mae_improvement := mae_unweighted - mae_weighted]
  sample_stats[, rmse_improvement := rmse_unweighted - rmse_weighted]
  sample_stats[, weighted_helped := mae_improvement > 0]
  
  # Merge with full_data to get additional metadata if provided
  if (!is.null(full_data)) {
    cat("Merging with full_data for harv_year and loc metadata...\n")
    
    # Convert full_data to data.table if needed
    if (!is.data.table(full_data)) {
      full_data <- as.data.table(full_data)
    }
    
    # Focus on key metadata columns, especially harv_year and loc
    key_metadata_cols <- c("ith_in_data_set", "harv_year", "loc")
    additional_cols <- c("year", "cultivar", "type", "in_ny", "clean_loc")
    
    metadata_cols <- intersect(names(full_data), c(key_metadata_cols, additional_cols))
    
    if (length(metadata_cols) > 0) {
      metadata <- full_data[, ..metadata_cols]
      
      # Merge with sample_stats using sample_id = ith_in_data_set
      if ("ith_in_data_set" %in% names(metadata)) {
        sample_stats <- merge(sample_stats, metadata, 
                              by.x = "sample_id", by.y = "ith_in_data_set", 
                              all.x = TRUE)
        
        # Check if we got the key fields
        if ("harv_year" %in% names(sample_stats)) {
          cat("✓ Successfully added harv_year data\n")
        } else {
          cat("✗ Warning: harv_year not found in full_data\n")
        }
        
        if ("loc" %in% names(sample_stats)) {
          cat("✓ Successfully added loc data\n")
        } else {
          cat("✗ Warning: loc not found in full_data\n")
        }
        
        cat("Added columns:", setdiff(names(sample_stats), c("sample_id", "mean_actual", "mean_pred_unweighted", "mean_pred_weighted", "mae_unweighted", "mae_weighted", "rmse_unweighted", "rmse_weighted", "bias_unweighted", "bias_weighted", "n_predictions", "clean_location", "primary_error", "secondary_error", "mae_improvement", "rmse_improvement", "weighted_helped")), "\n")
      } else {
        cat("✗ Error: Could not find 'ith_in_data_set' column in full_data for merging\n")
        cat("Available columns in full_data:", names(full_data)[1:10], "...\n")
      }
    } else {
      cat("✗ Error: No expected metadata columns found in full_data\n")
      cat("Looking for:", paste(c(key_metadata_cols, additional_cols), collapse = ", "), "\n")
      cat("Found in full_data:", names(full_data)[1:10], "...\n")
    }
  } else {
    cat("No full_data provided - using only weighted sampling results\n")
  }
  
  # Sort to find worst samples
  sample_stats_sorted <- sample_stats[order(-primary_error)]
  
  # Get worst and best samples
  worst_samples <- sample_stats_sorted[1:n_worst]
  best_samples <- sample_stats_sorted[(nrow(sample_stats_sorted) - n_worst + 1):nrow(sample_stats_sorted)]
  
  cat("\n=== WORST", n_worst, "SAMPLES ===\n")
  
  # Create display table with harv_year and loc prioritized
  base_cols <- c("sample_id", "clean_location", "mean_actual", "mae_unweighted", 
                 "mae_weighted", "mae_improvement", "weighted_helped")
  
  # Prioritize harv_year and loc, then other metadata
  priority_metadata <- c("harv_year", "loc")
  other_metadata <- c("cultivar", "type", "year")
  
  metadata_display_cols <- c(
    intersect(priority_metadata, names(worst_samples)),
    intersect(other_metadata, names(worst_samples))
  )
  
  display_cols <- c(base_cols, metadata_display_cols)
  available_cols <- intersect(display_cols, names(worst_samples))
  
  worst_table <- worst_samples[, ..available_cols]
  
  # Round numeric columns
  numeric_cols <- c("mean_actual", "mae_unweighted", "mae_weighted", "mae_improvement")
  for (col in intersect(numeric_cols, names(worst_table))) {
    worst_table[, (col) := round(get(col), 2)]
  }
  
  print(worst_table)
  
  cat("\n=== BEST", n_worst, "SAMPLES ===\n")
  best_table <- best_samples[, ..available_cols]
  
  # Round numeric columns  
  for (col in intersect(numeric_cols, names(best_table))) {
    best_table[, (col) := round(get(col), 2)]
  }
  
  print(best_table)
  
  # === FOCUSED ANALYSIS: HARV_YEAR AND LOC PATTERNS ===
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("HARVEST YEAR AND LOCATION ANALYSIS OF WORST SAMPLES\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")
  
  # Location analysis of worst samples (enhanced with original location if available)
  cat("\n=== LOCATION BREAKDOWN (CLEAN LOCATIONS) ===\n")
  worst_by_location <- worst_samples[, .N, by = clean_location][order(-N)]
  print(worst_by_location)
  
  # Original location codes breakdown - KEY ANALYSIS
  if ("loc" %in% names(worst_samples)) {
    cat("\n=== ORIGINAL LOCATION CODES (LOC) BREAKDOWN ===\n")
    worst_by_orig_location <- worst_samples[, .N, by = loc][order(-N)]
    print(worst_by_orig_location)
    
    # Compare to overall dataset proportions
    if ("loc" %in% names(sample_stats)) {
      cat("\nOverall dataset original location distribution:\n")
      overall_by_orig_location <- sample_stats[, .N, by = loc][order(-N)]
      print(overall_by_orig_location)
      
      # Calculate representation ratios
      loc_comparison <- merge(worst_by_orig_location, overall_by_orig_location, 
                              by = "loc", suffixes = c("_worst", "_total"))
      loc_comparison[, pct_of_worst := round(100 * N_worst / n_worst, 1)]
      loc_comparison[, pct_of_total := round(100 * N_total / nrow(sample_stats), 1)]
      loc_comparison[, representation_ratio := round(pct_of_worst / pct_of_total, 2)]
      
      cat("\nLocation representation analysis (values > 1.0 = overrepresented in worst samples):\n")
      print(loc_comparison[order(-representation_ratio)])
    }
  } else {
    cat("\n⚠️  Original location codes (loc) not available - check full_data merge\n")
  }
  
  # Harvest year analysis - KEY ANALYSIS
  if ("harv_year" %in% names(worst_samples)) {
    cat("\n=== HARVEST YEAR (HARV_YEAR) BREAKDOWN ===\n")
    worst_by_year <- worst_samples[, .(
      n_worst_samples = .N,
      mean_actual_cp = round(mean(mean_actual), 1),
      mean_mae_unweighted = round(mean(mae_unweighted), 2),
      mean_mae_weighted = round(mean(mae_weighted), 2),
      pct_helped_by_weighting = round(100 * mean(weighted_helped), 1)
    ), by = harv_year][order(harv_year)]
    print(worst_by_year)
    
    # Compare to overall dataset
    if ("harv_year" %in% names(sample_stats)) {
      cat("\nOverall dataset harvest year distribution:\n")
      overall_by_year <- sample_stats[, .N, by = harv_year][order(harv_year)]
      print(overall_by_year)
      
      # Calculate representation ratios
      year_comparison <- merge(worst_by_year[, .(harv_year, n_worst_samples)], 
                               overall_by_year, 
                               by = "harv_year", suffixes = c("_worst", "_total"))
      year_comparison[, pct_of_worst := round(100 * n_worst_samples / n_worst, 1)]
      year_comparison[, pct_of_total := round(100 * N / nrow(sample_stats), 1)]
      year_comparison[, representation_ratio := round(pct_of_worst / pct_of_total, 2)]
      
      cat("\nHarvest year representation analysis (values > 1.0 = overrepresented in worst samples):\n")
      print(year_comparison[order(-representation_ratio)])
    }
  } else {
    cat("\n⚠️  Harvest year (harv_year) not available - check full_data merge\n")
  }
  
  # Cross-tabulation of loc and harv_year for worst samples
  if ("loc" %in% names(worst_samples) && "harv_year" %in% names(worst_samples)) {
    cat("\n=== CROSS-TABULATION: LOC × HARV_YEAR FOR WORST SAMPLES ===\n")
    cross_tab <- worst_samples[, .N, by = .(loc, harv_year)][order(loc, harv_year)]
    print(cross_tab)
    
    # Create a wider table
    cross_tab_wide <- dcast(cross_tab, loc ~ harv_year, value.var = "N", fill = 0)
    cat("\nWide format (loc × harv_year):\n")
    print(cross_tab_wide)
  }
  
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  
  # Overall dataset location distribution for comparison
  overall_by_location <- sample_stats[, .N, by = clean_location][order(-N)]
  cat("\nOverall dataset location distribution:\n")
  print(overall_by_location)
  
  # Calculate location representation in worst samples
  comparison_table <- merge(worst_by_location, overall_by_location, by = "clean_location", suffixes = c("_worst", "_total"))
  comparison_table[, pct_of_worst := round(100 * N_worst / n_worst, 1)]
  comparison_table[, pct_of_total := round(100 * N_total / nrow(sample_stats), 1)]
  comparison_table[, overrepresented := pct_of_worst > pct_of_total]
  
  cat("\nLocation representation analysis:\n")
  print(comparison_table)
  
  # =============================================================================
  # PLOT 1: Worst Samples Comparison (Weighted vs Unweighted)
  # =============================================================================
  
  # Prepare worst samples for plotting
  worst_plot_data <- rbind(
    worst_samples[, .(sample_id, clean_location, mean_actual, 
                      error = mae_unweighted, approach = "Unweighted")],
    worst_samples[, .(sample_id, clean_location, mean_actual,
                      error = mae_weighted, approach = "Weighted")]
  )
  
  p1 <- ggplot(worst_plot_data, aes(x = reorder(sample_id, mean_actual), y = error, fill = approach)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("Unweighted" = "#A23B72", "Weighted" = "#2E86AB")) +
    labs(
      title = paste("Mean Absolute Error for", n_worst, "Worst Predicted Samples"),
      subtitle = "Samples ordered by actual CP concentration",
      x = "Sample ID (ordered by actual CP)",
      y = "Mean Absolute Error",
      fill = "Model Type"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "bottom"
    ) +
    facet_wrap(~clean_location, scales = "free_x", nrow = 1)
  
  # =============================================================================
  # PLOT 2: Worst Samples by Location
  # =============================================================================
  
  p2 <- ggplot(comparison_table, aes(x = reorder(clean_location, -N_total))) +
    geom_col(aes(y = pct_of_total), fill = "lightgray", alpha = 0.7, width = 0.6) +
    geom_col(aes(y = pct_of_worst), fill = "#A23B72", alpha = 0.8, width = 0.4) +
    geom_text(aes(y = pct_of_worst, label = paste0(N_worst, "/", N_total)), 
              vjust = -0.5, fontface = "bold") +
    labs(
      title = paste("Location Distribution: Worst", n_worst, "Samples vs Overall Dataset"),
      subtitle = "Dark bars = % of worst samples, Light bars = % of total dataset",
      x = "Location",
      y = "Percentage"
    ) +
    theme_classic()
  
  # =============================================================================
  # PLOT 3: Actual vs Predicted for Worst Samples
  # =============================================================================
  
  # Create actual vs predicted plot
  worst_scatter_data <- rbind(
    worst_samples[, .(sample_id, clean_location, actual = mean_actual, 
                      predicted = mean_pred_unweighted, approach = "Unweighted")],
    worst_samples[, .(sample_id, clean_location, actual = mean_actual,
                      predicted = mean_pred_weighted, approach = "Weighted")]
  )
  
  p3 <- ggplot(worst_scatter_data, aes(x = actual, y = predicted, color = clean_location, shape = approach)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    geom_point(size = 3, alpha = 0.8) +
    labs(
      title = paste("Actual vs Predicted CP for", n_worst, "Worst Samples"),
      subtitle = "Dashed line represents perfect prediction",
      x = "Actual CP Concentration",
      y = "Predicted CP Concentration",
      color = "Location",
      shape = "Model Type"
    ) +
    theme_classic() +
    theme(legend.position = "bottom")
  
  # =============================================================================
  # PLOT 4: Distribution of Errors
  # =============================================================================
  
  all_errors <- rbind(
    sample_stats[, .(sample_id, clean_location, error = mae_unweighted, 
                     sample_type = "All Samples", approach = "Unweighted")],
    worst_samples[, .(sample_id, clean_location, error = mae_unweighted,
                      sample_type = "Worst Samples", approach = "Unweighted")]
  )
  
  p4 <- ggplot(all_errors, aes(x = error, fill = sample_type)) +
    geom_histogram(alpha = 0.7, bins = 30, position = "identity") +
    geom_vline(data = worst_samples, aes(xintercept = min(mae_unweighted)), 
               linetype = "dashed", color = "red", size = 1) +
    scale_fill_manual(values = c("All Samples" = "lightblue", "Worst Samples" = "#A23B72")) +
    labs(
      title = "Distribution of Mean Absolute Errors",
      subtitle = paste("Red line shows threshold for worst", n_worst, "samples"),
      x = "Mean Absolute Error",
      y = "Count",
      fill = "Sample Type"
    ) +
    theme_classic() +
    theme(legend.position = "bottom")
  
  cat("\nAnalysis complete!\n")
  
  return(list(
    worst_samples = worst_samples,
    best_samples = best_samples,
    location_comparison = comparison_table,
    sample_stats = sample_stats,
    plot_worst_comparison = p1,
    plot_location_distribution = p2,
    plot_actual_vs_predicted = p3,
    plot_error_distribution = p4,
    combined_plots = grid.arrange(p1, p2, p3, p4, ncol = 2)
  ))
}

# Helper function to create a detailed table of worst samples with focus on harv_year and loc
create_worst_samples_table <- function(worst_analysis, include_metadata = TRUE) {
  library(kableExtra)
  
  worst_data <- worst_analysis$worst_samples
  
  # Basic columns
  base_table <- worst_data[, .(
    `Sample ID` = sample_id,
    Location = clean_location,
    `Actual CP` = round(mean_actual, 1),
    `Unweighted MAE` = round(mae_unweighted, 2),
    `Weighted MAE` = round(mae_weighted, 2),
    `Improvement` = round(mae_improvement, 2),
    `Weighted Better?` = ifelse(weighted_helped, "Yes", "No")
  )]
  
  # Add metadata columns if available and requested, prioritizing harv_year and loc
  if (include_metadata) {
    # Priority metadata columns
    if ("harv_year" %in% names(worst_data)) {
      base_table[["Harvest Year"]] <- worst_data$harv_year
    }
    
    if ("loc" %in% names(worst_data)) {
      base_table[["Original Location"]] <- worst_data$loc
    }
    
    # Additional metadata columns
    if ("cultivar" %in% names(worst_data)) {
      base_table[["Cultivar"]] <- worst_data$cultivar
    }
    
    if ("type" %in% names(worst_data)) {
      base_table[["Type"]] <- worst_data$type
    }
  }
  
  # Create header groupings
  n_base_cols <- 7  # Number of base columns
  n_metadata_cols <- ncol(base_table) - n_base_cols
  
  table_result <- kable(base_table,
                        caption = paste("Worst", nrow(worst_data), "Predicted Samples: Focus on Harvest Year and Location"),
                        align = rep("c", ncol(base_table))) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 10)
  
  if (n_metadata_cols > 0) {
    table_result <- table_result %>%
      add_header_above(c("Sample Info" = 2, "Model Performance" = 5, "Sample Metadata" = n_metadata_cols))
  } else {
    table_result <- table_result %>%
      add_header_above(c("Sample Info" = 2, "Model Performance" = 5))
  }
  
  table_result %>%
    footnote(general = c(
      "Improvement = Unweighted MAE - Weighted MAE. Positive values indicate weighted sampling performed better.",
      "Focus on Harvest Year and Original Location patterns for temporal and geographic analysis."
    ), general_title = "Notes:")
}

# Usage:
# full_data <- tar_read(full_data)  # or however you access your full dataset
# worst_analysis <- identify_worst_predicted_samples(weighted_sampling_results, full_data, n_worst = 15)
# worst_analysis$plot_worst_comparison  # Show specific plots
# worst_table <- create_worst_samples_table(worst_analysis, include_metadata = TRUE)
# Helper function to create a detailed table of worst samples with focus on harv_year and loc
create_worst_samples_table <- function(worst_analysis, include_metadata = TRUE) {
  library(kableExtra)
  
  worst_data <- worst_analysis$worst_samples
  
  # Basic columns
  base_table <- worst_data[, .(
    `Sample ID` = sample_id,
    Location = clean_location,
    `Actual CP` = round(mean_actual, 1),
    `Unweighted MAE` = round(mae_unweighted, 2),
    `Weighted MAE` = round(mae_weighted, 2),
    `Improvement` = round(mae_improvement, 2),
    `Weighted Better?` = ifelse(weighted_helped, "Yes", "No")
  )]
  
  # Add metadata columns if available and requested, prioritizing harv_year and loc
  if (include_metadata) {
    # Priority metadata columns
    if ("harv_year" %in% names(worst_data)) {
      base_table[["Harvest Year"]] <- worst_data$harv_year
    }
    
    if ("loc" %in% names(worst_data)) {
      base_table[["Original Location"]] <- worst_data$loc
    }
    
    # Additional metadata columns
    if ("cultivar" %in% names(worst_data)) {
      base_table[["Cultivar"]] <- worst_data$cultivar
    }
    
    if ("type" %in% names(worst_data)) {
      base_table[["Type"]] <- worst_data$type
    }
  }
  
  # Create header groupings
  n_base_cols <- 7  # Number of base columns
  n_metadata_cols <- ncol(base_table) - n_base_cols
  
  table_result <- kable(base_table,
                        caption = paste("Worst", nrow(worst_data), "Predicted Samples: Focus on Harvest Year and Location"),
                        align = rep("c", ncol(base_table))) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), font_size = 10)
  
  if (n_metadata_cols > 0) {
    table_result <- table_result %>%
      add_header_above(c("Sample Info" = 2, "Model Performance" = 5, "Sample Metadata" = n_metadata_cols))
  } else {
    table_result <- table_result %>%
      add_header_above(c("Sample Info" = 2, "Model Performance" = 5))
  }
  
  table_result %>%
    footnote(general = c(
      "Improvement = Unweighted MAE - Weighted MAE. Positive values indicate weighted sampling performed better.",
      "Focus on Harvest Year and Original Location patterns for temporal and geographic analysis."
    ), general_title = "Notes:")
}

