library(targets)
library(tarchetypes)

# Source all functions
source("R/data_functions.R")
source("R/preprocessing_functions.R") 
source("R/modeling_functions.R")
source("R/analysis_functions.R")
source("R/plotting_functions.R")
source("R/spectral_analysis_functions.R")  # NEW
source("R/multi_algorithm_functions.R")  # ADD THIS to your source() statements at top
source("R/location_weighting_functions.R")
source("R/fixed_location_weighting_function.R")

# Set target options
# 
# In your _targets.R file, update the packages list:
tar_option_set(
  packages = c(
    "tidyverse", "data.table", "caret", "prospectr", 
    "pls", "tidymodels", "nlme", "kableExtra", "skimr",
    "emmeans", "multcomp", "broom", "lme4", "readr", "gridExtra"  # ADDED gridExtra
  ),
  format = "rds"
)

list(
  # Data loading and preparation
  tar_target(
    full_data_raw,
    fread("./input_data/final_data_set/full_hemp_data.csv")
  ),
  
  tar_target(
    supp_tab_raw,
    fread("./input_data/final_data_set/cultivar_table_clean.csv")
  ),
  
  tar_target(
    preproc_key_raw,
    fread("./input_data/final_data_set/preprocessing_key.csv")
  ),
  
  # Clean and prepare data
  tar_target(
    full_data,
    prepare_full_data(full_data_raw)
  ),
  
  tar_target(
    supp_tab,
    prepare_supplemental_table(supp_tab_raw)
  ),
  
  tar_target(
    preproc_key,
    prepare_preprocessing_key(preproc_key_raw)
  ),
  
  # Summary statistics
  tar_target(
    protein_summary,
    calculate_protein_summary(full_data)
  ),
  
  tar_target(
    location_summary, 
    calculate_location_summary(full_data)
  ),
  
  # Preprocessing method comparison (100 iterations)
  tar_target(
    preprocessing_results,
    run_preprocessing_comparison(full_data, n_iterations = 100)
  ),
  
  tar_target(
    preprocessing_analysis,
    analyze_preprocessing_methods(preprocessing_results, preproc_key)
  ),
  
  tar_target(
    best_preprocessing_method,
    select_best_preprocessing(preprocessing_analysis)
  ),
  
  # Final model development (1000 iterations with best method)
  tar_target(
    final_model_results,
    run_single_method_modeling(
      full_data, 
      best_preprocessing_method, 
      n_iterations = 1000
    )
  ),
  tar_target(
    balanced_data,
    prepare_location_balanced_data(full_data)
  ),
  
  tar_target(
    weighting_comparison,
    compare_weighted_unweighted(
      balanced_data, 
      full_data, 
      best_preprocessing_method, 
      n_iterations = 1000
    )
  ),
  tar_target(
    # Fixed comparison that preserves location information
    weighting_comparison_fixed,
    compare_weighted_unweighted_fixed(
      balanced_data, 
      full_data, 
      best_preprocessing_method, 
      n_iterations = 1000
    )
  ),
  
  # Simple Geneva analysis that avoids the column selection issue
  # Replace your geneva_performance_analysis target with this:
  
  tar_target(
    geneva_performance_analysis,
    {
      cat("=== SIMPLE GENEVA ANALYSIS (NO COLUMN ISSUES) ===\n")
      
      library(dplyr)
      
      wc_fixed <- weighting_comparison_fixed
      bd <- balanced_data
      
      # Simple approach - calculate directly without complex pivoting
      
      # Weighted results by location
      weighted_results <- wc_fixed$weighted$predictions %>%
        group_by(location) %>%
        summarise(
          weighted_rmse = sqrt(mean((actual - predicted)^2)),
          weighted_bias = mean(predicted - actual),
          weighted_mae = mean(abs(actual - predicted)),
          weighted_rsq = cor(actual, predicted)^2,
          weighted_n = n(),
          .groups = 'drop'
        )
      
      # Unweighted results by location  
      unweighted_results <- wc_fixed$unweighted$predictions %>%
        group_by(location) %>%
        summarise(
          unweighted_rmse = sqrt(mean((actual - predicted)^2)),
          unweighted_bias = mean(predicted - actual),
          unweighted_mae = mean(abs(actual - predicted)),
          unweighted_rsq = cor(actual, predicted)^2,
          unweighted_n = n(),
          .groups = 'drop'
        )
      
      # Simple join and calculate improvements
      location_comparison <- weighted_results %>%
        left_join(unweighted_results, by = "location") %>%
        mutate(
          rmse_improvement = unweighted_rmse - weighted_rmse,
          rmse_pct_improvement = 100 * (unweighted_rmse - weighted_rmse) / unweighted_rmse,
          bias_improvement = abs(unweighted_bias) - abs(weighted_bias),
          rsq_improvement = weighted_rsq - unweighted_rsq
        ) %>%
        left_join(bd$location_summary, by = "location")
      
      cat("Location comparison results:\n")
      print(location_comparison %>% 
              select(location, weight, weighted_rmse, unweighted_rmse, 
                     rmse_improvement, rmse_pct_improvement))
      
      # Geneva-specific analysis
      geneva_results <- location_comparison %>% filter(location == "geneva")
      
      cat("\n=== GENEVA RESULTS ===\n")
      cat("Weight:", round(geneva_results$weight, 2), "x\n")
      cat("Weighted RMSE:", round(geneva_results$weighted_rmse, 3), "\n")
      cat("Unweighted RMSE:", round(geneva_results$unweighted_rmse, 3), "\n")
      cat("RMSE improvement:", round(geneva_results$rmse_improvement, 3), "\n")
      cat("% improvement:", round(geneva_results$rmse_pct_improvement, 1), "%\n")
      cat("R² improvement:", round(geneva_results$rsq_improvement, 4), "\n")
      
      # Statistical test for Geneva
      geneva_weighted_errors <- wc_fixed$weighted$predictions %>% 
        filter(location == "geneva") %>%
        mutate(abs_error = abs(actual - predicted)) %>%
        pull(abs_error)
      
      geneva_unweighted_errors <- wc_fixed$unweighted$predictions %>% 
        filter(location == "geneva") %>%
        mutate(abs_error = abs(actual - predicted)) %>%
        pull(abs_error)
      
      cat("\nGeneva sample sizes:\n")
      cat("Weighted:", length(geneva_weighted_errors), "\n")
      cat("Unweighted:", length(geneva_unweighted_errors), "\n")
      
      geneva_significance <- NULL
      
      if (length(geneva_weighted_errors) > 10 && length(geneva_unweighted_errors) > 10) {
        tryCatch({
          wilcox_test <- wilcox.test(geneva_weighted_errors, 
                                     geneva_unweighted_errors,
                                     alternative = "less")
          
          # Effect size
          mean_diff <- mean(geneva_unweighted_errors) - mean(geneva_weighted_errors)
          pooled_sd <- sqrt((var(geneva_weighted_errors) + var(geneva_unweighted_errors)) / 2)
          cohens_d <- mean_diff / pooled_sd
          
          cat("\nStatistical test results:\n")
          cat("p-value:", round(wilcox_test$p.value, 4), "\n")
          cat("Cohen's d:", round(cohens_d, 3), "\n")
          cat("Significant (p < 0.05):", wilcox_test$p.value < 0.05, "\n")
          cat("Meaningful effect (|d| > 0.2):", abs(cohens_d) > 0.2, "\n")
          
          geneva_significance <- list(
            p_value = wilcox_test$p.value,
            cohens_d = cohens_d,
            significant = wilcox_test$p.value < 0.05,
            meaningful = abs(cohens_d) > 0.2
          )
          
          # Final conclusion
          cat("\n=== GENEVA CONCLUSION ===\n")
          if (geneva_significance$significant && geneva_significance$meaningful) {
            cat("🎯 WEIGHTING SIGNIFICANTLY HELPS GENEVA!\n")
            cat("✓ Statistically significant AND meaningful effect size\n")
            conclusion <- "significant_improvement"
          } else if (geneva_results$rmse_improvement > 0.05) {
            cat("⚠️ WEIGHTING HELPS GENEVA BUT NOT SIGNIFICANTLY\n")
            cat("✓ Shows improvement but needs larger effect\n")
            conclusion <- "modest_improvement"
          } else {
            cat("❌ WEIGHTING DOES NOT MEANINGFULLY HELP GENEVA\n")
            conclusion <- "no_improvement"
          }
          
        }, error = function(e) {
          cat("Error in statistical test:", e$message, "\n")
          conclusion <- "test_failed"
        })
      } else {
        cat("Insufficient data for statistical test\n")
        conclusion <- "insufficient_data"
      }
      
      return(list(
        location_comparison = location_comparison,
        geneva_results = geneva_results,
        geneva_significance = geneva_significance,
        conclusion = if(exists("conclusion")) conclusion else "unknown"
      ))
    }
  ),
  
  tar_target(
    # Create Geneva-focused plot
    fig_geneva_improvement,
    create_geneva_improvement_plot(geneva_performance_analysis)
  ),
  
  tar_target(
    # Generate interpretation text
    geneva_interpretation,
    generate_geneva_interpretation(geneva_performance_analysis)
  ),
  
  tar_target(
    # Summary table of location-specific results
    table_location_specific_results,
    {
      library(knitr)
      library(dplyr)
      
      location_data <- geneva_performance_analysis$location_comparison %>%
        select(location, weight, rmse_unweighted, rmse_weighted, 
               rmse_improvement, rmse_pct_improvement) %>%
        mutate(
          significant = case_when(
            location == "geneva" & !is.null(geneva_performance_analysis$geneva_significance) ~ 
              ifelse(geneva_performance_analysis$geneva_significance$significant, "*", ""),
            TRUE ~ ""
          )
        )
      
      kable(location_data,
            col.names = c("Location", "Weight", "Unweighted RMSE", "Weighted RMSE", 
                          "Improvement", "% Improvement", "Significant"),
            digits = 3,
            caption = "Location-specific RMSE performance comparison (* p < 0.05)")
    }
  ),
  
  # Optional: Compare with current weighting_comparison for validation
  tar_target(
    location_comparison_validation,
    {
      cat("=== VALIDATING LOCATION PRESERVATION ===\n")
      
      # Check if location info is preserved in both approaches
      old_comparison <- weighting_comparison  # Your original
      new_comparison <- weighting_comparison_fixed  # Fixed version
      
      cat("Original weighted predictions columns:", names(old_comparison$weighted$predictions), "\n")
      cat("Fixed weighted predictions columns:", names(new_comparison$weighted$predictions), "\n")
      
      cat("Original unweighted predictions columns:", names(old_comparison$unweighted$predictions), "\n") 
      cat("Fixed unweighted predictions columns:", names(new_comparison$unweighted$predictions), "\n")
      
      # Check if we can now do location analysis
      has_location_old_weighted <- "location" %in% names(old_comparison$weighted$predictions)
      has_location_old_unweighted <- "location" %in% names(old_comparison$unweighted$predictions)
      has_location_new_weighted <- "location" %in% names(new_comparison$weighted$predictions)
      has_location_new_unweighted <- "location" %in% names(new_comparison$unweighted$predictions)
      
      cat("Location column availability:\n")
      cat("- Original weighted:", has_location_old_weighted, "\n")
      cat("- Original unweighted:", has_location_old_unweighted, "\n")
      cat("- Fixed weighted:", has_location_new_weighted, "\n")
      cat("- Fixed unweighted:", has_location_new_unweighted, "\n")
      
      validation_success <- has_location_new_weighted && has_location_new_unweighted
      cat("Geneva analysis possible:", validation_success, "\n")
      
      return(list(
        validation_success = validation_success,
        can_analyze_geneva = validation_success
      ))
    }
  ),
  tar_target(
    weighting_analysis,
    analyze_weighting_comparison(weighting_comparison)
  ),
  
  tar_target(
    weighted_model_results,
    run_weighted_modeling_ultrafast(
      balanced_data,
      best_preprocessing_method, 
      n_iterations = 1000
    )
  ),
  
  tar_target(
    location_performance_analysis,
    analyze_location_performance(weighted_model_results, balanced_data)
  ),
  tar_target(
    debug_weights,
    debug_weighting_data(weighting_comparison, weighting_analysis)
  ),
  tar_target(
    table_weighting_comparison,
    create_weighting_comparison_table(weighting_comparison)
  ),
  
  tar_target(
    fig_weighting_comparison,
    create_weighting_comparison_plot(weighting_comparison)  # Use weighting_comparison, NOT weighting_analysis
  ),
  
  tar_target(
    fig_location_performance,
    create_location_performance_plot(location_performance_analysis)
  ),
  
  tar_target(
    final_model_analysis,
    analyze_final_models(final_model_results)
  ),
  
  tar_target(
    error_analysis,
    analyze_prediction_errors(final_model_results, full_data)
  ),
  tar_target(
    multi_algorithm_results,
    run_multi_algorithm_comparison(
      full_data, 
      best_preprocessing_method, 
      n_iterations = 100,        # ← VERY SMALL for testing
      algorithms = c("pls", "svmRadial", "rf")  # ← Skip RF initially
    )
  ),
  
  tar_target(
    multi_algorithm_analysis,
    analyze_multi_algorithm_results(multi_algorithm_results)
  ),
  
  tar_target(
    fig_algorithm_comparison,
    create_algorithm_comparison_plot(multi_algorithm_analysis)
  ),
  
  tar_target(
    table_algorithm_comparison,
    create_algorithm_comparison_table(multi_algorithm_analysis)
  ),
  
  tar_target(
    algorithm_interpretation,
    generate_algorithm_interpretation(multi_algorithm_analysis)
  ),

  # NEW: Spectral feature analysis
  tar_target(
    spectral_analysis,
    run_spectral_analysis(full_data, best_preprocessing_method)
  ),
  # Protein-focused spectral analysis
  tar_target(
    protein_focused_analysis,
    run_protein_focused_analysis(full_data, best_preprocessing_method)
  ),
  
  # Comparison between full and protein-focused models
  tar_target(
    model_comparison,
    compare_full_vs_protein_models(spectral_analysis, protein_focused_analysis)
  ),
  
  # Protein-focused plots
  tar_target(
    fig_protein_coefficients,
    create_coefficient_plot(protein_focused_analysis)
  ),
  
  tar_target(
    fig_protein_vip,
    create_vip_plot(protein_focused_analysis)
  ),
  
  # Comparison plot (both models side by side)
  tar_target(
    fig_model_comparison,
    create_model_comparison_plot(spectral_analysis, protein_focused_analysis)
  ),
  
  # Protein-focused table
  tar_target(
    table_protein_wavelengths,
    create_wavelength_table(protein_focused_analysis)
  ),
  
  # Performance comparison table
  tar_target(
    table_model_comparison,
    create_performance_comparison_table(model_comparison)
  ),
  
  # Interpretation text
  tar_target(
    protein_interpretation,
    generate_protein_interpretation(protein_focused_analysis, model_comparison)
  ),
  
  tar_target(
    fig_coefficients,
    create_coefficient_plot(spectral_analysis)
  ),
  
  tar_target(
    fig_vip_scores,
    create_vip_plot(spectral_analysis)
  ),
  
  tar_target(
    table_wavelengths,
    create_wavelength_table(spectral_analysis)
  ),
  
  tar_target(
    wavelength_interpretation,
    generate_wavelength_interpretation(spectral_analysis)
  ),
  
  # Tables
  tar_target(
    table_hemp_provenance,
    create_hemp_provenance_table(supp_tab)
  ),
  
  tar_target(
    table_protein_summary,
    create_protein_summary_table(protein_summary)
  ),
  
  tar_target(
    table_preprocessing_comparison,
    create_preprocessing_table(preprocessing_analysis)
  ),
  
  # Figures
  tar_target(
    fig_model_calibration,
    create_calibration_plot(final_model_results)
  ),
  
  tar_target(
    fig_final_metrics,
    create_final_metrics_plot(final_model_analysis)
  ),
  
  tar_target(
    fig_validation_performance,
    create_validation_plot(error_analysis, full_data)
  ),
  
  # Export outputs
  tar_target(
    export_tables,
    export_all_tables(
      table_hemp_provenance,
      table_protein_summary, 
      table_preprocessing_comparison
    ),
    format = "file"
  ),
  
  tar_target(
    export_figures,
    export_all_figures(
      fig_model_calibration,
      fig_final_metrics,
      fig_validation_performance
    ),
    format = "file"
  )
)