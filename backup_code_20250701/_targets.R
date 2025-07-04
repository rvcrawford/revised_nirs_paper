# STREAMLINED _targets.R - Hemp NIR Paper Essentials Only
# Focus: NIR prediction of hemp protein with preprocessing comparison

library(targets)
library(tarchetypes)

# Source only essential functions  
source("R/config.R")           # Development/production modes
source("R/core_functions.R")   # All essential functions consolidated

tar_option_set(
  packages = c(
    "tidyverse", "data.table", "caret", "prospectr", 
    "pls", "kableExtra", "ggplot2", "quarto"
  ),
  format = "rds"
)

list(
  # =============================================================================
  # DATA PIPELINE
  # =============================================================================
  
  tar_target(
    hemp_data_raw,
    fread("./input_data/final_data_set/full_hemp_data.csv")
  ),
  tar_target(
    preproc_key_raw,
    fread("./input_data/final_data_set/preprocessing_key.csv")
  ),
  
  tar_target(
    hemp_data,
    prepare_hemp_data(hemp_data_raw)
  ),
  
  tar_target(
    preproc_key,
    prepare_preprocessing_key(preproc_key_raw)
  ),
  
  # =============================================================================
  # CORE ANALYSIS
  # =============================================================================
  
  # Compare 8 preprocessing methods
  tar_target(
    preprocessing_comparison,
    run_preprocessing_comparison(hemp_data, n_iterations = 3)
  ),
  
  # Analyze results and select best method
  tar_target(
    preprocessing_analysis,
    analyze_preprocessing_methods(preprocessing_comparison, preproc_key)
  ),
  
  tar_target(
    best_method,
    select_best_preprocessing_method(preprocessing_analysis)
  ),
  
  # ADD THIS TARGET (missing in your current _targets.R):
  tar_target(
    spectral_analysis,
    run_spectral_analysis(hemp_data, best_method)
  ),
  
  tar_target(
    multi_algorithm_comparison,
    run_multi_algorithm_comparison(
      data = hemp_data,
      best_method = preproc_key[method_name%in%best_method]$preproc,
      n_iterations = 100,                            # Perfect for stark differences
      algorithms = c("pls", "svmRadial", "rf"),     # All 3 algorithms
      training_mode = FALSE,                        # Full validation rigor
      fair_comparison = TRUE                        # Fair hyperparameter optimization
    )
  ),
  
  tar_target(
    multi_algorithm_analysis,
    analyze_multi_algorithm_results(multi_algorithm_comparison)
  ),
  
  # Create visualizations and outputs
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
  
  # Your existing targets:
  tar_target(
    protein_focused_analysis, 
    run_protein_focused_analysis(hemp_data, best_method)
  ),
  
  tar_target(
    model_comparison, 
    compare_full_vs_protein_models(spectral_analysis, protein_focused_analysis)
  ),
  
  
  # Final modeling with best method
  tar_target(
    final_model_results,
    run_final_modeling(hemp_data, best_method, n_iterations = 3)
  ),
  
  tar_target(
    final_model_analysis,
    analyze_final_model_performance(final_model_results)
  ),
  
  tar_target(
    error_analysis,
    analyze_prediction_errors(final_model_results, hemp_data)
  ),
  
  # =============================================================================
  # PAPER OUTPUTS
  # =============================================================================
  
  # Tables
  tar_target(
    table_sample_summary,
    create_sample_summary_table(hemp_data)
  ),
  
  tar_target(
    table_preprocessing,
    create_preprocessing_comparison_table(preprocessing_analysis)
  ),
  
  # Figures  
  tar_target(
    fig_model_calibration,
    create_calibration_plot(final_model_results)
  ),
  
  tar_target(
    fig_performance_boxplot,
    create_performance_boxplot(final_model_analysis)
  ),
  
  tar_target(
    fig_validation_errors,
    create_validation_error_plot(error_analysis)
  ),
  
  tar_target(fig_model_comparison, create_model_comparison_plot(spectral_analysis, protein_focused_analysis)),
  
  tar_target(table_model_comparison, create_performance_comparison_table(model_comparison))
  
  # =============================================================================
  # MANUSCRIPT GENERATION
  # =============================================================================
  
  # =============================================================================
  # MANUSCRIPT GENERATION - WITH EXPLICIT DEPENDENCIES
  # =============================================================================
  
  # =============================================================================
  # MANUSCRIPT GENERATION - SIMPLE VERSION
  # =============================================================================
  
  # tar_target(
  #   manuscript,
  #   {
  #     cat("Rendering manuscript with data:\n")
  #     cat("- Hemp data:", nrow(hemp_data), "samples\n")
  #     cat("- Best method:", best_method, "\n")
  #     
  #     # Just render in place - no moving files around
  #     quarto::quarto_render("./manuscript/hemp_nir_paper.qmd")
  #   },
  #   format = "file"
  # )
)
