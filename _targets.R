# _targets.R - COMPLETE Simplified Hemp NIR Analysis Pipeline
# All analyses included: preprocessing, algorithms, spectral analysis, protein-focused

library(targets)
library(tarchetypes)

# Source essential functions only
source("R/config.R")
source("R/core_functions.R")

tar_option_set(
  packages = c(
    "tidyverse", "data.table", "caret", "prospectr", 
    "pls", "kableExtra", "ggplot2"
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
  # PREPROCESSING COMPARISON
  # =============================================================================
  
  tar_target(
    preprocessing_comparison,
    run_preprocessing_comparison(hemp_data, n_iterations = 10)
  ),
  
  tar_target(
    preprocessing_analysis,
    analyze_preprocessing_methods(preprocessing_comparison, preproc_key)
  ),
  
  tar_target(
    best_method,
    select_best_preprocessing_method(preprocessing_analysis)
  ),
  
  # =============================================================================
  # FINAL MODELING WITH BEST METHOD
  # =============================================================================
  
  tar_target(
    final_model_results,
    run_final_modeling(hemp_data, best_method, n_iterations = 10)
  ),
  
  tar_target(
    final_model_analysis,
    analyze_final_model(final_model_results)
  ),
  
  # =============================================================================
  # MULTI-ALGORITHM COMPARISON
  # =============================================================================
  
  tar_target(
    multi_algorithm_results,
    run_multi_algorithm_comparison(
      hemp_data, 
      best_method, 
      n_iterations = 6,
      algorithms = c("pls", "svmRadial", "rf")
    )
  ),
  
  tar_target(
    multi_algorithm_analysis,
    analyze_multi_algorithm_results(multi_algorithm_results)
  ),
  
  # =============================================================================
  # SPECTRAL ANALYSIS (FULL SPECTRUM)
  # =============================================================================
  
  tar_target(
    spectral_analysis,
    run_spectral_analysis(hemp_data, best_method)
  ),
  
  # =============================================================================
  # PROTEIN-FOCUSED SPECTRAL ANALYSIS
  # =============================================================================
  
  tar_target(
    protein_focused_analysis,
    run_protein_focused_analysis(hemp_data, best_method)
  ),
  
  tar_target(
    model_comparison,
    compare_full_vs_protein_models(spectral_analysis, protein_focused_analysis)
  ),
  
  # =============================================================================
  # TABLES FOR MANUSCRIPT
  # =============================================================================
  
  tar_target(
    table_sample_summary,
    create_sample_summary_table(hemp_data)
  ),
  
  tar_target(
    table_preprocessing,
    create_preprocessing_comparison_table(preprocessing_analysis)
  ),
  
  tar_target(
    table_algorithm_comparison,
    create_algorithm_comparison_table(multi_algorithm_analysis)
  ),
  
  tar_target(
    table_model_comparison,
    create_model_comparison_table(model_comparison)
  ),
  tar_target(
    error_analysis,
    analyze_prediction_errors(final_model_results, hemp_data)
  ),
  
  # =============================================================================
  # FIGURES FOR MANUSCRIPT
  # =============================================================================
  
  tar_target(
    fig_calibration,
    create_calibration_plot(final_model_results)
  ),
  
  tar_target(
    fig_performance,
    create_performance_boxplot(final_model_analysis)
  ),
  
  tar_target(
    fig_algorithm_comparison,
    create_algorithm_comparison_plot(multi_algorithm_analysis)
  ),
  
  tar_target(
    fig_vip_scores,
    create_vip_plot(spectral_analysis)
  ),
  
  tar_target(
    fig_model_comparison,
    create_model_comparison_plot(spectral_analysis, protein_focused_analysis)
  ),
  
  tar_target(
    fig_validation_errors,
    create_validation_error_plot(error_analysis)
  )
  
  # =============================================================================
  # MANUSCRIPT RENDERING (OPTIONAL)
  # =============================================================================
  
  # Uncomment to auto-render manuscript when pipeline completes:
  # tar_target(
  #   manuscript,
  #   quarto::quarto_render("manuscript/hemp_nir_paper.qmd"),
  #   format = "file"
  # )
)