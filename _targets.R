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
  
  # Generate manuscript
  tar_quarto(
    manuscript,
    "manuscript/hemp_nir_paper.qmd"
  )
)