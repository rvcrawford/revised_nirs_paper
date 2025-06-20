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
    table_weighting_comparison,
    create_weighting_comparison_table(weighting_analysis)
  ),
  
  tar_target(
    fig_weighting_comparison,
    create_weighting_comparison_plot(weighting_analysis)
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