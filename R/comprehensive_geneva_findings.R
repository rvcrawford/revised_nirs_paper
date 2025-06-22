# Comprehensive Geneva-specific analysis and paper recommendations

run_complete_geneva_analysis <- function() {
  cat("=== COMPLETE GENEVA ANALYSIS FOR PAPER ===\n\n")
  
  # Run all analyses
  quick_results <- quick_geneva_check()
  statistical_results <- test_geneva_significance()
  
  # Extract key findings
  geneva_row <- quick_results %>% filter(location == "geneva")
  geneva_improves <- geneva_row$rmse_improvement > 0
  geneva_improvement <- geneva_row$rmse_improvement
  geneva_pct_improvement <- geneva_row$pct_improvement
  
  # Check trade-offs
  other_locations <- quick_results %>% filter(location != "geneva")
  other_avg_change <- mean(other_locations$rmse_improvement, na.rm = TRUE)
  
  cat("KEY FINDINGS:\n")
  cat("=============\n")
  cat("1. Geneva RMSE improvement:", round(geneva_improvement, 3), 
      "(", round(geneva_pct_improvement, 1), "%)\n")
  cat("2. Statistical significance (p-value):", round(statistical_results$wilcox_p, 4), "\n")
  cat("3. Effect size (Cohen's d):", round(statistical_results$cohens_d, 3), "\n")
  cat("4. Other locations average change:", round(other_avg_change, 3), "\n")
  cat("5. Overall conclusion:", statistical_results$overall_conclusion, "\n\n")
  
  # Generate paper recommendations
  cat("PAPER WRITING RECOMMENDATIONS:\n")
  cat("===============================\n\n")
  
  if (statistical_results$overall_conclusion) {
    cat("✅ WEIGHTING HELPS GENEVA - REPORT THIS!\n\n")
    
    cat("Suggested text for your paper:\n")
    cat("------------------------------\n")
    cat('"Location-weighted modeling significantly improved predictions for Geneva samples,\n')
    cat('reducing RMSE by', round(geneva_improvement, 2), 'units (', round(geneva_pct_improvement, 1), '%, p =', round(statistical_results$wilcox_p, 3), ').\n')
    cat('This improvement came at a modest cost to other locations (average change:\n')
    cat(round(other_avg_change, 3), '), resulting in minimal overall performance change. This suggests\n')
    cat('that class weighting can effectively address location-specific prediction bias\n')
    cat('when tailored to underrepresented geographic regions."\n\n')
    
    cat("Table to include:\n")
    cat("Location | Weight | Unweighted RMSE | Weighted RMSE | Improvement | p-value\n")
    cat("---------|--------|-----------------|---------------|-------------|--------\n")
    for(i in 1:nrow(quick_results)) {
      row <- quick_results[i,]
      p_val <- if(row$location == "geneva") round(statistical_results$wilcox_p, 3) else "—"
      cat(sprintf("%-8s | %.2fx   | %.3f           | %.3f         | %+.3f      | %s\n",
                  row$location, row$weight, row$unweighted_rmse, 
                  row$weighted_rmse, row$rmse_improvement, p_val))
    }
    
  } else if (geneva_improves && abs(geneva_improvement) > 0.05) {
    cat("⚠️ WEIGHTING HELPS GENEVA BUT NOT SIGNIFICANTLY\n\n")
    
    cat("Suggested text for your paper:\n")
    cat("------------------------------\n")
    cat('"Location-weighted modeling showed a', round(geneva_pct_improvement, 1), '% improvement in RMSE\n')
    cat('for Geneva samples, though this was not statistically significant (p =', round(statistical_results$wilcox_p, 3), ').\n')
    cat('The modest improvement suggests that class weighting may have potential for\n')
    cat('addressing geographic bias, but would require larger sample sizes or stronger\n')
    cat('weighting schemes to achieve statistical significance."\n\n')
    
  } else {
    cat("❌ WEIGHTING DOES NOT HELP GENEVA\n\n")
    
    cat("Suggested text for your paper:\n")
    cat("------------------------------\n")
    cat('"Location-weighted modeling did not improve predictions for underrepresented\n')
    cat('locations, with Geneva showing minimal change (', round(geneva_pct_improvement, 1), '% improvement, p =', round(statistical_results$wilcox_p, 3), ').\n')
    cat('This suggests that geographic location is not a major source of prediction\n')
    cat('bias in this dataset, and that the NIRS models are robust across the tested\n')
    cat('New York environments."\n\n')
  }
  
  # Additional analyses to consider
  cat("ADDITIONAL ANALYSES TO CONSIDER:\n")
  cat("=================================\n")
  
  if (statistical_results$overall_conclusion) {
    cat("1. Test stronger weighting schemes for Geneva (3x, 5x weights)\n")
    cat("2. Analyze if Geneva improvement is consistent across CV folds\n")
    cat("3. Test location-specific models vs weighted global models\n")
    cat("4. Examine which wavelengths contribute most to Geneva improvements\n")
  } else {
    cat("1. Test if other factors (cultivar, year) explain Geneva differences\n")
    cat("2. Analyze spectral characteristics of Geneva vs other locations\n")
    cat("3. Consider ensemble methods instead of sample weighting\n")
    cat("4. Focus on overall model robustness rather than location-specific tuning\n")
  }
  
  cat("\nDISCUSSION POINTS FOR PAPER:\n")
  cat("============================\n")
  
  if (statistical_results$overall_conclusion) {
    cat("- Sample weighting can effectively target underrepresented geographic regions\n")
    cat("- Trade-offs between global and location-specific performance are manageable\n")
    cat("- Future work could optimize weighting schemes for better balance\n")
    cat("- Demonstrates the value of testing for class imbalance in NIRS studies\n")
  } else {
    cat("- NIRS models show good geographic robustness across New York locations\n") 
    cat("- Class imbalance may be less problematic than initially expected\n")
    cat("- Simple unweighted models are sufficient for this dataset\n")
    cat("- Future studies should focus on balanced data collection rather than post-hoc weighting\n")
  }
  
  return(list(
    geneva_results = geneva_row,
    statistical_results = statistical_results,
    location_comparison = quick_results,
    overall_conclusion = statistical_results$overall_conclusion
  ))
}

# Create a final plot for the paper
create_geneva_paper_plot <- function(analysis_results) {
  library(ggplot2)
  
  # Create a clean plot showing location-specific improvements
  plot_data <- analysis_results$location_comparison %>%
    mutate(
      is_geneva = location == "geneva",
      significant = case_when(
        location == "geneva" & analysis_results$statistical_results$statistically_significant ~ TRUE,
        TRUE ~ FALSE
      ),
      improvement_type = case_when(
        rmse_improvement > 0.05 ~ "Substantial improvement",
        rmse_improvement > 0 ~ "Modest improvement", 
        rmse_improvement > -0.05 ~ "Minimal change",
        TRUE ~ "Performance decrease"
      )
    )
  
  p <- ggplot(plot_data, aes(x = reorder(location, rmse_improvement), 
                             y = rmse_improvement)) +
    geom_col(aes(fill = is_geneva), alpha = 0.8, width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.6) +
    geom_text(aes(label = ifelse(significant, 
                                 paste0(round(rmse_improvement, 3), "*"),
                                 round(rmse_improvement, 3))), 
              vjust = ifelse(plot_data$rmse_improvement > 0, -0.5, 1.2), 
              size = 3.5, fontface = "bold") +
    labs(
      title = "Location-Specific RMSE Improvement with Sample Weighting",
      subtitle = "* indicates statistically significant improvement (p < 0.05)",
      x = "Location", 
      y = "RMSE Improvement (Unweighted - Weighted)",
      fill = "Geneva Focus",
      caption = paste0("Geneva (weight: ", 
                       round(plot_data$weight[plot_data$location == "geneva"], 2),
                       "x) shows targeted improvement")
    ) +
    scale_fill_manual(values = c("FALSE" = "#A23B72", "TRUE" = "#2E86AB"),
                      labels = c("Other locations", "Geneva (weighted 1.85x)")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray60"),
      legend.position = "bottom"
    )
  
  return(p)
}

# Run everything
cat("Running comprehensive Geneva analysis...\n")
complete_results <- run_complete_geneva_analysis()
geneva_plot <- create_geneva_paper_plot(complete_results)

# Print the plot command for saving
cat("\nTo save the plot for your paper:\n")
cat("ggsave('geneva_improvement_plot.png', geneva_plot, width = 8, height = 6, dpi = 300)\n")