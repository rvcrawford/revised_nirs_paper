



#' Create table for weighted sampling comparison results
create_weighted_sampling_table <- function(analysis_results) {
  library(kableExtra)
  library(dplyr)
  
  # Extract the performance table from the analysis
  performance_table <- analysis_results$performance_table
  
  # Format for publication
  table_formatted <- performance_table %>%
    mutate(
      Metric = c("RMSE", "R²", "RPD", "RPIQ"),
      `Unweighted Mean` = sprintf("%.3f", Unweighted_Mean),
      `Weighted Mean` = sprintf("%.3f", Weighted_Mean), 
      `Mean Improvement` = sprintf("%.3f", Mean_Improvement),
      `% Iterations Better` = sprintf("%.1f%%", Pct_Better)
    ) %>%
    select(Metric, `Unweighted Mean`, `Weighted Mean`, `Mean Improvement`, `% Iterations Better`)
  
  # Create formatted table
  kable(table_formatted,
        caption = "Performance comparison: Weighted sampling vs Unweighted models (1000 iterations)",
        align = c("l", "c", "c", "c", "c")) %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
    add_header_above(c(" " = 1, "Mean Performance" = 2, "Improvement" = 2)) %>%
    footnote(general = "Positive improvements indicate weighted sampling performed better. RMSE improvement = Unweighted - Weighted (lower RMSE is better).",
             general_title = "Note:")
}

#' Create visualization of weighted sampling comparison
create_weighted_sampling_plot <- function(analysis_results) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # Extract metrics from full results
  metrics_data <- analysis_results$full_metrics
  
  # Prepare data for plotting - focus on improvement distributions
  improvement_data <- metrics_data %>%
    select(iteration, rmse_improvement, rsq_improvement, rpd_improvement) %>%
    pivot_longer(cols = ends_with("_improvement"), 
                 names_to = "metric", 
                 values_to = "improvement") %>%
    mutate(
      metric = case_when(
        metric == "rmse_improvement" ~ "RMSE Improvement",
        metric == "rsq_improvement" ~ "R² Improvement", 
        metric == "rpd_improvement" ~ "RPD Improvement"
      ),
      metric = factor(metric, levels = c("RMSE Improvement", "R² Improvement", "RPD Improvement"))
    )
  
  # Create the plot
  p <- ggplot(improvement_data, aes(x = improvement, fill = metric)) +
    geom_histogram(bins = 30, alpha = 0.7, color = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
    facet_wrap(~metric, scales = "free", ncol = 1) +
    labs(
      title = "Distribution of Performance Improvements\nWeighted Sampling vs Unweighted Models",
      subtitle = "1000 iterations with different train/test splits",
      x = "Improvement (Weighted - Unweighted)",
      y = "Frequency",
      caption = "Positive values indicate weighted sampling performed better.\nRed line at zero indicates no improvement."
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11)
    ) +
    scale_fill_manual(values = c("#2E86AB", "#A23B72", "#F18F01"))
  
  return(p)
}

#' Analyze location-specific performance from weighted sampling
analyze_location_specific_performance <- function(comparison_results, balanced_data) {
  library(dplyr)
  
  # Extract predictions with location info
  predictions <- comparison_results$predictions
  
  # Calculate location-specific metrics
  location_performance <- predictions %>%
    filter(!is.na(location)) %>%
    group_by(iteration, location) %>%
    summarise(
      n_samples = n(),
      unw_rmse = sqrt(mean((actual - pred_unweighted)^2, na.rm = TRUE)),
      w_rmse = sqrt(mean((actual - pred_weighted)^2, na.rm = TRUE)),
      rmse_improvement = unw_rmse - w_rmse,
      unw_rsq = cor(actual, pred_unweighted, use = "complete.obs")^2,
      w_rsq = cor(actual, pred_weighted, use = "complete.obs")^2,
      rsq_improvement = w_rsq - unw_rsq,
      .groups = 'drop'
    ) %>%
    group_by(location) %>%
    summarise(
      total_samples = mean(n_samples),  # Should be consistent across iterations
      mean_rmse_improvement = mean(rmse_improvement, na.rm = TRUE),
      sd_rmse_improvement = sd(rmse_improvement, na.rm = TRUE),
      pct_rmse_better = 100 * mean(rmse_improvement > 0, na.rm = TRUE),
      mean_rsq_improvement = mean(rsq_improvement, na.rm = TRUE),
      sd_rsq_improvement = sd(rsq_improvement, na.rm = TRUE),
      pct_rsq_better = 100 * mean(rsq_improvement > 0, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Add weight information from balanced_data
  if("location_summary" %in% names(balanced_data)) {
    location_performance <- location_performance %>%
      left_join(
        balanced_data$location_summary %>% 
          select(location, weight) %>%
          rename(location_weight = weight),
        by = "location"
      )
  }
  
  return(location_performance)
}

#' Create robustness plot showing confidence intervals
create_robustness_plot <- function(analysis_results) {
  library(ggplot2)
  library(dplyr)
  
  # Extract summary stats
  summary_stats <- analysis_results$summary_stats
  
  # Prepare data for plotting
  plot_data <- data.frame(
    Metric = c("RMSE", "R²", "RPD"),
    Mean_Improvement = c(
      summary_stats$mean_rmse_improvement,
      summary_stats$mean_rsq_improvement, 
      summary_stats$mean_rpd_improvement
    ),
    CI_Lower = c(
      summary_stats$rmse_improvement_ci_lower,
      summary_stats$rsq_improvement_ci_lower,
      quantile(analysis_results$full_metrics$rpd_improvement, 0.025, na.rm = TRUE)
    ),
    CI_Upper = c(
      summary_stats$rmse_improvement_ci_upper,
      summary_stats$rsq_improvement_ci_upper,
      quantile(analysis_results$full_metrics$rpd_improvement, 0.975, na.rm = TRUE)
    ),
    Pct_Better = c(
      summary_stats$pct_rmse_better,
      summary_stats$pct_rsq_better,
      summary_stats$pct_rpd_better
    )
  )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = Metric, y = Mean_Improvement)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), 
                  width = 0.2, size = 1, color = "#2E86AB") +
    geom_point(size = 4, color = "#2E86AB") +
    geom_text(aes(label = paste0(round(Pct_Better, 1), "%\nbetter")), 
              vjust = -0.5, hjust = 0.5, size = 3.5) +
    labs(
      title = "Weighted Sampling Improvement with 95% Confidence Intervals",
      subtitle = "Based on 1000 different train/test splits",
      x = "Performance Metric",
      y = "Mean Improvement (Weighted - Unweighted)",
      caption = "Error bars show 95% confidence intervals.\nPercentages show fraction of iterations where weighted performed better."
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11)
    )
  
  return(p)
}
