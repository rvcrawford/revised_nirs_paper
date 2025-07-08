# code emmeans comparison of these

preprocessing_comparison[,lapply(.SD, mean), .SDcols = c(4:7), by = preprocessing_method]
preproc_key
preprocessing_comparison[,c("preprocessing_method", "iteration"):=lapply(.SD, as.factor), .SDcols = c("preprocessing_method", "iteration")]

library(emmeans)
library(nlme)
library(multcomp)
mods <- preprocessing_comparison[,lapply(.SD, function(x) lme(x~preprocessing_method,
                                                            random = ~1|iteration)), .SDcols = c(4:7)]


method_comparison_marginal_means <- preprocessing_comparison |> 
  mutate(across(c(preprocessing_method, iteration), as.factor)) |> 
  pivot_longer(-c(1:3)) |> 
  nest(-name) |> 
  mutate(model = map(data, ~lme(value~preprocessing_method, random = ~1|iteration,
          weights = varIdent(form = ~ 1 | preprocessing_method), data = .x)),
         marginal_means = map(model, ~emmeans(.x, "preprocessing_method") |> multcomp::cld() |> 
                                as.data.frame()))


# =============================================================================
# FRESH FINAL MODEL ANALYSIS FUNCTION - RPD PRIMARY
# =============================================================================

analyze_final_model_rpd_primary <- function(final_model_results) {
  # Analyze final model performance using RPD-primary unified criteria
  # This function uses a completely different name to avoid conflicts
  
  cat("=== USING FRESH RPD-PRIMARY ANALYSIS ===\n")
  
  metrics <- final_model_results$metrics
  predictions <- final_model_results$predictions
  
  # DEBUG: Check data
  cat("Total models:", nrow(metrics), "\n")
  cat("RPD range:", paste(round(range(metrics$rpd, na.rm = TRUE), 2), collapse = " to "), "\n")
  cat("Mean RPD:", round(mean(metrics$rpd, na.rm = TRUE), 2), "\n")
  
  # Overall statistics
  overall_stats <- metrics[, .(
    mean_rmse = mean(rmse, na.rm = TRUE),
    sd_rmse = sd(rmse, na.rm = TRUE),
    mean_rsq = mean(rsq, na.rm = TRUE),
    sd_rsq = sd(rsq, na.rm = TRUE),
    mean_rpd = mean(rpd, na.rm = TRUE),
    sd_rpd = sd(rpd, na.rm = TRUE),
    mean_rpiq = mean(rpiq, na.rm = TRUE),
    sd_rpiq = sd(rpiq, na.rm = TRUE),
    n_successful = .N
  )]
  
  # Classification using RPD ONLY (primary approach)
  # Create classification vectors
  excellent_mask <- metrics$rpd > 3.0
  good_mask <- metrics$rpd >= 2.0 & metrics$rpd <= 3.0
  fair_mask <- metrics$rpd >= 1.4 & metrics$rpd < 2.0
  poor_mask <- metrics$rpd < 1.4
  
  # Count each category
  excellent_count <- sum(excellent_mask, na.rm = TRUE)
  good_count <- sum(good_mask, na.rm = TRUE)
  fair_count <- sum(fair_mask, na.rm = TRUE)
  poor_count <- sum(poor_mask, na.rm = TRUE)
  
  # DEBUG: Show classification counts
  cat("Classification counts:\n")
  cat("- Excellent (RPD > 3.0):", excellent_count, "\n")
  cat("- Good (RPD 2.0-3.0):", good_count, "\n")
  cat("- Fair (RPD 1.4-2.0):", fair_count, "\n")
  cat("- Poor (RPD < 1.4):", poor_count, "\n")
  cat("- Total:", excellent_count + good_count + fair_count + poor_count, "\n")
  
  # Supporting metric analysis (for those in each RPD category)
  excellent_with_support <- sum(excellent_mask & metrics$rpiq > 4.1 & metrics$rsq > 0.9, na.rm = TRUE)
  good_with_support <- sum(good_mask & metrics$rpiq >= 2.3 & metrics$rsq >= 0.8, na.rm = TRUE)
  fair_with_support <- sum(fair_mask & metrics$rpiq >= 1.5 & metrics$rsq >= 0.5, na.rm = TRUE)
  
  # Calculate percentages
  total_models <- nrow(metrics)
  excellent_pct <- round(100 * excellent_count / total_models, 1)
  good_pct <- round(100 * good_count / total_models, 1)
  fair_pct <- round(100 * fair_count / total_models, 1)
  poor_pct <- round(100 * poor_count / total_models, 1)
  
  # Supporting metric concordance rates
  excellent_support_pct <- ifelse(excellent_count > 0, round(100 * excellent_with_support / excellent_count, 1), 0)
  good_support_pct <- ifelse(good_count > 0, round(100 * good_with_support / good_count, 1), 0)
  fair_support_pct <- ifelse(fair_count > 0, round(100 * fair_with_support / fair_count, 1), 0)
  
  # Summary metrics
  quantitative_capable <- excellent_pct + good_pct
  total_acceptable <- excellent_pct + good_pct + fair_pct
  
  # Create results structure
  performance_summary <- list(
    total_models = total_models,
    excellent_models = excellent_count,
    good_models = good_count,
    fair_models = fair_count,
    poor_models = poor_count,
    excellent_percent = excellent_pct,
    good_percent = good_pct,
    fair_percent = fair_pct,
    poor_percent = poor_pct,
    excellent_support_rate = excellent_support_pct,
    good_support_rate = good_support_pct,
    fair_support_rate = fair_support_pct,
    quantitative_capable = quantitative_capable,
    total_acceptable = total_acceptable
  )
  
  # Print results
  cat("\n=== FINAL RESULTS ===\n")
  cat("Performance metrics (mean ± SD):\n")
  cat("- RMSE:", round(overall_stats$mean_rmse, 2), "±", round(overall_stats$sd_rmse, 2), "g/kg\n")
  cat("- R²:", round(overall_stats$mean_rsq, 3), "±", round(overall_stats$sd_rsq, 3), "\n")
  cat("- RPD:", round(overall_stats$mean_rpd, 2), "±", round(overall_stats$sd_rpd, 2), "\n")
  cat("- RPIQ:", round(overall_stats$mean_rpiq, 2), "±", round(overall_stats$sd_rpiq, 2), "\n\n")
  
  cat("Model classification (RPD-primary):\n")
  cat("- Excellent (RPD > 3.0):", excellent_count, "(", excellent_pct, "%)\n")
  cat("- Good (RPD 2.0-3.0):", good_count, "(", good_pct, "%)\n")
  cat("- Fair (RPD 1.4-2.0):", fair_count, "(", fair_pct, "%)\n")
  cat("- Poor (RPD < 1.4):", poor_count, "(", poor_pct, "%)\n\n")
  
  cat("Supporting metric concordance:\n")
  cat("- Excellent with full support:", excellent_support_pct, "%\n")
  cat("- Good with full support:", good_support_pct, "%\n")
  cat("- Fair with full support:", fair_support_pct, "%\n\n")
  
  cat("Summary:\n")
  cat("- Quantitative capable:", quantitative_capable, "%\n")
  cat("- Total acceptable:", total_acceptable, "%\n")
  
  return(list(
    overall_stats = overall_stats,
    performance_summary = performance_summary,
    metrics = metrics,
    predictions = predictions,
    classification_method = "rpd_primary_fresh"
  ))
}

analyze_final_model_rpd_primary(final_model_results)
# =============================================================================
# SIMPLE TEST FUNCTION
# =============================================================================

test_classification_logic <- function(rpd_values) {
  # Simple test to verify classification logic
  
  cat("Testing classification logic:\n")
  
  for (rpd in rpd_values) {
    category <- case_when(
      rpd > 3.0 ~ "Excellent",
      rpd >= 2.0 & rpd <= 3.0 ~ "Good",
      rpd >= 1.4 & rpd < 2.0 ~ "Fair",
      rpd < 1.4 ~ "Poor",
      TRUE ~ "Unknown"
    )
    cat("RPD =", rpd, "→", category, "\n")
  }
}

# Test with some example values
# test_classification_logic(c(1.2, 1.5, 2.0, 2.5, 3.0, 3.5))


# analyze_final_model <- function(final_model_results) {
#   # Analyze final model performance
#   
#   metrics <- final_model_results$metrics
#   predictions <- final_model_results$predictions
#   
#   # Overall statistics
#   overall_stats <- metrics[, .(
#     mean_rmse = mean(rmse, na.rm = TRUE),
#     sd_rmse = sd(rmse, na.rm = TRUE),
#     mean_rsq = mean(rsq, na.rm = TRUE),
#     sd_rsq = sd(rsq, na.rm = TRUE),
#     mean_rpd = mean(rpd, na.rm = TRUE),
#     sd_rpd = sd(rpd, na.rm = TRUE),
#     mean_rpiq = mean(rpiq, na.rm = TRUE),
#     sd_rpiq = sd(rpiq, na.rm = TRUE),
#     n_successful = .N
#   )]
#   
#   # Performance categories
#   good_models <- metrics[rpd >= 2.0 & rsq >= 0.7, .N]
#   excellent_models <- metrics[rpd >= 2.5 & rsq >= 0.8, .N]
#   total_models <- nrow(metrics)
#   
#   performance_summary <- list(
#     total_models = total_models,
#     good_models = good_models,
#     excellent_models = excellent_models,
#     good_percent = round(100 * good_models / total_models, 1),
#     excellent_percent = round(100 * excellent_models / total_models, 1)
#   )
#   
#   cat("Final model analysis:\n")
#   cat("- Total successful models:", total_models, "\n")
#   cat("- Good models (RPD≥2.0, R²≥0.7):", good_models, "(", performance_summary$good_percent, "%)\n")
#   cat("- Excellent models (RPD≥2.5, R²≥0.8):", excellent_models, "(", performance_summary$excellent_percent, "%)\n")
#   
#   return(list(
#     overall_stats = overall_stats,
#     performance_summary = performance_summary,
#     metrics = metrics,
#     predictions = predictions
#   ))
# }
