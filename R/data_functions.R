# Corrected data preparation functions based on actual data structure

prepare_full_data <- function(full_data_raw) {
  # Convert to data.table for easier manipulation
  full_data <- copy(full_data_raw)
  
  # Convert crude protein from % to g/kg (multiply by 10)
  full_data[, crude_protein := crude_protein * 10]
  
  # Create cleaner location names for analysis
  full_data[, clean_loc := fcase(
    loc %in% c("ithaca", "eith"), "Ithaca",
    loc %in% c("freev", "free"), "Freeville", 
    loc == "chazy", "Chazy",
    loc == "willsboro", "Willsboro",
    loc %in% c("rsgeneva", "rn041_gen", "rn028_gen"), "Geneva",
    loc == "cnoll", "Cornell",
    default = tools::toTitleCase(tolower(loc))
  )]
  
  cat("Data preparation complete:\n")
  cat("Samples:", nrow(full_data), "\n")
  cat("CP range:", round(min(full_data$crude_protein)), "-", round(max(full_data$crude_protein)), "g/kg\n")
  cat("Unique locations:", length(unique(full_data$clean_loc)), "\n")
  cat("Unique cultivars:", length(unique(full_data$cultivar)), "\n")
  
  return(full_data)
}

prepare_supplemental_table <- function(supp_tab_raw) {
  supp_tab2 <- copy(supp_tab_raw)
  supp_tab2[, Cultivar := toupper(cultivar2)]
  
  # Use base R function instead of stringr
  names(supp_tab2) <- tools::toTitleCase(tolower(names(supp_tab2)))
  
  # Return sortable version
  supp_tab2[, c(8, 2:7)]
}

prepare_preprocessing_key <- function(preproc_key_raw) {
  preproc_key <- copy(preproc_key_raw)
  
  # Match the 8 methods in your analysis
  preproc_key[, full_name := c(
    "Raw Spectra",                                    # raw
    "First Derivative",                               # first_derivative  
    "Savitzky-Golay",                                # sav_gol
    "Gap-segment Derivative",                         # gap_der
    "Standard Normal Variate",                        # snv
    "Standard Normal Variate following Savitzky-Golay", # snv_sg
    "Standard Normal Variate-Detrend",               # snv_detrend
    "Multiplicative Scatter Correction"              # msc
  )]
  
  preproc_key
}

calculate_protein_summary <- function(full_data) {
  # Calculate summary statistics for protein data
  cp_data <- full_data$crude_protein
  
  # Create summary - values should already be in g/kg from prepare_full_data
  data.frame(
    Mean = round(mean(cp_data, na.rm = TRUE), 0),
    SD = round(sd(cp_data, na.rm = TRUE), 0),
    Minimum = round(min(cp_data, na.rm = TRUE), 0),
    `First Quartile` = round(quantile(cp_data, 0.25, na.rm = TRUE), 0),
    Median = round(median(cp_data, na.rm = TRUE), 0),
    `Third Quartile` = round(quantile(cp_data, 0.75, na.rm = TRUE), 0),
    Maximum = round(max(cp_data, na.rm = TRUE), 0),
    check.names = FALSE
  )
}

calculate_location_summary <- function(full_data) {
  # Summary by location
  location_summary <- full_data[, .N, by = .(clean_loc, in_ny)]
  setorder(location_summary, -N)
  
  cat("Location summary:\n")
  print(location_summary)
  
  # Return just NY vs non-NY count for backward compatibility
  table(full_data$in_ny)
}

# Utility function for skewness calculation
calculate_skewness <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  sd_x <- sqrt(sum((x - mean_x)^2) / n)
  z <- (x - mean_x) / sd_x
  sum(z^3) / n
}