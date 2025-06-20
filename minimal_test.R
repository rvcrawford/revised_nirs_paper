library(tidyverse)
library(data.table) 
library(caret)
library(prospectr)
library(pls)

# Source your existing functions
source("R/data_functions.R")
source("R/preprocessing_functions.R") 
source("R/modeling_functions.R")
source("R/analysis_functions.R")

# Load data
full_data_raw <- fread("./input_data/final_data_set/full_hemp_data.csv")
full_data <- prepare_full_data(full_data_raw)

# Test just the preprocessing
spectral_cols <- grep("^x[0-9]+$", names(full_data), value = TRUE)
spectra_matrix <- as.matrix(full_data[, ..spectral_cols])
y <- full_data$crude_protein

inTrain <- split_spectra(y)
processed_spectra <- my_preprocess(spectra_matrix[inTrain, ], spectra_matrix[-inTrain, ])

cat("Preprocessing methods:", length(processed_spectra[[1]]), "\n")
cat("Method names:", names(processed_spectra[[1]]), "\n")