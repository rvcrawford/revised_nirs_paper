# Preprocessing functions based on your original implementation

my_preprocess <- function(train, test) {
  # Ensure inputs are matrices
  if (!is.matrix(train)) train <- as.matrix(train)
  if (!is.matrix(test)) test <- as.matrix(test)
  
  # Training set preprocessing with error handling
  raw_train <- train
  
  first_deriv_train <- tryCatch({
    t(diff(t(train), diff = 1))
  }, error = function(e) { train[, -ncol(train)] })  # Fallback if diff fails
  
  second_deriv_train <- tryCatch({
    t(diff(t(train), diff = 2))
  }, error = function(e) { train[, -c((ncol(train)-1):ncol(train))] })
  
  sg_train <- tryCatch({
    prospectr::savitzkyGolay(train, m = 1, p = 3, w = 5)
  }, error = function(e) { train })
  
  gap_der_train <- tryCatch({
    prospectr::gapDer(X = train, m = 1, w = 11, s = 5)
  }, error = function(e) { train })
  
  snv_train <- tryCatch({
    prospectr::standardNormalVariate(train)
  }, error = function(e) { train })
  
  snv_sg_train <- tryCatch({
    prospectr::standardNormalVariate(sg_train)
  }, error = function(e) { train })
  
  snv_detrend_train <- tryCatch({
    prospectr::detrend(train, wav = as.numeric(colnames(train) |> readr::parse_number()))
  }, error = function(e) { train })
  
  msc_train <- tryCatch({
    prospectr::msc(train)
  }, error = function(e) { train })
  
  # Apply same preprocessing to test set
  raw_test <- test
  first_deriv_test <- tryCatch({
    t(diff(t(test), diff = 1))
  }, error = function(e) { test[, -ncol(test)] })
  
  second_deriv_test <- tryCatch({
    t(diff(t(test), diff = 2))
  }, error = function(e) { test[, -c((ncol(test)-1):ncol(test))] })
  
  sg_test <- tryCatch({
    prospectr::savitzkyGolay(test, m = 1, p = 3, w = 5)
  }, error = function(e) { test })
  
  gap_der_test <- tryCatch({
    prospectr::gapDer(X = test, m = 1, w = 11, s = 5)
  }, error = function(e) { test })
  
  snv_test <- tryCatch({
    prospectr::standardNormalVariate(test)
  }, error = function(e) { test })
  
  snv_sg_test <- tryCatch({
    prospectr::standardNormalVariate(sg_test)
  }, error = function(e) { test })
  
  snv_detrend_test <- tryCatch({
    prospectr::detrend(test, wav = as.numeric(colnames(test) |> readr::parse_number()))
  }, error = function(e) { test })
  
  msc_test <- tryCatch({
    prospectr::msc(test, ref_spectrum = attr(msc_train, "Reference spectrum"))
  }, error = function(e) { test })
  
  # Package outputs
  output_train <- list(raw_train, first_deriv_train, second_deriv_train,
                       sg_train, gap_der_train, snv_train, snv_sg_train, 
                       snv_detrend_train, msc_train)
  names(output_train) <- c("raw_train", "first_derivative_train", "second_derivative_train",
                           "sav_gol_train", "gap_der_train", "snv_train", "snv_sg_train",
                           "snv_detrend_train", "msc_train")
  
  output_test <- list(raw_test, first_deriv_test, second_deriv_test,
                      sg_test, gap_der_test, snv_test, snv_sg_test, 
                      snv_detrend_test, msc_test)
  names(output_test) <- c("raw_test", "first_derivative_test", "second_derivative_test",
                          "sav_gol_test", "gap_der_test", "snv_test", "snv_sg_test",
                          "snv_detrend_test", "msc_test")
  
  return(list(output_train, output_test))
}

split_spectra <- function(y) {
  inTrain <- createDataPartition(
    y = y,
    p = .75,
    list = FALSE,
    groups = 3
  )
  return(inTrain)
}

get_preprocessing_methods <- function() {
  c("raw", "first_derivative", "sav_gol", "gap_der", 
    "snv", "snv_sg", "snv_detrend", "msc")
}

validate_preprocessing_method <- function(method) {
  valid_methods <- get_preprocessing_methods()
  if (!method %in% valid_methods) {
    stop("Invalid preprocessing method. Choose from: ", paste(valid_methods, collapse = ", "))
  }
  method
}