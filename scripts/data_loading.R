# data_loading.R

#' Load Excel Data with Error Handling
#'
#' @param file_path Character string of the full file path to the Excel file
#' @param sheet Optional sheet name or index to load
#' @return Loaded data frame
#' @importFrom readxl read_excel
#' @importFrom cli cli_alert_danger
#' @importFrom dplyr filter
load_excel_data <- function(file_path, sheet = NULL, filter_na = TRUE) {
  # Validate input file path
  if (!file.exists(file_path)) {
    cli_alert_danger(paste("Error: File not found at", file_path))
    stop("Invalid file path")
  }
  
  # Attempt to read Excel file
  tryCatch({
    # Load data with optional sheet specification
    data <- read_excel(file_path, sheet = sheet)
    
    # Optional filtering of NA values
    if (filter_na) {
      data <- data %>% 
        filter(!is.na(serogroup))
    }
    
    # Optional: Add data validation checks
    validate_data_structure(data)
    
    return(data)
  }, error = function(e) {
    cli_alert_danger(paste("Error loading Excel file:", e$message))
    stop(e)
  })
}

#' Validate Data Structure
#'
#' @param data Data frame to validate
#' @return Invisible TRUE if valid, throws error if not
validate_data_structure <- function(data) {
  # Check for required columns
  required_columns <- c("sample", "serogroup", "serovar", "subspecies")
  missing_columns <- setdiff(required_columns, names(data))
  
  if (length(missing_columns) > 0) {
    stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
  }
  
  # Optional: Additional data validation
  if (nrow(data) == 0) {
    warning("Loaded dataset is empty")
  }
  
  return(invisible(TRUE))
}

#' Load Multiple Related Datasets
#'
#' @return List of loaded datasets
load_saints_datasets <- function() {
  # Centralized file paths (consider moving to a config file)
  file_paths <- list(
    main_seq = "/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard/2024.11.17/work/2024.11.21/2024.11.22.SAINTS_master.xlsx",
    elisa = "/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard/2024.11.17/work/2024.11.21/2024.11.21.CQJvsCQQ_elisa_seq.xlsx",
    age_elisa = "/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard/2024.11.26/2024.11.26.SAINTS_merged_with_age.xlsx",
    wilson_data = "/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard/2024.11.26/cat_wilson_salmonella_cariage_serovar.xlsx"
  )
  
  # Load datasets with error handling
  datasets <- list(
    saints_seq_only = load_excel_data(file_paths$main_seq),
    elisa = load_excel_data(file_paths$elisa),
    age_elisa = load_excel_data(file_paths$age_elisa),
    wilson_data = load_excel_data(file_paths$wilson_data)
  )
  
  return(datasets)
}

# Example usage in main script
main <- function() {
  # Load all datasets
  datasets <- load_saints_datasets()
  
  # Now you can access datasets like:
  saints_seq_only <- datasets$saints_seq_only
  elisa <- datasets$elisa
  
  # Rest of your analysis...
}