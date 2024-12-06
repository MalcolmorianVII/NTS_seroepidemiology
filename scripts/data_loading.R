# data_loading.R

#' Load Excel Data with Error Handling
#'
#' @param file_path Character string of the full file path to the Excel file
#' @param sheet Optional sheet name or index to load
#' @return Loaded data frame
#' @importFrom readxl read_excel
#' @importFrom cli cli_alert_danger
#' @importFrom dplyr filter
load_excel_data <- function(file_path,required_columns,filter_cols=NULL) {
  # Validate input file path
  if (!file.exists(file_path)) {
    cli_alert_danger(paste("Error: File not found at", file_path))
    stop("Invalid file path")
  }
  
  # Attempt to read Excel file
  tryCatch({
    # Load data with optional sheet specification
    data <- read_excel(file_path)
    
    # Filtering samples with no serogroups
    if (!is.null(filter_cols) && length(filter_cols) > 0) {
      data <- data %>% 
        filter(!is.na(filter_cols))
    }
    
    # Optional: Add data validation checks ..required cols is a list
    validate_data(data, required_columns)
    
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
validate_data <- function(data, required_columns) {
  # Check for required columns
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
  required_seq = c("sample", "h1", "h2", "o_antigen",
                        "serogroup", "serovar", "serovar_cgmlst")
  required_elisa = c("sample", "day0STmELISA", "day0SEnELISa", "day90STmELISA",
                        "day90SEnELISA", "Antibody_fold_change", "Absolute_change") # add age as well
  required_age_elisa = c("sample", "day0STmELISA", "day0SEnELISa", "day90STmELISA",
                        "day90SEnELISA", "Antibody_fold_change", "ageyr", "agemon")
  required_wilson = c("serovar")

  datasets <- list(
    saints_seq_only = load_excel_data(file_paths$main_seq,required_seq,"serogroup"),
    elisa = load_excel_data(file_paths$elisa,required_elisa),
    age_elisa = load_excel_data(file_paths$age_elisa,required_age_elisa),
    wilson_data = load_excel_data(file_paths$wilson_data,required_wilson)
  )
  
  return(datasets)
}
