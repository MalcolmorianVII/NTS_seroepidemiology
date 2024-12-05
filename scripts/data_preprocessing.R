# data_preprocessing.R

#' Preprocess Salmonella Sequence Data
#'
#' @param data Raw input dataframe
#' @return Preprocessed dataframe
preprocess_saints_data <- function(data) {
  # Exclude problematic samples and NA entries
  data <- data %>%
    filter(
      !(sample == "CQJ142" | is.na(serogroup))
    )
  return(data)
}

#' Split Data by Subspecies
#'
#' @param data Preprocessed dataframe
#' @return List containing subspecies I and II dataframes
split_by_subspecies <- function(data) {
  data$subspecies <- ifelse(grepl("II", data$serovar_cgmlst), "II", "I")
  list(
    subspecies1 = data %>% filter(subspecies == "I"),
    subspecies2 = data %>% filter(subspecies == "II")
  )
}

#' Get Singletons in Dataset
#'
#' @param data Dataframe to analyze
#' @param classification Column to check for singletons
#' @return Dataframe of singleton entries
get_singletons <- function(data, classification = "serovar") {
  singletons <- data %>%
    group_by(!!sym(classification)) %>%
    filter(n() == 1) %>%
    ungroup()
  
  singletons$singletype <- classification
  return(singletons)
}


#' Remove Singletons from Dataset
#'
#' @param data Dataframe to process
#' @param classification Column to check for multiple entries
#' @return Dataframe with singletons removed
remove_singletons <- function(data, classification = "serovar") {
  data %>%
    group_by(!!sym(classification)) %>%
    filter(n() > 1) %>%
    ungroup()
  
}

#' Reorder Categorical Columns by Frequency
#'
#' @param data Dataframe to process
#' @param columns Columns to reorder
#' @return Dataframe with reordered factor columns
reorder_columns_by_frequency <- function(data, columns = c("serovar", "serogroup", "serovar_cgmlst")) {
  # Function to reorder a single column
  reorder_column <- function(data, column) {
    if (column %in% names(data)) {
      data[[column]] <- factor(
        data[[column]], 
        levels = names(sort(table(data[[column]]), decreasing = TRUE))
      )
    }
    return(data)
  }
  
  # Apply reordering to specified columns
  for (col in columns) {
    data <- reorder_column(data, col)
  }
  return(data)
}

#' Perform Comprehensive Data Preprocessing
#'
#' @param saints_seq_only Raw input dataframe
#' @return List of processed datasets and metadata
process_salmonella_data <- function(saints_seq_only) {
  # Preprocess main dataset
  preprocessed_data <- preprocess_saints_data(saints_seq_only)
  
  # Split by subspecies
  subspecies_data <- split_by_subspecies(preprocessed_data)
  
  # Identify and handle singletons
  singletons <- list(
    subspecies1_O_antigen = get_singletons(subspecies_data$subspecies1, "serovar"),
    subspecies1_cgMLST = get_singletons(subspecies_data$subspecies1, "serovar_cgmlst"),
    subspecies2_O_antigen = get_singletons(subspecies_data$subspecies2, "serovar"),
    subspecies2_cgMLST = get_singletons(subspecies_data$subspecies2, "serovar_cgmlst")
  )
  
  # Remove singletons from datasets
  cleaned_data <- list(
    subspecies1_serovar = remove_singletons(subspecies_data$subspecies1, "serovar"),
    subspecies1_cgMLST = remove_singletons(subspecies_data$subspecies1, "serovar_cgmlst"),
    subspecies2_serovar = remove_singletons(subspecies_data$subspecies2, "serovar"),
    subspecies2_cgMLST = remove_singletons(subspecies_data$subspecies2, "serovar_cgmlst")
  )
  
  # Reorder columns by frequency
  reordered_data <- list(
    subspecies1 = reorder_columns_by_frequency(subspecies_data$subspecies1),
    subspecies2 = reorder_columns_by_frequency(subspecies_data$subspecies2)
  )
  
  # Compute metadata
  metadata <- list(
    subspecies1_total = nrow(subspecies_data$subspecies1),
    subspecies2_total = nrow(subspecies_data$subspecies2)
  )
  
  # Write singletons to Excel (optional)
  write_singletons_to_excel(singletons)
  
  # Return comprehensive processing results
  return(list(
    raw_data = preprocessed_data,
    subspecies_data = subspecies_data,
    singletons = singletons,
    cleaned_data = cleaned_data,
    reordered_data = reordered_data,
    metadata = metadata
  ))
}

#' Write Singletons to Excel File
#'
#' @param singletons List of singleton dataframes
#' @return NULL (writes file to disk)
write_singletons_to_excel <- function(singletons) {
  # Combine singletons into a single dataframe
  singletons_df <- do.call(rbind, singletons)
  
  # Write to Excel
  output_path <- file.path(
    CONFIG$output_dir, 
    paste0("serovar_singletons_", format(Sys.Date(), "%Y.%m.%d"), ".xlsx")
  )
  
  writexl::write_xlsx(singletons_df, output_path)
  
  # Optional: Log file writing
  message(paste("Singletons written to:", output_path))
}
