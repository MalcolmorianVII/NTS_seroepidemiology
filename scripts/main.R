# Load required libraries
source("scripts/libraries.R")

# Load configuration
# source("R/config.R")
# Load data
source("scripts/data_loading.R")
# # Source utility functions
# source("R/utils.R")
source("scripts/data_preprocessing.R")
# source("R/statistical_analysis.R")
# source("R/visualization.R")

# Main analysis pipeline
main <- function() {
  # Load all datasets
  datasets <- load_saints_datasets()
  saints_seq_only <- datasets$saints_seq_only
  elisa <- datasets$elisa
  age_elisa <- datasets$age_elisa
  wilson_data <- datasets$wilson_data
  
  # Preprocess data
  processed_data <- process_salmonella_data(saints_seq_only)

  # Access processed components
  subspecies1 <- processed_data$reordered_data$subspecies1
  subspecies2 <- processed_data$reordered_data$subspecies2
  
  # Visualize serovar,serogroup distributions in subspecies 1,2

  # Perform analyses
  perform_statistical_analyses(processed_data)
  
  # Generate visualizations
  generate_all_visualizations(processed_data)
}

# Execute main pipeline
main()