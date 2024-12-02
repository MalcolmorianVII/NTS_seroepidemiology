# Load required libraries
source("R/libraries.R")

# Load configuration
source("R/config.R")

# Source utility functions
source("R/utils.R")
source("R/data_preprocessing.R")
source("R/statistical_analysis.R")
source("R/visualization.R")

# Main analysis pipeline
main <- function() {
  # Load all datasets
  datasets <- load_saints_datasets()
  saints_seq_only <- datasets$saints_seq_only
  elisa <- datasets$elisa
  age_elisa <- datasets$age_elisa
  wilson_data <- datasets$wilson_data
  
  # Preprocess data
  processed_data <- preprocess_salmonella_data(saints_seq_only)
  
  # Perform analyses
  perform_statistical_analyses(processed_data)
  
  # Generate visualizations
  generate_all_visualizations(processed_data)
}

# Execute main pipeline
main()