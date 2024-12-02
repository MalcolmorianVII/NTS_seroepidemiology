# config.R

# Project-level configuration settings

# File Paths Configuration
CONFIG <- list(
  # Base directory for the project
  base_dir = "/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard",
  
  # Input file paths
  input_files = list(
    main_seq = file.path(
      "/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard", 
      "2024.11.17/work/2024.11.21/2024.11.22.SAINTS_master.xlsx"
    ),
    elisa = file.path(
      "/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard", 
      "2024.11.17/work/2024.11.21/2024.11.21.CQJvsCQQ_elisa_seq.xlsx"
    ),
    age_elisa = file.path(
      "/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard", 
      "2024.11.26/2024.11.26.SAINTS_merged_with_age.xlsx"
    ),
    wilson_data = file.path(
      "/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard", 
      "2024.11.26/cat_wilson_salmonella_cariage_serovar.xlsx"
    )
  ),
  
  # Output directory configuration
  output_dir = file.path(
    "/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard", 
    "outputs"
  ),
  
  # Analysis parameters
  analysis_params = list(
    # Filtering parameters
    min_sample_size = 10,
    exclude_samples = c("CQJ142"),
    
    # Statistical analysis thresholds
    fold_change_threshold = 4,
    p_value_threshold = 0.05
  ),
  
  # Visualization settings
  plot_settings = list(
    theme = "minimal",
    color_palette = c("blue", "red"),
    font_size = 12
  )
)

# Optional: Validation function for configuration
validate_config <- function(config) {
  # Check if all required input files exist
  for (file_path in config$input_files) {
    if (!file.exists(file_path)) {
      stop(paste("Configuration error: File not found -", file_path))
    }
  }
  
  # Validate output directory
  if (!dir.exists(config$output_dir)) {
    warning(paste("Output directory does not exist. Creating:", config$output_dir))
    dir.create(config$output_dir, recursive = TRUE)
  }
  
  return(TRUE)
}

# Validate configuration on load
validate_config(CONFIG)