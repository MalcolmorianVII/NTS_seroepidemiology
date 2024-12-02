library(dplyr)
library(ggplot2)
library(stats)

# Statistical Analysis Module
salmonella_stats_analysis <- function(data, 
                                      group_variable, 
                                      value_variable, 
                                      fold_change_threshold = 4) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }
  
  if (!(group_variable %in% names(data))) {
    stop(paste("Group variable", group_variable, "not found in data"))
  }
  
  if (!(value_variable %in% names(data))) {
    stop(paste("Value variable", value_variable, "not found in data"))
  }
  
  # Perform t-test
  t_test_result <- tryCatch({
    t.test(as.formula(paste(value_variable, "~", group_variable)), data = data)
  }, error = function(e) {
    warning("T-test failed: ", e$message)
    return(NULL)
  })
  
  # Group-wise summary statistics
  summary_stats <- data %>%
    group_by(!!sym(group_variable)) %>%
    summarise(
      mean_value = mean(!!sym(value_variable), na.rm = TRUE),
      median_value = median(!!sym(value_variable), na.rm = TRUE),
      sd_value = sd(!!sym(value_variable), na.rm = TRUE),
      count_total = n(),
      count_above_threshold = sum(!!sym(value_variable) >= fold_change_threshold, na.rm = TRUE),
      percentage_above_threshold = (count_above_threshold / count_total) * 100
    )
  
  # Visualization: Boxplot with statistical comparison
  boxplot_visualization <- ggplot(data, aes_string(x = group_variable, y = value_variable)) +
    geom_boxplot(aes_string(fill = group_variable), alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    stat_compare_means(method = "t.test", label = "p.signif") +
    labs(
      title = paste("Distribution of", value_variable, "by", group_variable),
      x = group_variable,
      y = value_variable
    ) +
    theme_minimal()
  
  # Effect size calculation (Cohen's d)
  effect_size <- tryCatch({
    group_means <- data %>% 
      group_by(!!sym(group_variable)) %>% 
      summarise(group_mean = mean(!!sym(value_variable), na.rm = TRUE))
    
    pooled_sd <- data %>%
      group_by(!!sym(group_variable)) %>%
      summarise(group_sd = sd(!!sym(value_variable), na.rm = TRUE)) %>%
      summarise(pooled_sd = sqrt(mean(group_sd^2))) %>%
      pull(pooled_sd)
    
    diff_means <- diff(group_means$group_mean)
    cohen_d <- abs(diff_means / pooled_sd)
    
    cohen_d
  }, error = function(e) {
    warning("Effect size calculation failed: ", e$message)
    return(NA)
  })
  
  # Construct results object
  results <- list(
    t_test = t_test_result,
    summary_statistics = summary_stats,
    visualization = boxplot_visualization,
    effect_size = effect_size
  )
  
  return(results)
}

# Example usage
# elisa_analysis <- salmonella_stats_analysis(
#   data = elisa_data, 
#   group_variable = "serogroup", 
#   value_variable = "Antibody_fold_change"
# )