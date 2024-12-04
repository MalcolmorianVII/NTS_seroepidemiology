library(ggplot2)
library(dplyr)
library(ggpubr)
library(patchwork)
library(scales)

# Visualization Configuration
viz_config <- list(
  theme_base = theme_minimal(),
  color_palette = c(
    "#1F77B4",  # Blue
    "#FF7F0E",  # Orange
    "#2CA02C",  # Green
    "#D62728",  # Red
    "#9467BD"   # Purple
  ),
  font_family = "Arial",
  title_size = 14,
  subtitle_size = 10,
  axis_text_size = 10
)

# Custom ggplot theme
theme_salmonella <- function() {
  theme_minimal() +
  theme(
    text = element_text(family = viz_config$font_family),
    plot.title = element_text(
      size = viz_config$title_size, 
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      size = viz_config$subtitle_size, 
      hjust = 0.5
    ),
    axis.text = element_text(size = viz_config$axis_text_size),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )
}

# Violin Plot with Boxplot Overlay
create_distribution_plot <- function(
  data, 
  x_var, 
  y_var, 
  fill_var = NULL, 
  title = NULL
) {
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }
  
  # Create base plot
  plot <- ggplot(data, aes_string(x = x_var, y = y_var, fill = fill_var)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.2, position = position_dodge(0.9), alpha = 0.7) +
    scale_fill_manual(values = viz_config$color_palette) +
    theme_salmonella() +
    labs(
      title = title %||% paste("Distribution of", y_var),
      x = x_var,
      y = y_var
    )
  
  # Add statistical annotations if fill_var is provided
  if (!is.null(fill_var)) {
    plot <- plot + 
      stat_compare_means(
        aes_string(group = fill_var),
        label = "p.signif",
        method = "t.test"
      )
  }
  
  return(plot)
}

# Scatter Plot with Marginal Distributions
create_scatter_plot <- function(
  data, 
  x_var, 
  y_var, 
  color_var = NULL, 
  title = NULL
) {
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }
  
  # Base scatter plot
  scatter <- ggplot(data, aes_string(x = x_var, y = y_var, color = color_var)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
    scale_color_manual(values = viz_config$color_palette) +
    theme_salmonella() +
    labs(
      title = title %||% paste(y_var, "vs", x_var),
      x = x_var,
      y = y_var
    )
  
  # Marginal distributions
  x_density <- ggplot(data, aes_string(x = x_var)) +
    geom_density(fill = viz_config$color_palette[1], alpha = 0.5) +
    theme_void()
  
  y_density <- ggplot(data, aes_string(x = y_var)) +
    geom_density(fill = viz_config$color_palette[2], alpha = 0.5) +
    coord_flip() +
    theme_void()
  
  # Combine plots
  combined_plot <- x_density + 
    (scatter + theme(legend.position = "none")) + 
    y_density + 
    plot_layout(
      ncol = 2, 
      nrow = 2, 
      widths = c(4, 1), 
      heights = c(1, 4)
    )
  
  return(combined_plot)
}

# Serovar Distribution Bar Plot
create_serovar_distribution_plot <- function(
  data, 
  serovar_var, 
  count_type = "serovar",
  title = NULL
) {
  # Count serovars
  serovar_counts <- data %>%
    group_by(!!sym(serovar_var)) %>%
    summarise(count = n()) %>%
    arrange(desc(count))
  
  # Create bar plot
  plot <- ggplot(serovar_counts, aes(x = !!sym(serovar_var), y = count)) +
    geom_bar(stat = "identity", fill = viz_config$color_palette[1], alpha = 0.7) +
    theme_salmonella() +
    labs(
      title = title %||% paste(count_type, "Distribution"),
      x = serovar_var,
      y = "Count"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(breaks = pretty_breaks())
  
  return(plot)
}

# Antibody Fold Change Heatmap
create_fold_change_heatmap <- function(
  data, 
  group_var, 
  value_var, 
  threshold = 4
) {
  # Prepare data
  heatmap_data <- data %>%
    group_by(!!sym(group_var)) %>%
    summarise(
      mean_fold_change = mean(!!sym(value_var), na.rm = TRUE),
      percent_above_threshold = sum(!!sym(value_var) >= threshold, na.rm = TRUE) / n() * 100
    )
  
  # Create heatmap
  plot <- ggplot(heatmap_data, aes(x = !!sym(group_var), y = "Antibody Response")) +
    geom_tile(aes(fill = mean_fold_change), color = "white") +
    scale_fill_gradient(low = "white", high = viz_config$color_palette[1]) +
    geom_text(aes(label = sprintf("%.1f%%\n(Mean: %.2f)", percent_above_threshold, mean_fold_change))) +
    theme_salmonella() +
    labs(
      title = "Antibody Fold Change Heatmap",
      fill = "Mean Fold Change"
    )
  
  return(plot)
}

# Visualization Wrapper Function
generate_salmonella_visualizations <- function(
  data, 
  group_var = "serogroup", 
  value_var = "Antibody_fold_change"
) {
  # Generate multiple visualizations
  plots <- list(
    distribution = create_distribution_plot(
      data, 
      x_var = group_var, 
      y_var = value_var
    ),
    scatter = create_scatter_plot(
      data, 
      x_var = "age", 
      y_var = value_var
    ),
    serovar_dist = create_serovar_distribution_plot(
      data, 
      serovar_var = "serovar"
    ),
    fold_change_heatmap = create_fold_change_heatmap(
      data, 
      group_var = group_var, 
      value_var = value_var
    )
  )
  
  # Optional: Combine plots using patchwork
  combined_plot <- wrap_plots(plots) + 
    plot_annotation(
      title = "Comprehensive Salmonella Data Visualization",
      theme = theme_salmonella()
    )
  
  return(list(
    individual_plots = plots,
    combined_plot = combined_plot
  ))
}

# Example usage (commented out)
# visualizations <- generate_salmonella_visualizations(elisa_data)