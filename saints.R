library(tidyr)
library(dplyr)
library(ggplot2)
library(readxl)
library(ggpubr)
library(scales)
library(vegan)
library(purrr)

# -------------- GENOMIC DATASET ANALYSIS ------------------------------
saints_seq_only <- read_excel("/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard/2024.11.17/work/2024.11.21/2024.11.22.SAINTS_master.xlsx")
# Exclude CQJ142 ( all 3 antigens absent & serovar prediction extremely ambiguous)
saints_seq_only <- saints_seq_only %>%
  filter(!(sample == "CQJ142" | is.na(serogroup)))

# Split saint_seq_only into Subspecies I and II
subspecies1 <- saints_seq_only %>% filter(subspecies == "I")
subspecies2 <- saints_seq_only %>% filter(subspecies == "II")
subspecies1_total <- nrow(subspecies1)
subspecies2_total <- nrow(subspecies2)

#Get singletons
get_singletons <- function(data,classification) {
  data %>%
    group_by(classification) %>%
    summarise(count = n()) %>%
    filter(count == 1)
}
singletons1_O_antigen <- get_singletons(subspecies1,"serovar")
singletons1_cgMLST <- get_singletons(subspecies1,"serovar_cgmlst")
singletons2_O_antigen <- get_singletons(subspecies2)
singletons2_cgMLST <- get_singletons(subspecies2,"serovar_cgmlst")
# create a data frame of all singletons and write to an Excel file
singletons_df <- rbind(singletons1_O_antigen, singletons1_cgMLST, singletons2_O_antigen, singletons2_cgMLST)
write.xlsx(singletons_df, "/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard/2024.11.26/serovar_singletons.xlsx")

# --------- REMOVE SINGLETONS -------------
remove_singletons <- function(data,classification) {
  data %>%
    group_by(serovar) %>%
    filter(n() > 1)
}
subspecies1_serovar <- remove_singletons(subspecies1,"serovar")
subspecies1_cgMLST <- remove_singletons(subspecies1,"serovar_cgmlst")
subspecies2_serovar <- remove_singletons(subspecies2,"serovar")
subspecies2_cgMLST <- remove_singletons(subspecies2,"serovar_cgmlst")

# Reorder the 'serovar' & 'serogroup' columns by its frequency in descending order
order_cols <- function(data) {
  data %>%
    mutate(serovar = factor(serovar, levels = names(sort(table(serovar), decreasing = TRUE))))
}
subspecies1_serovar$serovar <- factor(subspecies1_serovar$serovar, 
                              levels = names(sort(table(subspecies1_serovar$serovar), decreasing = TRUE)))
subspecies2_serovar$serovar <- factor(subspecies2_serovar$serovar, 
                              levels = names(sort(table(subspecies2_serovar$serovar), decreasing = TRUE)))
subspecies1$serogroup <- factor(subspecies1$serogroup, 
                              levels = names(sort(table(subspecies1$serogroup), decreasing = TRUE)))

subspecies2$serogroup <- factor(subspecies2$serogroup, 
                                levels = names(sort(table(subspecies2$serogroup), decreasing = TRUE)))

subspecies1_cgMLST$serovar_cgmlst <- factor(subspecies1_cgMLST$serovar_cgmlst, 
                              levels = names(sort(table(subspecies1_cgMLST$serovar_cgmlst), decreasing = TRUE)))
subspecies2_cgMLST$serovar_cgmlst <- factor(subspecies2_cgMLST$serovar_cgmlst, 
                              levels = names(sort(table(subspecies2_cgMLST$serovar_cgmlst), decreasing = TRUE)))

# SEROGROUP DISTRIBUTION ACROSS SUBSPECIES I & II
ggplot(subspecies1, aes(x = serogroup)) +
  geom_bar(position = "dodge", fill = "blue") + # Set bar fill color to blue
  labs(
    title = paste("In silico serogroup distribution in Subspecies I (n=", subspecies1_total, ")"),
    x = "Subspecies I",
    y = "Count",
    fill = "Serogroup"
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + # Ensure integer breaks
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Adjust x-axis labels
    panel.grid = element_blank() # Remove grid lines
  )

ggplot(subspecies2, aes(x = serogroup)) +
  geom_bar(position = "dodge", fill = "blue") + # Set bar fill color to blue
  labs(
    title = paste("Serogroup distribution in Subspecies II (n=", subspecies2_total, ")"),
    x = "Subspecies II",
    y = "Count",
    fill = "Serogroup"
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + # Ensure integer breaks
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Adjust x-axis labels
    panel.grid = element_blank() # Remove grid lines
  )
# Serovar comparison _ Antibody profile SUBSPECIES I
ggplot(subspecies1, aes(x = serovar)) +
  geom_bar(position = "dodge", fill = "blue") + # Set bar fill color to blue
  labs(
    title = paste("Serovar distribution(O-antigen) in Subspecies I (n=", subspecies1_total, ")"),
    x = "Subspecies I",
    y = "Count",
    fill = "Serovar"
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + # Ensure integer breaks
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Adjust x-axis labels
    panel.grid = element_blank() # Remove grid lines
  )


# -----------
# Read elisa data
elisa <- read_excel("/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard/2024.11.17/work/2024.11.21/2024.11.21.CQJvsCQQ_elisa_seq.xlsx")
age_elisa <-  read_excel("/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard/2024.11.26/2024.11.26.SAINTS_merged_with_age.xlsx")

wilson_data <- read_excel("/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard/2024.11.26/cat_wilson_salmonella_cariage_serovar.xlsx")
# Jaccard data
# Calculate Jaccard similarity for the two sets
jaccard_similarity <- function(set1, set2) {
  intersection_size <- length(intersect(set1, set2))
  union_size <- length(union(set1, set2))
  return(intersection_size / union_size)
}

# Apply function to calculate Jaccard similarity
jaccard_similarity_value <- jaccard_similarity(saints_seq_only$serovar, wilson_data$serovar)

# Print the result
print(jaccard_similarity_value)

# ......
list1 <- intersect(saints_seq_only$serovar, wilson_data$serovar)

# List 2: Values only in saints_seq_only
list2 <- setdiff(saints_seq_only$serovar, wilson_data$serovar)

# List 3: Values only in wilson_data
list3 <- setdiff(wilson_data$serovar, saints_seq_only$serovar)

# Combine the lists into a data frame for Excel output
merged_data <- data.frame(
  Common_Serovars = c(list1, rep(NA, length(list2) + length(list3))),
  Saints_Only_Serovars = c(rep(NA, length(list1)), list2, rep(NA, length(list3))),
  Wilson_Only_Serovars = c(rep(NA, length(list1) + length(list2)), list3)
)

# Write the merged data to an Excel file
write.xlsx(merged_data, "/Users/malcolmorian/Documents/Bioinformatics/Projects2024/SAINTS/Richard/2024.11.26/serovar_lists.xlsx")

# AGE RELATED ANALYES
age_elisa$serogroup <- ifelse(age_elisa$serogroup == "B", "B", "Non-B")
age_elisa$classification <- ifelse(grepl("1,4,\\[5\\],12,\\[27\\]", age_elisa$serovar_cgmlst), "1,4,[5],12,[27]", "Others")
age_elisa$TyphimuriumsVsNon <- ifelse(grepl("Typhimurium", age_elisa$serovar_cgmlst), "Typhimurium", "Non-Typhimurium")
age_elisa$ageyr <- trunc(age_elisa$ageyr)
# Age vs Antibody_fold change
# age in months
ggplot(age_elisa, aes(x = agemon, y = Antibody_fold_change, color = classification)) +
  geom_point(alpha = 0.7) + # Scatterplot with semi-transparent points
  labs(
    title = "Age vs Antibody Fold Change",
    x = "Age (Months)",
    y = "Antibody Fold Change"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red")) # Optional: custom color palette for "TyphimuriumsVsNon"

# age in years
ggplot(age_elisa, aes(x = factor(ageyr), y = Antibody_fold_change)) +
  geom_boxplot(alpha = 0.1, fill = "blue", color = "black") +  # Very transparent boxplot
  geom_point(position = position_jitter(width = 0.2), alpha = 0.7, color = "blue") +  # Scatter plot with jittering
  labs(
    title = "Age vs Antibody Fold Change",
    x = "Age (Years)",
    y = "Antibody Fold Change"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for better readability
  )

# ------- LOG FOLD CHANGES -------------------------------------------#
elisa$serogroup <- ifelse(elisa$serogroup == "B", "B", "Non-B")
saints_seq_only$subspecies <- ifelse(grepl("II", saints_seq_only$serovar_cgmlst), "II", "I")
elisa$classification <- ifelse(grepl("1,4,\\[5\\],12,\\[27\\]", elisa$serovar_cgmlst), "1,4,[5],12,[27]", "Others")
elisa$TyphimuriumsVsNon <- ifelse(grepl("Typhimurium", elisa$serovar_cgmlst), "Typhimurium", "Non-Typhimurium")
#elisa$Antibody_fold_change <- log10(elisa$Antibody_fold_change)

# Typhimuriums vs Non-Typhimuriums i.e not applicable cuurently
ggplot(elisa, aes(x = TyphimuriumsVsNon, y = Antibody_fold_change, color = TyphimuriumsVsNon)) +
  geom_point(position = position_jitter(width = 0.2)) + # Jitter to spread out points
  #geom_hline(yintercept = 4, linetype = "dotted", color = "red", size = 1) + # Add dotted line at 4x
  labs(
    title = "Log(10) Antibody Fold Change by cgMLST",
    x = "Typhimuriums vs Non-Typhimuriums",
    y = "Antibody Fold Change"
  ) +
  theme_minimal() +
  stat_compare_means(
    method = "t.test",            # Specify the statistical test
    label = "p.signif",           # Use significance stars
    label.x = 1.5,                # Position of the label on the x-axis
    comparisons = list(c("Typhimurium", "Non-Typhimurium")) # Specify groups to compare
  )

# Plot histogram of Antibody Fold Change Typhimuriums vs Non-Typhimuriums
t_test_result_muriums <- t.test(Antibody_fold_change ~ TyphimuriumsVsNon, data = elisa)
ggplot(elisa, aes(x = Antibody_fold_change, fill = TyphimuriumsVsNon)) +
  geom_histogram(position = "dodge", bins = 30, alpha = 0.7) +
  labs(
    title = "Absolute Antibody Fold Change Typhimuriums vs Non-Typhimuriums(cgMLST)",
    x = "Antibody Fold Change",
    y = "Count"
  ) +
  theme_minimal() +
  annotate(
    "text", 
    x = max(elisa$Antibody_fold_change, na.rm = TRUE) * 0.9, # Place near the right edge of the x-axis
    y = max(table(cut(elisa$Antibody_fold_change, breaks = 30))) + 5, # Slightly above the max count
    label = paste0("p = ", signif(t_test_result_muriums$p.value, 3)), # Add the p-value
    color = "red",
    size = 5 # Adjust text size
  )


# Create the scatter plot serogroup B vs others
ggplot(elisa, aes(x = serogroup, y = Antibody_fold_change, color = serogroup)) +
  geom_point(position = position_jitter(width = 0.2)) + # Jitter to spread out points
  geom_hline(yintercept = 4, linetype = "dotted", color = "red", size = 1) + # Add dotted line at 4X
  labs(
    title = "Absolute Antibody Fold Change by Serogroup",
    x = "Serogroup (B vs Non-B)",
    y = "Antibody Fold Change"
  ) +
  theme_minimal() +
  stat_compare_means(
    method = "t.test",            # Specify the statistical test
    label = "p.signif",           # Use significance stars for labels
    label.x = 1.5,                # Position the label on the x-axis
    comparisons = list(c("B", "Non-B")) # Specify groups to compare
  )

# Histogram for serogroup B vs others
t_test_result_serogroup <- t.test(Antibody_fold_change ~ serogroup, data = elisa)
ggplot(elisa, aes(x = Antibody_fold_change, fill = serogroup)) +
  geom_histogram(position = "dodge", bins = 30, alpha = 0.7) +
  labs(
    title = "Absolute Antibody Fold Change by Serogroup B vs Others",
    x = "Antibody Fold Change",
    y = "Count"
  ) +
  theme_minimal() +
  # Dynamically calculate the y-position based on the max count
  annotate(
    "text", 
    x = max(elisa$Antibody_fold_change, na.rm = TRUE) * 0.8, # Adjust x position (80% of the max)
    y = max(table(cut(elisa$Antibody_fold_change, breaks = 30))) + 2, # Adjust y position slightly above the histogram
    label = paste0("p = ", signif(t_test_result_serogroup$p.value, 3)), # Display the p-value
    color = "red",
    size = 5 # Adjust text size
  )

#---------
ggplot(elisa, aes(x = Antibody_fold_change, fill = serogroup)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "dodge") +
  labs(
    title = "Antibody Fold Change Distribution by Serogroup",
    x = "Antibody Fold Change",
    y = "Count"
  ) +
  theme_minimal()

# Create the scatter plot 1,4,[5],12,[27] vs others
ggplot(elisa, aes(x = classification, y = Antibody_fold_change, color = classification)) +
  geom_point(position = position_jitter(width = 0.2)) + # Jitter to spread out points
  geom_hline(yintercept = 4, linetype = "dotted", color = "red", size = 1) + # Add dotted line at 4X
  labs(
    title = "Absolute Antibody Fold Change by 1,4,[5],12,[27] vs Others (cgMLST)",
    x = "1,4,[5],12,[27] vs Others",
    y = "Antibody Fold Change"
  ) +
  theme_minimal() +
  stat_compare_means(
    method = "t.test",            # Specify the statistical test
    label = "p.signif",           # Use significance stars for labels
    label.x = 1.5,                # Position the label on the x-axis
    comparisons = list(c("1,4,[5],12,[27]", "Others")) # Specify groups to compare
  )


# Perform t-test
t_test_result <- t.test(Antibody_fold_change ~ classification, data = elisa)

# Histogram
ggplot(elisa, aes(x = Antibody_fold_change, fill = classification)) +
  geom_histogram(bins = 30, position = "dodge", alpha = 0.7) + # Dodged histogram
  labs(
    title = "Absolute Antibody Fold Change 1,4,[5],12,[27] vs Others (cgMLST)",
    x = "Antibody Fold Change",
    y = "Count"
  ) +
  theme_minimal() +
  annotate(
    "text", 
    x = max(elisa$Antibody_fold_change, na.rm = TRUE) * 0.9, # Place near the right edge of the x-axis
    y = max(table(cut(elisa$Antibody_fold_change, breaks = 30))) + 5, # Slightly above the max count
    label = paste0("p = ", signif(t_test_result$p.value, 3)), # Add the p-value
    color = "red",
    size = 5 # Adjust text size for better visibility
  )

# Summarize fold changes by comparison groups and add the line where the X 4 FC
elisa_summary <- elisa %>%
  group_by(TyphimuriumsVsNon) %>%
  summarise(
    count_at_least_4 = sum(Antibody_fold_change >= 4, na.rm = TRUE), # Count values >= 4
    total = n(),                                                     # Total data points in each group
    percentage = (count_at_least_4 / total) * 100                    # Calculate percentage
  )
print(elisa_summary)

# Summarize fold changes by comparison groups and add the line where the X 4 FC
elisa_summary <- elisa %>%
  group_by(serogroup) %>%
  summarise(
    count_at_least_4 = sum(Antibody_fold_change >= 4, na.rm = TRUE), # Count values >= 4
    total = n(),                                                     # Total data points in each group
    percentage = (count_at_least_4 / total) * 100                    # Calculate percentage
  )
print(elisa_summary)


elisa_summary <- elisa %>%
  group_by(classification) %>%
  summarise(
    count_at_least_4 = sum(Antibody_fold_change >= 4, na.rm = TRUE), # Count values >= 4
    total = n(),                                                     # Total data points in each group
    percentage = (count_at_least_4 / total) * 100                    # Calculate percentage
  )
print(elisa_summary)



# Serovar comparison _ Antibody profile SUBSPECIES II
ggplot(subspecies2, aes(x = serovar)) +
  geom_bar(position = "dodge", fill = "blue") + # Set bar fill color to blue
  labs(
    title = paste("Serovar distribution(O-antigen) in Subspecies II (n=", subspecies2_total, ")"),
    x = "Subspecies II",
    y = "Count",
    fill = "Serovar"
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + # Ensure integer breaks
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Adjust x-axis labels
    panel.grid = element_blank() # Remove grid lines
  )

# Reorder the 'serovar' column by its frequency in descending order
subspecies1$serovar_cgmlst <- factor(subspecies1$serovar_cgmlst, 
                                     levels = names(sort(table(subspecies1$serovar_cgmlst), decreasing = TRUE)))

subspecies2$serovar_cgmlst <- factor(subspecies2$serovar_cgmlst, 
                                     levels = names(sort(table(subspecies2$serovar_cgmlst), decreasing = TRUE)))
# Serovar comparison _ cgMLST SUBSPECIES I
ggplot(subspecies1, aes(x = serovar_cgmlst)) +
  geom_bar(position = "dodge", fill = "blue") + # Set bar fill color to blue
  labs(
    title = paste("Serovar distribution(cgMLST) in Subspecies I (n=", subspecies1_total, ")"),
    x = "Subspecies I",
    y = "Count",
    fill = "Serovar"
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + # Ensure integer breaks
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Adjust x-axis labels
    panel.grid = element_blank() # Remove grid lines
  )

# Serovar comparison _ cgMLST SUBSPECIES II
ggplot(subspecies2, aes(x = serovar_cgmlst)) +
  geom_bar(position = "dodge", fill = "blue") + # Set bar fill color to blue
  labs(
    title = paste("Serovar distribution(cgMLST) in Subspecies II (n=", subspecies2_total, ")"),
    x = "Subspecies II",
    y = "Count",
    fill = "Serovar"
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + # Ensure integer breaks
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Adjust x-axis labels
    panel.grid = element_blank() # Remove grid lines
  )

# Function to calculate concordance
calculate_concordance <- function(data) {
  # Convert factors to character before comparison
  matches <- sum(as.character(data$serovar) == as.character(data$serovar_cgmlst)) # Count matches
  total <- nrow(data)                                                           # Total observations
  percentage_concordance <- (matches / total) * 100
  return(percentage_concordance)
}


# Calculate concordance for subspecies 1 and 2
concordance_subsp1 <- calculate_concordance(subspecies1)
concordance_subsp2 <- calculate_concordance(subspecies2)

# Print results
cat("Concordance for Subspecies 1:", concordance_subsp1, "%\n")
cat("Concordance for Subspecies 2:", concordance_subsp2, "%\n")

# VISUALIZE CONCORDANCE OPTIONAL 

# Combine concordance results
concordance_df <- data.frame(
  subspecies = c("Subspecies 1", "Subspecies 2"),
  concordance = c(concordance_subsp1, concordance_subsp2)
)

# Plot concordance
ggplot(concordance_df, aes(x = subspecies, y = concordance, fill = subspecies)) +
  geom_bar(stat = "identity", color = "black") +
  labs(title = "Percentage Concordance by Subspecies (Antibody gene Vs cgMLST)",
       x = "Subspecies",
       y = "Percentage Concordance (%)") +
  theme_minimal()
