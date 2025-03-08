library(openxlsx)
data <- read.xlsx("Data_hospital_patients.xlsx")

# Select necessary columns
data <- select(data, patient_SN, Drug.Name, 
               Pre_value_Platelet, Post_value_Platelet,
               Pre_value_Creatinine, Post_value_Creatinine,
               Pre_value_Bilirubin, Post_value_Bilirubin,
               Pre_value_MCV, Post_value_MCV,
               Pre_value_RBC, Post_value_RBC)

library(tidyr)
library(dplyr)
library(ggplot2)

# Create a function to plot violin plots for the specified indicator
plot_violin_for_indicator <- function(data, pre_column, post_column, indicator_name) {
  upper_limit_pre <- quantile(data[[pre_column]], 0.95, na.rm = TRUE)
  upper_limit_post <- quantile(data[[post_column]], 0.95, na.rm = TRUE)
  
  data_selected <- data %>%
    select(patient_SN, Drug.Name, all_of(pre_column), all_of(post_column)) %>%
    filter(!!sym(pre_column) <= upper_limit_pre, !!sym(post_column) <= upper_limit_post) %>%
    gather(key = "condition", value = "value", -patient_SN, -Drug.Name) %>%
    mutate(Treatment = ifelse(grepl("Pre", condition), "Pre-Treatment", "Post-Treatment")) %>%
    separate(condition, into = c("Indicator", "Temp"), sep = "_") %>%
    select(-Temp) %>%
    mutate(Condition = interaction(Indicator, Treatment, sep = "\n")) %>%
    mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
    group_by(Condition) %>%
    mutate(BP = value) %>%
    ungroup()
  
  colors <- c("Pre-Treatment" = "#B5CDE1", "Post-Treatment" = "#F6D2CC")
  
  sampling_rate <- 0.50
  
  sampled_data <- data_selected %>%
    group_by(Condition) %>%
    sample_frac(sampling_rate) %>%
    ungroup()
  
  p <- ggplot(data_selected, aes(x = Condition, y = BP, fill = Treatment)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.4, notch = TRUE, outlier.size = -1, color = "black", lwd = 0.8, alpha = 1) +
    geom_jitter(data = sampled_data, width = 0.15, size = 1, alpha = 0.6, shape = 21, color = "black") +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.margin = margin(b = 5),
      legend.box.margin = margin(b = 5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 8, color = "black"),
      axis.text.y = element_text(size = 7, color = "black"),
      legend.title = element_text(size = 10, color = "black"),
      legend.text = element_text(size = 8, color = "black"),
      plot.title = element_text(size = 10, hjust = 0.5),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      axis.ticks.margin = unit(3, "mm"),
      axis.line = element_line(color = "black")
    ) +
    labs(title = sprintf("Comparison of %s Before and After Treatment", indicator_name),
         x = "Condition",
         y = "Measured Value")
  
  print(p)
}

# Call the function to plot each indicator
indicators <- c("Platelet","Creatinine", "Bilirubin", "MCV", "RBC")
pre_cols <- c("Pre_value_Platelet","Pre_value_Creatinine", "Pre_value_Bilirubin", "Pre_value_MCV", "Pre_value_RBC")
post_cols <- c("Post_value_Platelet","Post_value_Creatinine", "Post_value_Bilirubin", "Post_value_MCV", "Post_value_RBC")

for (i in seq_along(indicators)) {
  plot_violin_for_indicator(data, pre_cols[i], post_cols[i], indicators[i])
}
