library(openxlsx)
data <- read.xlsx("Data_hospital_patients.xlsx")
library(dplyr)
library(tidyr)
library(broom)
data <- select(data, patient_SN, Drug.Name, 
               Pre_value_Platelet, Post_value_Platelet,
               Pre_value_Creatinine, Post_value_Creatinine,
               Pre_value_Bilirubin, Post_value_Bilirubin,
               Pre_value_MCV, Post_value_MCV,
               Pre_value_RBC, Post_value_RBC)
stats_results_df <- data.frame()
pairs <- list(
  Creatinine = c("Pre_value_Creatinine", "Post_value_Creatinine"),
  Bilirubin = c("Pre_value_Bilirubin", "Post_value_Bilirubin"),
  Platelet = c("Pre_value_Platelet", "Post_value_Platelet"),
  MCV = c("Pre_value_MCV", "Post_value_MCV"),
  RBC = c("Pre_value_RBC", "Post_value_RBC")
)

# Initialize an empty dataframe to store results
stats_results_df <- data.frame()

for (indicator in names(pairs)) {
  pair <- pairs[[indicator]]
  current_data <- data[, c("patient_SN", pair[1], pair[2]), drop = FALSE]
  current_data <- dplyr::filter(current_data, !is.na(current_data[[pair[1]]]) & !is.na(current_data[[pair[2]]]))
  
  # Skip if there's no valid data
  if (nrow(current_data) == 0) next
  
  # Perform the paired Wilcoxon test
  test_without_outliers <- wilcox.test(current_data[[pair[1]]], current_data[[pair[2]]], paired = TRUE)
  
  # Append results to the dataframe
  stats_results_df <- rbind(stats_results_df,
                            data.frame(
                              indicator = indicator,
                              number = nrow(current_data),
                              `average previous value` = mean(current_data[[pair[1]]], na.rm = TRUE),
                              `average post value` = mean(current_data[[pair[2]]], na.rm = TRUE),
                              p_value = test_without_outliers$p.value
                            ))
}

library(dplyr)
library(dplyr)
pre_means <- data %>%
  summarise(
    Avg_Pre_Platelets = mean(Pre_value_Platelet, na.rm = TRUE),
    Avg_Pre_Creatinine = mean(Pre_value_Creatinine, na.rm = TRUE),
    Avg_Pre_Bilirubin = mean(Pre_value_Bilirubin, na.rm = TRUE),
    Avg_Pre_MCV = mean(Pre_value_MCV, na.rm = TRUE),
    Avg_Pre_RBC = mean(Pre_value_RBC, na.rm = TRUE)
  )
data <- data %>%
  mutate(
    Z_Pre_Platelets = sign((Pre_value_Platelet - pre_means$Avg_Pre_Platelets)) * log2(abs((Pre_value_Platelet - pre_means$Avg_Pre_Platelets) / sd(Pre_value_Platelet, na.rm = TRUE)) + 1),
    Z_Post_Platelets = sign((Post_value_Platelet - pre_means$Avg_Pre_Platelets)) * log2(abs((Post_value_Platelet - pre_means$Avg_Pre_Platelets) / sd(Post_value_Platelet, na.rm = TRUE)) + 1),
    Z_Pre_Creatinine = sign((Pre_value_Creatinine - pre_means$Avg_Pre_Creatinine)) * log2(abs((Pre_value_Creatinine - pre_means$Avg_Pre_Creatinine) / sd(Pre_value_Creatinine, na.rm = TRUE)) + 1),
    Z_Post_Creatinine = sign((Post_value_Creatinine - pre_means$Avg_Pre_Creatinine)) * log2(abs((Post_value_Creatinine - pre_means$Avg_Pre_Creatinine) / sd(Post_value_Creatinine, na.rm = TRUE)) + 1),
    Z_Pre_Bilirubin = sign((Pre_value_Bilirubin - pre_means$Avg_Pre_Bilirubin)) * log2(abs((Pre_value_Bilirubin - pre_means$Avg_Pre_Bilirubin) / sd(Pre_value_Bilirubin, na.rm = TRUE)) + 1),
    Z_Post_Bilirubin = sign((Post_value_Bilirubin - pre_means$Avg_Pre_Bilirubin)) * log2(abs((Post_value_Bilirubin - pre_means$Avg_Pre_Bilirubin) / sd(Post_value_Bilirubin, na.rm = TRUE)) + 1),
    Z_Pre_MCV = sign((Pre_value_MCV - pre_means$Avg_Pre_MCV)) * log2(abs((Pre_value_MCV - pre_means$Avg_Pre_MCV) / sd(Pre_value_MCV, na.rm = TRUE)) + 1),
    Z_Post_MCV = sign((Post_value_MCV - pre_means$Avg_Pre_MCV)) * log2(abs((Post_value_MCV - pre_means$Avg_Pre_MCV) / sd(Post_value_MCV, na.rm = TRUE)) + 1),
    Z_Pre_RBC = sign((Pre_value_RBC - pre_means$Avg_Pre_RBC)) * log2(abs((Pre_value_RBC - pre_means$Avg_Pre_RBC) / sd(Pre_value_RBC, na.rm = TRUE)) + 1),
    Z_Post_RBC = sign((Post_value_RBC - pre_means$Avg_Pre_RBC)) * log2(abs((Post_value_RBC - pre_means$Avg_Pre_RBC) / sd(Post_value_RBC, na.rm = TRUE)) + 1)
  )
data <- data %>%
  mutate(
    Score_Pre_Platelets = ifelse((Pre_value_Platelet - pre_means$Avg_Pre_Platelets) < 0, abs(Z_Pre_Platelets), 0),
    Score_Post_Platelets = ifelse((Post_value_Platelet - pre_means$Avg_Pre_Platelets) < 0, abs(Z_Post_Platelets), 0),
    Score_Pre_Creatinine = ifelse((Pre_value_Creatinine - pre_means$Avg_Pre_Creatinine) > 0, Z_Pre_Creatinine, 0),
    Score_Post_Creatinine = ifelse((Post_value_Creatinine - pre_means$Avg_Pre_Creatinine) > 0, Z_Post_Creatinine, 0),
    Score_Pre_Bilirubin = ifelse((Pre_value_Bilirubin - pre_means$Avg_Pre_Bilirubin) > 0, Z_Pre_Bilirubin, 0),
    Score_Post_Bilirubin = ifelse((Post_value_Bilirubin - pre_means$Avg_Pre_Bilirubin) > 0, Z_Post_Bilirubin, 0),
    Score_Pre_MCV = ifelse((Pre_value_MCV - pre_means$Avg_Pre_MCV) > 0, abs(Z_Pre_MCV), 0),
    Score_Post_MCV = ifelse((Post_value_MCV - pre_means$Avg_Pre_MCV) > 0, abs(Z_Post_MCV), 0),
    Score_Pre_RBC = ifelse((Pre_value_RBC - pre_means$Avg_Pre_RBC) < 0, abs(Z_Pre_RBC), 0),
    Score_Post_RBC = ifelse((Post_value_RBC - pre_means$Avg_Pre_RBC) < 0, abs(Z_Post_RBC), 0)
  )
data <- data %>%
  mutate(
    Total_Score_Pre = Score_Pre_Platelets + Score_Pre_Creatinine + Score_Pre_Bilirubin + Score_Pre_MCV + Score_Pre_RBC,
    Total_Score_Post = Score_Post_Platelets + Score_Post_Creatinine + Score_Post_Bilirubin + Score_Post_MCV + Score_Post_RBC
  )
print(data)
library(dplyr)
data <- data %>%
  mutate(
    Risk_Category_Pre = case_when(
      Total_Score_Pre >= 0 & Total_Score_Pre <= 1   ~ "Low Risk",
      Total_Score_Pre > 1 & Total_Score_Pre <= 3 ~ "Moderate-Low Risk",
      Total_Score_Pre > 3 & Total_Score_Pre <= 5 ~ "Moderate Risk",
      Total_Score_Pre > 5 ~ "High Risk"
    ),
    Risk_Category_Post = case_when(
      Total_Score_Post >= 0 & Total_Score_Post <= 1   ~ "Low Risk",
      Total_Score_Post > 1 & Total_Score_Post <= 3 ~ "Moderate-Low Risk",
      Total_Score_Post > 3 & Total_Score_Post <= 5 ~ "Moderate Risk",
      Total_Score_Post > 5 ~ "High Risk"
    )
  )
risk_table <- table(data$Risk_Category_Pre, data$Risk_Category_Post)
print(risk_table)
library(dplyr)
risk_changes <- data %>%
  transmute(
    Pre = Risk_Category_Pre,
    Post = Risk_Category_Post
  )
risk_changes$Pre <- factor(risk_changes$Pre, levels = c("Low Risk", "Moderate-Low Risk", "Moderate Risk", "High Risk"), ordered = TRUE)
risk_changes$Post <- factor(risk_changes$Post, levels = c("Low Risk", "Moderate-Low Risk", "Moderate Risk", "High Risk"), ordered = TRUE)
risk_changes <- risk_changes %>%
  mutate(
    Type = case_when(
      Pre == Post ~ "No Change",
      Pre < Post ~ "Increase",
      TRUE ~ "Decrease"
    )
  )
print(risk_changes)
risk_change_summary <- risk_changes %>%
  count(Pre, Post, Type, sort = TRUE)
print(risk_change_summary)

library(DescTools)
g_test_result <- GTest(risk_table)
print(g_test_result)
merged_data_Crea <- data
merged_data_Crea <- select(merged_data_Crea, patient_SN, Drug.Name, Pre_value_Platelet, Post_value_Platelet,
                           Pre_value_Creatinine, Post_value_Creatinine)
merged_data_Crea <- merged_data_Crea %>%
  mutate(Pre_Match = ifelse(Pre_value_Platelet < 30 & Pre_value_Creatinine < 199, 1, 0),
         Post_Match = ifelse(Post_value_Platelet < 30 & Post_value_Creatinine < 199, 1, 0))
merged_data_table <- table(merged_data_Crea$Pre_Match, merged_data_Crea$Post_Match)
print(merged_data_table)
merged_data_chi_sq_test <- fisher.test(merged_data_table)
print(merged_data_chi_sq_test)
bevacizumab_data <- merged_data_Crea %>%
  filter(Drug.Name == "Bevacizumab")
bevacizumab_table <- table(bevacizumab_data$Pre_Match, bevacizumab_data$Post_Match)
print(bevacizumab_table)
bevacizumab_fisher_test <- fisher.test(bevacizumab_table)
print(bevacizumab_fisher_test)
other_drugs_data <- merged_data_Crea %>%
  filter(Drug.Name != "Bevacizumab")
other_drugs_table <- table(other_drugs_data$Pre_Match, other_drugs_data$Post_Match)
print(other_drugs_table)
other_drugs_fisher_test <- fisher.test(other_drugs_table)
print(other_drugs_fisher_test)

