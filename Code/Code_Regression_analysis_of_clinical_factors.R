load("Demo_VEGFi.Rdata")
library(dplyr)
data <- Demo_VEGFi
drug_list <- c("bevacizumab.all", "ranibizumab.all", "brolucizumab.all", "aflibercept.all", "conbercept.all", "pegaptanib.all", "ramucirumab.all", "nintedanib.all", "apatinib.all", "axitinib.all", "sunitinib.all", "sorafenib.all", "regorafenib.all", "cabozantinib.all", "pazopanib.all", "lenvatinib.all", "anlotinib.all", "fruquintinib.all", "tivozanib.all", "cediranib.all", "brivanib.all", "vandetanib.all")
classify_drug <- function(row) {
  for (drug in drug_list) {
    if (row[drug] == 1) {
      return(sub(".all", "", drug))
    }
  }
  return(NA)
}
data$Drug_Group <- apply(data, 1, classify_drug)
table(data$Drug_Group)
drug_type <- c("vegf.all", "vegfr.all")
classify_drug_type <- function(row) {
  for (drug in drug_type) {
    if (row[drug] == 1) {
      return(sub(".all", "", drug))
    }
  }
  return(NA)
}
data$Drug_Type <- apply(data, 1, classify_drug_type)
table(data$Drug_Type)
Demo <- data
blood_pt <- read.xlsx("PT_直接相关.xlsx")
target_pt <- unique(blood_pt$PT)
target_pt <- tolower(target_pt)
setDT(Demo)
Demo$TargetPT.all <- as.integer(sapply(Demo$PT_split, function(pt_array) {
  any(target_pt %in% pt_array)
}))
table(Demo$TargetPT.all)
library(dplyr)
Demo <- Demo %>% 
  select(-c("bevacizumab.all", "ranibizumab.all", "brolucizumab.all", "aflibercept.all",
            "conbercept.all", "pegaptanib.all", "ramucirumab.all", "nintedanib.all",
            "apatinib.all", "axitinib.all", "sunitinib.all", "sorafenib.all", 
            "regorafenib.all", "cabozantinib.all", "pazopanib.all", "lenvatinib.all",
            "anlotinib.all", "fruquintinib.all", "tivozanib.all", "cediranib.all",
            "brivanib.all", "vandetanib.all", "vegf.all", "vegfr.all", "all_drug.all", "chemo.all"))
names(Demo)
write.xlsx(Demo, "Demo_vegfi_clean.xlsx")
save(Demo, file = "Demo_vegfi_clean.Rdata")
load("Demo_vegfi_clean.Rdata")
data <- Demo
data <- Demo
filtered_data <- data[data$SEX %in% c("F", "M"), ]
filtered_data$group <- as.numeric(filtered_data$SEX == "F")
model_F_as_control <- glm(TargetPT.all ~ group, family = binomial(link = "logit"), data = filtered_data)
summary(model_F_as_control)
coef_est <- coef(summary(model_F_as_control))["group", "Estimate"]
coef_se <- coef(summary(model_F_as_control))["group", "Std. Error"]
odds_ratio <- exp(coef_est)
lower_ci <- exp(coef_est - 1.96 * coef_se)
upper_ci <- exp(coef_est + 1.96 * coef_se)
p_value <- coef(summary(model_F_as_control))["group", "Pr(>|z|)"]
sex_results_df <- data.frame(
  Term = "Male vs Female",
  OddsRatio = odds_ratio,
  LowerCI = lower_ci,
  UpperCI = upper_ci,
  PValue = p_value
)
print(sex_results_df)
write.xlsx(sex_results_df, "Sex.xlsx")
process_age_data <- function(data, age_column = "AGE", age_code_column = "AGE_COD", 
                             breaks = c(0, 45, 59, 70, Inf), labels = c("0-45", "46-59", "60-70", "70+")) {
  data <- data[!is.na(data[[age_column]]) & !is.na(data[[age_code_column]]), ]
  data[[age_column]] <- ifelse(data[[age_code_column]] == "YR", as.numeric(data[[age_column]]),
                               ifelse(data[[age_code_column]] == "DY", as.numeric(data[[age_column]]) / 365,
                                      ifelse(data[[age_code_column]] == "MON", as.numeric(data[[age_column]]) / 12,
                                             ifelse(data[[age_code_column]] == "DEC", as.numeric(data[[age_column]]) * 10, NA))))
  data <- data[data[[age_column]] > 1, ]
  mean_age <- mean(data[[age_column]], na.rm = TRUE)
  sd_age <- sd(data[[age_column]], na.rm = TRUE)
  median_age <- median(data[[age_column]], na.rm = TRUE)
  min_age <- min(data[[age_column]], na.rm = TRUE)
  max_age <- max(data[[age_column]], na.rm = TRUE)
  q1_age <- quantile(data[[age_column]], 0.25, na.rm = TRUE)
  q3_age <- quantile(data[[age_column]], 0.75, na.rm = TRUE)
  age_distribution <- table(cut(data[[age_column]], breaks = breaks, labels = labels, include.lowest = TRUE))
  age_distribution_percent <- prop.table(age_distribution) * 100
  age_stats <- data.frame(
    `Age` = c("Mean (SD)", "Median [Min, Max]", "Q1 [25th Percentile]", "Q3 [75th Percentile]"),
    `Value` = c(sprintf("%.1f (%.1f)", mean_age, sd_age), 
                sprintf("%.1f [%.1f, %.1f]", median_age, min_age, max_age),
                sprintf("%.1f", q1_age),
                sprintf("%.1f", q3_age))
  )
  age_group_distribution <- data.frame(
    `Age group` = labels,
    `Count` = as.integer(age_distribution),
    `Percentage` = sprintf("%.1f%%", age_distribution_percent)
  )
  list(Cleaned_Data = data, Age_Stats = age_stats, Age_Group_Distribution = age_group_distribution)
}
process_weight_data <- function(data, weight_column = "WT", weight_code_column = "WT_COD") {
  data[[weight_column]][data[[weight_column]] == ""] <- NA
  data[[weight_code_column]][data[[weight_code_column]] == ""] <- NA
  data_weight <- subset(data, !is.na(data[[weight_column]]) & !is.na(data[[weight_code_column]]))
  data_weight[[weight_column]][data_weight[[weight_code_column]] == "LBS"] <- as.numeric(data_weight[[weight_column]][data_weight[[weight_code_column]] == "LBS"]) * 0.453592
  data_weight[[weight_column]] <- as.numeric(data_weight[[weight_column]])
  quartiles <- quantile(data_weight[[weight_column]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  median <- median(data_weight[[weight_column]], na.rm = TRUE)
  mean <- mean(data_weight[[weight_column]], na.rm = TRUE)
  weight_stats <- data.frame(
    `Weight` = c("Q1 [25th Percentile]", "Median [50th Percentile]", "Q3 [75th Percentile]", "Mean"),
    `Value` = c(sprintf("%.1f", quartiles[1]), 
                sprintf("%.1f", quartiles[2]),
                sprintf("%.1f", quartiles[3]),
                sprintf("%.1f", mean))
  )
  list(Cleaned_Data = data_weight, Weight_Stats = weight_stats)
}
result_age <- process_age_data(data)
result_weight <- process_weight_data(data)
cleaned_age_data <- result_age$Cleaned_Data
cleaned_weight_data <- result_weight$Cleaned_Data
print(result_age$Age_Stats)
print(result_age$Age_Group_Distribution)
print(result_weight$Weight_Stats)
cleaned_age_data$Age_Group <- cut(cleaned_age_data$AGE, 
                                  breaks = c(-Inf, 45, 59, 70, Inf), 
                                  labels = c("0-45", "46-59", "60-70", ">70"))
cleaned_weight_data$Weight_Group <- cut(cleaned_weight_data$WT, 
                                        breaks = c(-Inf, 60, 70, 80, Inf), 
                                        labels = c("0-60", "60-70", "70-80", ">80"))
table(cleaned_age_data$Age_Group)
table(cleaned_weight_data$Weight_Group)
cleaned_age_data$Age_Group <- relevel(cleaned_age_data$Age_Group, ref = "0-45")
age_model <- glm(TargetPT.all ~ Age_Group, family = binomial(link = "logit"), data = cleaned_age_data)
age_coefs <- summary(age_model)$coefficients
age_or <- exp(age_coefs[, "Estimate"])
age_ci_lower <- exp(age_coefs[, "Estimate"] - 1.96 * age_coefs[, "Std. Error"])
age_ci_upper <- exp(age_coefs[, "Estimate"] + 1.96 * age_coefs[, "Std. Error"])
age_p_values <- age_coefs[, "Pr(>|z|)"]
age_results_df <- data.frame(
  Group = rownames(age_coefs),
  OddsRatio = age_or,
  LowerCI = age_ci_lower,
  UpperCI = age_ci_upper,
  PValue = age_p_values
)
print(age_results_df)
weight_model <- glm(TargetPT.all ~ Weight_Group, family = binomial(link = "logit"), data = cleaned_weight_data)
weight_coefs <- summary(weight_model)$coefficients
weight_or <- exp(weight_coefs[, "Estimate"])
weight_ci_lower <- exp(weight_coefs[, "Estimate"] - 1.96 * weight_coefs[, "Std. Error"])
weight_ci_upper <- exp(weight_coefs[, "Estimate"] + 1.96 * weight_coefs[, "Std. Error"])
weight_p_values <- weight_coefs[, "Pr(>|z|)"]
weight_results_df <- data.frame(
  Group = rownames(weight_coefs),
  OddsRatio = weight_or,
  LowerCI = weight_ci_lower,
  UpperCI = weight_ci_upper,
  PValue = weight_p_values
)
print(weight_results_df)
data <- data[data$SEX %in% c("F", "M"), ]
data$Sex_Group <- as.numeric(data$SEX == "F")
data <- process_age_data(data)$Cleaned_Data
data <- process_weight_data(data)$Cleaned_Data
data$Age_Group <- cut(data$AGE, 
                      breaks = c(-Inf, 39, 54, 69, Inf), 
                      labels = c("0-39", "40-54", "55-69", ">70"))
data$Weight_Group <- cut(data$WT, 
                         breaks = c(-Inf, 63.5, 76, 92, Inf), 
                         labels = c("0-63.5", "63.5-76.0", "76.0-92.0", ">92.0"))
table(data$Age_Group)
table(data$Weight_Group)
multivar_model <- glm(TargetPT.all ~ Sex_Group + Age_Group + Weight_Group, family = binomial(link = "logit"), data = data)
summary(multivar_model)
coefs <- summary(multivar_model)$coefficients
odds_ratio <- exp(coefs[, "Estimate"])
ci_lower <- exp(coefs[, "Estimate"] - 1.96 * coefs[, "Std. Error"])
ci_upper <- exp(coefs[, "Estimate"] + 1.96 * coefs[, "Std. Error"])
multivar_results <- data.frame(
  Term = rownames(coefs),
  OddsRatio = odds_ratio,
  LowerCI = ci_lower,
  UpperCI = ci_upper,
  PValue = coefs[, "Pr(>|z|)"]
)
print(multivar_results)
