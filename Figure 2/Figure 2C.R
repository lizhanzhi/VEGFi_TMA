load("Demo_VEGFi_TMA.Rdata")

data <- Demo_VEGFi_TMA

analyze_tto <- function(data, start_date_col, event_date_col) {
  convert_date_8digit <- function(date) {
    if (nchar(date) != 8) {
      return(NA)
    } else {
      tryCatch({
        as.Date(date, "%Y%m%d")
      }, error = function(e) NA)
    }
  }
  
  data[[start_date_col]] <- as.character(data[[start_date_col]])
  data[[event_date_col]] <- as.character(data[[event_date_col]])
  
  valid_rows_indices <- which(!is.na(data[[start_date_col]]) & !is.na(data[[event_date_col]]) &
                                nchar(data[[start_date_col]]) == 8 & nchar(data[[event_date_col]]) == 8)
  
  valid_data <- data[valid_rows_indices, ]
  
  valid_data[[start_date_col]] <- sapply(valid_data[[start_date_col]], convert_date_8digit)
  valid_data[[event_date_col]] <- sapply(valid_data[[event_date_col]], convert_date_8digit)
  
  valid_data$Days_to_Event <- as.numeric(valid_data[[event_date_col]] - valid_data[[start_date_col]])
  
  valid_data_positive <- valid_data[valid_data$Days_to_Event > 0, ]
  
  median_days_to_event <- median(valid_data_positive$Days_to_Event, na.rm = TRUE)
  min_days_to_event <- min(valid_data_positive$Days_to_Event, na.rm = TRUE)
  max_days_to_event <- max(valid_data_positive$Days_to_Event, na.rm = TRUE)
  mean_days_to_event <- mean(valid_data_positive$Days_to_Event, na.rm = TRUE)
  sd_days_to_event <- sd(valid_data_positive$Days_to_Event, na.rm = TRUE)
  missing_count <- nrow(data) - nrow(valid_data_positive)
  
  days_to_event_stats <- data.frame(
    Statistic = c("Mean (SD)", "Median [Min, Max]", "Missing"),
    Value = c(
      sprintf("%.1f (%.1f)", mean_days_to_event, sd_days_to_event),
      sprintf("%.1f [%.1f, %.1f]", median_days_to_event, min_days_to_event, max_days_to_event),
      sprintf("%d (%.1f%%)", missing_count, (missing_count / nrow(data)) * 100)
    )
  )
  
  list(
    Valid_Data = valid_data_positive,
    Statistics = days_to_event_stats
  )
}

result <- analyze_tto(data, "START_DT", "EVENT_DT")

data <- result$Valid_Data

library(survival)
library(survminer)

fit <- survfit(Surv(Days_to_Event, TargetPT.all) ~ vegf.all, data = data)

morandi_colors <- c("#E57373", "#6a4c93")

surv_plot <- ggsurvplot(
  fit,
  data = data,
  fun = "event",
  censor = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.2,
  risk.table.col = "strata",
  linetype = 1,
  ggtheme = theme_bw(),
  break.x.by = 200,
  legend.title = "Group",
  legend.labs = c("VEGF", "VEGFR"),
  pval = FALSE,
  risk.table.fontsize = 3,
  xlab = "Time to Event (days)",
  ylab = "Cumulative Percentage(%)",
  xlim = c(0, 1400),
  palette = morandi_colors,
  legend = "none"
)

surv_plot$plot <- surv_plot$plot + 
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) paste0(x * 100)
  ) +
  annotate("text", x = 800, y = 0.75, label = "TTO Median", color = "black", size = 3.5, vjust = -0.5) +
  annotate("text", x = 800, y = 0.65, label = paste("VEGF Median:", 201, "[70 - 415]"), color = "black", size = 3.5, vjust = -0.5) +
  annotate("text", x = 800, y = 0.55, label = paste("VEGFR Median:", 35, "[21 - 70]"), color = "black", size = 3.5, vjust = -0.5) +
  annotate("text", x = 800, y = 0.45, label = paste("Wilcoxon test: P <", 2.2e-16), color = "black", size = 3.5, vjust = -0.5) +
  theme(
    legend.justification = c(0.3, 0.5),
    legend.position = c(0.3, 0.5),
    legend.box.just = "right",
    legend.box.spacing = unit(1, "lines"),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.7, "cm"),
    legend.text = element_text(color = "black", size = 10),
    text = element_text(color = "black")
  )

print(surv_plot)

library(dplyr)

data <- result$Valid_Data
data <- data %>%
  filter(!is.na(SEX) & trimws(SEX) != "")

fit <- survfit(Surv(Days_to_Event, TargetPT.all) ~ SEX, data = data)
morandi_colors <- c("#f9b384", "#74c69d")

surv_plot <- ggsurvplot(
  fit,
  data = data,
  fun = "event",
  censor = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.2,
  risk.table.col = "strata",
  linetype = 1,
  ggtheme = theme_bw(),
  break.x.by = 200,
  legend.title = "Group",
  legend.labs = c("Female", "Male"),
  pval = FALSE,
  risk.table.fontsize = 3,
  xlab = "Time to Event (days)",
  ylab = "Cumulative Percentage(%)",
  xlim = c(0, 1400),
  palette = morandi_colors
)

surv_plot$plot <- surv_plot$plot + 
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) paste0(x * 100)
  ) +
  annotate("text", x = 800, y = 0.75, label = "TTO Median", color = "black", size = 3.5, vjust = -0.5) +
  annotate("text", x = 800, y = 0.65, label = paste("Female Median:", 25, "[17 - 159]"), color = "black", size = 3.5, vjust = -0.5) +
  annotate("text", x = 800, y = 0.55, label = paste("VEGFR Median:", 29, "[55 - 152]"), color = "black", size = 3.5, vjust = -0.5) +
  annotate("text", x = 800, y = 0.45, label = paste("Wilcoxon test: P =", 0.252), color = "black", size = 3.5, vjust = -0.5) +
  theme(
    legend.justification = c(0.3, 0.5),
    legend.position = c(0.3, 0.5),
    legend.box.just = "right",
    legend.box.spacing = unit(1, "lines"),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.7, "cm"),
    legend.text = element_text(color = "black", size = 10),
    text = element_text(color = "black")
  )

print(surv_plot)

library(openxlsx)

data <- read.xlsx("Data_TTO_Vi_TMA.xlsx")

library(dplyr)
library(tidyr)

cleaned_data <- data %>%
  filter(TimeToOnsetMin != "") %>%
  mutate(TimeToOnsetMin = strsplit(as.character(TimeToOnsetMin), ",\\s*")) %>%
  unnest(TimeToOnsetMin, keep = TRUE) %>%
  mutate(TimeToOnsetMin = as.numeric(TimeToOnsetMin)) %>%
  filter(TimeToOnsetMin > 0) %>%
  group_by(UMCReportId) %>%
  mutate(MinTimeToOnset = min(TimeToOnsetMin)) %>%
  filter(TimeToOnsetMin == MinTimeToOnset) %>%
  distinct(UMCReportId, .keep_all = TRUE)

data <- cleaned_data

fit <- survfit(Surv(MinTimeToOnset, any_pt) ~ Drug_Type, data = data)

morandi_colors <- c("#E57373", "#6a4c93")

surv_plot <- ggsurvplot(
  fit,
  data = data,
  fun = "event",
  censor = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.2,
  risk.table.col = "strata",
  linetype = 1,
  ggtheme = theme_bw(),
  break.x.by = 200,
  legend.title = "Group",
  legend.labs = c("VEGF", "VEGFR"),
  pval = FALSE,
  risk.table.fontsize = 3,
  xlab = "Time to Event (days)",
  ylab = "Cumulative Percentage(%)",
  xlim = c(0, 1400),
  palette = morandi_colors,
  legend = "none"
)

surv_plot$plot <- surv_plot$plot + 
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) paste0(x * 100)
  ) +
  annotate("text", x = 800, y = 0.75, label = "TTO Median", color = "black", size = 3.5, vjust = -0.5) +
  annotate("text", x = 800, y = 0.65, label = paste("VEGF Median:", 183, "[55 - 284]"), color = "black", size = 3.5, vjust = -0.5) +
  annotate("text", x = 800, y = 0.55, label = paste("VEGFR Median:", 27, "[20 - 75]"), color = "black", size = 3.5, vjust = -0.5) +
  annotate("text", x = 800, y = 0.45, label = paste("Wilcoxon test: P <", 0.001), color = "black", size = 3.5, vjust = -0.5) +
  theme(
    legend.justification = c(0.3, 0.5),
    legend.position = c(0.3, 0.5),
    legend.box.just = "right",
    legend.box.spacing = unit(1, "lines"),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.7, "cm"),
    legend.text = element_text(color = "black", size = 10),
    text = element_text(color = "black")
  )

print(surv_plot)

data <- data %>%
  filter(!is.na(Gender) & trimws(Gender) != "")

fit <- survfit(Surv(MinTimeToOnset, any_pt) ~ Gender, data = data)

morandi_colors <- c("#f9b384", "#74c69d")

surv_plot <- ggsurvplot(
  fit,
  data = data,
  fun = "event",
  censor = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.2,
  risk.table.col = "strata",
  linetype = 1,
  ggtheme = theme_bw(),
  break.x.by = 200,
  legend.title = "Group",
  legend.labs = c("Female", "Male"),
  pval = FALSE,
  risk.table.fontsize = 3,
  xlab = "Time to Event (days)",
  ylab = "Cumulative Percentage(%)",
  xlim = c(0, 1400),
  palette = morandi_colors
)

surv_plot$plot <- surv_plot$plot + 
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = function(x) paste0(x * 100)
  ) +
  annotate("text", x = 800, y = 0.75, label = "TTO Median", color = "black", size = 3.5, vjust = -0.5) +
  annotate("text", x = 800, y = 0.65, label = paste("Female Median:", 75, "[25 - 258]"), color = "black", size = 3.5, vjust = -0.5) +
  annotate("text", x = 800, y = 0.55, label = paste("Male Median:", 68, "[41 - 258]"), color = "black", size = 3.5, vjust = -0.5) +
  annotate("text", x = 800, y = 0.45, label = paste("Wilcoxon test: P =", 0.720), color = "black", size = 3.5, vjust = -0.5) +
  theme(
    legend.justification = c(0.3, 0.5),
    legend.position = c(0.3, 0.5),
    legend.box.just = "right",
    legend.box.spacing = unit(1, "lines"),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.7, "cm"),
    legend.text = element_text(color = "black", size = 10),
    text = element_text(color = "black")
  )

print(surv_plot)
