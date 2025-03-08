
library(openxlsx)
library(ggplot2)
library(forcats)
library(scales)

data1 <- read.xlsx("Data_Logistic_Regression_Vigibase.xlsx")
data2 <- read.xlsx("Data_Logistic_Regression_FAERS.xlsx")

data1$Group <- factor(data1$Group, levels = c(
  "Drug_VEGF",
  "Drug_VEGFR(ref)",
  "Age_Group>75",
  "Age_Group65-74",
  "Age_Group45-64",
  "Age_Group0-44(ref)",
  "Sex_Female",
  "Sex_Male(ref)"
))

data2$Group <- factor(data2$Group, levels = c(
  "Drug_VEGF",
  "Drug_VEGFR(ref)",
  "Age_Group>75",
  "Age_Group65-74",
  "Age_Group45-64",
  "Age_Group0-44(ref)",
  "Sex_Female",
  "Sex_Male(ref)"
))

data1 <- data1 %>%
  mutate(
    Color_ULR = ifelse(LowerCI_ULR <= 1 & UpperCI_ULR >= 1, "gray",
                       ifelse(OddsRatio_ULR > 1, "light_orange", "light_green")),
    Color_MLR = ifelse(LowerCI_MLR <= 1 & UpperCI_MLR >= 1, "gray",
                       ifelse(OddsRatio_MLR > 1, "dark_orange", "dark_green"))
  )

data2 <- data2 %>%
  mutate(
    Color_ULR = ifelse(LowerCI_ULR <= 1 & UpperCI_ULR >= 1, "gray",
                       ifelse(OddsRatio_ULR > 1, "light_orange", "light_green")),
    Color_MLR = ifelse(LowerCI_MLR <= 1 & UpperCI_MLR >= 1, "gray",
                       ifelse(OddsRatio_MLR > 1, "dark_orange", "dark_green"))
  )

a <- ggplot(data1, aes(y = Group)) +
  geom_point(aes(x = OddsRatio_ULR, color = Color_ULR), 
             shape = 16,  
             size = 4) + 
  geom_errorbarh(aes(xmin = LowerCI_ULR, xmax = UpperCI_ULR, color = Color_ULR), 
                 height = 0.4,  
                 size = 0.8) + 
  geom_point(aes(x = OddsRatio_MLR, color = Color_MLR), 
             shape = 15,  
             size = 4) + 
  geom_errorbarh(aes(xmin = LowerCI_MLR, xmax = UpperCI_MLR, color = Color_MLR), 
                 height = 0.4,  
                 size = 0.8) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
  scale_x_log10(breaks = c(0.0,0.5, 1,2, 3,4,5),labels = label_comma()) +
  scale_color_manual(values = c("gray" = "gray", "light_orange" = "#f9b384", "light_green" = "#74c69d",
                                "dark_orange" = "#f3722c", "dark_green" = "#40916c")) +
  labs(x = "Odds Ratio", y = "Group", title = "Univariate and multivariate logistic regression analysis of clinical factors \nassociated with TMA in patients using VEGFi inhibitors based on the Vigibase") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill =  "#E2E5D9", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 13),
    axis.title.x = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.title.y = element_blank()
  )

b <- ggplot(data2, aes(y = Group)) +
  geom_point(aes(x = OddsRatio_ULR, color = Color_ULR), 
             shape = 16,  
             size = 4) + 
  geom_errorbarh(aes(xmin = LowerCI_ULR, xmax = UpperCI_ULR, color = Color_ULR), 
                 height = 0.4,  
                 size = 0.8) + 
  geom_point(aes(x = OddsRatio_MLR, color = Color_MLR), 
             shape = 15,  
             size = 4) + 
  geom_errorbarh(aes(xmin = LowerCI_MLR, xmax = UpperCI_MLR, color = Color_MLR), 
                 height = 0.4,  
                 size = 0.8) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 0.5) +
  scale_x_log10(breaks = c(0.0,0.5, 1,2, 3,4,5),labels = label_comma()) +
  scale_color_manual(values = c("gray" = "gray", "light_orange" = "#f9b384", "light_green" = "#74c69d",
                                "dark_orange" = "#f3722c", "dark_green" = "#40916c")) +
  labs(x = "Odds Ratio", y = "Group", title = "Univariate and multivariate logistic regression analysis of clinical factors \nassociated with TMA in patients using VEGFi inhibitors based on the FAERS") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill =  "#E4E2DD", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.text.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 13),
    axis.title.x = element_blank(),
    axis.ticks.y = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_line(color = "black", linewidth = 0.5),
    axis.title.y = element_blank()
  )

library(patchwork)

b <- b + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

combined_plot <- a + b

print(combined_plot)
