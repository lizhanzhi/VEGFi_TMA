library(ggplot2)
library(dplyr)
library(scales)
library(readxl)
library(openxlsx)
library(forcats)

data <- read.xlsx("Data_FA_Vi_TMA.xlsx")

data <- na.omit(data)

data$Drug_Database <- interaction(data$Drug, data$Database, sep = "\n")

light_gray <- "#d9d9d9"
dark_gray <- "#404040"
light_red <- "#ffcccc"
dark_red <- "#990000"

gray_colors <- colorRampPalette(c(dark_gray, light_gray))(5)
red_colors <- colorRampPalette(c(light_red, dark_red))(75)

min_value <- min(data$ROR[data$ROR > 0], na.rm = TRUE)
max_value <- max(data$ROR, na.rm = TRUE)
break_points <- c(0, 1, max_value)
values_rescaled <- scales::rescale(break_points)

ggplot(data, aes(x = fct_reorder(Drug_Database, Drug, .desc = TRUE), y = PT, size = a, color = ROR)) +
  geom_point(alpha = 0.6) +
  scale_size(range = c(1, 20), guide = "legend") +
  scale_color_gradientn(colors = c(gray_colors, red_colors), values = c(0, 0.5, 1), limits = c(0, max_value)) +
  labs(title = "Risk bubble chart of AEs to TMA caused by VEGFi drugs",
       x = "DrugClass",
       y = "PT",
       size = "Report Count",
       color = "ROR") +
  theme_minimal(base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(color = "black", hjust = 0.5, size = 14),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", colour = "white"),
    panel.grid.minor = element_blank()
  ) +
  guides(size = guide_legend(override.aes = list(color = "black")))+
  scale_x_discrete(limits = c("All\nFAERS", "All\nVigibase", "VEGF\nFAERS", "VEGF\nVigibase", "VEGFR\nFAERS", "VEGFR\nVigibase",
                              "Bevacizumab\nFAERS", "Bevacizumab\nVigibase", "Ramucirumab\nFAERS", "Ramucirumab\nVigibase", "Sunitinib\nFAERS", "Sunitinib\nVigibase", "Aflibercept\nFAERS","Aflibercept\nVigibase"))
