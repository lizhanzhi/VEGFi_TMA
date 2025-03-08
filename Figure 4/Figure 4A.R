library(openxlsx)
library(ggplot2)
library(dplyr)
ror_fa <- read.xlsx("TCGA_TMA_ROR_FA.xlsx")
ror_vi <- read.xlsx("TCGA_TMA_ROR_Vi.xlsx")

ror_fa$Data <- 'FA'
ror_vi$Data <- 'VI'

combined_data <- rbind(ror_fa, ror_vi)

combined_data$TCGA_DATA <- interaction(combined_data$TCGA_TYPE, combined_data$Data)

base_colors <- c('#D0C1DE','#B3CC47','#9A80B5',
                 '#F2ACAF','#791A1C','#F6D0D6',
                 '#502A81','#DB2D86','#EEB167',
                 '#F8E1C5','#10A296','#0B77AA',
                 '#6C779C','#D17927','#90CBA4',
                 '#20A3DC','#F0E33D','#0C477B',
                 '#BAA233','#D7EEF8','#9BD6F1',
                 '#C8CBD9','#CC99C0','#A74D96',
                 '#0C8D42','#EF9128',"darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247")

if (length(base_colors) < length(unique(combined_data$TCGA_TYPE))) {
  base_colors <- rep(base_colors, length.out = length(unique(combined_data$TCGA_TYPE)))
}

light_colors <- sapply(base_colors[1:length(unique(combined_data$TCGA_TYPE))], adjustcolor, 1)
dark_colors <- sapply(base_colors[1:length(unique(combined_data$TCGA_TYPE))], adjustcolor, 1)

color_mapping <- setNames(c(light_colors, dark_colors), unique(combined_data$TCGA_DATA))

print(head(color_mapping))

vi_order <- combined_data %>%
  filter(Data == "VI") %>%
  arrange(desc(ROR)) %>%
  pull(TCGA_TYPE)

combined_data$TCGA_TYPE <- factor(combined_data$TCGA_TYPE, levels = rev(vi_order))

dodge_width <- 0.9

ggplot(data = combined_data, aes(x = TCGA_TYPE, y = ROR, fill = TCGA_DATA)) +
  geom_bar(stat = "identity", position = position_dodge(width = dodge_width)) +
  geom_text(aes(label = sprintf("%.2f", ROR), group = TCGA_DATA),
            position = position_dodge(width = dodge_width), vjust = -0.3, color = "black", size = 3.5, check_overlap = TRUE) +
  scale_fill_manual(values = color_mapping) +
  labs(title = "Pan-cancer analysis of FAERS and Vigibase",
       x = "TCGA Type",
       y = "ROR Value",
       fill = "Data Set") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.ticks.y = element_line(color = "black"),
    legend.position = "none"
  ) +
  theme(plot.title = element_text(hjust = 0.5, size = 14))
