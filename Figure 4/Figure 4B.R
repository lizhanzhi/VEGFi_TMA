library(openxlsx)
library(ggplot2)
library(ggrepel)

load("TCGA_TMA.Rdata")
ror_cancer <- read.xlsx("TCGA_TMA_ROR_Vi.xlsx")
target <- read.xlsx("TCGA_traget_pathway_TMA.xlsx")

selected_paths <- target$Data
medians_selected <- all_medians_df[selected_paths,]

medians_long <- t(medians_selected)
colnames(medians_long) <- target$Pathway
medians_long_df <- as.data.frame(medians_long)
medians_long_df$cancer <- rownames(medians_long_df)

ror_values <- ror_cancer[, c("TCGA_TYPE", "ROR")]
ror_values$TCGA_TYPE <- toupper(ror_values$TCGA_TYPE)

medians_with_ror_only <- merge(medians_long_df, ror_values, by.x="cancer", by.y="TCGA_TYPE")

pathway_of_interest <- target$Pathway[6]
plot_data <- medians_with_ror_only[, c('ROR', pathway_of_interest, 'cancer')]
colnames(plot_data) <- c('ROR', 'val', 'cancer')

want <- target$Pathway[6]
head(plot_data)
tmp <- plot_data

cor_test <- cor.test(tmp$val, tmp$ROR, method = "spearman")

txt <- paste0("Spearman's rho = ", round(cor_test$estimate, 3), 
              ", P-value = ", format(cor_test$p.value, digits = 3))

cols <- c('#D0C1DE','#B3CC47','#9A80B5','#F2ACAF','#791A1C','#F6D0D6','#502A81','#DB2D86','#EEB167',
          '#F8E1C5','#10A296','#0B77AA','#6C779C','#D17927','#90CBA4','#20A3DC','#F0E33D','#0C477B',
          '#BAA233','#D7EEF8','#9BD6F1','#C8CBD9','#CC99C0','#A74D96','#0C8D42','#EF9128',"darkgreen",
          "chocolate4","blueviolet","#223D6C","#D20A13","#088247")

ggplot(tmp, aes(ROR, val)) + 
  geom_point(aes(colour = cancer), size = 2.5) + 
  geom_smooth(span = 1.5, method = glm, color = "black", se = FALSE, linewidth = 0.8) + 
  geom_text_repel(aes(label = cancer), size = 3, point.padding = 0.5, min.segment.length = 0, max.time = 1, 
                  max.overlaps = 25, max.iter = 1e5, box.padding = 0.5) + 
  xlab("Platinum-drugs related inflammatory myopathy AEs ROR") + 
  ylab('ssGSEA scores') +
  labs(title = want) +
  theme_bw() +
  scale_color_manual(values = cols) +
  theme(
    plot.title = element_text(colour ='black', size = 10, lineheight = 1, hjust =0.5),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 10, color = "black"),
    panel.background = element_blank(),
    panel.grid = element_blank()) +
  annotate("text", x = 2, y = 0.46, hjust = 0, fontface = 3, size = 4, label = txt)

library(gridExtra)

plots <- list()

for (i in 1:6) {
  pathway_of_interest <- target$Pathway[i]
  plot_data <- medians_with_ror_only[, c('ROR', pathway_of_interest, 'cancer')]
  colnames(plot_data) <- c('ROR', 'val', 'cancer')
  
  cor_test <- cor.test(plot_data$val, plot_data$ROR, method = "spearman")
  txt <- paste0("Spearman's rho = ", round(cor_test$estimate, 3), 
                ", P-value = ", format(cor_test$p.value, digits = 3))
  
  p <- ggplot(plot_data, aes(x = ROR, y = val, colour = cancer)) +
    geom_point(size = 2.5) +
    geom_smooth(method = glm, se = FALSE, color = "black", linewidth = 0.8) +
    geom_text_repel(aes(label = cancer), size = 3, point.padding = 0.5, 
                    min.segment.length = 0, max.iter = 10000, box.padding = 0.5) +
    labs(title = pathway_of_interest, subtitle = txt) +
    xlab("VEGFi drugs related TMA AEs ROR") + 
    ylab('ssGSEA scores') +
    theme_bw() +
    scale_color_manual(values = cols) +
    theme(plot.title = element_text(size = 12, lineheight = 1.2, hjust = 0.5),
          axis.text = element_text(size = 10), 
          axis.title = element_text(size = 10),
          legend.position = if (i == 6) "none" else "none")
  
  plots[[i]] <- p
}

grid.arrange(grobs = plots, nrow = 2, ncol = 3)

ggplot(tmp, aes(ROR, val)) + 
  geom_point(aes(colour = cancer), size = 2.5) + 
  geom_smooth(span = 1.5, method = glm, color = "black", se = FALSE, linewidth = 0.8) + 
  geom_text_repel(aes(label = cancer), size = 3, point.padding = 0.5,
                  min.segment.length = 0, max.time = 1, 
                  max.overlaps = 25, max.iter = 1e5, box.padding = 0.5) +
  xlab("Platinum-drugs related inflammatory myopathy AEs ROR") + 
  ylab('ssGSEA scores') +
  labs(title = want) +
  theme_bw() +
  scale_color_manual(values = cols) +
  theme(
    plot.title = element_text(size = 10, lineheight = 1, hjust = 0.5),
    axis.text = element_text(size = 10, color = "black"), 
    axis.title = element_text(size = 10, color = "black"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.box = "horizontal"
  ) +
  annotate("text", x = 2, y = 0.46, hjust = 0, fontface = 3, size = 4, label = txt)

