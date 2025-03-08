library(data.table)
library(limma)
library(openxlsx)
gene_read_counts <- read.table("Kid_genes.tpm.txt", header = TRUE, sep = "\t", row.names = NULL)

if(any(duplicated(gene_read_counts[[1]]))) {
  gene_read_counts[[1]] <- make.unique(as.character(gene_read_counts[[1]]))
}

rownames(gene_read_counts) <- gene_read_counts[[1]]
gene_read_counts[[1]] <- NULL

bev_treatment_samples <- grep("^AB|^CB", colnames(gene_read_counts), value = TRUE)
bev_control_samples <- grep("^AP|^CP", colnames(gene_read_counts), value = TRUE)

vegfr_treatment_samples <- grep("^AS|^CS", colnames(gene_read_counts), value = TRUE)
vegfr_control_samples <- grep("^AP|^CP", colnames(gene_read_counts), value = TRUE)

selected_data_bev <- gene_read_counts[, c(bev_treatment_samples, bev_control_samples)]
group_bev <- factor(c(rep("Treatment_Bev", length(bev_treatment_samples)), rep("Control_Bev", length(bev_control_samples))))
design_bev <- model.matrix(~ group_bev)

selected_data_vegfr <- gene_read_counts[, c(vegfr_treatment_samples, vegfr_control_samples)]
group_vegfr <- factor(c(rep("Treatment_VEGFR", length(vegfr_treatment_samples)), rep("Control_VEGFR", length(vegfr_control_samples))))
design_vegfr <- model.matrix(~ group_vegfr)

log_selected_data_bev <- log2(selected_data_bev + 1)
log_selected_data_vegfr <- log2(selected_data_vegfr + 1)

gene_variances_bev <- apply(log_selected_data_bev, 1, var)
gene_variances_vegfr <- apply(log_selected_data_vegfr, 1, var)

threshold <- 0.1
log_selected_data_bev <- log_selected_data_bev[gene_variances_bev > threshold, ]
log_selected_data_vegfr <- log_selected_data_vegfr[gene_variances_vegfr > threshold, ]

fit_bev <- lmFit(log_selected_data_bev, design_bev)
fit_bev <- eBayes(fit_bev)

allDiff_bev <- topTable(fit_bev, coef=2, adjust="BH", number = Inf)
allDiff_bev$drug <- "Bevacizumab"
allDiff_bev$SYMBOL <- rownames(allDiff_bev)

fit_vegfr <- lmFit(log_selected_data_vegfr, design_vegfr)
fit_vegfr <- eBayes(fit_vegfr)

allDiff_vegfr <- topTable(fit_vegfr, coef=2, adjust="BH", number = Inf)
allDiff_vegfr$drug <- "VEGFR"
allDiff_vegfr$SYMBOL <- rownames(allDiff_vegfr)

library(dplyr)

filtered_bev <- allDiff_bev %>%
  filter(P.Value < 0.05)

filtered_vegfr <- allDiff_vegfr %>%
  filter(P.Value < 0.05)

genes_in_filtered_vegfr <- filtered_vegfr$SYMBOL
filtered_log_data_vegfr <- log_selected_data_vegfr[rownames(log_selected_data_vegfr) %in% genes_in_filtered_vegfr, ]

genes_in_filtered_bev <- filtered_bev$SYMBOL
filtered_log_data_bev <- log_selected_data_bev[rownames(log_selected_data_bev) %in% genes_in_filtered_bev, ]

pathway_genes <- c("Crk", "Pdpk1", "Cyfip2", "Nrp1", "Cdh5", "Flt1", "Actg1", 
                   "Hsp90aa1", "Nrp2", "Ncf4", "Ctnnb1", "Flt4", "Dock1", "Kdr",
                   "C3", "Cfd", "Igkv14-126", "C4b", "Ighv1-76", "C7", "Iglc2", "Ighg3", 
                   "Ighv11-2", "Igkv1-110", "Igkv1-135", "Igkv1-117", "Cfp", "Ighv1-64", 
                   "Ighv4-1", "C1qa", "Igkv15-103", "Ighv1-74", "Vtn", "Ighv3-6", 
                   "Igkv8-21", "C6", "Ighv1-18", "Igkv14-111", "Serping1", "Igkv13-85", 
                   "C1qb", "C5ar1", "Iglc1", "Ighv1-55", "Ighv1-69",
                   "Orm1", "Trf", "Cfd", "Fgg", "Fga", "Timp1", "Itih4", "F13a1", 
                   "Pf4", "Apoh", "Abcc4", "Vegfd", "Orm3", "Clec3b", "Cdc37l1", 
                   "Srgn", "Flna", "Mmrn1")

filtered_pathway_genes <- filtered_bev[filtered_bev$SYMBOL %in% pathway_genes, ]
print(head(filtered_pathway_genes))

filtered_pathway_genes_vegfr <- filtered_vegfr[filtered_vegfr$SYMBOL %in% pathway_genes, ]
print(head(filtered_pathway_genes_vegfr))

complement_genes <- c("C3", "Apoc1", "Kcnip2", "Timp1", "Cebpb", "Kynu", "Cp", "Pim1", 
                      "Tnfaip3", "Plat", "Gzmk", "Spock2", "Ccl5", "Maff", "Cblb", "Gzmb", 
                      "C1qa", "Ang", "Pcsk9", "Msrb1", "Plaur", "Csrp1", "Notch4", "Gngt2", 
                      "Pla2g7", "Serping1", "S100a13", "Kcnip3", "Mmp14", "Lck", "C1qc", 
                      "Lrp1", "Timp2", "Mmp15", "Prdm4", "C1s1", "Was")

platelet_genes <- c("Orm1", "Trf", "Cfd", "Fgg", "Fga", "Timp1", "Itih4", "F13a1", 
                    "Pf4", "Apoh", "Abcc4", "Vegfd", "Orm3", "Clec3b", "Cdc37l1", 
                    "Srgn", "Flna", "Mmrn1")

filtered_complement_genes <- filtered_vegfr[filtered_vegfr$SYMBOL %in% complement_genes, ]
filtered_platelet_genes <- filtered_vegfr[filtered_vegfr$SYMBOL %in% platelet_genes, ]

library(ComplexHeatmap)

for (i in 1:nrow(filtered_log_data_bev)) {
  filtered_log_data_bev[i, ] <- scale(log2(unlist(filtered_log_data_bev[i, ]) + 1))
}

mat <- as.matrix(filtered_log_data_bev)

samples <- ifelse(grepl('^AP|^CP', colnames(mat)), 'Control', 'Experiment')

heat <- Heatmap(mat, 
                col = colorRampPalette(c("navy", "white", "firebrick"))(100),
                heatmap_legend_param = list(grid_height = unit(10, 'mm')),
                show_row_names = FALSE,
                top_annotation = HeatmapAnnotation(Group = samples, 
                                                   simple_anno_size = unit(2, 'mm'), 
                                                   col = list(Group = c('Control' = '#00DAE0', 'Experiment' = '#FF9289')),
                                                   show_annotation_name = TRUE), 
                column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 6))

pathway_genes <- c("Orm1", "Fgg")

heat + rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% pathway_genes), 
                                      labels = pathway_genes, labels_gp = gpar(fontsize = 10)))

draw(heat)

filtered_log_data_vegfr$AS1_Kid.1 <- NULL

for (i in 1:nrow(filtered_log_data_vegfr)) {
  filtered_log_data_vegfr[i, ] <- scale(log2(unlist(filtered_log_data_vegfr[i, ]) + 1))
}

mat <- as.matrix(filtered_log_data_vegfr)

samples <- ifelse(grepl('^AP|^CP', colnames(mat)), 'Control', 'Experiment')

heat <- Heatmap(mat, 
                col = colorRampPalette(c("navy", "white", "firebrick"))(100),
                heatmap_legend_param = list(grid_height = unit(10, 'mm')),
                show_row_names = FALSE,
                top_annotation = HeatmapAnnotation(Group = samples, 
                                                   simple_anno_size = unit(2, 'mm'), 
                                                   col = list(Group = c('Control' = '#00DAE0', 'Experiment' = '#FF9289')),
                                                   show_annotation_name = TRUE), 
                column_names_gp = gpar(fontsize = 10), row_names_gp = gpar(fontsize = 6))

pathway_genes <- c("Crk", "Pdpk1", "Cyfip2", "Nrp1", "Cdh5", "Flt1", "Actg1", 
                   "Hsp90aa1", "Nrp2", "Ncf4", "Ctnnb1", "Flt4", "Dock1", "Kdr",
                   "C3", "Cfd", "Igkv14-126", "C4b", "Ighv1-76", "C7", "Iglc2", "Ighg3", 
                   "Ighv11-2", "Igkv1-110", "Igkv1-135", "Igkv1-117", "Cfp", "Ighv1-64", 
                   "Ighv4-1", "C1qa", "Igkv15-103", "Ighv1-74", "Vtn", "Ighv3-6", 
                   "Igkv8-21", "C6", "Ighv1-18", "Igkv14-111", "Serping1", "Igkv13-85", 
                   "C1qb", "C5ar1", "Iglc1", "Ighv1-55", "Ighv1-69",
                   "Orm1", "Trf", "Cfd", "Fgg", "Fga", "Timp1", "Itih4", "F13a1", 
                   "Pf4", "Apoh", "Abcc4", "Vegfd", "Orm3", "Clec3b", "Cdc37l1", 
                   "Srgn", "Flna", "Mmrn1")

heat <- heat + rowAnnotation(link = anno_mark(at = which(rownames(mat) %in% pathway_genes), 
                                              labels = pathway_genes, labels_gp = gpar(fontsize = 10)))

draw(heat)
