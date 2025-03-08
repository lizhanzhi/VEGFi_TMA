
res <- read.csv('Raw_Camp_input.csv', check.names = FALSE)

res_BA <- res[res$drug == "Bevacizumab", ]
res_VE <- res[res$drug == "VEGFR", ]

res_BA_unique <- res_BA[!duplicated(res_BA$SYMBOL), ]
res_VE_unique <- res_VE[!duplicated(res_VE$SYMBOL), ]

library(org.Mm.eg.db)
library(clusterProfiler)

convert_symbols_to_entrez <- function(df) {
  converted <- tryCatch({
    bitr(df$SYMBOL, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
  }, error = function(e) {
    message("Error converting symbols: ", e$message)
    return(NULL)
  })
  
  converted_unique <- converted[!duplicated(converted$SYMBOL), ]
  
  merge(df, converted_unique, by = "SYMBOL")
}

res_BA_converted <- convert_symbols_to_entrez(res_BA_unique)
res_VE_converted <- convert_symbols_to_entrez(res_VE_unique)

geneList_BA <- sort(setNames(res_BA_converted$logFC, res_BA_converted$SYMBOL), decreasing = TRUE)

hallmarks <- read.gmt("msigdb.v2023.2.Mm.symbols.gmt")

gsea_BA_results <- GSEA(geneList = geneList_BA, 
                        TERM2GENE = hallmarks, 
                        minGSSize=10, 
                        maxGSSize=1000,  
                        pAdjustMethod="BH", 
                        pvalueCutoff=0.05)

geneList_VE <- sort(setNames(res_VE_converted$logFC, res_VE_converted$SYMBOL), decreasing = TRUE)

gsea_VE_results <- GSEA(geneList = geneList_VE, 
                        TERM2GENE = hallmarks, 
                        minGSSize=10, 
                        maxGSSize=1000,  
                        pAdjustMethod="BH", 
                        pvalueCutoff=0.05)

library(ggplot2)
library(ggnewscale)

dotplot(gsea_VE_results, showCategory = 20, split = ".sign") +
  facet_grid(~.sign)

dotplot(gsea_VE_results, showCategory = 12, split = ".sign") +
  facet_grid(~.sign) +
  scale_y_discrete(labels = function(x) stringi::stri_sub(x, 1, 10))

if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}
cnetplot(gsea_VE_results, showCategory = 4, foldChange = geneList_VE, colorEdge = TRUE)

library(GseaVis)

gseaNb(object = gsea_BA_results,
       geneSetID = 'REACTOME_SIGNALING_BY_VEGF',
       newGsea = FALSE,
       subPlot = 3,
       addPval = TRUE,
       pvalX = 0.75,
       pvalY = 0.9,
       pCol = 'black',
       pHjust = 0)

gseaNb(object = gsea_VE_results,
       geneSetID = 'REACTOME_COMPLEMENT_CASCADE',
       newGsea = FALSE,
       subPlot = 3,
       addGene = c("C3", "C4b", "C1qa", "C1qb", "Serping1"),
       addPval = TRUE,
       pvalX = 0.75,
       pvalY = 0.9,
       pCol = 'black',
       pHjust = 0)

gseaNb(object = gsea_VE_results,
       geneSetID = 'REACTOME_RESPONSE_TO_ELEVATED_PLATELET_CYTOSOLIC_CA2',
       newGsea = FALSE,
       subPlot = 3,
       addGene = c("Pf4", "Fgg", "Fga", "Vegfd"),
       addPval = TRUE,
       pvalX = 0.75,
       pvalY = 0.9,
       pCol = 'black',
       pHjust = 0,
       rmSegment = TRUE)