

data <- as.data.frame(fread("data/input_withmissing_rat_skm_proteomics_022625.txt"))


cleanedData <- data[data$pct_missing==0,]
rownamesMeasurements <- cleanedData$Protein.Names

measurements <- cleanedData[,5:25]


rownames(measurements) <- rownamesMeasurements



sampleInfo <- as.data.frame(matrix(ncol=2, nrow = ncol(measurements)))
colnames(sampleInfo) <- c("SampleID","SampleGroup")
sampleInfo$SampleID <- colnames(measurements)
sampleInfo$SampleGroup <- c(rep("ZSF1_OB",8),rep("WKY",6),rep("ZSF1_LN",7))



ProteinInfo <- cleanedData[,c(2:4,26:27)]
rownames(ProteinInfo) <- ProteinInfo$Protein.Names


dge <- DGEList(counts = measurements, samples = sampleInfo, genes = ProteinInfo)

#Design a model for comparison of samples
design <- model.matrix(~ 0 + SampleGroup, data = dge$samples)
colnames(design) <- c("WKY","ZSF1_LN","ZSF1_OB")

#Design a contrast for pairwise comparison of samples
contrasts <- makeContrasts(
  levels = colnames(design),
  ZSF_O_Vs_ZSF_L = (ZSF1_OB - ZSF1_LN),
  WKY_Vs_ZSF_O = (WKY - ZSF1_OB),
  WKY_Vs_ZSF_L = (WKY - ZSF1_LN)
)
ContrastsUsed <- c("ZSF_O_Vs_ZSF_L","WKY_Vs_ZSF_O","WKY_Vs_ZSF_L")

fit <- lmFit(dge$counts, design)%>%
  contrasts.fit(contrasts)



fit <- eBayes(fit, trend=TRUE)

for (contrast in ContrastsUsed) {
  
  allDEProteins <- topTable(fit,coef = contrast,number = Inf, adjust.method = "fdr")
  allDEProteins <- allDEProteins %>%
    dplyr::mutate(isSignificant = case_when(
      adj.P.Val < 0.05 & abs(logFC) > 1 ~ TRUE,
      TRUE ~ FALSE
    ))
  allDEProteins$Protein <- rownames(allDEProteins)
  write.csv(allDEProteins,file = paste0("data/",contrast,".csv"),row.names = TRUE)
  
  allDEProteinsSig <- allDEProteins[allDEProteins$adj.P.Val < 0.05 & abs(allDEProteins$logFC) > 1,]
  allDEProteinsSig <- allDEProteinsSig[order(allDEProteinsSig$logFC, decreasing = TRUE),]
  
  write.csv(allDEProteinsSig,file = paste0("data/",contrast,"_Significant.csv"),row.names = TRUE)
  
  genesOfInterest <- cleanedData$Genes[cleanedData$Protein.Names %in% rownames(allDEProteinsSig)]
  
  write(paste0(genesOfInterest,"\n"), file = paste0("data/",contrast,"_SigGenes.txt"))
    
    
  top10 <- head(allDEProteinsSig,10)
      
  #Volcano plot to visualize the results
  volcanoPlot <- allDEProteins %>%
        ggplot(aes(x = logFC,
                   y = -log10(adj.P.Val),
                   colour = isSignificant)) +
        geom_vline(xintercept = 1, linetype = "dotted") +
        geom_vline(xintercept = -1, linetype = "dotted") +
        geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
        geom_point(size = 1, alpha = 0.5) +
        scale_colour_manual(values = c("grey", "red")) +
        #geom_text(data = top10 ,aes(label = GeneSymbol),check_overlap = TRUE) +
        geom_text_repel(data = top10 ,aes(label = Protein),size = 5) +
        ggtitle(paste0(contrast,", Thresholds: FDR<0.05, logFC>1")) +
        guides(color = guide_legend(title = "IsSignificant")) +
        theme(plot.title = element_text(hjust = 0.5,size = 15),
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 15),
              axis.title = element_text(size = 15),
              axis.text = element_text(size = 10)
        )
      
      
  ggsave(paste0("data/",contrast,"_VolcanoPlot.jpeg"), volcanoPlot, width = 35, height = 9,dpi = 600)
  
  
  MAPlot <- allDEProteins %>%
    ggplot(aes(x = AveExpr,
               y = logFC,
               colour = isSignificant)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_point(size = 1, alpha = 0.5) +
    scale_colour_manual(values = c("grey", "red")) +
    geom_text(data = top10 ,aes(label = Protein),check_overlap = TRUE) +
    ggtitle(paste0(contrast)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(paste0("data/",contrast,"_MAPlot.jpeg"), MAPlot, width = 35, height = 9,dpi = 600)
  
  
  if(contrast == "ZSF_O_Vs_ZSF_L"){
    
    requiredColumns <- grep("^zsf1_ob",colnames(measurements))
    requiredColumns <- c(requiredColumns,grep("^ZSF1_Lean",colnames(measurements)))
    requiredMeasurements <- measurements[,requiredColumns]
    
    sampleNames <- colnames(requiredMeasurements)
    
    zfs1_ob_SamplesCount <- length(grep(paste0("^", "zsf1_ob"), sampleNames))
    zsf1_lean_SamplesCount <- length(grep(paste0("^", "ZSF1_Lean"), sampleNames))
    
    sampleGroups <- factor(c(rep("ZSF1_OB",zfs1_ob_SamplesCount),
                             rep("ZSF1_LN",zsf1_lean_SamplesCount)
    ))
    
    groupColors <- c("ZSF1_OB" = "blue",
                     "ZSF1_LN" = "red")
    
    sigProteins <- allDEProteinsSig$Protein
    
    sigProteinIntensityData <- requiredMeasurements[rownames(requiredMeasurements) %in% sigProteins,]
    
    scaledSigExpData <- t(scale(t(sigProteinIntensityData)))
    
    topAnnotation <- HeatmapAnnotation(
      Group = sampleGroups,
      col = list(Group = groupColors)
    )
    
    jpeg(paste0("data/",contrast,"_TopHeatMap.jpeg"), width = 1200, height = 1000,quality = 100)
    
    topHeatmap <- Heatmap(scaledSigExpData,
                          name = "Expression",
                          row_names_side = "left",
                          column_names_side = "top",
                          clustering_distance_rows = "euclidean",
                          cluster_columns = FALSE,
                          top_annotation = topAnnotation,
                          column_title = paste0("Heatmap of significant proteins (FDR<0.05) for ",contrast," contrast"),
                          column_title_gp = gpar(fontsize = 20, fontface = "bold")
    )
    
    print(topHeatmap)
    
    dev.off()
    
    
  }else if (contrast == "WKY_Vs_ZSF_O"){
    
    requiredColumns <- grep("^wky",colnames(measurements))
    requiredColumns <- c(requiredColumns,grep("^zsf1_ob",colnames(measurements)))
    requiredMeasurements <- measurements[,requiredColumns]
    
    sampleNames <- colnames(requiredMeasurements)
    
    wky_SamplesCount <- length(grep(paste0("^", "wky"), sampleNames))
    zfs1_ob_SamplesCount <- length(grep(paste0("^", "zsf1_ob"), sampleNames))
    
    
    sampleGroups <- factor(c(rep("WKY",wky_SamplesCount),
                             rep("ZSF1_OB",zfs1_ob_SamplesCount)
    ))
    
    groupColors <- c("WKY" = "blue",
                     "ZSF1_OB" = "red")
    
    sigProteins <- allDEProteinsSig$Protein
    
    sigProteinIntensityData <- requiredMeasurements[rownames(requiredMeasurements) %in% sigProteins,]
    
    scaledSigExpData <- t(scale(t(sigProteinIntensityData)))
    
    topAnnotation <- HeatmapAnnotation(
      Group = sampleGroups,
      col = list(Group = groupColors)
    )
    
    jpeg(paste0("data/",contrast,"_TopHeatMap.jpeg"), width = 1200, height = 1000,quality = 100)
    
    topHeatmap <- Heatmap(scaledSigExpData,
                          name = "Expression",
                          row_names_side = "left",
                          column_names_side = "top",
                          clustering_distance_rows = "euclidean",
                          cluster_columns = FALSE,
                          top_annotation = topAnnotation,
                          column_title = paste0("Heatmap of significant proteins (FDR<0.05) for ",contrast," contrast"),
                          column_title_gp = gpar(fontsize = 20, fontface = "bold")
    )
    
    print(topHeatmap)
    
    dev.off()
    
  }else if (contrast == "WKY_Vs_ZSF_L"){
    
    requiredColumns <- grep("^wky",colnames(measurements))
    requiredColumns <- c(requiredColumns,grep("^ZSF1_Lean",colnames(measurements)))
    requiredMeasurements <- measurements[,requiredColumns]
    
    sampleNames <- colnames(requiredMeasurements)
    
    wky_SamplesCount <- length(grep(paste0("^", "wky"), sampleNames))
    zfs1_ob_SamplesCount <- length(grep(paste0("^", "ZSF1_Lean"), sampleNames))
    
    
    sampleGroups <- factor(c(rep("WKY",wky_SamplesCount),
                             rep("ZSF1_LEAN",zfs1_ob_SamplesCount)
    ))
    
    groupColors <- c("WKY" = "blue",
                     "ZSF1_LEAN" = "red")
    
    sigProteins <- allDEProteinsSig$Protein
    
    sigProteinIntensityData <- requiredMeasurements[rownames(requiredMeasurements) %in% sigProteins,]
    
    scaledSigExpData <- t(scale(t(sigProteinIntensityData)))
    
    topAnnotation <- HeatmapAnnotation(
      Group = sampleGroups,
      col = list(Group = groupColors)
    )
    
    jpeg(paste0("data/",contrast,"_TopHeatMap.jpeg"), width = 1200, height = 1000,quality = 100)
    
    topHeatmap <- Heatmap(scaledSigExpData,
                          name = "Expression",
                          row_names_side = "left",
                          column_names_side = "top",
                          clustering_distance_rows = "euclidean",
                          cluster_columns = FALSE,
                          top_annotation = topAnnotation,
                          column_title = paste0("Heatmap of significant proteins (FDR<0.05) for ",contrast," contrast"),
                          column_title_gp = gpar(fontsize = 20, fontface = "bold")
    )
    
    print(topHeatmap)
    
    dev.off()
    
  }
  
  DoEnrichRAnnotation(contrast = contrast)
  
}
  




