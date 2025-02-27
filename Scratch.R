

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
      
      
  ggsave(paste0(contrast,"_VolcanoPlot.jpeg"), volcanoPlot, width = 35, height = 9,dpi = 600)
  
}




