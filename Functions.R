DoEnrichRAnnotation <- function(contrast){
  
    
    genes <- read.table(paste0("data/",contrast,"_SigGenes.txt"))
    enrichmentResults <- enrichr(genes = genes$V1 , databases = "WikiPathways_2024_Mouse")
    res <- enrichmentResults$WikiPathways_2024_Mouse
    res <- res[res$`P.value` < 0.01,]
    
    #Visualize the bar plot of the pathways
    barPlot <- ggplot(res, aes(x = reorder(`Term`, `Combined.Score`), y = `Combined.Score`)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      coord_flip() +
      labs(x = "Pathway", y = "Combined Score", title = paste0("Pathway Enrichment for ",contrast," contrast")) +
      theme_minimal()
    
    
    #Visualize the dot plot of the pathways
    dotPlot <- ggplot(res, aes(x=`Overlap`, y=reorder(`Term`, -log10(`Adjusted.P.value`)), 
                               size=-log10(`Adjusted.P.value`), color=`Odds.Ratio`)) +
      geom_point() +
      scale_color_gradient(low="blue", high="red") +
      labs(x="Gene Overlap", y="Pathway", title= paste0("Pathway Enrichment Dot Plot for ",contrast," contrast")) +
      theme_minimal()
    
    ggsave(paste0("data/",contrast,"_BarPlot.jpeg"), barPlot,dpi = 600, width = 20)
    
    ggsave(paste0("data/",contrast,"_DotPlot.jpeg"), dotPlot,dpi = 600, width = 20)
    
    
  
  
}
