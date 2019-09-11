# UI elements for Downloads tab

output$downloads <- renderUI({
  
  tagList(
    
    verbatimTextOutput("downloadNetworkInfo"),
    downloadButton("downloadNodes", "Node Attributes"),
    downloadButton("downloadEdges", "Edge Attributes"),
    
    br(),
    br(),
    
    if (dim(computed_values()$nodes)[1] < 500) {
      
      tagList (
        verbatimTextOutput("downloadReportInfo"),
        downloadButton("downloadReport", "Report"),
        
        br(),
        br()
        
      )
    },
    
    verbatimTextOutput("downloadUniprotInfo"),
    downloadButton("downloadUniprot", "UniProt Annotated Background Gene List")
    
  )
  
})

# download network information ----
output$downloadNetworkInfo <- renderPrint(
  cat("Download text files (tab-separated) with inferred and provided
      protein and edge attributes. This can be immediately input into Cytoscape
      for further analysis. Attributes include:
      protein - localizations, GO terms, organism name, abundance, etc. and
      interaction (provided and STRING) - onset, terminus, confidence, etc.")
)

# output nodes ----
output$downloadNodes <- downloadHandler(
  filename= function() { 
    "interactome_node_attributes.txt"}
  ,
  content = function(file) {
    saveNodeAttributes(computed_values()$nodes, 
                       computed_values()$bait_ids,
                       computed_values()$timepoints,
                       file)
  }
)

# output edges ----
output$downloadEdges <- downloadHandler(
  filename= function() { 
    "interactome_edge_attributes.txt"}
  ,
  content = function(file) {
    saveEdgeAttributes(computed_values()$edges,
                       computed_values()$timepoints,
                       computed_values()$abund_prov,
                       computed_values()$norm_prot,
                       file)
  }
)


# download report information ----
output$downloadReportInfo <- renderPrint(
  cat("Download network analysis graphs containing: spells of protein activity
      (in each localization and associated with each GO term), plot of
      time of interaction formation, and histogram of interaction durations.")
)

# output report ----
output$downloadReport <- downloadHandler(
  filename= function() { 
    "interactome_analytical_graphs.pdf"}
  ,
  content = function(file) {
    pdf(file)
    analyticalTools(computed_values()$net3.dyn,
                    computed_values()$nodes,
                    computed_values()$edges,
                    computed_values()$timepoints,
                    computed_values()$localizations,
                    computed_values()$topGOterms)
    dev.off()
  }
) 

# download background gene list information ----
output$downloadUniprotInfo <- renderPrint(
  cat("Download annotated background gene list that contains UniProt-derived information,
      such as gene ontology terms for background genes. Using this list will 
      significantly reduce Inter-ViSTA computation time in the future.")
)


# download table
output$downloadUniprot <- downloadHandler(
  filename = function() {
    "annotated_background_gene_list.csv"}
  ,
  content = function(file) {
    write.table(computed_values()$background, file = file, sep = ",", 
                quote = F, row.names = F)
  }
)