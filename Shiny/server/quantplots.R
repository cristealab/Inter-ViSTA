# UI elements in quantplots main panel that need to be updated dynamically ----

output$quantplotsMain <- renderUI({
  
  v <- computed_values()
  
  if (v$abund_prov) {
    
    tagList(
      
      # interactive plots of abundance of single (or neighboring) proteins ----
      h4("Plot protein and interactor quantitative values"),
      plotOutput("abundPlot"),
      downloadButton("dAbund", "Download Quantitative Plot"),
      
      br(),
      br(),
      
      # plots of proteins clustered by abundance profiles ----
      verbatimTextOutput("clusterNumber"),
      br(),
      h4("Prey Proteins per Bait Clustered by Temporal Quantitative Profiles"),
      uiOutput("clusterPlots"),
      downloadButton("dCluster", "Download Cluster Plots"),
      
      br(),
      br(),
      
      if (v$norm_prot) {
        
        tagList(
          
          # heatmaps comparing normalized to not-normalized protein abundances ----
          h4("Heatmaps of Scaled Quantitative Profiles, Prior to and After Normalization to Proteome"),
          uiOutput("heatmapPlots"),
          downloadButton("dHeatmaps", "Download Heatmaps")
          
        )
        
      }
      
    )
  } else {
    
    verbatimTextOutput("No abundances provided")
    
  }
  
})


# server logic to build quantitative plots ----

# on button click, find neighbors for given proteins ----
confNeighbors <- eventReactive(input$abundGo, {
  
  v <- computed_values()
  calcNeighbors(input$confidence, input$checkLocalizations, 
                input$checkSpecies, input$checkGOTerms,
                v$nodes, v$edges, v$timepoints, v$localizations)
})

observeEvent(input$abundGo, {

  v <- computed_values()

  # render abundances plot ----
  output$abundPlot <- renderPlot({

    abundPlot(input$genesym, input$neighbors, confNeighbors(),
              v$nodes, v$plotabundances, v$timepoints)

  })
  
  addTooltip(session, id = "abundPlot", title = "Interactors are determined using parameters selected
                       in 'Interactome' tab, such as STRING confidence threshold, localizations, 
            GO terms, etc. Modify parameters to change displayed interactors", 
             placement = "bottom", trigger = "hover",
             options = NULL)

  # download abundances plot ----
  output$dAbund <- downloadHandler(

    filename= function() {
      "abundance_plot.pdf"}
    ,
    content = function(file) {
      pdf(file)
      abundPlot(input$genesym, input$neighbors, confNeighbors(),
                v$nodes, v$plotabundances, v$timepoints)
      dev.off()
    }
  )

  # output cluster number for selected protein ----
  output$clusterNumber <- renderPrint({

    printClusterNumber(input$genesym,
                       v$nodes, v$edges, v$bait_ids)

  })

})

## FOR EACH BAIT, PLOT CLUSTERING RESULTS ----

# function to make list of cluster plots for each bait ----
get_clusterplot_output_list <- function(NUM_BAITS) {
  
  # Insert plot output objects into the list
  
  plot_output_list <- lapply(1:NUM_BAITS, function(i) {
    
    plotname <- paste("plot", i, sep="")
    plot_output_object <- plotOutput(plotname)
    plot_output_object <- renderPlot({
    
      v <- computed_values()
      clusterPlot(v$bait_ids[i], v$nodes, v$abundances, v$timepoints)
        
    })
    
  })
  
  do.call(tagList, plot_output_list) # needed to display properly.
  
  return(plot_output_list)
  
}

# make UI for list of cluster plots ----
output$clusterPlots <- renderUI({ 
  
  NUM_BAITS <- length(computed_values()$bait_ids)
  get_clusterplot_output_list(NUM_BAITS) 
  
})


# download cluster plots ----
output$dCluster <- downloadHandler(

  filename= function() {
    "cluster_plots.pdf"}
  ,
  content = function(file) {
    pdf(file)
    for (bait in computed_values()$bait_ids) {
      clusterPlot(bait, computed_values()$nodes, computed_values()$abundances,
                  computed_values()$timepoints)
    }
    dev.off()
  }
)


## HEATMAPS IF NORMALIZATION TO THE PROTEOME IS REQUIRED

# function to make list of heatmap plots for each bait ----
get_heatmapplot_output_list <- function(NUM_BAITS) {
  
  # Insert plot output objects into the list
  
  plot_output_list <- lapply(1:NUM_BAITS, function(i) {
    
    plotname <- paste("plot", i, sep="")
    plot_output_object <- plotOutput(plotname)
    plot_output_object <- renderPlot({
      
      v <- computed_values()
      plotHeatmaps(v$bait_ids[i], v$nodes, v$normalized_abundances, v$timepoints)
      
    })
    
  })
  
  do.call(tagList, plot_output_list) # needed to display properly.
  
  return(plot_output_list)
  
}

# make UI for list of heatmap plots ----
output$heatmapPlots <- renderUI({ 
  
  NUM_BAITS <- length(computed_values()$bait_ids)
  get_heatmapplot_output_list(NUM_BAITS) 
  
})


# download heatmap plots ----
output$dHeatmap <- downloadHandler(
  
  filename= function() {
    "heatmap_plots.pdf"}
  ,
  content = function(file) {
    pdf(file)
    for (bait in computed_values()$bait_ids) {
      plotHeatmaps(bait, computed_values()$nodes, computed_values()$normalized_abundances,
                  computed_values()$timepoints)
    }
    dev.off()
  }
)



