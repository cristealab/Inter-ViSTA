
# UI elements in interactome sidebar that need to be updated dynamically ----

output$interactomeToggles <- renderUI({
  
  v <- computed_values()
  
  NUM_BAITS <- length(v$bait_ids)
  LOCS <- v$localizations
  TOP_GO_TERMS <- v$topGOterms
  ORGANISMS <- v$organisms
  
  tagList(
    
    # Radio buttons for number of shared proteins ----
    radioButtons("shared",
                 label="Show proteins shared between at least 'n' baits:",
                 choices=1:NUM_BAITS, selected=1, inline=T),
    
    # Checkboxes for which localizations to include ----
    checkboxGroupInput("checkLocalizations", label="Show proteins localized to:",
                       choices=c(LOCS, "All"), selected="All"),
    
    # Checkboxes for which GO terms to include ----
    if (!is.null(TOP_GO_TERMS)) {
      checkboxGroupInput("checkGOTerms", label="Show proteins associated with:",
                         choices=c(TOP_GO_TERMS, "All"),
                         selected="All")
    },
    
    # Checkboxes for which species to include ----
    checkboxGroupInput("checkSpecies", label="Show species:",
                       choices=ORGANISMS, selected=ORGANISMS)
    
  )
  
})

# UI elements in main panel (interactome) that need to be updated ----

output$interactomeMain <- renderUI({
  
  NUM_NODES <- dim(computed_values()$nodes)[1]
  
  # if there are fewer than 500 nodes, build network
  if (NUM_NODES < 500) {
    
    tagList(
      # stop boostrap css from messing up the tooltip in the widget
      tags$style(HTML(".tooltip {opacity: 1}")), 
      
      # animation widget for interactome
      ndtv:::ndtvAnimationWidgetOutput("netPlot", 
                                       width="100%", height="800px"),
      
      # legend for localization node color
      plotOutput("localizationLegend")
      
    )
    
  } else {
    
    h3("Your input has too many proteins (>= 500) to visualize in interactome network. However, other app functionalities are still available!")
    
  }
  
})


## server elements for network ----

# calculate inside a reactive to prevent multiple computations
inputNet <- eventReactive(input$renderNetwork, {
  
  v <- computed_values()
  
  if (dim(v$nodes)[1] < 500) {
    deactivateNE(input$shared, input$confidence, input$checkLocalizations, 
                 input$checkSpecies, input$checkGOTerms, v$net3.dyn,
                 v$nodes, v$edges, v$timepoints, v$bait_ids)
  }
  
})

# output network animation
output$netPlot <- ndtv:::renderNdtvAnimationWidget({
  
  v <- computed_values()
  
  withProgress(message="Updating network parameters", detail="This may take a while...", 
               value=0.7, {
                 plotCompNet(inputNet(), lab=input$labels, disp_tps = input$display_timepoints, v$bait_ids, v$abund_prov, v$named_timepoints)
               })
  
})

# output a legend for the interactome network
makeLegend <- eventReactive(input$renderNetwork, {
  
  v <- computed_values()
  
  locs <- as.list(get.vertex.attribute(v$net3.dyn, "localization_string_output"))
  colors <- as.list(get.vertex.attribute(v$net3.dyn, "loc_color"))
  
  node_colors <- colors
  names(node_colors) <- as.list(get.vertex.attribute(v$net3.dyn, "localization_string"))
  
  names(locs) <- colors
  legend_mapping <- locs[!duplicated(locs)]
  
  keep <- c()
  
  if (!grepl('All', input$checkLocalizations)){
    for (loc in input$checkLocalizations){
      for (l in legend_mapping){
        if (grepl(loc, l)){
          if (grepl(';', l)){
            keep <- append(keep, 'Multiple localizations')
            
          } else {
            keep <- append(keep, l)
          }
        }
      }
    }
  } else {
    
    keep <- names(node_colors)
  }
  
  node_colors <- unlist(node_colors[unique(keep)])
  
  labels <- names(node_colors)
  labels[labels == ''] <- 'Unspecified'
  
  names(node_colors) <- NULL
  
  plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
  legend("topleft", legend = labels, col = node_colors, pch=16, pt.cex=3, cex=1.5, bty='n', border = '#000000', ncol = 1)
  mtext("Uniprot annotated localization", cex=2)
  
})



# render lecalization legend ----
output$localizationLegend <- renderPlot({
  
  makeLegend()
  
})
