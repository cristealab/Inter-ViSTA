
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
                       choices=LOCS, selected=LOCS),
    
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
                                       width="100%", height="800px")
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
                 v$nodes, v$edges, v$timepoints, v$interval, v$bait_ids)
  }
  
})

# output network animation
output$netPlot <- ndtv:::renderNdtvAnimationWidget({
  
  v <- computed_values()
  
  withProgress(message="Updating network parameters", detail="This may take a while...", 
               value=0.7, {
                 plotCompNet(inputNet(), lab=input$labels, v$bait_ids, v$abund_prov)
               })
  
})
