# Run Inter-ViSTA ----

# switch to interactome tab and disable file inputs

observeEvent(input$run, {
  
  # disable run button to prevent multiple clicks
  shinyjs::disable("run")
  
  # disable inputs
  shinyjs::disable("input_sidebar")
  
  # switch to viewing interactome tab
  updateTabsetPanel(session, "navbar",
                    selected = "interactome")
  
})

# run Inter-ViSTA functions

computed_values <- eventReactive(input$run, {
  
  # NEW DATA ----
  # set of functions to run if input files are new data
  
  if (ncol(user_input$nodes) <= 4) {
    
    # Clean data ----
    
    withProgress(message = "Calculating parameters", value = 0.7, {
      
      cleaned_data <- cleanData(user_input$nodes, user_input$edges, 
                                user_input$timepoints, user_input$spec_threshold)
      
      # returns: nodes, edges, weight_threshold, timepoints, interval,
      # loc_prov, abund_prov, bait_ids, taxids
      
      computed_values <- cleaned_data
      
    })
    
    
    # Protein complexes ----
    withProgress(message = "Detecting protein complexes", value = 0.9, {
      
      complex_data <- calculateComplexes(computed_values$nodes)
      
      # returns: nodes, complexes
      computed_values$nodes <- complex_data$nodes
      computed_values$complexes <- complex_data$complexes
      
    })
    
    # Quantitative tools ----
    if (computed_values$abund_prov) {
      withProgress(message = "Performing quantitative analysis", value = 0.8, {
        
        # Cluster proteins ----
        cluster_data <- clusterProteins(computed_values$nodes, computed_values$edges,
                                        computed_values$bait_ids,
                                        computed_values$timepoints)
        # returns: edges, abundances
        computed_values$edges <- cluster_data$edges
        computed_values$abundances <- cluster_data$abundances
        
        # Normalize by proteome abundance ----
        # If there is a proteome abundance file provided, normalize by proteome abundance
        if (!is.null(user_input$proteome_abundance)) {
          normalized_data <- normalizeByProteomeAbund(user_input$proteome_abundance,
                                                      computed_values$nodes,
                                                      computed_values$edges,
                                                      computed_values$timepoints)
          # returns: edges, normalized abundances
          computed_values$edges <- normalized_data$edges
          computed_values$normalized_abundances <- normalized_data$normalized_abundances
          computed_values$norm_prot <- TRUE
        } else {
          computed_values$norm_prot <- FALSE
        }
        
        # Size nodes by abundance (single bait) ----
        # If there is only one bait and there are fewer than 500 proteins,
        # size interactome nodes by abundance
        if ( length(computed_values$bait_ids) == 1 && dim(computed_values$nodes)[1] < 500 ) {
          sized_data <- sizeByAbund(computed_values$nodes, computed_values$edges,
                                    computed_values$bait_ids,
                                    computed_values$timepoints)
          # returns: nodes
          computed_values$nodes <- sized_data
        }
      })
    }
    
    # Get Uniprot data ----
    withProgress(message = "Annotating UniProt data", detail="This may take a while...",
                 value = 0.6, {
                   
                   uniprot_data <- getUniprotData(user_input$background, computed_values$nodes,
                                                    computed_values$taxids, computed_values$loc_prov)
                   
                   # returns: nodes, organisms, uniprot, GOList
                   computed_values$nodes <- uniprot_data$nodes
                   computed_values$organisms <- uniprot_data$organisms
                   computed_values$localizations <- uniprot_data$localizations
                   computed_values$uniprot <- uniprot_data$uniprot
                   computed_values$GOList <- uniprot_data$GOList
                   
                 })
    
    # GO enrichment ----
    withProgress(message = "Performing gene ontology (GO) enrichment", 
                 detail="This may take a while...",
                 value = 0.8, {
                   
                   GO_data <- performGOEnrichment(computed_values$GOList, computed_values$nodes)
                   
                   # returns: nodes, GOTable, topGOterms
                   computed_values$nodes <- GO_data$nodes
                   computed_values$GOTable <- GO_data$GOTable
                   computed_values$topGOterms <- GO_data$topGOterms
                   
                 })
    
    # Calculate activity spells ----
    withProgress(message = "Calculating activity spells", value = 0.7, {
      
      activity_data <- calculateActivitySpells(computed_values$nodes, computed_values$edges,
                                               computed_values$timepoints, computed_values$interval,
                                               computed_values$weight_threshold,
                                               computed_values$bait_ids)
      # returns: nodes, edges, bait_ids
      computed_values$nodes <- activity_data$nodes
      computed_values$edges <- activity_data$edges
      computed_values$bait_ids <- activity_data$bait_ids
      
    })
    
    # Get STRING interactions ----
    withProgress(message = "Integrating STRING interactions", detail="This may take a while...",
                 value = 0.6, {
                   
                   string_data <- getSTRINGData(computed_values$nodes, computed_values$edges,
                                                computed_values$taxids, computed_values$timepoints)
                   
                   # returns: edges
                   computed_values$edges <- string_data
                   
                 })
    
    
  } else {
    
    # PREVIOUSLY CALCULATED DATA ----
    
    # user has input files that have already been calculated by Inter-ViSTA
    
    # Clean data ----
    
    withProgress(message = "Calculating parameters", value = 0.7, {
      
      cleaned_data <- cleanDataPREV(user_input$nodes, user_input$edges)
      
      # returns: nodes, edges, num_edges, timepoints, interval,
      # abund_prov, norm_prot, bait_ids
      computed_values <- cleaned_data
      
      # gather data that hasn't changed (uniprot list)
      computed_values$uniprot <- user_input$background
      computed_values$GOTable <- user_input$GOTable
      
    })
    
    # Protein complexes ----
    withProgress(message = "Detecting protein complexes", value = 0.9, {
      
      complex_data <- calculateComplexesPREV(computed_values$nodes)
      
      # returns: complexes
      computed_values$complexes <- complex_data
      
    })
    
    # Quantitative tools ----
    if (computed_values$abund_prov) {
      withProgress(message = "Performing quantitative analysis", value = 0.8, {
        
        # Cluster proteins ----
        cluster_data <- clusterProteinsPREV(computed_values$edges, computed_values$timepoints)
        # returns: abundances
        computed_values$abundances <- cluster_data
        
        # Normalize by proteome abundance ----
        # If there is a proteome abundance file provided, normalize by proteome abundance
        if ( computed_values$norm_prot ) {
          normalized_data <- normalizeByProteomeAbundPREV(computed_values$edges, 
                                                          computed_values$timepoints)
          # returns: normalized abundances
          computed_values$normalized_abundances <- normalized_data
        }
        
        # Size nodes by abundance (single bait) ----
        # If there is only one bait and there are fewer than 500 proteins,
        # size interactome nodes by abundance
        if ( length(computed_values$bait_ids) == 1 && dim(computed_values$nodes)[1] < 500 ) {
          sized_data <- sizeByAbundPREV(computed_values$nodes, computed_values$edges,
                                        computed_values$bait_ids,
                                        computed_values$timepoints)
          # returns: nodes
          computed_values$nodes <- sized_data
        }
      })
    }
    
    # GO enrichment ----
    if (is.null(computed_values$GOTable)) {
      
      withProgress(message = "Performing gene ontology (GO) enrichment", 
                   detail="This may take a while...",
                   value = 0.8, {
                     
                     GO_data <- performGOEnrichmentPREV(computed_values$uniprot, 
                                                        computed_values$nodes)
                     
                     # returns: GOTable
                     computed_values$GOTable <- GO_data
                     
                   })
      
    }
    
    # Calculate activity spells ----
    withProgress(message = "Calculating activity spells", value = 0.7, {
      
      activity_data <- calculateActivitySpellsPREV(computed_values$nodes, computed_values$edges,
                                               computed_values$timepoints, computed_values$interval,
                                               computed_values$bait_ids)
      # returns: nodes, edges, bait_ids
      computed_values$nodes <- activity_data$nodes
      computed_values$edges <- activity_data$edges
      computed_values$bait_ids <- activity_data$bait_ids
      
    })
  }
  
  # Build interactome ----
  if (dim(computed_values$nodes)[1] < 500) {
    withProgress(message = "Building interactome network", value = 0.9, {
      
      network_data <- buildNetwork(computed_values$nodes, computed_values$edges,
                                   computed_values$num_edges, computed_values$abund_prov,
                                   computed_values$norm_prot, computed_values$timepoints, 
                                   computed_values$interval, computed_values$bait_ids)
      
      # returns: nodes, edges, net3
      computed_values$nodes <- network_data$nodes
      computed_values$edges <- network_data$edges
      computed_values$net3.dyn <- network_data$net3.dyn
      
    })
  }
  
  # Tidy data for abundance plots ----
  if (computed_values$abund_prov) {
    withProgress(message = "Tidying data for plots", value = 0.7, {
      
      plot_data <- cleanAbundanceData(computed_values$edges, computed_values$timepoints)
      
      # returns: cleaned abundances to plot
      computed_values$plotabundances <- plot_data
      
    })
  }
  
  computed_values
  
})