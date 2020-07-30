# Store User Input ----

# store user input data as reactive values (variables) for further computation,
# then store computed values as variables for later use in server

user_input <- reactiveValues()

# User Chooses Sample Datasets ----
# when sample files are loaded

observeEvent(input$loadsamples, {
  sample_ds <- input$samples
  bg_name <- ""
  
  if (sample_ds != "none") {
    user_input$nodes <- read.table(paste("Sample Datasets/", sample_ds, "_nodes.txt", 
                                       sep=""), sep="\t", stringsAsFactors = F, header = T)
    user_input$edges <- read.table(paste("Sample Datasets/", sample_ds, "_edges.txt", 
                                       sep=""), sep="\t", comment.char="?", stringsAsFactors = F,
                                   header = T)
    user_input$background <- read.csv(paste("Sample Datasets/", sample_ds, "_uniprot.csv", 
                                       sep=""), stringsAsFactors = F)
    user_input$GOTable <- read.table(paste("Sample Datasets/", sample_ds, "_GOTable.txt", 
                                           sep=""), sep="\t", stringsAsFactors = F, header = T)
    
    # change the background gene list and specificity threshold to be used based on 
    # which sample dataset is selected 
    # (THIS HAS TO BE HARD CODED MANUALLY BASED ON OUR SELECTED THRESHOLDS)
    
    if (sample_ds == "apex") {
      user_input$spec_threshold <- 0.63
      user_input$timepoints <- c(1, 3, 10, 30)
      bg_name <- "hum_hela"
      
    } else if (sample_ds == "cyclins") {
      user_input$spec_threshold <- 0.05
      user_input$timepoints <- c(1, 2, 3, 4)
      bg_name <- "hum_hela"
      
    } else if (sample_ds == "ul13_ul37") {
      user_input$spec_threshold <- 0.9
      user_input$timepoints <- c(24, 48, 72, 96, 120)
      bg_name <- "hum_fib"
    } else if (sample_ds == "us9") {
      user_input$spec_threshold <- 1
      user_input$timepoints <- c(3, 8, 12, 18, 24)
      bg_name <- "rat_neuron"
    }
    

    # update timepoints textbox
    updateTextInput(session, "timepoints", value = paste(user_input$timepoints, collapse=", "))
    
    # update specificity filter slider
    updateSliderInput(session, "specificity", value = user_input$spec_threshold)
    
    # update sample background file dropdown
    updateSelectInput(session, "bgsamples", selected = bg_name)
  }
})

# when sample background files are loaded
observeEvent(input$loadbgsamples, {
  
  if (input$bgsamples != "none"){
  
    user_input$background <- read.csv(paste("Sample Datasets/Background Gene Lists/",
                                            input$bgsamples, ".csv", sep = ""),
                                      header = T, sep = ",", stringsAsFactors = F)
  }
})

# User Uploads Own Datasets ----

# when user uploads their own files (this rewrites whatever sample files were loaded)
observeEvent(input$nodes, {
  user_input$nodes <- read.csv(input$nodes$datapath, header = TRUE, sep = ",",
                               stringsAsFactors = F)
  user_input$new_data <- TRUE
  
})

observeEvent(input$edges, {
  user_input$edges <- read.csv(input$edges$datapath, header = TRUE, sep = ",",
                               stringsAsFactors = F)
})

observeEvent(input$background, {
  user_input$background <- read.csv(input$background$datapath, header = TRUE, sep = ",",
                                    stringsAsFactors = F)
})

observeEvent(input$specificity, {
  user_input$spec_threshold <- input$specificity
})

observeEvent(input$timepoints, {
  user_input$timepoints <- unlist(strsplit(input$timepoints, ", "))
})

observeEvent(input$proteome_abundance, {
  user_input$proteome_abundance <- read.csv(input$proteome_abundance$datapath, header = TRUE,
                                            sep = ",", stringsAsFactors = F)
})


# User Uploads Previously Calculated Datasets ----

# when user uploads previously calculated files (this rewrites whatever sample files were loaded)
observeEvent(input$prev_nodes, {
  user_input$nodes <- read.table(input$prev_nodes$datapath, header = TRUE, sep = "\t",
                               stringsAsFactors = F)
  user_input$new_data <- FALSE
  
})

observeEvent(input$prev_edges, {
  user_input$edges <- read.table(input$prev_edges$datapath, header = TRUE, sep = "\t",
                               stringsAsFactors = F, comment.char = "?")
})

observeEvent(input$prev_uniprot, {
  user_input$background <- read.csv(input$prev_uniprot$datapath, header = TRUE, sep = ",",
                                    stringsAsFactors = F)
})

observeEvent(input$prev_GOTable, {
  user_input$GOTable <- read.csv(input$prev_GOTable$datapath, header = TRUE, sep = "\t",
                                    stringsAsFactors = F)
})

# Print Loaded Data ----
# print uploaded data files for user to verify that data has been parsed correctly

output$nodescontents <- renderDataTable(user_input$nodes, options = list(pageLength = 5, 
                                                                         scrollX = TRUE))
addTooltip(session, "nodescontents", "File must contain, in the following order: 
           accession, gene name, taxid, (optional) localization",
           placement = "top", trigger = "hover",
           options = NULL)

output$edgescontents <- renderDataTable(user_input$edges, options = list(pageLength = 5,
                                                                         scrollX = TRUE))
addTooltip(session, "edgescontents", "File must contain, in the following order: 
           bait accession, prey accession, confidence_condition1, (optional) quant_condition1,
           confidence_condition2, (optional) quant_condition2, etc... where 
           'condition1', 'condition2', etc. are numeric",
           placement = "top", trigger = "hover",
           options = NULL)

output$bgcontents <- renderDataTable(user_input$background, options = list(pageLength = 5,
                                                                                 scrollX = TRUE))
addTooltip(session, "bgcontents", "If annotated file provided, file must contain:
           UniProt accession, gene name, GOID, localization, organism", 
           placement = "top", trigger = "hover", options = NULL)

output$proteomecontents <- renderDataTable(head(user_input$proteome_abundance),
                                           options = list(pageLength = 5, scrollX = TRUE))
addTooltip(session, "proteomecontents", "File must contain, in the following order: 
           UniProt accession, gene name, abund_condition1, abund_condition2, etc... where
           'condition1', 'condition2', etc. match exactly with the conditions provided in 
           edges file above", 
           placement = "top", trigger = "hover", options = NULL)

output$GOTablecontents <- renderDataTable(head(user_input$GOTable),
                                           options = list(pageLength = 5, scrollX = TRUE))

addTooltip(session, "GOTablecontents", "Gene ontology analysis file previously generated by Inter-ViSTA",
           placement = "top", trigger = "hover", options = NULL)

