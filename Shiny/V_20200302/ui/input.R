# Tab: Upload data  ----
tabPanel(title = "Upload Data", value = "input", fluid = TRUE,
         
         sidebarLayout(
           
           # Sidebar: Input parameters ----
           
           sidebarPanel(id = "input_sidebar",
             
             # Example Datasets Input ----
             
             h4("Load Sample Dataset"),
             
             selectInput("samples", "Choose Sample Dataset",
                         c("None" = "none",
                           "APEX (Lobingier et al. 2017)" = "apex",
                           "Cyclins (Pagliuca et al. 2011)" = "cyclins",
                           #"HCMV pUL13 and pUL37" = "ul13_ul37",
                           "PRV pUS9 (Kramer et al. 2012)" = "us9")
             ),
             
             actionButton("loadsamples", "Load Sample Dataset"),
             
             # Horizontal line ----
             tags$hr(),
             
             # Separate user input from example dataset input
             h4("Upload Your Own Data"),
             
             # Nodes Input: Select nodes file ----
             fileInput("nodes", "Upload Nodes CSV File",
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values",
                                  ".csv")),
             
             # Edges Input: Select edges file ----
             fileInput("edges", "Upload Edges CSV File",
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values",
                                  ".csv")),
             
             # Timepoint Input: textbox for timepoints / conditions ----
             textInput("timepoints", 'Enter Timepoints / Conditions separated by commas (","):',
                       value = "Ex: 24, 48, 72, 96, 120"),
             
             # Specificity filter ----
             sliderInput("specificity", label="Select interaction specificity threshold:",
                         min=0, max=1, value=0.7),
             bsTooltip("specificity", title = "Interactions with confidence scores below this 
                       threshold at every timepoint / condition will be discarded", 
                       placement = "bottom", trigger = "hover",
                       options = NULL),
             
             # Background gene list ----
             
             # use sample list
             selectInput("bgsamples", "Choose Sample Background Genes List", 
                         list("None" = "none",
                              "Human" = c("Fibroblast Tissue" = "hum_fib",
                                          "HeLa" = "hum_hela"),
                              "Rat" = c("Neuron (PC12)" = "rat_neuron"))
             ),
             
             actionButton("loadbgsamples", "Upload Sample Background Genes List"),
             
             br(),
             br(),
             
             # OR user provides
             fileInput("background", "OR Upload Background Genes CSV File",
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
             bsPopover(id = "background", title = "Provide either: a list of UniProt accessions 
                       OR a previously generated annotated file", 
                       placement = "bottom", trigger = "hover",
                       options = NULL),
             
             # (Optional) Proteome abundance file ----
             
             fileInput("proteome_abundance", "(Optional) Upload Proteome Abundance CSV File",
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
             
             bsPopover(id = "proteome_abundance", title = "Provide a file to normalize provided
                       interaction bait abundances to proteome abundances 
                       at each timepoint / condition", 
                       placement = "bottom", trigger = "hover",
                       options = NULL),
             
             # Horizontal line ----
             tags$hr(),
             
             # Previous Datasets Input ----
             
             h4("Load Datasets Previously Calculated by Inter-ViSTA"),
             
             # Previous Nodes Input: Select nodes file ----
             fileInput("prev_nodes", "Upload Previously Calculated Nodes File",
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values")),
             
             # Previous Edges Input: Select edges file ----
             fileInput("prev_edges", "Upload Previously Calculated Edges File",
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".txt")),
             
             # Previous UniProt Input: Select edges file ----
             fileInput("prev_uniprot", "Upload Previously Calculated UniProt Annotations File",
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
             
             # Previous GO enrichment table input: select GO Entrichment file ----
             fileInput("prev_GOTable", "Upload Previously Calculated Gene Ontology Enrichment Analysis",
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv"))
           ),
           
           # Main panel: Tables output ----
           mainPanel(
             div(style = "position:absolute;right:1em;", 
                 actionButton("run", "Run Inter-ViSTA", class = "btn-primary")),
    
             # Output: Nodes file ----
             h4("Nodes File"),
             dataTableOutput("nodescontents"),
             
             # Output: Edges file ----
             h4("Edges File"),
             dataTableOutput("edgescontents"),
             
             # Output: Background Gene List ----
             h4("Background Gene List"),
             dataTableOutput("bgcontents"),
             
             # Output: (Optional) Proteome Abundance File ----
             h4("(Optional) Proteome Abundance File"),
             dataTableOutput("proteomecontents"),
             
             # Output: (Optional) Gene Ontology Enrichment Analysis ----
             h4("(Optional) Gene Ontology Enrichment Analysis Table"),
             dataTableOutput("GOTablecontents")

           )
         )
)
