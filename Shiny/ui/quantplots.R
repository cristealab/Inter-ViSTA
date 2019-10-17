# Tab panel: Quantitative Plots ----

tabPanel(title = "Quantitative Plots", value = "quantplots", fluid = TRUE,
         
         sidebarLayout(
           
           # Sidebar for inputs to control quantitative plots ----
           sidebarPanel(
             
             # textbox for genes to plot
             textInput("genesym", label = "See quantitative values for:", 
                       value = "Ex: IFT27; CCL19; TRAV35"),
             
             # checkbox to include node neighbors in plot
             checkboxInput(inputId = "neighbors", "Display interactor quants", value=F),
             
             br(),
             
             p('Enter gene names (or partial gene names) above to show their quantitative abundance plot(s). Separate multiple gene names with a semicolon (;) character'),
             
             br(),
             
             p(em('Note: To search explicitly for a given gene or genes, enclose your query with quotes.')),
             
             p(em('Ex: "UL12" will only return the profile for gene UL12, while UL12 without quotes will return any gene containing "UL12" within its name')),
             
             br(),
             
             # action button to calculate plot and switch tabs to plots
             actionButton("abundGo", "Render Quantitative Plot", class = "btn-primary")
             
           ),
           
           # Main panel to display quantitative plots ----- 
           mainPanel(
             
             uiOutput("quantplotsMain")
             
           )
           
         )
         
)