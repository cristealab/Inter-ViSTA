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
            
             # action button to calculate plot and switch tabs to plots
             actionButton("abundGo", "Render Quantitative Plot", class = "btn-primary")
             
           ),
           
           # Main panel to display quantitative plots ----- 
           mainPanel(

             uiOutput("quantplotsMain")
             
           )
           
         )
         
)