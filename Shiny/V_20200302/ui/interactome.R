# Tab panel: Interactome Network ----

tabPanel(title = "Interactome Network", value = "interactome", fluid = TRUE, 
         
         sidebarLayout(
           
           sidebarPanel(
             
             # Checkbox to display labels ----
             checkboxInput(inputId = "labels", "Display node labels", value=T),
             
             # checkbox to display timepoints ----
             checkboxInput(inputId = "display_timepoints", "Display timepoints below network", value=T),
             
             # Slider for STRING confidence ----
             sliderInput(inputId = "confidence", label="STRING edges confidence threshold:",
                         min=0, max=1, value=0.5),
             
             # UI output for variables dependent on user input ----
             # 1. radio buttons for shared proteins
             # 2. checkboxes for localizations
             # 3. checkboxes for GO terms
             # 4. checkboxes for species
             
             uiOutput("interactomeToggles"),
             
             # Action button to refresh network ----
             actionButton("renderNetwork", "Render Interactome Network", class = "btn-primary")
             
           ),
           
           mainPanel(
             
             uiOutput("interactomeMain")
             
           )
         )
)