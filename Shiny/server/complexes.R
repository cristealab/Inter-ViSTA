# UI elements in Protein Complexes tab ----

output$complexes <- renderUI({
  
  if (is.null(computed_values()$complexes)) {
    
    verbatimTextOutput("No CORUM complexes had > 40% of their members detected in the input dataset.")
    
  } else {
    
    tagList(
      
      h2("Detected Protein Complexes from CORUM"),
      h4("CORUM protein complexes with > 40% of their members detected in the input dataset
         are displayed here."),
      br(),
      
      downloadButton("dComplexTable", "Download Protein Complexes Table"),
      
      br(),
      br(),
      
      dataTableOutput("complexesTable")
      
      )
    
  }
  
})


# server logic for complexes tab ----

## COMPLEXES

# download table
output$dComplexTable <- downloadHandler(
  filename = function() {
    "protein_complexes_table.csv"}
  ,
  content = function(file) {
    write.table(computed_values()$complexes, file = file, sep = ",", 
                quote = F, row.names = F)
  }
)

# print table output
output$complexesTable <- renderDataTable(computed_values()$complexes)
