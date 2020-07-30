# UI elements in Gene Ontology tab ----

output$GOTerms <- renderUI({
  
  if (is.null(computed_values()$topGOterms)) {
    
    verbatimTextOutput("No significant gene ontology terms (p-value < 0.05) were found.")
    
  } else {
    
    tagList(
      
      h2("Significant gene ontology terms"),
      h4("Terms with corrected p-value < 0.05 are displayed; if no such terms were found,
         terms with (uncorrected) p-value < 0.05 are displayed."),
      br(),
      
      downloadButton("dGOTable", "Download Gene Ontology Table"),
      
      br(),
      br(),
      
      DT::dataTableOutput("GOTable")
      
    )
    
  }
  
})


# server logic for gene ontology tab ----

## GENE ONTOLOGY

# download table
output$dGOTable <- downloadHandler(
  filename = function() {
    "GO_enrichment_table.csv"}
  ,
  content = function(file) {
    write.table(computed_values()$GOTable, file = file, sep = ",", 
                quote = F, row.names = F)
  }
)

# print table output
output$GOTable <- DT::renderDataTable(computed_values()$GOTable)

