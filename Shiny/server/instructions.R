
# server logic for instructions page

opts <- list(scrollX = TRUE, searching = FALSE, ordering = FALSE, info = FALSE)

output$nodesTable <- renderDataTable(head(read.csv("Sample Datasets/Sample Uploads/nodes.csv", header = T, sep = ",")), 
                                                  options = opts)

output$edgesTable <- renderDataTable(head(read.csv("Sample Datasets/Sample Uploads/edges.csv", header = T, sep = ",")), 
                                                  options = opts)

output$abundancesTable <- renderDataTable(head(read.csv("Sample Datasets/Sample Uploads/abundances.csv", header = T, sep = ",")), 
                                          options = opts)