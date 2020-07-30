
# server logic for instructions page

opts <- list(scrollX = TRUE, searching = FALSE, ordering = FALSE, info = FALSE)

output$nodesTable <- DT::renderDataTable(head(read.csv("Sample Datasets/Sample Uploads/nodes.csv", header = T, sep = ",")), 
                                                  options = opts)

output$edgesTable <- DT::renderDataTable(head(read.csv("Sample Datasets/Sample Uploads/edges.csv", header = T, sep = ",")), 
                                                  options = opts)

output$abundancesTable <- DT::renderDataTable(head(read.csv("Sample Datasets/Sample Uploads/abundances.csv", header = T, sep = ",")), 
                                          options = opts)