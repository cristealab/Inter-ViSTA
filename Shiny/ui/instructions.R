

# UI for instructions tab

## Author: Michelle Kennedy
## Date: September 23, 2019

tabPanel(title = "Instructions", value = "instructions", fluid = TRUE,
         
         navlistPanel(
           
           #### Intro panel ####
           
           'INTRODUCTION',
           
           tabPanel('About Inter-ViSTA',
                    
                    verticalLayout(
                      
                      br(),
                      
                      div(img(src = 'logo.png', align = 'center', width = '500px'),
                          style="text-align: center;"),
                      
                      h1('Welcome to Inter-ViSTA', align = 'center'),
                      
                      h4('Interaction Visualization in Space and Time Analysis', align = 'center'),
                      
                      p(strong("Inter-ViSTA"), em("(Interaction Visualization in Space and Time Analysis)"), "is a computational platform to integrate spatial and temporal proteome and protein interaction information from user derived experiments with the wealth of knowledge in existing databases. Inter-ViSTA utilizes cutting-edge graph visualization algorithms and automatic programmatic access to protein databases to deliver an intuitive and user-friendly data visualization platform. The platform enables users to discover underlying patterns in multiple interaction networks across conditions using an interactive interface with dynamic network visuals and exploratory quantitative analysis. Users can integrate interaction networks of one or many baits, include subcellular localization information and functional annotations to modify networks, analyze dynamic quantitative properties of the network, and identify properties of proteins shared between multiple baits - all with their own data, with few manual steps.", align = 'justify'),
                      
                      p("Read on to discover how to use Inter-ViSTA to analyze your own datasets, or load a sample dataset to explore Inter-ViSTA's user-friendly interface, visualization, and analysis capacities.", align = 'justify'),
                      
                      br(),
 
                      # strong('If you want to be able to access the capabilities of Inter-ViSTA offline, you can download a portable-version of intervista on github by clicking', a('this link', href = 'https://github.com/cristealab/Inter-ViSTA'), style = 'color: #525252', align = 'justify'),
                 
                      
                      br(),
                      br(),
                      br(),
                      br()
                      
                      )
                    
                    ),
             
            # tabPanel('Citing this tool',
            #          
            #          h2('License'),
            #          
            #          # replace with correct link
            #          p('Inter-ViSTA is provided under the,', a('MIT License', href = 'https://github.com/cristealab/VISTA/blob/v1/LICENSE') ,'as defined on the Inter-ViSTA github repository'),
            #          
            #          h2('Citations'),
            #          
            #          # replace with correct journal
            #          p('If you find Inter-ViSTA useful in your research, please cite it via our initial publication of Inter-ViSTA by Federspiel et al. in ???.'),
            #          
            #          # replace with correct DOI
            #          a('DOI # here', href = 'https://github.com/cristealab/Inter-ViSTA')
            #          
            #          ),
           
           'ANALYZING YOUR OWN DATA',
           
           #### User data panel ####
           
           tabPanel('Tutorial video',
                    
                    verticalLayout(
                      
                      h1('Inter-ViSTA tutorial video'),
                      p('The video below demonstrates how to format and upload your own data for visualization and analysis in Inter-ViSTA'),
                    
                      # video only loads and plays in browser for some reason (but not in Joel's browser??)
                      tags$video(src = "tutorial.mp4", type = "video/mp4", controls = "controls", width = "100%")
                      
                    )),
           
           #########################################################################################################
           
           tabPanel('Uploading data',
                    
                    verticalLayout(
                      h1('How to format and upload your own data for Inter-ViSTA analysis'),
                      
                      wellPanel(style = "background: #347AB6; padding: 5px 15px 5px 15px; border-radius: 5px; margin: 10px 0px 10px 0px;",
                        h4(strong('STEP 1:  Upload your nodes, edges, and conditions'), align = 'justify', style = 'color: white')
                      ),
                    
                      fluidRow(align = 'justify',
                        column(width = 12, 
                               br(),
                          
                               img(src = 'nodes_edges_upload.png', width = '65%', style = 'padding: 5px 15px; float: right;', align = 'right'),
                               
                               p('Inter-ViSTA requires two standard inputs: a nodes file and an edges file. Both of these must be saved as "comma-separated value" or ".csv" files.'),
                               
                               br(),
                               
                               h4("The nodes file"),
                               
                               p('The nodes file defines each potential node in the dataset and must contain the UniProt accession number, gene name, and species of each node. Optionally, you may also include curated localizations for each node--if this column is not prodived, Inter-ViSTA will automatically retrieve node localization annotations from UniProt. It is important that the columns in your file are provided explicity in this order: accession, gene name, species taxid, localizations (optional).'),
                               
                               br(),
                               
                               h4("The edges file"),
                               
                               p('The edges file defines interactions between nodes in the dataset and includes information regarding the strength and confidence of their interaction. For each edge, this file should include the gene names and UniProt accessions for a pair of nodes,  a confidence value for the interaction (e.g. as provided by the SAINT algorithm), and, optionally, the abundance of the interaction across conditions. If an interaction is not present or detected at a given condition, its abundance and confidence values should be zero.')
                        )
                      ),

                      br(),
                      
                      h4("Experiments with multiple timepoints and/or conditions"),
                      
                      fluidRow(align = 'justify',
                        column(width = 12,
                                      
                               img(src = 'conditions_confidence.png', style = 'padding: 5px 15px; float: right;', width = '65%'),
                                            
                               p('One of the main utilities of Inter-ViSTA is visualizating interaction dynamics across different conditions. If your dataset includes more than one condition and/or timepoint, you must label the column headers of your edges file accordingly. For each condition, there should be a corresponding edges file column named "abundances_" and "confidence_score_", where the condition name is appended after the "_" symbol. You may name your conditions as you wish, but please ensure that the conditions entered into Inter-ViSTA match those in your edges file and that they are ordered as you desire.'),
                              
                               br(),
                              
                               p('After entering in your conditions, you can also adjust the interaction specificity cutoff. Edges between nodes that do not meet this confidence score threshold will not be displayed.')
                              
                        )
                      ),
                      
                      br(),
                      
                      h4("Formatting your data"),
                      
                      # hyperlink to add later
                      # or use our data formatting tool (', a('Inter-ViSTA formatting GUI', href = 'https://github.com/cristealab/Inter-ViSTA'), ')
                      p('See the example tables below to properly format your files to work with Inter-ViSTA.', align = 'justify'),
                      
                      br(),
                      
                      fluidRow(column(width = 5, style="text-align: center;",
                                      strong('example nodes table'),
                                      br(),
                                      img(src = 'nodes_example.png', width = '100%')
                                      ), 
                               
                               column(width = 7, style="text-align: center;",
                                      strong('example edges table'),
                                      br(),
                                      img(src = 'edges_example.png', width = '100%')
                                      )
                               ),
                  
                      br(),
                      
                      wellPanel(style = "background: #347AB6; padding: 5px 15px 5px 15px; border-radius: 5px; margin: 10px 0px 10px 0px;",
                        h4(strong('STEP 2:  Upload or select a background gene list for gene ontology enrichment analysis'), align = 'justify', style = 'color: white')
                        ),
                      
                      br(),
                      
                      fluidRow(column(width = 12, align = 'justify',
                                      img(src = 'background_upload.png', style = 'padding: 5px 15px; float: right;', width = '50%'),
                                      
                                      p('The final piece of required user-provided information is a background gene list to perform gene ontology (GO) enrichment. This is a list of UniProt accessions and taxonomic IDs for all proteins you wish to include as background. We have already pre-loaded several potentially-useful background datasets including those for human fibroblast tissue, Hela cells, and rat primary neurons.'),
                                      
                                      br(),
                      
                                      p('However, if you would like to create a custom background gene list, all you have to do is upload a .csv file of all the UniProt accession numbers you would like included. Inter-ViSTA will automatically query UniProt and retrieve GO ID annotations for these proteins.')
                                      )
                               ),
                      
                      br(),
                      
                      wellPanel(style = "background: #347AB6; padding: 5px 15px 5px 15px; border-radius: 5px; margin: 10px 0px 10px 0px;",
                        h4(strong('STEP 3 (optional): Normalize abundances to proteome abundance'), align = 'justify', style = 'color: white')
                      ),
                      
                      br(),
                      
                      fluidRow(column(width = 5, align = 'justify',
                                      p('Changes in interactions may be driven by either functional or proteomic abundance changes. To account for the latter, Inter-ViSTA can normalize the provided interaction abundances to proteome abundances and produce heatmaps of interaction abundances prior to and after normalization to the proteome. To perform this scaling, you may upload a file that contains UniProt accession numbers, gene names, and protein abundance at each timepoint or condition. The abundance columns should be labeled similarly to the edges file, i.e. "abundance_condition" for each condition in your dataset (see example below)'),
                                      br()
                                      ), 
                               
                               column(width = 7, style="text-align: center;",
                                      img(src = 'abundances_upload.png', width = '100%'), br(),
                                      br(),
                                      img(src = 'abundance_example.png', width = '100%')
                                      )
                      ),
                      p(strong('Note:'), em('Inter-ViSTA will not automatically use the normalized interaction abundances across the platform. If you wish to use the normalized abundances for further analysis, download these values from the nodes and edges files output by VISTA and replace the provided abundances with the new normalized abundances before running VISTA again.')),
                      br(), br(), br()
                      
                      )
                    ),
           
           #################################################################################################
           
           tabPanel('Running Inter-ViSTA and viewing your analysis',
                    verticalLayout(
                      
                      h1('The anatomy of an Inter-ViSTA analysis'),
                    
                      img(src = 'anatomy_network.png', align = 'center'),
                      
                      br(),
                       
                      fluidRow(
                        column(width = 12,
                               br(),
                               img(src = 'network_ul1337.png', width = '60%', style = 'padding: 5px 15px; float: right;'),
                               p(em('A major hurdle in analyzing an interactions study is visualizating the results in a meaningful way.'), align = 'justify'),
                               p(strong('Inter-ViSTA makes this EASY!'), align = 'justify'),
                               p('The first analysis tab is titled', strong('"Interactome Network"'), '. Use the sidebar on the left to select various properties of the network-STRING edge confidence threshold, shared protein attributes, localization (provided or UniProt), functional annotation, etc. Then click "Render Network Plot" to build the interactome, which will appear in the main panel. Each node represents a protein, while each edge is an interaction. Nodes are colored by protein localization, with outline colors indicating the duration for which the protein is associated with a bait. Edges are colored by type - user-provided or STRING-inferred, and edge thickness indicates interaction confidence. Click the play and pause buttons to view dynamic changes to PPIs; click nodes or edges to view their properties (localization, duration of activity, confidence, etc.); double-click a node to highlight all its neighbors. To modify the network, choose different parameters from the sidebar and click "Render Network Plot" again.', align = 'justify')
                        ) 
                        
                      ),
                      
                      p(),
                    
                      img(src = 'anatomy_quantplots.png', align = 'center'), br(),
                      br(),
                      p('Not only does Inter-ViSTA help you visualize your network, but its second tab titled ', strong('"Quantitative Plots"'),' outputs the results of a quantitative analysis that includes clustering to help you visualize broad trends within your dataset. If you are interested in the trends of one or more specific proteins, you can use the input box on the left to plot quantitative profiles for these proteins and their neighbors. Finally, if proteome abundance data if provided, heatmaps comparing interaction abundances prior to and after normalization will also be displayed in this tab.', align = 'justify'),
                      br(),
                      p(strong('Note:'), em('  All the plots from this page can be downloaded as a pdf file.'), align = 'justify'),
                      
                      fluidRow(
                        column(width = 6, style="text-align: center;",
                               br(),
                               em('Inter-ViSTA query for interactions of the catenin complex with viral protein pUL13'),
                               p('(from HCMV pUL13 and pUL37 sample dataset)'),
                               br(),
                               img(src = 'query_CTN.png', width = '100%'), br()
                               
                        ), 
                        
                        column(width = 6, style="text-align: center;",
                               br(),
                               em('Inter-ViSTA clustering of HCMV pUL13 and pUL37 sample dataset'), br(),
                               br(),
                               img(src = 'clusters_ul1337.png', width = '100%'), br()
                        ), 
                        
                        column(width = 12, style="text-align: center;",
                               br(),
                               em('Interactions with viral protein pUL37 normalized by protein abundance'),
                               br(),
                               img(src = 'heatmaps.png', width = '100%'), br()
                               
                               )
                      ),
  
                      br(),
                      
                      img(src = 'anatomy_GO.png', align = 'center'), br(),
                      br(),
                      p('The third tab, titled ', strong('"Gene Ontology"'), ', outputs an interactive data table with significant GO terms found in the data. Significance is assessed as having a corrected p-value of < 0.05 after multiple-test correction is applied on Fisher test p-values for enrichment. Data in the table can be sorted, searched, and filtered.', align = 'justify'),
                      br(),
                      p(strong('Note:'), em('  This table can be downloaded as a tab-separated file.'), align = 'justify'),
                      
                      br(),
                      
                      img(src = 'anatomy_complexes.png', align = 'center'), br(),
                      br(),
                      p('The fourth tab, titled ', strong('"Protein Complexes"'), ' outputs an interactive data table with a list of mammalian protein complexes from CORUM in which over 40% of members are detected in the user-provided data. Inter-ViSTA provides the complex name, fraction of complex components that were detected, and lists the identified proteins in each complex. Data in the table can be sorted, searched, and filtered.', align = 'justify'),
                      br(),
                      p(strong('Note:'), em('  This table can be downloaded as a tab-separated file.'), align = 'justify'),
                      br(),
                      
                      img(src = 'anatomy_downloads.png', align = 'center'), br(),
                      br(),
                      p('The last tab is ', strong('"Downloads"'), ', which provides links to download the annotated nodes and edges files produced by Inter-ViSTA. These annotated files contain information on protein and interaction durations, localization, GO terms, cluster membership, etc. and can be directly exported to Cytoscape. Finally, an interactome "report" that provides exploratory analysis on the duration of activity for various proteins in the network, separated by their localizations and functional annotations, as well as plots of interaction onsets and durations that can be downloaded as a pdf.', align = 'justify'),
                      br(),
                      br(),
                      br()
                    )
            ),
           
           tabPanel('Saving and loading Inter-ViSTA analyses from memory',
                    verticalLayout(
                      br(),
                      h2('Want to come back to your analysis later?'),
                     
                      fluidRow(
                        column(width = 6, align = 'justify',
                               h3('Save some time!'),
                               br(),
                               p("Initially, running Inter-ViSTA may take some time (up to several minutes or longer depending on the speed of your internet connection). This is mostly because querying and downloading responses from API's like UniProt and STRING can be slow. However, once Inter-ViSTA has done this once for your dataset, it doesn't need to do it again!"),
                               br(),
                               p("If you are finished (or need to take a break) from your analysis, download the 'Node Attributes', 'Edge Attributes', and 'UniProt Annotated Gene List' files to your local machine."),
                               br(),
                               p("The next time you want to analyze the same data, simply upload these files in the appropriate fields, and Inter-ViSTA will render your analysis in less than 30 seconds!")
                            ),
                        
                        column(width = 5, style="text-align: center;", 
                               style = "background: #ffffff; border-style: solid; border-color: #6E84AD; border-width: 5px; border-radius: 5px; margin: 10px;",
                               
                               img(src = 'to_download_edited.png', width = '100%', style = 'padding: 5px 15px; float: right;'), br(),
                               img(src = 'arrow.png', width = '25%'), br(),
                               img(src = 'load_prev.png', width = '100%', style = 'padding: 5px 15px; float: right;'), br()

                             )
                      ),
                      br(),
                      p()
  
                    )
                    )
           
           #### examples panel ####
           
           # 'EXAMPLES',
           # 
           # tabPanel('Analyzing sample datasets',
           #          verticalLayout(
           # 
           # 
           # 
           #          )
           # )
      
      )
)
        