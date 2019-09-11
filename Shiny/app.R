library(shiny)
library(shinyBS)
library(shinyjs)
library(zip)

##### GET AND CLEAN DATA ####

cleanData <- function(nodes_in, edges_in, spec_threshold) {
  
  nodes_out <- nodes_in
  edges_out <- edges_in
  
  # ensure all values are filled
  edges_out[is.na(edges_out)] <- 0
  
  # get edge weight threshold
  # weight threshold is set lower than cutoff to prevent artificial on-off behavior
  WEIGHT_THRESHOLD <- NULL
  if (spec_threshold <= 0.1 ) { 
    WEIGHT_THRESHOLD <- spec_threshold 
  } else {
    WEIGHT_THRESHOLD <- spec_threshold - 0.1
  } 
  
  # get timepoint information from column names
  num_link_attr <- ncol(edges_out)
  timepoints <- unique(as.numeric(unlist(regmatches(colnames(edges_out)[3:num_link_attr],
                                                    gregexpr("[[:digit:]]+\\.*[[:digit:]]*", 
                                                             colnames(edges_out)[3:num_link_attr])))))
  NUM_TPS <- length(timepoints)
  
  # set interval length to smallest interval between timepoints
  interval <- {}
  for (i in 1:(NUM_TPS-1)) {
    interval[i] <- timepoints[i+1] - timepoints[i]
  }
  INTERVAL <- min(interval)
  
  # rename node columns
  LOC_PROV <- F
  if (ncol(nodes_out) == 3) { 
    colnames(nodes_out) <- c("accession", "gene_name", "taxid") 
  } else { 
    LOC_PROV <- T
    colnames(nodes_out) <- c("accession", "gene_name", "taxid", "localization") 
  }
  
  # rename link columns 
  colnames(edges_out)[c(1,2)] <- c("bait_accession", "prey_accession")
  ABUND_PROV <- F
  if (num_link_attr == NUM_TPS+2) {
    colnames(edges_out)[3:num_link_attr] <- paste("w", timepoints, sep="")
  } else { 
    ABUND_PROV <- T
    cnames <- c()
    for (i in 1:NUM_TPS) {
      cnames <- c(cnames, paste("w", timepoints[i], sep=""))
      cnames <- c(cnames, paste("a", timepoints[i], sep=""))
    }
    colnames(edges_out)[3:num_link_attr] <- cnames
  }
  
  # only keep links where at least one timepoint is over the cutoff threshold
  keep_links <- apply(edges_out[, paste("w", timepoints, sep="")], 1, 
                      function(x) { sum(x >= spec_threshold) > 0 } )
  edges_out <- edges_out[keep_links, ]
  
  # remove duplicate nodes and links
  nodes_out <- nodes_out[!duplicated(nodes_out$accession), ]
  nodes_out <- nodes_out[nodes_out$accession %in% edges_out$bait_accession | 
                           nodes_out$accession %in% edges_out$prey_accession, ]
  NUM_NODES <- dim(nodes_out)[1]
  
  edges_out <- edges_out[!duplicated(edges_out[, c("bait_accession", "prey_accession")]), ]
  NUM_EDGES <- dim(edges_out)[1]
  
  # assign unique ids to nodes
  nodes_out$id <- c(1:NUM_NODES)
  rownames(nodes_out) <- nodes_out$id
  
  # attach ids and gene names to links
  edges_out <- merge(edges_out, nodes_out[, c("accession", "gene_name", "id")], 
                     by.x="bait_accession", 
                     by.y="accession")
  colnames(edges_out)[c(ncol(edges_out)-1, ncol(edges_out))] <- c("bait_gene_name", "bait_id")
  edges_out <- merge(edges_out, nodes_out[, c("accession", "gene_name", "id")], 
                     by.x="prey_accession", 
                     by.y="accession")
  colnames(edges_out)[c(ncol(edges_out)-1, ncol(edges_out))] <- c("prey_gene_name", "prey_id")
  
  BAIT_IDS <- unique(edges_out$bait_id)
  NUM_BAITS <- length(BAIT_IDS)
  
  TAXIDS <- unique(nodes_out$taxid)
  
  # calculate shared nodes
  
  for (tp in timepoints) {
    shared_tp <- c(mode='numeric', length=NUM_NODES)
    shared_bait_tp <- c("", length=NUM_NODES)
    for (i in 1:NUM_NODES) {
      shared_tp[i] <- sum(edges_out[edges_out$prey_accession == nodes_out$accession[i], 
                                    paste('w', tp, sep="")] >= WEIGHT_THRESHOLD)
      if (shared_tp[i] > 0) {
        shared_baits <- edges_out[which(edges_out$prey_accession == nodes_out$accession[i] & 
                                          edges_out[, paste('w', tp, sep="")] >= WEIGHT_THRESHOLD), 
                                  "bait_gene_name"]
        shared_bait_tp[i] <- paste0(shared_baits, collapse = "; ")
      } else { shared_bait_tp[i] <- "" }
    }
    nodes_out <- cbind(nodes_out, shared_tp)
    colnames(nodes_out)[ncol(nodes_out)] <- paste("shared_", tp, sep="")
    nodes_out[, ncol(nodes_out)] <- as.character(nodes_out[, ncol(nodes_out)])
    
    nodes_out <- cbind(nodes_out, shared_bait_tp)
    colnames(nodes_out)[ncol(nodes_out)] <- paste("shared_baits_", tp, sep="")
    nodes_out[, ncol(nodes_out)] <- as.character(nodes_out[, ncol(nodes_out)])
  }
  
  
  # return list of: nodes, edges, weight threshold, timepoints, interval length,
  # loc_prov, abund_prov, bait IDs, taxids
  
  cleaned_data <- list("nodes" = nodes_out,
                       "edges" = edges_out,
                       "num_edges" = NUM_EDGES,
                       "weight_threshold" = WEIGHT_THRESHOLD,
                       "timepoints" = timepoints,
                       "interval" = INTERVAL,
                       "loc_prov" = LOC_PROV,
                       "abund_prov" = ABUND_PROV,
                       "bait_ids" = BAIT_IDS,
                       "taxids" = TAXIDS
  )
  return (cleaned_data)
}

##### CORUM COMPLEX INFORMATION ####  

calculateComplexes <- function(nodes_in) {
  
  nodes_out <- nodes_in
  NUM_NODES <- dim(nodes_out)[1]
  
  corum_complexes <- read.table("AppFiles/coreComplexes.txt", header = T, sep = "\t", quote = "", 
                                fill = T, stringsAsFactors = F)
  corum_complexes <- corum_complexes[, c("ComplexName", "subunits.UniProt.IDs.")]
  
  list_complexes <- strsplit(corum_complexes[, "subunits.UniProt.IDs."], ";")
  names(list_complexes) <- corum_complexes$ComplexName
  
  # only keep complexes with at least 3 members (no duplexes)
  list_complexes <- Filter(function(x) length(x) > 2, list_complexes)
  
  # only keep complexes with at least 40% of its members in the provided dataset
  keep_function <- function(x) {
    if (sum(x %in% nodes_out$accession) / length(x) >= 0.4 ) { 
      return (TRUE) 
    } else { return (FALSE) }
  }
  list_complexes <- Filter(keep_function, list_complexes)
  
  chr2string <- function(x) {
    gene_names <- nodes_out[nodes_out$accession %in% x, "gene_name"]
    return (paste0(gene_names, collapse="; "))
  }
  
  percent_detected <- function(x) {
    return ( sum(x %in% nodes_out$accession) / length(x) )
  }
  
  # output complexes table (for download)
  
  if (length(list_complexes) != 0) {
    
    detected_complexes <- cbind.data.frame(names(list_complexes), 
                                           sapply(list_complexes, percent_detected),
                                           sapply(list_complexes, chr2string))
    colnames(detected_complexes) <- c("Complex Name", "Fraction of Complex Detected",
                                      "Detected Complex Members")
    detected_complexes <- detected_complexes[order(-detected_complexes$`Fraction of Complex Detected`), ]
    
    # annotate nodes with complex membership
    complex_membership <- rep("", NUM_NODES)
    for (i in 1:NUM_NODES) {
      complexes <- which(grepl(nodes_out$gene_name[i], detected_complexes$`Detected Complex Members`))
      if (length(complexes) > 0) {
        complex_membership[i] <- paste0(names(list_complexes)[complexes], collapse="; ")
      }
    }
    nodes_out$complexes <- complex_membership
    
  } else {
    
    # no complexes found in CORUM
    
    detected_complexes <- NULL
    nodes_out$complexes <- NA 
    
  }
  
  # function returns: nodes annotated with complex membership and detected complexes
  
  return (list("nodes" = nodes_out,
               "complexes" = detected_complexes))
}

##### QUANTITATIVE FUNCTIONS ####

## CLUSTER PROTEINS BASED ON ABUNDANCES OVER TIME

# only use this function if abundances provided

clusterProteins <- function(nodes_in, edges_in, bait_ids, timepoints) {
  
  library(ggplot2)
  theme_set(theme_bw())
  
  edges_out <- edges_in
  
  edges_out$cluster = NA
  relative_abundances <- {}
  
  for (bait_id in bait_ids) {
    bait_gene_name <- nodes_in[which(nodes_in$id == bait_id), "gene_name"]
    
    # get abundances for all nodes associated with this bait (except bait itself)
    abundances <- edges_out[which(edges_out$bait_id == bait_id & edges_out$prey_id != bait_id), 
                            c("prey_gene_name", paste("a", timepoints, sep=""))]
    
    # scale values from 0 to 1
    zero2one <- function(x) {
      return ((x - min(x)) / (max(x) - min(x)) )
    }
    
    # scale values to max 1
    maxOne <- function(x) {
      return ( x / max(x) )
    }
    
    to_scale <- abundances[, c(2:ncol(abundances))]
    abundances[, c(2:ncol(abundances))] <- t(apply(to_scale, 1, maxOne))
    
    # calculate pca
    abundances[is.na(abundances)] <- 0
    # only calculate pca on columns with differential expresson
    comp_pca <- c()
    for (column in 2:ncol(abundances)) {
      if (dim(table(abundances[, column])) > 1) { 
        comp_pca <- c(comp_pca, column) 
      }
    }
    abundances_pca <- prcomp(scale(abundances[, comp_pca]))
    # find number of clusters as number of PCs that explain over 90% of the variance
    eigs <- abundances_pca$sdev^2 
    pve <- eigs / sum(eigs) # get percent of variance explained by each PC
    total_pve <- 0
    for (i in 1:length(pve)) {
      total_pve <- total_pve + pve[i]
      if (total_pve >= 0.9) { break }
    }
    num_clusters <- i
    
    # k-means clustering
    
    clusters <- kmeans(abundances[, c(2:ncol(abundances))], num_clusters)
    abundances <- cbind(abundances, clusters$cluster)
    colnames(abundances)[ncol(abundances)] <- "cluster"
    
    # add cluster numbers to links 
    for (i in 1:dim(abundances)[1]) {
      edges_out$cluster[which(edges_out$prey_gene_name == abundances$prey_gene_name[i] & 
                                edges_out$bait_gene_name == bait_gene_name)] <- 
        abundances$cluster[i]
    }
    
    # save relative abundances to plot
    abundances$bait_gene_name <- bait_gene_name
    relative_abundances <- rbind(relative_abundances, abundances)
  }
  
  # return list of: edges, relative abundances file
  return (list("edges" = edges_out,
               "abundances" = relative_abundances))
}

# function to plot individual protein profiles split by cluster

clusterPlot <- function(bait_id, nodes_in, relative_abundances, timepoints) {
  
  # plot profiles of proteins split by cluster
  melt_abund <- {}
  bait_gene_name <- nodes_in$gene_name[which(nodes_in$id == bait_id)]
  abundances <- relative_abundances[which(relative_abundances$bait_gene_name == 
                                            bait_gene_name), ]
  
  NUM_TPS <- length(timepoints)
  
  for (i in rownames(abundances)) {
    new <- cbind.data.frame(timepoints, rep(abundances[i, "cluster"], NUM_TPS))
    new <- cbind.data.frame(new, rep(abundances[i, "prey_gene_name"], NUM_TPS))
    new <- cbind.data.frame(new, as.vector(t(abundances[i, paste("a", timepoints, sep="")])))
    melt_abund <- rbind(melt_abund, new)
  }
  colnames(melt_abund) <- c("time", "cluster", "prey_gene_name", "abundance")
  rownames(melt_abund) <- 1:dim(melt_abund)[1]
  
  print(ggplot(data=melt_abund, aes(x=time, y=abundance, group=prey_gene_name)) +
          geom_line(aes(color=as.factor(cluster))) +
          geom_point(aes(color=as.factor(cluster))) +
          scale_x_continuous(breaks = timepoints) +
          labs(title=paste("Bait", bait_gene_name), x="Time", y="Scaled Relative Abundance", 
               color="Cluster") +
          facet_wrap(~as.factor(cluster)))
  
}

## NORMALIZE BY PROTEOME ABUNDANCE 

# only use this function if abundances provided and proteome abundances file provided

normalizeByProteomeAbund <- function(proteome_abundance, nodes_in, edges_in, timepoints) {
  
  library(ggplot2)
  theme_set(theme_bw())
  library(gridExtra)
  
  edges_out <- edges_in
  
  colnames(proteome_abundance) <- c("accession", "gene_name", 
                                    paste("prot_a", timepoints, sep=""))
  
  proteome_abundance <- proteome_abundance[proteome_abundance$accession %in% 
                                             nodes_in$accession, ]
  
  heatmap_abundance <- edges_out[, c("bait_gene_name", "prey_gene_name", "cluster",
                                     paste("a", timepoints, sep=""))]
  
  NUM_EDGES <- dim(edges_out)[1]
  
  for (tp in timepoints) {
    norm_abund <- vector(mode='numeric', length=NUM_EDGES)
    for (i in 1:NUM_EDGES) {
      prey <- edges_out$prey_accession[i]
      if (prey %in% proteome_abundance$accession) {
        norm_abund[i] <- edges_out[i, paste("a", tp, sep="")] / 
          proteome_abundance[proteome_abundance$accession == prey, paste("prot_a", tp, sep="")]
      } else {
        norm_abund[i] <- NA
      }
    }
    heatmap_abundance <- cbind.data.frame(heatmap_abundance, norm_abund, stringsAsFactors = F)
    colnames(heatmap_abundance)[ncol(heatmap_abundance)] <- paste("norm_a", tp, sep="")
    heatmap_abundance[, paste("norm_a", tp, sep="")] <- 
      as.numeric(heatmap_abundance[, paste("norm_a", tp, sep="")])
  }
  
  scaled_reg <- t(apply(heatmap_abundance[, c(paste("a", timepoints, sep=""))], 1, 
                        zero2one))
  scaled_norm <- t(apply(heatmap_abundance[, c(paste("norm_a", timepoints, sep=""))], 1, 
                         zero2one))
  
  heatmap_abundance_scaled <- cbind.data.frame(heatmap_abundance[, c("bait_gene_name",
                                                                     "prey_gene_name",
                                                                     "cluster")],
                                               scaled_reg, scaled_norm, stringsAsFactors = F)
  
  # append data to edge attributes file
  edges_out <- cbind.data.frame(edges_out, heatmap_abundance[, c(paste("norm_a", timepoints, 
                                                                       sep=""))])
  
  # this should return edges and heatmap dataframe
  return (list("edges" = edges_out,
               "normalized_abundances" = heatmap_abundance_scaled))
}

# function to plot heatmaps of abundances normalized to proteome and not-normalized

plotHeatmaps <- function(bait_id, nodes_in, normalized_abundances, timepoints) {
  bait_gene_name <- nodes_in$gene_name[which(nodes_in$id == bait_id)]
  
  bait_data <- normalized_abundances[normalized_abundances$bait_gene_name == bait_gene_name, ]
  
  clusters <- sort(unique(bait_data$cluster))
  grobList <- list()
  
  NUM_TPS <- length(timepoints)
  
  for (cluster in clusters) {
    cluster_data <- bait_data[bait_data$cluster == cluster, c("prey_gene_name", 
                                                              paste("a", timepoints, 
                                                                    sep=""),
                                                              paste("norm_a", 
                                                                    timepoints, sep=""))]
    prey_gene_name <- rep(cluster_data$prey_gene_name, each=NUM_TPS*2)
    type <- rep(c("Not Norm.", "Prot. Norm."), each=NUM_TPS, times=dim(cluster_data)[1])
    time <- rep(timepoints, times=dim(cluster_data)[1]*2)
    values <- c(t(cluster_data[, -1]))
    cluster_data_melt <- cbind.data.frame(prey_gene_name, type, time, values)
    
    grobList[[cluster]] <- ggplot(cluster_data_melt, aes(time, prey_gene_name)) + 
      geom_tile(aes(fill = values), colour = "white") + 
      scale_fill_gradient(low = "navyblue", high = "yellow3", 
                          name = "Scaled Abund.") + 
      facet_wrap(aes(type)) +
      xlab("Condition") + 
      ylab("Prey Gene Name") +
      scale_x_continuous(breaks = timepoints, labels = timepoints) +
      ggtitle(paste("Cluster", cluster))
  }
  do.call("grid.arrange", c(grobList, ncol=floor(sqrt(length(clusters))), 
                            top = paste("Bait", bait_gene_name)))
}

## SIZE NODES BY ABUNDANCE FOR SINGLE-BAIT 

# only use this function if abundances provided, there is only one bait, and fewer than 
# 500 nodes in all

sizeByAbund <- function(nodes_in, edges_in, bait_ids, timepoints) {
  
  nodes_out <- nodes_in
  
  NUM_TPS <- length(timepoints)
  
  # convert abundance to size and size by quartile
  
  overall_avg <- vector(mode='numeric', length=NUM_TPS)
  for (i in 1:NUM_TPS) {
    overall_avg[i] <- mean(edges_in[, paste('a', timepoints[i], sep="")])
  }
  
  # set 2.5 as the average node size
  mult_factor <- 2.5 / mean(overall_avg)
  
  returnSize <- function(x) {
    y <- max(1, mult_factor*x) # prevent nodes from being too small
    y <- min(4, y) # prevent nodes from being too big
    return (y)
  }
  
  for (i in 1:NUM_TPS) {
    size <- lapply(edges_in[, paste('a', timepoints[i], sep="")], returnSize)
    size <- cbind(edges_in$prey_id, unlist(size))
    size <- rbind(size, c(bait_ids, 2.5))
    colnames(size) <- c("id", paste("size_", timepoints[i], sep=""))
    nodes_out <- merge(nodes_out, size, by="id")
  }
  
  return (nodes_out)
}


##### GET NODE ANNOTATIONS FROM UNIPROT ####

# get UniProt data for proteins in each species and background list for enrichment
# and assign localizations to nodes

getUniprotData <- function(background_list, nodes_in, taxids, loc_prov) {
  
  library(reticulate)
  
  nodes_out <- nodes_in
  
  if (ncol(background_list) > 1) {
    # if the provided list already has Uniprot information, skip the extraction step for 
    # the background genes, but perform the extraction for nodes being visualized by Inter-ViSTA
    uniprot_list <- nodes_out$accession
    
  } else {
    # otherwise add any nodes whose accessions do not appear in the background list
    uniprot_list <- c(nodes_out$accession, background_list[,1])
    uniprot_list <- unlist(unique(uniprot_list))
  }
  
  source_python("AppFiles/UniprotQuery.py")
  source_python("AppFiles/ruq.py")
  
  # break into chunks of 2000 if larger than 2000 
  CHUNK_SIZE <- 2000

  NUM_CHUNKS <- (length(uniprot_list) %/% CHUNK_SIZE)+1
  
  uniprot_data <- list()
  
  for (i in 1:NUM_CHUNKS) {
    
    n_start <- CHUNK_SIZE*(i-1) + 1
    n_end <- min(CHUNK_SIZE*i, length(uniprot_list))
    
    uniprot_chunk <- uniprot_list[n_start:n_end]
    
    uniprot_data <- append(uniprot_data, run_uniprot_query(uniprot_chunk))
    
  }
  
  # merge organism name into nodes file
  given_uniprot <- uniprot_data[nodes_out$accession]
  
  nodes_out$organism <- lapply(names(given_uniprot), 
                               function(nm) { given_uniprot[[nm]][["organism"]] } )
  nodes_out[which(nodes_out$organism == "NULL"), "organism"] = NA
  
  # create mapping from taxid to organism
  txidToOrg <- list()
  for (tx in taxids) {
    names <- nodes_out[which(nodes_out$taxid == tx), "organism"]
    names <- names[!is.na(names)]
    txidToOrg[[toString(tx)]] <- names[1]
  }
  # if organism name was NA, replace with known 
  unspecified <- which(is.na(nodes_out$organism))
  for (u in unspecified) {
    nodes_out[u, "organism"] <- txidToOrg[toString(nodes_out[u, "taxid"])]
  }
  ORGANISMS <- unique(nodes_out$organism)
  
  # create node Uniprot annotations
  na_accession <- names(given_uniprot)
  na_genename <- sapply(names(given_uniprot), function(nm) { 
    paste(given_uniprot[[nm]][["genename"]], collapse = "; " ) } )
  na_GOIDs <- sapply(names(given_uniprot), function(nm) { 
    paste(given_uniprot[[nm]][["GO IDs"]], collapse = "; " ) } )
  na_localizations <- sapply(names(given_uniprot), function(nm) { 
    paste(given_uniprot[[nm]][["localizations"]], collapse = "; " ) } )
  na_organism <- sapply(names(given_uniprot), function(nm) { 
    paste(given_uniprot[[nm]][["organism"]], collapse = "; " ) } )
  
  node_annotations <- cbind(na_accession, na_genename, na_GOIDs, na_localizations, na_organism)
  colnames(node_annotations) <- c("accession", "genename", "GOIDs", "localizations", "organism")
  node_annotations <- as.data.frame(node_annotations, stringsAsFactors = F)
  
  # if this is the first time program was run on this background list, save the background
  # gene list with all Uniprot annotations
  if (ncol(background_list) == 1) {
    background_uniprot <- uniprot_data[background_list[, 1]]
    
    up_accession <- names(background_uniprot)
    up_genename <- sapply(names(background_uniprot), function(nm) { 
      paste(background_uniprot[[nm]][["genename"]], collapse = "; " ) } )
    up_GOIDs <- sapply(names(background_uniprot), function(nm) { 
      paste(background_uniprot[[nm]][["GO IDs"]], collapse = "; " ) } )
    up_localizations <- sapply(names(background_uniprot), function(nm) { 
      paste(background_uniprot[[nm]][["localizations"]], collapse = "; " ) } )
    up_organism <- sapply(names(background_uniprot), function(nm) { 
      paste(background_uniprot[[nm]][["organism"]], collapse = "; " ) } )
    
    uniprot_df <- cbind(up_accession, up_genename, up_GOIDs, up_localizations, 
                        up_organism)
    colnames(uniprot_df) <- c("accession", "genename", "GOIDs", "localizations", "organism")
    
    uniprot_df <- as.data.frame(uniprot_df, stringsAsFactors = F)
    
    # new bg list has all the annotations included
    background_list <- uniprot_df 
  }
  
  all_annotations <- rbind(node_annotations, background_list)
  # create functional GO annotations list to later pipe into topGO
  accessionToGO <- strsplit(all_annotations$GOIDs, "; ")
  names(accessionToGO) <- all_annotations$accession
  
  ## ASSIGN LOCALIZATIONS
  
  # localizations already provided
  if (loc_prov) { 
    # parse all given localizations along separator ";"
    locs <- unique(unlist(strsplit(nodes_out$localization, "; ")))
    NUM_LOCS <- length(locs)
    
    NUM_NODES <- dim(nodes_out)[1]
    
    for (loc in locs) {
      check <- vector(mode="logical", length=NUM_NODES)
      check <- grepl(loc, nodes_out$localization, ignore.case=T)
      nodes_out <- cbind(nodes_out, check)
      colnames(nodes_out)[ncol(nodes_out)] <- paste("loc_", loc, sep="")
    }
  } else { 
    locs <- c("Cell Membrane", "Cytoplasm", "Endoplasmic Reticulum", "Golgi", 
              "Mitochondrion", "Nucleus", "Peroxisome")
    NUM_LOCS <- length(locs)
    all_locs <- node_annotations$accession
    # assign locs based on whether 'subcellular location' comment in UniProt mentions location
    for (loc in locs) {
      x <- as.numeric(grepl(loc, node_annotations$localizations, ignore.case=T))
      all_locs <- cbind.data.frame(all_locs, x)
      colnames(all_locs)[ncol(all_locs)] <- paste("loc_", loc, sep="")
    }
    colnames(all_locs)[1] <- "accession"
    nodes_out <- merge(nodes_out, all_locs, by="accession", all.x=T)
    nodes_out <- nodes_out[!duplicated(nodes_out), ]
    nodes_out[is.na(nodes_out)] <- 0
  }
  
  NUM_NODES <- dim(nodes_out)[1]
  
  # convert localization information into single string
  localization_string <- vector(mode='character', length=NUM_NODES)
  localization_string_output <- vector(mode='character', length=NUM_NODES)
  
  for (loc in locs) {
    index <- which(nodes_out[, paste('loc_', loc, sep="")] == 1)
    if (length(index) > 0) {
      localization_string_output[index] <- 
        paste(localization_string_output[index], loc, sep="; ")
      for (j in index) {
        if (localization_string[j] == "" ) { localization_string[j] <- loc; }
        # multiple possiblities
        else { localization_string[j] <- "Multiple" }
      }
    }
  }
  nodes_out$localization_string <- localization_string
  nodes_out$localization_string_output <- localization_string_output
  
  # if node lacks localization information replace with organism name
  unspecified <- which(nodes_out$localization_string == "")
  for (u in unspecified) {
    nodes_out[u, "localization_string"] <- txidToOrg[toString(nodes_out[u, "taxid"])]
  }
  
  # function returns: nodes, organisms, localizations, background_list (out), GOList
  return (list("nodes" = nodes_out,
               "organisms" = ORGANISMS,
               "localizations" = locs,
               "background" = background_list,
               "GOList" = accessionToGO))
}

##### ASSIGN FUNCTIONAL ANNOTATIONS AFTER GO ENRICHMENT #####

performGOEnrichment <- function(GOList, nodes_in) {
  
  library(topGO)
  
  nodes_out <- nodes_in
  
  gene_universe <- names(GOList)
  interesting_genes <- nodes_out$accession
  
  # assign whether a gene is "interesting" - in the provided data
  gene_list <- factor(as.integer(gene_universe %in% interesting_genes))
  names(gene_list) <- gene_universe
  
  # create topGO object
  # Biological Process (BP) ontology
  GOdata <- new("topGOdata", ontology = "BP", allGenes = gene_list, 
                annot = annFUN.gene2GO, gene2GO = GOList)
  
  # use Fisher's test to determine enrichment
  test <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher <- getSigGroups(GOdata, test)
  
  # perform multiple test correction using fdr (BH) method
  fisher_p_values <- score(resultFisher)
  corr_p_values <- p.adjust(fisher_p_values, method = "BH")
  corrected_p_values <- as.data.frame(corr_p_values)
  corrected_p_values$GO.ID <- rownames(corrected_p_values)
  
  # generate results table
  allRes <- GenTable(GOdata, Fisher_p_value = resultFisher,
                     topNodes = dim(corrected_p_values)[1])
  allRes <- merge(allRes, corrected_p_values, by="GO.ID")
  # only keep terms with p-value < 0.05 
  # (if not more than 5 from multiple correction, then just Fisher)
  if (length(which(allRes$corr_p_values < 0.05)) > 5) {
    allRes <- allRes[which(allRes$corr_p_values < 0.05), ]
    # sort by corrected p_values
    allRes <- allRes[order(allRes$corr_p_values), ]
  } else {
    allRes <- allRes[which(allRes$Fisher_p_value < 0.05), ]
    # sort by Fisher p_values
    allRes <- allRes[order(allRes$Fisher_p_value), ]
  }
  
  colnames(allRes)[ncol(allRes)] <- "Corrected p-value"
  
  # if significant GO terms are found
  if (dim(allRes)[1] > 0) {
    
    # get list of top (upto 10) GO terms
    num_go_terms <- min(dim(allRes)[1], 10)
    TOP_GO_TERMS <- allRes$Term[1:num_go_terms]
    
    # list the genes for each of the significant GO terms
    annotated_genes <- lapply(allRes$GO.ID,
                              function(x) as.character(
                                unlist(genesInTerm(object = GOdata, whichGO = x))))
    names(annotated_genes) <- allRes$Term
    # get GO terms for the original genes in the provided list
    significant_genes <- lapply(annotated_genes, function(x) intersect(x, interesting_genes))
    
    # create columns for each GO term designating whether a protein is associated with the term
    # and string with all of top 10 significant GO terms
    
    NUM_NODES <- dim(nodes_out)[1]
    nodes_out$GO_term_string <- vector(mode = "character", length = NUM_NODES)
    for (got in TOP_GO_TERMS) {
      genes_in_got <- rep(0, NUM_NODES)
      i <- which(nodes_out$accession %in% unlist(significant_genes[got]))
      nodes_out$GO_term_string[i] <- paste(nodes_out$GO_term_string[i], got, sep="; ")
      genes_in_got[i] <- 1
      nodes_out <- cbind(nodes_out, genes_in_got)
      colnames(nodes_out)[ncol(nodes_out)] <- paste("got_", got, sep="")
    }
    
    # output GO Table (for download)
    GOTable <- allRes
    NUM_GO_TERMS <- dim(GOTable)[1]
    GOTable$Proteins <- vector(mode = "character", length = NUM_GO_TERMS)
    
    for (i in 1:NUM_GO_TERMS) {
      GOTable$Proteins[i] <- paste0(nodes_out$gene_name[which(nodes_out$accession %in%
                                                                unlist(significant_genes[GOTable$Term[i]]))],
                                    collapse = "; ")
    }
    
  } else {
    
    # no significant GO terms found
    TOP_GO_TERMS <- NULL
    nodes_out$GO_term_string <- NA
    GOTable <- NULL
  }
  
  # function returns: nodes, GO terms table, top GO terms
  return (list("nodes" = nodes_out,
               "GOTable" = GOTable,
               "topGOterms" = TOP_GO_TERMS))
  
}

##### CALCULATE ACTIVITY SPELLS ####

calculateActivitySpells <- function(nodes_in, edges_in, timepoints, interval, weight_threshold,
                                    bait_ids) {
  
  nodes_out <- nodes_in
  edges_out <- edges_in
  
  edges_out$type <- "Provided"
  
  # calculate edge onset and terminus information based on threshold
  # onset when weight >= threshold, terminus when weight <= threshold
  # duplicate link if there are non-overlapping spells of activity (ex: tp1-2, tp4-5)
  
  NUM_EDGES <- dim(edges_out)[1]
  
  f <- 0
  freq <- vector(mode='numeric', length=NUM_EDGES) 
  # count number of times edge must repeat: non-overlapping spells
  onset <- {}
  terminus <- {} 
  
  NUM_TPS <- length(timepoints)
  temp_timepoints <- c(timepoints, timepoints[NUM_TPS]+interval)
  for (n in (1:NUM_EDGES)) {
    j <- 1
    freq[n] <- 0
    while (j <= NUM_TPS) {
      count <- 0
      time <- temp_timepoints[j]
      w <- edges_out[n, paste('w', time, sep="")]
      if (w >= weight_threshold) {
        f <- f+1
        freq[n] <- freq[n] + 1
        onset[f] <- time
        while (w >= weight_threshold) {
          count <- count+1
          j <- j+1
          time <- temp_timepoints[j]
          if (time %in% timepoints) {
            w <- edges_out[n, paste('w', time, sep="")] 
          } else { w <- -1 } # set a negative weight to break out of while loop
        }
        terminus[f] <- time
      }
      j <- j+1
      time <- temp_timepoints[j]
    }
  }
  
  edges_out <- edges_out[rep(rownames(edges_out), freq), ]
  edges_out$onset <- onset
  edges_out$terminus <- terminus
  rownames(edges_out) <- 1:nrow(edges_out)
  edges_out <- edges_out[order(edges_out$prey_id), ]
  
  ## ASSIGN NODE SPELLS AND SHARED NODES
  
  # assign unique node onset and terminus times from edge info such that node activity encompasses 
  # incident edge activity - examples:
  # 1. edge tps=(1-2, 2-5) -> node tps=(1-5)
  # 2. edge tps=(1-5, 2-4) -> node tps=(1-5)
  # repeat until as condense as possible, then remove any remaining duplicates (on-off periods of edges) by
  # extending node activity for entire edge spell
  
  nodes_out <- nodes_in
  
  nodes_out <- nodes_out[!duplicated(nodes_out$id), ]
  nodes_out <- nodes_out[order(nodes_out$id), ]
  rownames(nodes_out) <- nodes_out$id
  
  # get the number of times to repeat each node (if it appears in multiple baits)
  freq <- as.data.frame(table(edges_out$prey_id), stringsAsFactors = F)
  keep_nodes <- nodes_out[rep(freq$Var1, freq$Freq), ]
  keep_nodes$onset <- edges_out$onset
  keep_nodes$terminus <- edges_out$terminus
  
  current_dim <- dim(keep_nodes)[1]
  row.names(keep_nodes) <- 1:current_dim
  new_dim <- 0
  
  count <- 1
  while (new_dim < current_dim) {
    count <- count+1
    remove_ids <- numeric()
    i <- 1
    while (i < dim(keep_nodes)[1]) {
      currentNode <- keep_nodes[i,"id"]
      while (i < dim(keep_nodes)[1] & keep_nodes[i+1, "id"] == currentNode) {
        # case 1
        if (keep_nodes[i, "onset"] == keep_nodes[i+1, "terminus"] | 
            keep_nodes[i, "terminus"] == keep_nodes[i+1, "onset"]) {
          keep_nodes[i, "onset"] <- min(keep_nodes[i, "onset"], keep_nodes[i+1, "onset"])
          keep_nodes[i, "terminus"] <- max(keep_nodes[i, "terminus"], 
                                           keep_nodes[i+1, "terminus"])
          remove_ids <- append(remove_ids, i+1)
        }
        # case 2
        else if ((keep_nodes[i, "onset"] <= keep_nodes[i+1, "onset"] & 
                  keep_nodes[i, "terminus"] >= keep_nodes[i+1, "terminus"]) |
                 (keep_nodes[i, "onset"] >= keep_nodes[i+1, "onset"] & 
                  keep_nodes[i, "terminus"] <= keep_nodes[i+1, "terminus"])) {
          keep_nodes[i, "onset"] <- min(keep_nodes[i, "onset"], keep_nodes[i+1, "onset"])
          keep_nodes[i, "terminus"] <- max(keep_nodes[i, "terminus"], 
                                           keep_nodes[i+1, "terminus"])
          remove_ids <- append(remove_ids, i+1)
        }
        i <- i+1
      }
      i <- i+1
    }
    current_dim <- dim(keep_nodes)[1]
    if (length(remove_ids) == 0) { break }
    keep_nodes <- keep_nodes[-remove_ids, ]
    new_dim <- dim(keep_nodes)[1]
    rownames(keep_nodes) <- 1:new_dim
  }
  
  # remove remaining duplicates by extending them for max possible time when edges active
  remove_ids <- numeric()
  duplicates <- names(which(table(keep_nodes$id) > 1))
  if (length(duplicates) > 0) {
    for (i in 1:length(duplicates)) {
      check <- which(keep_nodes$id == duplicates[i])
      onset <- min(keep_nodes[check, "onset"])
      terminus <- max(keep_nodes[check, "terminus"])
      keep_nodes[check[1], "onset"] <- onset
      keep_nodes[check[1], "terminus"] <- terminus
      remove_ids <- append(remove_ids, check[2:length(check)])
    }
    keep_nodes <- keep_nodes[-remove_ids, ]
  }
  
  # ensure that bait nodes are always present
  NUM_TPS <- length(timepoints)
  missing_baits <- bait_ids[!(bait_ids %in% keep_nodes$id)]
  if (length(missing_baits) > 0) {
    add_baits <- nodes_out[missing_baits, ]
    add_baits$onset <- timepoints[1]
    add_baits$terminus <- timepoints[NUM_TPS]
    
    keep_nodes <- rbind(keep_nodes, add_baits)
  }
  
  nodes_out <- keep_nodes 
  nodes_out$onset[which(nodes_out$id %in% bait_ids)] <- timepoints[1]
  nodes_out$terminus[which(nodes_out$id %in% bait_ids)] <- timepoints[NUM_TPS] + interval
  
  nodes_out$duration <- nodes_out$terminus - nodes_out$onset
  nodes_out <- nodes_out[order(nodes_out$id), ]
  nodes_out <- na.omit(nodes_out)
  
  NUM_NODES <- dim(nodes_out)[1]
  
  # convert old node ids to new ids
  old_to_new <- setNames(as.list(1:NUM_NODES), as.character(nodes_out$id))
  nodes_out$id <- 1:NUM_NODES
  
  edges_out$bait_id <- unlist(old_to_new[as.character(edges_out$bait_id)])
  edges_out$prey_id <- unlist(old_to_new[as.character(edges_out$prey_id)])
  
  bait_ids <- unlist(old_to_new[as.character(bait_ids)])
  
  # function returns: nodes, edges, bait_ids
  return (list("nodes" = nodes_out,
               "edges" = edges_out,
               "bait_ids" = bait_ids))
} 











##### GET STRING-DB INFORMATION ####

getSTRINGData <- function(nodes_in, edges_in, taxids, timepoints) {
  
  edges_out <- edges_in
  
  source_python("AppFiles/getStringInteractors.py")
  source_python("AppFiles/gsi.py")
  
  string_links <- {}
  
  for (species in taxids) {
    values <- nodes_in[which(nodes_in$taxid == species), "accession"]
    keys <- nodes_in[which(nodes_in$taxid == species), "gene_name"]
    node_list <- setNames(as.list(values), keys)
    interactions <- run_functions(dict(node_list), species)
    string_links <- rbind(string_links, interactions)
  }
  
  colnames(string_links) <- c("bait_accession", "prey_accession",
                              "bait_gene_name", "prey_gene_name", "weight")
  
  # only keep interactions that have been experimentally verified 
  string_links <- string_links[which(string_links$weight != 0), ]
  
  NUM_STRING_EDGES <- dim(string_links)[1]
  
  bait_id <- vector(mode = 'numeric', length = NUM_STRING_EDGES)
  prey_id <- vector(mode = 'numeric', length = NUM_STRING_EDGES)
  bait_gene_name <- vector(mode = 'character', length = NUM_STRING_EDGES)
  prey_gene_name <- vector(mode = 'character', length = NUM_STRING_EDGES)
  
  for (i in 1:NUM_STRING_EDGES) {
    bait_id[i] <- nodes_in[which(nodes_in$accession == string_links[i, "bait_accession"]), "id"]
    prey_id[i] <- nodes_in[which(nodes_in$accession == string_links[i, "prey_accession"]), "id"]
    bait_gene_name[i] <- nodes_in[which(nodes_in$accession == string_links[i, "bait_accession"]), 
                                  "gene_name"][1]
    prey_gene_name[i] <- nodes_in[which(nodes_in$accession == string_links[i, "prey_accession"]), 
                                  "gene_name"][1]
  }
  
  string_links$bait_id <- bait_id
  string_links$prey_id <- prey_id
  string_links$bait_gene_name <- bait_gene_name
  string_links$prey_gene_name <- prey_gene_name
  
  # give edge same weight and NA abundance at all time points
  for (tp in timepoints) {
    weight_var <- string_links$weight
    string_links <- cbind(string_links, weight_var)
    colnames(string_links)[length(string_links)] <- paste('w', tp, sep="")
    abundance_var <- rep(NA, dim(string_links)[1])
    string_links <- cbind(string_links, abundance_var)
    colnames(string_links)[length(string_links)] <- paste('a', tp, sep="")
    norm_abundance_var <- rep(NA, dim(string_links)[1])
    string_links <- cbind(string_links, norm_abundance_var)
    colnames(string_links)[length(string_links)] <- paste('norm_a', tp, sep="")
  }
  
  # make onset and terminus infinity
  string_links$onset <- Inf
  string_links$terminus <- Inf
  
  # assign type to be STRING inferred or both
  both <- merge(string_links, edges_out, by=c("bait_id", "prey_id"))
  if (dim(both)[1] > 0) {
    both$type <- "Both"
    string_links <- merge(string_links, both, by=c("bait_id", "prey_id"), all=T)
    string_links$type[is.na(string_links$type)] <- "STRING"
  } else {
    string_links$type <- "STRING"
  }
  
  # assign clusters to be NA 
  string_links$cluster <- NA
  
  # bind STRING links to provided links file
  string_links <- string_links[, colnames(edges_out)]
  edges_out <- rbind.data.frame(edges_out, string_links)
  
  # return: edges
  return (edges_out)
}

##### BUILD NETWORK ####

# only do this if there are fewer than 500 nodes in the data
buildNetwork <- function(nodes_in, edges_in, num_edges, abund_prov, timepoints, interval,
                         bait_ids) {
  
  library(networkDynamic)
  library(RColorBrewer)
  library(ndtv)
  library(tsna)
  
  edges_out <- edges_in
  nodes_out <- nodes_in
  
  # rearrange link dataframes to network
  if (abund_prov) {
    edges_out <- edges_out[, c("bait_id", "prey_id", "bait_accession", "prey_accession", 
                               "bait_gene_name", "prey_gene_name", "cluster", "onset", 
                               "terminus", "type", 
                               paste("w", timepoints, sep=""),
                               paste("a", timepoints, sep=""))]
  } else {
    edges_out <- edges_out[, c("bait_id", "prey_id", "bait_accession", "prey_accession", 
                               "bait_gene_name", "prey_gene_name", "onset", "terminus", "type", 
                               paste("w", timepoints, sep=""))]
  }
  
  colnames(edges_out)[c(1,2)] <- c("from", "to")
  
  nodes_out$localization_string <- unlist(nodes_out$localization_string)
  
  # give localization, duration and edge type a factor numeric value (to later assign color)
  nodes_out$locanum <- as.numeric(factor(nodes_out$localization_string))
  # give duration factor based on what quartile of expression it is in
  nodes_out$durfactor <- cut(nodes_out$duration, unique(quantile(nodes_out$duration, probs=0:4/4)),
                             include.lowest=T, labels=F)
  edges_out$edge_color <- c("black", "navyblue", "gray75")[as.numeric(factor(edges_out$type))]
  
  # create network object
  
  net3 <- network(edges_out, vertex.attr=nodes_out, matrix.type="edgelist", 
                  loops=F, multiple=F, directed=F, ignore.eval = F)
  
  # assign colors to localization, duration factors
  newPalette <- colorRampPalette(brewer.pal(8, "Set2"))
  loc_color <- newPalette(n = nlevels(factor(nodes_out$localization_string))+1)
  
  dur_color <- brewer.pal(n = nlevels(factor(nodes_out$durfactor)), name = "Reds")
  
  net3 %v% "loc_color" <- loc_color[net3 %v% "locanum"]
  net3 %v% "dur_color" <- dur_color[net3 %v% "durfactor"]
  
  # dynamic properties
  NUM_NODES <- dim(nodes_out)[1]
  
  vs <- data.frame(onset=nodes_out$onset, terminus=nodes_out$terminus, vertex.id=1:NUM_NODES)
  es <- data.frame(onset=edges_out$onset, terminus=edges_out$terminus, tail=edges_out$to, 
                   head=edges_out$from)
  
  # create dynamic network object
  ## NET3.DYN HAS ALL THE ORIGINAL DATA, NEVER MODIFY THIS BUT ONLY MAKE COPIES TO MODIFY
  
  net3.dyn <- networkDynamic(base.net=net3, vertex.spells=vs, edge.spells=es)
  
  net3.dyn <- reconcile.edge.activity(net3.dyn, mode="reduce.to.vertices")
  
  # set dynamic edge weights
  NUM_TPS <- length(timepoints)
  NUM_EDGES <- dim(edges_out)[1]
  
  for (i in 1:NUM_TPS) {
    temp_timepoints <- c(timepoints, timepoints[NUM_TPS]+interval)
    tp <- temp_timepoints[i]
    next_tp <- temp_timepoints[i+1]
    net3.dyn <- activate.edge.attribute(net3.dyn, 'edge_weight', 
                                        (edges_out[, paste('w', tp, sep="")])*2, 
                                        onset=tp, terminus=next_tp, e=1:NUM_EDGES)
    # set dynamic vertex sizes (if single bait)
    if (length(bait_ids) == 1 & abund_prov) {
      net3.dyn <- activate.vertex.attribute(net3.dyn, 'node_size',
                                            nodes_out[, paste('size_', tp, sep="")],
                                            onset=tp, terminus=next_tp, v=1:NUM_NODES)
    }
  }
  
  # function returns: nodes, edges, dynamic net3 object
  return (list("nodes" = nodes_out,
               "edges" = edges_out,
               "net3.dyn" = net3.dyn))
}

# output protein information
saveNodeAttributes <- function(nodes_in, bait_ids, timepoints, filename) {
  output_nodes <- nodes_in[, c("gene_name", "accession", "organism", 
                               "localization_string_output", "complexes", 
                               "GO_term_string")]
  output_nodes$GO_term_string <- sub("^; ", "", output_nodes$GO_term_string)
  output_nodes$localization_string_output <- sub("^; ", "",
                                                 output_nodes$localization_string_output)
  output_nodes$organism <- unlist(output_nodes$organism)
  colnames(output_nodes) <- c("gene_name", "accession", "species", "localizations",
                              "complexes", "GO_terms")

  NUM_BAITS <- length(bait_ids)

  if (NUM_BAITS > 1) {
    output_nodes <- cbind(output_nodes,
                          nodes_in[, paste("shared_baits_", timepoints, sep="")])
  }
  
  # replace NA values with empty string
  output_nodes[is.na(output_nodes)] <- ""
  
  write.table(output_nodes, file = filename, sep = "\t", 
              quote = F, row.names = F)
}

# output edge information
saveEdgeAttributes <- function(edges_in, timepoints, abund_prov, norm_prot, filename) {
  
  # get edge information to output
  output_edges <- edges_in[, c("bait_gene_name", "prey_gene_name",
                               "bait_accession", "prey_accession",
                               "onset", "terminus", "type", 
                               paste("w", timepoints, sep=""))]
  output_edges_colnames <- c("#node1", "node2", "node1_accession", "node2_accession", 
                             "onset", "terminus", "type", 
                             paste("confidence_score_", timepoints, sep=""))
  
  if (abund_prov) {
    output_edges <- cbind(output_edges, edges_in[, c("cluster", 
                                                     paste("a", timepoints, sep=""))])
    output_edges_colnames <- c(output_edges_colnames, "cluster",
                               paste("abundance_", timepoints, sep=""))
    if (norm_prot) {
      output_edges <- cbind(output_edges, edges_in[, paste("norm_a", timepoints, sep="")])
      output_edges_colnames <- c(output_edges_colnames, 
                                 paste("normalized_to_proteome_abundance_", 
                                       timepoints, 
                                       sep=""))
    }
  }
  
  output_edges[is.na(output_edges)] <- ""
  
  # rename abundance columns
  colnames(output_edges) <- output_edges_colnames
  
  write.table(output_edges, file = filename, sep = "\t", 
              quote = F, row.names = F)
}

##### MODIFY NETWORK BASED ON INPUTS ####

# only run these functions if network exists

deactivateNE <- function(min_shared, st, locns, species, go_terms, net3.dyn, nodes_in, 
                         edges_in, timepoints, interval, bait_ids) {
  
  # calculate shared nodes
  plot_net <- net3.dyn
  # deactivate nodes based on number of baits they share
  for (tp in timepoints) {
    keep_nodes <- nodes_in[nodes_in[, paste("shared_", tp, sep="")] >= min_shared, "id"]
    keep_nodes <- c(keep_nodes, bait_ids)
    keep_nodes <- unique(keep_nodes)
    deactivate_nodes <- nodes_in[!(nodes_in$id %in% keep_nodes), "id"]
    deactivate.vertices(plot_net, onset=tp, terminus=tp+interval, v=deactivate_nodes, 
                        deactivate.edges=T)
  }
  
  
  # deactivate edges (based on STRING confidence threshold)
  deactivate_link_ids <- which((edges_in$type == "STRING") &
                                 (edges_in[, paste("w", timepoints[1], sep="")] <= st))
  
  if (!is.null(deactivate_link_ids)) {
    deactivate.edges(plot_net, e=deactivate_link_ids)
  }
  
  # deactivate nodes (based on localizations, species, functions, etc.)
  
  # localizations
  keep_loc_vertex_ids <- {}
  for (loc in locns) {
    keep_loc_vertex_ids <- append(keep_loc_vertex_ids, 
                                  nodes_in[which(nodes_in[, paste("loc_", loc, sep="")] == 1), "id"])
  }
  
  # species 
  keep_spec_vertex_ids <- {}
  for (spec in species) {
    keep_spec_vertex_ids <- append(keep_spec_vertex_ids, 
                                   nodes_in[which(nodes_in$organism == spec), "id"])
  }
  
  # GO_terms
  if (!is.null(go_terms)) {
    
    keep_got_vertex_ids <- {}
    if ("All" %in% go_terms) {
      keep_got_vertex_ids <- nodes_in$id
    } else {
      for (got in go_terms) {
        keep_got_vertex_ids <- append(keep_got_vertex_ids, 
                                      nodes_in[which(nodes_in[, paste("got_", got, sep="")] == 1), "id"])
      }
    }
    
    keep_vertex_ids <- Reduce(intersect, list(keep_loc_vertex_ids, keep_spec_vertex_ids,
                                              keep_got_vertex_ids))
  } else {
    
    # there are no GOterms
    keep_vertex_ids <- Reduce(intersect, list(keep_loc_vertex_ids, keep_spec_vertex_ids))
    
  }
  
  
  deactivate_vertex_ids <- setdiff(nodes_in$id, keep_vertex_ids)
  
  NUM_TPS <- length(timepoints)
  
  if (!is.null(deactivate_vertex_ids)) {
    # don't deactivate baits
    deactivate_vertex_ids <- deactivate_vertex_ids[!deactivate_vertex_ids %in% bait_ids]
    
    deactivate.vertices(plot_net, onset=timepoints[1], terminus=(timepoints[NUM_TPS]+interval),
                        v=deactivate_vertex_ids, deactivate.edges=T)
  }
  
  # compute network animation
  compute.animation(plot_net, animation.mode = "kamadakawai", default.dist=2, 
                    slice.par=list(start=timepoints[1], end=timepoints[NUM_TPS], 
                                   interval=interval, 
                                   aggregate.dur=interval, rule='any'))
  return (plot_net)
}

##### PLOT DYNAMIC NETWORK ####

# only run this function if network exists

plotCompNet <- function(net3_network, lab, bait_ids, abund_prov) {
  
  if(length(bait_ids) == 1 & abund_prov) {
    # plot network animation
    render.d3movie(net3_network,
                   render.par=list(tween.frames = 10, show.time = F),
                   plot.par=list(bg='white'),
                   d3.options = list(animationDuration=6000,enterExitAnimationFactor=0.3),
                   usearrows = F,
                   displaylabels = lab, label=function(slice){slice%v%'gene_name'},
                   vertex.border=function(slice){slice%v%'dur_color'},
                   vertex.cex = function(slice){slice%v%'node_size'},
                   vertex.lwd = 3,
                   vertex.col = function(slice){slice%v%'loc_color'},
                   edge.lwd = 'edge_weight', 
                   edge.col = function(slice){slice%e%'edge_color'},
                   vertex.tooltip = function(slice){paste("<b>Gene Name:</b>", slice%v%'gene_name', 
                                                          "<br>", "<b>Localization:</b>", 
                                                          slice%v%'localization_string',
                                                          "<br>", "<b>Complexes:</b>", 
                                                          slice%v%'complexes',
                                                          "<br>", "<b>Time Active:</b>", 
                                                          slice%v%'duration')},
                   edge.tooltip = function(slice){paste('From:',slice%e%'bait_gene_name', ' To: ', 
                                                        slice%e%'prey_gene_name', '<br>',
                                                        'Weight:', (slice%e%'edge_weight')/2)},
                   output.mode = 'htmlWidget')
  } else {
    
    # plot network animation
    render.d3movie(net3_network,
                   render.par=list(tween.frames = 10, show.time = F),
                   plot.par=list(bg='white'),
                   d3.options = list(animationDuration=6000,enterExitAnimationFactor=0.3),
                   usearrows = F,
                   displaylabels = lab, label=function(slice){slice%v%'gene_name'},
                   vertex.border=function(slice){slice%v%'dur_color'},
                   vertex.cex = 1.5,
                   vertex.lwd = 3,
                   vertex.col = function(slice){slice%v%'loc_color'},
                   edge.lwd = 'edge_weight', 
                   edge.col = function(slice){slice%e%'edge_color'},
                   vertex.tooltip = function(slice){paste("<b>Gene Name:</b>", slice%v%'gene_name', 
                                                          "<br>", "<b>Localization:</b>", 
                                                          slice%v%'localization_string',
                                                          "<br>", "<b>Complexes:</b>", 
                                                          slice%v%'complexes',
                                                          "<br>", "<b>Time Active:</b>", 
                                                          slice%v%'duration')},
                   edge.tooltip = function(slice){paste('From:',slice%e%'bait_gene_name', ' To: ', 
                                                        slice%e%'prey_gene_name', '<br>',
                                                        'Weight:', (slice%e%'edge_weight')/2)},
                   output.mode = 'htmlWidget')
    
  }
}

##### ANALYTICAL TOOLS REPORT ####

analyticalTools <- function(net3.dyn, nodes_in, edges_in, timepoints, localizations, topGOterms) {
  
  # timeline of node activity separated by localization
  for (loc in localizations) {
    vs <- nodes_in[which(nodes_in[, paste("loc_", loc, sep="")] == 1), "id"]
    print(timeline(net3.dyn, v=vs, plot.edge.spells=F, displaylabels = F, 
                   xlab = "Time", ylab = paste(loc, "Proteins")))
  }
  
  # timeline of node activity separated by GO terms
  if (!is.null(topGOterms)) {
    for (got in topGOterms) {
      vs <- nodes_in[which(nodes_in[, paste("got_", got, sep="")] == 1), "id"]
      print(timeline(net3.dyn, v=vs, plot.edge.spells=F, displaylabels = F, 
                     xlab = "Time", ylab = paste(got, "proteins")))
    }
  }
  
  # edge formation
  data <- edges_in[edges_in$type == "Provided", "onset"]
  
  print(ggplot(as.data.frame(data), aes(x = data)) +
          geom_histogram() +
          scale_x_continuous(breaks = timepoints) +
          labs(x = "Time", y = "Number of Edges Formed", title = "Time of Edge Formation"))
  
  # edge duration
  provided_edges <- which(edges_in$type == "Provided")
  data <- data.frame(x=edgeDuration(net3.dyn, e=provided_edges))
  print(ggplot(data) +
          geom_histogram(mapping = aes(x = data$x)) +
          scale_x_continuous(breaks = timepoints) +
          labs(x = "Duration", y = "Number of Edges", title = "Duration of Edge Activity"))
}

##### ABUNDANCE PLOTS #####

# only perform this function if abundances provided
cleanAbundanceData <- function(edges_in, timepoints) {
  
  # re-orient data frame to make it suitable for ggplot2 plotting
  abund_info <- edges_in[which(edges_in$type == "Provided"), ]
  abund_info <- abund_info[, c("bait_gene_name", "prey_gene_name", 
                               paste("a", timepoints, sep=""))]
  
  # melt abundances dataframe down
  abundances <- {}
  NUM_TPS <- length(timepoints)
  
  for (i in 1:dim(abund_info)[1]) {
    new <- cbind.data.frame(timepoints, rep(abund_info[i, "bait_gene_name"], NUM_TPS))
    new <- cbind.data.frame(new, rep(abund_info[i, "prey_gene_name"], NUM_TPS))
    new <- cbind.data.frame(new, as.vector(t(abund_info[i, paste("a", timepoints, sep="")])))
    abundances <- rbind(abundances, new)
  }
  colnames(abundances) <- c("time", "bait_gene_name", "prey_gene_name", "abundance")
  rownames(abundances) <- 1:dim(abundances)[1]
  
  return (abundances)
}

# create list of neighbors for every node
calcNeighbors <- function(confidence, locns, species, go_terms, 
                          nodes_in, edges_in, timepoints, localizations) {
  
  # localizations
  keep_loc_nodes <- {}
  if (all(locns == localizations)) { 
    keep_loc_nodes <- nodes_in$accession 
  } else {
    for (loc in locns) {
      keep_loc_nodes <- append(keep_loc_nodes, 
                               nodes_in[which(nodes_in[, paste("loc_", loc, sep="")] == 1), 
                                        "accession"])
    }
  }
  
  # species 
  keep_spec_nodes <- {}
  for (spec in species) {
    keep_spec_nodes <- append(keep_spec_nodes, 
                              nodes_in[which(nodes_in$organism == spec), "accession"])
  }
  
  # GO_terms
  
  if (!is.null(go_terms)) {
    
    keep_got_nodes <- {}
    if ("All" %in% go_terms) {
      keep_got_nodes <- nodes_in$accession
    } else {
      for (got in go_terms) {
        keep_got_nodes <- append(keep_got_nodes, 
                                 nodes_in[which(nodes_in[, paste("got_", got, sep="")] == 1), 
                                          "accession"])
      }
    }
    
    keep_nodes <- Reduce(intersect, list(keep_loc_nodes, keep_spec_nodes,
                                         keep_got_nodes))
    
  } else {
    
    keep_nodes <- Reduce(intersect, list(keep_loc_nodes, keep_spec_nodes))
    
  }
  
  
  keep_links <- edges_in[which(edges_in$bait_accession %in% keep_nodes &
                                 edges_in$prey_accession %in% keep_nodes &
                                 edges_in[, paste('w', timepoints[1], sep="")] >= confidence &
                                 edges_in$type == "STRING"), ]
  neighbors <- list()
  NUM_NODES <- dim(nodes_in)[1]
  
  for (i in 1:NUM_NODES) {
    add <- keep_links[which(keep_links$bait_accession == nodes_in[i, "accession"]), 
                      "prey_accession"]
    add <- c(add, keep_links[which(keep_links$prey_accession == nodes_in[i, "accession"]), 
                             "bait_accession"])
    neighbors[[nodes_in[i, "accession"]]] <- add
  }
  return (neighbors)
}

# plot abundances
abundPlot <- function(genesyms, nhbr=F, neighborList,
                      nodes_in, plotabundances, timepoints) {
  
  # decimal places to show in label
  scalefunction <- function(x) sprintf("%.2f", x)
  
  genes <- unlist(strsplit(genesyms, "; "))
  # display neighboring nodes abundances as well
  if (nhbr == T) {
    accessions <- {}
    for (gene in genes) {
      accessions <- c(accessions, 
                      nodes_in[grepl(gene, nodes_in$gene_name, ignore.case = T), "accession"])
    }
    neighbor_genes <- nodes_in[nodes_in$accession %in% unlist(neighborList[accessions]), 
                               "gene_name"]
    genes <- as.vector(c(genes, neighbor_genes))
  }
  current_genes <- {}
  for (gene in genes) {
    current_genes <- rbind(current_genes, 
                           plotabundances[grepl(gene, plotabundances$prey_gene_name, 
                                                ignore.case = T), ])
  }
  print(ggplot(data = current_genes) +
          geom_line(mapping = aes(x=time, y=abundance, color=prey_gene_name, 
                                  linetype=bait_gene_name), size=1.5) +
          scale_y_continuous(trans='log', labels=scalefunction) +
          scale_x_continuous(breaks = timepoints) + 
          labs(x="Time", y="Abundance", color="Search Genes", 
               linetype="Bait"))
}

printClusterNumber <- function(genesyms, nodes_in, edges_in, bait_ids) {
  genes <- unlist(strsplit(genesyms, ", "))
  for (gene in genes) {
    for (bait in bait_ids) {
      bait_gene_name <- nodes_in$gene_name[which(nodes_in$id == bait)]
      clusterNumber <- edges_in[which(edges_in$from == bait & 
                                        grepl(gene, edges_in$prey_gene_name, ignore.case = T)),
                                "cluster"]
      if (length(clusterNumber) != 0) {
        out <- cat(paste(gene, " is in cluster ", clusterNumber, " for bait ", 
                         bait_gene_name, "\n", collapse=""))
      }
    }
  }
}

##### BUILD UI ####

## Code to display generic error message when there is an error ----

options(shiny.sanitize.errors = TRUE)

ui = fluidPage(
  
  # Enable shinyjs
  shinyjs::useShinyjs(),
  
  # Aesthetics ----
  tags$head(
    tags$style(HTML("
                    hr {border-top: 4px solid #000000;}
                    .shiny-notification {
                    height: 100px;
                    width: 800px;
                    position:fixed;
                    top: calc(50% - 50px);;
                    left: calc(50% - 400px);;
                    }
                    "))
  ),
  
  # Create tabs ----
  tabsetPanel(
    id = "navbar",
    
    # source each tab UI from separate file
    source(file.path("ui", "input.R"), local = TRUE)$value,
    source(file.path("ui", "interactome.R"), local = TRUE)$value,
    source(file.path("ui", "quantplots.R"), local = TRUE)$value,
    source(file.path("ui", "GOterms.R"), local = TRUE)$value,
    source(file.path("ui", "complexes.R"), local = TRUE)$value,
    source(file.path("ui", "downloads.R"), local = TRUE)$value
    
  )
)

##### DEFINE SERVER LOGIC ####

# Define server logic to read selected file
server <- function(input, output, session) {
  
  # source each tab server from separate file
  source(file.path("server", "input.R"), local = TRUE)$value
  source(file.path("server", "programfunctions.R"), local = TRUE)$value
  source(file.path("server", "interactome.R"), local = TRUE)$value
  source(file.path("server", "quantplots.R"), local = TRUE)$value
  source(file.path("server", "GOterms.R"), local = TRUE)$value
  source(file.path("server", "complexes.R"), local = TRUE)$value
  source(file.path("server", "downloads.R"), local = TRUE)$value
  
  # Stop App ----
  # stop the app when the browser window is closed
  
  session$onSessionEnded(stopApp)
}

##### RUN APP ####

shinyApp(ui, server)