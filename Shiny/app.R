library(shiny)
library(reticulate)
library(shinyBS)
library(shinyjs)
library(zip)
library(DT)
library(ggplot2)
library(gridExtra)
library(topGO)
library(networkDynamic)
library(RColorBrewer)
library(tsna)


# make a virtualenv
reticulate::virtualenv_create(envname = 'intervista', python = '/usr/bin/python3')

# install Python packages
reticulate::virtualenv_install(envname= 'intervista', packages=c('requests', 'pandas', 'biopython==1.76'))

# use intervista environment
reticulate::use_virtualenv("intervista", required=TRUE)

source_python("AppFiles/getStringInteractors.py")
source_python("AppFiles/gsi.py")
source_python("AppFiles/UniprotQuery.py")
source_python("AppFiles/ruq.py")

library(ndtv)

##### GET AND CLEAN DATA ####

cleanData <- function(nodes_in, edges_in, timepoints, spec_threshold) {
  
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
  
  # get timepoint information 
  NUM_TPS <- length(timepoints)
  
  # save original timepoints for network rendering
  named_timepoints <- timepoints
  names(named_timepoints) <- 1:NUM_TPS
  
  # if timepoints aren't numeric, assign them integers
  if (sum(grepl("[A-Za-z]", timepoints)) != 0) {
    timepoints <- 1:NUM_TPS
  } else {
    timepoints <- as.numeric(timepoints)
  }
  
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
  num_link_attr <- ncol(edges_out)
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
  edge_accessions <- unique(c(edges_out$bait_accession, edges_out$prey_accession))
  nodes_out <- nodes_out[nodes_out$accession %in% edge_accessions, ]
  NUM_NODES <- dim(nodes_out)[1]
  
  edges_out <- edges_out[!duplicated(edges_out[, c("bait_accession", "prey_accession")]), ]
  NUM_EDGES <- dim(edges_out)[1]
  
  # assign unique ids to nodes
  nodes_out$id <- 1:NUM_NODES
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
  
  
  # return list of: nodes, edges, weight threshold, timepoints
  # loc_prov, abund_prov, bait IDs, taxids
  
  cleaned_data <- list("nodes" = nodes_out,
                       "edges" = edges_out,
                       "num_edges" = NUM_EDGES,
                       "weight_threshold" = WEIGHT_THRESHOLD,
                       "timepoints" = timepoints,
                       "loc_prov" = LOC_PROV,
                       "abund_prov" = ABUND_PROV,
                       "bait_ids" = BAIT_IDS,
                       "taxids" = TAXIDS,
                       "named_timepoints" = named_timepoints
  )
  return (cleaned_data)
}

cleanDataPREV <- function(nodes_in, edges_in) {
  
  nodes_out <- nodes_in
  edges_out <- edges_in
  
  # rename edge columns
  colnames(edges_out)[1:4] <- c("bait_gene_name", "prey_gene_name", 
                                "bait_accession", "prey_accession")
  colnames(edges_out)[colnames(edges_out) %in% c("onset_condition", "terminus_condition")] <- 
    c("orig_onset", "orig_terminus")
  
  # get number of edges (provided, no repeats) for building interactome network
  provided_edges <- edges_out[edges_out$type == "Provided", c("bait_accession", "prey_accession")]
  NUM_EDGES <- dim(unique(provided_edges))[1]
  
  # get named timepoints
  # get column names that will have timepoints
  named_timepoints <- colnames(edges_out)[grep("^confidence", colnames(edges_out))]
  named_timepoints <- sub(".*_", "", named_timepoints)
  
  # create timepoints
  NUM_TPS <- length(named_timepoints)
  timepoints <- 1:NUM_TPS
  names(named_timepoints) <- timepoints
  named_timepoints <- as.list(named_timepoints)
  
  # rename edge columns
  colnames(edges_out)[grep("^confidence", colnames(edges_out))] <- paste("w", timepoints, sep="")
  
  # check if abundances were provided
  ABUND_PROV <- F
  if ( sum(grepl("^abundance", colnames(edges_out))) != 0 ) { ABUND_PROV <- T }
  
  if (ABUND_PROV) {
    # rename edge columns
    colnames(edges_out)[grep("^abundance", colnames(edges_out))] <- paste("a", timepoints, sep="")
    
  }
  
  # check if abundances were normalized
  NORM_PROT <- F
  if ( sum(grepl("^normalized", colnames(edges_out))) != 0 ) { NORM_PROT <- T }
  
  if (NORM_PROT) {
    # rename edge columns
    colnames(edges_out)[grep("^normalized", colnames(edges_out))] <- paste("norm_a", 
                                                                           timepoints, sep="")
    
  }
  
  # assign unique ids to nodes
  NUM_NODES <- nrow(nodes_out)
  nodes_out$id <- c(1:NUM_NODES)
  rownames(nodes_out) <- nodes_out$id
  
  # get shared nodes (if they exist)
  if ( sum(grepl("^shared", colnames(nodes_out))) != 0 ) {
    i <- 1
    for (tp in named_timepoints) {
      nodes_out <- cbind(nodes_out, lengths(strsplit(nodes_out[, paste("shared_baits_", tp, sep="")], 
                                                     ";")))
      colnames(nodes_out)[ncol(nodes_out)] <- paste("shared_", timepoints[i], sep="")
      nodes_out[, ncol(nodes_out)] <- as.character(nodes_out[, ncol(nodes_out)])
      i <- i+1
    }
  }
  
  # get list of localizations
  LOCS <- unique(unlist(strsplit(nodes_out$localizations, "; ")))
  
  # assign localization string (replace multiple localizations with string "Multiple")				
  nodes_out$localization_string <- nodes_out$localizations				
  nodes_out[grep(";", nodes_out$localization_string), "localization_string"] <- "Multiple localizations"				
  nodes_out[nodes_out$localization_string == "", "localization_string"] <- "Unspecified"
  
  # create loc_ columns with T/F for each localization
  for (loc in LOCS) {
    check <- vector(mode="logical", length=NUM_NODES)
    check <- grepl(loc, nodes_out$localizations, ignore.case=T)
    nodes_out <- cbind(nodes_out, check)
    colnames(nodes_out)[ncol(nodes_out)] <- paste("loc_", loc, sep="")
  }
  
  # get top GO terms
  topGOterms <- unique(unlist(strsplit(nodes_out$GO_terms, "; ")))
  
  # create topGO columns with T/F for each GOterm
  for (got in topGOterms) {
    check <- vector(mode="logical", length=NUM_NODES)
    check <- grepl(got, nodes_out$GO_terms, ignore.case=T)
    nodes_out <- cbind(nodes_out, check)
    colnames(nodes_out)[ncol(nodes_out)] <- paste("got_", got, sep="")
  }
  
  # get organism names
  colnames(nodes_out)[which(colnames(nodes_out) == "species")] <- "organism"
  ORGANISMS <- unique(unlist(strsplit(nodes_out$organism, "; ")))
  
  # attach ids and gene names to links
  edges_out <- merge(edges_out, nodes_out[, c("accession", "id")], 
                     by.x="bait_accession", 
                     by.y="accession")
  colnames(edges_out)[ncol(edges_out)] <- "bait_id"
  edges_out <- merge(edges_out, nodes_out[, c("accession", "id")], 
                     by.x="prey_accession", 
                     by.y="accession")
  colnames(edges_out)[ncol(edges_out)] <- "prey_id"
  
  # get node durations
  nodes_out$duration <- nodes_out$terminus - nodes_out$onset
  
  # get provided bait ids
  BAIT_IDS <- unique(edges_out[edges_out$type == "Provided", "bait_id"])
  
  # rename nodes file columns
  colnames(nodes_out)[which(colnames(nodes_out) == "localizations")] <- "localization_string_output"
  colnames(nodes_out)[which(colnames(nodes_out) == "GO_terms")] <- "GO_term_string"
  colnames(nodes_out)[which(colnames(nodes_out) %in% c("onset_condition", "terminus_condition"))] <- 
    c("orig_onset", "orig_terminus")
  
  # returns list of: nodes, edges, num_edges, timepoints, localizations, topGOterms,
  # organisms, abund_prov, norm_prot, bait_ids
  return (list("nodes" = nodes_out, 
               "edges" = edges_out,
               "num_edges" = NUM_EDGES,
               "localizations" = LOCS,
               "topGOterms" = topGOterms,
               "organisms" = ORGANISMS,
               "timepoints" = timepoints,
               "abund_prov" = ABUND_PROV,
               "norm_prot" = NORM_PROT,
               "bait_ids" = BAIT_IDS, 
               "named_timepoints" = named_timepoints))
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

calculateComplexesPREV <- function(nodes_in) {
  
  NUM_NODES <- dim(nodes_in)[1]
  
  corum_complexes <- read.table("AppFiles/coreComplexes.txt", header = T, sep = "\t", quote = "", 
                                fill = T, stringsAsFactors = F)
  corum_complexes <- corum_complexes[, c("ComplexName", "subunits.UniProt.IDs.")]
  
  list_complexes <- strsplit(corum_complexes[, "subunits.UniProt.IDs."], ";")
  names(list_complexes) <- corum_complexes$ComplexName
  
  # only keep complexes with at least 3 members (no duplexes)
  list_complexes <- Filter(function(x) length(x) > 2, list_complexes)
  
  # only keep complexes with at least 40% of its members in the provided dataset
  keep_function <- function(x) {
    if (sum(x %in% nodes_in$accession) / length(x) >= 0.4 ) { 
      return (TRUE) 
    } else { return (FALSE) }
  }
  list_complexes <- Filter(keep_function, list_complexes)
  
  chr2string <- function(x) {
    gene_names <- nodes_in[nodes_in$accession %in% x, "gene_name"]
    return (paste0(gene_names, collapse="; "))
  }
  
  percent_detected <- function(x) {
    return ( sum(x %in% nodes_in$accession) / length(x) )
  }
  
  # output complexes table (for download)
  
  if (length(list_complexes) != 0) {
    
    detected_complexes <- cbind.data.frame(names(list_complexes), 
                                           sapply(list_complexes, percent_detected),
                                           sapply(list_complexes, chr2string))
    colnames(detected_complexes) <- c("Complex Name", "Fraction of Complex Detected",
                                      "Detected Complex Members")
    detected_complexes <- detected_complexes[order(-detected_complexes$`Fraction of Complex Detected`), ]
    rownames(detected_complexes) <- 1:nrow(detected_complexes)
    
  } else {
    
    # no complexes found in CORUM
    
    detected_complexes <- NULL
    
  }
  
  # function returns: detected complexes
  
  return (detected_complexes)
  
}

##### QUANTITATIVE FUNCTIONS ####

## CLUSTER PROTEINS BASED ON ABUNDANCES OVER TIME

# only use this function if abundances provided

clusterProteins <- function(nodes_in, edges_in, bait_ids, timepoints) {
  
  theme_set(theme_bw())
  
  edges_out <- edges_in
  
  edges_out$cluster <- NA
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

clusterProteinsPREV <- function(edges_in, timepoints) {
  
  theme_set(theme_bw())
  
  edges_in <- edges_in[edges_in$type == "Provided", ]
  abundances <- edges_in[!duplicated(edges_in[3:4]), 
                         c("bait_gene_name", "prey_gene_name", "cluster", 
                           paste("a", timepoints, sep=""))]
  
  # scale values to max 1
  maxOne <- function(x) {
    return ( x / max(x) )
  }
  
  to_scale <- abundances[, c(4:ncol(abundances))]
  abundances[, c(4:ncol(abundances))] <- t(apply(to_scale, 1, maxOne))
  
  abundances <- abundances[!is.na(abundances$cluster), ]
  
  # return relative abundances file
  return (abundances)
}

# function to plot individual protein profiles split by cluster

clusterPlot <- function(bait_id, nodes_in, relative_abundances, timepoints, named_timepoints) {
  
  custom_colors <- c('#a14949', '#69802e', '#5c86b4', '#896dc1', '#329a90', '#954f72', '#6b696a')
  
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
          scale_x_continuous(breaks = timepoints, labels=named_timepoints) +
          scale_color_manual(values = custom_colors) +
          labs(title=paste("Bait", bait_gene_name), x="Time", y="Scaled Relative Abundance", 
               color="Cluster") +
          facet_wrap(~as.factor(cluster))+ 
          theme(text = element_text(size=20)))
  
}

## NORMALIZE BY PROTEOME ABUNDANCE 

# only use this function if abundances provided and proteome abundances file provided

normalizeByProteomeAbund <- function(proteome_abundance, nodes_in, edges_in, timepoints) {
  
  theme_set(theme_bw())

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
  
  # scale values from 0 to 1
  zero2one <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)) )
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

normalizeByProteomeAbundPREV <- function(edges_in, timepoints) {
  
  theme_set(theme_bw())
  
  edges_in <- edges_in[edges_in$type == "Provided", ]
  
  abundances <- edges_in[!duplicated(edges_in[3:4]), 
                         c("bait_gene_name", "prey_gene_name", "cluster", 
                           paste("a", timepoints, sep=""),
                           paste("norm_a", timepoints, sep=""))]
  
  # scale values from 0 to 1
  zero2one <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)) )
  }
  
  scaled_reg <- t(apply(abundances[, c(paste("a", timepoints, sep=""))], 1, 
                        zero2one))
  scaled_norm <- t(apply(abundances[, c(paste("norm_a", timepoints, sep=""))], 1, 
                         zero2one))
  
  abundances <- cbind.data.frame(abundances[, c("bait_gene_name", "prey_gene_name",
                                                "cluster")],
                                 scaled_reg, scaled_norm, stringsAsFactors = F)
  
  # this should return heatmap dataframe
  return ("normalized_abundances" = abundances)
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
      facet_wrap(vars(type)) +
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

sizeByAbundPREV <- function(nodes_in, edges_in, bait_ids, timepoints) {
  
  nodes_out <- nodes_in
  edges_in <- edges_in[edges_in$type == "Provided", ]
  
  NUM_TPS <- length(timepoints)
  
  # convert abundance to size and size by quartile
  
  overall_avg <- vector(mode='numeric', length=NUM_TPS)
  for (i in 1:NUM_TPS) {
    overall_avg[i] <- mean(edges_in[, c(paste("a", timepoints[i], sep=""))])
  }
  
  # set 2.5 as the average node size
  mult_factor <- 2.5 / mean(overall_avg)
  
  returnSize <- function(x) {
    y <- max(1, mult_factor*x) # prevent nodes from being too small
    y <- min(4, y) # prevent nodes from being too big
    return (y)
  }
  
  for (i in 1:NUM_TPS) {
    size <- lapply(edges_in[, c(paste("a", timepoints[i], sep=""))], returnSize)
    size <- cbind(edges_in$prey_id, unlist(size))
    size <- rbind(size, c(bait_ids, 2.5))
    colnames(size) <- c("id", paste("size_", timepoints[i], sep=""))
    nodes_out <- merge(nodes_out, size, by="id")
  }
  
  nodes_out <- nodes_out[!duplicated(nodes_out[ , "accession"]), ]
  
  return (nodes_out)
}

##### GET NODE ANNOTATIONS FROM UNIPROT ####

# get UniProt data for proteins in each species and background list for enrichment
# and assign localizations to nodes

getUniprotData <- function(background_list, nodes_in, taxids, loc_prov) {
   
  nodes_out <- nodes_in
  
  if (ncol(background_list) > 1) {
    
    # if the provided list already has Uniprot information, skip the extraction step for 
    # all genes that already have the information provided
    uniprot_list <- nodes_out[!(nodes_out$accession %in% background_list$accession), "accession"]
    
  } else {
    # otherwise add any nodes whose accessions do not appear in the background list
    uniprot_list <- c(nodes_out$accession, background_list[,1])
    uniprot_list <- unlist(unique(uniprot_list))
  }
  
  # if there is more data to be collected, collect the data
  if (length(uniprot_list) > 0) {
    
    updateUniprot <- function(uniprot_chunk) {
      
      tryCatch({
        run_uniprot_query(uniprot_chunk)
      }, error = function (err) {
        return (list())
      })
      
    }
    
    # break into chunks of 2000 if larger than 2000 
    CHUNK_SIZE <- 2000
    
    NUM_CHUNKS <- (length(uniprot_list) %/% CHUNK_SIZE)+1
    
    uniprot_data <- list()
    
    for (i in 1:NUM_CHUNKS) {
      
      n_start <- CHUNK_SIZE*(i-1) + 1
      n_end <- min(CHUNK_SIZE*i, length(uniprot_list))
      
      uniprot_chunk <- uniprot_list[n_start:n_end]
      
      uniprot_data <- append(uniprot_data, updateUniprot(uniprot_chunk))
      
    }
    
    if (length(uniprot_data) > 0) {
      # if new uniprot data has been collected
      # create dataframe of newly acquired uniprot data
      up_accession <- names(uniprot_data)
      up_genename <- sapply(names(uniprot_data), function(nm) { 
        paste(uniprot_data[[nm]][["genename"]], collapse = "; " ) } )
      up_GOIDs <- sapply(names(uniprot_data), function(nm) { 
        paste(uniprot_data[[nm]][["GO IDs"]], collapse = "; " ) } )
      up_localizations <- sapply(names(uniprot_data), function(nm) { 
        paste(uniprot_data[[nm]][["localizations"]], collapse = "; " ) } )
      up_organism <- sapply(names(uniprot_data), function(nm) { 
        paste(uniprot_data[[nm]][["organism"]], collapse = "; " ) } )
      
      uniprot_df <- cbind(up_accession, up_genename, up_GOIDs, up_localizations, 
                          up_organism)
      colnames(uniprot_df) <- c("accession", "genename", "GOIDs", "localizations", "organism")
      rownames(uniprot_df) <- 1:nrow(uniprot_df)
      
      uniprot_df <- as.data.frame(uniprot_df, stringsAsFactors = F)
      
      # add a column that specifies whether the protein is a node in the dataset,
      # part of the background, or both
      
      uniprot_df$type <- "Background Only"
      uniprot_df[uniprot_df$accession %in% nodes_out$accession, "type"] <- "In NodeList Only"
      
      if (ncol(background_list) > 1) {
        bg_terms <- background_list[background_list$type != "In NodeList Only", "accession"]
        all_accessions <- c(background_list$accession, uniprot_df$accession)
        
      } else {
        bg_terms <- background_list[,1]
        all_accessions <- c(bg_terms, uniprot_df$accession)
        
      }
      
      # determine which proteins are in both the background and nodes files
      both_accessions <- all_accessions[(all_accessions %in% bg_terms) & (all_accessions %in% nodes_out$accession)]
      
      # combine newly collected uniprot data with previously collected
      if (ncol(background_list) > 1) {
        
        uniprot_df <- rbind.data.frame(uniprot_df, background_list)
        uniprot_df <- uniprot_df[!duplicated(uniprot_df[, "accession"]), ]
        
      # assign proteins in background and nodes as type "both"
      uniprot_df[uniprot_df$accession %in% both_accessions, "type"] <- "Both"
        
      }
      
    } else {
      
      # no new uniprot data collected
      
      uniprot_df <- background_list
      
    }
    
    
  } else {
    
    # otherwise simply retain the background list
    
    uniprot_df <- background_list
    
  }
  
  # merge organism name into nodes file
  nodes_out <- merge(nodes_out, uniprot_df[, c("accession", "organism")], by="accession", all.x=T)
  
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
  
  
  # create nodes Uniprot annotations
  node_annotations <- uniprot_df[uniprot_df$type != "Background Only", ]
  
  # create functional GO annotations list to later pipe into topGO
  accessionToGO <- strsplit(uniprot_df$GOIDs, "; ")
  names(accessionToGO) <- uniprot_df$accession
  
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
    locs <- c("Cell Membrane", "Cytoplasm", "Endoplasmic Reticulum", "Golgi", "Lysosome",
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
  localization_string <- vector(mode='character', length=NUM_NODES) # for node colors
  localization_string_output <- vector(mode='character', length=NUM_NODES) # for node tooltip
  
  for (loc in locs) {
    index <- which(nodes_out[, paste('loc_', loc, sep="")] == 1)
    if (length(index) > 0) {
      localization_string_output[index] <- 
        paste(localization_string_output[index], loc, sep="; ")
      for (j in index) {
        if (localization_string[j] == "" ) { localization_string[j] <- loc; }
        # multiple possiblities
        else { localization_string[j] <- "Multiple localizations" }
      }
    }
  }
  
  nodes_out$localization_string <- localization_string
  nodes_out$localization_string_output <- sub("^; ", "", localization_string_output)
  
  # if node lacks localization information, say "Unspecified"
  nodes_out[nodes_out$localization_string == "", "localization_string"] <- "Unspecified"
  nodes_out[nodes_out$localization_string_output == "", "localization_string_output"] <- "Unspecified"
  
  # function returns: nodes, organisms, localizations, uniprot, GOList
  return (list("nodes" = nodes_out,
               "organisms" = ORGANISMS,
               "localizations" = locs,
               "uniprot" = uniprot_df,
               "GOList" = accessionToGO))
}

##### ASSIGN FUNCTIONAL ANNOTATIONS AFTER GO ENRICHMENT #####

performGOEnrichment <- function(GOList, nodes_in) {
  
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

performGOEnrichmentPREV <- function(background_list, nodes_in) {
  
  # create functional GO annotations list to later pipe into topGO
  accessionToGO <- strsplit(background_list$GOIDs, "; ")
  names(accessionToGO) <- background_list$accession
  
  gene_universe <- background_list$accession
  interesting_genes <- nodes_in$accession
  
  # assign whether a gene is "interesting" - in the provided data
  gene_list <- factor(as.integer(gene_universe %in% interesting_genes))
  names(gene_list) <- gene_universe
  
  # create topGO object
  # Biological Process (BP) ontology
  GOdata <- new("topGOdata", ontology = "BP", allGenes = gene_list, 
                annot = annFUN.gene2GO, gene2GO = accessionToGO)
  
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
    
    # list the genes for each of the significant GO terms
    annotated_genes <- lapply(allRes$GO.ID,
                              function(x) as.character(
                                unlist(genesInTerm(object = GOdata, whichGO = x))))
    names(annotated_genes) <- allRes$Term
    # get GO terms for the original genes in the provided list
    significant_genes <- lapply(annotated_genes, function(x) intersect(x, interesting_genes))
    
    # output GO Table (for download)
    GOTable <- allRes
    NUM_GO_TERMS <- dim(GOTable)[1]
    GOTable$Proteins <- vector(mode = "character", length = NUM_GO_TERMS)
    
    for (i in 1:NUM_GO_TERMS) {
      GOTable$Proteins[i] <- paste0(nodes_in$gene_name[which(nodes_in$accession %in%
                                                               unlist(significant_genes[GOTable$Term[i]]))],
                                    collapse = "; ")
    }
    
  } else {
    
    # no significant GO terms found
    GOTable <- NULL
  }
  
  # function returns: GO terms table
  return (GOTable)
  
}

##### CALCULATE ACTIVITY SPELLS ####

calculateActivitySpells <- function(nodes_in, edges_in, timepoints, named_timepoints, weight_threshold,
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
  temp_timepoints <- 1:NUM_TPS
  weights <- edges_out[,grep("^w", colnames(edges_out))]
  for (n in (1:NUM_EDGES)) {
    j <- 1
    freq[n] <- 0
    while (j <= NUM_TPS) {
      count <- 0
      w <- weights[n, j]
      if (w >= weight_threshold) {
        f <- f+1
        freq[n] <- freq[n] + 1
        onset[f] <- j
        while (w >= weight_threshold) {
          count <- count+1
          j <- j+1
          if (j %in% temp_timepoints) {
            w <- weights[n, j] 
          } else { w <- -1 } # set a negative weight to break out of while loop
        }
        terminus[f] <- j
      }
      j <- j+1
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
    add_baits$onset <- 1
    add_baits$terminus <- NUM_TPS+1
    
    keep_nodes <- rbind(keep_nodes, add_baits)
  }
  
  nodes_out <- keep_nodes 
  nodes_out$onset[which(nodes_out$id %in% bait_ids)] <- 1
  nodes_out$terminus[which(nodes_out$id %in% bait_ids)] <- NUM_TPS+1
  
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
  
  # also save the original timepoints / conditions for onset and terminus information
  edges_out$orig_onset <- named_timepoints[edges_out$onset]
  edges_out$orig_terminus <- named_timepoints[edges_out$terminus]
  nodes_out$orig_onset <- named_timepoints[nodes_out$onset]
  nodes_out$orig_terminus <- named_timepoints[nodes_out$terminus]
  
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
  string_links$orig_onset <- NA
  string_links$orig_terminus <- NA
  
  # assign type to be STRING inferred or both
  both <- merge(edges_out, string_links, by=c("bait_id", "prey_id"))
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
buildNetwork <- function(nodes_in, edges_in, num_edges, abund_prov, norm_prot, timepoints,
                         bait_ids) {
  
  edges_out <- edges_in
  nodes_out <- nodes_in
  
  # rearrange link dataframes to network
  if (abund_prov && norm_prot) {
    edges_out <- edges_out[, c("bait_id", "prey_id", "bait_accession", "prey_accession", 
                               "bait_gene_name", "prey_gene_name", "cluster", "onset", 
                               "terminus", "orig_onset", "orig_terminus", "type", 
                               paste("w", timepoints, sep=""),
                               paste("a", timepoints, sep=""),
                               paste("norm_a", timepoints, sep=""))]
  } else if (abund_prov && !norm_prot) {
    edges_out <- edges_out[, c("bait_id", "prey_id", "bait_accession", "prey_accession", 
                               "bait_gene_name", "prey_gene_name", "cluster", "onset", 
                               "terminus", "orig_onset", "orig_terminus", "type", 
                               paste("w", timepoints, sep=""),
                               paste("a", timepoints, sep=""))]
    
  } else {
    edges_out <- edges_out[, c("bait_id", "prey_id", "bait_accession", "prey_accession", 
                               "bait_gene_name", "prey_gene_name", "onset", "terminus", 
                               "orig_onset", "orig_terminus", "type", 
                               paste("w", timepoints, sep=""))]
  }
  
  colnames(edges_out)[c(1,2)] <- c("from", "to")
  
  temp_locs <- unlist(nodes_out$localization_string)
  
  # give localization, duration and edge type a factor numeric value (to later assign color)
  nodes_out$locanum <- as.numeric(factor(temp_locs))
  # give duration factor based on what quartile of expression it is in
  nodes_out$durfactor <- cut(nodes_out$duration, unique(quantile(nodes_out$duration, probs=0:4/4)),
                             include.lowest=T, labels=F)
  edges_out$edge_color <- c("black", "navyblue", "gray75")[as.numeric(factor(edges_out$type))]
  
  # create network object
  
  net3 <- network(edges_out, vertex.attr=nodes_out, matrix.type="edgelist", 
                  loops=F, multiple=F, directed=F, ignore.eval = F)
  
  # assign colors to localization, duration factors
  custom_colors <- c('#c69191', '#a5b281', '#9db6d2', '#b8a7d9', '#84c2bc', '#bf95aa', '#a6a5a5')
  newPalette <- colorRampPalette(custom_colors)
  loc_color <- newPalette(n = nlevels(factor(temp_locs))+1)
  
  newPalette_1 = colorRampPalette(c('#A14949', '#ecdada'))
  dur_color <- newPalette_1(n = nlevels(factor(temp_locs))+1)
  
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
  weights <- edges_out[,grep("^w", colnames(edges_out))]
  sizes <- nodes_out[,grep("^size", colnames(nodes_out))]
  
  for (i in 1:NUM_TPS) {
    net3.dyn <- activate.edge.attribute(net3.dyn, 'edge_weight', 
                                        (weights[, i])*2, 
                                        onset=i, terminus=i+1, e=1:NUM_EDGES)
    # set dynamic vertex sizes (if single bait)
    if (length(bait_ids) == 1 & abund_prov) {
      net3.dyn <- activate.vertex.attribute(net3.dyn, 'node_size',
                                            sizes[, i],
                                            onset=i, terminus=i+1, v=1:NUM_NODES)
    }
  }
  
  # function returns: nodes, edges, dynamic net3 object
  return (list("nodes" = nodes_out,
               "edges" = edges_out,
               "net3.dyn" = net3.dyn))
}

# output protein information
saveNodeAttributes <- function(nodes_in, bait_ids, timepoints, named_timepoints, filename) {
  
  output_nodes <- nodes_in[, c("gene_name", "accession", "organism", 
                               "localization_string_output", "complexes", 
                               "GO_term_string", "orig_onset", "orig_terminus", 
                               "onset", "terminus")]
  output_nodes$GO_term_string <- sub("^; ", "", output_nodes$GO_term_string)
  output_nodes$localization_string_output <- sub("^; ", "",
                                                 output_nodes$localization_string_output)
  output_nodes$organism <- unlist(output_nodes$organism)
  colnames(output_nodes) <- c("gene_name", "accession", "species", "localizations",
                              "complexes", "GO_terms", "onset_condition", "terminus_condition",
                              "onset", "terminus")
  
  NUM_BAITS <- length(bait_ids)
  
  
  if (NUM_BAITS > 1) {
    output_nodes <- cbind(output_nodes,
                          nodes_in[, paste("shared_baits_", named_timepoints, sep="")])
  }
  
  # replace NA values with empty string
  output_nodes[is.na(output_nodes)] <- ""
  
  write.table(output_nodes, file = filename, sep = "\t", 
              quote = F, row.names = F)
}

# output edge information
saveEdgeAttributes <- function(edges_in, timepoints, named_timepoints, abund_prov, norm_prot, 
                               filename) {
  
  # get edge information to output
  output_edges <- edges_in[, c("bait_gene_name", "prey_gene_name",
                               "bait_accession", "prey_accession",
                               "orig_onset", "orig_terminus",
                               "onset", "terminus", "type", 
                               paste("w", timepoints, sep=""))]
  output_edges_colnames <- c("#node1", "node2", "node1_accession", "node2_accession", 
                             "onset_condition", "terminus_condition",
                             "onset", "terminus", "type", 
                             paste("confidence_score_", named_timepoints, sep=""))
  
  if (abund_prov) {
    output_edges <- cbind(output_edges, edges_in[, c("cluster", 
                                                     paste("a", timepoints, sep=""))])
    output_edges_colnames <- c(output_edges_colnames, "cluster",
                               paste("abundance_", named_timepoints, sep=""))
    if (norm_prot) {
      output_edges <- cbind(output_edges, edges_in[, paste("norm_a", timepoints, sep="")])
      output_edges_colnames <- c(output_edges_colnames, 
                                 paste("normalized_to_proteome_abundance_", 
                                       named_timepoints, 
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
                         edges_in, timepoints, bait_ids) {
  
  # calculate shared nodes
  plot_net <- net3.dyn
  # deactivate nodes based on number of baits they share
  if (length(bait_ids) > 1) {
    
    for (i in 1:length(timepoints)) {
      tp <- timepoints[i]
      keep_nodes <- nodes_in[nodes_in[, paste("shared_", tp, sep="")] >= min_shared, "id"]
      keep_nodes <- c(keep_nodes, bait_ids)
      keep_nodes <- unique(keep_nodes)
      deactivate_nodes <- nodes_in[!(nodes_in$id %in% keep_nodes), "id"]
      deactivate.vertices(plot_net, onset=i, terminus=i+1, v=deactivate_nodes, 
                          deactivate.edges=T)
    }
    
  }
  
  # deactivate edges (based on STRING confidence threshold)
  deactivate_link_ids <- which((edges_in$type != "Provided") &
                                 (edges_in[, paste("w", timepoints[1], sep="")] <= st))
  
  if (!is.null(deactivate_link_ids)) {
    deactivate.edges(plot_net, e=deactivate_link_ids)
  }
  
  # deactivate nodes (based on localizations, species, functions, etc.)
  
  # localizations
  keep_loc_vertex_ids <- {}
  if ("All" %in% locns) {
    keep_loc_vertex_ids <- nodes_in$id
  } else {
    for (loc in locns) {
      keep_loc_vertex_ids <- append(keep_loc_vertex_ids, 
                                    nodes_in[which(nodes_in[, paste("loc_", loc, sep="")] == 1), "id"])
    }
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
    
    deactivate.vertices(plot_net, onset=1, terminus=NUM_TPS+1,
                        v=deactivate_vertex_ids, deactivate.edges=T)
  }
  
  # compute network animation
  compute.animation(plot_net, animation.mode = "kamadakawai", default.dist=2, 
                    slice.par=list(start=1, end=NUM_TPS, 
                                   interval=1, 
                                   aggregate.dur=1, rule='earliest'))
  return (plot_net)
}

##### PLOT DYNAMIC NETWORK ####

# only run this function if network exists

plotCompNet <- function(net3_network, lab, disp_tps, bait_ids, abund_prov, named_timepoints) {
  
  if (disp_tps){
    xlab <- function(onset,terminus){ifelse(onset==length(named_timepoints),paste(named_timepoints[onset],sep=''),
                                            paste(named_timepoints[onset]," to ", named_timepoints[terminus],sep=''))}
  } else {
    xlab <- ''
  }
  
  if(length(bait_ids) == 1 & abund_prov) {
    # plot network animation
    render.d3movie(net3_network,
                   render.par=list(tween.frames = 10, show.time = F),
                   plot.par=list(bg='white', xlab = xlab),
                   slice.par=list(rule="latest"),
                   d3.options = list(animationDuration=4500,enterExitAnimationFactor=0.3),
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
                                                          slice%v%'localization_string_output',
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
                   plot.par=list(bg='white', xlab = xlab),
                   d3.options = list(animationDuration=4500,enterExitAnimationFactor=0.3),
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
                                                          slice%v%'localization_string_output',
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
          scale_x_continuous(breaks = 1:length(timepoints)) +
          labs(x = "Time", y = "Number of Edges Formed", title = "Time of Edge Formation"))
  
  # edge duration
  provided_edges <- which(edges_in$type == "Provided")
  data <- data.frame(x=edgeDuration(net3.dyn, e=provided_edges))
  print(ggplot(data) +
          geom_histogram(mapping = aes(x = data$x)) +
          scale_x_continuous(breaks = 1:length(timepoints)) +
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
  
  if ("All" %in% locns) {
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
                                 edges_in$type != "Provided"), ]
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
                      nodes_in, plotabundances, timepoints, named_timepoints) {
  
  custom_colors <- c('#a14949', '#69802e', '#5c86b4', '#896dc1', '#329a90', '#954f72', '#6b696a')

  # decimal places to show in label
  scalefunction <- function(x) sprintf("%.2f", x)
  
  # determine whether to search for exact gene name matches or not
  exact_match <- FALSE
  
  if (grepl('"', genesyms)){
    exact_match <- TRUE
  }
  
  genes <- gsub('"', '', genesyms) # remove quotes if present
  genes <- unlist(strsplit(genes, ";")) # split on semicolon operator
  genes <- mapply(trimws, genes) # trim leading or trailing whitespace from each gene
  names(genes) <- NULL
  
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
    if (exact_match){
      current_genes <- rbind(current_genes, 
                             plotabundances[which(plotabundances$prey_gene_name == gene), ])
    } else{
      current_genes <- rbind(current_genes, 
                             plotabundances[grepl(gene, plotabundances$prey_gene_name, ignore.case = T), ])
    }
  }
  
  if (nrow(current_genes) != 0){
  
    if (length(unique(current_genes$prey_gene_name))<= length(custom_colors)){
      
      print(ggplot(data = current_genes) +
            geom_line(mapping = aes(x=time, y=abundance, color=prey_gene_name, linetype=bait_gene_name), size=1.5) +
            scale_y_continuous(trans='log', labels=scalefunction) +
            scale_x_continuous(breaks = timepoints, labels=named_timepoints) + 
            labs(x="Time", y="Abundance", color="Search Genes", 
                 linetype="Bait") +
            scale_color_manual(values = custom_colors) +
            theme(text = element_text(size=20))
    )
      
    } else {
      
      print(ggplot(data = current_genes) +
            geom_line(mapping = aes(x=time, y=abundance, color=prey_gene_name, linetype=bait_gene_name), size=1.5) +
            scale_y_continuous(trans='log', labels=scalefunction) +
            scale_x_continuous(breaks = timepoints, labels=named_timepoints) + 
            labs(x="Time", y="Abundance", color="Search Genes", 
                 linetype="Bait") +
            theme(text = element_text(size=20))
      )
      
    }
  }
  
}

printClusterNumber <- function(genesyms, nodes_in, edges_in, bait_ids) {
  
  # determine whether to search for exact gene name matches or not
  exact_match <- FALSE
  
  if (grepl('"', genesyms)){
    exact_match <- TRUE
  }
  
  genes <- gsub('"', '', genesyms) # remove quotes if present
  genes <- unlist(strsplit(genes, ";")) # split on semicolon operator
  genes <- mapply(trimws, genes) # trim leading or trailing whitespace from each gene
  names(genes) <- NULL
  
  for (gene in genes) {
    for (bait in bait_ids) {
      bait_gene_name <- nodes_in$gene_name[which(nodes_in$id == bait)]
      prey_gene_name <- edges_in[which(edges_in$from == bait & 
                                         grepl(gene, edges_in$prey_gene_name, ignore.case = T)), 
                                 'prey_gene_name']
      clusterNumber <- edges_in[which(edges_in$from == bait & 
                                        grepl(gene, edges_in$prey_gene_name, ignore.case = T)),
                                "cluster"]
      if (length(clusterNumber) != 0) {
        out <- cat(paste(prey_gene_name, " is in cluster ", clusterNumber, " for bait ", 
                         bait_gene_name, "\n", collapse=""))
      }
    }
  }
}

##### BUILD UI ####

## Code to display generic error message when there is an error ----

options(shiny.sanitize.errors = TRUE)

ui = fluidPage(
  
  list(tags$head(HTML('<link rel="icon", href="intervista_icon_1.png", type="image/png" />'))),
  
  titlePanel(title=div(img(src="logo_intervista_new.png", width = 800)), windowTitle = 'Inter-ViSTA'),
  
  # Enable shinyjs
  shinyjs::useShinyjs(),
  
  # Aesthetics ----
  tags$head(
    tags$style(HTML("
                    hr {border-top: 4px solid #000000;}
                    .shiny-notification {
                    height: 150px;
                    width: 600px;
                    position:fixed;
                    top: calc(50%);;
                    left: calc(33.33333%);;
                    }
                    "))
  ),
  
  # Create tabs ----
  tabsetPanel(
    id = "navbar",
    
    # source each tab UI from separate file
    source(file.path("ui", "instructions.R"), local = TRUE)$value,
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
  source(file.path("server", "instructions.R"), local = TRUE)$value
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
