#' Determine taxa labels of tips on each side of a node\cr
#'
#' Takes in a tree, a table of taxa labels per tip, and either a node label or index.\cr
#' 
#' **The format of the taxa label table is very important to note: it must have the tips as the rownames and taxonomic levels (ranging from highest to lowest) as the column names**.
#' 
#' For each side of the specified node separately, this function returns the lowest possible taxon label shared by at
#' least the specified proportion of tips (set by the "threshold" variable). Will return "Unclear" if there is no applicable taxon.\cr
#' 
#' To clarify, the taxon that meets the threshold at the lowest possible taxonomic level will be used as the representative label.
#' For example, if all the tips on one side of the node are members of the Pseudomonas genus, but only 60% are members of the Pseudomonas aeruginosa species specifically,
#' then Pseudomonas will be used as the representative label (based on a threshold of 0.75 or higher and assuming that species are the last column in the table).
#'
#' @param in_tree Phylo object
#' 
#' @param taxon_labels Dataframe of taxa labels for tips in tree.
#' All tips underlying the specified node must be present, although typically all tips in the tree would be present.
#' Row names must be the tip labels. The column names correspond to each taxonomic level, such as Kingdom, Phylum, etc.
#' The actual column names do not matter: it is just important that the order of the taxonomic levels goes from the highest taxonomic level present (e.g., Kingdom), to the lowest taxonomic level present (e.g., Species).
#' 
#' @param node_label Optional label of node for which the representative taxon label will be determined. Either this option or the node_index option must be specified, but not both.
#' 
#' @param node_index As above for the node_label option, but to specify a node by index rather than by label.
#'
#' @param threshold Float > 0.5 and <= 1.0 specifying the proportion of tips that must share a taxon label for it to be considered representative.
#'
#' @param combine_labels Boolean flag for whether taxon labels should be combined, so that all higher taxonomic labels are included.
#' Specifically, when TRUE, all higher labels are concatenated and delimited by "; ".
#' E.g., rather than just the genus "Odoribacter" the label would be "Bacteria; Bacteroidetes; Bacteroidia; Bacteroidales; Porphyromonadaceae; Odoribacter",
#' given that those were the labels of the higher taxonomic levels of that genus.
#'
#' @return Character vector of size two with the representative taxon for tips on each side of the specified node.
#' 
#' @export
node_taxa <- function(in_tree, taxon_labels, node_label=NULL, node_index=NULL, threshold=0.75, combine_labels = TRUE) {

  if (threshold <= 0.5 | threshold > 1) {
    stop("Stopping - the set threshold needs to be > 0.5 and <= 1.")
  }
  
  if (! is.null(node_label) & is.null(node_index)) {
    underlying_tips <- lhs_rhs_tips(in_tree, node_label, get_node_index=TRUE)
  } else if (is.null(node_label) & ! is.null(node_label)) {
    underlying_tips <- lhs_rhs_tips(in_tree, node_index, get_node_index=FALSE)
  } else if (! is.null(node_index) & ! is.null(node_label)) {
    stop("Stopping - only one of node_index and node_label can be specified.")
  } else {
    stop("Stopping - one of node_index or node_label must be specified.")
  }
  
  if (length(which(! c(underlying_tips$lhs, underlying_tips$rhs) %in% rownames(taxon_labels))) > 0) {
    message("The following tips are missing as row names from the taxa label table.")
    message(c(underlying_tips$lhs, underlying_tips$rhs)[which(! c(underlying_tips$lhs, underlying_tips$rhs) %in% rownames(taxon_labels))])
    stop("Stopping - please make sure the row names of the input taxa label table are the tip labels of the tree.")
  }

  # Determine representative taxon.
  taxon_lhs <- feature_consensus_taxon(taxa=taxon_labels, features=underlying_tips$lhs, threshold=threshold, combine_labels=combine_labels)
  taxon_rhs <- feature_consensus_taxon(taxa=taxon_labels, features=underlying_tips$rhs, threshold=threshold, combine_labels=combine_labels)
  
  return(c(taxon_lhs, taxon_rhs))

}


#' Compute isometric log ratio based on abundance of feature sets
#'
#' Computes isometric log ratio between two sets of feature abundances, for each sample separately. Requires an abundance table,
#' with two sets of features for which the ratio will be computed. 
#'
#' @param abun_table Abundance table, e.g., read counts or relative abundance.
#' Should be dataframe with column names correspond to sample names and row names corresponding to the feature ids.
#' No 0's are permitted unless the "pseudocount" option is set.
#'
#' @param set1_features Features (rows of abundance table) that make up one side of the ratio to be computed (numerator).
#' 
#' @param set2_features Same as "set1_features", but for the other side of the ratio (denominator).
#' 
#' @param pseudocount Constant to add to all abundance values, to ensure that there are only non-zero values. For read count data this would typically be 1.
#' 
#' @return Numeric vector of the computed isometric log ratio for each sample (where samples are taken to be each column in the input table).
#' 
#' @export
abun_isometric_log_ratios <- function(abun_table, set1_features, set2_features, pseudocount=NULL) {
  
  num_set1 <- length(set1_features)
  num_set2 <- length(set2_features)
  
  if (num_set1 == 0 | num_set2 == 0) {
    stop("At least one feature must be specified in each set.")
  }
  
  # Check that features are in table.
  if (length(which(! c(set1_features, set2_features) %in% rownames(abun_table))) > 0) {
    stop("Stopping - at least one feature in the specified sets is not present as a row name in the abundance table.") 
  }
  
  # Check that no features intersect between the sets.
  if (length(which(set1_features %in% set2_features)) > 0) {
    stop("Stopping - at least one feature overlaps between the input sets.") 
  }
  
  abun_table <- abun_table[c(set1_features, set2_features), ]
  
  if (! is.null(pseudocount)) {
    abun_table <- abun_table + pseudocount
  }
  
  if (length(which(abun_table == 0)) > 0) {
    stop("At least one 0 is present in the abundance table, which means that at least some isometric log ratios cannot be computed. Fix this by setting (or changing) the \"pseudocount\" option.")
  }
  
  half1 <- sqrt((num_set1 * num_set2) / (num_set1 + num_set2))
  
  sample_names <- colnames(abun_table)
  
  half2 <- c()
  
  for(sample in sample_names) {
    sample_set1_values <- abun_table[set1_features, sample]
    sample_set2_values <- abun_table[set2_features, sample]
    
    gmean_set1 <- exp(mean(log(sample_set1_values)))
    gmean_set2 <- exp(mean(log(sample_set2_values)))
    
    half2 <- c(half2, log(gmean_set1 / gmean_set2))
  }
  
  ilr_out <- half1 * half2
  
  names(ilr_out) <- sample_names
  
  return(ilr_out)

}


#' Compute balances at tree nodes. 
#' 
#' Computes balances (i.e., isometric log ratios, for each sample separately) of feature abundances at each non-negligible node in the tree.
#' 
#' @param tree Phylo object with tip labels matching row names of input abundance table. Note that node labels are required.
#'
#' @param abun_table Abundance table, e.g., read counts or relative abundance.
#' Should be dataframe with column names correspond to sample names and row names corresponding to the tips of the tree.
#' No 0's are permitted unless the "pseudocount" option is set.
#' 
#' @param min_num_tips Minimum number of tips that must be found on each side of a node for it to be included (i.e., to be considered non-negligible).
#'
#' @param ncores Number of cores to use for steps of function that can be run in parallel.
#'
#' @param pseudocount Optional constant to add to all abundance values, to ensure that there are only non-zero values. For read count data this would typically be 1.
#'
#' @param derep_nodes Boolean setting to specify whether nodes should be dereplicated based on the Jaccard similarity of the underlying tips.
#' When TRUE, nodes with pairwise Jaccard similarity >= jaccard_cutoff will be collapsed into the same cluster. A node will be added to a cluster if it is adequately similar to any nodes in a cluster.
#' One representative per cluster will be retained, which will correspond to the node with the fewest underlying tips. Note that this step is performed after the step involving the min_num_tips screening.
#'
#' @param jaccard_cutoff Numeric vector of length 1. Must be between 0 and 1 (inclusive). Corresponds to the Jaccard cut-off used for clustering nodes based on similar sets of underlying tips.
#'
#' @param subset_to_test Optional vector of node labels (*not indices*) that correspond to the subset of nodes that should be considered.
#' Note that balances will still only be computed at each of these nodes if they have a sufficient number of underlying tips (as specified by the "min_num_tips" argument).
#' If this argument is not specified then all nodes will be considered.
#'
#' @return List containing three objects:
#' 
#' "tips_underlying_nodes": the tips on the left-hand side (lhs; the numerator) and right-hand side (rhs; the denominator) of each node.
#' Note that which side of the node is denoted as the left-hand or right-hand side is arbitrary.
#' 
#' "balances": list with each non-negligible node as a separate element. The sample balances for each node are provided as a numeric vector within each of these elements.
#'  
#' "negligible_nodes": character vector of node labels considered negligible. This is defined as those with fewer tips on either side of the node than specified by the "min_num_tips" argument.
#' 
#' 
#' When derep_nodes = TRUE, additional elements will also be returned:
#' 
#' "ignored_redundant_nodes": character vector of (non-negligible) node labels ignored due to being in sharing high Jaccard similarity with at least one other node.
#' 
#' "node_pairwise_jaccard": dataframe of pairwise Jaccard similarity for all non-negligible nodes.
#' 
#' "node_clusters": list with the node labels clustered into each unique cluster of nodes based on Jaccard similarities.
#' Each list element is a separate cluster for which only one node was selected as a representative (whichever one had the fewest underlying tips).
#' 
#' @export
compute_node_balances <- function(tree,
                                  abun_table,
                                  min_num_tips = 10,
                                  ncores = 1,
                                  pseudocount = NULL,
                                  derep_nodes = FALSE,
                                  jaccard_cutoff = 0.75,
                                  subset_to_test = NULL) {
  
  if (is.null(tree$node.label)) {
    stop("Stopping - input tree does not have any node labels.") 
  }
  
  if (! all(tree$tip.label %in% rownames(abun_table))) {
    stop("Stopping - not all tips are found as row names in the abundance table.")
  }
  
  # Test all nodes unless subset specified.  
  if(! is.null(subset_to_test)) {
    if (length(which(! subset_to_test %in% tree$node.label)) > 0) {
      stop("Stopping - some labels in subset_to_test do not match node labels in the tree.")  
    }
    nodes_to_test <- subset_to_test
  } else {
    nodes_to_test <- tree$node.label
  }
  
  # Get tips on either side of each node.
  node_features <- parallel::mclapply(nodes_to_test,
                                      lhs_rhs_tips,
                                      tree = tree,
                                      get_node_index = TRUE,
                                      mc.cores = ncores)
  
  names(node_features) <- nodes_to_test
  
  # Identify negligible nodes based on having insufficient underlying tips.
  negligible_nodes <- sapply(names(node_features),
                             function(x) {
                               lhs_feat_num <- length(node_features[[x]]$lhs)
                               rhs_feat_num <- length(node_features[[x]]$rhs)
                               if((lhs_feat_num < min_num_tips) || (rhs_feat_num < min_num_tips)) {
                                 return(TRUE)
                               } else {
                                 return(FALSE) 
                               }
                             })
  
  negligible_nodes_i <- which(negligible_nodes)
  if(length(negligible_nodes_i) > 0) {
    negligible_nodes <- names(negligible_nodes[negligible_nodes_i])
  } else {
    negligible_nodes <- as.character()
  }
  
  nodes_to_consider <- names(node_features)[which(! names(node_features) %in% negligible_nodes)]
  
  if (length(nodes_to_consider) == 0) {
    stop("Stopping - no non-negligible nodes remain after filtering based on mininum number of tips on left and right-hand side of each node.")
  }
  
  if (derep_nodes) {
    jaccard_derep_out <- jaccard_derep_nodes(tree = tree,
                                             node_features = node_features[nodes_to_consider],
                                             cutoff = jaccard_cutoff,
                                             ncores = ncores)

    nodes_to_consider <- nodes_to_consider[which(! nodes_to_consider %in% jaccard_derep_out$ignored_redundant_nodes)]
      
  }
  
  
  
  # Calculate balances at each node.
  balance_calc <- parallel::mclapply(nodes_to_consider,
                                       function(x) {
                                         return(abun_isometric_log_ratios(abun_table = abun_table,
                                                                          set1_features = node_features[[x]]$lhs,
                                                                          set2_features = node_features[[x]]$rhs,
                                                                          pseudocount = pseudocount))
                                       },
                                       mc.cores=ncores)
    
  names(balance_calc) <- nodes_to_consider
  
  if (derep_nodes) {
    
    return(list(tips_underlying_nodes = node_features,
                balances = balance_calc,
                negligible_nodes = negligible_nodes,
                ignored_redundant_nodes = jaccard_derep_out$ignored_redundant_nodes,
                node_pairwise_jaccard = jaccard_derep_out$node_jaccard,
                node_clusters = jaccard_derep_out$node_clusters))
    
  } else {
  
    return(list(tips_underlying_nodes = node_features,
                balances = balance_calc,
                negligible_nodes = negligible_nodes))
  }
  
}


# Get representative taxon label for set of features.
feature_consensus_taxon <- function(taxa, features, threshold, combine_labels = FALSE) {

  taxa_subset <- taxa[features, , drop = FALSE]
  
  for (i in ncol(taxa_subset):1) {
    
    if (combine_labels & i > 1) {
      labels <- apply(taxa_subset[ , 1:i], 1, paste, collapse = "; ")
    } else {
      labels <- taxa_subset[, i, drop=TRUE]
    }
    
    labels_table <- table(labels)
    
    if (length(labels_table) == 0) { next }
    
    if ((max(labels_table) / length(features)) >= threshold) {

      return(paste(names(labels_table)[which(labels_table == max(labels_table))], " ", "(", colnames(taxa_subset)[i], ")", sep=""))
    
    }
  }
  
  return("Unclear")
}


# Get tips on left-hand and right-hand side of node.
lhs_rhs_tips <- function(tree, node, get_node_index=FALSE) {
  
  if (get_node_index) {
    node <- which(tree$node.label == node) + length(tree$tip.label)
  }
  
  if (inherits(node, "character")) { stop("Node is a character and not numeric - set get_node_index=TRUE if you want to input node label rather than node number.") }
  
  node_children <- phangorn::Children(tree, node)
  
  node_lhs_descendants <- tree$tip.label[phangorn::Descendants(x = tree, node = node_children[1], type = "tips")[[1]]]
  
  node_rhs_descendants <- tree$tip.label[phangorn::Descendants(x = tree, node = node_children[2], type = "tips")[[1]]]
  
  return(list(lhs=node_lhs_descendants,
              rhs=node_rhs_descendants,
              node_i=node))
}


# Identify clusters of nodes based on how similar their underlying tips are.
# Keep a single node per cluster, which will correspond to the node with the smallest number of underlying tips.
jaccard_derep_nodes <- function(tree, node_features, cutoff = 0.75, ncores = 1) {
  
  if (length(cutoff) != 1) {
    stop("Jaccard cut-off argument must be of length one.")
  } else if ((! inherits(cutoff, "integer")) & (! inherits(cutoff, "numeric"))) {
    stop("Jaccard cut-off argument must be a number.")
  } else if (cutoff < 0 | cutoff > 1) {
    stop("Jaccard cut-off must be between 0 and 1 (inclusively).")
  }
  
  node_features_underlying <- parallel::mclapply(node_features, function(x) { c(x$lhs, x$rhs) },
                                                 mc.cores = ncores)
  
  node_combos_i <- utils::combn(length(node_features_underlying), 2, simplify = FALSE)
  
  node_jaccard <- unlist(parallel::mclapply(node_combos_i,
                                            function(x) {
                                              length(intersect(node_features_underlying[[x[1]]], node_features_underlying[[x[2]]])) / 
                                                length(union(node_features_underlying[[x[1]]], node_features_underlying[[x[2]]]))
                                            },
                                            mc.cores = ncores))
  
  node_jaccard <- data.frame(cbind(t(utils::combn(length(node_features_underlying), 2)), node_jaccard))
  
  colnames(node_jaccard) <- c("NodeA", "NodeB", "Jaccard_Sim")
  node_jaccard$NodeA <- names(node_features_underlying)[node_jaccard$NodeA]
  node_jaccard$NodeB <- names(node_features_underlying)[node_jaccard$NodeB]
  
  node_jaccard_outliers <- node_jaccard[which(node_jaccard$Jaccard_Sim >= cutoff), ]
  
  if (nrow(node_jaccard_outliers) == 0) {
    
    # No redundant nodes identified.
    return(list(node_jaccard = node_jaccard,
                ignored_redundant_nodes = as.character(),
                node_clusters = list()))

  }
  
  node_clusters <- list()
  
  unique_nodeA <- unique(node_jaccard_outliers$NodeA)
  
  cluster_num <- 1
  
  for (node_label in unique_nodeA) {
    
    matching_A <- which(node_jaccard_outliers$NodeA == node_label)
    
    if (length(matching_A) == 0) { next }  
    
    node_clusters[[cluster_num]] <- c(node_label, node_jaccard_outliers[matching_A, "NodeB"])
    node_jaccard_outliers <- node_jaccard_outliers[-matching_A, ]
    
    matching_A <- which(node_jaccard_outliers$NodeA %in% node_clusters[[cluster_num]])
    matching_B <- which(node_jaccard_outliers$NodeB %in% node_clusters[[cluster_num]])
    all_match_indices <- c(matching_A, matching_B)
    
    while(length(all_match_indices) > 0) {
      
      node_clusters[[cluster_num]] <- unique(c(node_clusters[[cluster_num]],
                                               node_jaccard_outliers[matching_A, "NodeB"],
                                               node_jaccard_outliers[matching_B, "NodeA"]))
      
      node_jaccard_outliers <- node_jaccard_outliers[-all_match_indices, ]
      
      matching_A <- which(node_jaccard_outliers$NodeA %in% node_clusters[[cluster_num]])
      matching_B <- which(node_jaccard_outliers$NodeB %in% node_clusters[[cluster_num]])
      all_match_indices <- c(matching_A, matching_B)
      
    }
    
    cluster_num <- cluster_num + 1
  }
  
  ignored_redundant_nodes <- as.character()
  
  for (cluster_i in 1:length(node_clusters)) {
    
    num_underlying_tips <- sapply(node_features_underlying[node_clusters[[cluster_i]]], length)
    
    smallest_node_i <- which.min(num_underlying_tips)
    smallest_node <- names(num_underlying_tips)[smallest_node_i]
    
    ignored_redundant_nodes <- c(ignored_redundant_nodes, node_clusters[[cluster_i]][-smallest_node_i])
    
  }

  return(list(node_jaccard = node_jaccard,
              ignored_redundant_nodes = ignored_redundant_nodes,
              node_clusters = node_clusters))

}