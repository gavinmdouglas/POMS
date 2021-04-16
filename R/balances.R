#' @export
calc_balances <- function(abun_table, lhs_features, rhs_features, pseudocount=NULL) {
  
  if(pseudocount) {
    abun_table <- abun_table + pseudocount
  }
  
  num_lhs <- length(lhs_features)
  num_rhs <- length(rhs_features)
  
  half1 <- sqrt((num_lhs * num_rhs) / (num_lhs + num_rhs))
  
  sample_col <- colnames(abun_table)
  
  half2 <- c()
  
  for(sample in sample_col) {
    sample_lhs_values <- abun_table[lhs_features, sample]
    sample_rhs_values <- abun_table[rhs_features, sample]
    
    gmean_lhs <- exp(mean(log(sample_lhs_values)))
    gmean_rhs <- exp(mean(log(sample_rhs_values)))
    
    half2 <- c(half2, log(gmean_lhs / gmean_rhs))
  }
  
  balances_out <- half1 * half2
  
  names(balances_out) <- sample_col
  
  return(balances_out)
}

#' @export
compute_tree_node_balances <- function(phylogeny, abun, min_num_tips, ncores=1, pseudocount=1, subset2test=NULL) {
  
  ### Function to perform isomatric log-ratio transformation of feature abundances at each node in the tree.
  ### Will return a list containing the features on the left-hand side (lhs) and right-hand side (rhs) of each
  ### node in a tree ("features") and also the computed balances ("balances").
  
  if(is.null(phylogeny$node.label)) {
    stop("Stopping - input tree does not have any node labels.") 
  }
  
  # Test all nodes unless subset specified.  
  if(! is.null(subset2test)) {
    nodes2test <- subset2test
  } else {
    nodes2test <- phylogeny$node.label
  }
  
  # Get ASVs on either side of each node.
  node_features <- parallel::mclapply(nodes2test,
                                      lhs_rhs_asvs,
                                      tree=phylogeny,
                                      get_node_index=TRUE,
                                      mc.cores=ncores)
  
  names(node_features) <- nodes2test
  
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
    negligible_nodes <- c()
  }
  
  nonnegligible_nodes_i <- which(! names(node_features) %in% negligible_nodes)
  
  if(length(nonnegligible_nodes_i) > 0) {
    nonnegligible_nodes <- names(node_features)[nonnegligible_nodes_i]
    
    # Calculate balances at each node.
    balance_calc <- mclapply(nonnegligible_nodes,
                             function(x) {
                               return(calc_balances(abun_table=abun,
                                                    lhs_features=node_features[[x]]$lhs,
                                                    rhs_features=node_features[[x]]$rhs,
                                                    pseudocount=pseudocount))
                             },
                             mc.cores=ncores)
    
    names(balance_calc) <- nonnegligible_nodes

  } else {
    balance_calc <- NULL
    message("Skipping balance calculation, because no non-negligible nodes remain after filtering based on mininum number of tips of left and right-hand side of each node.")
  }
  
  return(list(features=node_features, balances=balance_calc, negligible_nodes=negligible_nodes))
  
}


#' @export
node_taxa <- function(lhs_features, rhs_features, taxa, threshold=0.75) {
  
  if(threshold <= 0.5 | threshold > 1) {
    stop("The set threshold needs to be > 0.5 and <= 1.")
  }
  
  taxon_lhs <- feature_consensus_taxon(taxa, lhs_features, threshold)
  taxon_rhs <- feature_consensus_taxon(taxa, rhs_features, threshold)
  
  return(paste(taxon_lhs, taxon_rhs, sep=" / "))
}

feature_consensus_taxon <- function(taxa, features, threshold) {
  
  taxa_subset <- taxa[features, ]
  
  for(i in ncol(taxa_subset):1) {
    
    taxa_subset_table <- table(taxa_subset[, i, drop=TRUE])
    
    if(length(taxa_subset_table) == 0) { next }
    
    if((max(taxa_subset_table) / length(features)) >= threshold) {
      return(paste(names(taxa_subset_table)[which(taxa_subset_table == max(taxa_subset_table))], " ", "(", colnames(taxa_subset)[i], ")", sep=""))
    }
  }
  
  return("Unclear")
}

#' @export
internode_mean_max_dist <- function(phy, dist_matrix, node_labels) {
  
  # If fewer than 2 nodes in set then return NA.
  if(length(node_labels) <= 1) {
    missing_val <- c(NA, NA)
    names(missing_val) <- c("mean", "max")
    return(missing_val)
  }
  
  nodes_i <- which(phy$node.label %in% node_labels) + length(phy$tip.label)
  
  if(length(nodes_i) != length(node_labels)) {
    stop(paste("Wrong number of node labels match input set. Looking for ",
               paste(node_labels, collapse = ", "),
               " and found indices ",
               paste(as.character(nodes_i), collapse=", "),
               ".",
               sep=""))
  }
  
  dist_matrix_subset <- as.dist(dist_matrix[nodes_i, nodes_i])
  
  internode_dist_metrics <- c(mean(dist_matrix_subset), max(dist_matrix_subset))
  names(internode_dist_metrics) <- c("mean", "max")
  
  return(internode_dist_metrics)
  
}


#' @export
lhs_rhs_asvs <- function(tree, node, get_node_index=FALSE) {
  
  if(get_node_index) {
    node <- which(tree$node.label == node) + length(tree$tip.label)
  }
  
  if(class(node) == "character") { stop("Node is a character and not numeric - set get_node_index=TRUE if you want to input node label rather than node number.") }
  
  node_children <- phangorn::Children(tree, node)
  
  node_lhs_descendants <- tree$tip.label[phangorn::Descendants(x = tree, node = node_children[1], type = "tips")[[1]]]
  
  node_rhs_descendants <- tree$tip.label[phangorn::Descendants(x = tree, node = node_children[2], type = "tips")[[1]]]
  
  return(list(lhs=node_lhs_descendants,
              rhs=node_rhs_descendants,
              count=c(length(node_lhs_descendants),
                      length(node_rhs_descendants)),
              node_i=node))
}

