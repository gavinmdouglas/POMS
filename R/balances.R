#' Determine representative taxa labels of tips on each side of a specific node.\cr
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
#' then Pseudomonas will be used as the representative label (based on a threshold of 0.75 or lower and assuming that species are the last column in the table).
#'
#' @param in_tree Phylo object
#' 
#' @param taxon_labels Dataframe of taxa labels for tips in tree.
#' All tips underlying the specified node must be present, although typically all tips in the tree would be present.
#' Rownames must be the tip labels. The column names correspond to each taxonomic level, such as Kingdom, Phylum, etc.
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
#' E.g., rather than just the genus "Odoribacter" the label would be "Bacteria; Bacteroidetes; Bacteroidia; Bacteroidales; Porphyromonadaceae; Odoribacter".
#'
#' @return Character string of representative taxon of tips underlying the specified node, delimited by " / ".
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
    stop("Stopping - one of node_index and node_label must be specified.")
  }
  
  if (length(which(! c(underlying_tips$lhs, underlying_tips$rhs) %in% rownames(taxon_labels))) > 0) {
    message("The following tips are missing as rownames from the taxa label table.")
    message(c(underlying_tips$lhs, underlying_tips$rhs)[which(! c(underlying_tips$lhs, underlying_tips$rhs) %in% rownames(taxon_labels))])
    stop("Stopping - please make sure the rownames of the input taxa label table are the tip labels of the tree.")
  }

  # Determine actual feature consensus step:
  taxon_lhs <- feature_consensus_taxon(taxon_labels, underlying_tips$lhs, threshold, combine_labels)
  taxon_rhs <- feature_consensus_taxon(taxon_labels, underlying_tips$rhs, threshold, combine_labels)
  
  return(paste(taxon_lhs, taxon_rhs, sep=" / "))

}


#' @export
calc_balances <- function(abun_table, lhs_features, rhs_features, pseudocount=NULL) {
  
  if (pseudocount) {
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


feature_consensus_taxon <- function(taxa, features, threshold, combine_labels = FALSE) {

  taxa_subset <- taxa[features, ]
  
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
  
  # Get tips on either side of each node.
  node_features <- parallel::mclapply(nodes2test,
                                      lhs_rhs_tips,
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
    balance_calc <- parallel::mclapply(nonnegligible_nodes,
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


feature_consensus_taxon_OLD <- function(taxa, features, threshold) {
  
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


lhs_rhs_tips <- function(tree, node, get_node_index=FALSE) {
  
  if(get_node_index) {
    node <- which(tree$node.label == node) + length(tree$tip.label)
  }
  
  if(class(node) == "character") { stop("Node is a character and not numeric - set get_node_index=TRUE if you want to input node label rather than node number.") }
  
  node_children <- phangorn::Children(tree, node)
  
  node_lhs_descendants <- tree$tip.label[phangorn::Descendants(x = tree, node = node_children[1], type = "tips")[[1]]]
  
  node_rhs_descendants <- tree$tip.label[phangorn::Descendants(x = tree, node = node_children[2], type = "tips")[[1]]]
  
  return(list(lhs=node_lhs_descendants,
              rhs=node_rhs_descendants,
              node_i=node))
}

