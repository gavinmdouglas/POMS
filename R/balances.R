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
  
  return(c(taxon_lhs, taxon_rhs))

}


#' Computes isometric log ratio between two sets of feature abundances, for each sample separately. Requires an abundance table,
#' with two sets of features for which the ratio will be computed. 
#'
#' @param abun_table Abundance table, e.g., read counts or relative abundance.
#' Should be dataframe with column names correspond to sample names and row names corresponding to the feature ids.
#' No 0 values are permitted unless the "pseudocount" option is set.
#'
#' @param set1_features Features (rows of abundance table) that make up one side of the ratio to be computed (numerator).
#' 
#' @param set2_features Same as "set1_features", but for the other side of the ratio (denominator).
#' 
#' @param pseudocount Constant to add to all abundance values, to ensure that there are only non-zero values. For read count data this would typically be 1.
#' 
#' @return Numeric vector of the computed isometric log ratio for each sample (taken to be each column in the input table).
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
  
  # Check that no features intersect.
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


#' Computes balances (i.e., isometric log ratios, for each sample separately) of feature abundances at each node in the tree.
#' 
#' @param phylogeny Phylo object with tip labels matching row names of input abundance table. Note that node labels are required.
#'
#' @param abun_table Abundance table, e.g., read counts or relative abundance.
#' Should be dataframe with column names correspond to sample names and row names corresponding to the tips of the tree.
#' No 0 values are permitted unless the "pseudocount" option is set.
#' 
#' @param min_num_tips Minimum number of tips that must be found on each side of a node for it to be included.
#'
#' @param ncores Number of cores to use for steps of function that can be run in parallel.
#'
#' @param pseudocount Constant to add to all abundance values, to ensure that there are only non-zero values. For read count data this would typically be 1.
#'
#' @param subset_to_test Vector of node labels (*not indices*) that correspond to the subset of nodes that should be considered.
#' Note that balances will still only be computed at each of these nodes if they have a sufficient number of underlying tips (as specified by the "min_num_tips" argument).
#'
#' @return List containing three items:
#' 
#' "tips_underlying_nodes": the tips on the left-hand side (lhs; the numerator) and right-hand side (rhs; the denominator) of each node.
#' Note that which side of the node is denoted as the left-hand or right-hand side is arbitrary.
#' 
#' "balances": list with each non-negligible nodes as a separate element. The sample balances for each node are provided as a numeric vector within each of these elements.
#'  
#' "negligible_nodes": character vector of node labels where there are fewer tips on either side of the node than specified by the "min_num_tips" argument.
#' 
#' @export
compute_node_balances <- function(phylogeny, abun_table, min_num_tips=5, ncores=1, pseudocount=NULL, subset_to_test=NULL) {
  
  if (is.null(phylogeny$node.label)) {
    stop("Stopping - input tree does not have any node labels.") 
  }
  
  if (! all(phylogeny$tip.label %in% rownames(abun_table))) {
    stop("Stopping - not all tips are found as row names in the abundance table.")
  }
  
  # Test all nodes unless subset specified.  
  if(! is.null(subset_to_test)) {
    if (length(which(! subset_to_test %in% phylogeny$node.label)) > 0) {
      stop("Stopping - some labels in subset_to_test do not match node labels in the tree.")  
    }
    nodes_to_test <- subset_to_test
  } else {
    nodes_to_test <- phylogeny$node.label
  }
  
  # Get tips on either side of each node.
  node_features <- parallel::mclapply(nodes_to_test,
                                      lhs_rhs_tips,
                                      tree=phylogeny,
                                      get_node_index=TRUE,
                                      mc.cores=ncores)
  
  names(node_features) <- nodes_to_test
  
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
  
  nonnegligible_nodes_i <- which(! names(node_features) %in% negligible_nodes)
  
  if (length(nonnegligible_nodes_i) > 0) {
    nonnegligible_nodes <- names(node_features)[nonnegligible_nodes_i]
    
    # Calculate balances at each node.
    balance_calc <- parallel::mclapply(nonnegligible_nodes,
                                       function(x) {
                                         return(abun_isometric_log_ratios(abun_table=abun_table,
                                                                          set1_features=node_features[[x]]$lhs,
                                                                          set2_features=node_features[[x]]$rhs,
                                                                          pseudocount=pseudocount))
                                       },
                                       mc.cores=ncores)
    
    names(balance_calc) <- nonnegligible_nodes
    
  } else {
    stop("Stopping - no non-negligible nodes remain after filtering based on mininum number of tips of left and right-hand side of each node.")
  }
  
  return(list(tips_underlying_nodes=node_features, balances=balance_calc, negligible_nodes=negligible_nodes))
  
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

