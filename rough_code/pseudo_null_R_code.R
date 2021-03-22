#' @export
POM_pseudo_null <- function(num_sig_nodes,
                            funcs_per_node,
                            subset_to_assess=NULL,
                            ncores=1,
                            num_null_rep=1000,
                            p_corr_cutoff=0.05,
                            ran_seed=6276) {
  
  set.seed(ran_seed)
  
  if(is.null(subset_to_assess)) {
    repeated_func <- as.character(do.call(c, sapply(funcs_per_node, rownames)))
    subset_to_assess <- repeated_func[-which(duplicated(repeated_func))]
  }
  
  # Generate profile of enrichments per node to sample from.
  node_enrichments <- mclapply(funcs_per_node,
                               compute_node_enrichments,
                               functions_to_include=subset_to_assess,
                               p_corr_cutoff=p_corr_cutoff,
                               mc.cores = ncores)
  
  names(node_enrichments) <- names(funcs_per_node)
  
  # For the specificed no. replicates, sum up the no. enrichments of each type per gene given a sampled set of nodes.
  abs_enrich_out <- mclapply(1:num_null_rep,
                             function(x) { return(random_nodes_enrichment(possible_enrichments=node_enrichments,
                                                                          num_sig_nodes=num_sig_nodes)) },
                             mc.cores = ncores)
  
  abs_enrich_out <- as.data.frame(do.call(cbind, abs_enrich_out))
  
  rownames(abs_enrich_out) <- subset_to_assess
  
  return(abs_enrich_out)
}



#' @export
pseudo_null_pvalues <- function(pseudo_null_out, actual_df) {
  
  func_pseudo_pvalues <- c()
  
  for(func in rownames(pseudo_null_out)) {
    
    abs_enrich <- abs(actual_df[func, "num_sig_nodes_pos_enrich"] - actual_df[func, "num_sig_nodes_neg_enrich"])
    
    func_pseudo_pvalues <- c(func_pseudo_pvalues, (length(which(pseudo_null_out[func, ] >= abs_enrich)) / ncol(pseudo_null_out)))
    
  }
  
  names(func_pseudo_pvalues) <- rownames(pseudo_null_out)
  
  return(func_pseudo_pvalues)
}

#' @export
compute_node_enrichments <- function(node_fishers_out, functions_to_include, p_corr_cutoff=0.05) {
  
  non_included_genes <- which(! rownames(node_fishers_out) %in% functions_to_include)
  if(length(non_included_genes) > 0) {
    node_fishers_out <- node_fishers_out[-non_included_genes, ]
  }
  
  node_fishers_out$pos_enrich <- 0
  node_fishers_out$neg_enrich <- 0
  
  for(gene in rownames(node_fishers_out)) {
    if(node_fishers_out[gene, "P_corr"] < p_corr_cutoff) {
      if(node_fishers_out[gene, "OR"] > 1) {
        node_fishers_out[gene, "pos_enrich"] <- 1
      }  else if(node_fishers_out[gene, "OR"] < 1) {
        node_fishers_out[gene, "neg_enrich"] <- 1
      }
    }
  }
  
  missing_genes <- functions_to_include[which(! functions_to_include %in% rownames(node_fishers_out))]
  
  if(length(missing_genes) > 0) {
    node_fishers_out[missing_genes, ] <- 0
  }
  
  return(node_fishers_out[functions_to_include, c("pos_enrich", "neg_enrich"), drop=FALSE])
  
}

#' @export
random_nodes_enrichment <- function(possible_enrichments, num_sig_nodes) {
  
  observed_enrichments <- possible_enrichments[sample(names(possible_enrichments), num_sig_nodes)]
  
  observed_enrichments <- base::Reduce("+", observed_enrichments)
  
  observed_enrichments$abs_enrich <- abs(observed_enrichments$pos_enrich - observed_enrichments$neg_enrich)
  
  return(observed_enrichments$abs_enrich)
  
}
