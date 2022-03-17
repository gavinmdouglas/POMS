#' Two-group POMS pipeline
#'
#' Key function to run POMS pipeline.\cr
#' 
#' This function will identify significant nodes based on sample balances using a Wilcoxon test by default, or significant nodes can be specified.
#' Significant nodes are referred to as Balance-Significant Nodes (BSNs).\cr
#' 
#' Fisher's exact tests are run at each node in the tree with sufficient numbers of underlying tips on each side.
#' Significant nodes based on this test are referred to as Function-Significant Nodes (FSNs).\cr
#' 
#' The key output is the tally of the intersecting nodes based on the sets of BSNs and FSNs.\cr
#' 
#' Each FSN can be categorized in one of three ways:\cr
#' (1) It does not intersect with any BSN.\cr
#' (2) It intersects with a BSN and the functional enrichment is within the taxa that are relatively more abundant in group 1 samples.\cr
#' (3) Same as #2, but enriched within taxa that are relatively more abundant in group 2 samples.\cr
#' 
#' A multinomial test is run to see if the tallies the FSNs in these three categories is significantly different from the random expectation.
#'
#' @param abun Dataframe of abundances of taxa which are at the tips of the input phylogeny, which would usually be individual genomes.
#' The taxa should be the rows and the samples the columns.
#' 
#' @param func Dataframe of the number of copies of each function that are encoded by each input taxon.
#' This pipeline only considers the presence/absence of functions across taxa.
#' Taxa (with rownames intersecting with the "abun" table) should be the rows and the functions should be the columns.
#' 
#' @param phylogeny Phylo object with tip labels that match the rownames of the "abun" and "func" tables.
#' This object is usually a newick tree that has been read into R with the ape R package.
#' 
#' @param group1_samples Character vector of column names of "abun" table that correspond to the first sample group.
#' This grouping is used for testing for significant sample balances at each node.
#' 
#' @param group2_samples Same as "group1_samples", but corresponding to the second sample group.
#' 
#' @param ncores Integer specifying how many cores to run sections of pipeline that are parallelized.
#'
#' @param pseudocount Number added to all cells of "abun" table to avoid 0 values.
#' Set this to be 0 if this is not desired, although there will be issues with the balance tree approach if any 0's are present.
#'
#' @param significant_nodes Optional vector of node names that match node labels of input phylogeny.
#' These nodes will be considered the significant balance tree set, and the Wilcoxon tests will not be run.
#' The group means of the balances at each node will still be used to determine which group has higher values.
#' Note this requires that "tested_nodes" is also specified.
#' 
#' @param tested_nodes Optional vector of node names which represent all tested nodes that resulted in the input to the significant_nodes vector.
#' I.e., this vector must include all names in the significant_nodes vector, but also all non-significant tested nodes as well.
#' 
#' @param min_num_tips The minimum number of tips on each side of a node that is required that node to be retained for the analysis.
#' Ignored if significant nodes are specified manually.
#' 
#' @param min_func_instances The minimum number of taxa that must encode the function for it to be retained for the analysis.
#'
#' @param min_func_prop The minimum proportion of taxa that must encode the function for it to be retained for the analysis.
#' 
#' @param multinomial_min_FSNs The minimum number of FSNs required to run a multinomial test for a given function.

#' @param BSN_p_cutoff Significance cut-off for identifying BSNs.
#' 
#' @param BSN_correction Multiple-test correction to use on Wilcoxon test p-values when identifying BSNs.
#' Must be a p.adjust option.
#' 
#' @param FSN_p_cutoff Significance cut-off for identifying FSNs.
#' 
#' @param FSN_correction Multiple-test correction to use on Fisher's exact test p-values when identifying FSNs.
#' Must be a p.adjust option.
#' 
#' @param func_descrip_infile Optional path to mapfile of function ids (column 1) to descriptions (column 2).
#' This should be tab-delimited with no header and one function per line.
#' If this option is specified then an additional description column will be added to the output table.
#' 
#' @param multinomial_correction Multiple-test correction to use on raw multinomial test p-values.
#' Must be a p.adjust option.
#'  
#' @param detailed_output Boolean flag to indicate that several intermediate objects should be included in the final output.
#' This is useful when troubleshooting issues, but is not expected to be useful for most users.
#'
#' @param verbose Boolean flag to indicate that log information should be written to the console, to help keep track of the pipeline's progress.
#' 
#' @return A list containing at minimum three elements:\cr\cr
#' "summary_df" - a dataframe with each tested function as a row and the numbers of FSNs of each type as columns, as well as the multinomial test output.\cr\cr
#' "balances_info" - a list of the sample balances at each node.\cr\cr
#' "sig_nodes" - the labels of BSNs.\cr\cr
#' 
#' @export
POMS_pipeline <- function(abun,
                          func,
                          phylogeny,
                          group1_samples=NULL,
                          group2_samples=NULL,
                          ncores=1,
                          pseudocount=1,
                          significant_nodes=NULL,
                          tested_nodes=NULL,
                          min_num_tips=10,
                          min_func_instances=10,
                          min_func_prop=0.001,
                          multinomial_min_FSNs=5,
                          BSN_p_cutoff = 0.05,
                          BSN_correction = "none",
                          FSN_p_cutoff = 0.05,
                          FSN_correction = "none",
                          func_descrip_infile = NULL,
                          multinomial_correction="BH",
                          detailed_output=FALSE,
                          verbose=FALSE) {
  
  if(verbose) { message("Checking input arguments.") }
  input_param <- check_POMS_pipeline_args(abun=abun,
                                          func=func,
                                          phylogeny=phylogeny,
                                          group1_samples=group1_samples,
                                          group2_samples=group2_samples,
                                          ncores=ncores,
                                          pseudocount=pseudocount,
                                          significant_nodes=significant_nodes,
                                          tested_nodes=tested_nodes,
                                          min_num_tips=min_num_tips,
                                          min_func_instances=min_func_instances,
                                          min_func_prop=min_func_prop,
                                          multinomial_min_FSNs=multinomial_min_FSNs,
                                          BSN_p_cutoff=BSN_p_cutoff,
                                          BSN_correction=BSN_correction,
                                          FSN_p_cutoff=FSN_p_cutoff,
                                          FSN_correction=FSN_correction,
                                          func_descrip_infile=func_descrip_infile,
                                          multinomial_correction=multinomial_correction,
                                          detailed_output=detailed_output,
                                          verbose=verbose)
  
  if(verbose) { message("Prepping input phylogeny.") }
  phylogeny <- prep_tree(phy=phylogeny, tips2keep=rownames(abun))
  
  if((min_func_instances > 0) || (min_func_prop > 0)) {
    func <- filter_rare_table_cols(func,  min_func_instances, min_func_prop, verbose)
  }
  
  if (is.null(significant_nodes)) {
    
    if (verbose) { message("Calculating balances.") }
    
    # Check whether pseudocount setting is reasonable.
    # This is especially important if the pseudocount is set to 1, but the data is proportional.
    min_nonzero_abun <- min(abun[abun > 0])
    if (min_nonzero_abun < pseudocount) {
       message("WARNING: specified pseudocount is larger than the smallest non-zero abundance value.")
       message("MAKE SURE TO USE A LOWER PSEUDOCOUNT VALUE IF THE INPUT TABLE IS PROPORTIONAL/PERCENTAGE DATA.")
    }
    
    calculated_balances <- compute_tree_node_balances(abun=abun,
                                                      phylogeny=phylogeny,
                                                      ncores=ncores,
                                                      min_num_tips = min_num_tips,
                                                      pseudocount=pseudocount)
    
    if (length(calculated_balances$balances) == 0) {
      message("Cannot run POMS workflow because no nodes are non-negligible based on specified settings.")
      return(list(error = "Cannot run POMS workflow because no nodes are non-negligible based on specified settings."))
    }
    
    if(verbose) { message("Identified ", length(calculated_balances$balances), " of ", length(phylogeny$node.label), " nodes as non-negligible and will be used for analyses.") }
    
    if(verbose) { message("Running Wilcoxon tests to test for differences in balances between groups at each non-negligible node.") }
    
    pairwise_node_out <- pairwise_mean_direction_and_wilcoxon(calculated_balances$balances,
                                                              group1_samples,
                                                              group2_samples,
                                                              corr_method=BSN_correction,
                                                              skip_wilcoxon=FALSE)
    
    sig_nodes <- names(calculated_balances$balances)[which(pairwise_node_out$wilcox_corrected_p < BSN_p_cutoff)]
    
  } else {
    
    sig_nodes <- significant_nodes
    
    if (any(! sig_nodes %in% phylogeny$node.label)) { stop("Not all sig. nodes are not found in phylogeny.")}
    
    calculated_balances <- list()
    calculated_balances$balances <- tested_nodes
    calculated_balances$features <- parallel::mclapply(names(tested_nodes),
                                                       lhs_rhs_asvs,
                                                       tree=phylogeny,
                                                       get_node_index=TRUE,
                                                       mc.cores=ncores)
    
    names(calculated_balances$features) <- names(tested_nodes)
    
    negligible_nodes_i <- which(! names(phylogeny$node.label) %in% names(tested_nodes))
    if (length(negligible_nodes_i) > 0) {
      calculated_balances$negligible_nodes <- phylogeny$node.label[negligible_nodes_i]
    } else {
      calculated_balances$negligible_nodes <- c()
    }
    
    pairwise_node_out <- pairwise_mean_direction_and_wilcoxon(calculated_balances$balances,
                                                              group1_samples,
                                                              group2_samples,
                                                              skip_wilcoxon=TRUE)
  }
  
  if(verbose) { message("Identifying enriched functions at all non-negligible nodes.") }
  
  all_balances_enriched_funcs <- parallel::mclapply(names(calculated_balances$balances),
                                          function(x) {
                                            return(node_func_fisher(node = x,
                                                                    in_tree = phylogeny,
                                                                    in_func = func,
                                                                    higher_group=pairwise_node_out$mean_direction[x],
                                                                    add_pseudocount=TRUE,
                                                                    multiple_test_corr=FSN_correction))
                                          },
                                          mc.cores = ncores)
  
  names(all_balances_enriched_funcs) <- names(calculated_balances$balances)
  
  if (verbose) { message("Summarizing significant functions across nodes.") }
  
  func_summaries <- summarize_node_enrichment(all_balances_enriched_funcs, sig_nodes, FSN_p_cutoff)
  
  # Get single DF summarizing the key metrics and print this out.
  all_func_id <- c()
  for (balance in names(all_balances_enriched_funcs)) {
    all_func_id <- c(all_func_id, rownames(all_balances_enriched_funcs[[balance]]))
  }
  
  if (length(which(duplicated(all_func_id))) > 0) {
    all_func_id <- all_func_id[-which(duplicated(all_func_id))]
  }
  
  if (verbose) { message("Creating results dataframe.") }

  summary_df <- data.frame(matrix(NA, nrow=length(all_func_id), ncol=4))
  
  rownames(summary_df) <- all_func_id
  
  colnames(summary_df) <- c("num_FSNs",
                            "num_FSNs_group1_enrich",
                            "num_FSNs_group2_enrich",
                            "num_FSNs_at_nonBSNs")
  
  rownames(summary_df) <- all_func_id
  
  if (verbose) { message("Creating results dataframe.") }

  if (verbose) { message("Running multinomial test on every function (that meets the multinomial_min_FSNs cut-off).") } 
  
  prop_sig_node_balances <- length(sig_nodes) / length(calculated_balances$balances)
  
  multinomial_exp_prop <- c(prop_sig_node_balances * 0.5, prop_sig_node_balances * 0.5, 1 - prop_sig_node_balances)
  
  names(multinomial_exp_prop) <- c("exp_prop_FSNs_group1_enrich",
                                   "exp_prop_FSNs_group2_enrich",
                                   "exp_prop_FSNs_at_nonBSNs")
  
  summary_df$multinomial_p <- NA
  
  for (func_id in all_func_id) {
    
      summary_df[func_id, c("num_FSNs_group1_enrich",
                            "num_FSNs_group2_enrich",
                            "num_FSNs_at_nonBSNs")] <- c(length(func_summaries[[func_id]]$positive_nodes),
                                                         length(func_summaries[[func_id]]$negative_nodes),
                                                         length(func_summaries[[func_id]]$enriched_nonsig_nodes))
      
      
      summary_df[func_id, "num_FSNs"] <- sum(as.numeric(summary_df[func_id, c("num_FSNs_group1_enrich",
                                                                              "num_FSNs_group2_enrich",
                                                                              "num_FSNs_at_nonBSNs")]))
          
      all_nodes_present <- c(func_summaries[[func_id]]$nonenriched_sig_nodes,
                             func_summaries[[func_id]]$positive_nodes,
                             func_summaries[[func_id]]$negative_nodes,
                             func_summaries[[func_id]]$nonenriched_nonsig_nodes,
                             func_summaries[[func_id]]$enriched_nonsig_nodes)
      
      if (max(table(all_nodes_present)) > 1) {
        stop("Node categorized into at least 2 mutually exclusive groups.")
      }
       
      observed_counts <- as.numeric(summary_df[func_id, c("num_FSNs_group1_enrich",
                                                          "num_FSNs_group2_enrich",
                                                          "num_FSNs_at_nonBSNs")])
      
      if ((length(sig_nodes) > 0) && (summary_df[func_id, "num_FSNs"] >= multinomial_min_FSNs) && (prop_sig_node_balances != 1)) {
        summary_df[func_id, "multinomial_p"] <- XNomial::xmulti(obs=observed_counts,
                                                                expr=multinomial_exp_prop, detail=0)$pProb 
      }
    }

    if (multinomial_correction != "none") {
      summary_df$multinomial_corr <- p.adjust(summary_df$multinomial_p, multinomial_correction)
    }
  
    if (! is.null(func_descrip_infile)) {
    if (verbose) { message("Adding function descriptions to output.") }
    func_descrip <- read.table(func_descrip_infile,
                               header=FALSE, sep="\t", row.names=1, stringsAsFactors = FALSE, quote="")
    summary_df$description <- func_descrip[rownames(summary_df), "V2"]
  } else {
    if (verbose) { message("Function description mapfile not specified (func_descrip_infile argument), so no descriptions will be added.") } 
  }
  
  results <- list(balances_info=calculated_balances,
                  sig_nodes=sig_nodes,
                  df=summary_df)
  
  results[["multinomial_exp_prop"]] <- multinomial_exp_prop
  
  if (detailed_output) {
      # Restrict vector of mean directions to significant nodes only to avoid confusion.
      pairwise_node_out$mean_direction <- pairwise_node_out$mean_direction[sig_nodes]

      results[["balance_comparisons"]] <- pairwise_node_out
      results[["funcs_per_node"]] <- all_balances_enriched_funcs
      results[["out_list"]] <- func_summaries
      results[["tree"]] <- phylogeny
      results[["input_param"]] <- input_param
  }
  
  return(results)

}


check_POMS_pipeline_args <- function(abun,
                                     func,
                                     phylogeny,
                                     group1_samples,
                                     group2_samples,
                                     ncores,
                                     pseudocount,
                                     significant_nodes,
                                     tested_nodes,
                                     min_num_tips,
                                     min_func_instances,
                                     min_func_prop,
                                     multinomial_min_FSNs,
                                     BSN_p_cutoff,
                                     BSN_correction,
                                     FSN_p_cutoff,
                                     FSN_correction,
                                     func_descrip_infile,
                                     multinomial_correction,
                                     detailed_output,
                                     verbose) {

  if (((is.null(significant_nodes)) && (! is.null(tested_nodes))) || ((! is.null(significant_nodes)) && (is.null(tested_nodes)))) {
    stop("Stopping - arguments significant_nodes and tested_nodes either both need to be given or neither should be specified.")
  }

  if ((! is.null(significant_nodes)) && (length(significant_nodes) == 0)) { stop("Stopping - vector specified for significant_nodes argument is empty.") }

  if ((! is.null(significant_nodes)) && (! is.null(tested_nodes))) {
    if (length(which(significant_nodes %in% tested_nodes)) != length(significant_nodes)) {
      stop("Stopping - not all nodes in significant_nodes vector are in in tested_nodes") 
    }
  }
  
  if (class(abun) != "data.frame") { stop("Stopping - argument abun needs to be of the class data.frame.") }
  if (class(func) != "data.frame") { stop("Stopping - argument func needs to be of the class data.frame.") }
  if (class(phylogeny) != "phylo") { stop("Stopping - argument phylo needs to be of the class phylo.") }

  if (class(group1_samples) != "character") { stop("Stopping - argument group1_samples needs to be of the class character.") }
  if (class(group2_samples) != "character") { stop("Stopping - argument group2_samples needs to be of the class character.") }
  if (length(group1_samples) == 0) { stop("Stopping - argument group1_samples is of length 0.") }
  if (length(group2_samples) == 0) { stop("Stopping - argument group2_samples is of length 0.") }
  if (length(which(group1_samples %in% colnames(abun))) != length(group1_samples)) { stop("Stopping - not all group1_samples match columns of abun argument.") }
  if (length(which(group2_samples %in% colnames(abun))) != length(group2_samples)) { stop("Stopping - not all group2_samples match columns of abun argument.") }
  if (length(which(group1_samples %in% group2_samples)) > 0) { stop("Stopping - at least one sample overlaps between group1 and group2.") }
  
  if ((class(ncores) != "integer") && (class(ncores) != "numeric")) { stop("Stopping - ncores argument needs to be of class numeric or integer.") }
  if (ncores <= 0) { stop("Stopping - ncores argument needs to be higher than 0.") }

  if ((class(pseudocount) != "integer") && (class(pseudocount) != "numeric")) { stop("Stopping - pseudocount argument needs to be of class numeric or integer.") }
  if (pseudocount < 0) { stop("Stopping - pseudocount argument cannot be lower than 0.") }

  if ((class(multinomial_min_FSNs) != "integer") && (class(multinomial_min_FSNs) != "numeric")) { stop("Stopping - multinomial_min_FSNs argument needs to be of class numeric or integer.") }
  if (multinomial_min_FSNs < 0) { stop("Stopping - multinomial_min_FSNs argument cannot be lower than 0.") }
  
  if (min_num_tips > length(phylogeny$tip.label) / 2) { stop("Stopping - the min_num_tips argument cannot be higher than half of the total number of tips.") }

  if ((min_func_prop < 0) || (min_func_prop > 1)) { stop("Stopping - the min_func_prop argument must be between 0 and 1.") }
  if ((BSN_p_cutoff < 0) || (BSN_p_cutoff > 1)) { stop("Stopping - the BSN_p_cutoff argument must be between 0 and 1.") }
  if ((FSN_p_cutoff < 0) || (FSN_p_cutoff > 1)) { stop("Stopping - the FSN_p_cutoff argument must be between 0 and 1.") }
  
  if ((min_func_instances < 0) || (min_func_instances > ncol(func))) { stop("Stopping - the min_func_instances argument must be between 0 and 1.") }

  if (! FSN_correction %in% p.adjust.methods) { stop("Stopping - FSN_correction argument needs to be found in p.adjust.methods.") }
  if (! BSN_correction %in% p.adjust.methods) { stop("Stopping - BSN_correction argument needs to be found in p.adjust.methods.") }
  if (! multinomial_correction %in% p.adjust.methods) { stop("Stopping - multinomial_correction argument needs to be found in p.adjust.methods.") }
  
  if (! is.null(func_descrip_infile) && (! file.exists(func_descrip_infile))) { stop("Stopping - func_descrip_infile is non-NULL, but the specified file was not found.") }

  if (! is.logical(detailed_output)) { stop("Stopping - detailed_output argument needs to be TRUE or FALSE.") }
  
  if (! is.logical(verbose)) { stop("Stopping - verbose argument needs to be TRUE or FALSE.") }

  return(list(group1_samples=group1_samples,
              group2_samples=group2_samples,
              ncores=ncores,
              pseudocount=pseudocount,
              significant_nodes=significant_nodes,
              min_num_tips=min_num_tips,
              min_func_instances=min_func_instances,
              min_func_prop=min_func_prop,
              multinomial_min_FSNs=multinomial_min_FSNs,
              BSN_p_cutoff=BSN_p_cutoff,
              BSN_correction=BSN_correction,
              FSN_p_cutoff=FSN_p_cutoff,
              FSN_correction=FSN_correction,
              func_descrip_infile=func_descrip_infile,
              multinomial_correction=multinomial_correction,
              detailed_output=detailed_output,
              verbose=verbose))
}


summarize_node_enrichment <- function(enriched_funcs, sig_nodes, func_p_cutoff) {
  # For each function that is significant at least once get:
  #   - Names of significant nodes where the function is present.
  #   - Names of positively and negatively enriched significant nodes
  #   - Names of nonsignificant nodes where the function would have been significant.
  #   - Names of nonsignificant nodes where the function was present.

  # First get all unique functions.
  all_func_id <- c()
  for (node in names(enriched_funcs)) {
    all_func_id <- c(all_func_id, rownames(enriched_funcs[[node]]))
  }

  if (length(which(duplicated(all_func_id))) > 0) {
    all_func_id <- all_func_id[-which(duplicated(all_func_id))]
  }

  func_summaries <- list()

  # Loop through each function and get breakdown of contributing node names.
  for (func_id in all_func_id) {

    func_summaries[[func_id]]$positive_nodes <- c()
    func_summaries[[func_id]]$negative_nodes <- c()
    func_summaries[[func_id]]$nonenriched_sig_nodes <- c()
    func_summaries[[func_id]]$nonenriched_nonsig_nodes <- c()
    func_summaries[[func_id]]$enriched_nonsig_nodes <- c()
    func_summaries[[func_id]]$missing_sig_nodes <- c()
    func_summaries[[func_id]]$missing_nonsig_nodes <- c()

    for (node in names(enriched_funcs)) {
      if (func_id %in% rownames(enriched_funcs[[node]])) {

        if (node %in% sig_nodes) {
          # If this node was significant then categorize it as either positive, negative, or nonenriched.
          if (as.numeric(enriched_funcs[[node]][func_id, "P_corr"]) < func_p_cutoff) {

            if (as.numeric(enriched_funcs[[node]][func_id, "OR"]) > 1) {
              func_summaries[[func_id]]$positive_nodes <- c(func_summaries[[func_id]]$positive_nodes, node)
            } else if (as.numeric(enriched_funcs[[node]][func_id, "OR"]) < 1) {
              func_summaries[[func_id]]$negative_nodes <- c(func_summaries[[func_id]]$negative_nodes, node)
            } else {
              print(enriched_funcs[[node]][func_id, ])
              stop("Significant function but OR for above info not different from 1?!")
            }
          } else {
            func_summaries[[func_id]]$nonenriched_sig_nodes <- c(func_summaries[[func_id]]$nonenriched_sig_nodes, node)
          }
        } else {
          # Since this node was NOT significant then categorize it as either enriched or nonenriched.
          if (as.numeric(enriched_funcs[[node]][func_id, "P_corr"]) < func_p_cutoff) {
            func_summaries[[func_id]]$enriched_nonsig_nodes <- c(func_summaries[[func_id]]$enriched_nonsig_nodes, node)
          } else {
            func_summaries[[func_id]]$nonenriched_nonsig_nodes <- c(func_summaries[[func_id]]$nonenriched_nonsig_nodes, node)
          }
        }
      } else {
        if (node %in% sig_nodes) {
          func_summaries[[func_id]]$missing_sig_nodes <- c(func_summaries[[func_id]]$missing_sig_nodes, node)
        } else {
          func_summaries[[func_id]]$missing_nonsig_nodes <- c(func_summaries[[func_id]]$missing_nonsig_nodes, node)
        }
      }
    }
  }
  return(func_summaries)
}
