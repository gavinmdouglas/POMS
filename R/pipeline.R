#' Key function to run POMS pipeline.
#'
#' See details below.\cr
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
#' @param abun Dataframe of abundances of taxa which are at the tips of the input tree, which would usually be individual genomes.
#' The taxa should be the rows and the samples the columns.
#' 
#' @param func Dataframe of the number of copies of each function that are encoded by each input taxon.
#' This pipeline only considers the presence/absence of functions across taxa.
#' Taxa (with rownames intersecting with the "abun" table) should be the rows and the functions should be the columns.
#' 
#' @param tree Phylo object with tip labels that match the rownames of the "abun" and "func" tables.
#' This object is usually a newick tree that has been read into R with the ape R package.
#' 
#' @param group1_samples Character vector of column names of "abun" table that correspond to the first sample group.
#' This grouping is used for testing for significant sample balances at each node.
#' Required unless the "manual_BSN_dir" argument is set (i.e., if the binary directions of BSNs are specified manually).
#' 
#' @param group2_samples Same as "group1_samples", but corresponding to the second sample group.
#' 
#' @param ncores Integer specifying how many cores to run sections of pipeline that are parallelized.
#'
#' @param pseudocount Number added to all cells of "abun" table to avoid 0 values.
#' Set this to be 0 if this is not desired, although there will be issues with the balance tree approach if any 0's are present.
#'
#' @param manual_BSNs Optional vector of node names that match node labels of input tree.
#' These nodes will be considered the significant balance tree set, and the Wilcoxon tests will not be run.
#' The group means of the b/alances at each node will still be used to determine which group has higher values.
#' Note this requires that/ "manual_balances" is also specified.
#' 
#' @param manual_balances Optional list of balance values which represent the balances at all tested nodes that resulted in the input to the manual_BSNs vector.
#' This list must include balances for all nodes in the manual_BSNs vector, but also all non-significant tested nodes as well.
#' These node labels must all be present in the input tree.
#' The required list format is the "balances" object in the output of compute_node_balances (where balance values could be replaced by another approach if needed).
#' 
#' @param manual_BSN_dir Optional character vector specifying "group1" or "group2", depending on the direction of the BSN difference.
#' This must be a named vector, with all names matching the set of nodes specified by the manual_BSNs argument.
#' Although this requires that the exact labels "group1" or "group2" are specified, these categories could represent different binary divisions rather than strict sample groups. 
#' For instance, "group1" could be used to represent nodes where sample balances are positively associated with a continuous variable (rather than a discrete grouping),
#' whereas "group2" could represent nodes where sample balances are negatively associated.
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
#' The additional results include "balance_comparisons" (summary of Wilcoxon tests on balances),
#' "func_enrichments" (Fisher's exact test output for all functions at each node), and
#' "input_param" (a list containing the input parameters specified).
#'
#' @param verbose Boolean flag to indicate that log information should be written to the console, to help keep track of the pipeline's progress.
#' 
#' @return A list containing at minimum these elements:\cr\cr
#' - results: dataframe with each tested function as a row and the numbers of FSNs of each type as columns, as well as the multinomial test output.\cr\cr
#' - balances: list of the sample balances at each tested node (including non-significant nodes).\cr\cr
#' - BSNs: character vector with BSNs as names and values of \"group1\" and \"group2\" to indicate for which sample group (or other binary division) the sample balances were higher\cr\cr
#' - FSNs_summary: list containing each tested function as a separate element. For functions with FSNs, will provide the node labels for nodes in each category of the multinomial test.
#' - tree: the prepped tree used by the pipeline, including the added node labels if a tree lacking labels was provided.

#' @export
POMS_pipeline <- function(abun,
                          func,
                          tree,
                          group1_samples=NULL,
                          group2_samples=NULL,
                          ncores=1,
                          pseudocount=1,
                          manual_BSNs=NULL,
                          manual_balances=NULL,
                          manual_BSN_dir=NULL,
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
                                          tree=tree,
                                          group1_samples=group1_samples,
                                          group2_samples=group2_samples,
                                          ncores=ncores,
                                          pseudocount=pseudocount,
                                          manual_BSNs=manual_BSNs,
                                          manual_balances=manual_balances,
                                          manual_BSN_dir=manual_BSN_dir,
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
  
  if(verbose) { message("Prepping input tree.") }
  tree <- prep_tree(phy=tree, tips2keep=rownames(abun))
  
  if((min_func_instances > 0) || (min_func_prop > 0)) {

    func <- filter_rare_table_cols(in_tab = func,
                                   min_nonzero_count = min_func_instances,
                                   min_nonzero_prop = min_func_prop,
                                   drop_missing_rows = FALSE,
                                   verbose = verbose)
  }
  
  if (is.null(manual_BSNs)) {
    
    if (verbose) { message("Calculating balances.") }
    
    # Check whether pseudocount setting is reasonable.
    # This is especially important if the pseudocount is set to 1, but the data is proportional.
    min_nonzero_abun <- min(abun[abun > 0])
    if (min_nonzero_abun < pseudocount) {
       message("WARNING: specified pseudocount is larger than the smallest non-zero abundance value.")
       message("MAKE SURE TO USE A LOWER PSEUDOCOUNT VALUE IF THE INPUT TABLE IS PROPORTIONAL/PERCENTAGE DATA.")
    }
    
    calculated_balances <- compute_node_balances(abun=abun,
                                                 tree=tree,
                                                 ncores=ncores,
                                                 min_num_tips = min_num_tips,
                                                 pseudocount=pseudocount)
    
    if (length(calculated_balances$balances) == 0) {
      message("Cannot run POMS workflow because no nodes are non-negligible based on specified settings.")
      return(list(error = "Cannot run POMS workflow because no nodes are non-negligible based on specified settings."))
    }
    
    if(verbose) { message("Identified ", length(calculated_balances$balances), " of ", length(tree$node.label), " nodes as non-negligible and will be used for analyses.") }
    
    if(verbose) { message("Running Wilcoxon tests to test for differences in balances between groups at each non-negligible node.") }
    
    pairwise_node_out <- pairwise_mean_direction_and_wilcoxon(calculated_balances$balances,
                                                              group1_samples,
                                                              group2_samples,
                                                              corr_method=BSN_correction,
                                                              skip_wilcoxon=FALSE)
    
    BSNs <- names(calculated_balances$balances)[which(pairwise_node_out$wilcox_corrected_p < BSN_p_cutoff)]
    
  } else {
    
    BSNs <- manual_BSNs
    
    if (any(! BSNs %in% tree$node.label)) { stop("Not all sig. nodes are not found in tree.")}
    
    calculated_balances <- list()
    calculated_balances$balances <- manual_balances
    calculated_balances$features <- parallel::mclapply(names(manual_balances),
                                                       lhs_rhs_tips,
                                                       tree=tree,
                                                       get_node_index=TRUE,
                                                       mc.cores=ncores)
    
    names(calculated_balances$features) <- names(manual_balances)
    
    negligible_nodes_i <- which(! names(tree$node.label) %in% names(manual_balances))
    if (length(negligible_nodes_i) > 0) {
      calculated_balances$negligible_nodes <- tree$node.label[negligible_nodes_i]
    } else {
      calculated_balances$negligible_nodes <- c()
    }
    
    if (! is.null(manual_BSN_dir)) {
      pairwise_node_out <- list(mean_direction = manual_BSN_dir)
    } else {
      pairwise_node_out <- pairwise_mean_direction_and_wilcoxon(calculated_balances$balances,
                                                                group1_samples,
                                                                group2_samples,
                                                                corr_method=BSN_correction,
                                                                skip_wilcoxon=TRUE)
    }

  }
  
  if (verbose) { message("Identifying enriched functions at all non-negligible nodes.") }
  
  all_node_enriched_funcs <- parallel::mclapply(names(calculated_balances$balances),
                                          function(x) {
                                            return(node_func_fisher(node = x,
                                                                    in_tree = tree,
                                                                    in_func = func,
                                                                    higher_group=pairwise_node_out$mean_direction[x],
                                                                    add_pseudocount=FALSE,
                                                                    multiple_test_corr=FSN_correction))
                                          },
                                          mc.cores = ncores)
  
  names(all_node_enriched_funcs) <- names(calculated_balances$balances)
  
  if (verbose) { message("Summarizing significant functions across nodes.") }
  
  func_summaries <- summarize_node_enrichment(all_node_enriched_funcs, BSNs, FSN_p_cutoff)
  
  # Get single DF summarizing the key metrics and print this out.
  all_func_id <- c()
  for (balance in names(all_node_enriched_funcs)) {
    all_func_id <- c(all_func_id, rownames(all_node_enriched_funcs[[balance]]))
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

  if (verbose) { message("Running multinomial test on every function (that meets the multinomial_min_FSNs cut-off).") } 
  
  prop_tested_BSNs <- length(BSNs) / length(calculated_balances$balances)
  
  multinomial_exp_prop <- c(prop_tested_BSNs * 0.5, prop_tested_BSNs * 0.5, 1 - prop_tested_BSNs)
  
  names(multinomial_exp_prop) <- c("exp_prop_FSNs_group1_enrich",
                                   "exp_prop_FSNs_group2_enrich",
                                   "exp_prop_FSNs_at_nonBSNs")
  
  summary_df$multinomial_p <- NA
  
  for (func_id in all_func_id) {
    
      summary_df[func_id, c("num_FSNs_group1_enrich",
                            "num_FSNs_group2_enrich",
                            "num_FSNs_at_nonBSNs")] <- c(length(func_summaries[[func_id]]$FSNs_group1_enrich),
                                                         length(func_summaries[[func_id]]$FSNs_group2_enrich),
                                                         length(func_summaries[[func_id]]$FSNs_at_nonBSNs))
      
      
      summary_df[func_id, "num_FSNs"] <- sum(as.numeric(summary_df[func_id, c("num_FSNs_group1_enrich",
                                                                              "num_FSNs_group2_enrich",
                                                                              "num_FSNs_at_nonBSNs")]))
       
      observed_counts <- as.numeric(summary_df[func_id, c("num_FSNs_group1_enrich",
                                                          "num_FSNs_group2_enrich",
                                                          "num_FSNs_at_nonBSNs")])
      
      if ((length(BSNs) > 0) && (summary_df[func_id, "num_FSNs"] >= multinomial_min_FSNs) && (prop_tested_BSNs != 1)) {
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
  
  results <- list(results=summary_df,
                  BSNs=pairwise_node_out$mean_direction[BSNs],
                  balances=calculated_balances)
  
  results[["multinomial_exp_prop"]] <- multinomial_exp_prop
  results[["FSNs_summary"]] <- func_summaries
  results[["tree"]] <- tree
  
  if (detailed_output) {
      # Restrict vector of mean directions to significant nodes only to avoid confusion.
      pairwise_node_out$mean_direction <- pairwise_node_out$mean_direction[BSNs]

      results[["balance_comparisons"]] <- pairwise_node_out
      results[["func_enrichments"]] <- all_node_enriched_funcs
      results[["input_param"]] <- input_param
  }
  
  return(results)

}


check_POMS_pipeline_args <- function(abun,
                                     func,
                                     tree,
                                     group1_samples,
                                     group2_samples,
                                     ncores,
                                     pseudocount,
                                     manual_BSNs,
                                     manual_balances,
                                     manual_BSN_dir,
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

  if (((is.null(manual_BSNs)) && (! is.null(manual_balances))) || ((! is.null(manual_BSNs)) && (is.null(manual_balances)))) {
    stop("Stopping - arguments manual_BSNs and manual_balances either both need to be given or neither should be specified.")
  }

  if (! is.null(manual_BSN_dir) && is.null(manual_BSNs)) {
    stop("Stopping - the manual_BSN_dir argument can only be set if the manual_BSNs argument is set as well.") 
  }
  
  if ((! is.null(manual_BSNs)) && (length(manual_BSNs) == 0)) { stop("Stopping - vector specified for manual_BSNs argument is empty.") }

  if ((! is.null(manual_BSNs)) && (! is.null(manual_balances))) {
    
    if (! is.null(manual_BSN_dir)) {

      if (length(group1_samples) > 0 || length(group2_samples) > 0) {
        stop("Stopping - group1_samples and group2_samples arguments should not be set when the manual_BSN_dir argument is specified.") 
      }
      
      if (class(manual_BSN_dir) != "character") {
        stop("Stopping - input to manual_BSN_dir argument needs to be a character vector.") 
      }
  
      if (length(manual_BSN_dir) != length(manual_BSNs)) {
        stop("Stopping - values for manual_BSN_dir and manual_BSNs should be vectors of the same length, but currently are not.")
      }
      
      if (length(which(manual_BSNs %in% names(manual_BSN_dir))) < length(manual_BSNs)) {
        stop("Stopping - all nodes in the manual_BSNs input must be present as names in manual_BSN_dir vector (when the manual_BSN_dir argument is used), but currently are not.")
      }

    }
    
    if (length(which(manual_BSNs %in% names(manual_balances))) != length(manual_BSNs)) {
      stop("Stopping - not all nodes in manual_BSNs vector are in in manual_balances") 
    }
    
    if (! "node.label" %in% names(tree)) {
      stop("Stopping - node labels must be present in tree if manual balances (i.e., manual_balances argument) are specificed.") 
    }
    
    if (length(which(! names(manual_balances) %in% tree$node.label)) > 0) {
      stop("Stopping - some balance labels (in manually input manual_balances argument) are missing from tree node labels.") 
    }
    
  }
  
  if (class(abun) != "data.frame") { stop("Stopping - argument abun needs to be of the class data.frame.") }
  if (class(func) != "data.frame") { stop("Stopping - argument func needs to be of the class data.frame.") }
  if (class(tree) != "phylo") { stop("Stopping - argument phylo needs to be of the class phylo.") }

  if (is.null(manual_BSN_dir)) {
    if (class(group1_samples) != "character") { stop("Stopping - argument group1_samples needs to be of the class character.") }
    if (class(group2_samples) != "character") { stop("Stopping - argument group2_samples needs to be of the class character.") }
    if (length(group1_samples) == 0) { stop("Stopping - argument group1_samples is of length 0.") }
    if (length(group2_samples) == 0) { stop("Stopping - argument group2_samples is of length 0.") }
    if (length(which(group1_samples %in% colnames(abun))) != length(group1_samples)) { stop("Stopping - not all group1_samples match columns of abun argument.") }
    if (length(which(group2_samples %in% colnames(abun))) != length(group2_samples)) { stop("Stopping - not all group2_samples match columns of abun argument.") }
    if (length(which(group1_samples %in% group2_samples)) > 0) { stop("Stopping - at least one sample overlaps between group1_samples and group2_samples, but these should be non-overlapping sets.") }
  }

  if ((class(ncores) != "integer") && (class(ncores) != "numeric")) { stop("Stopping - ncores argument needs to be of class numeric or integer.") }
  if (ncores <= 0) { stop("Stopping - ncores argument needs to be higher than 0.") }

  if ((class(pseudocount) != "integer") && (class(pseudocount) != "numeric")) { stop("Stopping - pseudocount argument needs to be of class numeric or integer.") }
  if (pseudocount < 0) { stop("Stopping - pseudocount argument cannot be lower than 0.") }

  if ((class(multinomial_min_FSNs) != "integer") && (class(multinomial_min_FSNs) != "numeric")) { stop("Stopping - multinomial_min_FSNs argument needs to be of class numeric or integer.") }
  if (multinomial_min_FSNs < 0) { stop("Stopping - multinomial_min_FSNs argument cannot be lower than 0.") }
  
  if (min_num_tips > length(tree$tip.label) / 2) { stop("Stopping - the min_num_tips argument cannot be higher than half of the total number of tips.") }

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

  return(list(abun=abun,
              func=func,
              tree=tree,
              group1_samples=group1_samples,
              group2_samples=group2_samples,
              ncores=ncores,
              pseudocount=pseudocount,
              manual_BSNs=manual_BSNs,
              manual_balances=manual_balances,
              manual_BSN_dir=manual_BSN_dir,
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


summarize_node_enrichment <- function(enriched_funcs, identified_BSNs, func_p_cutoff) {
  # For each function that is significant at least once get:
  #   - Names of BSNs where the function is present.
  #   - Names of positively and negatively enriched significant nodes (with respect to BSN direction)
  #   - Names of non-BSNs where the function was significant.

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

    func_summaries[[func_id]]$FSNs_group1_enrich <- as.character()
    func_summaries[[func_id]]$FSNs_group2_enrich <- as.character()
    func_summaries[[func_id]]$FSNs_at_nonBSNs <- as.character()
    
    for (node in names(enriched_funcs)) {
      if (func_id %in% rownames(enriched_funcs[[node]])) {

        if (node %in% identified_BSNs) {
          # If this node was significant then categorize it as either positive, negative, or nonenriched.
          if (as.numeric(enriched_funcs[[node]][func_id, "P_corr"]) < func_p_cutoff) {

            if (as.numeric(enriched_funcs[[node]][func_id, "OR"]) > 1) {
              func_summaries[[func_id]]$FSNs_group1_enrich <- c(func_summaries[[func_id]]$FSNs_group1_enrich, node)
            } else if (as.numeric(enriched_funcs[[node]][func_id, "OR"]) < 1) {
              func_summaries[[func_id]]$FSNs_group2_enrich <- c(func_summaries[[func_id]]$FSNs_group2_enrich, node)
            } else {
              print(enriched_funcs[[node]][func_id, ])
              stop("Significant function but OR for above info not different from 1?!")
            }
          }
    
        } else {
          # Since this node was NOT significant then categorize it as either enriched or nonenriched.
          if (as.numeric(enriched_funcs[[node]][func_id, "P_corr"]) < func_p_cutoff) {
            func_summaries[[func_id]]$FSNs_at_nonBSNs <- c(func_summaries[[func_id]]$FSNs_at_nonBSNs, node)
          }
        }
      }
    }
    
  }

  return(func_summaries)
}
