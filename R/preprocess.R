#' Filters out columns of dataframe based on number of proportion of non-zero cells
#'
#' Filters dataframe columns with either a low absolute count of non-zero values or a low proportion of rows with non-zero counts.
#' Note that this function is intended for positively-bounded data only (e.g., the function or taxon abundance tables), and will not work properly
#' if the table contains negative values. Included in package simply to make running workflow easier.
#'  
#' @param in_tab input dataframe
#'
#' @param min_nonzero_count minimum number of cells in column that must be non-zero for column to be retained.
#'
#' @param min_nonzero_prop minimum proportion of cells in column that must be non-zero for column to be retained.
#'
#' @param drop_missing_rows boolean flag to indicate whether rows with all zero values (after dropping columns based on specified cut-offs) should be removed.
#'
#' @param verbose boolean flag to indicate that the number of columns removed should be written to the console.
#'
#' @return dataframe with columns that did not meet the `min_nonzero_count` and/or `min_nonzero_prop` options removed (and potentially rows dropped too if drop_missing_rows=TRUE).
#'
#' @export
filter_rare_table_cols <- function(in_tab, min_nonzero_count, min_nonzero_prop, drop_missing_rows=TRUE, verbose=TRUE) {

  if (class(in_tab) != "data.frame") {
    stop("Error - argument \"in_tab\" not class data.frame") 
  }
  
  nonzero_counts <- colSums(in_tab > 0)

  col2remove <- which((nonzero_counts < min_nonzero_count) | ((nonzero_counts / nrow(in_tab)) < min_nonzero_prop))
  
  if(length(col2remove) > 0) {
    
    if (verbose) { message("Filtering out ", length(col2remove), " rare functions from input function table.") }

    in_tab <- in_tab[, -col2remove, drop = FALSE]
    
  } else {
    
    if (verbose) { message("No columns removed.") }
    
  }
  
  if (drop_missing_rows) {
   
    missing_rows <- which(rowSums(in_tab > 0) == 0)
     
    if (length(missing_rows) > 0) {

      if (verbose) { message("Filtering out ", length(missing_rows), " rows that contain no non-zero values.") }
  
      in_tab <- in_tab[-missing_rows, , drop = FALSE]
    
    }
    
  }
  
  return(in_tab)
  
}



#' Subset dataframe by column names and then post-filter
#'
#' Subset table by set of column names. After doing this, it will remove any rows and columns that are all 0's.
#'  
#' @param in_tab input dataframe
#' 
#' @param col2keep column names to retain in output (as long as they have at least one non-zero value).
#' 
#' @param verbose flag to indicate that the final number of rows and columns (as well as the number removed) should be reported.
#'
#' @return dataframe with subset of specified columns (if they have at least one non-zero value), also with rows that only contain 0's removed.
#' 
#' @export
subset_by_col_and_filt <- function(in_tab, col2keep, verbose = TRUE) {
  
  if (class(in_tab) != "data.frame") {
    stop("Error - argument \"in_tab\" not class data.frame") 
  }
  
  if (class(col2keep) != "character") {
    stop("Error - argument \"col2keep\" not class character") 
  }
  
  if (length(col2keep) == 0) {
    stop("Error - at least one column name to keep must be specified.")
  }
  
  if (length(which(! col2keep %in% colnames(in_tab))) > 0) {
    stop("Error - not all specified columns in \"col2keep\" are present in table. Please adjust.")
  }

  in_tab <- in_tab[, col2keep, drop = FALSE]
  
  if (verbose) {
    message("Subset dataframe to ", length(col2keep), " specified columns.") 
  }
  
  missing_rows <- which(rowSums(in_tab) == 0)
  missing_cols <- which(colSums(in_tab) == 0)
  
  if (length(missing_rows) > 0) {
    in_tab <- in_tab[-missing_rows, ]
    
    if (verbose) {
      message("Removed ", length(missing_rows), " rows with no non-zero values.")
    }
  } else if (verbose) {
    message("All rows contain at least one non-zero value.")
  }
  
  if (length(missing_cols) > 0) {
    in_tab <- in_tab[, -missing_cols]
    if (verbose) {
      message("Removed ", length(missing_cols), " columns with no non-zero values (of the set remaining after subsetting by column name).")
    }
  } else if (verbose) {
      message("All remaining columns after subsetting by column name contain at least one non-zero value.")
  }
  
  if (verbose) {
    message("Returning dataframe with ", nrow(in_tab), " rows and ", ncol(in_tab), " columns.")
  }
  
  return(in_tab)

}

#' Get node indices of FSN and BSN categories across tree for a given function
#'
#' Parse POMS_pipeline output to look at FSNs for a specific function (e.g., a specific gene family). Will also parse BSN information (which is not dependent on a particular function).
#' This is convienient to do before plotting the distribution of FSNs and BSNs across the tree with the ggtree R package for instance.
#' When a taxa label table is specified, labels of tested nodes in the tree (found in the POMS_pipeline output object) will be renamed to be the representative taxa on each side.
#'  
#' @param POMS_output output object from POMS_pipeline function.
#'
#' @param func_id label of function for which should FSNs should be parsed. Must be present in POMS_output$FSNs_summary.
#'
#' @param taxa_table optional dataframe containing taxa labels for each tip of tree. Must be in same format as expected for node_taxa function.
#'
#' @param taxa_threshold float > 0.5 and <= 1.0 specifying the proportion of tips that must share a taxon label for it to be considered representative. Only relevant if taxa_table specified.
#'
#' @param full_taxon_label boolean flag for whether taxon labels should be combined, so that all higher taxonomic labels are included. Specifically, when TRUE, all higher labels are concatenated and delimited by "; ".
#' E.g., rather than just the genus "Odoribacter" the label would be "Bacteria; Bacteroidetes; Bacteroidia; Bacteroidales; Porphyromonadaceae; Odoribacter". Only relevant if taxa_table specified.
#'
#' @return List containing final tree as well as indices of nodes corresponding to different FSN and BSN categories.
#' If taxa_table was specified, then node labels in tree will correspond to representative taxa on each side of the
#' nodes that were tested (i.e., those that were non-negligible).
#'
#' @export
prep_func_node_info <- function(POMS_output, func_id, taxa_table = NULL, taxa_threshold = 0.75, full_taxon_label = FALSE) {
  
  num_tips <- length(POMS_output$tree$tip.label)
  
  FSNs_group1_enrich_i <- which(POMS_output$tree$node.label %in% POMS_output$FSNs_summary[[func_id]]$FSNs_group1_enrich) + num_tips
  FSNs_group2_enrich_i <- which(POMS_output$tree$node.label %in% POMS_output$FSNs_summary[[func_id]]$FSNs_group2_enrich) + num_tips
  FSNs_at_nonBSNs_i <- which(POMS_output$tree$node.label %in% POMS_output$FSNs_summary[[func_id]]$FSNs_at_nonBSNs) + num_tips
  
  all_FSNs_i <- c(FSNs_group1_enrich_i, FSNs_group2_enrich_i, FSNs_at_nonBSNs_i)
  
  all_BSNs_i <- which(POMS_output$tree$node.label %in% names(POMS_output$BSNs)) + num_tips
  BSNs_at_nonFSNs_i <- all_BSNs_i[which(! all_BSNs_i %in% all_FSNs_i)]
  
  tested_nodes_i <- which(POMS_output$tree$node.label %in% names(POMS_output$balance_info$balances)) + num_tips
  
  # Specify taxa names at all tested nodes if taxa table provided.
  if (! is.null(taxa_table)) {
    tested_nodes <- names(POMS_output$balance_info$balances)
    tested_node_taxa <- list()
    for (node in tested_nodes) {
      tested_node_taxa[[node]] <- node_taxa(in_tree = POMS_output$tree,
                                            taxon_labels = taxa_table,
                                            node_label = node,
                                            threshold = taxa_threshold,
                                            combine_labels = full_taxon_label)
    }
    
    prepped_tree <- POMS_output$tree
    
    nodes2ignore <- which(! prepped_tree$node.label %in% tested_nodes)
    prepped_tree$node.label[nodes2ignore] <- ""
    
    for (node in tested_nodes) {
      prepped_tree$node.label[which(prepped_tree$node.label == node)] <- paste(tested_node_taxa[[node]], collapse = " / ")
    }
    
  } else {
    prepped_tree <- POMS_output$tree
  }
  
  return(list(prepped_tree = prepped_tree,
              tested_nodes_i = tested_nodes_i,
              all_FSNs_i = all_FSNs_i,
              all_BSNs_i = all_BSNs_i,
              FSNs_group1_enrich_i = FSNs_group1_enrich_i,
              FSNs_group2_enrich_i = FSNs_group2_enrich_i,
              FSNs_at_nonBSNs_i = FSNs_at_nonBSNs_i,
              BSNs_at_nonFSNs_i = BSNs_at_nonFSNs_i))

}


prep_tree <- function(phy, tips2keep) {

  # Function to subset tree to specified tips and run sanity checks.
  # Will midpoint root tree if necessary.
  # Will also add node labels to tree before returning (unless they are already present).
  
  tips2remove <- phy$tip.label[which(! phy$tip.label %in% tips2keep)]
  
  phy <- ape::drop.tip(phy = phy, tip = tips2remove)
  
  if(! ape::is.binary(phy)) {
    stop("Stopping - input tree is non-binary.")
  }
  
  if(! ape::is.rooted(phy)) {
    phy <- phangorn::midpoint(phy)
  }
  
  if (! "node.label" %in% names(phy)) {
    phy$node.label <- NULL
    phy <- ape::makeNodeLabel(phy, method="number", prefix='n')
  }
  
  return(phy)

}
