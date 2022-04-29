#' Filters out columns of dataframe based on number of proportion of non-zero cells
#'
#' Filters dataframe columns with either a low absolute count of non-zero values or a low proportion of rows with non-zero counts.
#' Note that this function is intended for positively-bounded data only (e.g., the function or taxon abundance tables), and will not work properly
#' if the table contains negative values. Included in package simply to make running workflow easier.
#'  
#' @param in_tab Input dataframe
#'
#' @param min_nonzero_count Min number of cells in column that must be non-zero for column to be retained.
#'
#' @param min_nonzero_prop Min proportion of cells in column that must be non-zero for column to be retained.
#'
#' @param verbose Boolean flag to indicate whether rows with all zero values (after dropping columns based on specified cut-offs) should be removed.
#'
#' @param verbose Boolean flag to indicate that the number of columns removed should be written to the console.
#'
#' @return Dataframe with columns that did not meet the `min_nonzero_count` and/or `min_nonzero_prop` options removed.
#'
#' @export
filter_rare_table_cols <- function(in_tab, min_nonzero_count, min_nonzero_prop, drop_missing_rows=TRUE, verbose=TRUE) {

  if (class(in_tab) != "data.frame") {
    stop("Error - argument \"in_tab\" not class data.frame") 
  }
  
  nonzero_counts <- colSums(in_tab > 0)

  col2remove <- which((nonzero_counts < min_nonzero_count) | ((nonzero_counts / nrow(in_tab)) < min_nonzero_prop))
  
  if(length(col2remove) > 0) {
    
    if (verbose) { message(paste("Filtering out", as.character(length(col2remove)), "rare functions from input function table.", sep = " ")) }

    in_tab <- in_tab[, -col2remove, drop = FALSE]
    
  } else {
    
    if (verbose) { message("No columns removed.") }
    
  }
  
  if (drop_missing_rows) {
   
    missing_rows <- which(rowSums(in_tab > 0) == 0)
     
    if (length(missing_rows) > 0) {

      if (verbose) { message(paste("Filtering out", as.character(length(missing_rows)), "rows that contain no non-zero values.", sep = " ")) }
  
      in_tab <- in_tab[-missing_rows, , drop = FALSE]
    
    }
    
  }
  
  return(in_tab)
  
}

#' Subset dataframe by column names and then post-filter
#'
#' Subset table by set of column names. After doing this, it will remove any rows and columns that are all 0's.
#'  
#' @param in_tab Input dataframe
#' 
#' @param col2keep Column names to retain in output (as long as they have at least one non-zero value).
#' 
#' @param verbose Flag to indicate that the final number of rows and columns (as well as the number removed) should be reported.
#'
#' @return Dataframe with subset of specified columns (if they have at least one non-zero value), also with all rows with no non-zero values removed.
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
    message(paste("Subset dataframe to", as.character(length(col2keep)), "specified columns.", sep = " ")) 
  }
  
  missing_rows <- which(rowSums(in_tab) == 0)
  missing_cols <- which(colSums(in_tab) == 0)
  
  if (length(missing_rows) > 0) {
    in_tab <- in_tab[-missing_rows, ]
    
    if (verbose) {
      message(paste("Removed", as.character(length(missing_rows)), "rows with no non-zero values.", sep = " "))
    }
  } else if (verbose) {
    message("All rows contain at least one non-zero value.")
  }
  
  if (length(missing_cols) > 0) {
    in_tab <- in_tab[, -missing_cols]
    if (verbose) {
      message(paste("Removed", as.character(length(missing_cols)), "columns with no non-zero values (of the set remaining after subsetting by column name).", sep = " "))
    }
  } else if (verbose) {
      message("All remaining columns after subsetting by column name contain at least one non-zero value.")
  }
  
  if (verbose) {
    message("Returning dataframe with ", as.character(nrow(in_tab)), "rows and", as.character(ncol(in_tab)), "columns.", sep = " ")
  }
  
  return(in_tab)

}


prep_tree <- function(phy, tips2keep) {

  # Function to subset tree to specified tips only.
  # and run sanity checks.
  # Will midpoint root tree if necessary.
  # Add node labels to tree before returning (unless they are already present).
  
  tips2remove <- phy$tip.label[which(! phy$tip.label %in% tips2keep)]
  
  phy <- ape::drop.tip(phy = phy, tip = tips2remove)
  
  if(! ape::is.binary(phy)) {
    stop("Tree is non-binary.")
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

