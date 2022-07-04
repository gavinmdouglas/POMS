#' Wrapper for running phylogenetic regression with phylolm\cr
#'
#' Runs basic case of single x and y variables (dummy or continuous). Note that the ordering of the input vectors
#' and the tree tip labels needs to be checked by the user beforehand:
#' this script does not require that the y and x variables are named, and so no name check is performed.
#'
#' @param y variable to use for y component of model.
#' 
#' @param x variable to use for x component of model.
#' 
#' @param in_tree phylo object. Tip label order is assumed to match the y and x variables.
#' 
#' @param model_type length-one character vector specifying which phylogenetic model to use (must be a possible setting of the model argument to the phylolm function).
#'
#' @return Numeric vector of length three, providing the estimated coefficients for the intercept and slope, along with the p-value.
#' 
#' @export
phylolm_summary <- function(y, x, in_tree, model_type = "BM") {
  
  if (length(y) != length(x)) {
    stop("Stopping - input vectors are not of equal length.")  
  } else if (length(y) != length(in_tree$tip.label)) {
    stop("Stopping - input vectors need to be the same length as tip labels.") 
  } 
  
  in_data <- data.frame(x = x, y = y)
  rownames(in_data) <- in_tree$tip.label
  
  phylolm_out <- summary(phylolm::phylolm(y ~ x, data = in_data, phy = in_tree, model = model_type))
  
  coef_out <- c(phylolm_out$coefficients[, "Estimate"], phylolm_out$coefficients[2, "p.value"])
  names(coef_out) <- c("intercept", "slope", "p")

  return(coef_out)
  
}


#' Phylogenetic regression of input vector against function presence/absence.
#'
#' Runs phylogenetic regression with phylolm on each function (or trait) in the specified function table. 
#'
#' @param y variable to use for y component of model. Typically would be either a binary vector indicating which taxa are significantly different, or the normalized specicity or normalized prevalence values.
#' Must be a named numeric vector with names matching the rows of the func dataframe.
#' These names also must match the tree tip labels, but they can be a subset and any missing tips will be dropped.
#' 
#' @param func dataframe of the number of copies of each function that are encoded by each input taxon.
#' This pipeline only considers the presence/absence of functions across taxa.
#' Taxa (with row names intersecting with the "abun" table) should be the rows and the functions should be the columns.
#' 
#' @param in_tree phylo object. Tip labels must include the row names of the func dataframe and the names of the y input vector.
#' 
#' @param ncores integer specifying how many cores to use for parallelized sections of pipeline.
#'
#' @param model_type length-one character vector specifying which phylogenetic model to use (must be a possible setting of the model argument to the phylolm function).
#'
#' @return Dataframe summarizing the phylolm coefficients and model p-values for each y ~ function comparison.
#' Will include the intercept, slope, and p-value for each case. Row names will be function ids.
#' 
#' @export
genome_content_phylo_regress <- function(y, func, in_tree, ncores = 1, model_type = "BM") {
  
  if (length(y) != nrow(func)) {
    stop("Stopping - y length and number of rows in func are mismatched.")  
  } else if (length(which(! names(y) %in% rownames(func))) > 0) {
    stop("Stopping - names of elements in y and row names of func dataframe are mismatched.") 
  } else if (length(which(! names(y) %in% in_tree$tip.label)) > 0) {
    stop("Stopping - not all taxa in func and y are present as tips in tree.")
  }

  if ((class(ncores) != "integer") && (class(ncores) != "numeric")) {
    stop("Stopping - ncores argument needs to be of class numeric or integer.")
  } else if (length(ncores) != 1) {
    stop("Stopping - ncores argument must be of length 1.")
  } else if (ncores <= 0) {
    stop("Stopping - ncores argument must be 1 or higher.")
  }
  
  if (! model_type %in% c("BM", "OUrandomRoot", "OUfixedRoot", "lambda", "kappa", "delta", "EB", "trend")) {
     stop("Stopping - specified model_type not a model option for phylolm.")
  }
  
  raw_out <- parallel::mclapply(X = 1:ncol(func), FUN = function(func_i) { phylolm_summary(y = y,
                                                                                           x = func[, func_i, drop = TRUE],
                                                                                           in_tree = in_tree,
                                                                                           model_type = model_type)
                                                                          }, mc.cores = ncores)

  output <- as.data.frame(data.table::transpose(raw_out),
                          col.names = c("intersect", "slope", "p"),
                          row.names = colnames(func))

  return(output)
  
}
