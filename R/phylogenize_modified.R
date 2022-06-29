#' Calculate logit of a value or vector of values.
#' This function was taken unchanged from the phylogenize R codebase
#' (except it is no longer exported).
#'
logit <- function(x) (log(x / (1 - x)))

#' Calculate inverse-logit of a value or vector of values.
#' This function was taken unchanged from the phylogenize R codebase.
#' (except it is no longer exported).
logistic <- function(x) exp(x) / (1 + exp(x))

#' Apply a function to a vector of names, with the returned list having those
#' names.
#' This function was taken unchanged from the phylogenize R codebase.
# (except the description lines were removed).
lapply.across.names <- function(X, FUN, ...) {
  r <- lapply(X, FUN, ...)
  names(r) <- X
  r
}

#' Compute additive smoothed prevalence of features (e.g, taxa), restricted to samples of a particular metadata category.\cr
#'
#' This code replicates the prevalence score introduced in phylogenize. The code here is modified from the phylogenize code base 
#' (https://bitbucket.org/pbradz/phylogenize/src/master/package/phylogenize/R/; commit 6f1bdba9c5a9ff04e90a8ad77bcee8ec9281730d).
#' 
#' This algorithm is descibed in detail in Bradley et al. 2018. Phylogeny-corrected identification of microbial gene families relevant to human gut colonization. PLOS Computational Biology.
# 10.1371/journal.pcbi.1006242. Note that this version of the algorithm does not allow for multiple datasets.
#'
#' @param abun_table Abundance table to use for computing prevalence. Features must be rows and samples columns. All values greater than 0 will be interpreted as present.
#' 
#' @param meta_table Dataframe object containing metadata for all samples. Must include at least one column corresponding to the sample ids and one column containing the metadata of interest that will be focused on when computing prevalence.
#' 
#' @param focal_var_level length-one character vector specifying the variable value to restrict inferences of prevalence to. In other words, prevalence will be computed based on the sample set that contain this value of the variable of interest in the metadata table. 
#' 
#' @param var_colname length-one character vector specifying the name of column in the metadata table that contains the metadata of interest (e.g., where focal_var_level can be found).
#' 
#' @param sample_colname length-one character vector specifying the name of column in the metadata table that contains the sample ids.
#' 
#' @param silence_citation length-one Boolean vector specifying whether to silence message notifying user about phylogenize package and paper.
#'
#' @return Numeric vector with the normalized prevalence score for each input feature (i.e., for each row of abun_table).
#' 
#' @export
prevalence_norm_logit <- function(abun_table,
                                  meta_table,
                                  focal_var_level,
                                  var_colname,
                                  sample_colname,
                                  silence_citation = FALSE) {

  if (! silence_citation) {
    message("The prevalence calculations performed here are closely based on code distributed with the phylogenize package, which is under an MIT license.")
    message("If you use these normalized prevalence values, please cite the original authors: 10.1371/journal.pcbi.1006242.")
  }

  if (! var_colname %in% colnames(meta_table)) {
    stop("Stopping - specified variable column name ", var_colname, " not found in metadata table.")
  }
  
  if (! sample_colname %in% colnames(meta_table)) {
    stop("Stopping - specified sample column name ", sample_colname, " not found in metadata table.")
  }
  
  if (! focal_var_level %in% meta_table[, var_colname]) {
    stop("Stopping - variable value ", focal_var_level, " not found in specified metadata column.")
  }

  focal_var_rows <- which(meta_table[, var_colname] == focal_var_level)

  intersecting_samples <- as.character(intersect(colnames(abun_table), meta_table[focal_var_rows, sample_colname]))
  
  if (length(intersecting_samples) == 0) {
    stop("Stopping - no sample ids intersect between abundance table column names and specified sample column of metadata table.")
  }
  
  num_samples <- length(intersecting_samples)

  taxa_raw_prevalence <- rowSums(1 * (abun_table[, intersecting_samples, drop = FALSE] > 0))
  
  taxa_norm_prevalence <- (1 + taxa_raw_prevalence) / (2 + num_samples)

  # Compute logit of these normalize prevalences and return.
  return(logit(taxa_norm_prevalence))

}



sim_presence_absence <- function(effect.size = 2,
                                baseline.distro = c(shape1 = 0.66,
                                                    shape2 = 2.62),
                                samples = c(H = 38, D = 13),
                                focal_level = "D",
                                taxa = 2000,
                                tpr = 0.25,
                                sign.pos.prob = 0.5) {
  tp.taxa <- round(tpr * taxa)
  neg.taxa <- taxa - (tp.taxa)
  n.classes <- length(samples)
  pT <- sapply(1:taxa, function(.) {
    rbeta(1, baseline.distro[1], baseline.distro[2])
  })
  fx <- c(((2 * rbinom(n = tp.taxa, size = 1, sign.pos.prob)) - 1),
          rep(0, neg.taxa))
  pTbs <- sapply(1:taxa, function(i) {
    pTb <- logistic(logit(pT) + (fx * (effect.size)))
  })
  if (is.null(names(samples))) names(samples) <- 1:length(samples)
  sim.mtx <- t(sapply(1:taxa, function(i) {
    Reduce(c, lapply.across.names(names(samples), function(smp) {
      if (smp == focal_level) { p <- pTbs[i] } else { p <- pT[i] }
      rbinom(samples[smp], size = 1, p)
    }))
  }))
  ids <- Reduce(c, lapply(1:length(samples), function(i) rep(i, samples[i])))
  return(list(mtx = sim.mtx,
              pT = pT,
              pTbs = pTbs,
              fx = fx,
              ids = ids,
              input.params = list(effect.size = effect.size,
                                  baseline = baseline.distro,
                                  samples = samples,
                                  taxa = taxa,
                                  tpr = tpr,
                                  sign.pos.prob = sign.pos.prob)))
}


#' Score a simulated regularization by how well it recapitulates the ground truth.
#' This function was taken unchanged from the phylogenize R codebase.
# (except the description lines were removed).
score.regularization <- function(mtx,
                                 ids,
                                 real.fx,
                                 which.env = 2,
                                 prior = 0.002,
                                 b = 0.1,
                                 add.pc = FALSE,
                                 tol = prior * 0.005,
                                 ...) {
  regularized <- apply(mtx, 1, function(x) {
    regularize.pET(x,
                   ids,
                   which.env = which.env,
                   prior = prior,
                   b = b,
                   add.pc = add.pc)
  })
  posteriors <- regularized[2, , drop=TRUE]
  signif <- 1 * (abs(posteriors - prior) > tol)
  predicted.signs <- signif * sign(posteriors - prior)
  fpr <- (signif[real.fx == 0] %>% mean)
  pwr.hi <- mean(predicted.signs[real.fx == 1] == 1)
  pwr.lo <- mean(predicted.signs[real.fx == -1] == -1)
  return(c(fpr = fpr, pwr.hi = pwr.hi, pwr.lo = pwr.lo))
}

#' Get a distribution of environment prevalences from a matrix.
#' This function was taken unchanged from the phylogenize R codebase.
#' 
#' This is used to optimize the value of the free parameter $b$ in
#' \code{regularize.pET}.
#'
#' @param mtx A matrix of presence/absence values.
#' @param ids A named factor assigning samples (matrix columns) to environments.
#' @param fallback A two-element numeric vector, giving the beta parameters to
#'     use if fitting fails.
#' @return A best-fit of prevalences to a beta distribution.
#' @keywords internal
fit.beta.list <-  function(mtx, ids, fallback = c(NA, NA)) {
  lapply(unique(ids), function(i) {
    tryCatch(
      MASS::fitdistr(densfun = "beta",
                     start = list(shape1 = 1, shape2 = 1),
                     apply(mtx[, which(ids == i), drop=FALSE], 1, function(x) {
                       mean(c(x, 0, 1) > 0)
                     }))$estimate,
      error = function(e) fallback)
  })
}

#' Based on a real presence/absence matrix, optimize the value of the
#' regularization parameter $b$.
optimize_b_wrapper <- function(real_abun_table,
                               real_sample_values,
                               focal_var_level,
                               effect.size = 2,
                               prior,
                               tol,
                               bounds = c(0.01, 5),
                               add.pc = TRUE,
                               a = 0.05,
                               pos.prop = 0.5) {
  
  shape.n <- which(sort(unique(real_sample_values)) == focal_var_level)
  
  # Fall back to the overall beta if not enough information to fit a
  # distribution (or too wacky)
  overall.shape <- suppressWarnings(MASS::fitdistr(densfun = "beta",
                                                   start = list(shape1 = 1, shape2 = 1),
                                                   apply(real_abun_table, 1, function(x) {
                                                     mean(c(x, 0, 1) > 0)
                                                   }))$estimate)

  shapes <- suppressWarnings(fit.beta.list(real_abun_table,
                                           real_sample_values,
                                           fallback = overall.shape))

  level_tallies <- table(real_sample_values)
  var_order <- names(level_tallies)
  level_tallies <- as.numeric(level_tallies)
  names(level_tallies) <- var_order
  
  level_tallies_focal_i <- which(names(level_tallies) == focal_var_level)
  
  get.optim <- function(b) {

    sim <- sim_presence_absence(effect.size = effect.size,
                                baseline.distro = shapes[[shape.n]],
                                samples = level_tallies,
                                taxa = nrow(real_abun_table),
                                tpr = 0.25,
                                sign.pos.prob = pos.prop,
                                focal_level = focal_var_level)

    s <- score.regularization(sim$mtx,
                              sim$ids,
                              sim$fx,
                              level_tallies_focal_i,
                              prior,
                              b,
                              add.pc,
                              tol)

    # Summarize the statistics from score.regularization into a single metric.
    # If proportion of false positives is a value, return 1 - FPR; otherwise,
    # return 1 + the average (geometric mean) power on positive and negative effect sizes.
    if (s["fpr"] <= a) {
      return(1 + exp(sum(log(s[c("pwr.hi", "pwr.lo")])) / length(s[c("pwr.hi", "pwr.lo")])))
    } else {
      return((1 - s["fpr"]))
    }
  }
  
  suppressWarnings(optimize(get.optim, bounds, maximum = TRUE))

}


#' Obtain a regularized estimate of specificity.
#' This function was taken unchanged from the phylogenize R codebase.
# (except the description lines were removed).
regularize.pET <- function(vec,
                           env.ids,
                           which.env = 1,
                           prior = 0.05,
                           b = 1,
                           add.pc = FALSE,
                           min.limit = -10,
                           max.limit = 10) {

  # n and k in binomial are #(E) and #(E,T), respectively (used on P(T|E))
  n <- sum(env.ids == which.env) # #(E)
  k <- sum(vec[env.ids == which.env] > 0) # #(E & T)
  if (add.pc) {
    n <- n + 1
    k <- k + 1
  }
  if (!add.pc) {
    pT.E <- mean(c(vec[which(env.ids == which.env)] > 0))
    pT.nE <- mean(c(vec[which(env.ids != which.env)]) > 0)
  } else {
    pT.E <- mean(c(0, 1, vec[which(env.ids == which.env)] > 0))
    pT.nE <- mean(c(0, 1, vec[which(env.ids != which.env)]) > 0)
  }
  pT <- pT.E * prior + pT.nE * (1 - prior)
  map <- function(p) {
    x <- (p * pT) / prior # P(T|C)
    dbinom(k, n, x) * ((1/2*b)) * exp(-(abs(logit(p)-logit(prior))/b))
  }
  map.logit <- function(logit.p) {
    p <- logistic(logit.p)
    x <- (p * pT) / prior # P(T|C)
    # return log probability also, better numerical stability
    dbinom(k, n, x, log = TRUE) +
      log(1 / (2*b)) -
      (abs(logit(p)-logit(prior))/b)
  }
  max.p <- prior / pT
  max.p <- min(max.p, logistic(max.limit))
  results <- optimize(map.logit, c(min.limit, logit(max.p)), maximum = TRUE)
  initial.x <- pT.E
  initial.p <- (pT.E * prior) / pT
  final.x <- logistic(results$maximum) * pT / prior
  final.p <- logistic(results$maximum)
  return(c(x = final.x,
           p = final.p,
           x.init = initial.x,
           p.init = initial.p,
           pT = pT))
}


#' Test whether a value is between two other values (non-inclusive).
#' This function was taken unchanged from the phylogenize R codebase.
# (except the description lines were removed).
`%btwn%` <- function(x, y) { (x > min(y)) & (x < max(y)) }


#' Compute shrunken specificity score of a feature, which represents how the presence of a feature is associated with a given sample grouping.
#'
#' This code replicates the environmental specificity score introduced in phylogenize. The code here is modified from the phylogenize code base 
#' (https://bitbucket.org/pbradz/phylogenize/src/master/package/phylogenize/R/; commit 6f1bdba9c5a9ff04e90a8ad77bcee8ec9281730d).
#' 
#' This algorithm is descibed in detail in Bradley et al. 2018. Phylogeny-corrected identification of microbial gene families relevant to human gut colonization. PLOS Computational Biology.
# 10.1371/journal.pcbi.1006242.
#'
#' Note thee can be some random fluctuations between re-runs of this function. The differences are usually minor, but users are strongly suggested to set a random seed before use to ensure their workflow is reproducible.
#'
#' @param abun_table Abundance table to use for computing specificity Features must be rows and samples columns. All values greater than 0 will be interpreted as present.
#' 
#' @param meta_table Dataframe object containing metadata for all samples. Must include at least one column corresponding to the sample ids and one column containing the metadata of interest that will be focused on.
#' 
#' @param focal_var_level length-one character vector specifying the variable value to restrict inferences of prevalence to. In other words, prevalence will be computed based on the sample set that contain this value of the variable of interest in the metadata table. 
#' 
#' @param var_colname length-one character vector specifying the name of column in the metadata table that contains the metadata of interest (e.g., where focal_var_level can be found).
#' 
#' @param sample_colname length-one character vector specifying the name of column in the metadata table that contains the sample ids.
#' 
#' @param silence_citation length-one Boolean vector specifying whether to silence message notifying user about phylogenize package and paper.
#'
#' @return Numeric vector with the specificity score for each input feature (i.e., for each row of abun_table).
#' 
#' @export
specificity_scores <- function(abun_table,
                               meta_table,
                               focal_var_level,
                               var_colname,
                               sample_colname,
                               silence_citation = FALSE) {
  if (! silence_citation) {
    message("The specificity score calculations performed here are closely based on code distributed with the phylogenize package, which is under an MIT license.")
    message("If you use these values, please cite the original authors: 10.1371/journal.pcbi.1006242.")
  }

  if (! var_colname %in% colnames(meta_table)) {
    stop("Stopping - specified variable column name ", var_colname, " not found in metadata table.")
  }
  
  if (! sample_colname %in% colnames(meta_table)) {
    stop("Stopping - specified sample column name ", sample_colname, " not found in metadata table.")
  }
  
  if (! focal_var_level %in% meta_table[, var_colname]) {
    stop("Stopping - variable value ", focal_var_level, " not found in specified metadata column.")
  }
  
  meta_table[, sample_colname] <- as.character(meta_table[, sample_colname])
  meta_table[, var_colname] <- as.character(meta_table[, var_colname])
  
  focal_var_rows <- which(meta_table[, var_colname] == focal_var_level)

  intersecting_samples_i <- which(meta_table[, sample_colname] %in% colnames(abun_table))
  
  if (length(intersecting_samples_i) == 0) {
    stop("Stopping - no sample ids intersect between abundance table column names and specified sample column of metadata table.")
  }
  
  meta_table_present <- meta_table[intersecting_samples_i, ]

  # Compute uninformative prior.
  prior <- 1 / length(unique(meta_table_present[, var_colname]))

  sample_values <- sapply(colnames(abun_table),
                          function (sampleid) { meta_table_present[which(meta_table_present[, sample_colname] == sampleid), var_colname] })
  
  if (length(unique(sample_values)) < 2) {
    stop("Stopping - only one unique variable value found after restricting to samples in the abundance table.")
  }
  
  tolerance <- prior * 0.01

  b.optim <- optimize_b_wrapper(real_abun_table = abun_table,
                                real_sample_values = sample_values,
                                focal_var_level = focal_var_level,
                                effect.size = 2,
                                prior = prior,
                                tol = tolerance,
                                bounds = c(0.01, 5),
                                a = 0.05,
                                pos.prop = 0.5,
                                add.pc = TRUE)$maximum

  regularized <- apply(abun_table, 1, function(x) {
    regularize.pET(x,
                   sample_values,
                   which.env = focal_var_level,
                   b = b.optim,
                   prior = prior,
                   add.pc = TRUE)
  })

  logist.pheno <- regularized[2, , drop=TRUE]
  
  # "hard-shrink" anything shrunk almost
  # to the prior to prevent these tiny
  # differences from affecting the result
  # in the absence of a strong change
  logist.pheno[which(logist.pheno %btwn%
                       c(prior - tolerance,
                         prior + tolerance))] <- prior
  
  phenoP <- logit(prior)
  
  return(list(b.optim = b.optim,
              ess = logit(logist.pheno),
              regularized = regularized,
              prior = prior,
              phenoP = phenoP))
}
