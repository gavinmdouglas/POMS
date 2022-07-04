
library(POMS)

ex_taxa_abun <- read.table("../../example_files/ex_taxa_abun.tsv.gz", header = TRUE, sep = "\t", row.names = 1)
ex_func <- read.table("../../example_files/ex_func.tsv.gz", header = TRUE, sep = "\t", row.names = 1)

ex_group1 <- read.table("../../example_files/ex_group1.txt.gz", stringsAsFactors = FALSE)$V1
ex_group2 <- read.table("../../example_files/ex_group2.txt.gz", stringsAsFactors = FALSE)$V1

ex_tree <- ape::read.tree("../../example_files/ex_tree.newick")

ex_taxa_abun <- ex_taxa_abun[ex_tree$tip.label, ]
ex_func <- ex_func[ex_tree$tip.label, ]

ex_meta <- data.frame(Sample = colnames(ex_taxa_abun), Grouping = NA)
ex_meta[which(ex_meta$Sample %in% ex_group1), "Grouping"] <- "Group1"
ex_meta[which(ex_meta$Sample %in% ex_group2), "Grouping"] <- "Group2"

test_that("phylolm_summary returns expected values.", {
  
  # Note that just using two columns of function column as a quick test.
  observed <- phylolm_summary(y = ex_func[, 1, drop = TRUE],
                              x = ex_func[, 2, drop = TRUE],
                              in_tree = ex_tree)

  # Can compare with:
  # summary(phylolm::phylolm(formula = K07106 ~ K02036, data = ex_func, phy = ex_tree))$coefficients
  expected <- c(0.463879166, 0.004026632, 0.9746031)
  names(expected) <- c("intercept", "slope", "p")
  
  expect_equal(observed, expected, tolerance = 1e-7)

})


test_that("prevalence_norm_logit returns expected values.", {
  
  observed <- prevalence_norm_logit(abun_table = ex_taxa_abun[c(1, 4, 10, 14, 30, 41), ],
                                    meta_table = ex_meta,
                                    focal_var_level = "Group1",
                                    var_colname = "Grouping",
                                    sample_colname = "Sample",
                                    silence_citation = TRUE)
  
  expected <- c(5.170484, 5.866468, 5.866468, 5.866468, -4.762174, -4.471639 )
  names(expected) <- c("SRR5558051_bin.9", "GCF_000243215", "ERR1190644_bin.42",
                       "ERR866592_bin.16", "ERR414271_bin.13", "ERR1913101_bin.17")
  expect_equal(observed, expected, tolerance = 1e-7)
  
})


test_that("specificity_scores returns expected values.", {
  
  # Set seed to make reproducible.
  set.seed(1231)
  
  observed <- specificity_scores(abun_table = ex_taxa_abun[c(2, 20, 30, 40, 41), ],
                                    meta_table = ex_meta,
                                    focal_var_level = "Group1",
                                    var_colname = "Grouping",
                                    sample_colname = "Sample",
                                    silence_citation = TRUE)
  
  # Can compare with output from calc.ess function of phylogenize.
  # Note that this function is not 100% reproducible.
  expected <- c(4.4242361, 0.0000000, -0.3790120, 0.1635954, 0.0000000)
  names(expected) <- c("ERR321543_bin.17", "ERR1305897_bin.4",
                       "ERR414271_bin.13", "ERR688539_bin.27",
                       "ERR1913101_bin.17")
  expect_equal(observed$ess, expected, tolerance = 1e-7)
  
})


test_that("genome_content_phylo_regress returns expected values.", {
  
  # Use mean abundance as test values.
  mean_abun <- rowMeans(ex_taxa_abun)
  
  observed <- genome_content_phylo_regress(y = mean_abun,
                                           func = ex_func,
                                           in_tree = ex_tree)
  
  expected <- data.frame(intersect = c(0.7691243, 0.5128196),
                         slope = c(-0.01247047, 0.27961290),
                         p = c(0.9554795, 0.1885335),
                         row.names = c("K07106", "K02036"))
  
  expect_equal(observed, expected, tolerance = 1e-7)
  
})

