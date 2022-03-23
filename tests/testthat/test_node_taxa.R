
library(POMS)

ex1_taxa_labels <- read.table("../../example_files/ex_taxa_labels.tsv.gz", header = TRUE, sep = "\t", row.names = 1)
ex1_tree <- ape::read.tree("../../example_files/ex_tree.newick")


test_that("Check taxa on each side of the node - combined labels and threshold of 0.51", {
 
  exp_output <- c("Bacteria; Bacteroidetes; Bacteroidia; Bacteroidales; Porphyromonadaceae (Family)",
                  "Bacteria; Bacteroidetes; Bacteroidia; Bacteroidales; Porphyromonadaceae; NA (Genus)")

  obs_output <- node_taxa(in_tree = ex1_tree, taxon_labels = ex1_taxa_labels,
                          node_label = "n1", combine_labels = TRUE, threshold = 0.51)
  
  expect_equal(obs_output, exp_output)

})


test_that("Check taxa on each side of the node - simple (non-combined) labels and threshold of 0.95", {
  
  exp_output <- c("Porphyromonadaceae (Family)", "Bacteroidales (Order)")
  
  obs_output <- node_taxa(in_tree = ex1_tree, taxon_labels = ex1_taxa_labels,
                          node_label = "n1", combine_labels = FALSE, threshold = 0.95)
  
  expect_equal(obs_output, exp_output)
  
})
