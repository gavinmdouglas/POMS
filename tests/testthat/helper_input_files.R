
cat(getwd(), file = "/home/gavin/tmp/test")

ex_taxa_abun <- read.table("../example_files/ex_taxa_abun.tsv.gz", header = TRUE, sep = "\t", row.names = 1)
ex_func <- read.table("../example_files/ex_func.tsv.gz", header = TRUE, sep = "\t", row.names = 1)

ex_group1 <- read.table("../example_files/ex_group1.txt.gz", stringsAsFactors = FALSE)$V1
ex_group2 <- read.table("../example_files/ex_group2.txt.gz", stringsAsFactors = FALSE)$V1

ex_tree <- ape::read.tree("../example_files/ex_tree.newick")

ex_taxa_labels <- read.table("../example_files/ex_taxa_labels.tsv.gz",
                             header = TRUE,
                             sep = "\t",
                             stringsAsFactors = FALSE,
                             row.names = 1)
