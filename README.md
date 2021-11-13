# Phylogenetic Organization of Metagenomic Signals (POMS)

This repository contains the code used for the joint taxa-function bioinformatics analysis tool described in Douglas et al. (_In Prep_).

The available POMS code is currently an **alpha release**, meaning that it is subject to major changes and is still being tested for bugs.

### Installation

```R
library(devtools)
install_github("gavinmdouglas/POMS", ref = "main")
```

### Example usage

The key POMS function is `POMS_pipeline`, which requires tables of taxa and function abundances, a tree of taxa, and sample names split into group1 and group2. The below code will run the pipeline on example input files that are part of this repository. These example files are highly simplified to run quickly - for instance, there are only two functions tested and only 60 tips in the tree.

_Note: the options `min_num_tips`, `multinomial_min_sig`, and `min_func_instances` are set to non-default values to work with this small example, but normally you could leave them to be default._

```R
setwd("path/to/POMS/example_files/")

library(ape)
library(POMS)

ex_taxa_abun <- read.table("ex_taxa_abun.tsv.gz", header = TRUE, sep = "\t", row.names = 1)
ex_func <- read.table("ex_func.tsv.gz", header = TRUE, sep = "\t", row.names = 1)
ex_tree <- read.tree("ex_tree.newick")
ex_group1 <- read.table("ex_group1.txt.gz", stringsAsFactors = FALSE)$V1
ex_group2 <- read.table("ex_group2.txt.gz", stringsAsFactors = FALSE)$V1


# Example of how to run main POMS function. 
POMS_out <- POMS_pipeline(abun = ex_taxa_abun,
                          func = ex_func,
                          phylogeny = ex_tree,
                          group1_samples = ex_group1,
                          group2_samples = ex_group2,
                          ncores = 1,
                          min_num_tips = 4,
                          multinomial_min_sig = 3,
                          min_func_instances = 0,
                          verbose = TRUE)
```

### Example output

The above usage example will produce a list called `POMS_out`. This can be a very large list full of intermediate objects when `detailed=TRUE`. However, most users will just be interested in the output results dataframe. Usually there will be many rows to this dataframe, but in this case there are only two because that's how many functions were tested. Accordingly, we can look at the transposed version of this dataframe to see it better.

```R
t(POMS_out$df)

#                             K07106 K02036
# num_nodes_enriched               3 5.0000
# num_sig_nodes_group1_enrich      2 5.0000
# num_sig_nodes_group2_enrich      1 0.0000
# num_nonsig_nodes_enrich          0 0.0000
# multinomial_p                    1 0.0272
# multinomial_corr                 1 0.0544
```

* `num_nodes_enriched`: Number of nodes with differential enrichment of the function in one subtree (the next three categories sum to this count)
* `num_sig_nodes_group1_enrich`: Number of nodes where the function is enriched in the subtree that is relatively higher in group1 samples.
* `num_sig_nodes_group2_enrich`: Number of nodes where the function is enriched in the subtree that is relatively higher in group2 samples.
* `num_nonsig_nodes_enrich`: Number of nodes where the function is enriched and there is no difference in sample balances.
* `multinomial_p` and `multinomial_corr` P and corrected P-value (`BH` by default) based on multinomial test.

The other default output elements in `POMS_out` are:
* `balances_info`: Summary of the balance tree information, including the balances at each node and lists of the tip names within each subtree below each node.
* `sig_nodes`: Nodes that are significantly different based on sample balances.
* `multinomial_exp_prop`: Expected proportions of the three function enrichment categories based on the counts of significant nodes based on functions and balances.
