rm(list=ls(all=TRUE))

setwd("/home/gavin/projects/functional_reference_frames/MAGs/Almeida2019/")
source("../../scripts/balance_tree_functions.R")

# Read in input files.
almeida_func <- read.table(file = "functional_analyses/kegg_summary.csv", sep=",", stringsAsFactors = FALSE, quote="", comment.char = "",
                           header=TRUE, check.names = FALSE, row.names=1)

almeida_func <- data.frame(t(almeida_func), check.names = FALSE)

almeida_tree <- read.tree(file = "phylogenies/raxml_hgr-umgs_phylogeny.nwk")
almeida_abun <- read.table(file = "mapping_results/modified/bwa_depth_min25coverage.tsv", header=TRUE, sep="\t", check.names=FALSE,
                           row.names=1, quote="", comment.char="")

almeida_sample_info <- read.table("MGS_samples_info_SuppTable1.txt", header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "")

almeida_taxa_umgs <- read.table("taxonomy/taxonomy_umgs.tab", header=TRUE, sep="\t", stringsAsFactors = FALSE)
rownames(almeida_taxa_umgs) <- almeida_taxa_umgs$MAG_ID
almeida_taxa_umgs <- almeida_taxa_umgs[, -which(colnames(almeida_taxa_umgs) %in% c("UMGS_ID", "MAG_ID"))]

almeida_taxa_hgr <- read.table("taxonomy/taxonomy_hgr.tab", header=TRUE, sep="\t", stringsAsFactors = FALSE)
rownames(almeida_taxa_hgr) <- almeida_taxa_hgr$Genome
almeida_taxa_hgr <- almeida_taxa_hgr[, -which(colnames(almeida_taxa_hgr) == "Genome")]

almeida_taxa <- rbind(almeida_taxa_hgr, almeida_taxa_umgs)


# Obesity
ERP002061_almeida_sample_info <- almeida_sample_info[which(almeida_sample_info$Study == "ERP002061"), ]


ERP002061_almeida_abun <- subset_abun_table(in_abun = almeida_abun, col2keep = ERP002061_almeida_sample_info$Run)

ERP002061_almeida_sample_info <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Run %in% colnames(ERP002061_almeida_abun)), ]
ERP002061_group1_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state == "Diseased"), "Run"]
ERP002061_group2_samples <- ERP002061_almeida_sample_info[which(ERP002061_almeida_sample_info$Health.state == "Healthy"), "Run"]

ERP002061_almeida_func <- almeida_func[rownames(ERP002061_almeida_abun), ]

ERP002061_almeida_out <- two_group_balance_tree_pipeline(abun=ERP002061_almeida_abun,
                                                         func=ERP002061_almeida_func,
                                                         phylogeny=almeida_tree,
                                                         taxa=almeida_taxa,
                                                         group1_samples = ERP002061_group1_samples,
                                                         group2_samples = ERP002061_group2_samples,
                                                         ncores=30)



