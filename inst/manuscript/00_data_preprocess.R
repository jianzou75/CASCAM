library(CASCAM)
load("../20220713_Software/CASCAM_demo/BRCA_tumor_count.RData")
load("../20220713_Software/CASCAM_demo/BRCA_tumor_cell_aligned_DE_filter.RData")

gene_info <- create_InformativeGenes(brca_ct, brca_label2, "ILC") ## obtain the DE and GSEA information from the tumor read counts.
CASCAM_ILC <- create_CASCAM(tumor_aligned, tumor_label, cell_aligned, gene_info) ## create an CASCAM object for downstream analysis
