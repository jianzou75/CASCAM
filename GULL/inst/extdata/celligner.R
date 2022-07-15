library(tidyverse)
library(magrittr)

# Parameters
global <- list(
  n_genes = 'all', # set to 'all' to use all protein coding genes found in both datasets 
  umap_n_neighbors = 10, # num nearest neighbors used to create UMAP plot
  umap_min_dist = 0.5, # min distance used to create UMAP plot
  mnn_k_CL = 5, # number of nearest neighbors of tumors in the cell line data
  mnn_k_tumor = 50, # number of nearest neighbors of cell lines in the tumor data
  top_DE_genes_per = 1000, # differentially expressed genes with a rank better than this is in the cell line or tumor data
  # are used to identify mutual nearest neighbors in the MNN alignment step
  remove_cPCA_dims = c(1,2,3,4), # which cPCA dimensions to regress out of the data 
  distance_metric = 'euclidean', # distance metric used for the UMAP projection
  mod_clust_res = 5, # resolution parameter used for clustering the data
  mnn_ndist = 3, # ndist parameter used for MNN
  n_PC_dims = 70, # number of PCs to use for dimensionality reduction
  reduction.use = 'umap', # 2D projection used for plotting
  fast_cPCA = 10 # to run fast cPCA (approximate the cPCA eigenvectors instead of calculating all) set this to a value >= 4
)

create_Seurat_object <- function(exp_mat, ann, type = NULL) {
  seu_obj <- Seurat::CreateSeuratObject(t(exp_mat),
                                        min.cells = 0,
                                        min.features = 0,
                                        meta.data = ann %>%
                                          magrittr::set_rownames(ann$sampleID))
  if(!is.null(type)) {
    seu_obj@meta.data$type <- type
  }
  # mean center the data, important for PCA
  seu_obj <- Seurat::ScaleData(seu_obj, features = rownames(Seurat::GetAssayData(seu_obj)), do.scale = F)
  
  seu_obj %<>% Seurat::RunPCA(assay='RNA',
                              features = rownames(Seurat::GetAssayData(seu_obj)),
                              npcs = global$n_PC_dims, verbose = F)

  
  return(seu_obj)
}


cluster_data <- function(seu_obj) {
  seu_obj <- Seurat::FindNeighbors(seu_obj, reduction = 'pca',
                                   dims = 1:global$n_PC_dims,
                                   k.param = 20, 
                                   force.recalc = TRUE,
                                   verbose = FALSE)
  
  seu_obj %<>% Seurat::FindClusters(reduction = 'pca',
                                    resolution = global$mod_clust_res)
  
  seu_obj@meta.data$cluster <- seu_obj@meta.data$seurat_clusters
  
  return(seu_obj)
  
}

find_differentially_expressed_genes <- function(seu_obj) {
  n_clusts <- nlevels(seu_obj@meta.data$seurat_clusters)
  if (n_clusts > 2) {
    cur_DE_genes <- run_lm_stats_limma_group(
      t(Seurat::GetAssayData(seu_obj, assay='RNA', slot='scale.data')),
      seu_obj@meta.data %>% dplyr::select(seurat_clusters),
      limma_trend = TRUE) %>%
      dplyr::select(Gene, gene_stat = F_stat)
  } else if (n_clusts == 2) {
    cur_DE_genes <- run_lm_stats_limma(t(Seurat::GetAssayData(seu_obj, assay='RNA', slot='scale.data')),
                                       seu_obj@meta.data$cluster,
                                       limma_trend = TRUE) %>%
      dplyr::mutate(gene_stat = abs(t_stat)) %>%
      dplyr::select(Gene, gene_stat)
  } else {
    cur_DE_genes <- data.frame(Gene = colnames(seu_obj), gene_stat = NA)
  }
  
  return(cur_DE_genes)
  
}

run_cPCA <- function(TCGA_obj, CCLE_obj, pc_dims = NULL) {
  cov_diff_eig <- run_cPCA_analysis(t(Seurat::GetAssayData(TCGA_obj, assay='RNA', slot='scale.data')), 
                                    t(Seurat::GetAssayData(CCLE_obj, assay='RNA', slot='scale.data')), 
                                    TCGA_obj@meta.data, CCLE_obj@meta.data, pc_dims=pc_dims)
  return(cov_diff_eig) 
}


run_MNN <- function(CCLE_cor, TCGA_cor,  k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist, 
                    subset_genes) {
  mnn_res <- modified_mnnCorrect(CCLE_cor, TCGA_cor, k1 = k1, k2 = k2, ndist = ndist, 
                                 subset_genes = subset_genes)
  
  return(mnn_res)
}



Celligner <- function(tumor, cell){
  tumor_tr <- t(tumor); cell_tr <- t(cell)
  ann <- data.frame(sampleID = c(rownames(tumor_tr), rownames(cell_tr)),
                    type = c(rep('tumor', nrow(tumor_tr)), rep('CL', nrow(cell_tr))))
  
  tumor_ann <- dplyr::filter(ann, type=='tumor')
  cell_ann <- dplyr::filter(ann, type=='CL')
  
  tumor_obj <- create_Seurat_object(tumor_tr, tumor_ann, type = 'tumor')
  cell_obj <- create_Seurat_object(cell_tr, cell_ann, type = 'CL')
  
  tumor_obj <- cluster_data(tumor_obj)
  cell_obj <- cluster_data(cell_obj)
  
  tumor_DE_genes <- find_differentially_expressed_genes(tumor_obj)
  cell_DE_genes <- find_differentially_expressed_genes(cell_obj)
  
  DE_genes <- full_join(tumor_DE_genes, cell_DE_genes, by = 'Gene', suffix = c('_tumor', '_CL')) %>%
    mutate(
      tumor_rank = dplyr::dense_rank(-gene_stat_tumor),
      CL_rank = dplyr::dense_rank(-gene_stat_CL),
      best_rank = pmin(tumor_rank, CL_rank, na.rm=T))
  
  DE_gene_set <- DE_genes %>%
    dplyr::filter(best_rank < global$top_DE_genes_per) %>% .[['Gene']]
  
  cov_diff_eig <- run_cPCA(tumor_obj, cell_obj, global$fast_cPCA)
  
  if(is.null(global$fast_cPCA)) {
    cur_vecs <- cov_diff_eig$vectors[, global$remove_cPCA_dims, drop = FALSE]
  } else {
    cur_vecs <- cov_diff_eig$rotation[, global$remove_cPCA_dims, drop = FALSE]
  }
  
  rownames(cur_vecs) <- colnames(tumor_tr)
  tumor_cor <- resid(lm(t(tumor_tr) ~ 0 + cur_vecs)) %>% t()
  cell_cor <- resid(lm(t(cell_tr) ~ 0 + cur_vecs)) %>% t()
  
  mnn_res <- run_MNN(cell_cor, tumor_cor,  k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist, 
                     subset_genes = DE_gene_set)
  combined_mat <- rbind(mnn_res$corrected, cell_cor)
  
  comb_obj <- create_Seurat_object(combined_mat, ann)
  comb_obj <- cluster_data(comb_obj)
  
  return(comb_obj)
  
}