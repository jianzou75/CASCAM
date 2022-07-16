# 01: Celligner
# The data alignment step of analysis, and the normalized data are normalized and aligned to be comparable

#' Celligner: tumor and cell line model normalized RNASeq data alignment
#'
#' A two-step machine learning method to make normalized RNASeq data (log2 TPM preferred) from cell lines and tumors comparable.
#' To gain the benefits from the large sample size, it is recommended to merge the interested tumor samples with the public datasets such as TCGA.
#' The parameters are set as default and consistent with the Celligner paper, and the users can refer to the origional paper for more details.
#'
#' @param tumor a matrix or data frame of normalized tumor gene expression on log scale (log2-TPM preferred) with rows for genes and columns for samples.
#' @param cell a matrix or data frame of normalized cell line gene expression on log scale (log2-TPM preferred) with rows for genes and columns for samples.
#' @import Seurat
#'
#' @return a Seurat object with aligned data included. \code{Seurat::GetAssayData()} can be used to extract the aligned dataset.
#' @export
#' @references
#' Warren, Allison, Yejia Chen, Andrew Jones, Tsukasa Shibue, William C. Hahn, Jesse S. Boehm, Francisca Vazquez, Aviad Tsherniak, and James M. McFarland. 2021. “Global Computational Alignment of Tumor and Cell Line Transcriptional Profiles.” Nature Communications 12 (1): 22. https://doi.org/10.1038/s41467-020-20294-x.
Celligner <- function(tumor, cell){
  source(system.file(package = 'GULL', 'extdata/Celligner_helpers.R'))
  source(system.file(package = 'GULL', 'extdata/analysis_helpers.R'))
  source(system.file(package = 'GULL', 'extdata/celligner.R'))

  tumor <- data.frame(tumor); cell <- data.frame(cell)

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
