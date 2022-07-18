#' @import fgsea
#' @import DESeq2
#' @import sparseLDA
NULL


#' The InformativeGenes class
#'
#' @slot DEA A data frame of the results of differential gene expression analysis.
#' @slot GSEA A data frame of the results of gene set enrichment analysis.
#' @slot interested_subtype A string indicating the interested tumor subtype.
#'
#' @return
#' @exportClass InformativeGenes
#'
#' @examples
InformativeGenes <- methods::setClass(
  "InformativeGenes",
  slots = c(
    DEA = "data.frame",
    GSEA = "data.frame",
    interested_subtype = "character"
    )
  )


#' Create an InformativeGenes object
#'
#' This function runs differential expression analysis (DESeq2) and gene set enrichment analysis (GSEA) to generate a InformativeGenes object,
#' which is necessary for creating the CASCAM object.
#'
#' @param tumor_ct a data frame or matrix of the un-normalized (estimated) counts of sequencing reads or fragments from the tumor data with
#'        Rows for genes with row names of HGNC gene symbols and columns for samples.
#' @param tumor_lab a vector of strings indicating the labels of the tumor samples.
#' @param interested_subtype a string indicating the interested tumor subtype, and it should be included in the tumor labels.
#'
#' @import DESeq2
#' @import fgsea
#' @import qusage
#'
#' @return An \code{InformativeGenes} object.
#'
#' @export
#'
#' @references
#' Love, Michael I., Wolfgang Huber, and Simon Anders. "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome biology 15, no. 12 (2014): 1-21.
create_InformativeGenes <- function(tumor_ct, tumor_lab, interested_subtype) {
  if (!interested_subtype %in% unique(tumor_lab)) {
    stop("Interested subtype should be included in the tumor labels.")
  }

  interested_level <- which(interested_subtype == unique(tumor_lab))
  uninterested_level <- which(interested_subtype != unique(tumor_lab))
  uninterested_subtype <- setdiff(unique(tumor_lab), interested_subtype)
  tumor_clin <- data.frame(tumor_lab = factor(tumor_lab, levels = c(uninterested_subtype, interested_subtype)))

  dds <- DESeq2::DESeqDataSetFromMatrix(countData = tumor_ct,
                                        colData = tumor_clin,
                                        design = ~tumor_lab)
  dds <- DESeq2::DESeq(dds)
  de <- data.frame(DESeq2::results(dds)[, c("padj", "log2FoldChange")])
  de$gene_symbol <- rownames(de)
  de$diffexpressed <- "NO"
  de$diffexpressed[de$log2FoldChange > log2(1.5) & de$padj < 0.05] <- "UP"
  de$diffexpressed[de$log2FoldChange < -log2(1.5) & de$padj < 0.05] <- "DOWN"
  de$delabel <- NA
  de$delabel[de$diffexpressed != "NO"] <- de$gene_symbol[de$diffexpressed != "NO"]

  pathways <- c(qusage::read.gmt(system.file(package = "CASCAM", "extdata/c2.cp.kegg.v7.4.symbols.gmt")),
                qusage::read.gmt(system.file(package = "CASCAM", "extdata/h.all.v7.4.symbols.gmt")))
  gene_rank <- de$log2FoldChange
  names(gene_rank) <- de$gene_symbol
  fgseaRes <- data.frame(fgsea::fgsea(pathways, gene_rank, minSize = 15, maxSize = 500))

  object <- methods::new(Class = "InformativeGenes",
                         DEA = de, GSEA = fgseaRes, interested_subtype = as.character(interested_subtype))

  return(object)
}


#' The CASCAM class
#'
#' @slot tumor_aligned_data A data frame of aligned tumor gene expression with rows for genes and columns for samples.
#' @slot tumor_label A vector of strings indicating the labels of the tumor samples.
#' @slot camod_aligned_data A data.frame of aligned cancer model gene expression with rows for genes and columns for samples.
#' @slot sda_model A list for the SDA model trained by the tumor data.
#' @slot tumor_norm_data A data frame of normalized \emph{tumor_aligned_data} with rows for samples and columns for genes.
#' @slot camod_norm_data A data frame of normalized \emph{camod_aligned_data} with rows for samples and columns for genes.
#' @slot tumor_sda_project A one column matrix showing the SDA projected vector of tumor gene expression.
#' @slot camod_sda_project A one column matrix showing the SDA projected vector of cancer model gene expression.
#' @slot sda_ds A data frame for SDA deviance score with the rows for camods and columns for tumor subtypes.
#' @slot sda_lds A data frame for SDA log deviance score with the rows for camods and columns for tumor subtypes.
#' @slot sda_predict_prob A two column matrix showing the SDA assignment probabilities for camods.
#' @slot sda_ds_pval A two column matrix showing the p-values of SDA deviance scores for camods.
#' @slot sda_lds_ci A list of two matrices showing the 95\% confidence intervals of \emph{sda_lds}.
#' @slot genome_figure A list of two ggplot2 figures for the genome-wide pre-selection (distribution and confidence interval).
#' @slot selected_camods The vector of genome-wide selected cancer models representing the interested tumor subtype by \emph{sda_predict_prob} and \emph{sda_ds_pval}.
#' @slot gene_ds A list of two matrices showing the gene-specific deviance score.
#' @slot pathway_ds A list of two matrices showing the pathway-specific deviance score.
#' @slot available_pathways A vector of available pathways.
#' @slot pathway_congruence_heatmap A ggplot2 object for the pathway congruence heatmap.
#' @slot pathway_gene_heatmap A ggplot2 object for the pathway specific heatmap.
#' @slot pathway_gene_ridgeline A ggplot2 object for the pathway specific gene expression ridgeline.
#' @slot pathview_figure A pathview figure for one specific pathway and one specific cancer model.
#'
#' @return
#' @export
#'
#' @examples
CASCAM <- methods::setClass(
  "CASCAM",
  slots = c(
    tumor_aligned_data = "data.frame",
    tumor_label = "character",
    camod_aligned_data = "data.frame",
    sda_model = "ANY",
    tumor_norm_data = "data.frame",
    camod_norm_data = "data.frame",
    tumor_sda_project = "matrix",
    camod_sda_project = "matrix",
    sda_ds = "data.frame",
    sda_lds = "data.frame",
    sda_predict_prob = "matrix",
    sda_ds_pval = "matrix",
    sda_lds_ci = "list",
    genome_figure = "list",
    selected_camods = "character",
    gene_ds = "list",
    pathway_ds = "list",
    available_pathways = "character",
    pathway_congruence_heatmap = "ANY",
    pathway_gene_heatmap = "ANY",
    pathway_gene_ridgeline = "ANY",
    pathview_figure = "ANY"
  ),
  contains = c("InformativeGenes")
)


#' Create a CASCAM object
#'
#' This function creates a \code{CASCAM} object for the downstream analysis. The data sets of tumor aligned gene expression data with the tumor subtype labels and
#' cancer model aligned gene expression data are needed as the input. Celligner is highly recommended to align the tumor and cancer model data.
#'
#' @param tumor_aligned_data A data frame of aligned tumor gene expression with rows for genes and columns for samples.
#' @param tumor_label A vector of strings indicating the labels of the tumor samples.
#' @param camod_aligned_data A data.frame of aligned cancer model gene expression with rows for genes and columns for samples.
#' @param info_object An \code{InformativeGenes} object for the study.
#'
#' @return A \code{CASCAM} object.
#' @export
#'
#' @examples
create_CASCAM <- function(tumor_aligned_data, tumor_label, camod_aligned_data, info_object) {
  DEG <- c(na.omit(info_object@DEA$delabel))
  DEG <- Reduce(intersect, list(DEG, rownames(tumor_aligned_data), rownames(camod_aligned_data)))
  tumor_aligned_data <- tumor_aligned_data[DEG,]
  camod_aligned_data <- camod_aligned_data[DEG,]

  object <- methods::new(Class = "CASCAM",
                         tumor_aligned_data = tumor_aligned_data, tumor_label = tumor_label,
                         camod_aligned_data = camod_aligned_data, info_object)

  return(object)
}
