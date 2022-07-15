##################################################################################################################################
# 04: Pathway selection
# The genome-wide pre-selected cell lines are used for further pathway-specific analysis

#' Pathway specific selection analysis
#'
#' @param object A \code{GULL} object.
#'
#' @return
#' @export
#'
#' @examples
pathway_analysis <- function(object){
  geometric.mean <- function(x) {exp(mean(log(x)))}

  subtype_levels <- levels(factor(object@tumor_label))
  tumor_expr <- t(object@tumor_norm_data)
  cell_expr <- t(object@cell_norm_data)[,object@selected_cells]
  deg <- rownames(tumor_expr)

  ## gene specific deviance score
  center_location <- sapply(subtype_levels, function(t) apply(tumor_expr[,object@tumor_label == t], 1, function(x) median(x)))
  pool_sd <- apply(tumor_expr, 1, function(y) mad(unlist(sapply(subtype_levels, function(t) y[object@tumor_label == t] - median(y[object@tumor_label == t])))))
  gene_ds <- lapply(subtype_levels, function(t) apply(cell_expr, 2, function(x) (x - center_location[,t])/pool_sd)); names(gene_ds) = subtype_levels

  ## pathway specific deviance score
  pathways <- c(qusage::read.gmt(system.file(package = 'GULL', 'extdata/c2.cp.kegg.v7.4.symbols.gmt')),
                qusage::read.gmt(system.file(package = 'GULL', 'extdata/h.all.v7.4.symbols.gmt')))
  pathway_ds_subtype1 <- pathway_ds_subtype2 <- c();                       ## PDS : pathway deviance score
  path_row_names <- c()
  for (i in 1:length(pathways)){
    if (sum(deg %in% pathways[[i]]) > 20){
      gene_ds_subtype1 <- gene_ds[[subtype_levels[[1]]]][intersect(deg, pathways[[i]]), ]
      gene_ds_subtype2 <- gene_ds[[subtype_levels[[2]]]][intersect(deg, pathways[[i]]), ]
      pathway_ds_subtype1 <- rbind(pathway_ds_subtype1, c(apply(abs(gene_ds_subtype1), 2, geometric.mean)))
      pathway_ds_subtype2 <- rbind(pathway_ds_subtype2, c(apply(abs(gene_ds_subtype2), 2, geometric.mean)))
      path_row_names<- c(path_row_names, names(pathways)[i])
    }
  }
  rownames(pathway_ds_subtype1) <- rownames(pathway_ds_subtype2) <- path_row_names
  colnames(pathway_ds_subtype1) <- colnames(pathway_ds_subtype2) <- colnames(cell_expr)
  pathway_ds <- list(pathway_ds_subtype1, pathway_ds_subtype2)
  names(pathway_ds) <- subtype_levels

  object@gene_ds <- gene_ds
  object@pathway_ds <- pathway_ds
  object@available_pathways <- path_row_names

  return(object)
}


#' Pathway congruence heatmap
#'
#' A heatmap showing the congruence between the cell lines and tumors for detailed pathways.
#'
#' @param object A \code{GULL} object.
#'
#' @return
#' @export
#'
#' @examples
pathway_congruence_heatmap <- function(object){
  fgseaRes <- object@GSEA[match(object@available_pathways, object@GSEA$pathway), ]
  pathway_annotate <- data.frame(NES =  cut(as.numeric(fgseaRes$NES), seq(-2.5, 2.5,length.out = 5), labels = c(-2, -1, 1, 2)),
                                 pathway_size = fgseaRes$size); row.names(pathway_annotate) = fgseaRes$pathway

  sample_annotate <- data.frame(SDA_DS = object@sda_ds[colnames(object@pathway_ds[[object@interested_subtype]]), object@interested_subtype])
  row.names(sample_annotate) <- colnames(object@pathway_ds[[object@interested_subtype]])

  pathway_congruence_heatmap <- pheatmap::pheatmap(object@pathway_ds[[object@interested_subtype]],
                                                   color = c("#DA2C43", "#E15566", "#E97E88", "#F0A8AB", "#F8D1CD", "#FFFAF0"),
                                                   breaks = seq(0, 2, length.out = 7),
                                                   annotation_row = pathway_annotate,
                                                   annotation_col = sample_annotate,
                                                   annotation_colors = list(NES = c("2" = "#0025F6", "1" ="grey", "-1" = "grey", "-2" ="#E60000"),
                                                                            pathway_size = c("#FAFA31", "#FFB131"),
                                                                            SDA_DS = c("#388230", "#EFED8F")),
                                                   clustering_method = "ward.D2", main = "Pathway congruence heatmap")
  object@pathway_congruence_heatmap <- pathway_congruence_heatmap

  return(object)
}


#' Pathway specific heatmap
#'
#' @param object A \code{GULL} object.
#' @param pathway_name The interested pathway, which should be within the KEGG and Hallmark pathways.
#'
#' @import reshape2
#' @import ggridges
#'
#' @return
#' @export
#'
#' @examples
pathway_specific_heatmap <- function(object, pathway_name){
  allcolour = c("#3AC9B0","#F2C935","#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
                "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
                "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

  pathways <- c(qusage::read.gmt(system.file(package = 'GULL', 'extdata/c2.cp.kegg.v7.4.symbols.gmt')),
                qusage::read.gmt(system.file(package = 'GULL', 'extdata/h.all.v7.4.symbols.gmt')))
  gene_pathway <- intersect(rownames(object@gene_ds[[object@interested_subtype]]),
                            pathways[[pathway_name]])
  gene_pathway <- gene_pathway[order(object@DEA[gene_pathway, 1])]
  uninterested_subtype <- setdiff(levels(factor(object@tumor_label)), object@interested_subtype)

  pathway_annotate <- data.frame(DS_pathway = object@pathway_ds[[object@interested_subtype]][pathway_name,])
  pathway_cell <- object@gene_ds[[object@interested_subtype]][gene_pathway, ]
  pathway_gene_heatmap <- pheatmap::pheatmap(abs(pathway_cell),
                                             color = c("#DA2C43", "#E15566", "#E97E88", "#F0A8AB", "#F8D1CD", "#FFFAF0"),
                                             breaks = seq(0, 3, length.out = 7),
                                             annotation_col = pathway_annotate,
                                             annotation_colors = list(DS_pathway = c("#388230", "#EFED8F")),
                                             clustering_method = "ward.D2",
                                             main = paste0("Heatmap for absolute value of DS_Gene in ", pathway_name))

  object@pathway_gene_heatmap <- pathway_gene_heatmap

  return(object)
}


#' Pathway specific gene expression ridgeline
#'
#' @param object A \code{GULL} object.
#' @param interested_cells A vector of interested cell lines.
#' @param pathway_name The interested pathway, which should be within the KEGG and Hallmark pathways.
#'
#' @importFrom reshape2 melt
#' @import ggridges
#'
#' @return
#' @export
#'
#' @examples
pathway_specific_ridgeline <- function(object, interested_cells, pathway_name){
  if(length(interested_cells) > 5){
    warning("It is highly recommended to list less than 5 cells, otherwise the figure will be too busy to read.")
  }

  allcolour = c("#3AC9B0","#F2C935","#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
                "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
                "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
                "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

  pathways <- c(qusage::read.gmt(system.file(package = 'GULL', 'extdata/c2.cp.kegg.v7.4.symbols.gmt')),
                qusage::read.gmt(system.file(package = 'GULL', 'extdata/h.all.v7.4.symbols.gmt')))
  gene_pathway <- intersect(rownames(object@gene_ds[[object@interested_subtype]]),
                            pathways[[pathway_name]])
  gene_pathway <- gene_pathway[order(object@DEA[gene_pathway, 1])]
  uninterested_subtype <- setdiff(levels(factor(object@tumor_label)), object@interested_subtype)

  tumor_pathway <- cbind.data.frame(label = object@tumor_label, object@tumor_norm_data[,gene_pathway], stringsAsFactors = FALSE)
  cell_pathway  <- object@cell_norm_data[, gene_pathway]
  tumor_pathway2 <- data.frame(reshape2::melt(tumor_pathway, id.vars = "label"))
  cell_pathway2 <- reshape2::melt(as.matrix(cell_pathway))
  pathway_gene_ridgeline <- ggplot() +
    geom_density_ridges(data =  tumor_pathway2, aes(x = value, y = variable, color = label), scale = 1,  alpha = 0, size = 1, quantile_lines = TRUE, quantiles = 2) +
    scale_colour_manual(name = "tumor distribution", breaks = c(uninterested_subtype, object@interested_subtype), values = c("#BBE7FE", "#FFAEBC")) +
    theme_bw() +
    geom_point(position = position_jitter(seed = 1, width = 0, height = 0), data = cell_pathway2[cell_pathway2$Var1 %in% interested_cells,], aes(x = value, y = Var2, fill = Var1), size = 3, shape=21, color = "white") +
    scale_fill_manual(name = "cell line", values = allcolour[1:length(interested_cells)]) +
    ggtitle(paste0(pathway_name, " Gene Distribution")) +
    labs(x = "Expression", y = "Gene") +
    xlim(round(range(cell_pathway2[cell_pathway2$Var1 %in% interested_cells,]$value), 2)) +
    ggtitle(paste0(pathway_name, " gene distribution"))

  object@pathway_gene_ridgeline <- pathway_gene_ridgeline

  return(object)
}


#' Pathview while not saving as files
#'
#' @param ...
#' @param save_image
#'
#' @import grid
#' @import pathview
#' @import png
#' @import ggplotify
#'
#' @return
#' @export
#'
#' @keywords internal
see_pathview <- function(..., save_image = FALSE){
  msg <- capture.output(pathview::pathview(...), type = "message")
  msg <- grep("image file", msg, value = T)
  filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
  img <- png::readPNG(filename)
  p <- as.ggplot(function() grid::grid.raster(img))
  if(!save_image) invisible(file.remove(filename))
  return(p)
}


#' Pathview of the interested cell in interested pathway
#'
#' @param object A \code{GULL} object.
#' @param interested_cell A string of one interested cell line.
#' @param pathway_name The interested pathway, which should be within the KEGG and Hallmark pathways.
#'
#' @import pathview
#'
#' @return
#' @export
pathview_analysis <- function(object, interested_cell, pathway_name){
  load(system.file(package = 'GULL', 'extdata/KEGG_name_ID_match.RData'))
  kegg_id = KEGG_name_ID_match$gs_exact_source[KEGG_name_ID_match$gs_name == pathway_name]

  pathview_input <- object@gene_ds[[object@interested_subtype]][,interested_cell]
  pathview_figure <- see_pathview(gene.data = pathview_input, pathway.id = kegg_id,
                                  species = "hsa", out.suffix = "", gene.idtype = "SYMBOL",
                                  limit=list(gene=c(-3,3)), bins = list(gene = 16),
                                  low = list(gene = "#B2E6F0"), mid = list(gene = "#FF0018"), high = list(gene = "#FFF485"),
                                  node.sum = "mean")
  invisible(file.remove(paste0("hsa",kegg_id,".png")))
  invisible(file.remove(paste0("hsa",kegg_id,".xml")))

  object@pathview_figure <- pathview_figure

  return(object)
}


