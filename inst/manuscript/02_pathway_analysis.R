set.seed(12345)
library(aricode);library(glmnet); library(caret); library(xtable); library(ggrepel);
library(tidyverse); library(ggpubr); library(patchwork); library(boot);
library(DESeq2); library(openxlsx); library(RColorBrewer); library(pathview); library(fgsea)
library(ggridges); library(data.table)
library(grid);library(gridExtra)

tumor_aligned <- CASCAM_ILC@tumor_aligned_data 
tumor_label <- CASCAM_ILC@tumor_label
cell_aligned <- CASCAM_ILC@camod_aligned_data
train_norm = sparseLDA::normalize(t(tumor_aligned))
test_norm  = sparseLDA::normalizetest(t(cell_aligned), train_norm)
train_data = t(train_norm$Xc)
train_label = tumor_label
label_level = c("IDC", "ILC")
test_data = t(test_norm)
deg <- rownames(train_data)
de <- CASCAM_ILC@DEA


## Supplementary figure 1 for complete heatmap
center_location <- sapply(label_level, function(t) apply(train_data[,train_label == t], 1, function(x) median(x)))
pool_sd <- apply(train_data, 1, function(y) mad(unlist(sapply(label_level, function(t) y[train_label == t] - median(y[train_label == t])))))

gene_devscore <- lapply(label_level, function(t) apply(test_data, 2, function(x) (x - center_location[,t])/pool_sd)); names(gene_devscore) = label_level

geometric.mean <- function(x) {exp(mean(log(x)))}

kegg_pathways <- c(qusage::read.gmt(system.file("extdata", "c2.cp.kegg.v7.4.symbols.gmt", package="CASCAM")),
                   qusage::read.gmt(system.file("extdata", "h.all.v7.4.symbols.gmt", package="CASCAM")))
pathway_pds_ILC <- c();                       ## PDS : Aggregation deviance score
path.row.names <- c()

for (i in 1:length(kegg_pathways)){
  if (sum(deg %in% kegg_pathways[[i]]) > 20){
    gds_ILC <- gene_devscore$ILC[intersect(deg, kegg_pathways[[i]]), ]
    pathway_pds_ILC <- rbind(pathway_pds_ILC, c(apply(abs(gds_ILC), 2, geometric.mean)))
    
    path.row.names<- c(path.row.names, names(kegg_pathways)[i])
  }
}
rownames(pathway_pds_ILC)  <- path.row.names
colnames(pathway_pds_ILC)  <- colnames(test_data)

gene_rank <- de$log2FoldChange
names(gene_rank) <- de$gene_symbol
fgseaRes <- fgsea(kegg_pathways[path.row.names], na.omit(gene_rank), minSize=15, maxSize = 500)

pathway_annot <- data.frame(NES =  cut(fgseaRes$NES, seq(-2.5,2.5,length.out = 5), labels = c(-2, -1, 1, 2)),
                            pathway_size = fgseaRes$size); row.names(pathway_annot) = fgseaRes$pathway
pathway_annot <- pathway_annot[path.row.names,]
sample_annot  <- data.frame(SDA_DS = genome.ci.table[(genome.ci.table$classification_P_LDS == "ILC" & genome.ci.table$SDA_ILC_P > 0.8)|(genome.ci.table$cell == "MDAMB134VI"), ]$SDA_ILC_DS);
row.names(sample_annot) = genome.ci.table[(genome.ci.table$classification_P_LDS == "ILC" & genome.ci.table$SDA_ILC_P > 0.8)|(genome.ci.table$cell == "MDAMB134VI"),]$cell

supp_figure1 <- pheatmap::pheatmap(pathway_pds_ILC[,row.names(sample_annot)],
                                   color = c("#d73027", "#fc8d59", "#fee090", "#ffffbf", "#e0f3f8", "#91bfdb", "#4575b4"),
                                   fontsize = 12, fontsize_row = 10, fontsize_col = 10,
                                   breaks = seq(0, 2, length.out = 7),
                                   annotation_row = pathway_annot,
                                   annotation_col = sample_annot,
                                   annotation_colors = list(NES = c("2" = "#0025F6", "1" ="grey", "-1" = "grey", "-2" ="#E60000"),
                                                            pathway_size = c("#FAFA31", "#FFB131"),
                                                            SDA_DS = c("#388230", "#EFED8F")),
                                   clustering_method = "ward.D2",  main = "")


## Figure 5A for the heatmap of selected pathways
path.row.names2 <- fgseaRes$pathway[abs(fgseaRes$NES) > 1.5 & fgseaRes$size < 200]
path.row.names2 <- path.row.names2[order(fgseaRes$NES[abs(fgseaRes$NES) > 1.5 & fgseaRes$size < 200], decreasing = T)]
path.row.names2 <- c(path.row.names2, "KEGG_CELL_ADHESION_MOLECULES_CAMS")
pathway_annot <- data.frame(NES =  cut(fgseaRes$NES, seq(-2.5,2.5,length.out = 5), labels = c(-2, -1, 1, 2)),
                            pathway_size = fgseaRes$size); row.names(pathway_annot) = fgseaRes$pathway
pathway_annot <- pathway_annot[path.row.names2,]
pathway_heatmap = pathway_pds_ILC[row.names(pathway_annot),row.names(sample_annot)]
pathway_heatmap = rbind(pathway_heatmap, combined_pathway = colMeans(pathway_heatmap[1:14,]))

figure5A <- pheatmap::pheatmap(pathway_heatmap[,rownames(sample_annot)[order(sample_annot$SDA_DS, decreasing = F)]],
                               color = c("#d73027", "#fc8d59", "#fee090", "#ffffbf", "#e0f3f8", "#91bfdb", "#4575b4"),
                               fontsize = 12, fontsize_row = 14, fontsize_col = 14,
                               breaks = seq(0, 1.6, length.out = 7),
                               annotation_row = rbind(pathway_annot, combined_pathway = c(NA, NA)),
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               annotation_col = sample_annot,
                               annotation_colors = list(NES = c("2" = "#0025F6", "1" ="grey", "-1" = "grey", "-2" ="#E60000"),
                                                        pathway_size = c("#FAFA31", "#FFB131"),
                                                        SDA_DS = c("#388230", "#EFED8F")),
                               clustering_method = "ward.D2",  main = "")

gene_devscore_tumor <- lapply(label_level, function(t) apply(train_data[,train_label == "ILC"], 2, function(x) (x - center_location[,t])/pool_sd)); names(gene_devscore_tumor) = label_level
geometric.mean <- function(x) {exp(mean(log(x)))}
pathway_pds_ILC_tumor <- c();                       
path.row.names <- c()

for (i in 1:length(kegg_pathways)){
  if (sum(deg %in% kegg_pathways[[i]]) > 20){
    gds_ILC <- gene_devscore_tumor$ILC[intersect(deg, kegg_pathways[[i]]), ]
    pathway_pds_ILC_tumor <- rbind(pathway_pds_ILC_tumor, c(apply(abs(gds_ILC), 2, geometric.mean)))
    
    path.row.names<- c(path.row.names, names(kegg_pathways)[i])
  }
}
rownames(pathway_pds_ILC_tumor)  <- path.row.names
colnames(pathway_pds_ILC_tumor)  <- colnames(train_data)[train_label == "ILC"]

pathway_pds_ILC_tumor <- pathway_pds_ILC_tumor[rownames(pathway_heatmap)[1:15], ]
pathway_pds_ILC_tumor = rbind(pathway_pds_ILC_tumor, combined_pathway = colMeans(pathway_pds_ILC_tumor[1:14,]))
tumor_path_ds_ecdf <- apply(pathway_pds_ILC_tumor, 1, ecdf)
pathway_heatmap_pvalue <- sapply(1:16, function(i) tumor_path_ds_ecdf[[i]](pathway_heatmap[i,]))
pathway_heatmap_pvalue <- t(1 - pathway_heatmap_pvalue)
rownames(pathway_heatmap_pvalue) <- rownames(pathway_heatmap) ## Table of p-values for the heatmap of above selected pathways
colnames(pathway_heatmap_pvalue) <- colnames(pathway_heatmap)


## Supplementary Figure 2 for enrichment score
E2F = plotEnrichment(kegg_pathways[["HALLMARK_E2F_TARGETS"]], na.omit(gene_rank)) + labs(title="Cell cycle related targets of E2F")
PPAR = plotEnrichment(kegg_pathways[["KEGG_PPAR_SIGNALING_PATHWAY"]], na.omit(gene_rank)) + labs(title="PPAR signaling pathway")
supp_figure2 <- plot_grid(E2F, PPAR, nrow = 2, rel_heights = c(1, 1))


## Figure 5B for cell adhesion heatmap
adhesion_annot = data.frame(DS_pathway = pathway_pds_ILC["KEGG_CELL_ADHESION_MOLECULES_CAMS",row.names(sample_annot)]); row.names(adhesion_annot) = row.names(sample_annot)
adhesion_cell <- gene_devscore$ILC[intersect(deg, kegg_pathways$KEGG_CELL_ADHESION_MOLECULES_CAMS), row.names(sample_annot)]
sequence_gene <- c("CDH15", "CADM3", "SELP", "CLDN19", "NRXN2", "CLDN9", "CNTNAP2", "CLDN5", "CLDN6",
                   "PVR", "CDH1", "CDH3", "NRCAM", "CADM1", "L1CAM", "CLDN15", "CDH2", "CLDN1",
                   "CDH4", "CLDN16", "CLDN11", "JAM2")

figure5B <- pheatmap::pheatmap(abs(adhesion_cell)[sequence_gene,rownames(sample_annot)[order(sample_annot$SDA_DS, decreasing = F)]],
                               fontsize = 12, fontsize_row = 14, fontsize_col = 14,
                               color = c("#DA2C43", "#E15566", "#E97E88", "#F0A8AB", "#F8D1CD", "#FFFAF0"),
                               breaks = seq(0, 2, length.out = 7),
                               border_color = NA, cluster_cols = F, cluster_rows = F,
                               annotation_col = adhesion_annot,
                               annotation_colors = list(DS_pathway = c("#388230", "#EFED8F")),
                               clustering_method = "ward.D2",  main = "") 


## Figure 5CD for cell adhesion PathView
figure5C <- pathview(gene.data = gene_devscore$ILC[,"BCK4"], pathway.id = "04514",
                     species = "hsa", out.suffix = "KEGG_CELL_ADHESION_MOLECULES_CAMS (BCK4)", gene.idtype = "SYMBOL",
                     limit=list(gene=c(-3,3)), bins = list(gene = 16),res=400,
                     low = list(gene = "#B2E6F0"), mid = list(gene = "#FF0018"), high = list(gene = "#FFF485"), node.sum = "mean")



## Supplementary Figure 3 for cell adhesion violin plot
gene_adhesion <- intersect(deg, kegg_pathways$KEGG_CELL_ADHESION_MOLECULES_CAMS)
tumor_adhesion <- cbind.data.frame(label = tumor_label, train_norm$Xc[,gene_adhesion], stringsAsFactors = FALSE)
cell_adhesion  <- test_norm[row.names(sample_annot), gene_adhesion]
tumor_adhesion2 <- data.frame(melt(tumor_adhesion, id.vars = "label"))
cell_adhesion2 <- melt(cell_adhesion)
adhesion_genelist <- colnames(cell_adhesion)

adhesion_gene_violin_p1 <- ggplot() +
  geom_violin(data = tumor_adhesion2 %>% filter(label == "ILC" & variable %in% adhesion_genelist[1:11]),
              aes(x = variable, y = value),
              color = "#FFAEBC") +
  theme_bw() + 
  geom_point(position = position_jitter(seed = 1, width = 0, height = 0),
             data = cell_adhesion2[cell_adhesion2$Var2 %in% adhesion_genelist[1:11] & cell_adhesion2$Var1 %in% c("HCC2218", "BCK4"),],
             aes(y = value, x = Var2, fill = Var1), size = 3, shape=23, color = "white") +  
  scale_fill_manual(name = "Cell line", values = c("#3AC9B0","#F2C935")) +
  labs(x = "", y = "") +
  theme(legend.position = "none") + 
  coord_flip() +
  theme(plot.title = element_text(size = 20, face = "bold"), axis.title = element_text(size = 15), axis.text = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15))

adhesion_gene_violin_p2 <- ggplot() +
  geom_violin(data = tumor_adhesion2 %>% filter(label == "ILC" & variable %in% adhesion_genelist[12:22]),
              aes(x = variable, y = value),
              color = "#FFAEBC") +
  theme_bw() + 
  geom_point(position = position_jitter(seed = 1, width = 0, height = 0),
             data = cell_adhesion2[cell_adhesion2$Var2 %in% adhesion_genelist[12:22] & cell_adhesion2$Var1 %in% c("HCC2218", "BCK4"),],
             aes(y = value, x = Var2, fill = Var1), size = 3, shape=23, color = "white") +  
  scale_fill_manual(name = "Cell line", values = c("#3AC9B0","#F2C935")) +
  labs(x = "", y = "") +
  coord_flip() +
  theme(plot.title = element_text(size = 20, face = "bold"), axis.title = element_text(size = 15), axis.text = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 17))

adhesion_gene_dist <- plot_grid(adhesion_gene_violin_p1, adhesion_gene_violin_p2, ncol = 2, rel_widths = c(1, 1.3))
y.grob <- textGrob("Gene", 
                   gp=gpar(fontsize=17), rot=90)
x.grob <- textGrob("Expression", 
                   gp=gpar(fontsize=17))
supp_figure3 <- grid.arrange(arrangeGrob(adhesion_gene_dist, left = y.grob, bottom = x.grob),top = textGrob("Cell adhesion gene distribution", gp=gpar(fontsize=20)))


