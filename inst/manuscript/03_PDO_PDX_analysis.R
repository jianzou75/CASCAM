library(CASCAM)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(patchwork)
allcolour = c("#74AF74", "#5278CC", "#FFAEBC")
load("./data/raw_data/BRCA_PDO_PDX_CASCAM.RData")

## Figure 6A for PDO PDX position
brca_pdo_pdx <- genome_selection(brca_pdo_pdx)
subtype_levels <- levels(factor(brca_pdo_pdx@tumor_label))
y <- brca_pdo_pdx@tumor_sda_project; train_lab <- brca_pdo_pdx@tumor_label
center <- sapply(subtype_levels, function(t) median(brca_pdo_pdx@tumor_sda_project[brca_pdo_pdx@tumor_label == t]))
pool_sd <- mad(unlist(sapply(subtype_levels, function(t) brca_pdo_pdx@tumor_sda_project[brca_pdo_pdx@tumor_label == t] - median(brca_pdo_pdx@tumor_sda_project[brca_pdo_pdx@tumor_label == t]))))

distance_train <- sapply(1:2, function(r) abs(y[train_lab == subtype_levels[r]] - center[r])/pool_sd)
distribution_SDA <- sapply(distance_train,  ecdf)
SDA_DS_pval = sapply(1:2, function(r) 1-distribution_SDA[[r]](brca_pdo_pdx@sda_ds[,r]))
colnames(SDA_DS_pval) = levels(factor(train_lab))
rownames(SDA_DS_pval) = rownames(brca_pdo_pdx@sda_ds)
brca_pdo_pdx@sda_ds_pval = SDA_DS_pval

uninterested_subtype <- setdiff(brca_pdo_pdx@sda_model$classes, brca_pdo_pdx@interested_subtype)
position = data.frame(place = "pos", val = brca_pdo_pdx@camod_sda_project,
                      classification_predict_prob = ifelse((brca_pdo_pdx@sda_predict_prob[,brca_pdo_pdx@interested_subtype] > 0.5),
                                                           brca_pdo_pdx@interested_subtype, uninterested_subtype),
                      classification_sda_ds_pval = ifelse(brca_pdo_pdx@sda_ds_pval[,brca_pdo_pdx@interested_subtype] > 0.05,
                                                          brca_pdo_pdx@interested_subtype, uninterested_subtype),
                      classification_combine = ifelse((brca_pdo_pdx@sda_predict_prob[,brca_pdo_pdx@interested_subtype] > 0.5 & brca_pdo_pdx@sda_ds_pval[,brca_pdo_pdx@interested_subtype] > 0.05),
                                                      brca_pdo_pdx@interested_subtype, uninterested_subtype),
                      sda_predict_prob_interested = brca_pdo_pdx@sda_predict_prob[,brca_pdo_pdx@interested_subtype],
                      sda_predict_prob_uninterested = brca_pdo_pdx@sda_predict_prob[,uninterested_subtype],
                      sda_lds_interested = brca_pdo_pdx@sda_lds[,brca_pdo_pdx@interested_subtype],
                      sda_lds_uninterested = brca_pdo_pdx@sda_lds[,uninterested_subtype],
                      sda_ds_interested = brca_pdo_pdx@sda_ds[,brca_pdo_pdx@interested_subtype],
                      sda_ds_uninterested = brca_pdo_pdx@sda_ds[,uninterested_subtype])
position$camod = rownames(position)
position$camod = gsub("\\_|\\.", "-", position$camod)
position$camod = gsub("X", "", position$camod)
tumor_density = data.frame(val = brca_pdo_pdx@tumor_sda_project, grp = brca_pdo_pdx@tumor_label)
position$camod = ifelse(grepl("171881", position$camod), position$camod, NA)
position$camod = gsub("171881-019-R-", "", position$camod)
position$camod[position$camod == "V1-organoid"] = "PDO.1"
position$camod[position$camod == "APW-DS2"] = "PDX.0"
position$camod[position$camod == "APYF68"] = "PDX.1A"
position$camod[position$camod == "APWG05"] = "PDX.1B"
position$camod[position$camod == "APWG05PF7"] = "PDX.2A"
position$camod[position$camod == "APVG40-RG-G15"] = "PDX.2B"

subtype_levels <- levels(factor(brca_pdo_pdx@tumor_label))
center <- sapply(subtype_levels, function(t) median(brca_pdo_pdx@tumor_sda_project[brca_pdo_pdx@tumor_label == t]))
scatter = ggplot(data = position, aes(x = place, y = val)) +
  geom_jitter(aes(color = classification_predict_prob, shape = classification_sda_ds_pval), size = 2, position = position_jitter(seed = 1, width = 0.1)) +
  guides(size = FALSE) +
  scale_colour_manual(name = "SDA_class selection",  breaks = c(brca_pdo_pdx@interested_subtype, uninterested_subtype), values = c("#FF5334", "#4F9BFA")) +
  scale_shape_manual(name = "SDA_DS selection",  breaks = c(brca_pdo_pdx@interested_subtype, uninterested_subtype), values =  c(19, 2)) +
  geom_hline(aes(yintercept = center[uninterested_subtype]), color = "#4F9BFA") +
  geom_text(aes(0.6, center[uninterested_subtype], label = paste0(uninterested_subtype, " tumor center"), vjust = -1) , color = "#4F9BFA") +
  geom_hline(aes(yintercept = center[brca_pdo_pdx@interested_subtype]), color = "#FF5334") +
  geom_text(aes(0.6, center[brca_pdo_pdx@interested_subtype], label = paste0(brca_pdo_pdx@interested_subtype, " tumor center"), vjust = -1) , color = "#FF5334") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none") +
  ylab("SDA Projection") +
  ylim(c(min(c(position$val, tumor_density$val)), max(c(position$val, tumor_density$val)))) +
  geom_text_repel(aes(label= camod), alpha = 0.5, position = position_jitter(seed = 1, width = 0.1),
                  min.segment.length = unit(0, 'lines')) +
  ggtitle("")
dens =  ggplot(tumor_density, aes(x = val, color = grp)) +
  geom_density(size = 1.5) +
  xlim(c(min(c(position$val, tumor_density$val)), max(c(position$val, tumor_density$val)))) +
  theme_void() +
  theme(legend.position = "none") +
  scale_colour_manual(breaks = c(brca_pdo_pdx@interested_subtype, uninterested_subtype), values = c("#FF5334", "#4F9BFA")) +
  coord_flip()
figure6A <- scatter + dens + plot_layout(ncol = 2, nrow = 1, widths = c(4, 1), heights = c(1, 4))


## Figure 6B for PDO PDX pathway heatmap
geometric.mean <- function(x) {exp(mean(log(x)))}

subtype_levels <- levels(factor(brca_pdo_pdx@tumor_label))
tumor_expr <- t(brca_pdo_pdx@tumor_norm_data)
camod_expr <- t(brca_pdo_pdx@camod_norm_data)[,]
deg <- rownames(tumor_expr)

center_location <- sapply(subtype_levels, function(t) apply(tumor_expr[,brca_pdo_pdx@tumor_label == t], 1, function(x) median(x)))
pool_sd <- apply(tumor_expr, 1, function(y) mad(unlist(sapply(subtype_levels, function(t) y[brca_pdo_pdx@tumor_label == t] - median(y[brca_pdo_pdx@tumor_label == t])))))
gene_ds <- lapply(subtype_levels, function(t) apply(camod_expr, 2, function(x) (x - center_location[,t])/pool_sd)); names(gene_ds) = subtype_levels

pathways <- c(qusage::read.gmt(system.file(package = 'CASCAM', 'extdata/c2.cp.kegg.v7.4.symbols.gmt')),
              qusage::read.gmt(system.file(package = 'CASCAM', 'extdata/h.all.v7.4.symbols.gmt')))
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
colnames(pathway_ds_subtype1) <- colnames(pathway_ds_subtype2) <- colnames(camod_expr)
pathway_ds <- list(pathway_ds_subtype1, pathway_ds_subtype2)
names(pathway_ds) <- subtype_levels

brca_pdo_pdx@gene_ds <- gene_ds
brca_pdo_pdx@pathway_ds <- pathway_ds
brca_pdo_pdx@available_pathways <- path_row_names

fgseaRes <- brca_pdo_pdx@GSEA[match(brca_pdo_pdx@available_pathways, brca_pdo_pdx@GSEA$pathway), ]
path.row.names <- fgseaRes$pathway[abs(fgseaRes$NES) > 1.5 & fgseaRes$size < 200 & fgseaRes$size > 30]
path.row.names <- path.row.names[order(fgseaRes$NES[abs(fgseaRes$NES) > 1.5 & fgseaRes$size < 200 & fgseaRes$size > 30], decreasing = T)]
path.row.names <- c(path.row.names, "KEGG_CELL_ADHESION_MOLECULES_CAMS")

pathway_annotate <- data.frame(NES =  cut(as.numeric(fgseaRes$NES), seq(-2.5, 2.5,length.out = 5), labels = c(-2, -1, 1, 2)),
                               pathway_size = fgseaRes$size); row.names(pathway_annotate) = fgseaRes$pathway
pathway_annotate <- pathway_annotate[path.row.names,]

sample_annotate <- data.frame(SDA_DS = brca_pdo_pdx@sda_ds[colnames(brca_pdo_pdx@pathway_ds[[brca_pdo_pdx@interested_subtype]]), brca_pdo_pdx@interested_subtype])
row.names(sample_annotate) <- colnames(brca_pdo_pdx@pathway_ds[[brca_pdo_pdx@interested_subtype]])

heatmap_mat <- brca_pdo_pdx@pathway_ds[[brca_pdo_pdx@interested_subtype]][, grepl("171881", colnames(brca_pdo_pdx@pathway_ds[[brca_pdo_pdx@interested_subtype]]))]
heatmap_mat <- heatmap_mat[, c("X171881_019.R_APW.DS2", "X171881_019.R_APYF68",
                               "X171881_019.R_APWG05", "X171881_019.R_V1.organoid",
                               "X171881_019.R_APWG05PF7", "X171881_019.R_APVG40_RG.G15")]
heatmap_mat = heatmap_mat[rownames(pathway_annotate),]
heatmap_mat = rbind(heatmap_mat, combined_pathway = colMeans(heatmap_mat[1:14,]))

figure6B <- pheatmap::pheatmap(heatmap_mat,
                               color =  c("#d73027", "#fc8d59", "#fee090", "#ffffbf", "#e0f3f8", "#91bfdb", "#4575b4"),#c("#DA2C43", "#E15566", "#E97E88", "#F0A8AB", "#F8D1CD", "#FFFAF0"),
                               cluster_cols = F, cluster_rows = F,
                               breaks = seq(0, 2, length.out = 7),
                               annotation_row = rbind(pathway_annotate, combined_pathway = c(NA, NA)),
                               annotation_col = sample_annotate,
                               annotation_colors = list(NES = c("2" = "#0025F6", "1" ="grey", "-1" = "grey", "-2" ="#E60000"),
                                                        pathway_size = c("#FAFA31", "#FFB131"),
                                                        SDA_DS = c("#388230", "#EFED8F")),
                               labels_col = c("PDX.0", "PDX.1A",
                                              "PDX.1B", "PDO",
                                              "PDX.2A", "PDX.2B"),
                               clustering_method = "ward.D2", main = "")

label_level <- colnames(center_location)
gene_devscore_tumor <- lapply(label_level, function(t) apply(t(brca_pdo_pdx@tumor_norm_data)[,brca_pdo_pdx@tumor_label == "ILC"], 2,
                                                             function(x) (x - center_location[,t])/pool_sd)); names(gene_devscore_tumor) = label_level
geometric.mean <- function(x) {exp(mean(log(x)))}
pathway_pds_ILC_tumor <- c();                       
path.row.names <- c()

kegg_pathways <- c(qusage::read.gmt(system.file(package = 'CASCAM', 'extdata/c2.cp.kegg.v7.4.symbols.gmt')),
                   qusage::read.gmt(system.file(package = 'CASCAM', 'extdata/h.all.v7.4.symbols.gmt')))
deg <- rownames(gene_devscore_tumor$ILC)
for (i in 1:length(kegg_pathways)){
  if (sum(deg %in% kegg_pathways[[i]]) > 20){
    gds_ILC <- gene_devscore_tumor$ILC[intersect(deg, kegg_pathways[[i]]), ]
    pathway_pds_ILC_tumor <- rbind(pathway_pds_ILC_tumor, c(apply(abs(gds_ILC), 2, geometric.mean)))
    
    path.row.names<- c(path.row.names, names(kegg_pathways)[i])
  }
}
rownames(pathway_pds_ILC_tumor)  <- path.row.names
colnames(pathway_pds_ILC_tumor)  <- colnames(brca_pdo_pdx@tumor_aligned_data)[brca_pdo_pdx@tumor_label == "ILC"]

pathway_pds_ILC_tumor <- pathway_pds_ILC_tumor[rownames(heatmap_mat)[1:15], ]
pathway_pds_ILC_tumor = rbind(pathway_pds_ILC_tumor, combined_pathway = colMeans(pathway_pds_ILC_tumor[1:14,]))
tumor_path_ds_ecdf <- apply(pathway_pds_ILC_tumor, 1, ecdf)
pathway_heatmap_pvalue <- sapply(1:16, function(i) tumor_path_ds_ecdf[[i]](heatmap_mat[i,]))
pathway_heatmap_pvalue <- t(1 - pathway_heatmap_pvalue)
rownames(pathway_heatmap_pvalue) <- rownames(heatmap_mat) ## p-values for the heatmap
colnames(pathway_heatmap_pvalue) <- colnames(heatmap_mat)


## Figure 6C for cell adhesion pathway violin plot
allcolour = c("#3AC9B0","#F2C935","#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
              "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
              "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
              "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
pathway_name = "KEGG_CELL_ADHESION_MOLECULES_CAMS"
interested_camods = c("X171881_019.R_V1.organoid", "X171881_019.R_APWG05")

pathways <- c(qusage::read.gmt(system.file(package = 'CASCAM', 'extdata/c2.cp.kegg.v7.4.symbols.gmt')),
              qusage::read.gmt(system.file(package = 'CASCAM', 'extdata/h.all.v7.4.symbols.gmt')))
gene_pathway <- intersect(c(na.omit(brca_pdo_pdx@DEA$delabel)), # rownames(brca_pdo_pdx@gene_ds[[brca_pdo_pdx@interested_subtype]]),
                          pathways[[pathway_name]])
gene_pathway <- gene_pathway[order(brca_pdo_pdx@DEA[gene_pathway, 1])]
uninterested_subtype <- setdiff(levels(factor(brca_pdo_pdx@tumor_label)), brca_pdo_pdx@interested_subtype)

tumor_pathway <- cbind.data.frame(label = brca_pdo_pdx@tumor_label, brca_pdo_pdx@tumor_norm_data[,gene_pathway], stringsAsFactors = FALSE)
camod_pathway  <- brca_pdo_pdx@camod_norm_data[, gene_pathway]
tumor_pathway2 <- data.frame(reshape2::melt(tumor_pathway, id.vars = "label"))
camod_pathway2 <- reshape2::melt(as.matrix(camod_pathway))
tumor_pathway2 <- tumor_pathway2 %>% filter(variable %in% c("NRXN2", "CDH15", "L1CAM", "CADM1", "CADM3", "CDH2"))
camod_pathway2 <- camod_pathway2 %>% filter(Var2 %in% c("NRXN2", "CDH15", "L1CAM", "CADM1", "CADM3", "CDH2"))

figure6C <- ggplot() +
  geom_violin(data =  tumor_pathway2 %>% filter(label == "ILC"), aes(x = value, y = variable, color = label), scale = 1,  alpha = 0, size = 1) +
  geom_violin(data =  tumor_pathway2 %>% filter(label == "IDC"), aes(x = value, y = variable, color = label), scale = 1,  alpha = 0, size = 1) +
  scale_colour_manual(name = "tumor distribution", values = c("IDC" = "#BBE7FE", "ILC" = "#FFAEBC")) +
  theme_bw() +
  geom_point(position = position_jitter(seed = 1, width = 0, height = 0), data = camod_pathway2[camod_pathway2$Var1 %in% interested_camods,], aes(x = value, y = Var2, fill = Var1), size = 3, shape=21, color = "white") +
  scale_fill_manual(name = "cancer model", values = c("#9E8964", "#76AB97"), labels = c("PDO.1", "PDX.1B")) +
  labs(x = "Expression", y = "Gene") +
  xlim(c(-0.1, 0.2)) +
  ggtitle("")


## Supplementary figure 5 for complete cell adhesion pathway violin plot
pathway_name = "KEGG_CELL_ADHESION_MOLECULES_CAMS"
interested_camods = c("X171881_019.R_V1.organoid", "X171881_019.R_APWG05")

pathways <- c(qusage::read.gmt(system.file(package = 'CASCAM', 'extdata/c2.cp.kegg.v7.4.symbols.gmt')),
              qusage::read.gmt(system.file(package = 'CASCAM', 'extdata/h.all.v7.4.symbols.gmt')))
gene_pathway <- intersect(c(na.omit(brca_pdo_pdx@DEA$delabel)), # rownames(brca_pdo_pdx@gene_ds[[brca_pdo_pdx@interested_subtype]]),
                          pathways[[pathway_name]])
gene_pathway <- gene_pathway[order(brca_pdo_pdx@DEA[gene_pathway, 1])]
uninterested_subtype <- setdiff(levels(factor(brca_pdo_pdx@tumor_label)), brca_pdo_pdx@interested_subtype)

tumor_pathway <- cbind.data.frame(label = brca_pdo_pdx@tumor_label, brca_pdo_pdx@tumor_norm_data[,gene_pathway], stringsAsFactors = FALSE)
camod_pathway  <- brca_pdo_pdx@camod_norm_data[, gene_pathway]
tumor_pathway2 <- data.frame(reshape2::melt(tumor_pathway, id.vars = "label"))
camod_pathway2 <- reshape2::melt(as.matrix(camod_pathway))
tumor_pathway2 <- tumor_pathway2 # %>% filter(variable %in% c("NRXN2", "CDH15", "L1CAM", "CADM1", "CADM3", "CDH2"))
camod_pathway2 <- camod_pathway2 # %>% filter(Var2 %in% c("NRXN2", "CDH15", "L1CAM", "CADM1", "CADM3", "CDH2"))

ecaherin_violin_p1 <- ggplot() +
  geom_violin(data =  tumor_pathway2 %>% filter(label == "ILC" & variable %in% gene_pathway[1:11]), aes(x = value, y = variable, color = label), scale = 1,  alpha = 0, size = 1) +
  geom_violin(data =  tumor_pathway2 %>% filter(label == "IDC" & variable %in% gene_pathway[1:11]), aes(x = value, y = variable, color = label), scale = 1,  alpha = 0, size = 1) +
  scale_colour_manual(name = "tumor distribution", values = c("IDC" = "#BBE7FE", "ILC" = "#FFAEBC")) +
  theme_bw() +
  geom_point(position = position_jitter(seed = 1, width = 0, height = 0), data = camod_pathway2[camod_pathway2$Var1 %in% interested_camods & camod_pathway2$Var2 %in% gene_pathway[1:11],], aes(x = value, y = Var2, fill = Var1), size = 3, shape=21, color = "white") +
  scale_fill_manual(name = "cancer model", values = c("#9E8964", "#76AB97"), labels = c("PDO.1", "PDX.1B")) +
  labs(x = "Expression", y = "Gene") +
  xlim(c(-0.1, 0.2)) +
  labs(x = "", y = "") +
  theme(legend.position = "none") + 
  ggtitle("") +
  theme(plot.title = element_text(size = 20, face = "bold"), axis.title = element_text(size = 15), axis.text = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15))

ecaherin_violin_p2 <- ggplot() +
  geom_violin(data =  tumor_pathway2 %>% filter(label == "ILC" & variable %in% gene_pathway[12:22]), aes(x = value, y = variable, color = label), scale = 1,  alpha = 0, size = 1) +
  geom_violin(data =  tumor_pathway2 %>% filter(label == "IDC" & variable %in% gene_pathway[12:22]), aes(x = value, y = variable, color = label), scale = 1,  alpha = 0, size = 1) +
  scale_colour_manual(name = "tumor distribution", values = c("IDC" = "#BBE7FE", "ILC" = "#FFAEBC")) +
  theme_bw() +
  geom_point(position = position_jitter(seed = 1, width = 0, height = 0), data = camod_pathway2[camod_pathway2$Var1 %in% interested_camods & camod_pathway2$Var2 %in% gene_pathway[12:22],], aes(x = value, y = Var2, fill = Var1), size = 3, shape=21, color = "white") +
  scale_fill_manual(name = "cancer model", values = c("#9E8964", "#76AB97"), labels = c("PDO.1", "PDX.1B")) +
  labs(x = "Expression", y = "Gene") +
  xlim(c(-0.1, 0.2)) +
  labs(x = "", y = "") +
  ggtitle("") +
  theme(plot.title = element_text(size = 20, face = "bold"), axis.title = element_text(size = 15), axis.text = element_text(size = 15), legend.text = element_text(size = 15), legend.title = element_text(size = 15))

adhesion_gene_dist <- plot_grid(ecaherin_violin_p1, ecaherin_violin_p2, ncol = 2, rel_widths = c(1, 1.3))
y.grob <- textGrob("Gene", 
                   gp=gpar(fontsize=17), rot=90)
x.grob <- textGrob("Expression", 
                   gp=gpar(fontsize=17))
supp_figure5 <- grid.arrange(arrangeGrob(adhesion_gene_dist, left = y.grob, bottom = x.grob))



