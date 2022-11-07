library(aricode);library(glmnet); library(caret); library(xtable); library(ggrepel); library(tidyverse); library(ggpubr); library(patchwork);
library(boot); library(cowplot); library(latex2exp)
tumor_aligned <- CASCAM_ILC@tumor_aligned_data 
tumor_label <- CASCAM_ILC@tumor_label
cell_aligned <- CASCAM_ILC@camod_aligned_data


genome_analysis <- function(train_dat, train_lab, test_dat){
  subtype_levels = levels(factor(train_lab))
  
  train_norm = sparseLDA::normalize(t(train_dat))
  test_norm  = sparseLDA::normalizetest(t(test_dat), train_norm)
  train_lab_mat = mclust::unmap(train_lab)
  colnames(train_lab_mat) = subtype_levels
  splda_model <- sparseLDA::sda(train_norm$Xc, train_lab_mat, stop = -153, lambda = 1e-04)
  
  x <- test_norm[,splda_model$varNames] %*% splda_model$beta	
  y <- train_norm$Xc[,splda_model$varNames] %*% splda_model$beta	
  
  center <- sapply(subtype_levels, function(t) median(y[train_lab == t]))
  pool_sd <- mad(unlist(sapply(subtype_levels, function(t) y[train_lab == t] - median(y[train_lab == t]))))
  
  ## SDA_DS & SDA_LDS
  dist <- sapply(1:nrow(splda_model$fit$means), function(r) abs(x - center[r])/pool_sd)
  colnames(dist) = levels(factor(train_lab))
  dist = data.frame(dist)
  rownames(dist) = rownames(x)
  SDA_DS = dist
  SDA_LDS = log2(SDA_DS)
  
  ## SDA_P
  SDA_P <- predict(splda_model, test_norm, prior = c(0.5, 0.5))$posterior
  
  ## SDA_DS_pval
  distance_train <- sapply(1:nrow(splda_model$fit$means), function(r) abs(y[train_lab == subtype_levels[r]] - center[r])/pool_sd)
  distribution_SDA <- sapply(distance_train,  ecdf)
  SDA_DS_pval = sapply(1:nrow(splda_model$fit$means), function(r) 1-distribution_SDA[[r]](dist[,r]))
  colnames(SDA_DS_pval) = levels(factor(train_lab))
  rownames(SDA_DS_pval) = rownames(dist)
  
  
  return(list(splda_model = splda_model,
              SDA_DS = SDA_DS,
              SDA_P = SDA_P,
              SDA_LDS = SDA_LDS,
              SDA_center = center,
              project_test = x,
              project_train = y,
              SDA_DS_pval = SDA_DS_pval))
}

bootstrap_ci <- function(train_dat, train_lab, test_dat){ ## bootstrapping with 1000 replications
  subtype_levels = levels(factor(train_lab))
  train_norm = sparseLDA::normalize(t(train_dat))
  test_norm  = sparseLDA::normalizetest(t(test_dat), train_norm)
  train_lab_mat = mclust::unmap(train_lab)
  colnames(train_lab_mat) = subtype_levels
  model <- sparseLDA::sda(train_norm$Xc, train_lab_mat, stop = -153, lambda = 1e-04)
  
  boot.stat.sda <- function(data, i){
    
    x <- test_norm[,model$varNames] %*% model$beta	
    y <- data[,model$varNames] %*% model$beta	
    
    center <- sapply(subtype_levels, function(t) median(y[i][train_lab[i] == t]))
    pool_sd <- mad(unlist(sapply(subtype_levels, function(t) y[i][train_lab[i] == t] - median(y[i][train_lab[i] == t]))))
    logILC <- log2(abs(x - center[2])/pool_sd)
    
    
    return(c(logILC))
  }
  
  results <- boot(data=train_norm$Xc, 
                  statistic = boot.stat.sda,
                  R = 1000)
  correct.ci <- sapply(1:nrow(test_norm), function(s) c(boot.ci(results, type = "norm", index = s)$normal[2:3]) + c(colMeans(results$t)[s] - results$t0[s]))
  logILC_ci = t(correct.ci); colnames(logILC_ci) = c("logILC_lower", "logILC_upper")
  
  return(list(LDS_DSA_ILC_ci = logILC_ci))
}


## Figure 4A for projected values
set.seed(123)
genome <- genome_analysis(tumor_aligned, tumor_label, cell_aligned)
position = data.frame(place = "pos", val = genome$project_test,
                      classification_P = ifelse((genome$SDA_P[,2] > 0.5), "ILC", "IDC"),
                      SDA_LDS_threshold = ifelse(genome$SDA_DS_pval[,2] > 0.05 , "ILC", "No ILC"),
                      classification_P_LDS = ifelse((genome$SDA_P[,2] > 0.5 & genome$SDA_DS_pval[,2] > 0.05 ), "ILC", "No ILC"),
                      SDA_ILC_P = genome$SDA_P[,2],
                      SDA_IDC_P = genome$SDA_P[,1],
                      SDA_ILC_LDS = genome$SDA_LDS[,2],
                      SDA_IDC_LDS = genome$SDA_LDS[,1],
                      SDA_ILC_DS  = genome$SDA_DS[,2],
                      SDA_IDC_DS  = genome$SDA_DS[,1])
position$cell = rownames(position)
tumor.density = data.frame(val = genome$project_train,
                           grp = tumor_label)
position$cell = ifelse(position$cell %in% c("SUM44PE", "UACC812", "ZR751", "DU4475", "MDAMB453", "ZR7530","MDAMB134VI", "OCUBM", "UACC893",
                                            rownames(position)[position$classification_P_LDS == "ILC"]), position$cell, NA)

p1 = ggplot(data = position, aes(x = place, y = val)) +
  geom_jitter(aes(color = classification_P, shape = SDA_LDS_threshold), size = 2, position = position_jitter(seed = 1, width = 0.2)) +
  guides(size = FALSE) +
  scale_colour_manual(name = "P_SDA selection",  labels = c("IDC", "ILC"), values = c("#4F9BFA", "#FF5334")) +
  scale_shape_manual(name = "DS_SDA selection",  labels = c("ILC", "No ILC"), values =  c(19, 2)) +
  geom_hline(aes(yintercept = genome$SDA_center["IDC"]), color = "#4F9BFA") +
  geom_text(aes(0.6, genome$SDA_center["IDC"], label = "IDC tumor center", vjust = -1) , color = "#4F9BFA") +
  geom_hline(aes(yintercept = genome$SDA_center["ILC"]), color = "#FF5334") +
  geom_text(aes(0.6, genome$SDA_center["ILC"], label = "ILC tumor center", vjust = -1) , color = "#FF5334") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none") +               
  ylab("SDA Projection") +
  ylim(c(-1.8, 2.1)) +
  geom_text_repel(aes(label= cell), alpha = 0.5, position = position_jitter(seed = 1, width = 0.2),
                  min.segment.length = unit(0, 'lines')) +
  ggtitle("") +
  theme(legend.position="none", plot.title = element_text(size = 15, face = "bold"), axis.title = element_text(size = 15), axis.text = element_text(size = 15))

dens =  ggplot(tumor.density, aes(x = val, color = grp)) + 
  geom_density(size = 1.5) + 
  xlim(c(-1.8, 2.1)) +
  theme_void() + 
  theme(legend.position = "none") + 
  scale_colour_manual(values = c("IDC" = "#4F9BFA", "ILC" = "#FF5334")) +
  coord_flip()
figure4A = p1 + dens + plot_layout(ncol = 2, nrow = 1, widths = c(2, 1), heights = c(1, 4))


## Figure 4B for confidence interval
position = data.frame(place = "pos", val = genome$project_test,
                      classification_P = ifelse((genome$SDA_P[,2] > 0.8), "ILC", "IDC"),
                      SDA_LDS_threshold = ifelse(genome$SDA_DS_pval[,2] > 0.1, "ILC", "No ILC"),
                      classification_P_LDS = ifelse((genome$SDA_P[,2] > 0.8 & genome$SDA_DS_pval[,2] > 0.1), "ILC", "No ILC"),
                      SDA_ILC_P = genome$SDA_P[,2],
                      SDA_IDC_P = genome$SDA_P[,1],
                      SDA_ILC_LDS = genome$SDA_LDS[,2],
                      SDA_IDC_LDS = genome$SDA_LDS[,1],
                      SDA_ILC_DS  = genome$SDA_DS[,2],
                      SDA_IDC_DS  = genome$SDA_DS[,1])
position$cell = rownames(position)
genome.ci <- bootstrap_ci(tumor_aligned, tumor_label, cell_aligned)
genome.ci.table <- cbind(position, genome.ci$LDS_DSA_ILC_ci)

figure4B <- ggplot(data=genome.ci.table[position$classification_P_LDS == "ILC"|genome.ci.table$cell == "MDAMB134VI", ] %>%
                     mutate(ILC_name = factor(cell, levels = cell[order(SDA_ILC_LDS, decreasing = T)])),
                   aes(y=ILC_name, x=SDA_ILC_LDS, xmin=logILC_lower, xmax=logILC_upper))+ 
  geom_point()+ 
  geom_errorbarh(height=.1) +
  theme_bw() +
  ylab("") +
  xlab("") +
  labs(title = "") + 
  scale_x_continuous(breaks = c(-7.5, -5, -2.5, 0, 1, log2(3)), 
                     labels = round(2^c(-7.5, -5, -2.5, 0, 1, log2(3)), 2)) +
  theme(legend.position="none", plot.title = element_text(size = 20, face = "bold"), axis.title = element_text(size = 15), axis.text = element_text(size = 15))

