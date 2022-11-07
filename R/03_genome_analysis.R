##################################################################################################################################
# 03: Genome-wide pre-selection analysis
# The step of genome-wide analysis, and all the cancer models are fed into the framework for pre-selection.

#' SDA model trained by tumor gene expression
#'
#' This function trains a sparse discriminant analysis model using tumor gene expression data. The parameters can be self-determined or
#' determined by the cross-validation. If cross-validation is prefered, users can refer to the function \link{sda_model_cv}.
#'
#' @param object A \code{CASCAM} object.
#' @param stop A negative value with its absolute value representing the number of genes selected for model training.
#' @param lambda A parameter for the L2-norm for elastic net regression.
#'
#' @import sparseLDA
#' @importFrom mclust unmap
#'
#' @return A \code{CASCAM} object with sda_model, camod_norm_data, tumor_norm_data, tumor_sda_project, and camod_sda_project slots.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' gene_info <- create_InformativeGenes(brca_ct, brca_label2, "ILC")
#' CASCAM_eg <- create_CASCAM(tumor_aligned, tumor_label, cell_aligned, gene_info)
#' CASCAM_eg <- sda_model(CASCAM_eg)
#' }
#'
#' @references
#' Clemmensen, Line, Trevor Hastie, Daniela Witten, and Bjarne Ersbøll. 2011. “Sparse Discriminant Analysis.”
#'  Technometrics: A Journal of Statistics for the Physical, Chemical, and Engineering Sciences 53 (4): 406–13. https://doi.org/10.1198/tech.2011.08118.
sda_model <- function(object, stop = NULL, lambda = NULL){
  set.seed(123)
  if(is.null(stop)| is.null(lambda)){
    stop <- -round(nrow(object@tumor_aligned_data)/20)
    lambda <- 1e-4
  }

  subtype_levels <- levels(factor(object@tumor_label))
  tumor_norm <- sparseLDA::normalize(t(object@tumor_aligned_data))
  camod_norm <- sparseLDA::normalizetest(t(object@camod_aligned_data), tumor_norm)
  tumor_label_mat <- mclust::unmap(object@tumor_label)
  colnames(tumor_label_mat) = subtype_levels

  sda_model_trained <- sparseLDA::sda(tumor_norm$Xc, tumor_label_mat, stop = stop, lambda = lambda)
  object@sda_model <- sda_model_trained
  object@tumor_norm_data <- data.frame(tumor_norm$Xc, check.names = F)
  object@camod_norm_data <- data.frame(camod_norm, check.names = F)
  object@tumor_sda_project <- as.matrix(object@tumor_norm_data[, object@sda_model$varNames]) %*% object@sda_model$beta
  object@camod_sda_project <- as.matrix(object@camod_norm_data[, object@sda_model$varNames]) %*% object@sda_model$beta

  return(object)
}


#' SDA model trained by tumor gene expression
#'
#' This function trains a sparse discriminant analysis model using tumor gene expression data with cross-validation for parameter selection.
#'  Users can refer to the function \link{sda_model} if parameters are known.
#'
#' @param object A \code{CASCAM} object.
#' @param stop_vector A vector of negative values with its absolute value representing the number of genes selected for model training.
#' @param lambda_vector A vector of values for the L2-norm for elastic net regression.
#' @param parallel_cores Number of cores for the parallel computation.
#'
#' @import caret
#' @import doParallel
#' @import sparseLDA
#'
#' @return A \code{CASCAM} object with sda_model, camod_norm_data, tumor_norm_data, tumor_sda_project, and camod_sda_project slots.
#' @export
#'
#' @examples
#' \dontrun{
#' gene_info <- create_InformativeGenes(tumor_ct, tumor_label2, "ILC")
#' CASCAM_eg <- create_CASCAM(tumor_aligned, tumor_label, camod_aligned, gene_info)
#' CASCAM_eg <- sda_model_cv(CASCAM_eg)
#' }
#'
#' @references
#' Clemmensen, Line, Trevor Hastie, Daniela Witten, and Bjarne Ersbøll. 2011. “Sparse Discriminant Analysis.”
#'  Technometrics: A Journal of Statistics for the Physical, Chemical, and Engineering Sciences 53 (4): 406–13. https://doi.org/10.1198/tech.2011.08118.
sda_model_cv <- function(object, stop_vector = NULL, lambda_vector = NULL, parallel_cores = 1){
  if(is.null(stop)| is.null(lambda)){
    searchgrid <- expand.grid(NumVars = round(c(nrow(object@tumor_aligned_data)/20, nrow(object@tumor_aligned_data)/10,
                                                nrow(object@tumor_aligned_data)/5,  nrow(object@tumor_aligned_data)/2)),
                              lambda = c(0, 0.1, 1e-4))
  } else{
    searchgrid <- expand.grid(NumVars = -round(stop_vector),
                              lambda = lambda_vector)
  }

  tumor_norm <- sparseLDA::normalize(t(object@tumor_aligned_data))
  camod_norm <- sparseLDA::normalizetest(t(object@camod_aligned_data), tumor_norm)

  cl <- makePSOCKcluster(parallel_cores)
  registerDoParallel(cl)
  tune <- caret::train(tumor_norm$Xc, train_lab, method = 'sparseLDA', tuneGrid = searchgrid)
  stopCluster(cl)

  object@sda_model <- tune$finalModel
  object@tumor_norm_data <- data.frame(tumor_norm$Xc, check.names = F)
  object@camod_norm_data <- data.frame(camod_norm, check.names = F)
  object@tumor_sda_project <- as.matrix(object@tumor_norm_data[, object@sda_model$varNames]) %*% object@sda_model$beta
  object@camod_sda_project <- as.matrix(object@camod_norm_data[, object@sda_model$varNames]) %*% object@sda_model$beta

  return(object)
}


#' Genome-wide pre-selection analysis
#'
#' @param object A \code{CASCAM} object.
#' @param R The number of bootstrap replicates. R = 1000 by default.
#' @param assignment_prob_cutoff The cutoff point for the mimimum SDA based assignment probability.
#' @param min_sda_ds_pval The cutoff point for the minimum SDA based deviance score p-value.
#' @param max_sda_ds_pval The cutoff point for the maximum SDA based deviance score p-value.
#'
#'
#' @import boot
#'
#' @return A \code{CASCAM} object with genome-wide pre-selection related slots.
#' @export
#'
#' @examples
#' \dontrun{
#' gene_info <- create_InformativeGenes(tumor_ct, tumor_label2, "ILC")
#' CASCAM_eg <- create_CASCAM(tumor_aligned, tumor_label, camod_aligned, gene_info)
#' CASCAM_eg <- sda_model(CASCAM_eg)
#' CASCAM_eg <- genome_selection(CASCAM_eg)
#' }
genome_selection <- function(object, R = 1000, assignment_prob_cutoff = 0.8, min_sda_ds_pval = 0.05){
  object@genome_selection_criteria <- c(assignment_prob_cutoff = assignment_prob_cutoff,
                                        min_sda_ds_pval = min_sda_ds_pval)

  subtype_levels <- levels(factor(object@tumor_label))
  center <- sapply(subtype_levels, function(t) median(object@tumor_sda_project[object@tumor_label == t]))
  pool_sd <- mad(unlist(sapply(subtype_levels, function(t) object@tumor_sda_project[object@tumor_label == t] - median(object@tumor_sda_project[object@tumor_label == t]))))

  ## SDA deviance score (sda_ds) and SDA log deviance score (sda_lds)
  dist <- sapply(1:length(subtype_levels), function(r) abs(object@camod_sda_project - center[r])/pool_sd)
  colnames(dist) <- subtype_levels
  dist <- data.frame(dist)
  rownames(dist) <- rownames(object@camod_sda_project)
  object@sda_ds <- dist
  object@sda_lds <- log2(object@sda_ds)

  ## SDA prediction probability (sda_predict_prob)
  object@sda_predict_prob <- sparseLDA::predict.sda(object@sda_model, object@camod_norm_data, prior = c(0.5, 0.5))$posterior

  ## p-value for SDA deviance scores (sda_ds_pval)
  distance_train <- sapply(1:length(subtype_levels), function(r) abs(object@tumor_sda_project[object@tumor_label == subtype_levels[r]] - center[r])/pool_sd)
  distribution_SDA <- sapply(distance_train,  ecdf)
  SDA_DS_pval = sapply(1:length(subtype_levels), function(r) 1-distribution_SDA[[r]](dist[,r]))
  colnames(SDA_DS_pval) = levels(factor(object@tumor_label))
  rownames(SDA_DS_pval) = rownames(dist)
  object@sda_ds_pval <- SDA_DS_pval

  ## bootstrap for confidence interval
  boot_stat_sda_cat1 <- function(data, i){
    center <- sapply(subtype_levels, function(t) median(data[i][object@tumor_label[i] == t]))
    pool_sd <- mad(unlist(sapply(subtype_levels, function(t) data[i][object@tumor_label[i] == t] - median(data[i][object@tumor_label[i] == t]))))
    log_cat1 <- log2(abs(object@camod_sda_project - center[1])/pool_sd)
    return(c(log_cat1))
  }
  results_cat1 <- boot::boot(data=object@tumor_sda_project, statistic =  boot_stat_sda_cat1, R = R)
  correct_ci_cat1 <- sapply(1:nrow(object@camod_norm_data), function(s) c(boot::boot.ci(results_cat1, type = "norm", index = s)$normal[2:3]) + c(colMeans(results_cat1$t)[s] - results_cat1$t0[s]))
  sda_lds_ci_cat1 <- t(correct_ci_cat1); colnames(sda_lds_ci_cat1) = c("lower", "upper"); rownames(sda_lds_ci_cat1) = rownames(object@sda_lds)
  boot_stat_sda_cat2 <- function(data, i){
    center <- sapply(subtype_levels, function(t) median(data[i][object@tumor_label[i] == t]))
    pool_sd <- mad(unlist(sapply(subtype_levels, function(t) data[i][object@tumor_label[i] == t] - median(data[i][object@tumor_label[i] == t]))))
    log_cat2 <- log2(abs(object@camod_sda_project - center[2])/pool_sd)
    return(c(log_cat2))
  }
  results_cat2 <- boot::boot(data=object@tumor_sda_project, statistic =  boot_stat_sda_cat2, R = R)
  correct_ci_cat2 <- sapply(1:nrow(object@camod_norm_data), function(s) c(boot::boot.ci(results_cat2, type = "norm", index = s)$normal[2:3]) + c(colMeans(results_cat2$t)[s] - results_cat2$t0[s]))
  sda_lds_ci_cat2 <- t(correct_ci_cat2); colnames(sda_lds_ci_cat2) = c("lower", "upper"); rownames(sda_lds_ci_cat2) = rownames(object@sda_lds)
  sda_lds_ci <- list(sda_lds_ci_cat1, sda_lds_ci_cat2)
  names(sda_lds_ci) <- names(center)
  object@sda_lds_ci <- sda_lds_ci

  uninterested_subtype <- setdiff(object@sda_model$classes, object@interested_subtype)
  classification_combine = ifelse((object@sda_predict_prob[,object@interested_subtype] > assignment_prob_cutoff & object@sda_ds_pval[,object@interested_subtype] > min_sda_ds_pval),
                                   object@interested_subtype, uninterested_subtype)
  object@selected_camods <- rownames(object@sda_predict_prob)[classification_combine == object@interested_subtype]

  return(object)
}


#' Visualization of genome-wide pre-selection analysis
#'
#' @param object
#'
#' @import magrittr
#' @import ggplot2
#' @import ggrepel
#' @import patchwork
#' @importFrom dplyr filter arrange
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' gene_info <- create_InformativeGenes(tumor_ct, tumor_label2, "ILC")
#' CASCAM_eg <- create_CASCAM(tumor_aligned, tumor_label, camod_aligned, gene_info)
#' CASCAM_eg <- sda_model(CASCAM_eg)
#' CASCAM_eg <- genome_selection(CASCAM_eg)
#' CASCAM_eg <- genome_selection_visualize(CASCAM_eg)
#' }
genome_selection_visualize <- function(object){
  assignment_prob_cutoff = object@genome_selection_criteria["assignment_prob_cutoff"]
  min_sda_ds_pval = object@genome_selection_criteria["min_sda_ds_pval"]

  uninterested_subtype <- setdiff(object@sda_model$classes, object@interested_subtype)
  position = data.frame(place = "pos", val = object@camod_sda_project,
                        classification_predict_prob = ifelse((object@sda_predict_prob[,object@interested_subtype] > assignment_prob_cutoff),
                                                              object@interested_subtype, uninterested_subtype),
                        classification_sda_ds_pval = ifelse(object@sda_ds_pval[,object@interested_subtype] > min_sda_ds_pval,
                                                            object@interested_subtype, uninterested_subtype),
                        classification_combine = ifelse((object@sda_predict_prob[,object@interested_subtype] > assignment_prob_cutoff & object@sda_ds_pval[,object@interested_subtype] > min_sda_ds_pval),
                                                         object@interested_subtype, uninterested_subtype),
                        sda_predict_prob_interested = object@sda_predict_prob[,object@interested_subtype],
                        sda_predict_prob_uninterested = object@sda_predict_prob[,uninterested_subtype],
                        sda_lds_interested = object@sda_lds[,object@interested_subtype],
                        sda_lds_uninterested = object@sda_lds[,uninterested_subtype],
                        sda_ds_interested = object@sda_ds[,object@interested_subtype],
                        sda_ds_uninterested = object@sda_ds[,uninterested_subtype])
  position$camod = rownames(position)
  tumor_density = data.frame(val = object@tumor_sda_project, grp = object@tumor_label)
  top_camod = position %>% dplyr::filter(classification_combine == object@interested_subtype) %>% dplyr::arrange(sda_ds_interested)
  if(nrow(top_camod) > 3){
    top_camod = top_camod$camod[1:3]
  } else{
    top_camod = top_camod$camod}
  position$camod = ifelse(position$camod %in% top_camod, position$camod, NA)

  ## Distribution figure (sda_project_position)
  subtype_levels <- levels(factor(object@tumor_label))
  center <- sapply(subtype_levels, function(t) median(object@tumor_sda_project[object@tumor_label == t]))
  scatter = ggplot(data = position, aes(x = place, y = val)) +
    geom_jitter(aes(color = classification_predict_prob, shape = classification_sda_ds_pval), size = 2, position = position_jitter(seed = 1, width = 0.1)) +
    guides(size = FALSE) +
    scale_colour_manual(name = "SDA_class selection",  breaks = c(object@interested_subtype, uninterested_subtype), values = c("#FF5334", "#4F9BFA")) +
    scale_shape_manual(name = "SDA_DS selection",  breaks = c(object@interested_subtype, uninterested_subtype), values =  c(19, 2)) +
    geom_hline(aes(yintercept = center[uninterested_subtype]), color = "#4F9BFA") +
    geom_text(aes(0.6, center[uninterested_subtype], label = paste0(uninterested_subtype, " tumor center"), vjust = -1) , color = "#4F9BFA") +
    geom_hline(aes(yintercept = center[object@interested_subtype]), color = "#FF5334") +
    geom_text(aes(0.6, center[object@interested_subtype], label = paste0(object@interested_subtype, " tumor center"), vjust = -1) , color = "#FF5334") +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="none") +
    ylab("SDA Projection") +
    ylim(c(min(c(position$val, tumor_density$val)), max(c(position$val, tumor_density$val)))) +
    geom_text_repel(aes(label= camod), alpha = 0.5, position = position_jitter(seed = 1, width = 0.1),
                    min.segment.length = unit(0, 'lines')) +
    ggtitle("Position of the cancer models")
  dens =  ggplot(tumor_density, aes(x = val, color = grp)) +
    geom_density(size = 1.5) +
    xlim(c(min(c(position$val, tumor_density$val)), max(c(position$val, tumor_density$val)))) +
    theme_void() +
    theme(legend.position = "none") +
    scale_colour_manual(breaks = c(object@interested_subtype, uninterested_subtype), values = c("#FF5334", "#4F9BFA")) +
    coord_flip()
  sda_project_position = scatter + dens + plot_layout(ncol = 2, nrow = 1, widths = c(4, 1), heights = c(1, 4))

  ## Confidence interval
  genome_ci <- cbind(position, object@sda_lds_ci[[object@interested_subtype]])
  genome_ci$camod <- rownames(position)
  sda_ds_ci_rank <-  ggplot(data = genome_ci[position$classification_combine == object@interested_subtype, ] %>%
                              dplyr::mutate(name = factor(camod, levels = camod[order(sda_lds_interested, decreasing = T)])),
                            aes(y = name, x = sda_lds_interested, xmin = lower, xmax = upper)) +
    geom_point() +
    geom_errorbarh(height=.1) +
    theme_bw() +
    ylab("") +
    xlab("") +
    labs(title = paste0("DS_SDA to ", object@interested_subtype, " center")) +
    scale_x_continuous(breaks = c(seq(from = floor(min(genome_ci$lower)), to = log2(1), length.out = 8), log2(2), log2(3)),
                       labels = round(2^c(seq(from = floor(min(genome_ci$lower)), to = log2(1), length.out = 8), log2(2), log2(3)), 2))

  object@genome_figure <- list(sda_project_position = sda_project_position,
                               sda_ds_ci_rank = sda_ds_ci_rank)

  return(object)
}


