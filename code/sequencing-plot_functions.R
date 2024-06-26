# These consolidate all common graphs so they call all be changed simultaneously

### 01 Preprocess ###

## Knee plot ##
# Made de novo from ggplot, keeps all features consistent
# Takes as input an array of cells, along with reads and rank
knee_plot <- function(x){

  plot <- ggplot(x$Metrics, aes(x = Rank, y = UMIs)) +
                 geom_point(aes(colour = FDR, shape = Retain)) +
                 geom_hline(yintercept = x$Knee_threshold, colour = "dodgerblue", linetype = "dashed") +
                 annotate("text", x = 0, y = x$Knee_threshold,
                          label = paste("Knee =", x$Knee_threshold),
                          colour = "dodgerblue", hjust = 0, vjust = -1) +
                 geom_hline(yintercept = x$Inflection_threshold,
                            colour = "forestgreen", linetype = "dashed") +
                 annotate("text", x = 0, y = x$Inflection_threshold,
                          label = paste("Inflection =", x$Inflection_threshold),
                          colour = "forestgreen", hjust = 0, vjust = -1) +
                  geom_hline(yintercept = 100, colour = "darkorchid", linetype = "dashed") +
                  annotate("text", x = 0, y = 100, label = "Empty threshold = 100",
                           colour = "darkorchid", hjust = 0, vjust = -1) +
                  scale_colour_gradientn(colours = heatmap_colors) +
                  scale_x_log10(labels = scales::number) + scale_y_log10(labels = scales::number) +
                  ylab("Total counts") +
                  theme_minimal()
  return(plot)
}

knee_plot_atac <- function(x){

  plot <- ggplot(x$Metrics, aes(x = Rank, y = Frags)) +
    geom_point(aes(colour = Retain)) +
    geom_hline(yintercept = x$Min_frags, colour = "darkorchid", linetype = "dashed") +
    annotate("text", x = 0, y = x$Min_frags, label = paste("Minimum peak fragments =", x$Min_frags),
             colour = "darkorchid", hjust = 0, vjust = -1) +
    scale_colour_discrete(type = filter_colors) +
    scale_x_log10(labels = scales::number) + scale_y_log10(labels = scales::number) +
    ylab("Total fragments in peaks") +
    theme_minimal()
  return(plot)
}

## Filter plot ##
# Made de novo from ggplot, keeps all features consistent
# Takes as input an array of cells, with filter status and filter characteristics
filter_plot <- function(x, x_col, y_col, colors){
  # do this to avoid having all show up as red if none get filtered out
  if (length(levels(as.factor(x[[colors]]))) == 1) {
    filter_color <- "dodgerblue"
  } else if (length(levels(as.factor(x[[colors]])))== 2) {
    filter_color <- filter_colors
  } else {

  }

  # tailors the labels and filters to data type
  if (colors == ">1000 Genes"){
    line_type <- geom_vline(xintercept = 1000, colour = "black", linetype = "dashed")
  } else if (colors == "<median+5xMAD"){
    line_type <- geom_hline(yintercept = median(x$UMIs)+5*mad(x$UMIs),
                            colour = "black", linetype = "dashed")
  } else if (colors == "<Q3+5xIQR"){
    line_type <- geom_hline(yintercept = quantile(x[[y_col]])[4]+
                              5*IQR(x[[y_col]]),
                            colour = "black", linetype = "dashed")
  } else if (colors == "<Med+5MAD Frags"){
    line_type <- geom_hline(yintercept = median(x$Fragments)+5*mad(x$Fragments),
                            colour = "black", linetype = "dashed")
  } else if (colors == ">2000 Peak Fragments"){
    line_type <- geom_hline(yintercept = 2000, colour = "black", linetype = "dashed")
  } else if (colors == "FRIP>40"){
    line_type <- geom_vline(xintercept = 40, colour = "black", linetype = "dashed")
  } else if (colors == "NS<4"){
    line_type <- geom_vline(xintercept = 4, colour = "black", linetype = "dashed")
  } else if (colors == "TSSe>3"){
    line_type <- geom_vline(xintercept = 3, colour = "black", linetype = "dashed")
  } else if (colors == "Top 10%"){
    line_type <- geom_vline(xintercept = .9, colour = "black", linetype = "dashed")
  } else {
    line_type <- geom_blank()
  }

  plot <- ggplot(x, aes(.data[[x_col]], .data[[y_col]], colour = .data[[colors]])) +
    line_type +
    scale_y_log10(labels = scales::number) +
    scale_x_continuous(labels = scales::number) +
    scale_colour_manual(values = filter_color) +
    geom_point(alpha = 0.3, size = 0.3) +
    guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    theme_minimal() + theme(legend.title=element_text(size=9))
  return(plot)
}

## Filter histogram ##
# Made de novo from ggplot, keeps all features consistent
# Takes as input an array of cells, with filter status and filter characteristics
filter_histogram <- function(x, attribute, filter_type, bin_size, filter_point){
  # do this to avoid having all show up as red if none get filtered out
  if (length(levels(as.factor(x[[filter_type]]))) == 1) {
    filter_color <- "dodgerblue"
  } else if (length(levels(as.factor(x[[filter_type]]))) == 2) {
    filter_color <- filter_colors
  }

  plot<-ggplot(x, aes(.data[[attribute]], fill = .data[[filter_type]])) +
    geom_histogram(binwidth = bin_size) +
    geom_vline(xintercept = filter_point, colour = "black", linetype = "dashed") +
    scale_y_log10(labels = scales::number, expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    scale_fill_manual(values = filter_color) +
    ylab("Count") +
    theme_minimal() + theme(legend.title=element_text(size=9))
  return(plot)
}


## Filter progression ##
# Made de novo from ggplot, keeps all features consistent
# Takes as input an array of cells with filter status for each
filter_progression <- function(x, filter_type){

  if (filter_type == "Total_nuclei"){
    y_label <- "Nuclei passing filter"
  } else if (filter_type == "Mean_counts"){
    y_label <- "Mean UMIs of nuclei passing filter"
  } else if (filter_type == "Mean_features"){
    y_label <- "Mean features of nuclei passing filter"
  } else if (filter_type == "Mean_mito"){
    y_label <- "Mean mitochondrial proportion of nuclei passing filter"
  } else if (filter_type == "Mean_ribo"){
    y_label <- "Mean ribosomal proportion of nuclei passing filter"
  } else if (filter_type == "Mean_frags"){
    y_label <- "Mean high-quality fragments of nuclei passing filter"
  } else if (filter_type == "Mean_peak_frags"){
    y_label <- "Mean fragments in peaks of nuclei passing filter"
  } else if (filter_type == "Mean_FRIP"){
    y_label <- "Mean FRIP of nuclei passing filter"
  } else if (filter_type == "Mean_NS"){
    y_label <- "Mean nucleosome signal of nuclei passing filter"
  } else if (filter_type == "Mean_TSS"){
    y_label <- "Mean TSS enrichment of nuclei passing filter"
  }

  plot <- ggplot(x, aes(Filter_step, .data[[filter_type]], group = 1)) +
                 geom_line(linetype = "dashed", colour = "forestgreen") +
                 geom_point(colour = "dodgerblue") +
                 geom_text(aes(label= .data[[filter_type]]),hjust=-.1, vjust=-.3, colour = "darkorchid") +
                 ylab(y_label) + xlab("Filter step") +
                 theme_minimal()
  return(plot)
}

### 02/04 Process ###

## Total nuclei plot ##
# Made de novo from ggplot, keeps all features consistent
# Takes as input an array of genes in terms of mean, variance, amount of cells expressing

total_nuclei_plot <- function(object){

  plot <- ggplot(object, aes(Region)) +
                 ylab("Total nuclei") +
                 geom_bar(aes(fill = Batch)) + scale_y_continuous(expand = c(0,0)) +
                 scale_fill_manual(values = batch_colors_rna) +
                 theme_classic() +  theme(axis.ticks.x = element_blank())
  return(plot)
}

total_nuclei_plot_atac <- function(object){

  plot <- ggplot(object, aes(Region)) +
    ylab("Total nuclei") +
    geom_bar(aes(fill = Batch)) + scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = batch_colors_atac) +
    theme_classic() +  theme(axis.ticks.x = element_blank())
  return(plot)
}

## Scatter plot for genes vs UMIs ##
# Made de novo from ggplot, keeps all features consistent
# Takes as input an array of genes in terms of mean, variance, amount of cells expressing

count_feature_plot <- function(object, group_id, colors){

  # needs aes string to specify correctly
  plot <- ggplot(object, aes(UMIs, Genes, colour = .data[[group_id]])) +
                 geom_point(alpha = 0.5, size = 0.3) +
                 geom_hline(yintercept = mean(object$Genes), colour = 'darkorchid', linetype = 2) +
                 geom_vline(xintercept = mean(object$UMIs), colour = 'forestgreen', linetype = 2) +
                 scale_colour_manual(values = colors) +
                 guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
                 labs(colour = group_id) + theme_minimal()
  return(plot)
}

count_feature_plot_atac <- function(object, group_id, colors){

  # needs aes string to specify correctly
  plot <- ggplot(object, aes(`Peak fragments`, Peaks, colour = .data[[group_id]])) +
    geom_point(alpha = 0.5, size = 0.3) +
    geom_hline(yintercept = mean(object$Peaks), colour = 'darkorchid', linetype = 2) +
    geom_vline(xintercept = mean(object$`Peak fragments`), colour = 'forestgreen', linetype = 2) +
    scale_colour_manual(values = colors) +
    guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    labs(colour = group_id) + theme_minimal()
  return(plot)
}

## Violin plot for metrics ##
# Formats color schemes together
# Takes seurat array as input

violin_metrics_formatted <- function(x, group, metric, label, colors){

  plot <- VlnPlot(x, features = metric, group.by = group, pt.size = 0) +
    ggtitle('') + xlab('') + ylab(label) +
    scale_fill_manual(values = colors) +
    NoLegend()
  return(plot)
}

## Scatter plot for gene filter in cells ##
# Made de novo from ggplot, keeps all features consistent
# Takes as input an array of genes in terms of mean, variance, amount of cells expressing

gene_filter_formatted <- function(object){

  plot <- ggplot(object, aes(x = Mean, y = Total_cells, colour = `>3 Nuclei`)) +
                 geom_point(alpha = 0.3, size = 0.3) +
                 geom_hline(yintercept = 3,
                            colour = "black", linetype = "dashed") +
                 scale_x_log10() + scale_y_log10() +
                 guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
                 xlab("Mean Expression") + ylab(paste("Nuclei expressing")) +
                 scale_colour_manual(values = filter_colors) +
                 theme_minimal() + theme(legend.title=element_text(size=10))
  return(plot)
}

## SCT Gene distribution plot ##
# Made de novo from ggplot, keeps all features consistent
# Takes as input an array of genes in terms of mean, variance, HVG, SCTransform model

sct_plot <- function(object, x_col, y_col){

  if (x_col == "Mean"){
    x_label <- "Mean Expression"
  } else if (x_col == "gmean"){
    x_label <- "Geometric Mean"
  }

  if (y_col == "CV2"){
    y_label <- "CV2"
  } else if (y_col == "Variance_raw"){
    y_label <- "Variance"
  } else if (y_col == "variance"){
    y_label <- "Variance"
  } else if (y_col == "residual_variance"){
    y_label <- "Residual Variance"
  }

  plot <- ggplot(object, aes(x = .data[[x_col]], y = .data[[y_col]], colour = HVG)) +
                 geom_point(alpha = 0.3, size = 0.3) +
                 scale_colour_manual(values = filter_colors) +
                 guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
                 scale_x_log10() + scale_y_log10() +
                 xlab(x_label) + ylab(y_label) +
                 theme_minimal()
  return(plot)
}

## Variance plot for PCs ##
# Made de novo from ggplot, keeps all features consistent
# Takes PCs and variance explained as an input

variance_explained_plot <- function(object){

  plot <- ggplot(object, aes(x = PC)) +
                 geom_point(aes(y = Proportion_variance, colour = "Individual")) +
                 geom_line(aes(y = Cumulative_variance/10, colour = "Cumulative")) +
                 scale_y_continuous(sec.axis = sec_axis(~.*10, name = "Cumulative Variance Explained"),
                                    name = "Individual Variance Explained") +
                 scale_colour_manual(name = "Variance", values = filter_colors) +
                 theme_bw()
  return(plot)
}

# adapted from Harvard Chan Bioinformatics Core
improved_elbow_plot <- function(object){

  plot_df <- data.frame(pct = object$Percent_variance,
                        cumu = object$Cumulative_variance,
                        rank = 1:length(object$Percent_variance))

  # Elbow plot to visualize
  plot <- ggplot(plot_df, aes(cumu, pct, label = rank, colour = rank <= object$Optimal_PCs)) +
    geom_text() +
    geom_vline(xintercept = 90, color = "black", linetype = "dashed") +
    geom_hline(yintercept = min(object$Percent_variance[object$Percent_variance > 5]), color = "black", linetype = "dashed") +
    scale_colour_manual(values = filter_colors, name = "Included") +
    xlab("Cumulative Variance") + ylab("Percent Variance") +
    theme_classic()
  return(plot)

}

### 03/05 Classify ###

## Resolution plot ##
# Formats labels and legends
# Takes seurat array as input

resolution_plot_formatted <- function(x, dimred, res, metric){

  if (dimred == 'pca'){
    x_label <- "PC 1"
    y_label <- "PC 2"
  } else if (dimred == "tsne"){
    x_label <- "tSNE 1"
    y_label <- "tSNE 2"
  } else if (dimred == "umap"){
    x_label <- "UMAP 1"
    y_label <- "UMAP 2"
  } else {
    stop("Incorrect reduction")
  }

  plot <- DimPlot(x, reduction = dimred, group.by = paste0(metric, '_snn_res.', res)) +
                  xlab(x_label) + ylab(y_label) + labs(colour = 'Cluster') +
                  theme_classic() +
                  theme(plot.title = element_blank(), legend.text=element_text(size=8),
                        axis.ticks = element_blank(), axis.text = element_blank())
  return(plot)
}

## Cluster plot ##
# Formats labels and legends
# Takes seurat array as input

cluster_plot_formatted <- function(x, dimred, clusters, colors){

  if (dimred == 'pca'){
    x_label <- "PC 1"
    y_label <- "PC 2"
  } else if (dimred == "tsne"){
    x_label <- "tSNE 1"
    y_label <- "tSNE 2"
  } else if (dimred == "umap"){
    x_label <- "UMAP 1"
    y_label <- "UMAP 2"
  } else {
    stop("Incorrect reduction")
  }

  plot <- DimPlot_scCustom(x, reduction = dimred, group.by = clusters) +
                  xlab(x_label) + ylab(y_label) +
                  scale_colour_manual(values = colors) +
                  theme_classic() +
                  theme(plot.title = element_blank(), legend.text=element_text(size=8))
  return(plot)
}

## Feature plot ##
# Builds in legend and title formatting, as well as color schemes
# Takes seurat array as input

feature_plot_formatted <- function(x, dimred, marker, colorset, na_cutoff){

  if (dimred == 'pca'){
    x_label <- "PC 1"
    y_label <- "PC 2"
  } else if (dimred == "tsne"){
    x_label <- "tSNE 1"
    y_label <- "tSNE 2"
  } else if (dimred == "umap"){
    x_label <- "UMAP 1"
    y_label <- "UMAP 2"
  } else {
    stop("Incorrect reduction")
  }

  if (marker == 'nCount_RNA'){
    marker_lab <- "UMIs"
  } else if (marker == 'nFeature_RNA'){
    marker_lab <- "Genes"
  } else if (marker ==  'Mito_proportion'){
    marker_lab <- "% Mito"
  } else if (marker == 'Ribo_proportion'){
    marker_lab <- "% Ribo"
  } else {
    marker_lab <- marker
  }

  plot <- FeaturePlot_scCustom(x, reduction = dimred, features = marker, colors_use = colorset, na_cutoff = na_cutoff, order = FALSE) +
                      xlab(x_label) + ylab(y_label) + labs(colour = marker_lab) +
                      theme_classic() +
                      theme(plot.title = element_blank())
  return(plot)
}

## Compositional bar plots showing absolute numbers of cells ##
# Formats color schemes together
# Takes arrays of cell identities as input

absolute_bar_plot_formatted <- function(x, bar_group, stack_division, colors){

  plot <- ggplot(x, aes(.data[[bar_group]])) +
                 geom_bar(aes(fill = .data[[stack_division]])) +
                 scale_y_continuous(expand = c(0,0)) +
                 scale_fill_manual(values = colors) + labs(fill = '') +
                 theme_classic() +
                 theme(axis.ticks.x = element_blank(),
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  return(plot)
}

## Compositional bar plots showing relative number in division compared to averages ##
# Formats color schemes together
# Takes arrays of cell identities as input
# Hard-coded to be based on region

relative_bar_plot_formatted <- function(x, bar_group, stack_group, colors){

  plot_prop_array <- table(x[[bar_group]], x[[stack_group]])
  plot_prop_array <- as.data.frame(plot_prop_array / rowSums(plot_prop_array))

  plot <- ggplot(plot_prop_array, aes(x = Var1, y = Freq, fill = Var2)) +
                 geom_bar(position="stack", stat="identity") + scale_y_continuous(expand = c(0,0)) +
                 scale_fill_manual(values = colors) +
                 labs(fill = '') + xlab('') + ylab('Proportion in region') +
                 theme_classic() +
                 theme(axis.ticks.x = element_blank(),
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  return(plot)
}

relative_bar_plot_formatted_ranked <- function(x, bar_group, stack_group, colors){

  plot_prop_array <- table(x[[bar_group]], x[[stack_group]])
  plot_prop_array <- as.data.frame(plot_prop_array / rowSums(plot_prop_array))

  plot_prop_array %>%
    filter(Var2 == "aplCoA") %>%
    mutate(x = factor(Var1, levels = Var1[order(Freq)])) %>%
    pull(x) %>%
    levels -> x_levels

  plot_prop_array <- plot_prop_array %>%
                     mutate(Var1 = factor(Var1, levels = x_levels))

  plot <- ggplot(plot_prop_array, aes(x = Var1, y = Freq, fill = Var2)) +
    geom_bar(position="stack", stat="identity") + scale_y_continuous(expand = c(0,0)) +
    scale_fill_manual(values = colors) +
    labs(fill = '') + xlab('') + ylab('Proportion in region') +
    theme_classic() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  return(plot)
}

## Cluster distance plot ##
# Takes median distance between clusters and the dendrogram to show relatedness between different clusters

dendrogram_heatmap <- function(cluster_dend_list, method_id){

  # create distance matrix based on the cluster medians at most variable genes
  if (method_id == "distance"){
    dist_medians <- dist(t(cluster_dend_list[[1]]))
    color_order <- rev(scatter_intensity_colors)
  }
  else if (method_id == "correlation"){
    dist_medians <- cor(cluster_dend_list[[1]], method = "spearman")
    color_order <- scatter_intensity_colors
  }

  plot <- heatmap3(as.matrix(dist_medians),
                   Rowv = cluster_dend_list[[2]], Colv = cluster_dend_list[[2]],
                   trace = "none", col = color_order, cexRow=0.8, cexCol=0.8, scale = "none",
                   seq(0.2, 1, length.out = 99))

  return(plot)
}

## DEG abundance plot ##
# Lists number of DEGs per region for the level of organization of interest
deg_abundance_plot <- function(deg_list, color_list, initial_array, axis_limits){

  deg_array <- unlist(lapply(deg_list, function(x) length(rownames(x))))

  deg_array <- cbind.data.frame(deg_array,
                                names(deg_array))
  colnames(deg_array) <- c("DEGs", "Class")

  deg_array <- deg_array[match(levels(initial_array), deg_array$Class),]
  deg_array$Class <- factor(deg_array$Class, levels = unique(deg_array$Class))

  plot <- ggplot(deg_array, aes(x = Class, y = DEGs, fill = Class)) +
                 geom_bar(stat="identity", position=position_dodge()) +
                 scale_y_continuous(expand = c(0,0), limits = axis_limits) +
                 scale_fill_manual(values = color_list) +
                 xlab("") +
                 theme_classic() +
                 theme(axis.ticks.x = element_blank(),
                       axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                       legend.position = "none")

  return(plot)
}

deg_abundance_plot_magnitude <- function(deg_list, color_list, axis_limits, region){

  if (region == "aplCoA"){
    deg_array <- cbind.data.frame(c(rep(1.5, times = length(names(deg_list))),
                                    rep(2.5, times = length(names(deg_list))),
                                    rep(5, times = length(names(deg_list)))),
                                  rep(names(deg_list), times = 3),
                                  c(unlist(lapply(deg_list,
                                                  function(x) length(rownames({x[x$p_val_adj < .05
                                                                                 & x$avg_logFC > 0.585
                                                                                 & x$avg_logFC < 1.322,]})))),
                                    unlist(lapply(deg_list,
                                                  function(x) length(rownames({x[x$p_val_adj < .05
                                                                                 & x$avg_logFC > 1.322
                                                                                 & x$avg_logFC < 2.322,]})))),
                                    unlist(lapply(deg_list,
                                                  function(x) length(rownames({x[x$p_val_adj < .05
                                                                                 & x$avg_logFC > 2.322,]})))))
    )
  } else if (region == "pplCoA"){
    deg_array <- cbind.data.frame(c(rep(1.5, times = length(names(deg_list))),
                                    rep(2.5, times = length(names(deg_list))),
                                    rep(5, times = length(names(deg_list)))),
                                  rep(names(deg_list), times = 3),
                                  c(unlist(lapply(deg_list,
                                                  function(x) length(rownames({x[x$p_val_adj < .05
                                                                                 & x$avg_logFC < -0.585
                                                                                 & x$avg_logFC > -1.322,]})))),
                                    unlist(lapply(deg_list,
                                                  function(x) length(rownames({x[x$p_val_adj < .05
                                                                                 & x$avg_logFC < -1.322
                                                                                 & x$avg_logFC > -2.322,]})))),
                                    unlist(lapply(deg_list,
                                                  function(x) length(rownames({x[x$p_val_adj < .05
                                                                                 & x$avg_logFC < -2.322,]})))))
    )
  }

  colnames(deg_array) <- c("Magnitude", "Class", "DEGs")
  deg_array$Magnitude <- factor(deg_array$Magnitude, levels = c("1.5", "2.5", "5"))

  plot <- ggplot(deg_array, aes(x = Class, y = DEGs, fill = Magnitude)) +
    geom_bar(position="stack", stat="identity") +
    scale_y_continuous(expand = c(0,0), limits = axis_limits) +
    scale_fill_manual(values = color_list) +
    xlab("") +
    theme_classic() +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position = "none")
  return(plot)
}

## Stability plot ##
# Made de novo from ggplot, keeps all features consistent

stability_plot <- function(clustree_stability_object, metric_name){

  stable_cluster <- clustree_stability_object$data
  stable_plot <- aggregate(stable_cluster[,7], list(stable_cluster[[metric_name]]), mean)
  stable_plot <- cbind.data.frame(stable_plot, aggregate(res_sil, list(stable_cluster[[metric_name]]), mean)[,2])
  colnames(stable_plot) <- c("resolution", "Stability", "Silhouette")
  stable_plot <- stable_plot %>% pivot_longer(!resolution)

  plot <- ggplot(stable_plot, aes(x = resolution, y = value, group = name, colour = name)) +
    geom_line() + geom_point() +
    scale_colour_manual(values = filter_colors) +
    guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    ylab("Mean cluster metric") + xlab("Leiden resolution") + labs(colour = 'Metric') +
    theme_classic() + theme(legend.title=element_text(size=10))
  return(plot)
}

silhouette_plot <- function(x, type, res){

  mean_sil <- mean(x@meta.data[[paste0("silhouette.", res)]])
  plot <- x@meta.data %>%
    mutate(barcode = rownames(.)) %>%
    arrange(.data[[paste0(type, "_snn_res.", res)]],-.data[[paste0("silhouette.", res)]]) %>%
    mutate(barcode = factor(barcode, levels = barcode)) %>%

    ggplot() + geom_col(aes(barcode, .data[[paste0("silhouette.", res)]], fill = .data[[paste0(type, "_snn_res.", res)]]),
                        show.legend = FALSE) +
    geom_hline(yintercept = mean_sil, color = 'black', linetype = 'dashed') +
    scale_x_discrete(name = 'Cells') +
    scale_y_continuous(name = 'Silhouette Score') +
    theme_classic() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  return(plot)
}

## Volcano plot ##
# Modified EnhancedVolcano

enhanced_volcano_plot <- function(x, title_id, label_size){

  if (class(x$avg_logFC) == "numeric"){
    plot <- EnhancedVolcano(x,
                            title = title_id,
                            subtitle = 'P-value < .01, Magnitude change > log2(0.585)-fold (1.5 or 50%)',
                            pointSize = 2,
                            labSize = label_size,
                            lab = rownames(x),
                            pCutoff = 1e-02,
                            FCcutoff = 0.585,
                            x = 'avg_logFC',
                            y = 'p_val_adj',
                            gridlines.major = FALSE,
                            gridlines.minor = FALSE)
  } else {
    plot <- NULL
  }

  return(plot)
}

tree_plot <- function(x, label, colors) {
  Idents(x) <- "Tissue"
  x <- BuildClusterTree(x,
                        dims = 1:50,
                        verbose = FALSE)

  tree <- x@tools$BuildClusterTree

  plot <- ggtree(tree, aes(x, y)) +
    scale_y_reverse() +
    geom_tree() +
    theme_tree() +
    geom_tiplab(offset = 1) +
    geom_tippoint(color = colors, shape = 16, size = 5) +
    coord_cartesian(clip = 'off') +
    theme(plot.margin = unit(c(0,2.5,0,0), 'cm'))

  return(plot)
}

## VG1/VG2 coexpression plot ##

vg1vg2_coexpression_plot_region <- function(x){

  subset_array <- FetchData(array, slot = "counts", vars = c("Slc17a6", "Slc17a7", "Region"))

  plot <- ggplot(subset_array, aes(x = Slc17a6, y = Slc17a7, colour = Region)) +
    geom_point(size = 0.1) + geom_jitter(size = 0.1) +
    scale_colour_manual(values = region_colors) +
    guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    xlab("Slc16a6 UMIs") + ylab("Slc17a7 UMIs") +
    theme_minimal()

  return(plot)
}

vg1vg2_coexpression_plot <- function(x, region){

  subset_array <- x[, as.factor(x$Region) %in% region]
  subset_array <- FetchData(subset_array, slot = "counts", vars = c("Slc17a6", "Slc17a7", "Type"))

  plot <- ggplot(subset_array, aes(x = Slc17a6, y = Slc17a7, colour = Class)) +
    geom_point(size = 0.1) + geom_jitter(size = 0.1) +
    scale_colour_manual(values = class_colors) +
    guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    xlab("Drd1 UMIs") + ylab("Drd2 UMIs") +
    theme_minimal()

  return(plot)
}

vg1vg2_coexpression_plot_batch <- function(x, batch){

  subset_array <- x[, as.factor(x$orig.ident) %in% batch]
  subset_array <- FetchData(subset_array, slot = "counts", vars = c("Drd1", "Drd2", "Type"))

  plot <- ggplot(subset_array, aes(x = Slc17a6, y = Slc17a7, colour = Class)) +
    geom_point(size = 0.1) + geom_jitter(size = 0.1) +
    scale_colour_manual(values = class_colors) +
    guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    xlab("Drd1 UMIs") + ylab("Drd2 UMIs") +
    theme_minimal()

  return(plot)
}
