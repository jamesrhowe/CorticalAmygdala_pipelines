# contains helper functions for analysis notebooks

# this used to collect summary stats for every dataset
count_all <- function(step, object){
  cat(count_stats(paste("counts per", step), object@meta.data$nCount_RNA),
      count_stats(paste("features per", step), object@meta.data$nFeature_RNA),
      count_stats(paste("mitochondrial proportion per", step), object@meta.data$Mito_proportion),
      count_stats(paste("ribosomal proportion per", step), object@meta.data$Ribo_proportion), sep = "\n")
}

count_all_atac <- function(step, object){
  cat(count_stats(paste("fragments per", step), object@meta.data$passed_filters),
      count_stats(paste("fragments in peaks per", step), object@meta.data$peak_region_fragments),
      count_stats(paste("fraction reads in peaks per", step), object@meta.data$fraction_reads_in_peaks),
      count_stats(paste("nucleosome signal per", step), object@meta.data$nucleosome_signal),
      count_stats(paste("TSS enrichment per", step), object@meta.data$TSS.enrichment), sep = "\n")
}

count_stats <- function(x, y){
  cat(paste("Min", x, round(min(y), 2)),
      paste("Max", x, round(max(y), 2)),
      paste("Mean", x, round(mean(y), 2)),
      paste("Median", x, round(median(y), 2)),
      paste("Std. dev.", x, round(sd(y), 2)),
      paste("Std. error", x, round(sd(y)/sqrt(length(y)), 2)),
      paste(),
      sep = "\n")
}

# Find markers

## Cluster markers (wilcoxon, single-cell)
markers_sc_wilcox <- function(x, label, verbosity = TRUE){
  Idents(x) <- label
  marker_list <- FindAllMarkers(object = x, logfc.threshold = 0.585, only.pos = TRUE)
  marker_list <- cbind.data.frame(marker_list$gene, marker_list$p_val_adj, marker_list$cluster)
  colnames(marker_list) <- c("gene", "combined_p_val", "cluster")
  marker_list <- split(marker_list, marker_list$cluster)
  marker_list <- lapply(marker_list, function(x) as.data.frame(x[,1:2]))
  return(marker_list)
}

## Region markers (edgeR-LRT, pseudobulk)
markers_pb_edgeRLRT <- function(x, label, pct_cutoff){
  metadata_col <- ifelse(x$Region == levels(as.factor(x$Region))[1], levels(as.factor(x$Region))[1], "others")
  x <- AddMetaData(x, metadata = metadata_col, col.name = "test_regions")
  marker_list <- run_de(x, cell_type_col = label, replicate_col = "orig.ident", label_col = "test_regions")
  marker_list <- marker_list[,1:5]
  marker_list <- split(marker_list, marker_list$cell_type)
  marker_list <- lapply(marker_list, function(x) as.data.frame(x[,2:5]))
  for (i in 1:length(marker_list)){
    rownames(marker_list[[i]]) <- marker_list[[i]]$gene
    marker_list[[i]] <- marker_list[[i]][,2:4]
    pct_expressed <- pct_expressed[rownames(pct_expressed) %in% rownames(marker_list[[i]]),]
    marker_list[[i]] <- marker_list[[i]][pct_expressed[, i] >= pct_cutoff,]
  }

  return(marker_list)
}

# adapted from Harvard Chan Bioinformatics Core instructions
find_optimum_pcs <- function(pc_array){

  pct <- pc_array[["pca"]]@stdev / sum(pc_array[["pca"]]@stdev) * 100

  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)

  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  metric_1 <- which(cumu > 90 & pct < 5)[1]

  # Determine the difference between variation of PC and subsequent PC, last point where change of % of variation is more than 0.1%.
  metric_2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

  pcs <- min(metric_1, metric_2)

  pc_metrics <- list(pct, cumu, metric_1, metric_2, pcs)
  names(pc_metrics) <- c("Percent_variance", "Cumulative_variance", "CumVar>90", "Diff<0.1", "Optimal_PCs")

  return(pc_metrics)
}

# adapted from Harvard Chan Bioinformatics Core instructions
find_optimum_pcs_atac <- function(pc_array){

  pct <- pc_array[["lsi_bin"]]@stdev / sum(pc_array[["lsi_bin"]]@stdev) * 100

  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)

  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  metric_1 <- which(cumu > 90 & pct < 5)[1]

  # Determine the difference between variation of PC and subsequent PC, last point where change of % of variation is more than 0.1%.
  metric_2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

  pcs <- min(metric_1, metric_2)

  pc_metrics <- list(pct, cumu, metric_1, metric_2, pcs)
  names(pc_metrics) <- c("Percent_variance", "Cumulative_variance", "CumVar>90", "Diff<0.1", "Optimal_PCs")

  return(pc_metrics)
}

targeted_comp_markers <- function(x, i1, i2, gene_vector){

  gene_list <- FindMarkers(x, test.use = "t", only.pos = TRUE, logfc.threshold = 0.5, verbose = FALSE,
                           ident.1 = i1, ident.2 = i2)
  gene_list <- gene_list[gene_list$p_val_adj < 0.05,]
  gene_list <- gene_list[gene_list$pct.1 > 2*gene_list$pct.2,]
  gene_list <- gene_list[!rownames(gene_list) %in% gene_vector,]
}

pc_list <- c(1, 2, 3, 50)

# Dimensionality reduction IDs
dr_params <- cbind.data.frame(c("PCA", "tSNE", "UMAP"),
                                  c("pca", "tsne", "umap"))
colnames(dr_params) <- c("name", "type")


vln_metrics <- cbind.data.frame(c("nCount_RNA", "nFeature_RNA", "Mito_proportion", "Ribo_proportion"),
                                    c("UMIs", "Genes", "% Mito", "% Ribo"))
colnames(vln_metrics) <- c("id", "label")

batch_region <- c("Batch", "Region")
identity_region <- c("orig.ident", "Region")
tissue_class <- c("Tissue", "Class")# contains helper functions for analysis notebooks

# this used to collect summary stats for every dataset
count_all <- function(step, object){
  cat(count_stats(paste("counts per", step), object@meta.data$nCount_RNA),
      count_stats(paste("features per", step), object@meta.data$nFeature_RNA),
      count_stats(paste("mitochondrial proportion per", step), object@meta.data$Mito_proportion),
      count_stats(paste("ribosomal proportion per", step), object@meta.data$Ribo_proportion), sep = "\n")
}

count_stats <- function(x, y){
  cat(paste("Min", x, round(min(y), 2)),
      paste("Max", x, round(max(y), 2)),
      paste("Mean", x, round(mean(y), 2)),
      paste("Median", x, round(median(y), 2)),
      paste("Std. dev.", x, round(sd(y), 2)),
      paste("Std. error", x, round(sd(y)/sqrt(length(y)), 2)),
      paste(),
      sep = "\n")
}

create_dendrogram <- function(array_id, label_id, color_id){

  # figure out cluster medians
  medians <- do.call("cbind",
                     tapply(names(label_id),
                            label_id,
                            function(x){matrixStats::rowMedians(as.matrix(GetAssayData(object = array_id, slot = "data")[,x]))
                            }))

  pvclust_matrix <- pvclust(data = medians, method.dist = "cor", method.hclust = "average", nboot = 100, parallel = TRUE)

  dend <- as.dendrogram(pvclust_matrix$hclust)

  color_labels <- setNames(color_id, levels(label_id))

  dend <- dend %>% set("labels_cex", 0.7)
  dend <- dend %>% set("leaves_pch", 19) %>% set("leaves_cex", 0.5)
  dend <- dend %>% set("labels_col", color_labels[labels(dend)])
  dend <- dend %>% set("leaves_col", color_labels[labels(dend)])

  output <- list(medians, dend)
  names(output) <- c("Medians", "Dendrogram")

  return(output)
}

test_region_specific_expression <- function(x, gene){

  gene_expression <- FetchData(x, vars = c(gene, "orig.ident", "Region"))
  gene_expression <- aggregate(gene_expression[[gene]] ~ Region + orig.ident, data = gene_expression, FUN = "mean")

  colnames(gene_expression) <- c("Region", "Batch", gene)

  anova_gene <- aov(gene_expression[[gene]] ~ Region, data = gene_expression)

  gene_test_list <- list(gene_expression,
                         summary(anova_gene),
                         TukeyHSD(anova_gene))
  names(gene_test_list) <- c("Summary", "ANOVA", "Tukey test")

  return(gene_test_list)
}

pseudobulk_pca <- function(x, identifier, label, colors){

  pseudobulk_array <- AggregateExpression(x, group.by = identifier, assays = "RNA", slot = "counts", return.seurat = TRUE)
  pseudobulk_array <- AddMetaData(pseudobulk_array, metadata = as.factor(names(pseudobulk_array$orig.ident)), col.name = label)
  pseudobulk_array <- NormalizeData(pseudobulk_array, verbose = FALSE) # results are slightly different
  pseudobulk_array <- FindVariableFeatures(pseudobulk_array, nfeatures = 3000)
  pseudobulk_array <- ScaleData(pseudobulk_array, verbose = FALSE)
  pseudobulk_array <- RunPCA(pseudobulk_array, verbose = FALSE, npcs = 2)

  var_explained <- pseudobulk_array@reductions$pca@stdev ^ 2 /
    sum(matrixStats::rowVars(GetAssayData(pseudobulk_array, assay = "RNA", slot = "scale.data")))
  var_explained <- 100 * round(var_explained, 3)

  plot <- DimPlot_scCustom(pseudobulk_array, reduction = 'pca', group.by = label, pt.size = 5) +
    xlab(paste0("PC 1 (", var_explained[1], "%)")) + ylab(paste0("PC 2 (", var_explained[2], "%)")) +
    scale_colour_manual(values = colors) +
    theme_classic() +
    theme(plot.title = element_blank(), legend.text=element_text(size=8))
  return(plot)
}

vg1vg2_coexpression_stats <- function(x, region){

  subset_array <- x[, as.factor(x$Region) %in% region]
  subset_array <- FetchData(subset_array, slot = "counts", vars = c("Slc17a6", "Slc17a7"))

  print(table(subset_array$Slc17a6, subset_array$Slc17a7))
  print(table(subset_array$Slc17a6 >= 1, subset_array$Slc17a7 >= 1))

  print(paste(region, "VG2+/VG1-:",
              table(subset_array$Slc17a6 >= 1, subset_array$Slc17a7 == 0)[2,2],
              round(table(subset_array$Slc17a6 >= 1, subset_array$Slc17a7 == 0)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Slc17a6 == 0, subset_array$Slc17a7 == 0)[2,2]), digits = 4)))
  print(paste(region, "VG1+/VG2-:",
              table(subset_array$Slc17a6 == 0, subset_array$Slc17a7 >= 1)[2,2],
              round(table(subset_array$Slc17a6 == 0, subset_array$Slc17a7 >= 1)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Slc17a6 == 0, subset_array$Slc17a7 == 0)[2,2]), digits = 4)))
  print(paste(region, "VG2+/VG1+:",
              table(subset_array$Slc17a6 >= 1, subset_array$Slc17a7 >= 1)[2,2],
              round(table(subset_array$Slc17a6 >= 1, subset_array$Slc17a7 >= 1)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Slc17a6 == 0, subset_array$Slc17a7 == 0)[2,2]), digits = 4)))
}

d1d2_coexpression_stats_batch <- function(x, batch){

  subset_array <- x[, as.factor(x$orig.ident) %in% batch]
  subset_array <- FetchData(subset_array, slot = "counts", vars = c("Drd1", "Drd2"))

  print(table(subset_array$Drd1, subset_array$Drd2))
  print(table(subset_array$Drd1 >= 1, subset_array$Drd2 >= 1))

  print(paste(batch, "D1+/D2-:",
              table(subset_array$Drd1 >= 1, subset_array$Drd2 == 0)[2,2],
              round(table(subset_array$Drd1 >= 1, subset_array$Drd2 == 0)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Drd1 == 0, subset_array$Drd2 == 0)[2,2]), digits = 4)))
  print(paste(batch, "D1-/D2+:",
              table(subset_array$Drd1 == 0, subset_array$Drd2 >= 1)[2,2],
              round(table(subset_array$Drd1 == 0, subset_array$Drd2 >= 1)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Drd1 == 0, subset_array$Drd2 == 0)[2,2]), digits = 4)))
  print(paste(batch, "D1+/D2+:",
              table(subset_array$Drd1 >= 1, subset_array$Drd2 >= 1)[2,2],
              round(table(subset_array$Drd1 >= 1, subset_array$Drd2 >= 1)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Drd1 == 0, subset_array$Drd2 == 0)[2,2]), digits = 4)))
}

d1d2_coexpression_stats_celltype <- function(x, region, grouping, cell_type){

  subset_array <- x[, as.factor(x$Region) %in% region]
  subset_array <- FetchData(subset_array, slot = "counts", vars = c("Drd1", "Drd2", grouping))
  colnames(subset_array) <- c("Drd1", "Drd2", "group")
  subset_array <- subset_array[subset_array$group %in% cell_type,]
  subset_array <- subset_array[,1:2]

  print(table(subset_array$Drd1, subset_array$Drd2))
  print(table(subset_array$Drd1 >= 1, subset_array$Drd2 >= 1))

  print(paste(region, cell_type, "D1+/D2-:",
              table(subset_array$Drd1 >= 1, subset_array$Drd2 == 0)[2,2],
              round(table(subset_array$Drd1 >= 1, subset_array$Drd2 == 0)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Drd1 == 0, subset_array$Drd2 == 0)[2,2]), digits = 4)))
  print(paste(region, cell_type, "D1-/D2+:",
              table(subset_array$Drd1 == 0, subset_array$Drd2 >= 1)[2,2],
              round(table(subset_array$Drd1 == 0, subset_array$Drd2 >= 1)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Drd1 == 0, subset_array$Drd2 == 0)[2,2]), digits = 4)))
  print(paste(region, cell_type, "D1+/D2+:",
              table(subset_array$Drd1 >= 1, subset_array$Drd2 >= 1)[2,2],
              round(table(subset_array$Drd1 >= 1, subset_array$Drd2 >= 1)[2,2] /
                      (length(rownames(subset_array)) - table(subset_array$Drd1 == 0, subset_array$Drd2 == 0)[2,2]), digits = 4)))
}

d1d2_coexpression_stats_4group <- function(x, region_id){

  subset_array <- x[, as.factor(x$Region) %in% region_id]
  subset_array <- FetchData(subset_array, slot = "counts", vars = c("Drd1", "Drd2", "Class"))
  colnames(subset_array) <- c("Drd1", "Drd2", "group")
  subset_array <- subset_array[subset_array$group %in% c("GABA.D1", "GABA.D2"),]

  d1d2_vector <- NULL
  for (i in 1:length(rownames(subset_array))){
    if (subset_array[i,1] >= 1) {
      if (subset_array[i,2] >= 1){
        d1d2_vector[i] <- "dual D1/D2"
      } else {
        d1d2_vector[i] <- "D1"
      }
    } else {
      if (subset_array[i,2] >= 1){
        d1d2_vector[i] <- "D2"
      } else {
        d1d2_vector[i] <- "no D1/D2"
      }
    }
  }
  stat_table <- table(d1d2_vector, subset_array$group)
  stat_table <- stat_table[,1:2]

return(stat_table)
}

pc_list <- c(1, 2, 3, 30)

# Dimensionality reduction IDs
dr_params <- cbind.data.frame(c("PCA", "tSNE", "UMAP"),
                              c("pca", "tsne", "umap"))
colnames(dr_params) <- c("name", "type")

vln_metrics <- cbind.data.frame(c("nCount_RNA", "nFeature_RNA", "Mito_proportion", "Ribo_proportion"),
                                c("UMIs", "Genes", "% Mito", "% Ribo"))
colnames(vln_metrics) <- c("id", "label")

batch_region <- c("Batch", "Region")
identity_region <- c("orig.ident", "Region")

marker_features <- c("Syp", "Eno2", # pan-neuronal
                     "Slc17a7", "Slc17a6", # glutamatergic neuron
                     "Gad1", "Gad2", # GABAergic neuron
                     "Aldh1l1", "Aqp4", # astrocyte
                     "Tmem119", "Ptprc", # microglia
                     "Pdgfra", "Cspg4", #OPC
                     "Bmp4", "Enpp6", # NFOL
                     "Mog", "Mal", # oligodendrocyte/OLG
                     "Cldn5", "Adgrl4", # endothelial cell
                     "Kcnj8", "Abcc9", # mural cell
                     "Slc47a1", "Rspo3", # ABC
                     "Slc6a13", "Slc13a4") # VLMC

marker_features_half <- c("Syp",
                          "Slc17a7", "Slc17a6", # glutamatergic neuron
                          "Gad2", # GABAergic neuron
                          "Aqp4", # Astrocytes
                          "Tmem119", "Ptprc", # microglia/macrophages
                          "Pdgfra", # OPC
                          "Bmp4", # NFOL
                          "Mog", # MOL
                          "Cldn5", # Endo
                          "Kcnj8", # Mural
                          "Slc47a1", # ABC
                          "Slc6a13") # VLMC

marker_features_neuron <- c("Syp", "Eno2", # pan-neuronal
                            "Slc17a7", "Slc17a6", # glutamatergic neuron
                            "Gad1", "Gad2") # GABAergic neuron

marker_features_glia <- c("Aldh1l1", "Aqp4", # astrocyte
                          "Tmem119", "Ptprc", # microglia
                          "Pdgfra", "Cspg4", #OPC
                          "Bmp4", "Enpp6", # NFOL
                          "Mog", "Mal", # oligodendrocyte/OLG
                          "Cldn5", "Adgrl4", # endothelial cell
                          "Kcnj8", "Abcc9", # mural cell
                          "Slc47a1", "Rspo3", # ABC
                          "Slc6a13", "Slc13a4") # VLMC

marker_features_glia_half <- c("Aqp4", # Astrocytes
                               "Tmem119", "Ptprc", # microglia/macrophages
                               "Pdgfra", # OPC
                               "Bmp4", # NFOL
                               "Mog", # MOL
                               "Cldn5", # Endo
                               "Kcnj8", # Mural
                               "Slc47a1", # ABC
                               "Slc6a13") # VLMC

micro_markers <- c("Tmem119", "Ptprc", # general markers
                   "Cd83", "Csf1", # dendritic cells
                   "F13a1", "Aoah") # macrophages
