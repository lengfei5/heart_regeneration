
########################################################
########################################################
# Section : subclustering and annotation
# 
########################################################
########################################################
subclustering_manual.annotation = function(aa)
{
  ##########################################
  # transfer the cell labels from spliced data
  ##########################################
  Transferring_cellLabels_splicedData = FALSE
  if(Transferring_cellLabels_splicedData){
    RdataDir_spliced = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/results/Rdata/'
    bb = readRDS(file = paste0(RdataDir_spliced, 'seuratObject_axloltl_scRNAseq_splicedOnly_manualAnnotation.rds'))
    mm = match(colnames(aa), colnames(bb))
    aa$celltypes = bb$celltypes[mm]
    aa$subtypes = bb$subtypes[mm]
    
    rm(bb)
    
    aa$clusters = aa$seurat_clusters
    
    saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_lognormamlized_pca_umap_manualAnnot.rds'))
    
    ##########################################
    # some cleaning for each cluster ??
    ##########################################
    
  }
  
  DimPlot(aa, label = TRUE, repel = TRUE) + ggtitle("snRNAseq (multiome)")
  p1 = FeaturePlot(aa, features = 'nFeature_RNA')
  p2 = FeaturePlot(aa, features = 'percent.mt')
  
  p1 | p2
  
  ggsave(filename = paste0(resDir, '/FeaturesPlot_nFeatures_pertMT.pdf'), width = 14, height = 6)
  
  # annotation for all cells
  Idents(aa) = aa$subtypes
  
  DimPlot(aa, label = TRUE, group.by = 'subtypes', repel = TRUE) + ggtitle("first manual annotation")
  ggsave(filename = paste0(resDir, '/first_test_umap_v2_manualAnnot.pdf'), width = 10, height = 8)
  
  aa$condition = factor(aa$condition, levels = paste0('Amex_scRNA_d', c(0, 1, 4, 7, 14)))
  aa$subtypes[which(is.na(aa$subtypes))] = 'unknow'
  
  DimPlot(aa, label = TRUE, repel = TRUE, split.by = 'condition') + NoLegend() +  
    ggtitle("first manual annotation")
  
  ggsave(filename = paste0(resDir, '/first_test_umap_v2_manualAnnot_byCondition.pdf'), width = 25, height = 8)
  
  # aa$subtypes = aa$celltypes
  # aa$clusters = aa$seurat_clusters
  # markers.coarse = readRDS(file = paste0(RdataDir_spliced, 'top10_markerGenes_coarseCluster.rds'))
  
  ##########################################
  # ## subclustering CMs
  ##########################################
  celltype.sels = 'CM'
  sub.obj = subset(aa, cells = colnames(aa)[!is.na(match(aa$celltypes, celltype.sels))])
  
  #sub.obj = SCTransform(aa, ncells = 3000, assay = "RNA", verbose = FALSE, 
  #            variable.features.n = 3000, return.only.var.genes = TRUE, vst.flavor = "v2")
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
  sub.obj = ScaleData(sub.obj)
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 30)
  
  nb.pcs = 20 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 20; min.dist = 0.05;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, 
                     n.neighbors = n.neighbors,
                     min.dist = min.dist)
  
  # sub.obj$clusters = sub.obj$seurat_clusters
  p1 = DimPlot(sub.obj, group.by = 'clusters', reduction = 'umap', label = TRUE, label.size = 5) +
    ggtitle(celltype.sels)
  
  features = rownames(aa)[grep('MYH6|ACTN2|NPPA|TNN|GATA4', rownames(aa))]
  FeaturePlot(aa, features = features, cols = c('gray', 'red'))
  
  sub.obj <- FindNeighbors(sub.obj, dims = 1:20)
  sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 0.7)
  
  p2 = DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 6) +
    ggtitle(paste0(celltype.sels, ' -- subclusters'))
  
  markers = FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3)
  
  markers %>%
    filter(!str_detect(gene, '^(AMEX|LOC|N/A)')) %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC) -> top10
  
  xx = subset(sub.obj, downsample = 500)
  #tops.common = c(markers.coarse$gene[!is.na(match(markers.coarse$cluster, c(1, 2, 3, 6, 13, 14, 15)))])
  
  p3 = DoHeatmap(xx, features = top10$gene) + NoLegend()
  p3
  #ggsave(filename = paste0(resDir, '/first_test_clusterMarkers_v2.pdf'), width = 10, height = 30)
  
  pdfname = paste0(resDir, '/subclustering_associatedMarkerGenes_', celltype.sels, '_v2.pdf')
  pdf(pdfname, width=16, height = 16)
  
  p1
  p2
  p3 
  dev.off()
  
  ## save the subclustering labels 
  cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 2)]
  mm = which(!is.na(match(colnames(aa), cell.sels)))
  cat(length(cell.sels), '--', length(mm), '\n')
  aa$subtypes[mm] = "CM.Endo.doublet"
  
  cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters != 2)]
  mm = which(!is.na(match(colnames(aa), cell.sels)))
  cat(length(cell.sels), '--', length(mm), '\n')
  
  
  ##########################################
  # subclustering FBs
  ##########################################
  celltype.sels = 'FB'
  sub.obj = subset(aa, cells = colnames(aa)[!is.na(match(aa$celltypes, celltype.sels))])
  
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
  
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 30)
  
  nb.pcs = 20 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 20; min.dist = 0.05;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, 
                     n.neighbors = n.neighbors,
                     min.dist = min.dist)
  
  
  # sub.obj$clusters = sub.obj$seurat_clusters
  p1 = DimPlot(sub.obj, group.by = 'clusters', reduction = 'umap', label = TRUE, label.size = 5) +
    ggtitle(celltype.sels)
  
  features = rownames(aa)[grep('VIM|COL1A2|FSTL1|POSTN', rownames(aa))]
  FeaturePlot(aa, features = features, cols = c('gray', 'red'))
  
  
  sub.obj <- FindNeighbors(sub.obj, dims = 1:20)
  sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 0.7)
  
  p2 = DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 5) +
    ggtitle(paste0(celltype.sels, ' -- subclusters'))
  
  p1 + p2
  
  markers = FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3)
  
  markers %>%
    filter(!str_detect(gene, '^(AMEX|LOC|N/A)')) %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC) -> top10
  
  xx = subset(sub.obj, downsample = 500)
  #tops.common = c(markers.coarse$gene[!is.na(match(markers.coarse$cluster, c(1, 2, 3, 6, 13, 14, 15)))])
  
  p3 = DoHeatmap(xx, features = top10$gene) + NoLegend()
  p3
  
  pdfname = paste0(resDir, '/subclustering_associatedMarkerGenes_', celltype.sels, '_v2.pdf')
  pdf(pdfname, width=18, height = 16)
  
  p1
  p2
  p3 
  dev.off()
  
  ## save the subclustering labels 
  cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 1)]
  mm = which(!is.na(match(colnames(aa), cell.sels)))
  cat(length(cell.sels), '--', length(mm), '\n')
  aa$subtypes[mm] = "FB/Endo"
  
  cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 8)]
  mm = which(!is.na(match(colnames(aa), cell.sels)))
  cat(length(cell.sels), '--', length(mm), '\n')
  aa$subtypes[mm] = "FB.CM.doublet"
  
  
  ##########################################
  # subclustering immune cells 
  ##########################################
  #celltype.sels = c('immune', 'blood', 'Macrophage', 'B')
  cluster.sels = c(24, 22, 18, 17, 9, 11, 20, 21, 13)
  sub.obj = subset(aa, cells = colnames(aa)[!is.na(match(aa$clusters, cluster.sels))])
  
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 5000)
  sub.obj = ScaleData(sub.obj)
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 30)
  
  nb.pcs = 20 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 20; min.dist = 0.1;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, 
                     n.neighbors = n.neighbors,
                     min.dist = min.dist)
  
  # sub.obj$clusters = sub.obj$seurat_clusters
  p1 = DimPlot(sub.obj, group.by = 'clusters', reduction = 'umap', label = TRUE, label.size = 5)
  p11 = DimPlot(sub.obj, group.by = 'subtypes', reduction = 'umap', label = TRUE, repel = TRUE, label.size = 5) + 
    NoLegend()
  
  p1
  p11
  
  features = rownames(sub.obj)[grep('PTPRC|CD68|CD8A|CD74|CSF1R|ITGAM', rownames(sub.obj))]
  FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))
  
  sub.obj <- FindNeighbors(sub.obj, dims = 1:20)
  sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 0.3)
  
  p2 = DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 5)
  
  p1 + p2
  
  markers = FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3)
  
  markers %>%
    filter(!str_detect(gene, '^(AMEX|LOC|N/A)')) %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC) -> top10
  
  xx = subset(sub.obj, downsample = 500)
  #tops.common = c(markers.coarse$gene[!is.na(match(markers.coarse$cluster, c(1, 2, 3, 6, 13, 14, 15)))])
  
  p3 = DoHeatmap(xx, features = top10$gene) + NoLegend()
  p3
  
  pdfname = paste0(resDir, '/subclustering_associatedMarkerGenes_', celltype.sels, '_v2.pdf')
  pdf(pdfname, width=18, height = 16)
  
  p1
  p11
  p2
  p3 
  dev.off()
  
  ## save the subclustering labels
  Save.subtype.Annot = FALSE
  if(Save.subtype.Annot){
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 0)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = "Macrophage"
    
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 1)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = 'blood'
    
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 2)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = 'Thrombocytes/Megakaryocytes'
    
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 3)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = 'Bcell'
    
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 4)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = 'notSure.immune.CM.doublet'
    
    
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 5)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = 'immune/MF'
    
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 6)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = 'blood/neutrophil'
    
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 7)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = 'Tcell'
    
  }
  
  ##########################################
  # subclustering Endo
  ##########################################
  celltype.sels = c('Endo')
  sub.obj = subset(aa, cells = colnames(aa)[!is.na(match(aa$celltypes, celltype.sels))])
  
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
  
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 30)
  
  nb.pcs = 20 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 20; min.dist = 0.05;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, 
                     n.neighbors = n.neighbors,
                     min.dist = min.dist)
  
  # sub.obj$clusters = sub.obj$seurat_clusters
  p1 = DimPlot(sub.obj, group.by = 'clusters', reduction = 'umap', label = TRUE, label.size = 5) +
    ggtitle(celltype.sels)
  
  #features = rownames(sub.obj)[grep('PTPRC|CD68|CD8A|CD74|CSF1R|ITGAM', rownames(sub.obj))]
  #FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))
  
  sub.obj <- FindNeighbors(sub.obj, dims = 1:10)
  sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 0.5)
  
  p2 = DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 5) +
    ggtitle(paste0(celltype.sels, ' -- subclusters'))
  
  p1 + p2
  
  markers = FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
  
  markers %>%
    filter(!str_detect(gene, '^(AMEX|LOC|N/A)')) %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC) -> top10
  
  xx = subset(sub.obj, downsample = 500)
  #tops.common = c(markers.coarse$gene[!is.na(match(markers.coarse$cluster, c(1, 2, 3, 6, 13, 14, 15)))])
  
  p3 = DoHeatmap(xx, features = top10$gene) + NoLegend()
  p3
  
  pdfname = paste0(resDir, '/subclustering_associatedMarkerGenes_', celltype.sels, '_v2.pdf')
  pdf(pdfname, width=18, height = 16)
  
  p1
  p2
  p3 
  dev.off()
  
  
}
