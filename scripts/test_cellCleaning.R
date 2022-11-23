##########################################################################
##########################################################################
# Project: heart regeneration
# Script purpose: to clean the cells (but keep interesting subpopulation) and 
# to identify the time-specific population 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Nov 21 16:54:43 2022
##########################################################################
##########################################################################
## here is the starting point where there are 48867 cells
aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                           '_lognormamlized_pca_umap_v3.rds'))

# Elad's 37609 cells with 42 subtypes
xx = readRDS(file = paste0('/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_subtypes_final_20221117.rds'))

missed = which(is.na(match(colnames(aa), colnames(xx))))

aa$keep = 'keep'
aa$keep[missed] = 'missed'
aa$keep = as.factor(aa$keep)

# transfer the subtype labels
aa$subtypes = NA
mm = match(colnames(xx), colnames(aa))
aa$subtypes[mm] = as.character(xx$subtypes)

rm(xx)

# aggregate subtypes into major cell types 
aa$celltypes = as.character(aa$subtypes)

aa$celltypes[grep('CM_|CMs_|_CM|_CM_', aa$subtypes)] = 'CM'
aa$celltypes[grep('EC_|_EC', aa$subtypes)] = 'EC'
aa$celltypes[grep('FB_', aa$subtypes)] = 'FB'
aa$celltypes[grep('B_cells', aa$subtypes)] = 'B'
aa$celltypes[grep('T_cells', aa$subtypes)] = 'T'

aa$celltypes[grep('Macrophages|_MF|Macs', aa$subtypes)] = 'Macs'
aa$celltypes[grep('Megakeryocytes', aa$subtypes)] = 'Megakeryocytes'
aa$celltypes[grep('RBC', aa$subtypes)] = 'RBC'
aa$celltypes[grep('Neu_', aa$subtypes)] = 'Neutrophile'


p1 = VlnPlot(aa, features = 'nFeature_RNA', group.by = 'keep',
        y.max = 10000) +
  geom_hline(yintercept = c(500, 2500, 3000))

p2 = VlnPlot(aa, features = 'nCount_RNA', group.by = 'keep', y.max = 50000)
p3 = VlnPlot(aa, features = 'percent.mt', group.by = 'keep', y.max = 100)

p1|p2|p3

ggsave(filename = paste0(resDir, '/compare_keep_discarded_cells.pdf'), width = 20, height = 8)


p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'keep') 
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes') 
p1 | p2

p3 = FeaturePlot(object = aa, features = 'nFeature_RNA')
p4 = FeaturePlot(aa, features = 'percent.mt') 
p3 | p4

(p1 | p2)/(p3 | p4)

ggsave(filename = paste0(resDir, '/compare_celltypes_nFeatures_pctMT_keep.vs.discarded_cells.pdf'), 
       width = 20, height = 16)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                          '_lognormamlized_pca_umap_keep.missed_subtypes.rds'))

########################################################
########################################################
# Section : we will keep the previous cell filtering (pct.MT < 60% & nFeatures.RNA > 200)  
# preditc doublets
########################################################
########################################################
library(DoubletFinder)
require(Seurat)

aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                        '_lognormamlized_pca_umap_keep.missed_subtypes.rds'))

cc = unique(aa$condition)

for(n in 1:length(cc))
{
  # n = 1]
  cat(n, ' -- ', cc[n], '\n')
  subs <- subset(aa, condition == cc[n])
  
  subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 5000)
  subs <- ScaleData(subs)
  
  subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = TRUE)
  
  subs <- FindNeighbors(subs, dims = 1:30)
  subs <- FindClusters(subs, resolution = 1)
  
  subs <- RunUMAP(subs, dims = 1:30)
  
  sweep.res.list_nsclc <- paramSweep_v3(subs)
  
  sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = FALSE)
  bcmvn_nsclc <- find.pK(sweep.stats_nsclc)
  
  pK <- bcmvn_nsclc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  
  pK <- as.numeric(as.character(pK[[1]]))
  annotations <- subs@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  
  nExp_poi <- round(0.076*nrow(subs@meta.data))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  subs <- doubletFinder_v3(subs, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj,  
                           reuse.pANN = FALSE, sct = FALSE)
  
  df_out = subs@meta.data
  subs$DF_out = df_out[, grep('DF.classification', colnames(df_out))]
  
  DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'DF_out',
          raster=FALSE)
  ggsave(filename = paste0(resDir, '/subs_doubletFinder_out_', cc[n], '.pdf'), 
         width = 12, height = 8)
  
  saveRDS(subs, file = paste0(RdataDir, 'subs_doubletFinder_out_', cc[n], '.rds'))
  
}


##########################################
# save the doubletFinder in the main table  
##########################################
cc = unique(aa$condition)
aa$DF_out = NA

for(n in 1:length(cc))
{
  # n = 1
  cat(n, '--', cc[n], '\n')
  subs = readRDS(file = paste0(RdataDir, 'subs_doubletFinder_out_', cc[n], '.rds'))
  aa$DF_out[match(colnames(subs), colnames(aa))] = subs$DF_out
  
}

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                          '_lognormamlized_pca_umap_keep.missed_subtypes_DFinderOut.rds'))

##########################################
# Visulize the doublet 
##########################################
# aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)
# aa <- ScaleData(aa)
# aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)
# 
# ElbowPlot(aa, ndims = 30)
# 
# Idents(aa) = aa$condition
# aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'DF_out', raster=FALSE)
ggsave(filename = paste0(resDir, '/umap_doubletFinder_results.pdf'), width = 12, height = 8)

VlnPlot(aa, features = 'nCount_RNA', group.by = 'DF_out', pt.size = 0.1) +
  geom_hline(yintercept = c(25000, 12500), col = 'red')

ggsave(filename = paste0(resDir, '/nCounts_RNA_double.vs.singlet_doubletFinder_results.pdf'), 
       width = 12, height = 8)

as_tibble(data.frame(condition = aa$condition, group= aa$DF_out)) %>%
  group_by(condition, group) %>% tally() 

pcts = c()
for(n in 1:length(cc))
{
  # n =1
  pcts = c(pcts, length(which(aa$DF_out== 'Doublet' & aa$condition == cc[n]))/length(which(aa$condition == cc[n])))
  
}

data.frame(condition = cc, pct = pcts) %>%
  ggplot(aes(x = condition, y = pct, fill = condition)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none")  + 
  ggtitle('pct of doublets by DF ') + 
  theme(axis.text.x = element_text(angle = 90)) 

ggsave(filename = paste0(resDir, '/Percentages_doublet.vs.total_doubletFinder_results.pdf'), 
       width = 12, height = 8)


aa = subset(aa, DF_out == 'Singlet')
saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                          '_lognormamlized_pca_umap_keep.missed_subtypes_DFinderFiltered.rds'))


## check the missed cells after doublet filtering
p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'keep') 
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes') 
p1 | p2

p3 = FeaturePlot(object = aa, features = 'nFeature_RNA')
p4 = FeaturePlot(aa, features = 'percent.mt') 
p3 | p4

(p1 | p2)/(p3 | p4)

ggsave(filename = paste0(resDir, '/DoubletFiltered_compare_celltypes_nFeatures_pctMT_keep.vs.discarded_cells.pdf'), 
       width = 20, height = 16)


aa$celltypes[is.na(aa$celltypes)] = 'missed'

aa$condition = factor(aa$condition, levels = c('Amex_scRNA_d0', 'Amex_scRNA_d1', 
                                               'Amex_scRNA_d4', 'Amex_scRNA_d7',
                                               'Amex_scRNA_d14'))

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', split.by = 'condition')


########################################################
########################################################
# Section : test if the time-specific subpopulations are robust to batch correction
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
              '_lognormamlized_pca_umap_keep.missed_subtypes_DFinderFiltered.rds'))

aa$condition = factor(aa$condition, levels = c('Amex_scRNA_d0', 'Amex_scRNA_d1', 
                                               'Amex_scRNA_d4', 'Amex_scRNA_d7',
                                               'Amex_scRNA_d14'))
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes')


## redo the umap and HVGs and clustering
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000)
all.genes <- rownames(aa)
aa <- ScaleData(aa, features = all.genes)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)

ElbowPlot(aa, ndims = 30)

aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 1.0)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

p1 = DimPlot(aa, label = TRUE, repel = TRUE) + ggtitle("doublet filtered - compare Elad's subtypes")
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes')

p1 + p2

ggsave(filename = paste0(resDir, '/snRNAseq_doubletFiltered_clusters_vs_celltypes.pdf'), 
       width = 20, height = 8)


# manaul 
DimPlot(aa, label = TRUE, repel = TRUE, cells.highlight = colnames(aa)[which(aa$seurat_clusters == '23')])

aa$labels = NA
cluster_sels = c(0, 1, 3, 4, 6, 8, 11, 14, 28, 31)
aa$labels[!is.na(match(aa$seurat_clusters, cluster_sels))] = 'EC'

cluster_sels = c(2, 13, 9, 15, 36, 22, 30, 27, 19, 21, 7, 17)
aa$labels[!is.na(match(aa$seurat_clusters, cluster_sels))] = 'CM'

cluster_sels = c(29, 20, 10)
aa$labels[!is.na(match(aa$seurat_clusters, cluster_sels))] = 'FB'

cluster_sels = c(12)
aa$labels[!is.na(match(aa$seurat_clusters, cluster_sels))] = 'Megakeryocytes'

cluster_sels = c(33, 5)
aa$labels[!is.na(match(aa$seurat_clusters, cluster_sels))] = 'Macs.Neurophil'


DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'labels')


features = rownames(aa)[grep('PECAM|TEK-|VIM|ACTA2', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('COL1A1|COL1A2|COL3A1|LUM-|POSTN', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

##########################################
# subset FB 
##########################################
sub.obj = subset(aa, labels == 'FB')

DimPlot(sub.obj, label = TRUE, repel = TRUE, group.by = 'labels')

sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
sub.obj = ScaleData(sub.obj)
sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
ElbowPlot(sub.obj, ndims = 30)

nb.pcs = 30 # nb of pcs depends on the considered clusters or ids
n.neighbors = 30; min.dist = 0.3;
sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, 
                   n.neighbors = n.neighbors,
                   min.dist = min.dist)

sub.obj <- FindNeighbors(sub.obj, dims = 1:30)
sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 0.5)


# sub.obj$clusters = sub.obj$seurat_clusters
p1 = DimPlot(sub.obj, group.by = 'clusters', reduction = 'umap', label = TRUE, label.size = 5) +
  ggtitle(celltype.sels)

p2 = DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 6) 

DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 6, group.by = 'keep') 

DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 4, group.by = 'keep', split.by = 'condition') 

DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 4,  split.by = 'condition') 


##########################################
# test the batch correction (RPCA)
##########################################
sub.list <- SplitObject(sub.obj, split.by = "condition")

# normalize and identify variable features for each dataset independently
sub.list <- lapply(X = sub.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = sub.list, nfeatures = 3000)
features.common = rownames(sub.obj)
sub.list <- lapply(X = sub.list, FUN = function(x) {
  x <- ScaleData(x, features = features.common, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  
})

sub.anchors <- FindIntegrationAnchors(object.list = sub.list, 
                                      anchor.features = features, 
                                      reduction = "rpca", 
                                      k.anchor = 5)

rm(sub.list)

# this command creates an 'integrated' data assay
sub.combined <- IntegrateData(anchorset = sub.anchors, features.to.integrate = features.common)

rm(sub.anchors)
# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(sub.combined) <- "integrated"

xx = DietSeurat(sub.combined, counts = TRUE, data = TRUE, scale.data = FALSE, assays = 'integrated')
xx@assays$integrated@counts = sub.combined@assays$RNA@counts

#saveRDS(xx, file = paste0(RdataDir, 
#'Seurat.obj_adultMiceHeart_Forte2020_Ren2020_subCombined_logNormalize_counts_v3.rds'))


# Run the standard workflow for visualization and clustering
sub.combined = xx;
rm(xx)
sub.combined$condition = factor(sub.combined$condition, 
                                levels = c('Amex_scRNA_d0', 'Amex_scRNA_d1', 
                                'Amex_scRNA_d4', 'Amex_scRNA_d7',
                                'Amex_scRNA_d14'))
sub.combined <- ScaleData(sub.combined, verbose = FALSE)
sub.combined <- RunPCA(sub.combined, npcs = 30, verbose = FALSE)

ElbowPlot(sub.combined, ndims = 30)

sub.combined <- FindNeighbors(sub.combined, reduction = "pca", dims = 1:20)
sub.combined <- FindClusters(sub.combined, algorithm = 3, resolution = 0.5)

sub.combined <- RunUMAP(sub.combined, reduction = "pca", dims = 1:30, n.neighbors = 30, min.dist = 0.3) 


p0 = DimPlot(sub.obj, reduction = "umap")
p1 = DimPlot(sub.combined, reduction = "umap", group.by = 'RNA_snn_res.0.5')

p0 | p1

DimPlot(sub.combined, reduction = 'umap', label = TRUE, label.size = 4,  split.by = 'condition') 

ggsave(paste0(resDir, '/Forte2020_Ren2020_IntegrationRPCA_', Normalization, '.pdf'), 
       width = 24, height = 10)


