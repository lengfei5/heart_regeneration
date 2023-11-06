##########################################################################
##########################################################################
# Project: Heart regeneration 
# Script purpose: process and analyze the single cell/nucleus RNA-seq data for the cell type reference
# for adult heart there are tww references: cardiomyocyte and non-cardiomycyte. 
# because the ST-specific deconvolution method requires UMI counts as input, at the end the deconvolution has to be done 
# separately with two reference and integration post-hoc 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Nov 11 14:14:12 2021
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_reference_20211111'

resDir = paste0("../results/scRNAseq_adultMouse", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

# required libraries
library(data.table)
require(Seurat)
library(SeuratObject)
require(sctransform)
require(ggplot2)
library(dplyr)
library(patchwork)
require(tictoc)
library(pryr) # monitor the memory usage
mem_used()

Normalization = 'lognormal' # ('lognormal or SCT')

########################################################
########################################################
# Section I : First process the adult cardiomyocyte from Ren et al., 2020
# 
########################################################
########################################################
Process.Cardiomyocyte.Ren.2020 = FALSE
if(Process.Cardiomyocyte.Ren.2020){
  dataDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/published_dataset/Ren_2020'
  metadata = read.delim(file = paste0(dataDir, '/GSE120064_TAC_clean_cell_info_summary.txt'), sep = '\t', header = TRUE)
  counts = fread(file = paste0(dataDir, '/GSE120064_TAC_raw_umi_matrix.csv'), header = TRUE, nThread = 6)
  #rownames(counts) = counts$V1
  
  # align the cellID in metadata and in the count table
  mm = match(metadata$CellID, colnames(counts))
  counts = counts[, c(1, mm), with=FALSE]
  
  # keep only the week 0 data
  jj = which(metadata$condition == '0w' | metadata$condition == '2w')
  
  metadata = metadata[jj, ]
  counts = counts[, c(1, jj+1), with = FALSE]
  counts = as.data.frame(counts)
  rownames(counts) = counts$V1
  counts = counts[, -1]
  
  # remove genes with 0 umi counts
  ss = apply(as.matrix(counts), 1, sum)
  counts = counts[which(ss>0), ]
  
  metadata = data.frame(metadata)
  rownames(metadata) = metadata$CellID
  
  save(metadata, counts, file = paste0(RdataDir, 'metadata_counts_week0.week2.Rdata'))
  
}

##########################################
# make Seurat object with metadata and counts  
##########################################
load(file = paste0(RdataDir, 'metadata_counts_week0.week2.Rdata'))

myo <- CreateSeuratObject(counts = counts, project = "adult", min.cells = 50, min.features = 500)
myo = AddMetaData(myo, metadata, col.name = NULL) 

myo[["percent.mt"]] <- PercentageFeatureSet(myo, pattern = "^mt-")

plot1 <- FeatureScatter(myo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(myo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

aa = subset(myo, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA < 5*10^5)
rm(myo)
rm(counts)

if(SCT.normalization){
  aa <- SCTransform(aa, assay = "RNA", verbose = FALSE, variable.features.n = 3000, return.only.var.genes = FALSE, 
                    min_cells=5) 
  
}else{
  aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
  
  plot1 <- VariableFeaturePlot(aa)
  
  top10 <- head(VariableFeatures(aa), 10)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2
  
  all.genes <- rownames(aa)
  aa <- ScaleData(aa, features = all.genes)
  
}

aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 30)

aa <- FindNeighbors(aa, dims = 1:10)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)

### KEEP the umap configuration for only week0: nfeatures = 2000, dims = 1:20, n.neighbors = 30 and min.dist = 0.05 (sctransform)
### umap configuraiton for week0 and week2: nfeatures = 3000, dims = 1:20, n.neighbors = 30/50, min.dist = 0.1 (sctransform)
### umap configuration for week0 and week2: nfeatures = 3000, dims = 1:20, n.neighbors = 20, min.dist = 0.05 (seurat.norm)

if(Normalization == 'SCT'){
  aa = RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
  #saveRDS(aa, file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_SCT_umap.rds'))
}else{
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 20, min.dist = 0.05)
  #saveRDS(aa, file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization_umap.rds'))
}

##########################################
# double check the cardiomyocyte subtypes  
##########################################
Double.check.adult.cardiomyocyte.major.celltypes.subtypes(aa)

########################################################
########################################################
# Section II : adult non-caridomyocyte single cell dataset 
# the dataset is from Forte et al. 2020 
# processed dataset from the collaborate Shoval
# in which the lognormal normalization was used and no batch correction or data integration were performed
# because they think it is not required for this dataset, or the batch difference is minor
########################################################
########################################################
dataDir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/data/'

aa = readRDS(file = paste0(dataDir, 'Forte_et_al_forEladJingkui.rds'))

# add some more information to the metadata
aa$timepoints = NA
aa$timepoints[which(aa$orig.ident == 'MF17010')] = 'd0'
aa$timepoints[which(aa$orig.ident == 'MF17013')] = 'd1'
aa$timepoints[which(aa$orig.ident == 'MF17014')] = 'd3'
aa$timepoints[which(aa$orig.ident == 'MF17015')] = 'd5'
aa$timepoints[which(aa$orig.ident == 'MF17016')] = 'd7'
aa$timepoints[which(aa$orig.ident == 'MF17017')] = 'd14'
aa$timepoints[which(aa$orig.ident == 'MF17018')] = 'd28'

meta = aa@meta.data

p1 = DimPlot(aa, reduction = 'umap_0.05', group.by = 'my_annot')
p2 = DimPlot(aa, reduction = "umap_0.05", group.by = c("timepoints"))

p1 + p2 
ggsave(paste0(resDir, '/Forte2020_Umap_clusters._cellType.original.pdf'), 
                 width = 16, height = 8)

##########################################
# explore the normalization and intergration 
##########################################
if(Normalization == 'SCT'){
  
  tic()
  aa <- SCTransform(aa, assay = "RNA", verbose = FALSE, variable.features.n = 3000, 
                    return.only.var.genes = FALSE, ncells = 5000, conserve.memory = FALSE,
                    min_cells=5) 
  toc()
  
  saveRDS(aa, file = paste0(RdataDir, 'Forte2020_SCTnorm_allgenes.rds'))
  
}else{
  
  aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
  
  plot1 <- VariableFeaturePlot(aa)
  
  top10 <- head(VariableFeatures(aa), 10)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2
  
  all.genes <- rownames(aa)
  aa <- ScaleData(aa, features = all.genes)
  
  saveRDS(aa, file = paste0(RdataDir, 'Forte2020_logNormalize_allgenes.rds'))
  
} 

aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 30)

aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 50, min.dist = 0.05)

p1 = DimPlot(aa, reduction = 'umap', group.by = 'my_annot') + ggtitle(paste0(Normalization, ' umap'))


p0 = DimPlot(aa, reduction = 'umap_0.05', group.by = 'my_annot') + ggtitle('original umap')
p0 + p1 
ggsave(paste0(resDir, '/Forte2020_myumap_cellType.Shoval_', Normalization, '.pdf'), 
               width = 16, height = 8)

p2 = DimPlot(aa, reduction = "umap", group.by = c("timepoints"))

p1 + p2 
ggsave(paste0(resDir, '/Forte2020_Umap_myclusters_cellType.Shoval_', Normalization, '.pdf'), 
                 width = 16, height = 8)

##########################################
# double check Shoval's cluster annotation
##########################################
Double.check.adult.non.cardiomyocyte.major.celltypes.subtypes(aa)


########################################################
########################################################
# Section II: # integrate Ren2020 and Forte2020 to have one reference using SCTransform and RPCA from Seurat
# original code from https://satijalab.org/seurat/articles/integration_rpca.html
# 
########################################################
########################################################
Merge.adult.mice.cardiomyocyte.noncardiomyocyte = FALSE

if(Merge.adult.mice.cardiomyocyte.noncardiomyocyte){
  
  aa = readRDS(file = paste0(RdataDir, 'Forte2020_logNormalize_allgenes_majorCellTypes_subtypes.rds'))
  cms = readRDS(file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization_umap_subtypes.rds'))
  
  aa$dataset = 'Forte2020'
  cms$dataset = 'Ren2020'
  #aa$annot.ref = aa$my_annot
  #cms$annot.ref = cms$CellType
  features.common = intersect(rownames(aa), rownames(cms))
  
  refs.merged = merge(aa, y = cms, add.cell.ids = c("Forte2020", "Ren2020"), project = "adultHeart")
  
  ref.list <- SplitObject(refs.merged, split.by = "dataset")
  
  rm(list = c('aa', 'cms', 'refs.merged')) # remove big seurat objects to clear memory
  
  # normalize and identify variable features for each dataset independently
  ref.list <- lapply(X = ref.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  })
  
  # select features that are repeatedly variable across datasets for integration run PCA on each
  # dataset using these features
  features <- SelectIntegrationFeatures(object.list = ref.list, nfeatures = 3000)
  
  ref.list <- lapply(X = ref.list, FUN = function(x) {
    x <- ScaleData(x, features = features.common, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
    
  })
  
  ref.anchors <- FindIntegrationAnchors(object.list = ref.list, anchor.features = features, reduction = "rpca", 
                                        k.anchor = 5)
  
  rm(ref.list)
  
  # this command creates an 'integrated' data assay
  ref.combined <- IntegrateData(anchorset = ref.anchors, features.to.integrate = features.common)
  
  rm(ref.anchors)
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(ref.combined) <- "integrated"
  
  xx = DietSeurat(ref.combined, counts = TRUE, data = TRUE, scale.data = FALSE, assays = 'integrated')
  xx@assays$integrated@counts = ref.combined@assays$RNA@counts
  
  saveRDS(xx, file = paste0(RdataDir, 'Seurat.obj_adultMiceHeart_Forte2020_Ren2020_refCombined_logNormalize_counts_v3.rds'))
  
  
  # Run the standard workflow for visualization and clustering
  ref.combined = readRDS(file =paste0(RdataDir, 
                                      'Seurat.obj_adultMiceHeart_Forte2020_Ren2020_refCombined_logNormalize_counts_v3.rds'))
  
  ref.combined <- ScaleData(ref.combined, verbose = FALSE)
  ref.combined <- RunPCA(ref.combined, npcs = 30, verbose = FALSE)
  
  ElbowPlot(ref.combined, ndims = 30)
  
  ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:20)
  ref.combined <- FindClusters(ref.combined, resolution = 0.2)
  
  ref.combined <- RunUMAP(ref.combined, reduction = "pca", dims = 1:30, n.neighbors = 50, min.dist = 0.05) 
  
  DimPlot(ref.combined, reduction = "umap")
  
  kk = which(ref.combined$dataset == 'Ren2020' & ref.combined$celltype != 'CM')
  ref.combined$celltype[kk] = paste0(ref.combined$celltype[kk], '_Ren2020')
  
  # Visualization
  p1 <- DimPlot(ref.combined, reduction = "umap", group.by = "dataset")
  p2 <- DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE,
                repel = TRUE)
  p1 + p2 
  
  ggsave(paste0(resDir, '/Forte2020_Ren2020_IntegrationRPCA_', Normalization, '.pdf'), 
                   width = 24, height = 10)
  
  
  ##########################################
  # clean the reference, i.e. remove the non-cardiomyocyte from Ren2020
  # change the confusing annotation names from Shoval
  ##########################################
  kk = which(ref.combined$dataset == 'Forte2020'| (ref.combined$dataset == 'Ren2020' & ref.combined$celltype == 'CM'))
  refs = ref.combined[,kk]
  
  p1 <- DimPlot(refs, reduction = "umap", group.by = "dataset")
  p2 <- DimPlot(refs, reduction = "umap", group.by = "celltype", label = TRUE,
                repel = TRUE)
  p1 + p2 
  
  ggsave(paste0(resDir, '/Forte2020_Ren2020onlyCM_IntegrationRPCA_', Normalization, '.pdf'), 
         width = 24, height = 10)
  
  p2 
  ggsave(paste0(resDir, '/Forte2020_Ren2020onlyCM_IntegrationRPCA_celltypes.in.Refs_', Normalization, '.pdf'), 
         width = 12, height = 10)
  
  DimPlot(refs, reduction = "umap", group.by = "subtype", label = TRUE,
          repel = TRUE)
  
  ggsave(paste0(resDir, '/Forte2020_Ren2020onlyCM_IntegrationRPCA_subtypes.in.Refs_', Normalization, '.pdf'), 
         width = 12, height = 10)
  
  # xx = subset(refs, cells = colnames(refs)[which(refs$dataset == 'Ren2020')])
  # 
  # xx <- ScaleData(xx, verbose = FALSE)
  # xx <- RunPCA(xx, npcs = 30, verbose = FALSE)
  # 
  # ElbowPlot(xx, ndims = 30)
  # 
  # xx <- FindNeighbors(xx, reduction = "pca", dims = 1:20)
  # xx <- FindClusters(xx, resolution = 0.2)
  # 
  # xx <- RunUMAP(xx, reduction = "pca", dims = 1:10, n.neighbors = 10, min.dist = 0.05) 
  # 
  # DimPlot(xx, reduction = "umap")
  # 
  
  DimPlot(refs, reduction = "umap", group.by = "celltype", split.by = 'dataset',  label = TRUE, repel = TRUE)
  # DimPlot(refs, reduction = 'umap', group.by = 'integrated_snn_res.0.5')
  
  
  
  rm(ref.combined)
  
  saveRDS(refs, file = paste0(RdataDir, 
                              'Seurat.obj_adultMiceHeart_Forte2020.nonCM_Ren2020CM_refCombined_',
                              'cleanAnnot_logNormalize_v4.rds'))
  saveRDS(refs, file = paste0(RdataDir, 
                              'SeuratObj_adultMiceHeart_refCombine_Forte2020.nonCM_Ren2020CM_',
                              'cleanAnnot_logNormalize_v4.rds'))
  
}



########################################################
########################################################
# Section III: post-integration: double check the celltypes  
# 
########################################################
aa = readRDS(file = paste0(RdataDir, 
                           'Seurat.obj_adultMiceHeart_Forte2020.nonCM_Ren2020CM_refCombined_',
                           'cleanAnnot_logNormalize_v4.rds'))

##########################################
# prepare reference data and double check the main cell types and subtypes 
##########################################
Double.check.adult.cardiomyocyte.major.celltypes.subtypes = FALSE
if(Double.check.adult.cardiomyocyte.major.celltypes.subtypes){
  
  p1 = DimPlot(aa, reduction = 'umap', group.by = 'seurat_clusters')
  p2 = DimPlot(aa, reduction = "umap", group.by = c("CellType"))
  p1 + p2 
  
  p3 = DimPlot(aa, reduction = 'umap', group.by = 'condition')
  p4 = DimPlot(aa, reduction = 'umap', group.by = 'sample')
  (p1 + p2)/(p3 + p4)
  ggsave(paste0(resDir, '/Umap_Ren2020_week0.week2_newClusters_vs_cellType.original_conditions_samples_seuratNorm.pdf'), 
         width = 12, height = 8)
  
  # double check the cardiomyocyte markers from Elad 
  gg.examples =  c('Nppa', 'Nppb', 'Myh6', 'Tnnc1', 'Tnni3', 
                   'Tnnt2', 'Actn2', 'Gata4', 'Nkx2-5', 'Gja1', 
                   'Myl2', 'Tpm1', 'Ryr2', 'Atp2a2', 'Acta1')
  
  p5 = FeaturePlot(aa, reduction = 'umap', features = gg.examples)
  
  p2 / p5 + ggsave(paste0(resDir, '/Umap_cellType.original_FeatuerPlot_cardiomyocyoteMarkers_seuratNorm.pdf'), 
                   width = 14, height = 20)
  
  p6 = VlnPlot(aa, features = gg.examples, group.by = 'CellType')
  
  p6 + ggsave(paste0(resDir, '/Umap_VlnPlot_cardiomyocyoteMarkers_seuratNorm.pdf'), 
              width = 14, height = 10)
  
  p1 + p2 + ggsave(paste0(resDir, '/Umap_newClusters_vs_cellType.original_seuratNorm.pdf'), 
                   width = 12, height = 8)
  
  #DimPlot(aa, reduction = "umap", group.by = c("SubCluster"))
  
  FeaturePlot(aa, reduction = 'umap', features = c('Csf1r', 'Cd163'))
  FeaturePlot(aa, reduction = 'umap', features = c('S100a3', 'S100a9'))
  
  ##########################################
  # double check the subtypes of CM 
  ##########################################
  aa$celltype = aa$CellType
  aa$subtype = aa$SubCluster
  
  mcells = 'CM'
  ax = subset(aa, cells = colnames(aa)[which(aa$CellType == mcells)])
  table(ax$celltype)
  table(ax$subtype)
  
  ax <- FindVariableFeatures(ax, selection.method = "vst", nfeatures = 5000)
  ax <- ScaleData(ax, features = rownames(ax))
  
  ax <- RunPCA(ax, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(ax, ndims = 30)
  
  # UMAP to visualize subtypes
  ax <- RunUMAP(ax, dims = 1:10, n.neighbors = 20, min.dist = 0.05, n_threads = 6)
  
  DimPlot(ax, reduction = 'umap', group.by = 'subtype') + 
    ggtitle(paste0(mcells, '-', ' cells UMAP (', Normalization, ' nfeature = 5000, ndim=10, neighbors=20, mdist=0.05)'))
  
  ggsave(paste0(resDir, '/Ref_Ren2020_UMAP_', mcells, '_subcelltypes.pdf'), 
         width = 10, height = 8)
  
  ax <- FindNeighbors(ax, dims = 1:10)
  
  ax <- FindClusters(ax, verbose = FALSE, algorithm = 3, resolution = 0.5)
  p1 = DimPlot(ax, reduction = 'umap', group.by = 'seurat_clusters')
  p0 = DimPlot(ax, reduction = 'umap', group.by = 'subtype')
  
  p0 + p1
  
  ax$subtype = paste0('CM', ax$seurat_clusters)
  Idents(ax) = ax$subtype
  
  FeaturePlot(ax, features = c('Sln', 'Myl4', 'Myl7', 'Nppa', 'Myl1', # Atrial CM markers
                               'Fabp3', 'Pln', 'Myl2', 'Myl3', 'Pin', 'Mb' # Ventricular CMs
  ))
  
  
  VlnPlot(ax, features = c('Sln', 'Myl4', 'Myl7', 'Nppa', 'Myl1', # Atrial CM markers
                           'Fabp3', 'Pln', 'Myl2', 'Myl3', 'Pin', 'Mb' # Ventricular CMs
  ))
  
  ggsave(paste0(resDir, '/VlnPLot_Atrial_Ventricular_markerGenes_', mcells, '_subtypes.pdf'), width = 12, height = 10)
  
  ax.markers <- FindAllMarkers(ax, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.3)
  # saveRDS(ax.markers, file = paste0(RdataDir, 'Forte2020_logNormalize_allgenes_majorCellTypes_markerGenes.rds'))
  
  ax.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top10
  
  DoHeatmap(ax, features = top10$gene)
  
  ggsave(paste0(resDir, '/heatmap_markerGenes_', mcells, '_subtypes.pdf'), width = 12, height = 26)
  
  aa$subtype[match(colnames(ax), colnames(aa))] = ax$subtype
  saveRDS(aa, file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization_umap_subtypes.rds'))
  
  rm(aa)
  rm(ax)
  
}

Double.check.adult.non.cardiomyocyte.major.celltypes.subtypes = FALSE
if(Double.check.adult.cardiomyocyte.major.celltypes.subtypes){
  
  # modify the cell annotations
  aa$celltype = as.character(aa$my_annot)
  
  # major cell type FB
  aa$celltype[which(aa$celltype == '0 - Fibro-I')] = 'FB1'
  aa$celltype[which(aa$celltype == '9 - Fibro-II')] = 'FB2'
  aa$celltype[which(aa$celltype == '11 - Fibro-III')] = 'FB3'
  aa$celltype[which(aa$celltype == '15 - Fibro-IV')] = 'FB4'
  aa$celltype[which(aa$celltype == '3 - MyoF')] = 'myoFB'
  
  # major cell type EC
  aa$celltype[which(aa$celltype == '2 - EC-I')] = 'EC1'
  aa$celltype[which(aa$celltype == '18 - EC-II')] = 'EC2'
  aa$celltype[which(aa$celltype == '19 - EC-III')] = 'EC3'
  aa$celltype[which(aa$celltype == '17 - Lymph-EC')] = 'Lymph.EC'
  
  # major cell type smooth muscle cells
  aa$celltype[which(aa$celltype == '16 - SMC')] = 'SMC'
  
  # immune cells
  aa$celltype[which(aa$celltype == '6 - B cell')] = 'B'
  aa$celltype[which(aa$celltype == '4 - GRN')] = 'GN'
  aa$celltype[which(aa$celltype == '12 - NK/T')] = 'NK.T'
  aa$celltype[which(aa$celltype == '13 - Monon/DC')] = 'MCT.DC'
  aa$celltype[which(aa$celltype == '5 - Chil3 Mono')] = 'MCT.Chil3'
  
  aa$celltype[which(aa$celltype == '1 - Trem2 Macs')] = 'Mphage.Trem2'
  aa$celltype[which(aa$celltype == '10 - Arg1 Macs')] = 'Mphage.Argl'
  aa$celltype[which(aa$celltype == '7 - MHCII Macs')] = 'MHCII.Mphage'
  aa$celltype[which(aa$celltype == '14 - Proliferating Macs')] = 'prolife.Mphage'
  
  # merge major cell types
  aa$subtype = aa$celltype
  
  aa$celltype[grep('FB', aa$subtype)]  = 'FB'
  aa$celltype[grep('EC1|EC2|EC3|Lymph.EC', aa$subtype)]  = 'EC'
  aa$celltype[grep('SMC', aa$subtype)] = 'SMC'
  
  aa$celltype[grep('GN|MCT|Mphage|NK', aa$subtype)]  = 'immune'
  aa$celltype[which(aa$subtype == 'B')] = 'immune'
  
  p0 = DimPlot(aa, reduction = 'umap', group.by = 'celltype') + ggtitle('Shoval UMAP')
  
  p1 = DimPlot(aa, reduction = 'umap', group.by = 'celltype') + ggtitle(paste0(Normalization, ' - Elad umap'))
  
  p0 + p1
  
  ggsave(paste0(resDir, '/Ref_Forte2020_majorCelltypes_Shoval.umap_vs_Elad.umap.', Normalization, '.pdf'), 
         width = 18, height = 8)
  
  
  for(mcells in c('all', 'FB', 'EC', 'immune'))
  {
    # mcells = 'immune.others'
    
    if(mcells == 'all'){
      ax = aa
      
      Idents(ax) = ax$celltype
      
    }else{
      ax = subset(aa, cells = colnames(aa)[which(aa$celltype == mcells & aa$subtype != 'B' & aa$subtype != 'GN' &
                                                   aa$subtype != 'NK.T')] )
      table(ax$celltype)
      table(ax$subtype)
      
      ax <- FindVariableFeatures(ax, selection.method = "vst", nfeatures = 3000)
      ax <- ScaleData(ax, features = rownames(ax))
      
      ax <- RunPCA(ax, verbose = FALSE, weight.by.var = TRUE)
      ElbowPlot(ax, ndims = 30)
      
      # UMAP to visualize subtypes
      ax <- RunUMAP(ax, dims = 1:30, n.neighbors = 30, min.dist = 0.05, n_threads = 6)
      
      DimPlot(ax, reduction = 'umap', group.by = 'subtype') + 
        ggtitle(paste0(mcells, '-', ' cells UMAP (', Normalization, ' nfeature = 3000, ndim=30, neighbors=30, mdist=0.05)'))
      
      ggsave(paste0(resDir, '/Ref_Forte2020_UMAP_', mcells, '_excluding B GN NK.T_subcelltypes.pdf'), 
             width = 10, height = 8)
      
      ax <- FindNeighbors(ax, dims = 1:20)
      
      ax <- FindClusters(ax, verbose = FALSE, algorithm = 3, resolution = 0.4)
      p1 = DimPlot(ax, reduction = 'umap', group.by = 'seurat_clusters')
      p0 = DimPlot(ax, reduction = 'umap', group.by = 'subtype')
      
      p0 + p1
      # marker genes and heatmaps
      Idents(ax) = ax$subtype
      
      gg.examples = c('Col1a2', 'Vim', 'Fstl1', 'Ddr2', 'Acta2', 'Postn', 'Tcf21', 'Pdgfra', 'Col3a1', 'Col1a1', 
                      'Gsn', 'Fbln2', 'Sparc', 'Mmp2', 'Msln', 'Rspo1', 'Lum', 'Col8a1')
      
      gg.examples = c('Tek', 'Pecam1', 'Emcn', 'Cdh5', 'Kdr', 'Vwf', 'Fabp4', 'Tie1', 
                      'Flt1', 'Epas1', 'Ednrb', 'Gpihbp1', 'Npr3')
      VlnPlot(ax, features = gg.examples)
      ggsave(paste0(resDir, '/VlnPlot_markerGenes_', mcells, '_subtypes.pdf'), width = 20, height = 16)
      
      
    }
    
    ax.markers <- FindAllMarkers(ax, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
    # saveRDS(ax.markers, file = paste0(RdataDir, 'Forte2020_logNormalize_allgenes_majorCellTypes_markerGenes.rds'))
    
    ax.markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC) -> top10
    
    DoHeatmap(ax, features = top10$gene)
    
    ggsave(paste0(resDir, '/heatmap_markerGenes_', mcells, '_subtypes.pdf'), width = 12, height = 20)
    
  }
  
  ##########################################
  # add distinct immune cells as major cell types  
  ##########################################
  jj = which(aa$subtype == 'B')
  aa$celltype[jj] = aa$subtype[jj] 
  
  jj = which(aa$subtype == 'GN')
  aa$celltype[jj] = aa$subtype[jj] 
  
  jj = which(aa$subtype == 'NK.T')
  aa$celltype[jj] = aa$subtype[jj] 
  
  jj = which(aa$subtype == 'prolife.Mphage')
  aa$celltype[jj] = aa$subtype[jj] 
  
  jj = which(aa$subtype == 'MHCII.Mphage')
  aa$celltype[jj] = aa$subtype[jj] 
  
  jj = which(aa$subtype == 'MCT.DC')
  aa$celltype[jj] = 'immune.others'
  
  jj = which(aa$celltype == 'immune')
  aa$celltype[jj] = 'immune.others'
  
  p0 = DimPlot(aa, reduction = 'umap', group.by = 'celltype')
  p1 = DimPlot(aa, reduction = 'umap', group.by = 'subtype')
  
  p0 + p1
  ggsave(paste0(resDir, '/heatmap_markerGenes_', mcells, '_majoyCelltypes_subtypes.pdf'), width = 12, height = 20)
  
  rm(ax)
  
  xx = DietSeurat(aa, counts = TRUE, data = TRUE, scale.data = FALSE, assays = 'RNA', dimreducs = 'umap')
  rm(aa)
  saveRDS(xx, file = paste0(RdataDir, 'Forte2020_logNormalize_allgenes_majorCellTypes_subtypes.rds'))
  
}

##########################################
# actually used the reference 
##########################################
refs = readRDS(file = paste0('../results/Rdata/', 
                             'Seurat.obj_adultMiceHeart_Forte2020.nonCM_Ren2020CM_refCombined_cleanAnnot_',
                             'logNormalize_v4.rds'))

refs$celltype[which(refs$celltype == 'immune.others')] = 'Mphage.MCT'
refs = subset(refs, cells = colnames(refs)[which(refs$celltype != 'SMC')])

DimPlot(refs, reduction = 'umap', group.by = 'celltype')
DimPlot(refs, reduction = 'umap', group.by = 'subtype')

saveRDS(refs, file = paste0('data_examples/ref_scRNAseq.rds'))

