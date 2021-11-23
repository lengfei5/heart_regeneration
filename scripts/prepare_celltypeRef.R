##########################################################################
##########################################################################
# Project: Heart regeneration 
# Script purpose: prepare the single cell/nucleus RNA-seq data for the cell type reference
# for adult heart there are tww references: cardiomyocyte and non-cardiomycyte. 
# because the ST-specific deconvolution method requires UMI counts as input, at the end the deconvolution has to be done 
# separately with two reference and integration post-hoc 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Nov 11 14:14:12 2021
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_scRNAseq_reference_20211111'

resDir = paste0("../results/scRNAseq_mouse", version.analysis)
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

Normalization = 'lognormal' # ('lognormal or SCT')

########################################################
########################################################
# Section : First process the adult cardiomyocyte from Ren et al., 2020
# 
########################################################
########################################################
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

myo = subset(myo, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA < 5*10^5)

aa = myo

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

saveRDS(aa, file = paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization.rds'))

### KEEP the umap configuration for only week0: nfeatures = 2000, dims = 1:20, n.neighbors = 30 and min.dist = 0.05 (sctransform)
### umap configuraiton for week0 and week2: nfeatures = 3000, dims = 1:20, n.neighbors = 30/50, min.dist = 0.1 (sctransform)
### umap configuration for week0 and week2: nfeatures = 3000, dims = 1:20, n.neighbors = 20, min.dist = 0.05 (seurat.norm)
aa = readRDS(file = paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization.rds'))

if(Normalization == 'SCT'){
  aa = RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
  #saveRDS(aa, file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_SCT_umap.rds'))
}else{
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 20, min.dist = 0.05)
  #saveRDS(aa, file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization_umap.rds'))
}

p1 = DimPlot(aa, reduction = 'umap', group.by = 'seurat_clusters')
p2 = DimPlot(aa, reduction = "umap", group.by = c("CellType"))
p1 + p2 

p3 = DimPlot(aa, reduction = 'umap', group.by = 'condition')
p4 = DimPlot(aa, reduction = 'umap', group.by = 'sample')
(p1 + p2)/(p3 + p4) + 
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
# define the subclusters of CM (to do)
##########################################


##########################################
# identify marker genes and heatmap 
##########################################
Refine.clustering.annotaiton = FALSE
if(Refine.clustering.annotaiton){
  # remove the mitochonio marker genes
  #jj = grep('^mt-', rownames(aa)) 
  #aa = aa[-jj, ]
  
  # remove non-annotated cells
  aa = aa[ ,!is.na(aa$CellType)]
  Idents(aa) = aa$CellType
  
  aa.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.3, logfc.threshold = 0.5)
  
  aa.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
  DoHeatmap(aa, features = top10$gene) + NoLegend() + ggsave(paste0(resDir, '/heatmap_markerGenes_rmMt.pdf'), 
                                                             width = 12, height = 16)
    
}

########################################################
########################################################
# Section : adult non-caridomyocyte single cell dataset 
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

p1 + p2 + ggsave(paste0(resDir, '/Forte2020_Umap_clusters._cellType.original.pdf'), 
                 width = 16, height = 8)

saveRDS(aa, file = paste0(RdataDir, 'Forte2020_processedSeurat.obj.from.Shoval.rds'))

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
p0 + p1 +  ggsave(paste0(resDir, '/Forte2020_myumap_cellType.Shoval_', Normalization, '.pdf'), 
               width = 16, height = 8)

p2 = DimPlot(aa, reduction = "umap", group.by = c("timepoints"))

p1 + p2 + ggsave(paste0(resDir, '/Forte2020_Umap_myclusters_cellType.Shoval_', Normalization, '.pdf'), 
                 width = 16, height = 8)

##########################################
# integrate Ren2020 and Forte2020 to have one reference using SCTransform and RPCA from Seurat
# original code from https://satijalab.org/seurat/articles/integration_rpca.html
# 
##########################################
if(Normalization == 'SCT'){
  cms = readRDS(file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_SCT_umap.rds'))
  
  ref.list = list(aa, cms)
  
  features <- SelectIntegrationFeatures(object.list = ref.list, nfeatures = 3000, assay = c('SCT', 'SCT'))
  #features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
  ref.list <- PrepSCTIntegration(object.list = ref.list, anchor.features = features)
  
  ref.list <- lapply(X = ref.list, FUN = RunPCA, features = features)
  
  # this command creates an 'integrated' data assay
  ref.anchors <- FindIntegrationAnchors(object.list = ref.list, normalization.method = "SCT",
                                        anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
  
  ref.sct <- IntegrateData(anchorset = ref.anchors, normalization.method = "SCT", dims = 1:30)
  
  ref.sct <- RunPCA(ref.sct, verbose = FALSE)
  ref.sct <- RunUMAP(ref.sct, reduction = "pca", dims = 1:30)
  
  DefaultAssay(ref.sct) <- "integrated"
  
  # Visualization
  p1 <- DimPlot(ref.sct, reduction = "umap", group.by = "stim")
  p2 <- DimPlot(ref.sct, reduction = "umap", group.by = "seurat_annotations", label = TRUE,
                repel = TRUE)
  p1 + p2
  
}else{
  
  aa = readRDS(file = paste0(RdataDir, 'Forte2020_logNormalize_allgenes.rds'))
  cms = readRDS(file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization_umap.rds'))
  aa$dataset = 'Forte2020'
  cms$dataset = 'Ren2020'
  aa$annot.ref = aa$my_annot
  cms$annot.ref = cms$CellType
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
  
  saveRDS(ref.combined, file = paste0(RdataDir, 'Seurat.obj_adultMiceHeart_Forte2020_Ren2020_refCombined_logNormalize_v1.rds'))
  
  # Run the standard workflow for visualization and clustering
  ref.combined = readRDS(file =paste0(RdataDir, 'Seurat.obj_adultMiceHeart_Forte2020_Ren2020_refCombined_logNormalize_v1.rds'))
  
  ref.combined <- ScaleData(ref.combined, verbose = FALSE)
  ref.combined <- RunPCA(ref.combined, npcs = 30, verbose = FALSE)
  
  ElbowPlot(ref.combined, ndims = 30)
  
  ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:20)
  ref.combined <- FindClusters(ref.combined, resolution = 0.5)
  
  ref.combined <- RunUMAP(ref.combined, reduction = "pca", dims = 1:30, n.neighbors = 50, min.dist = 0.05) 
  
  
  # Visualization
  p1 <- DimPlot(ref.combined, reduction = "umap", group.by = "dataset")
  p2 <- DimPlot(ref.combined, reduction = "umap", group.by = "annot.ref", label = TRUE,
                repel = TRUE)
  p1 + p2 + ggsave(paste0(resDir, '/Forte2020_Ren2020_IntegrationRPCA_', Normalization, '.pdf'), 
                   width = 24, height = 10)
  
  
  ##########################################
  # clean the reference, i.e. remove the non-cardiomyocyte from Ren2020
  # change the confusing annotation names from Shoval
  ##########################################
  kk = which(ref.combined$dataset == 'Forte2020'| (ref.combined$dataset == 'Ren2020' & ref.combined$annot.ref == 'CM'))
  refs = ref.combined[,kk]
  
  rm(ref.combined)
  
  saveRDS(refs, file = paste0(RdataDir, 
                              'Seurat.obj_adultMiceHeart_Forte2020.nonCM_Ren2020CM_refCombined_cleanAnnot_logNormalize_v1.rds'))
  
  # clean a bit the annotation
  refs$celltype = refs$annot.ref
  refs$celltype[which(refs$celltype == '0 - Fibro-I')] = 'FB1'
  refs$celltype[which(refs$celltype == '9 - Fibro-II')] = 'FB2'
  refs$celltype[which(refs$celltype == '11 - Fibro-III')] = 'FB3'
  refs$celltype[which(refs$celltype == '15 - Fibro-IV')] = 'FB4'
  
  refs$celltype[which(refs$celltype == '2 - EC-I')] = 'EC1'
  refs$celltype[which(refs$celltype == '18 - EC-II')] = 'EC2'
  refs$celltype[which(refs$celltype == '19 - EC-III')] = 'EC3'
  refs$celltype[which(refs$celltype == '17 - Lymph-EC')] = 'Lymph.EC'
  
  refs$celltype[which(refs$celltype == '6 - B cell')] = 'B'
  refs$celltype[which(refs$celltype == '3 - MyoF')] = 'myoFB'
  refs$celltype[which(refs$celltype == '16 - SMC')] = 'sMC'
  
  refs$celltype[which(refs$celltype == '4 - GRN')] = 'GN'
  refs$celltype[which(refs$celltype == '12 - NK/T')] = 'NK.T'
  refs$celltype[which(refs$celltype == '13 - Monon/DC')] = 'MCT.DC'
  refs$celltype[which(refs$celltype == '5 - Chil3 Mono')] = 'MCT.Chil3'
  
  refs$celltype[which(refs$celltype == '1 - Trem2 Macs')] = 'Mphage.Trem2'
  refs$celltype[which(refs$celltype == '10 - Arg1 Macs')] = 'Mphage.Argl'
  refs$celltype[which(refs$celltype == '7 - MHCII Macs')] = 'MHCII.Mphage'
  refs$celltype[which(refs$celltype == '14 - Proliferating Macs')] = 'prolife.Mphage'
  
  saveRDS(refs, file = paste0(RdataDir, 
                              'SeuratObj_adultMiceHeart_refCombine_Forte2020.nonCM_Ren2020CM_cleanAnnot_logNormalize_v1.rds'))
  
  
  # test refs without integration
  Test_refs_withoutIntegration = FALSE
  if(Test_refs_){
    DefaultAssay(refs) <- "RNA"
    refs = FindVariableFeatures(refs, selection.method = "vst", nfeatures = 3000)
    
    # Run the standard workflow for visualization and clustering
    refs <- ScaleData(refs, verbose = FALSE)
    refs <- RunPCA(refs, npcs = 30, verbose = FALSE)
    
    ElbowPlot(refs, ndims = 30)
    
    refs <- FindNeighbors(refs, reduction = "pca", dims = 1:20)
    refs <- FindClusters(refs, resolution = 0.5)
    
    refs <- RunUMAP(refs, reduction = "pca", dims = 1:30, n.neighbors = 50, min.dist = 0.05) 
    
    # Visualization
    p1 <- DimPlot(refs, reduction = "umap", group.by = "dataset")
    p2 <- DimPlot(refs, reduction = "umap", group.by = "annot.ref", label = TRUE,
                  repel = TRUE)
    p1 + p2 + ggsave(paste0(resDir, '/Forte2020_Ren2020_noCorrection_', Normalization, '.pdf'), 
                     width = 24, height = 10)
    
  }
  
  ##########################################
  # try to reversely calculated batch-corrected UMI counts using corrected gene expression matrix from Seurat 
  ##########################################
  refs = readRDS(file = paste0(RdataDir, 
                               'SeuratObj_adultMiceHeart_refCombine_Forte2020.nonCM_Ren2020CM_cleanAnnot_logNormalize_v1.rds'))
  
  p1 <- DimPlot(refs, reduction = "umap", group.by = "dataset")
  p2 <- DimPlot(refs, reduction = "umap", group.by = "celltype", label = TRUE,
                repel = TRUE)
  p1 + p2
  
  metadata = refs@meta.data   
  Ec = refs@assays$integrated@data 
  counts = refs@assays$RNA@counts
  counts = counts[match(rownames(Ec), rownames(counts)), ]
  
  
  
}




