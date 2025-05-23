##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: process the zebrafish scRNA-seq data from the Junker's paper
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Feb 27 16:11:33 2025
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_scRNAseq_20250227'

resDir = paste0("../results/scRNAseq_zebrafish", version.analysis,'/')
RdataDir = paste0(resDir, 'Rdata/')

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
library(Seurat)
library(Matrix)
library(ggplot2)
library(cowplot)
library(dplyr)
library(gridExtra)
library(ggrepel)
library(reshape2)


source('functions_scRNAseq.R')
source('functions_Visium.R')
source('utility_zebrafish.R')
mem_used()

dataDir = "../published_dataset/zebrafish/Hu_Junker_2022/zebrafish_heart_processing/"


########################################################
########################################################
# Section I: import the data and prepare the scRNA-seq data
# 
########################################################
########################################################
## load and integrate data in one seurat object, calculate mito reads, filter cells, 
## annotate cell types. ####


#mito.genes names
mito.genes <- read.table("Data/mito.genes.vs.txt",sep = ",")
mito.genes <- mito.genes$V3
mito.genes <- as.character(mito.genes)

# load all data
counts_folder <- "../GSE159032_RAW"
file_list = list.files(path = counts_folder, pattern = '*.h5', full.names = TRUE)

bcs = sapply(basename(file_list), function(x){unlist(strsplit(as.character(x), '_'))[2]})
batches <- c("H5","H6","H7","H8a","H8v","Hr1","Hr2a","Hr2b","Hr3","Hr4","Hr6a","Hr6v","Hr7a","Hr7v","Hr8",
             "Hr9","Hr10","Hr11","Hr12","Hr13","Hr14","Hr15","Hr19","Hr20","Hr21","Hr22","Hr23",
             "Hr24","Hr25","Hr26","Hr27","Hr28","Hr29","Hr30","Hr31","Hr32","Hr33","Hr34","Hr35")

bcs[which(is.na(match(bcs, batches)))]
batches[which(is.na(match(batches, bcs)))]

gather.data <- list()
data.stats <- numeric()

for (i in batches) {
  kk = grep(paste0('_', i, '_filtered_matrix.h5'), file_list)
  cat(i, ' -- ', file_list[kk], '\n')
  if(length(kk) == 1){
    gather.data[[i]] <- Read10X_h5(filename = file_list[kk])
    # Select row indices and not ERCC names 
    RFP.index <- grep(pattern = "^RFP", x = rownames(gather.data[[i]]), value = FALSE) 
    gather.data[[i]] <- gather.data[[i]][-RFP.index, ]
    gather.data[[i]] <- CreateSeuratObject(counts = gather.data[[i]],
                                           min.cells = 3, min.features = 150,
                                           project = i)
    mito.genes.use <- setdiff(mito.genes,setdiff(mito.genes,rownames(gather.data[[i]][["RNA"]])))
    gather.data[[i]][["percent.mito"]] <- PercentageFeatureSet(object = gather.data[[i]], 
                                                               features = mito.genes.use)
    gather.data[[i]] <- subset(x = gather.data[[i]], subset = percent.mito < 25 & nFeature_RNA < 4100)
    data.stats[i] <- length(colnames(x = gather.data[[i]]))
    
  }else{
    stop('no file found \n')
  }
}

# stats
print(data.stats)
# merge all data
all.hearts <- merge(gather.data[[1]], y = gather.data[-1], add.cell.ids = batches, project = "allheart")

saveRDS(all.hearts, file = 'Rdata/all_heart_beforeAnnotation.rds')

#remove raw data to free space
rm(gather.data)

##########################################
# cell annotation
##########################################
all.hearts = readRDS(file = 'Rdata/all_heart_beforeAnnotation.rds')
all.hearts$batch = all.hearts$orig.ident

#Annotate metadata
all.hearts@meta.data$time <- NA
all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("H5","H6","H7","H8a","H8v"),]$time <- "Ctrl"

all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("Hr10","Hr11","Hr12","Hr22","Hr23",
                                                            "Hr24","Hr25","Hr26","Hr27","Hr28",
                                                            "Hr29","Hr34","Hr35"),]$time <- "3dpi"

all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("Hr1","Hr2a","Hr2b",
                                                            "Hr8","Hr9","Hr13","Hr14","Hr15",
                                                            "Hr6a","Hr6v","Hr7a","Hr7v","Hr30","Hr31",
                                                            "Hr32","Hr33"),]$time <- "7dpi"

all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("Hr3","Hr4","Hr19","Hr20","Hr21"),]$time <- 
  "30dpi"

all.hearts@meta.data$time <- factor(x = all.hearts@meta.data$time, 
                                    levels = c("Ctrl","3dpi","7dpi","30dpi"))

all.hearts@meta.data$AV <- "Wholeheart"
all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("H8a","Hr6a","Hr7a"),]$AV <- "Atrium"
all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("H8v","Hr6v","Hr7v","Hr25"),]$AV <- "Ventricle"

all.hearts@meta.data$inhib <- "NULL"
all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("Hr28","Hr30"),]$inhib <- "DMSO"
all.hearts@meta.data[all.hearts@meta.data$orig.ident %in% c("Hr29","Hr31","Hr32",
                                                            "Hr33","Hr34","Hr35"),]$inhib <- "IWR1"

ncol(all.hearts)


# add annotated cell type names with filtering out erythrocytes and dead cells
cell.id.table <- read.csv("Data/cell.id.table.csv", row.names = "X")
row.names(cell.id.table) <- paste0(row.names(cell.id.table), "-1")

Correct_cellids_batch.Hr22.Hr23 = TRUE
if(Correct_cellids_batch.Hr22.Hr23){
  
  total = sapply(rownames(cell.id.table), function(x){unlist(strsplit(as.character(x), '_'))[1]})
  table(total)
  
  
  cellid = rownames(cell.id.table)
  jj1 = grep('Hr22_', cellid)
  jj2 = grep('Hr23_', cellid)
  
  cellid[jj1] = gsub('Hr22_', 'Hr23_', cellid[jj1])
  cellid[jj2] = gsub('Hr23_', 'Hr22_', cellid[jj2])
  
  rownames(cell.id.table) = cellid
  
  missed_index = which(is.na(match(rownames(cell.id.table), colnames(all.hearts))))
  missed_index
  
  missed = rownames(cell.id.table)[missed_index]
  missed = sapply(missed, function(x){unlist(strsplit(as.character(x), '_'))[1]})
  
  table(missed)
  
  head(cell.id.table[grep('Hr22', rownames(cell.id.table)), ])
  
}

all.hearts <- AddMetaData(all.hearts, cell.id.table)
all.hearts <- subset(all.hearts, 
                     cells = rownames(all.hearts@meta.data)[!is.na(all.hearts@meta.data$celltypes)], )

saveRDS(all.hearts, file = 'Rdata/all_heart_annotation.rds')


## normalize, pca, cluster and plot umap
all.hearts = readRDS(file = 'Rdata/all_heart_annotation.rds')
all.hearts$subtypes = all.hearts$subclustered.celltypes

all.hearts <- NormalizeData(object = all.hearts)

all.hearts <- FindVariableFeatures(object = all.hearts, selection.method = "vst",nfeatures = 3500)
all.hearts <- ScaleData(object = all.hearts, 
                        vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mito"))

all.hearts <- RunPCA(object = all.hearts, features = VariableFeatures(object = all.hearts),npcs = 100)
ElbowPlot(object = all.hearts,ndims = 100)

all.hearts <- FindNeighbors(object = all.hearts, dims = 1:50)
all.hearts <- FindClusters(object = all.hearts, resolution = 4)
all.hearts <- RunUMAP(all.hearts, dims = 1:15, verbose = F)

saveRDS(all.hearts, file = 'Rdata/all_heart_annotation_normalized.pca.umap.rds')

#plot umap
DimPlot(all.hearts, group.by = "celltypes",label = F, raster=FALSE)
ggsave(filename = 'output/UMAP_celltypes.pdf', height = 8, width = 16)

subtypes = unique(all.hearts$subtypes)

cell2keep = colnames(all.hearts)[grep('Ery|duplex', all.hearts$subtypes, invert = TRUE)]

all.hearts = subset(all.hearts, cells = cell2keep)

all.hearts = subset(all.hearts, 
                    cells = colnames(all.hearts)[grep('Dead', all.hearts$subtypes, invert = TRUE)])

saveRDS(all.hearts, file = 'Rdata/all_heart_annotation_normalized.pca.umap_filtered.duplex.dead.rds')


Annotate_from_seurat.clusters = FALSE
if(Annotate_from_seurat.clusters){
  # annotate cell types, annotation is based on differental gene expression of each cluster
  #annotate cell types without ery ####
  ct <- data.frame(cluster = 0:102, Cell.type = 0:102)
  ct[ct$cluster %in% c(87,101,18,94,83,78,18,42,9,45,82,44,11,77,46,20,
                       59,36,93,40,58,62,95,97,43,38,16),]$Cell.type <- "Immune Cells"
  ct[ct$cluster %in% c(57),]$Cell.type <- "B-cells"
  ct[ct$cluster %in% c(89,79,64),]$Cell.type <- "Neutrophils"
  ct[ct$cluster %in% c(44,46,11,77,20,36,93,40,59),]$Cell.type <- "Macrophages"
  ct[ct$cluster %in% c(94,83,87,78,101,18,42,9,82,45),]$Cell.type <- "Macrophages"
  ct[ct$cluster %in% c(62,95),]$Cell.type <- "Macrophages"
  ct[ct$cluster %in% c(58,43,38,16,97),]$Cell.type <- "T-cells"
  ct[ct$cluster %in% c(34,32,29,75,73,23,24,3,84,53,54,60,22,86,71),]$Cell.type <-"Fibroblasts"
  ct[ct$cluster %in% c(17,74,68,1,19,69,81,31,92,98),]$Cell.type <-"Smooth muscle cells"
  ct[ct$cluster %in% c(50,49,99),]$Cell.type <-"Bl.ves.EC (plvapb)"
  ct[ct$cluster %in% c(66),]$Cell.type <-"Bl.ves.EC (lyve1)"
  ct[ct$cluster %in% c(100,48,47),]$Cell.type <-"Bl.ves.EC (apnln)"
  ct[ct$cluster %in% c(0,6,27,51,72,13,15,70,37),]$Cell.type <-"Endocardium 1 (A)"
  ct[ct$cluster %in% c(2,5,28,67,33,7,21,91),]$Cell.type <-"Endocardium 1 (V)"
  ct[ct$cluster %in% c(70,37),]$Cell.type <-"Endocardium 2 (A)"
  ct[ct$cluster %in% c(21,91),]$Cell.type <-"Endocardium 2 (V)"
  ct[ct$cluster %in% c(15),]$Cell.type <-"Endocardium frzb (A)"
  ct[ct$cluster %in% c(55),]$Cell.type <-"Endocardium frzb (V)"
  ct[ct$cluster %in% c(41,61,10,26),]$Cell.type <-"Cardiomyocytes A"
  ct[ct$cluster %in% c(30,56,12,14,8,63,52,76),]$Cell.type <-"Cardiomyocytes V"
  ct[ct$cluster %in% c(4,39,35),]$Cell.type <-"Cardiomyocytes (ttn.2) A"
  ct[ct$cluster %in% c(80),]$Cell.type <-"Cardiomyocytes (ttn.2) V"
  ct[ct$cluster %in% c(90),]$Cell.type <-"Cardiomyocytes (proliferating)"
  ct[ct$cluster %in% c(65),]$Cell.type <-"Perivascular cells"
  ct[ct$cluster %in% c(25,96),]$Cell.type <-"Fibroblast-like cells"
  ct[ct$cluster %in% c(102),]$Cell.type <-"Neuronal cells"
  ct[ct$cluster %in% c(88),]$Cell.type <-"Myelin cells"
  ct[ct$cluster %in% c(85),]$Cell.type <-"Proliferating cells"
  
  final.all.hearts$"first.line.annotation" <- final.all.hearts$"seurat_clusters"
  final.all.hearts@meta.data$first.line.annotation <- 
    plyr::mapvalues(final.all.hearts@meta.data$first.line.annotation, from =ct$cluster, 
                    to = ct$Cell.type)
  
  final.all.hearts <- SetIdent(final.all.hearts,value = "first.line.annotation")
  
  #final.all.hearts <- SetIdent(final.all.hearts, 
  # cells = rownames(niche@meta.data),value = niche@meta.data$work.ident )
  #final.all.hearts <- SetIdent(final.all.hearts,
  # cells = rownames(immune@meta.data),value = immune@meta.data$work.ident )
  
}


########################################################
########################################################
# Section II: reanalyze the processed data
# 
########################################################
########################################################
aa = readRDS(file = paste0(dataDir,
                           'Rdata/all_heart_annotation_normalized.pca.umap_filtered.duplex.dead.rds'))

#plot umap
DimPlot(aa, group.by = "celltypes", label = TRUE, repel = TRUE, raster=FALSE)

ggsave(filename = paste0(resDir, 'UMAP_celltypes.pdf'), height = 8, width = 16)


DimPlot(aa, group.by = "batch", label = FALSE, raster=FALSE)

ggsave(filename = paste0(resDir, 'UMAP_batch.pdf'), height = 8, width = 16)

saveRDS(aa, file = paste0(RdataDir, 'dr_scRNAseq_39Batches_noCorrection.rds'))


##########################################
# Test batch correction  
# from the previous analysis, the fastMNN reduction seems to be good
# and also the batch-corrected expression from rpca were used 
##########################################
source('functions_dataIntegration.R')

outDir = paste0(resDir, '/batch_correction/')
if(!dir.exists(outDir)) dir.create(outDir)

aa = readRDS(file = paste0(RdataDir, 'dr_scRNAseq_39Batches_noCorrection.rds'))

p1 = DimPlot(aa, group.by = "celltypes", label = TRUE, repel = TRUE, raster=FALSE)
p2 = DimPlot(aa, group.by = "batch", label = FALSE, raster=FALSE)

p1 / p2

ggsave(filename = paste0(resDir, 'UMAP_celltypes.original_39batch.pdf'), height = 16, width = 12)



## RPCA integration
method = "Seurat_RPCA"
ref.combined = IntegrateData_Seurat_RPCA(aa, 
                                         group.by = 'batch', 
                                         nfeatures = 3000,
                                         #merge.order = 
                                         redo.normalization.scaling = TRUE,
                                         correct.all = TRUE)

DefaultAssay(ref.combined) = 'integrated'
saveRDS(ref.combined, file = paste0(outDir, 
                                    'dr_scRNAseq_39Batches_correction_SeuratRPCA.rds'))


source('functions_dataIntegration.R')
ref.combined = IntegrateData_runFastMNN(aa, group.by = 'batch', 
                                        nfeatures = 3000,
                                        #merge.order = list(list(4,3,1,2), list(8,7,5,6)),
                                        correct.all = TRUE)
DefaultAssay(ref.combined) = 'mnn.reconstructed'

saveRDS(ref.combined, file = paste0(outDir, 
                                    'dr_scRNAseq_39Batches_correction_fastMNN.rds'))





p1 = DimPlot(ref.combined, group.by = 'FineID', label = TRUE, repel = TRUE, raster=FALSE) + 
  ggtitle(method)
p2 = DimPlot(ref.combined, group.by = 'condition', label = TRUE, repel = TRUE) +
  ggtitle(method)
p3 = DimPlot(ref.combined, group.by = 'BroadID', label = TRUE, repel = TRUE, raster=FALSE) + 
  ggtitle(method)

p1 /(p3 + p2)

ggsave(filename = paste0(outDir, '/CM_Cui2020_noCM_Wang2020_P1_merged_dataIntegration_', method, '.pdf'), 
       width = 16, height = 12)

DimPlot(ref.combined, group.by = 'FineID', split.by = 'dataset', label = TRUE, repel = TRUE)
ggsave(filename = paste0(outDir, '/CM_Cui2020_noCM_Wang2020_P1_merged_subtypes.by.dataset', 
                         '_dataIntegration_', method, '.pdf'), 
       width = 16, height = 6)

DimPlot(ref.combined, group.by = 'BroadID', split.by = 'dataset', label = TRUE, repel = TRUE)
ggsave(filename = paste0(outDir, '/CM_Cui2020_noCM_Wang2020_P1_merged_celltypes.by.dataset', 
                         '_dataIntegration_', method, '.pdf'), 
       width = 16, height = 6)

DimPlot(ref.combined, group.by = 'FineID', split.by = 'condition', label = TRUE, repel = TRUE, ncol = 2)
ggsave(filename = paste0(outDir, '/CM_Cui2020_noCM_Wang2020_P1_merged_subtypes.by.conditions', 
                         '_dataIntegration_', method, '.pdf'),  width = 16, height = 12)

DimPlot(ref.combined, group.by = 'BroadID', split.by = 'condition', label = TRUE, repel = TRUE, ncol = 2)
ggsave(filename = paste0(outDir, '/CM_Cui2020_noCM_Wang2020_P1_merged_celltypes.by.conditions', 
                         '_dataIntegration_', method, '.pdf'),  width = 16, height = 12)





