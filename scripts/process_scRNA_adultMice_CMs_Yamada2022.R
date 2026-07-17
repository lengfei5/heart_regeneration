##########################################################################
##########################################################################
# Project: Heart regeneration 
# Script purpose: process the CMs scRNA-seq data from Yamada et al., 2022
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Jul 14 16:20:54 2026
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_scRNAseq_Yamada_20260714'

resDir = paste0("../results/scRNAseq_adultMice_CMs", version.analysis,'/')
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

library(dplyr)
library(Seurat)
library(patchwork)
library(DropletUtils)

sessionInfo()

source('functions_scRNAseq.R')
source('functions_Visium.R')
#source('utility_zebrafish.R')
mem_used()

dataDir = "../published_dataset/adult_mice/Yamada_2022/GSE176092_RAW/"


########################################################
########################################################
# Section I: import the processed tables
# original code from https://github.com/firstheart123/spatiotemporal_heart/blob/main/snRNA_seq.R
########################################################
########################################################
file_list = list.files(path = dataDir, all.files = TRUE, full.names = FALSE, no.. = TRUE)


samples = gsub('_barcodes.tsv.gz', '', file_list)
samples = gsub('_features.tsv.gz', '', samples)
samples = gsub('_matrix.mtx.gz|.tif.gz|_tissue_positions_list.csv.gz', '', samples)
samples = unique(samples)
samples = samples[grep('_CM_', samples)]


# import data from cellranger output
for(n in 1:length(samples))
{
  # n = 1
  cat(n, ' : ', samples[n], '\n')
  
  topdir = paste0(dataDir)
  
  exp = Matrix::readMM(paste0(topdir, samples[n], "_matrix.mtx.gz")) #read matrix
  bc = read.csv(paste0(topdir, samples[n], "_barcodes.tsv.gz"), header = F, stringsAsFactors = F)
  g = read.csv(paste0(topdir, samples[n], "_features.tsv.gz"), header = F, stringsAsFactors = F, sep = '\t')
  
  ## make unique gene names
  g$name = g$V2
  gg.counts = table(g$V2)
  gg.dup = names(gg.counts)[which(gg.counts>1)]
  index.dup = which(!is.na(match(g$V2, gg.dup)))
  g$name[index.dup] = paste0(g$V2[index.dup], '_', g$V1[index.dup])
  
  colnames(exp) = bc$V1
  rownames(exp) = g$name
  
  count.data = exp
  rm(exp);
  
  cat('get empty drops with UMI rank \n')
  
  # get emptyDrops and default cutoff cell estimates
  iscell_dd = defaultDrops(count.data, expected = 8000) # default cell estimate, similar to 10x cellranger
  sum(iscell_dd, na.rm=TRUE)
  
  ## not used the emptyDrops too slow 
  # eout = emptyDrops(count.data, lower = 200)
  # eout$FDR[is.na(eout$FDR)] = 1
  # iscell_ed = eout$FDR<=0.01
  # sum(iscell_ed, na.rm=TRUE)
  
  meta = data.frame(row.names = colnames(count.data), condition = samples[n],
                    iscell_dd = iscell_dd)
  
  # use defaultDrop to select cells.
  aa = CreateSeuratObject(counts = count.data[, iscell_dd],
                          meta.data = meta[iscell_dd, ], 
                          min.cells = 0, min.features = 0)
  aa$cell.id = paste0(samples[n], '_', colnames(aa))
  
  if(n == 1) {
    mnt = aa
  }else{
    mnt = merge(mnt, aa)
  }
}

#mnt[["percent.mt"]] <- PercentageFeatureSet(mnt, pattern = "^mt-")

saveRDS(mnt, file = paste0(RdataDir, 'seuratObject_',  version.analysis, '.rds'))


##########################################
# process the data
##########################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_',  version.analysis, '.rds'))

#aa = mnt 
#rm(mnt)

#images = file_list[grep('tif', file_list)]
aa$condition[grep('CM_Sham_1', aa$condition)] = 'Sham_1'
aa$condition[grep('CM_Sham_2', aa$condition)] = 'Sham_2'

aa$condition[grep('CM_MI_day1_IZ\\+BZ_1', aa$condition)] = 'MI_day1_IR_1'
aa$condition[grep('CM_MI_day1_IZ\\+BZ_2', aa$condition)] = 'MI_day1_IR_2'
aa$condition[grep('CM_MI_day1_RZ_1', aa$condition)] = 'MI_day1_RR_1'
aa$condition[grep('CM_MI_day1_RZ_2', aa$condition)] = 'MI_day1_RR_2'

aa$condition[grep('CM_MI_day7_IZ\\+BZ_1', aa$condition)] = 'MI_day7_IR_1'
aa$condition[grep('CM_day7_IZ\\+BZ_2', aa$condition)] = 'MI_day7_IR_2'
aa$condition[grep('CM_MI_day7_RZ_1', aa$condition)] = 'MI_day7_RR_1'
aa$condition[grep('CM_MI_day7_RZ_2', aa$condition)] = 'MI_day7_RR_2'

aa$condition[grep('CM_MI_day14_IZ\\+BZ_1', aa$condition)] = 'MI_day14_IR_1'
aa$condition[grep('CM_MI_day14_IZ\\+BZ_2', aa$condition)] = 'MI_day14_IR_2'
aa$condition[grep('CM_MI_day14_RZ_1', aa$condition)] = 'MI_day14_RR_1'
aa$condition[grep('CM_MI_day14_RZ_2', aa$condition)] = 'MI_day14_RR_2'


aa[["percent.mt"]] <- PercentageFeatureSet(aa, pattern = "^mt-")

mt.index <- grep(pattern = "^mt-", x = rownames(aa[["RNA"]]), value = FALSE)
aa_matrix <- aa[["RNA"]][-mt.index, ]
aa <- CreateSeuratObject(counts = aa_matrix, meta.data = aa@meta.data)

aa_list <- SplitObject(aa, split.by = "condition")

aa_list <- lapply(X = aa_list, FUN = function(x) {
  x <- subset(x, subset = nFeature_RNA > 500 & percent.mt < 60)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = aa_list)

aa_anchors <- FindIntegrationAnchors(object.list = aa_list, anchor.features = features)

combined <- IntegrateData(anchorset = aa_anchors)

combined <- ScaleData(combined, features = rownames(aa)) %>% 
  RunPCA(npcs = 30, features = VariableFeatures(combined)) %>% 
  RunUMAP(reduction = "pca", dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) 

combined <- FindClusters(combined, resolution = 0.25)

DimPlot(combined, reduction = "umap", pt.size = 1., label = TRUE)


saveRDS(combined, file = paste0(RdataDir, 'seuratObject_CCAIntegrated',  version.analysis, '.rds'))


##########################################
# subset the CMs
##########################################
combined = readRDS(file = paste0(RdataDir, 'seuratObject_CCAIntegrated',  version.analysis, '.rds'))

combined <- FindClusters(combined, resolution = 0.2)

DefaultAssay(combined) <- "RNA"

p1 = DimPlot(combined, reduction = "umap", pt.size = 1., label = TRUE)

p2 = FeaturePlot(combined, features= c('Myh6', 'Nppa', 'Nppb', 'Tnni3','Tnnt2', 'Actc1', 'Ttn',
                                       'Acta1', 'Myl2', 'Tnnc1', 'Actn2'))

p1 + p2


ggsave(paste0(resDir, 'UMAP_clusters_CMmarkers.pdf'), 
       width = 20, height = 8)


CM <- subset(combined, idents = c("0", "4"))

CM <- FindVariableFeatures(CM, selection.method = "vst", nfeatures = 1000) %>% 
  ScaleData(features = rownames(CM))

CM <- RunPCA(CM, npcs = 30, features = VariableFeatures(CM))
ElbowPlot(aa, ndims = 30)

CM =  RunUMAP(reduction = "pca", dims = 1:20) %>% FindNeighbors(reduction = "pca", dims = 1:20)

CM <- FindClusters(CM, resolution = 0.2)

DimPlot(CM, reduction = "umap", pt.size = 1.5)

cluster0_marker <- FindMarkers(CM, ident.1 = 0, only.pos = TRUE, logfc.threshold = 0.25) %>% 
  dplyr::filter(p_val_adj < 0.05)
cluster1_marker <- FindMarkers(CM, ident.1 = 1, only.pos = TRUE, logfc.threshold = 0.25) %>% 
  dplyr::filter(p_val_adj < 0.05)
cluster2_marker <- FindMarkers(CM, ident.1 = 2, only.pos = TRUE, logfc.threshold = 0.25) %>% 
  dplyr::filter(p_val_adj < 0.05)


