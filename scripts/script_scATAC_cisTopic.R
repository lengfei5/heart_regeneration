##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: test cisTopic 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Jan 20 16:42:18 2023
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_R13591_atac_reseq_20221115'

resDir = paste0("../results/sc_multiome", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')


if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R14353_ax_snATAC_reseq'

#source('functions_scATAC.R')
#source('functions_scRNAseq.R')
#source('functions_Visium.R')

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(patchwork)
require(SeuratObject)

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
require(tictoc)
options(future.globals.maxSize = 80000 * 1024^2)
set.seed(1234)
mem_used()

library(data.table)
require(cisTopic)
##########################################
# start the main function 
##########################################
srat_cr = readRDS(file = paste0(RdataDir, 
                                'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                'scATAC.merged.peaks.cr.',
                                '584K.annot_38280cells.rds'))
# normalize ATAC and UMAP
DefaultAssay(srat_cr) <- "ATAC"

Filter.cells.with.scATAC = FALSE
if(Filter.cells.with.scATAC){
  # quick filtering 
  srat_cr <- subset(
    x = srat_cr,
    subset = nCount_ATAC < 100000 &
      nCount_RNA < 25000 &
      nCount_ATAC > 200 &
      nCount_RNA > 1000 &
      nucleosome_signal < 6 &
      TSS.enrichment > 1
  )
  
}

Idents(srat_cr) = as.factor(srat_cr$condition)
#srat_cr = subset(srat_cr,  downsample = 1000)

run_dim_reduction = function(atac_matrix, cell_embeddings, dims, metadata=NULL, reduction='pca.l2') {
  if (is.null(metadata)) {
    seurat_obj = Seurat::CreateSeuratObject(atac_matrix, assay = 'peaks')
  } else {
    seurat_obj = Seurat::CreateSeuratObject(atac_matrix, meta.data = metadata, assay = 'peaks')
  }
  
  seurat_obj[['pca']] = Seurat::CreateDimReducObject(embeddings=cell_embeddings, key='PC_', 
                                                     assay='peaks')
  seurat_obj = seurat_obj %>%
    Seurat::L2Dim(reduction='pca') %>%
    Seurat::RunUMAP(reduction = reduction, dims = dims) %>%
    Seurat::RunTSNE(reduction = reduction, dims = dims) %>%
    Seurat::FindNeighbors(reduction=reduction, nn.eps=0.25, dims=dims)
  return(seurat_obj)
}


#data(counts_mel) 
bmat = srat_cr@assays$ATAC@counts
coords = str_split_fixed(rownames(bmat), '-', 3)
new_coords = paste0(coords[, 1], ':', coords[, 2], '-', coords[, 3])
rownames(bmat) = new_coords

cisTopicObject = cisTopic::createcisTopicObject(bmat, project.name='ax_heart')

rm(bmat)
#topic=seq(20, 120, by=10)
#topic = c(10, 20)
topic=c(seq(2, 18, by = 2), seq(20, 160, by=10))

cisTopicObject = cisTopic::runWarpLDAModels(cisTopicObject, 
                                            topic = topic, 
                                            seed=2023, 
                                            nCores=32, 
                                            tmp = "../results/tmp/",
                                            iterations = 500, 
                                            addModels = FALSE)


saveRDS(cisTopicObject, file = paste0(RdataDir, 'test_cisTopic_saved_v2.rds'))

# cisTopicObject = readRDS(file = paste0(RdataDir, 'test_cisTopic_saved.rds'))
# 
# par(mfrow=c(3,3))
# cisTopicObject <- selectModel(cisTopicObject, type='maximum')
# cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
# cisTopicObject <- selectModel(cisTopicObject, type='derivative')
# 
# cistopicObject = cisTopic::selectModel(cisTopicObject, type = 'maximum')
# 
# method='Z-score'
# #method = 'Probability'
# cistopicObject.reduced_space = t(cisTopic::modelMatSelection(cistopicObject,
#                                                              target='cell',
#                                                              method=method))
# 
# colnames(cistopicObject.reduced_space) = paste0('PC_', 1:ncol(cistopicObject.reduced_space))
# 
# dimensions = ncol(cistopicObject.reduced_space)
# 
# #cistopicObject.seurat = run_dim_reduction(,
# #                                          cistopicObject.reduced_space,
# #                                          dims=1:120,
# #                                          reduction='pca')
# 
# seurat_obj = Seurat::CreateSeuratObject(cistopicObject@binary.count.matrix, assay = 'peaks')
# seurat_obj[['pca']] = Seurat::CreateDimReducObject(embeddings=cistopicObject.reduced_space, 
#                                                    key='PC_', 
#                                                    assay='peaks')
# 
# seurat_obj$subtypes = srat_cr$subtypes[match(colnames(seurat_obj), colnames(srat_cr))]
# seurat_obj$celltypes = srat_cr$celltypes[match(colnames(seurat_obj), colnames(srat_cr))]
# 
# # seurat_obj = Seurat::L2Dim(seurat_obj, reduction='pca') 
# seurat_obj = Seurat::RunUMAP(seurat_obj, reduction = "pca", dims = 1:120, n.neighbors = 30, 
#                              min.dist = 0.05)
# 
# DimPlot(seurat_obj, reduction = 'umap', group.by = 'celltypes', label = TRUE, repel = TRUE) + 
#   NoLegend()
# 
# ggsave(paste0(resDir, '/UMAP_DM_cisTopic_test.pdf'), width = 8, height = 6)
# 
# DimPlot(seurat_obj, reduction = 'umap', group.by = 'subtypes', label = TRUE, repel = TRUE) + NoLegend()

#cistopicObject.seurat = cistopicObject.seurat %>%
#  Seurat::FindNeighbors(reduction=reduction, nn.eps=0.25, dims=1:dimensions) %>%
#  Seurat::FindClusters(reduction=reduction, n.start=20, resolution=resolution)

#srat_cr[['Topic']] = Seurat::CreateDimReducObject(embeddings=cistopicObject.reduced_space, 
#                                                key='Topic_', 
#                                                   assay='ATAC')

