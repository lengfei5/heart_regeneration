##########################################################################
##########################################################################
# Project:
# Script purpose: regress out the nCount_RNA in axolotl snRNA-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Sep  7 11:35:08 2022
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_R13591_intron.exon.20220729'

resDir = paste0("../results/sc_multiome", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')


if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R13591_axolotl_multiome'

source('functions_scRNAseq.R')
source('functions_Visium.R')

require(Seurat)
#require(sctransform)
library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
require(tictoc)
options(future.globals.maxSize = 64000 * 1024^2)

mem_used()

species = 'axloltl_scRNAseq'

##########################################
# main script 
##########################################
aa = readRDS(file = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets.rds')

Idents(aa) = aa$subtypes
#DimPlot(aa, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()

#VlnPlot(aa, features = 'nCount_RNA', group.by = 'condition')

aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000)
all.genes <- rownames(aa)

aa <- ScaleData(aa, features = all.genes, vars.to.regress = 'nCount_RNA')

saveRDS(aa, file = paste0(RdataDir, 'aa_annotated_no_doublets_Elad_regressed.nCountRNA.rds'))

