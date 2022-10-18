##########################################################################
##########################################################################
# Project: Heart regeneration 
# Script purpose: for Elad to check marker genes and annotate cell types and subtypes
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Aug  9 14:55:17 2022
##########################################################################
##########################################################################
rm(list = ls())
require(Seurat)
#require(sctransform)
#library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)

SeuratObj = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/Rdata_spliced/'
resDir = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/res_Elad'

aa = readRDS(file = paste0(SeuratObj, 'seuratObject_axloltl_scRNAseq_R13591_20220720_lognormamlized_pca_umap.rds'))

DimPlot(aa, label = TRUE, repel = TRUE) + ggtitle("scNuc (multiome)")

##########################################
# start to annotate the cell types
##########################################
features = rownames(aa)[grep('LY6', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))
VlnPlot(aa, features = features)

features = rownames(aa)[grep('PTPRC', rownames(aa))]

FeaturePlot(aa, features = features, cols = c('gray', 'red'))
VlnPlot(aa, features = features)

#ggsave(filename = paste0(resDir, '/FeaturePlot_CD45.pdf'), width = 8, height = 6)

features = rownames(aa)[grep('VIM|COL1A2|FSTL1|POSTN', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('MYH6|ACTN2|NPPA|TNNT2|GATA4', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('CD68|CD8A|CD74|CSF1R|ITGAM', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('MKI67|CCNB2|PCNA-|CDK1-', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

