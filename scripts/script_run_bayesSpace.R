##########################################################################
##########################################################################
# Project:
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Oct  4 10:00:26 2022
##########################################################################
##########################################################################
rm(list = ls())

species = 'axolotl'
version.analysis = '_R12830_resequenced_20220308'
dataDir = '../R12830_visium_reseqenced/nf_out'
resDir = paste0("../results/visium_axolotl", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

source('functions_Visium.R')
library(pryr) # monitor the memory usage
require(ggplot2)
mem_used()

options(future.globals.maxSize = 64000 * 1024^2)

##########################################
# main script
##########################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered', species, '.Rdata'))
st$condition = factor(st$condition, levels = design$condition)

cat('visium conditions :\n')
print(table(st$condition))
cc = design$condition

source('functions_Visium.R')
run_bayesSpace(st, outDir = paste0(resDir, '/bayesSpace/'))

