##########################################################################
##########################################################################
# Project: heart regeneration project
# Script purpose: First script to analyze the processed Visium data by spaceranger of 10x
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Aug 27 11:39:22 2021
##########################################################################
##########################################################################
# setup for data import and sequencing QCs
version.analysis = '_R11934_20210827'

resDir = paste0("../results/visium_mouse", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R11934_visium'

require(Seurat)
require(SeuratObject)
require(ggplot2)
require(tibble)
require(dplyr)

########################################################
########################################################
# Section : 
# 
########################################################
########################################################
