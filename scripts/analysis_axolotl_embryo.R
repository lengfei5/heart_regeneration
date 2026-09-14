##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: analyze the axolotl cardiac_lineage data from Sebastian
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Sep 14 13:33:26 2026
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_20260914'

resDir = paste0("../results/axoltol_embryo_cardiac_lineage", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/groups/tanaka/Collaborations/Jingkui-Elad/From_Sebestian_embryo'

source('functions_scRNAseq.R')
source('functions_Visium.R')
species = 'axloltl_scRNAseq'

require(Seurat)
#require(sctransform)3
library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
require(tictoc)
options(future.globals.maxSize = 80000 * 1024^2)

mem_used()


########################################################
########################################################
# Section I: double check the QC and overview of the data
# 
########################################################
########################################################

##########################################
# this option did not work
##########################################
library(Seurat)
library(SeuratDisk)
# Convert h5ad -> h5seurat first
Convert(paste0(dataDir, "/adata_cardiac_lineage_subset.h5ad"), dest = "h5seurat", 
        overwrite = TRUE)

# Load the h5seurat file
seurat_obj <- LoadH5Seurat("data.h5seurat")


##########################################
# try sceasy option
##########################################
library(sceasy)
library(reticulate)

use_condaenv('scvi-env')

sceasy::convertFormat(paste0(dataDir, "/adata_cardiac_lineage_subset.h5ad"), 
                      from = "anndata", to = "seurat",
                      outFile = "data.rds")



