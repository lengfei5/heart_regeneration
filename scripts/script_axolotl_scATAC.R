##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: analyze scATAC-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Aug 25 14:24:06 2022
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_R13591_atac.20220825'

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

########################################################
########################################################
# Section I: import the scRNAseq data by kalisto
# check QCs 
########################################################
########################################################
design = data.frame(sampleID = seq(197249, 197253), 
                    condition = c(paste0('Amex_scRNA_d', c(0, 1, 4, 7, 14))), stringsAsFactors = FALSE)

source('functions_scRNAseq.R')

for(n in 1:nrow(design))
{
  # n = 1
  cat('-----------', design$condition[n], '-------------\n')
  
  # load nf output and process
  topdir = paste0(dataDir, '/', design$condition[n], '/', design$condition[n], '/')
  
  aa = make_SeuratObj_scRNAseq(topdir = topdir,
                               saveDir = paste0(resDir, '/', design$condition[n], '_', design$sampleID[n], '/'), 
                               changeGeneName.axolotl = TRUE, 
                               defaultDrops.only = TRUE)
  
  aa$condition = design$condition[n]
  
  if(n == 1) {
    scn = aa
  }else{
    scn = merge(scn, aa)
  }
  
}

rm(aa)

save(design, scn, 
     file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, version.analysis, '.Rdata'))