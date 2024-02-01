##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: run the cell type deconvolution of visium data using RCTD
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Feb  1 09:45:07 2024
##########################################################################
##########################################################################
rm(list = ls())

library(pryr) # monitor the memory usage
require(ggplot2)
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
library(circlize)
library(RColorBrewer)
require(scran)
require(scater)
library(pryr) # monitor the memory usage
require(ggplot2)

source('functions_scRNAseq.R')
source('functions_Visium.R')


Run_RCTD_adultMice = TRUE
if(Run_RCTD_adultMice){
  dataPath_nichenet = '../data/NicheNet/'
  version.analysis = '_R11934_20240131'
  
  resDir = paste0("../results/visium_adultMice", version.analysis)
  RdataDir = paste0(resDir, '/Rdata/')
  
  if(!dir.exists(resDir)) dir.create(resDir)
  if(!dir.exists(RdataDir)) dir.create(RdataDir)
  
  dataDir = '../R11934_visium'
  
  
  mem_used()
  species = 'mouse_adult'
  
  
  load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, 
                     '_umap.clustered.Rdata'))
  
  refs = readRDS(file = paste0('../data/data_examples/ref_scRNAseq_adultMice_clean.v1.rds'))
  
  #jj = which(refs$dataset == 'Ren2020')
  #refs$timepoints[jj] = refs$condition[jj]
  #refs$condition = refs$timepoints
  
  p1 = DimPlot(refs, reduction = 'umap', group.by = 'dataset',  label = TRUE, repel = TRUE)
  p2 = DimPlot(refs, reduction = 'umap', group.by = 'celltype',  label = TRUE, repel = TRUE)
  p3 = DimPlot(refs, reduction = 'umap', group.by = 'subtype',  label = TRUE, repel = TRUE)
  p4 = DimPlot(refs, reduction = 'umap', group.by = 'timepoints',  label = TRUE, repel = TRUE)
  
  (p1 + p2)/(p3 + p4)
  
  ggsave(filename = paste0(resDir, '/UMAP_scRNAseq_refrence_dataset_timepoints_celltypes.pdf'), 
         width = 24, height = 12)
  
  
  refs$celltype_toUse = as.character(refs$subtype)
  length(table(refs$celltype_toUse))
  table(refs$celltype_toUse)
  #DimPlot(refs, reduction = 'umap', group.by = 'celltype_toUse')
  
  ## prepare the celltype to use and also specify the time-specific subtypes
  table(refs$condition)
  
  table(refs$celltype_toUse)
  length(table(refs$celltype_toUse))
  st$condition = factor(st$condition)
  table(st$condition)
  
  ## preapre the paramters for RCTD subtypes
  DefaultAssay(refs) = 'integrated'
  DefaultAssay(st) = 'Spatial'
  require_int_SpatialRNA = FALSE
  
  condition.specific.ref = FALSE
  
  outDir = paste0(resDir, '/celltype_deconvolution')
  RCTD_out = paste0(outDir, '/RCTD_', length(table(refs$celltype_toUse)), 'Subtype_ref_v1.0')
  
  max_cores = 16
  
  # st = subset(st, condition == 'adult.day7'); st$condition = droplevels(st$condition)
  Run.celltype.deconvolution.RCTD(st, refs, 
                                  condition.specific.ref = condition.specific.ref,
                                  #condition.specific_celltypes = condition.specific_celltypes,
                                  require_int_SpatialRNA = require_int_SpatialRNA,
                                  max_cores = max_cores,
                                  RCTD_out = RCTD_out,
                                  plot.RCTD.summary = FALSE, 
                                  PLOT.scatterpie = FALSE
                                  
  )
  
  plot.RCTD.results(st = st, 
                    RCTD_out = RCTD_out,
                    species = species,
                    plot.RCTD.summary = FALSE)
  
  
  
}
