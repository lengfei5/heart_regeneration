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
#library(nichenetr)
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


Run_RCTD_axolotl_allVisium = TRUE

if(Run_RCTD_axolotl_allVisium){
  #rm(list = ls())
  
  species = 'axolotl'
  version.analysis = '_R17246_R12830_allVisium_20240905'
  dataDir = '../R17246_visium_axolotl/nf_out'
  
  resDir = paste0("../results/visium_axolotl", version.analysis)
  RdataDir = paste0(resDir, '/Rdata/')
  
  if(!dir.exists(resDir)) dir.create(resDir)
  if(!dir.exists(RdataDir)) dir.create(RdataDir)
  
  source('functions_Visium.R')
  library(pryr) # monitor the memory usage
  require(ggplot2)
  options(future.globals.maxSize = 120000 * 1024^2)
  
  mem_used()
  
  st = readRDS(file = paste(resDir, 'RCTD_st_axolotl_allVisium.rds'))
  # st = subset(st, condition == 'Amex_d4')
  refs = readRDS(file = paste0(resDir, '/RCTD_refs_subtypes_final_20221117.rds'))
  
  condition.specific_celltypes = readRDS(file = paste0(resDir, '/RCTD_refs_condition_specificity_v1.rds'))
  
  source('functions_Visium.R')
  
  RCTD_out = paste0(resDir, '/RCTD_allVisium_subtype_out_41subtypes_ref.time.specific_v1.0')
  max_cores = 16
  
  ## preapre the paramters for RCTD subtypes
  DefaultAssay(refs) = 'RNA'
  DefaultAssay(st) = 'Spatial'
  require_int_SpatialRNA = FALSE
  
  condition.specific.ref = TRUE
  
  cc = names(table(st$condition))
  
  source('functions_Visium.R')
  Run.celltype.deconvolution.RCTD(st = st,
                                  refs, 
                                  mapping_ref = 'time',
                                  condition.specific.ref = condition.specific.ref,
                                  condition.specific_celltypes = condition.specific_celltypes,
                                  require_int_SpatialRNA = require_int_SpatialRNA,
                                  max_cores = max_cores,
                                  RCTD_out = RCTD_out,
                                  plot.RCTD.summary = FALSE, 
                                  PLOT.scatterpie = FALSE
                                  
  )
  
  
  source('functions_Visium.R')
  plot.RCTD.results(st = st,
                    RCTD_out = RCTD_out,
                    plot.RCTD.summary = FALSE)
  
  
}


Run_RCTD_adultMice = FALSE
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
  st = readRDS(file = paste0("../results/visium_adultMice_R11934_20210827/Rdata/",
                             'seuratObject_mouse_adult_cell.gene.filtered_umap.clustered.rds'))
  
  #load(file = paste0("../results/visium_adultMice_R11934_20210827/Rdata/", 
  #                   'seuratObject_design_variableGenes_mouse_adult_umap.clustered.Rdata'))
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
  RCTD_out = paste0(outDir, '/RCTD_', length(table(refs$celltype_toUse)), 'Subtype_ref_v1.2')
  
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
