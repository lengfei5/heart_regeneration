##########################################################################
##########################################################################
# Project: heart regeneration
# Script purpose: ligand-receptor anlaysis
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Oct 11 13:43:45 2022
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


library(pryr) # monitor the memory usage
require(ggplot2)
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
library(circlize)
library(RColorBrewer)
require(scran)
require(scater)
source('functions_scRNAseq.R')
source('functions_Visium.R')
dataPath_nichenet = '../data/NicheNet/'

mem_used()

##########################################
# load processed scRNA-seq and visium data
##########################################
# load ST data with additional region annotations
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered_manualSegmentation', 
                   species, '.Rdata')) # visium data

table(st$segmentation, st$condition)

## snRNA-seq reference  
refs = readRDS(file = paste0(RdataDir, 'RCTD_refs_subtypes_final_20221117.rds'))
refs$subtypes = refs$celltype_toUse # clean the special symbols
refs$celltypes = refs$celltype_toUse
#refs_file = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets_20221004_2.rds'
#refs = readRDS(file = refs_file)

table(refs$subtypes)
length(table(refs$subtypes))

# subtype time-specificity 
condition.specific_celltypes = readRDS(paste0(RdataDir, 'RCTD_refs_condition_specificity.rds'))

########################################################
########################################################
# Section : cell neighborhood analysis
# 
########################################################
########################################################
Run_Neighborhood_Enrichment_Analysis = FALSE
if(Run_Neighborhood_Enrichment_Analysis){
  source('functions_Visium.R')
  outDir = paste0(resDir, '/neighborhood_test/')
  RCTD_out = '../results/visium_axolotl_R12830_resequenced_20220308/RCTD_subtype_out_v3.5'
  
  run_neighborhood_analysis(st, 
                            outDir = outDir,
                            RCTD_out = RCTD_out)
  
  
}

########################################################
########################################################
# Section : ligand-receptor prediction analysis with 
# LIANA and NicheNet
# time-specifc and space-specific niches for nichenet
########################################################
########################################################

##########################################
# specific sub-populations to compare
# sender cells, receiver cells
# BZ-specific and Remote-specific populations (Nichenet specific)
##########################################
timepoint_specific = TRUE
# define a list of cell type for each time point, either manual defined or from neighborhood enrichment analysis
#celltypes = c('EC', 'EC_NOS3', 'EC_IS_IARS1', 'FB_IS_TFPI2', 'Mo.Macs_SNX22', 'Neu_IL1R1', 
#              'CM_IS', "RBC")

refs$celltypes = gsub('CM_ven_Robo2', 'CM_Robo2', refs$celltypes)

celltypes_BZ_timeSepcific = list(day1 = c('EC', 'EC_NOS3', 'EC_IS_IARS1', 'FB_IS_TFPI2', 'Mo.Macs_SNX22', 'Neu_IL1R1', 
                                       'CM_IS', "RBC"),
                              day4 = c('EC_IS_LOX', 'EC_IS_Prol', 'Mo.Macs_SNX22', 'Neu_DYSF', 'CM_IS', 
                                       'CM_Prol_IS', 'RBC'),
                              day7 = c('EC_IS_LOX', 'EC_IS_Prol', 'Mo.Macs_FAXDC2', 'Neu_DYSF', 'Neu_IL1R1', 'CM_IS', 
                                       'CM_Prol_IS', 'RBC'),
                              day14 = c('EC_IS_LOX', 'EC_IS_Prol', 'FB_PKD1', 'Neu_DYSF', 'CM_IS', 'Megakeryocytes', 
                                        'RBC')
)

celltypes_RZ_timeSepcific = list(day1 = c('EC', 'EC_NOS3', 'FB_PKD1', 'FB_IS_TFPI2', 'Mo.Macs_SNX22', 
                                          'CM_Robo2', 'CM_IS'),
                                 day4 = c('EC', 'EC_NOS3', 'EC_WNT4', 'FB_PKD1', 'CM_Robo2'),
                                 dya7 = c('EC', 'EC_NOS3', 'EC_IS_Prol', 'FB_PKD1', 'Neu_IL1R1', 'CM_Robo2'),
                                 day14 = c('EC', 'EC_NOS3', 'EC_IS_Prol', 'FB_PKD1', 'CM_Robo2')
                                 )

receiver_cells_BZ = 'CM_IS'
receiver_cells_RZ = 'CM_Robo2'

##########################################
# run LIANA 
##########################################
# set parameter for ligand-receptor analysis
outDir = paste0(resDir, '/Ligand_Receptor_analysis/LIANA_v3.2_41subtypes_receiverCells.CM_IS')
ntop = 100

source('functions_cccInference.R')
run_LIANA(refs, 
          timepoint_specific = timepoint_specific,
          celltypes_timeSpecific = celltypes_BZ_timeSepcific,
          receiver_cells = receiver_cells_BZ,
          outDir = outDir, 
          ntop = ntop)

########################################################
# diff Nichenet for ligand-receptor analysis
# original code from https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md
########################################################
outDir = paste0(resDir, '/Ligand_Receptor_analysis/DiffNicheNet_v3')
system(paste0('mkdir -p ', outDir))

source('functions_cccInference.R')

run_Diff_NicheNet(refs = refs, 
                  timepoint_specific = timepoint_specific,
                  celltypes_BZ_timeSepcific = celltypes_BZ_timeSepcific,
                  celltypes_RZ_timeSepcific = celltypes_RZ_timeSepcific,
                  receiver_cells_BZ = receiver_cells_BZ,
                  receiver_cells_RZ = receiver_cells_RZ, 
                  outDir = outDir
                  )

