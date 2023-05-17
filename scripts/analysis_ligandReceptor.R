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

refs$subtypes = as.factor(refs$subtypes) 

refs$celltypes = gsub('CM_ven_Robo2', 'CM_Robo2', refs$celltypes)
refs$celltypes = gsub('CM_ven_Cav3_1', 'CM_Cav3.1', refs$celltypes)


########################################################
########################################################
# Section : cell neighborhood analysis
# 
########################################################
########################################################
Run_Neighborhood_Enrichment_Analysis = FALSE
if(Run_Neighborhood_Enrichment_Analysis){
  
  outDir = paste0(resDir, '/neighborhood_test/Run_misty_v1.8_test/')
  
  RCTD_out = paste0('../results/visium_axolotl_R12830_resequenced_20220308/',
                    'RCTD_subtype_out_42subtypes_ref.time.specific_v4.3')
  
  levels(refs$subtypes)
  # condition-specific subtypes selected
  #condSpec_celltypes = readxl::read_xlsx("../data/neighbourhood_analysis_list_short.xlsx")
  #condSpec_celltypes = as.data.frame(condSpec_celltypes)
  # condSpec_celltypes = list(d1 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", 'EC_IS_IARS1', 
  #                                    "FB_TNXB",'FB_IS_TFPI2',
  #                                    'Mo.Macs_SNX22', "Neu_DYSF", 
  #                                     "CM_Cav3.1", "CM_Robo2", 'CM_IS',
  #                                    "Megakeryocytes","RBC"),
  #                           d4 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", 'EC_IS_IARS1', 
  # "EC_IS_LOX",
  #                                    "FB_PKD1", "FB_TNXB", 
  #                                    "Mo.Macs_resident",  "Mo.Macs_FAXDC2", 'Mo.Macs_SNX22', 'Neu_DYSF', 
  #                                    "CM_Cav3.1", "CM_Robo2", 'CM_IS',  'CM_Prol_IS', 
  #                                    "Megakeryocytes", 'RBC'),
  #                           d7 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", "EC_IS_LOX",
  #                                    "FB_PKD1", "FB_TNXB", "FB_VWA2", 
  #                                    "Mo.Macs_resident", "Mo.Macs_FAXDC2", 'Neu_DYSF', 
  #                                    "CM_Robo2", "CM_Cav3.1",  'CM_IS',  'CM_Prol_IS', 
  #                                    "Megakeryocytes", 'RBC'),
  #                           d14 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", "EC_IS_LOX",
  #                                     "FB_PKD1", "FB_TNXB", "FB_VWA2", 
  #                                     "Mo.Macs_resident", 'Neu_DYSF', 
  #                                     "CM_Robo2", "CM_Cav3.1",  'CM_IS',   
  #                                     "Megakeryocytes", 'RBC')
  # )
  # 
  condSpec_celltypes = list(d1 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", 'EC_IS_IARS1', "EC_Prol",
                                   "FB_PKD1", "FB_TNXB",'FB_IS_TFPI2',
                                   'Mo.Macs_SNX22', "Neu_DYSF", "Neu_IL1R1", 
                                   "CM_Robo2", "CM_Cav3.1",  'CM_IS', "CM_Prol_1", "CM_Prol_3",
                                   "Megakeryocytes", "Proliferating_Megakeryocytes", "RBC", "Proliferating_RBC"),
                            
                            d4 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", 'EC_IS_IARS1', "EC_IS_LOX", 
                                   "EC_IS_Prol", "EC_Prol",
                                   "FB_PKD1", "FB_TNXB",
                                   "Mo.Macs_Prol", "Mo.Macs_resident", "Mo.Macs_FAXDC2", 'Mo.Macs_SNX22', "Neu_DYSF", 
                                   "Neu_IL1R1",
                                   "CM_Robo2", "CM_Cav3.1",  'CM_IS', "CM_Prol_IS", "CM_Prol_1", "CM_Prol_3",
                                   "Megakeryocytes", "Proliferating_Megakeryocytes", "RBC", "Proliferating_RBC"),
                            
                            d7 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", "EC_IS_LOX", "EC_IS_Prol", 
                                   "EC_Prol",
                                   "FB_PKD1", "FB_TNXB",
                                   "Mo.Macs_Prol", "Mo.Macs_resident", "Mo.Macs_FAXDC2", "Neu_DYSF", "Neu_IL1R1",
                                   "CM_Robo2", "CM_Cav3.1",  'CM_IS', "CM_Prol_IS", "CM_Prol_1", "CM_Prol_3",
                                   "Megakeryocytes", "Proliferating_Megakeryocytes", "RBC", "Proliferating_RBC"),
                            
                            d14 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", "EC_IS_LOX", "EC_IS_Prol",
                                    "EC_Prol",
                                    "FB_PKD1", "FB_TNXB",
                                    "Mo.Macs_Prol", "Mo.Macs_resident", "Neu_DYSF", 
                                    "CM_Robo2", "CM_Cav3.1",  'CM_IS', "CM_Prol_1", "CM_Prol_3",
                                    "Megakeryocytes", "Proliferating_Megakeryocytes", "RBC", "Proliferating_RBC")
  )
  
  
  source('functions_Visium.R')
  # run_neighborhood_analysis(st, 
  #                           outDir = outDir,
  #                           RCTD_out = RCTD_out)
  run_misty_colocalization_analysis(st, 
                                    outDir = outDir,
                                    RCTD_out = RCTD_out,
                                    condSpec_celltypes = condSpec_celltypes
                                    )
  
  
}

########################################################
########################################################
# Section : ligand-receptor prediction analysis with 
# LIANA and NicheNet
# time-specifc and space-specific niches for nichenet
########################################################
########################################################

# subtype time-specificity 
# condition.specific_celltypes = readRDS(paste0(RdataDir, 'RCTD_refs_condition_specificity.rds'))

##########################################
# specific sub-populations to compare
# sender cells, receiver cells
# BZ-specific and Remote-specific populations (Nichenet specific)
##########################################
timepoint_specific = TRUE

# define a list of cell type for each time point, either manual defined or from neighborhood enrichment analysis
#celltypes = c('EC', 'EC_NOS3', 'EC_IS_IARS1', 'FB_IS_TFPI2', 'Mo.Macs_SNX22', 'Neu_IL1R1', 
#              'CM_IS', "RBC")



# celltypes_BZ_timeSpecific = list(day1 = c('EC', 'EC_NOS3', 'EC_IS_IARS1', 'FB_IS_TFPI2', 'Mo.Macs_SNX22', 
# 'Neu_IL1R1', 
#                                        'CM_IS', "RBC"),
#                               day4 = c('EC_IS_LOX', 'EC_IS_Prol', 'Mo.Macs_SNX22', 'Neu_DYSF', 'CM_IS', 
#                                        'CM_Prol_IS', 'RBC'),
#                               day7 = c('EC_IS_LOX', 'EC_IS_Prol', 'Mo.Macs_FAXDC2', 'Neu_DYSF', 'Neu_IL1R1', 
# 'CM_IS', 
#                                        'CM_Prol_IS', 'RBC'),
#                               day14 = c('EC_IS_LOX', 'EC_IS_Prol', 'FB_PKD1', 'Neu_DYSF', 'CM_IS', 'Megakeryocytes', 
#                                         'RBC')
# )

celltypes_BZ_timeSpecific = list(day1 = c('EC', 'EC_NOS3', "FB_TNXB", 
                                          "CM_Cav3.1", "CM_Robo2", "CM_IS", 
                                          "Megakeryocytes", "RBC"),
                                 day4 = c('EC', 'EC_IS_LOX', "FB_PKD1", 
                                          "Mo.Macs_FAXDC2", 'Mo.Macs_SNX22', 
                                          'Neu_DYSF', "CM_Robo2", 'CM_Prol_IS', "CM_IS",
                                          "Megakeryocytes" ,'RBC'),
                                 day7 = c('EC_NOS3', 'EC_IS_LOX', "FB_PKD1","FB_TNXB", 
                                          "Mo.Macs_FAXDC2", 'Neu_DYSF',
                                          "CM_Robo2", 'CM_Prol_IS', "CM_IS",
                                          "Megakeryocytes" ,'RBC'),
                                 day14 = c('EC_NOS3', "FB_PKD1", "FB_TNXB", 
                                           "Mo.Macs_resident", 'Neu_DYSF', "CM_IS", 'CM_Prol_IS',
                                           "CM_Robo2", 'RBC')
)

# celltypes_RZ_timeSpecific = list(day1 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2'),
#                                  day4 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2'),
#                                  dya7 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2'),
#                                  day14 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2')
# )

celltypes_RZ_timeSpecific = list(day1 = c('EC', 'EC_NOS3', "FB_TNXB", 'Mo.Macs_SNX22',
                                          "CM_Cav3.1", "CM_Robo2", 'RBC'),
                                 day4 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2', 'CM_Prol_IS', "RBC"),
                                 dya7 = c('EC', 'EC_NOS3', 'FB_PKD1', "FB_TNXB", 'CM_Robo2'),
                                 day14 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2', "RBC")
                                 )

receivers_BZ_timeSpecific = list(day1 = c("CM_IS"),
                         day4 = c('CM_Prol_IS', "CM_IS"),
                         dya7 = c('CM_Prol_IS', "CM_IS"),
                         day14 = c('CM_Prol_IS', "CM_IS")
                         )

receivers_RZ_timeSpecific = list(day1 = c("CM_Robo2"),
                         day4 = c("CM_Robo2"),
                         dya7 = c("CM_Robo2"),
                         day14 = c("CM_Robo2")
                         )

##########################################
# run LIANA 
##########################################
# set parameter for ligand-receptor analysis
outDir = paste0(resDir, '/Ligand_Receptor_analysis/LIANA_v4.0_receiverCells_CM.IS.CM.IS.PROL')

source('functions_cccInference.R')
run_LIANA(refs, 
          timepoint_specific = timepoint_specific,
          celltypes_timeSpecific = celltypes_BZ_timeSpecific,
          receiver_cells = receivers_BZ_timeSpecific,
          outDir = outDir)

## double check the ligand and receptor expression distribution
FeaturePlot(refs, features = rownames(refs)[grep('EGFC|VIPR2', rownames(refs))])

########################################################
# diff Nichenet for ligand-receptor analysis
# original code from https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md
########################################################
outDir = paste0(resDir, '/Ligand_Receptor_analysis/DiffNicheNet_v4.0_test')
system(paste0('mkdir -p ', outDir))

source('functions_cccInference.R')

run_Diff_NicheNet(refs = refs, 
                  timepoint_specific = timepoint_specific,
                  celltypes_BZ_timeSpecific = celltypes_BZ_timeSpecific,
                  celltypes_RZ_timeSpecific = celltypes_RZ_timeSpecific,
                  receivers_BZ_timeSpecific = receivers_BZ_timeSpecific,
                  receivers_RZ_timeSpecific = receivers_RZ_timeSpecific, 
                  outDir = outDir
                  )

extract_tables_from_res_Diff_NicheNet(outDir)


## quick construction of GAS6-AXL to targets 
library(nichenetr)
library(tidyverse)

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
ligand_tf_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_tf_matrix.rds"))

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
sig_network = readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
gr_network = readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))


ligands_all = "GAS6" # this can be a list of multiple ligands if required
targets_all = unique(table_targets$target[which(table_targets$ligand == ligands_all)])

targets_all = c("EGFR", "ETS1","LEF1", 'MYH9', 'SMAD3', 'STAT1')


active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, 
                                                     ligands_all = ligands_all, 
                                                     targets_all = targets_all, 
                                                     weighted_networks = weighted_networks)

# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and 
# gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% 
  mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% 
  mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, 
                                                  ligands_all = ligands_all, 
                                                  targets_all = targets_all, 
                                                  sig_color = "indianred", 
                                                  gr_color = "steelblue")

# To render the graph: uncomment following line of code
DiagrammeR::render_graph(graph_min_max, layout = "kk")


data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network, 
                                                   lr_network = lr_network, 
                                                   sig_network = sig_network, 
                                                   gr_network = gr_network)
head(data_source_network) 