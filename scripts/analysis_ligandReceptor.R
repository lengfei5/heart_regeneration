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
version.analysis = '_R17246_R12830_allVisium_20240905'
dataDir = '../R17246_visium_axolotl/nf_out'

resDir = paste0("../results/visium_axolotl", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

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
library(nichenetr)
library(tidyverse)
library(circlize)
library(randomcoloR)

source('functions_scRNAseq.R')
source('functions_Visium.R')

options(future.globals.maxSize = 120000 * 1024^2)


dataPath_nichenet = '../data/NicheNet/'

mem_used()

##########################################
# load processed scRNA-seq and visium data
##########################################
# load ST data with additional region annotations
st = readRDS(file = paste0(RdataDir, 'seuratObject_allVisiusmst_',
              'filtered.spots_time_conditions_manualSegmentation', 
              version.analysis, '.rds'))

table(st$segmentation, st$condition)

SpatialDimPlot(st, ncol = 4)

## snRNA-seq reference  
refs = readRDS(file = paste0(resDir, '/RCTD_refs_subtypes_final_20221117.rds'))
refs$subtypes = refs$celltype_toUse # clean the special symbols
refs$celltypes = refs$celltype_toUse


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
  
  outDir = paste0(resDir, '/neighborhood_test/Run_misty_v1.0/')
  #condition.specific_celltypes = readRDS(file = 
  #                                         paste0(resDir, '/RCTD_refs_condition_specificity_v3.rds'))
  
  
  RCTD_out = paste0(resDir, '/RCTD_allVisium_subtype_out_41subtypes_ref.time.specific_v2.0')
  
  levels(refs$subtypes)
  # condition-specific subtypes selected
  #condSpec_celltypes = readxl::read_xlsx("../data/neighbourhood_analysis_list_short.xlsx")
  #condSpec_celltypes = as.data.frame(condSpec_celltypes)
  
  condSpec_celltypes = list(d1 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", 'EC_IS_IARS1',
                                     "FB_TNXB",'FB_IS_TFPI2',
                                     'Mo.Macs_SNX22', "Neu_DYSF",
                                      "CM_Cav3.1", "CM_Robo2", 'CM_IS',
                                     "Megakeryocytes","RBC"),
                            
                            d4 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", 'EC_IS_IARS1',
  "EC_IS_LOX",
                                     "FB_PKD1", "FB_TNXB",
                                     "Mo.Macs_resident",  "Mo.Macs_FAXDC2", 'Mo.Macs_SNX22', 'Neu_DYSF',
                                     "CM_Cav3.1", "CM_Robo2", 'CM_IS',  'CM_Prol_IS',
                                     "Megakeryocytes", 'RBC'),
                            d7 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", "EC_IS_LOX",
                                     "FB_PKD1", "FB_TNXB", "FB_VWA2",
                                     "Mo.Macs_resident", "Mo.Macs_FAXDC2", 'Neu_DYSF',
                                     "CM_Robo2", "CM_Cav3.1",  'CM_IS',  'CM_Prol_IS',
                                     "Megakeryocytes", 'RBC'),
                            d14 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", "EC_IS_LOX",
                                      "FB_PKD1", "FB_TNXB", "FB_VWA2",
                                      "Mo.Macs_resident", 'Neu_DYSF',
                                      "CM_Robo2", "CM_Cav3.1",  'CM_IS',
                                      "Megakeryocytes", 'RBC')
  )
  
  celltypes_interest = c()
  for(n in 1:length(condSpec_celltypes))
  {
    celltypes_interest = unique(c(celltypes_interest, condSpec_celltypes[[n]]))
  }
  
  condSpec_celltypes$d0 = celltypes_interest
  
  
  st$segmentation = as.character(st$segmentation)
  st$segmentation[which(st$segmentation == 'RZ1')] = 'RZ'
  st$segmentation[which(st$segmentation == 'RZ2')] = 'RZ'
  st$segmentation[which(st$segmentation == 'Intact1')] = 'Intact'
  st$segmentation[which(st$segmentation == 'Intact2')] = 'Intact'
  
  st$segmentation = as.factor(st$segmentation)
  SpatialDimPlot(st, group.by = 'segmentation', ncol = 4)
  table(st$segmentation, st$condition)
  
  
  source('functions_Visium.R')
  run_misty_colocalization_analysis(st, 
                                    outDir = outDir,
                                    RCTD_out = RCTD_out,
                                    condSpec_celltypes = NULL,
                                    segmentation_annots = c('all', 'BZ', 'RZ', 'Intact')
                                    )
  
  
  
  source('functions_Visium.R')
  
  run_significanceTest_misty(st, 
                             outDir = outDir, 
                             segmentation_annots = c('all', 'BZ', 'RZ', 'Intact'))
  
  
  
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
# define a list of cell type for each time point, either manual defined or from neighborhood enrichment analysis
#celltypes = c('EC', 'EC_NOS3', 'EC_IS_IARS1', 'FB_IS_TFPI2', 'Mo.Macs_SNX22', 'Neu_IL1R1', 
#              'CM_IS', "RBC")

version_testing_short = FALSE
if(version_testing){
  timepoint_specific = TRUE
  
  celltypes_BZ_timeSpecific = list(day1 = c('EC', 'EC_NOS3', 'EC_IS_IARS1', 'FB_IS_TFPI2', 'Mo.Macs_SNX22',
  'Neu_IL1R1',
                                         'CM_IS', "RBC"),
                                day4 = c('EC_IS_LOX', 'EC_IS_Prol', 'Mo.Macs_SNX22', 'Neu_DYSF', 'CM_IS',
                                         'CM_Prol_IS', 'RBC'),
                                day7 = c('EC_IS_LOX', 'EC_IS_Prol', 'Mo.Macs_FAXDC2', 'Neu_DYSF', 'Neu_IL1R1',
  'CM_IS',
                                         'CM_Prol_IS', 'RBC'),
                                day14 = c('EC_IS_LOX', 'EC_IS_Prol', 'FB_PKD1', 'Neu_DYSF', 'CM_IS', 'Megakeryocytes',
                                          'RBC')
  )
  celltypes_RZ_timeSpecific = list(day1 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2'),
                                   day4 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2'),
                                   day7 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2'),
                                   day14 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2')
  )
  receivers_BZ_timeSpecific = list(day1 = c("CM_IS"),
                                   day4 = c('CM_Prol_IS', "CM_IS"),
                                   day7 = c('CM_Prol_IS', "CM_IS"),
                                   day14 = c('CM_Prol_IS', "CM_IS")
  )
  
  receivers_RZ_timeSpecific = list(day1 = c("CM_Robo2"),
                                   day4 = c("CM_Robo2"),
                                   day7 = c("CM_Robo2"),
                                   day14 = c("CM_Robo2")
  )
  
}


version_testing_long = FALSE
if(version_testing_long){
  timepoint_specific = TRUE
  
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
  
  celltypes_RZ_timeSpecific = list(day1 = c('EC', 'EC_NOS3', "FB_TNXB", 'Mo.Macs_SNX22',
                                            "CM_Cav3.1", "CM_Robo2", 'RBC'),
                                   day4 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2', 'CM_Prol_IS', "RBC"),
                                   day7 = c('EC', 'EC_NOS3', 'FB_PKD1', "FB_TNXB", 'CM_Robo2'),
                                   day14 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2', "RBC")
  )
  
  
  receivers_BZ_timeSpecific = list(day1 = c("CM_IS"),
                                   day4 = c('CM_Prol_IS', "CM_IS"),
                                   day7 = c('CM_Prol_IS', "CM_IS"),
                                   day14 = c('CM_Prol_IS', "CM_IS")
  )
  
  receivers_RZ_timeSpecific = list(day1 = c("CM_Robo2"),
                                   day4 = c("CM_Robo2"),
                                   day7 = c("CM_Robo2"),
                                   day14 = c("CM_Robo2")
  )
  
}

### version for Prateek
OnlyFB_version_for_Prateek = FALSE
if(OnlyFB_version_for_Prateek)
{
  timepoint_specific = TRUE
  celltypes_BZ_timeSpecific = list(day1 = c('FB_TNXB'),
                                   day4 = c('FB_TNXB', 'FB_PKD1'),
                                   day7 = c('FB_TNXB', 'FB_PKD1'),
                                   day14 = c('FB_TNXB', 'FB_PKD1')
  )
  receivers_BZ_timeSpecific = list(day1 = c("CM_IS"),
                                   day4 = c('CM_Prol_IS'),
                                   day7 = c('CM_Prol_IS'),
                                   day14 = c("CM_IS")
  )
  
}


### all pairs of subtypes
version_testing_all.subtype.pairs = FALSE
if(version_testing_all.subtype.pairs){
  
  # define cell subtype pairs in the border zone 
  celltypes_BZ_timeSpecific = list(day1 = list(CM_IS = c("CM_Cav3.1", "CM_Robo2", "EC_WNT4", 
                                                         'Mo.Macs_SNX22','Neu_IL1R1', "RBC")
                                               #CM_Prol_IS = c("CM_Cav3.1", "CM_Robo2", "EC_WNT4", 
                                               #           'Mo.Macs_SNX22','Neu_IL1R1', "RBC")
                                               ),
                                   
                                   day4 = list(CM_Prol_IS = c('CM_IS',"CM_Robo2", 'EC', 'FB_TNXB',  
                                                              'Mo.Macs_SNX22', "Mo.Macs_FAXDC2")
                                               
                                               ),
                                            
                                   day7 = list(CM_Prol_IS = c("CM_Robo2",'EC_IS_LOX', "EC_NOS3",
                                                              "FB_PKD1", "FB_TNXB",
                                                              'Mo.Macs_FAXDC2','Neu_DYSF','RBC')),
                                   
                                   day14 = list(CM_IS = c("CM_Prol_3", "CM_Robo2", 'EC_IS_LOX', 
                                                          'FB_PKD1', "FB_TNXB", "Mo.Macs_resident", 'RBC'))
                                   
  )
  
  # define cell subtype pairs in the remote zone as control for NicheNet
  celltypes_RZ_timeSpecific = list(day1 = list(CM_Robo2 = c("CM_Cav3.1", "EC_IS_IARS1", "EC_WNT4")),
                                   day4 = list(CM_Robo2 = c('CM_Robo2',"FB_TNXB", "Megakeryocytes",
                                                            'Mo.Macs_resident', 'RBC')),
                                   day7 = list(CM_Robo2 = c("EC_NOS3", "EC_WNT4", "Mo.Macs_resident", "RBC")),
                                   day14 = list(CM_Robo2 = c("CM_Cav3.1", "EC_IS_IARS1", "FB_TNXB", 
                                                             "Mo.Macs_resident"))
  )
  
}


##########################################
# run LIANA for all pairs
##########################################
# set parameter for ligand-receptor analysis
outDir_version = paste0(resDir, '/Ligand_Receptor_analysis/LIANA_v5.5_FB.Immue_intraOnly')
if(!dir.exists(outDir_version)) dir.create(outDir_version)

out_misty = paste0('../results/visium_axolotl_R12830_resequenced_20220308/neighborhood_test/',
                   'Run_misty_v1.8_short/Plots_RCTD_density')

misty_cutoff = 0.2

# run LIANA day by day
timepoint_specific = TRUE

#times_slice = c('d1', 'd4', 'd7', 'd14')
times_slice = c('d7')

subtypes = unique(refs$celltypes)

for(n in 1:length(times_slice))
{
  # n = 1
  source('functions_cccInference.R')
  time = times_slice[n]
  cat(' run LIANA for time -- ', time, '\n')
  
  outDir = paste(outDir_version, '/', time, collapse = '')
  outDir = gsub(' ', '', outDir)
  
  ## select the interacting subtype pairs  
  intra = read.csv2(paste0(out_misty, '/Amex_', time, '_all_summary_table_intra.csv'), row.names = 1)
  juxta = read.csv2(paste0(out_misty, '/Amex_', time, '_all_summary_table_juxta5.csv'), row.names = 1)
    
  intra = intra > misty_cutoff
  juxta = juxta > misty_cutoff
  
  intra[which(is.na(intra))] = FALSE
  juxta[which(is.na(juxta))] = FALSE
  
  pairs = intra + juxta > 0
  pairs = pairs[which(rownames(pairs) != 'RBC'), which(colnames(pairs) != 'RBC')]
  ss_row = apply(pairs, 1, sum)
  ss_col = apply(pairs, 2, sum)
  
  Select_specificPairs = TRUE
  if(Select_specificPairs){
    
    FB_Immune = c('FB.TNXB', 'FB.PKD1', 'Mo.Macs,FAXDC2', 'Mo.Macs.SNX22', 'Neu.DYSF')
    ii1 = which(!is.na(match(colnames(pairs), FB_Immune)))
    jj1 = which(!is.na(match(rownames(pairs), FB_Immune)))
    
    pairs = pairs[jj1, ii1] 
    
    ss_col = apply(pairs, 2, sum)
    ss_row = apply(pairs, 1, sum)
    
    pairs = pairs[which(ss_row>=1), which(ss_col >= 1)] # at least interacting with 1 receivers
    
    
  }else{
    
    pairs = pairs[ ,which(ss_col >= 1)] # at least interacting with 1 receivers
    pairs = pairs[which(ss_row >= 3), ] # at least have 3 senders 
    
  }
  
  
  colnames(pairs) = gsub("Cav3_1", "Cav3.1", gsub('Mo_Macs', 'Mo.Macs', gsub('[.]','_', colnames(pairs))))
  rownames(pairs) = gsub("Cav3_1", "Cav3.1", gsub('Mo_Macs', 'Mo.Macs', gsub('[.]','_', rownames(pairs))))
  
  cat(match(colnames(pairs), subtypes), '\n')
  cat(match(rownames(pairs), subtypes), '\n')
  
  celltypes_BZ_timeSpecific = vector("list", nrow(pairs))
  
  for(m in 1:nrow(pairs))
  {
    # m = 16
    celltypes_BZ_timeSpecific[[m]] = colnames(pairs)[which(pairs[m, ] == TRUE)]
    #x = pairs[which(rownames(pairs) == 'CM.Prol.IS'), ]
    names(celltypes_BZ_timeSpecific)[m] = rownames(pairs)[m]
    
  }
  
  run_LIANA(refs, 
            timepoint_specific = TRUE,
            include_autocrine = TRUE,
            celltypes_timeSpecific = celltypes_BZ_timeSpecific,
            outDir = outDir
  )
  
  source("functions_cccInference.R")
  res = aggregate_output_LIANA(liana_out = paste(outDir))
  
  require(cellcall)
  
  subtypes = unique(c(res$sender, res$receiver))
  cell_color = data.frame(color = distinctColorPalette(length(subtypes)), stringsAsFactors = FALSE)
  rownames(cell_color) <- subtypes
  
  # ViewInterCircos(object = mt, font = 2, cellColor = cell_color, 
  #                 lrColor = c("#F16B6F", "#84B1ED"),
  #                 arr.type = "big.arrow",arr.length = 0.04,
  #                 trackhight1 = 0.05, 
  #                 slot="expr_l_r_log2_scale",
  #                 linkcolor.from.sender = TRUE,
  #                 linkcolor = NULL, 
  #                 gap.degree = 2,
  #                 order.vector=c('ST', "SSC", "SPGing", "SPGed"),
  #                 trackhight2 = 0.032, 
  #                 track.margin2 = c(0.01,0.12), 
  #                 DIY = FALSE)
  
  ## test celltalker
  library(celltalker)
  #xx = readRDS(file = paste0(resDir, '/test_LianaOut_for_circoPlot.rds'))
  xx = res
  xx$sender = gsub('_', '.', xx$sender)
  xx$receiver = gsub('_', '.', xx$receiver)
  xx$interaction_pairs = paste0(xx$sender, '_', xx$receiver)
  xx$interaction = paste0(xx$ligand, '_', xx$receptor)
  xx = xx[,c(8, 9, 7, 1:5)]
  colnames(xx)[3:5] = c('value', 'cell_type1', 'cell_type2')
  
  top_stats_xx <- xx %>% as_tibble() %>%
    #mutate(fdr=p.adjust(p_val,method="fdr")) %>%
    #filter(fdr<0.05) %>%
    group_by(cell_type2) %>%
    top_n(-50, aggregate_rank) %>%
    ungroup()
  
  cat(nrow(top_stats_xx), ' top pairs \n')
  #top_stats_xx = as_tibble(xx)
  subtypes = unique(c(top_stats_xx$cell_type1, top_stats_xx$cell_type2))
  colors_use <- distinctColorPalette(length(subtypes))
  
  svg(paste0(outDir, "/ligand_target_circos_", time, "_test.svg"), width = 6, height = 6)
  
  source("functions_cccInference.R")
  circos_plot_customized(ligand_receptor_frame=top_stats_xx,
                         cell_group_colors=colors_use,
                         ligand_color="#84B1ED",
                         receptor_color="#F16B6F",
                         cex_outer=0.5,
                         cex_inner=0.2,
                         link.lwd=0.1, 
                         arr.length=0., 
                         arr.width=0.05
  )
  
  dev.off()
  
  
  
}


## double check the ligand and receptor expression distribution
#FeaturePlot(refs, features = rownames(refs)[grep('EGFC|VIPR2', rownames(refs))])

########################################################
# diff Nichenet for ligand-receptor analysis
# original code from https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md
########################################################
outDir_version = paste0(resDir, '/Ligand_Receptor_analysis/DiffNicheNet_v5.1_allpairs_intraOnly')

for(n in 1:length(celltypes_BZ_timeSpecific))
{
  # n = 2
  source('functions_cccInference.R')
  time = names(celltypes_BZ_timeSpecific)[n]
  cat(' run DiffNicheNet for time -- ', time, '\n')
  outDir = paste(outDir_version, '/', time, collapse = '')
  outDir = gsub(' ', '', outDir)
  
  system(paste0('mkdir -p ', outDir))
  
  source('functions_cccInference.R')
  
  run_Diff_NicheNet(refs = refs, 
                    timepoint_specific = TRUE,
                    include_autocrine = TRUE,
                    celltypes_BZ_specificDay = celltypes_BZ_timeSpecific[[n]],
                    celltypes_RZ_specificDay = celltypes_RZ_timeSpecific[[n]],
                    outDir = outDir
  )
  
  # extract_tables_from_res_Diff_NicheNet(outDir)
  
}

##########################################
# ## quick construction of GAS6-AXL to targets 
##########################################
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