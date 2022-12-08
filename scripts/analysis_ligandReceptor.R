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

# subtype time-specificity 
condition.specific_celltypes = readRDS(paste0(RdataDir, 'RCTD_refs_condition_specificity.rds'))

########################################################
########################################################
# Section : cell neighborhood analysis
# 
########################################################
########################################################
source('functions_Visium.R')
outDir = paste0(resDir, '/neighborhood_test/')
RCTD_out = '../results/visium_axolotl_R12830_resequenced_20220308/RCTD_subtype_out_v3.5'

run_neighborhood_analysis(st, 
                          outDir = outDir,
                          RCTD_out = RCTD_out)

########################################################
########################################################
# Section : ligand-receptor prediction analysis with LIANA 
# and NicheNet
########################################################
########################################################

##########################################
# run LIANA 
##########################################
# set parameter for ligand-receptor analysis
outDir = paste0(resDir, '/Ligand_Receptor_analysis/LIANA_v3_42subtypes_receiverCells.CM_IS')



timepoint_specific = TRUE
# define a list of cell type for each time point, either manual defined or from neighborhood enrichment analysis
celltypes = c('EC', 'EC_NOS3', 'EC_IS_IARS1', 'FB_IS_TFPI2', 'Mo.Macs_SNX22', 'Neu_IL1R1', 
              'CM_IS', "RBC")
celltypes_timeSpecific = list(day1 = c('EC', 'EC_NOS3', 'EC_IS_IARS1', 'FB_IS_TFPI2', 'Mo.Macs_SNX22', 'Neu_IL1R1', 
                                   'CM_IS', "RBC"),
                          day4 = c('EC_IS_LOX', 'EC_IS_Prol', 'Mo.Macs_SNX22', 'Neu_DYSF', 'CM_IS', 
                                   'CM_Prol_IS', 'RBC'),
                          day7 = c('EC_IS_LOX', 'EC_IS_Prol', 'Mo.Macs_FAXDC2', 'Neu_DYSF', 'Neu_IL1R1', 'CM_IS', 
                                   'CM_Prol_IS', 'RBC'),
                          day14 = c('EC_IS_LOX', 'EC_IS_Prol', 'FB_PKD1', 'Neu_DYSF', 'CM_IS', 'Megakeryocytes', 
                                    'RBC')
                          )
ntop = 100

source('functions_cccInference.R')

run_LIANA(refs, 
          celltypes_timeSpecific = celltypes_timeSpecific,
          timepoint_specific = timepoint_specific,
          receiver_cells = 'CM_IS',
          outDir = outDir, 
          ntop = ntop)


########################################################
########################################################
# Section : Nichenet for ligand-receptor analysis
# original code from https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md
########################################################
########################################################
## snRNA-seq as reference
refs_file = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets_20221004_2.rds'
refs = readRDS(file = refs_file)
table(refs$subtypes)

outDir = paste0(resDir, '/Ligand_Receptor_analysis/NicheNet_v2')
system(paste0('mkdir -p ', outDir))

##########################################
# preprocess seurat data with selected cell types
##########################################
celltypes_BZ = c('Mono_Macrophages', 'Proliferating_CM', 'Injury_specific_EC', 'FB_1')
celltypes_sel = c('CM', 'EC', 'FB', 'MP')


## select first CM, EC, FB, MP cell types
refs$celltypes = as.character(refs$subtypes)
refs$celltypes[grep('CM_|CMs_|_CM|_CM_', refs$subtypes)] = 'CM'
refs$celltypes[grep('EC_|_EC', refs$subtypes)] = 'EC'
refs$celltypes[grep('FB_', refs$subtypes)] = 'FB'
refs$celltypes[grep('Macrophages|_MF', refs$subtypes)] = 'MP'

table(refs$celltypes)

refs$celltypes[which(refs$celltypes == 'CM' & refs$subtypes != 'Proliferating_CM' & 
                       refs$subtypes != 'Ventricular_CM_ROBO2+')] = NA
table(refs$celltypes)

refs$celltypes[grep('CM|EC|FB|MP', refs$celltypes, invert = TRUE)] = NA
table(refs$celltypes)

## subset seurat object and change names
Idents(refs) = as.factor(refs$celltypes)
subref = subset(refs, cells = colnames(refs)[!is.na(match(refs$celltypes, celltypes_sel))])
subref$celltypes = droplevels(as.factor(subref$celltypes))
table(subref$celltypes)

## define BZ subpopution and non-BZ ones
subref$celltypes = as.character(subref$celltypes)
subref$celltypes[which(subref$subtypes == 'Mono_Macrophages')] = 'MP_BZ'
subref$celltypes[which(subref$subtypes == 'Proliferating_CM')] = 'CM_BZ'
subref$celltypes[which(subref$subtypes == 'Injury_specific_EC')] = 'EC_BZ'
subref$celltypes[which(subref$subtypes == 'FB_1')] = 'FB_BZ'


cat('celltype to consider -- ', names(table(subref$celltypes)), '\n')
Idents(subref) = as.factor(subref$celltypes)
table(subref$celltypes)

subref = subset(x = subref, downsample = 3000) # downsample the CM and EC for the sake of speed
table(subref$celltypes)

##########################################
# Gene filtering, in particular filtering gene without UMI counts or lowly expressed genes
##########################################
Gene.filtering.preprocessing = TRUE
if(Gene.filtering.preprocessing){
  sce <- as.SingleCellExperiment(subref)
  ggs = rownames(sce)
  
  colLabels(sce) = as.factor(sce$celltypes)
  rownames(sce) = get_geneName(rownames(sce))
  
  ave.counts <- calculateAverage(sce, assay.type = "counts")
  
  hist(log10(ave.counts), breaks=100, main="", col="grey80",
       xlab=expression(Log[10]~"average count"))
  
  num.cells <- nexprs(sce, byrow=TRUE)
  smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
                xlab=expression(Log[10]~"average count"))
  
  # detected in >= 10 cells, ave.counts cutoff are both very lose
  genes.to.keep <- num.cells > 10 & ave.counts >= 10^-5  & ave.counts <10^3  
  summary(genes.to.keep)
  sce <- sce[genes.to.keep, ]
  ggs = ggs[genes.to.keep]
  
  rownames(sce) = make.names(rownames(sce), unique = TRUE)
  
  geneDup = data.frame(gene = ggs, gg.uniq = rownames(sce), stringsAsFactors = FALSE)
  geneDup$geneSymbol = get_geneName(geneDup$gene)
  
  kk = which(geneDup$geneSymbol != geneDup$gg.uniq)
  
  subref = as.Seurat(sce, counts = "counts", data = "logcounts")
  rm(sce)
  Idents(subref) = as.factor(subref$celltypes)
  
  saveRDS(geneDup, paste0(outDir, '/geneSymbol_duplication_inLRanalysis.rds'))
  
}

##########################################
# overview of scRNA-seq data 
##########################################
rm(refs)
subref$celltype = subref$celltypes

# rerun the processing and umap
subref <- FindVariableFeatures(subref, selection.method = "vst", nfeatures = 3000)
subref <- ScaleData(subref)

subref <- RunPCA(subref, verbose = FALSE, features = VariableFeatures(object = subref), weight.by.var = TRUE,
                 reduction.key = 'pca_')
ElbowPlot(subref, ndims = 30)

subref <- RunUMAP(subref, dims = 1:30, n.neighbors = 30, min.dist = 0.1, reduction = 'pca', 
                  reduction.key = 'umap_')
DimPlot(subref, group.by = "celltype", label = TRUE, reduction = 'umap') 


geneDup = readRDS(paste0(outDir, '/geneSymbol_duplication_inLRanalysis.rds'))
kk = which(geneDup$geneSymbol != geneDup$gg.uniq)
ggs = unique(geneDup$geneSymbol[kk])
cat(length(ggs), ' genes have symbols issue \n')
ligand_dup = ggs[ggs %in% lr_network$ligand]
receptor_dup = ggs[ggs %in% lr_network$receptor]

ligand_dup = geneDup[!is.na(match(geneDup$geneSymbol, ligand_dup)), ]
ligand_dup = ligand_dup[order(ligand_dup$geneSymbol), ]

receptor_dup = geneDup[!is.na(match(geneDup$geneSymbol, receptor_dup)), ]
receptor_dup = receptor_dup[order(receptor_dup$geneSymbol), ]

write.csv2(ligand_dup, file = paste0(outDir, '/geneDuplication_issue_ligands.csv'))
write.csv2(receptor_dup, file = paste0(outDir, '/geneDuplication_issue_receptors.csv'))

saveRDS(subref, file = paste0(outDir, '/seuratObject_snRNAseq_subset_for_NicheNet.rds')) 

##########################################
# # reload the processed seurat object
##########################################
subref = readRDS(file = paste0(outDir, '/seuratObject_snRNAseq_subset_for_NicheNet.rds'))

#subref = subset(x = subref, downsample = 1000) 
table(subref$celltypes)

seurat_obj = SetIdent(subref, value = "celltype")
DimPlot(seurat_obj, group.by = "celltype", label = TRUE, reduction = 'umap') 

table(seurat_obj$celltype)

##########################################
# step 0):  load the Nichenet data 
##########################################
# NicheNet’s ligand-target prior model
ligand_target_matrix = readRDS(paste0(dataPath_nichenet,  "ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

# ligand-receptor network, Putative ligand-receptor links from NicheNet
lr_network = readRDS(paste0(dataPath_nichenet, "lr_network.rds"))
lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% 
  distinct(ligand, receptor, bonafide)

# If wanted, users can remove ligand-receptor interactions that were predicted based on 
# protein-protein interactions and 
# only keep ligand-receptor interactions that are described in curated databases.
# lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
head(lr_network)

# ## weighted integrated networks 
weighted_networks = readRDS(paste0(dataPath_nichenet,  "weighted_networks.rds"))
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

weighted_networks_lr = weighted_networks$lr_sig %>% 
  dplyr::rename(ligand = from, receptor = to) %>%
  inner_join(lr_network %>% 
               distinct(ligand,receptor), by = c("ligand","receptor"))

##########################################
# Step 1:  Define the niches/microenvironments of interest
##########################################
seurat_obj$celltype = as.factor(seurat_obj$celltype)
table(seurat_obj$celltype)
celltypes_all = levels(seurat_obj$celltype)
celltypes_BZ = celltypes_all[grep('_BZ', celltypes_all)]
celltypes_noBZ = setdiff(celltypes_all, celltypes_BZ)
rm(celltypes_all)

cat('cell types in BZ : ', celltypes_BZ, '\n')

cat('cell types in other regions : ', celltypes_noBZ, '\n')

for(receiver_cells in celltypes_BZ)
{
  
  # receiver_cells = 'CM_BZ'
  cat('--- receiver cells : ', receiver_cells, ' ---\n')
  
  # define the niches
  niches = list(
    "BZ_niche" = list(
      #"sender" = setdiff(celltypes_BZ, receiver_cells),
      "sender" = celltypes_BZ, # including autocrine
      "receiver" = receiver_cells),
    "noBZ_niche" = list(
      #"sender" = setdiff(celltypes_noBZ, gsub('_BZ',  '', receiver_cells)),
      "sender" = celltypes_noBZ, 
      "receiver" = gsub('_BZ',  '', receiver_cells))
  ) 
  
  ##########################################
  # step 2. Calculate differential expression between the niches 
  # issue solved because find markers for all genes, should run only for receptors and ligands
  ##########################################
  assay_oi = "RNA" 
  
  # only ligands important for sender cell types
  DE_sender = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), 
                                 niches = niches, 
                                 type = "sender", 
                                 assay_oi = assay_oi)
  
  # only receptors now, later on: DE analysis to find targets
  DE_receiver = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), 
                                   niches = niches, 
                                   type = "receiver", 
                                   assay_oi = assay_oi)
  
  DE_sender = DE_sender %>% 
    mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), 
                               ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
  DE_receiver = DE_receiver %>% 
    mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), 
                               ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
  
  # expected percentage of cells expressing the genes 
  expression_pct = 0.10 
  DE_sender_processed = process_niche_de(DE_table = DE_sender, 
                                         niches = niches, 
                                         expression_pct = expression_pct, 
                                         type = "sender")
  DE_receiver_processed = process_niche_de(DE_table = DE_receiver, 
                                           niches = niches, 
                                           expression_pct = expression_pct, 
                                           type = "receiver")
  
  # Combine sender-receiver DE based on L-R pairs:
  # As mentioned above, DE of ligands from one sender cell type is determined be calculating DE 
  # between that cell type, and all the sender cell types of the other niche. 
  # To summarize the DE of ligands of that cell type we have several options: 
  # we could take the average LFC, but also the minimum LFC compared to the other niche. 
  # We recommend using the minimum LFC, because this is the strongest specificity measure of ligand expression, 
  # because a high min LFC means that a ligand is more strongly expressed in the cell type of niche 1 
  # compared to all cell types of niche 2 (in contrast to a high average LFC, 
  # which does not exclude that one or more cell types in niche 2 also strongly express that ligand).
  specificity_score_LR_pairs = "min_lfc"
  DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, 
                                                  DE_receiver_processed, 
                                                  lr_network, 
                                                  specificity_score = specificity_score_LR_pairs)
  
  ##########################################
  # Step 3) Optional: Calculate differential expression between the different spatial regions 
  ##########################################
  # if not spatial info to include: put this to false 
  include_spatial_info_sender = FALSE 
  include_spatial_info_receiver = FALSE
  
  # user adaptation required on own dataset
  spatial_info = tibble(celltype_region_oi = "CAF_High", 
                        celltype_other_region = "myofibroblast_High", 
                        niche =  "pEMT_High_niche", celltype_type = "sender") 
  specificity_score_spatial = "lfc"
  
  # this is how this should be defined if you don't have spatial info (mock spatial info)
  if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
    spatial_info = tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% 
      mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
  } 
  
  
  if(include_spatial_info_sender == TRUE){
    sender_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% 
                                               subset(features = lr_network$ligand %>% unique()), 
                                             spatial_info = spatial_info %>% filter(celltype_type == "sender"))
    sender_spatial_DE_processed = process_spatial_de(DE_table = sender_spatial_DE, 
                                                     type = "sender", 
                                                     lr_network = lr_network, 
                                                     expression_pct = expression_pct, 
                                                     specificity_score = specificity_score_spatial)
    
    # add a neutral spatial score for sender celltypes in which the spatial is not known / not of importance
    sender_spatial_DE_others = get_non_spatial_de(niches = niches, 
                                                  spatial_info = spatial_info, 
                                                  type = "sender", 
                                                  lr_network = lr_network)
    sender_spatial_DE_processed = sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)
    
    sender_spatial_DE_processed = sender_spatial_DE_processed %>% 
      mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
    
  }else{
    #add a neutral spatial score for all sender celltypes (for none of them, spatial is relevant in this case)
    sender_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, 
                                                     type = "sender", lr_network = lr_network)
    sender_spatial_DE_processed = sender_spatial_DE_processed %>% 
      mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  
    
  }
  
  if(include_spatial_info_receiver == TRUE){
    receiver_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% 
                                                 subset(features = lr_network$receptor %>% unique()), 
                                               spatial_info = spatial_info %>% 
                                                 filter(celltype_type == "receiver"))
    receiver_spatial_DE_processed = process_spatial_de(DE_table = receiver_spatial_DE, 
                                                       type = "receiver", 
                                                       lr_network = lr_network, 
                                                       expression_pct = expression_pct, 
                                                       specificity_score = specificity_score_spatial)
    
    # add a neutral spatial score for receiver celltypes in which the spatial is not known / not of importance
    receiver_spatial_DE_others = get_non_spatial_de(niches = niches, 
                                                    spatial_info = spatial_info, 
                                                    type = "receiver", 
                                                    lr_network = lr_network)
    receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)
    
    receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% 
      mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
    
  }else{
    # # add a neutral spatial score for all receiver celltypes (for none of them, spatial is relevant in this case)
    receiver_spatial_DE_processed = get_non_spatial_de(niches = niches, 
                                                       spatial_info = spatial_info, 
                                                       type = "receiver", 
                                                       lr_network = lr_network)
    receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% 
      mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  }
  
  ##########################################
  # Step 4). Calculate ligand activities and infer active ligand-target links
  ##########################################
  ## first define the target genes by DE
  # It is always useful to check the number of genes in the geneset before doing the ligand activity analysis. 
  # We recommend having between 20 and 1000 genes in the geneset of interest, 
  # and a background of at least 5000 genes for a proper ligand activity analysis. 
  # If you retrieve too many DE genes, it is recommended to use a higher lfc_cutoff threshold. 
  # We recommend using a cutoff of 0.15 if you have > 2 receiver cells/niches to compare and 
  # use the min_lfc as specificity score. 
  # If you have only 2 receivers/niche, we recommend using a higher threshold (such as using 0.25). 
  # If you have single-cell data like Smart-seq2 with high sequencing depth, 
  #we recommend to also use higher threshold.
  # Smart-seq2 data and only 2 niches to compare, 
  # so we will use a stronger LFC threshold to keep less DE genes, but more trustworthy ones.
  
  lfc_cutoff = 0.15 # recommended for 10x as min_lfc cutoff. 
  specificity_score_targets = "min_lfc"
  
  # here DE all genes in the receiver cells, not only for ligand and receptors as before
  DE_receiver_targets = calculate_niche_de_targets(seurat_obj = seurat_obj, 
                                                   niches = niches, 
                                                   lfc_cutoff = lfc_cutoff, 
                                                   expression_pct = expression_pct, 
                                                   assay_oi = assay_oi)
  
  DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, 
                                                             niches = niches, 
                                                             expression_pct = expression_pct, 
                                                             specificity_score = specificity_score_targets)
  
  background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
  geneset_niche1 = DE_receiver_processed_targets %>% 
    filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & 
             target_present == 1) %>% pull(target) %>% unique()
  geneset_niche2 = DE_receiver_processed_targets %>% 
    filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & 
             target_present == 1) %>% pull(target) %>% unique()
  
  # Good idea to check which genes will be left out of the ligand activity analysis 
  # (=when not present in the rownames of the ligand-target matrix).
  # If many genes are left out, this might point to some issue in the gene naming 
  # (eg gene aliases and old gene symbols, bad human-mouse mapping)
  geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
  geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
  
  length(geneset_niche1)
  length(geneset_niche2)
  
  top_n_target = 250
  niche_geneset_list = list(
    "BZ_niche" = list(
      "receiver" = niches[[1]]$receiver,
      "geneset" = geneset_niche1,
      "background" = background),
    "noBZ_niche" = list(
      "receiver" = niches[[2]]$receiver,
      "geneset" = geneset_niche2 ,
      "background" = background)
  )
  
  # get_ligand_activites_targets calling the function predict_ligand_activities for each niche
  ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, 
                                                            ligand_target_matrix = ligand_target_matrix, 
                                                            top_n_target = top_n_target)
  
  head(ligand_activities_targets)
  ligand_activities_targets %>% arrange(-activity) %>% filter(receiver %in% niches$BZ_niche$receiver)
  
  
  ##########################################
  # step 5. Calculate (scaled) expression of ligands, receptors and targets 
  # across cell types of interest (log expression values and expression fractions)
  # we will calculate average (scaled) expression, and fraction of expression, of 
  # ligands, receptors, and target genes across all cell types of interest
  ##########################################
  features_oi = union(lr_network$ligand, lr_network$receptor) %>% 
    union(ligand_activities_targets$target) %>% setdiff(NA)
  
  # save dotplot 
  dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()), 
                                             features = features_oi, assay = assay_oi))
  
  exprs_tbl = dotplot$data %>% as_tibble()
  exprs_tbl = exprs_tbl %>% 
    dplyr::rename(celltype = id, gene = features.plot, expression = avg.exp, 
                  expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
    mutate(fraction = fraction/100) %>% as_tibble() %>% 
    select(celltype, gene, expression, expression_scaled, fraction) %>% 
    distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))
  
  exprs_tbl_ligand = exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% 
    dplyr::rename(sender = celltype, ligand = gene, ligand_expression = expression, 
                  ligand_expression_scaled = expression_scaled, ligand_fraction = fraction)
  
  exprs_tbl_receptor = exprs_tbl %>% 
    filter(gene %in% lr_network$receptor) %>% 
    dplyr::rename(receiver = celltype, receptor = gene, receptor_expression = expression, 
                  receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
  exprs_tbl_target = exprs_tbl %>% 
    filter(gene %in% ligand_activities_targets$target) %>% 
    dplyr::rename(receiver = celltype, target = gene, target_expression = expression, 
                  target_expression_scaled = expression_scaled, target_fraction = fraction)
  
  exprs_tbl_ligand = exprs_tbl_ligand %>% 
    mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% 
    mutate(ligand_fraction_adapted = ligand_fraction) %>% 
    mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% 
    mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))
  
  exprs_tbl_receptor = exprs_tbl_receptor %>% 
    mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled)) %>% 
    mutate(receptor_fraction_adapted = receptor_fraction) %>% 
    mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct) %>% 
    mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))
  
  
  ##########################################
  # step 6. Expression fraction and receptor
  # score ligand-receptor interactions based on expression strength of the receptor, 
  # in such a way that we give higher scores to the most strongly expressed receptor of a certain ligand, 
  # in a certain celltype. 
  # This will not effect the rank of individual ligands later on, 
  # but will help in prioritizing the most important receptors per ligand 
  # (next to other factors regarding the receptor - see later).
  ##########################################
  exprs_sender_receiver = lr_network %>% 
    inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
    inner_join(exprs_tbl_receptor, by = c("receptor")) %>% 
    inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))
  
  ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% 
    group_by(ligand, receiver) %>% 
    mutate(rank_receptor_expression = dense_rank(receptor_expression), 
           rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% 
    mutate(ligand_scaled_receptor_expression_fraction = 0.5* ((rank_receptor_fraction / max(rank_receptor_fraction)) + 
                                                                ((rank_receptor_expression / max(rank_receptor_expression))) ))  %>% 
    distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% 
    distinct() %>% ungroup() 
  
  
  ##########################################
  # step 7. Prioritization of ligand-receptor and ligand-target links 
  ##########################################
  prioritizing_weights = c(
    "scaled_ligand_score" = 5, # # niche-specific expression of the ligand: Recommended 5 (between 1-5)
    # scaled ligand expression in one sender compared to the other cell types in the dataset
    "scaled_ligand_expression_scaled" = 1,
    "ligand_fraction" = 1, # Ligands expressed in a smaller fraction of cells of cell type than cutoff (default: 0.10)
    "scaled_ligand_score_spatial" = 0, 
    # receptor DE score: niche-specific expression, Recommended 0.5 (>=0.5 and lower than "scaled_ligand_score")
    "scaled_receptor_score" = 1, 
    "scaled_receptor_expression_scaled" = 0.5, # Recommended weight: 0.5
    # Receptors that are expressed in a smaller fraction of cells of a cell type than exprs_cutoff(default: 0.10) 
    # will get a lower ranking, proportional to their fraction, Recommended weight: 1. 
    "receptor_fraction" = 1, 
    # this factor let us give higher weights to the most highly expressed receptor of a ligand in the receiver.
    # Recommended value: 1 (minimum: 0.5)
    "ligand_scaled_receptor_expression_fraction" = 1,
    "scaled_receptor_score_spatial" = 0,
    # Absolute ligand activity: to further prioritize ligand-receptor pairs based on their predicted effect of 
    # the ligand-receptor interaction on the gene expression in the receiver cell type - 
    # prioritizing_weights argument: "scaled_activity". Recommended weight: 0, 
    # unless absolute enrichment of target genes is of specific interest.
    "scaled_activity" = 0, 
    # Normalized ligand activity: to further prioritize ligand-receptor pairs based on their predicted effect of
    # the ligand-receptor interaction on the gene expression in the receiver cell type, Recommended weight: >=1.
    "scaled_activity_normalized" = 4,
    "bona_fide" = 1)
  
  output = list(DE_sender_receiver = DE_sender_receiver, 
                ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, 
                sender_spatial_DE_processed = sender_spatial_DE_processed, 
                receiver_spatial_DE_processed = receiver_spatial_DE_processed,
                ligand_activities_targets = ligand_activities_targets, 
                DE_receiver_processed_targets = DE_receiver_processed_targets, 
                exprs_tbl_ligand = exprs_tbl_ligand,  
                exprs_tbl_receptor = exprs_tbl_receptor, 
                exprs_tbl_target = exprs_tbl_target
  )
  
  prioritization_tables = get_prioritization_tables(output, prioritizing_weights)
  
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
    filter(receiver == niches[[1]]$receiver) %>% 
    head(10)
  
  prioritization_tables$prioritization_tbl_ligand_target %>% 
    filter(receiver == niches[[1]]$receiver) %>% 
    head(10)
  
  saveRDS(prioritization_tables, 
          file = paste0(outDir, '/nichenet_prioritization_tables_receiver_', receiver_cells, '.rds'))
  ##########################################
  # 8). Visualization of the Differential NicheNet output 
  ##########################################
  top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
    select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
    group_by(ligand) %>% 
    top_n(1, prioritization_score) %>% 
    ungroup() %>% 
    select(ligand, receptor, niche) %>% 
    dplyr::rename(top_niche = niche)
  
  top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
    select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
    group_by(ligand, receptor) %>% 
    top_n(1, prioritization_score) %>% 
    ungroup() %>% 
    select(ligand, receptor, niche) %>% 
    dplyr::rename(top_niche = niche)
  
  for(ntop in c(50, 100, 150, 200)){
    
    cat('ntop -- ', ntop, '\n')
    ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
      select(niche, sender, receiver, ligand, prioritization_score) %>% 
      group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% 
      ungroup() %>% distinct() %>% 
      inner_join(top_ligand_niche_df) %>% 
      filter(niche == top_niche) %>% 
      group_by(niche) %>% 
      top_n(n = ntop, prioritization_score) %>% 
      ungroup() # get the top50 ligands per niche
    
    receiver_oi = receiver_cells
    
    filtered_ligands = ligand_prioritized_tbl_oi %>% 
      filter(receiver == receiver_oi) %>% 
      pull(ligand) %>% unique()
    
    prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
      filter(ligand %in% filtered_ligands) %>% 
      select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
      distinct() %>% 
      inner_join(top_ligand_receptor_niche_df) %>% 
      group_by(ligand) %>% 
      filter(receiver == receiver_oi) %>% 
      top_n(2, prioritization_score) %>% ungroup() 
    
    lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, 
                                             prioritized_tbl_oi, 
                                             prioritization_tables$prioritization_tbl_ligand_receptor, 
                                             plot_legend = FALSE, 
                                             heights = NULL, widths = NULL)
    lfc_plot
    
    ggsave(paste0(outDir, '/Ligand_receptors_LFC_receiver.cells_', receiver_cells,  '_ntop', ntop, '.pdf'), 
           width = 12, height = 12*ntop/50, limitsize = FALSE)
    
    
    exprs_activity_target_plot = 
      make_ligand_activity_target_exprs_plot(receiver_oi, 
                                             prioritized_tbl_oi,  
                                             prioritization_tables$prioritization_tbl_ligand_receptor,  
                                             prioritization_tables$prioritization_tbl_ligand_target, 
                                             output$exprs_tbl_ligand,  
                                             output$exprs_tbl_target, 
                                             lfc_cutoff, 
                                             ligand_target_matrix, 
                                             plot_legend = FALSE, 
                                             heights = NULL, widths = NULL)
    
    exprs_activity_target_plot$combined_plot
    ggsave(paste0(outDir, '/Combined_plots_ligand_noFiltering_receiverCells_', receiver_cells,  '_ntop', ntop,  
                  '.pdf'), 
           width = 60, height = 12*ntop/50, limitsize = FALSE)
    
  }
  
}

Combine_NicheNet_Plots = FALSE
if(Combine_NicheNet_Plots){
   prioritization_tables = readRDS(file = paste0(outDir, 
                                                 '/nichenet_prioritization_tables_receiver_', 
                                                 receiver_cells, '.rds'))
  
   top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
     select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
     group_by(ligand) %>% 
     top_n(1, prioritization_score) %>% 
     ungroup() %>% 
     select(ligand, receptor, niche) %>% 
     dplyr::rename(top_niche = niche)
   
   top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
     select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
     group_by(ligand, receptor) %>% 
     top_n(1, prioritization_score) %>% 
     ungroup() %>% 
     select(ligand, receptor, niche) %>% 
     dplyr::rename(top_niche = niche)
   
   for(ntop in c(50, 100, 150, 200)){
     
     # ntop = 100
     cat('ntop -- ', ntop, '\n')
     ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
       select(niche, sender, receiver, ligand, prioritization_score) %>% 
       group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% 
       ungroup() %>% distinct() %>% 
       inner_join(top_ligand_niche_df) %>% 
       filter(niche == top_niche) %>% 
       group_by(niche) %>% 
       top_n(n = ntop, prioritization_score) %>% 
       ungroup() # get the top50 ligands per niche
     
     receiver_oi = receiver_cells
     
     filtered_ligands = ligand_prioritized_tbl_oi %>% 
       filter(receiver == receiver_oi) %>% 
       pull(ligand) %>% unique()
     
     prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
       filter(ligand %in% filtered_ligands) %>% 
       select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
       distinct() %>% 
       inner_join(top_ligand_receptor_niche_df) %>% 
       group_by(ligand) %>% 
       filter(receiver == receiver_oi) %>% 
       top_n(2, prioritization_score) %>% ungroup() 
     
     lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, 
                                              prioritized_tbl_oi, 
                                              prioritization_tables$prioritization_tbl_ligand_receptor, 
                                              plot_legend = FALSE, 
                                              heights = NULL, widths = NULL)
     lfc_plot
     
     ggsave(paste0(outDir, '/Ligand_receptors_LFC_receiver.cells_', receiver_cells,  '_ntop', ntop, '.pdf'), 
            width = 12, height = 12*ntop/50, limitsize = FALSE)
     
     
     exprs_activity_target_plot = 
       make_ligand_activity_target_exprs_plot(receiver_oi, 
                                              prioritized_tbl_oi,  
                                              prioritization_tables$prioritization_tbl_ligand_receptor,  
                                              prioritization_tables$prioritization_tbl_ligand_target, 
                                              output$exprs_tbl_ligand,  
                                              output$exprs_tbl_target, 
                                              lfc_cutoff, 
                                              ligand_target_matrix, 
                                              plot_legend = FALSE, 
                                              heights = NULL, widths = NULL)
     
     exprs_activity_target_plot$combined_plot
     ggsave(paste0(outDir, '/Combined_plots_ligand_noFiltering_receiverCells_', receiver_cells,  '_ntop', ntop,  
                   '.pdf'), 
            width = 80, height = 8*ntop/50, limitsize = FALSE)
     
   }
   
   
  # ## select only top 20 ligands
  filtered_ligands = ligand_prioritized_tbl_oi %>%
    filter(receiver == receiver_oi) %>%
    top_n(50, prioritization_score) %>%
    pull(ligand) %>% unique()

  prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>%
    filter(ligand %in% filtered_ligands) %>%
    select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>%
    distinct() %>% inner_join(top_ligand_receptor_niche_df) %>%
    group_by(ligand) %>% filter(receiver == receiver_oi) %>%
    top_n(2, prioritization_score) %>% ungroup()

  exprs_activity_target_plot =
    make_ligand_activity_target_exprs_plot(receiver_oi,
                                           prioritized_tbl_oi,
                                           prioritization_tables$prioritization_tbl_ligand_receptor,
                                           prioritization_tables$prioritization_tbl_ligand_target,
                                           output$exprs_tbl_ligand,
                                           output$exprs_tbl_target,
                                           lfc_cutoff,
                                           ligand_target_matrix,
                                           plot_legend = FALSE,
                                           heights = NULL, widths = NULL)

  exprs_activity_target_plot$combined_plot

  ggsave(paste0(outDir, '/Combined_plots_ligand_receiverCells_', receiver_cells, '_top50.pdf'),
         width = 50, height = 14, limitsize = FALSE)
}


## Circos plot of prioritized ligand-receptor pairs
Make.Circos.plot = FALSE
if(Make.Circos.plot){
  filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% 
    top_n(15, prioritization_score) %>% pull(ligand) %>% unique()
  
  prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
    filter(ligand %in% filtered_ligands) %>% 
    select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
    distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% 
    group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 
  
  colors_sender = brewer.pal(n = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% 
                               length(), name = 'Spectral') %>% 
    magrittr::set_names(prioritized_tbl_oi$sender %>% unique() %>% sort())
  colors_receiver = c("lavender")  %>% magrittr::set_names(prioritized_tbl_oi$receiver %>% unique() %>% sort())
  
  circos_output = make_circos_lr(prioritized_tbl_oi, colors_sender, colors_receiver)
  
}


