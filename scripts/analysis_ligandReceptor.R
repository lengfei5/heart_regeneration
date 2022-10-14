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
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered_manualSegmentation', 
                   species, '.Rdata')) # visium data

## snRNA-seq 
refs_file = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets_20221004_2.rds'
refs = readRDS(file = refs_file)
table(refs$subtypes)

##########################################
# run LIANA 
##########################################
# set parameter for ligand-receptor analysis
outDir = paste0(resDir, '/Ligand_Receptor_analysis/LIANA')

refs$celltypes = refs$subtypes
celltypes = c('Mono_Macrophages', 'Proliferating_CM', 'FB_1', 'Injury_specific_EC')

source('functions_Visium.R')
run_LIANA(refs, celltypes = celltypes, outDir = outDir)


########################################################
########################################################
# Section : Nichenet for ligand-receptor analysis
# original code from https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md
########################################################
########################################################
outDir = paste0(resDir, '/Ligand_Receptor_analysis/NicheNet')
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
refs$celltypes[grep('CM|EC|FB|MP', refs$celltypes, invert = TRUE)] = NA

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
subref = subset(x = subref, downsample = 3000) # downsample the CM and EC for the sake of speed
table(subref$celltypes)

##########################################
# Gene filtering, in particular filtering gene without UMI counts or lowly expressed genes
##########################################
Gene.filtering.preprocessing = TRUE
if(Gene.filtering.preprocessing){
  sce <- as.SingleCellExperiment(subref)
  colLabels(sce) = as.factor(sce$celltypes)
  rownames(sce) = get_geneName(rownames(sce))
  
  ave.counts <- calculateAverage(sce, assay.type = "counts")
  
  hist(log10(ave.counts), breaks=100, main="", col="grey80",
       xlab=expression(Log[10]~"average count"))
  
  num.cells <- nexprs(sce, byrow=TRUE)
  smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
                xlab=expression(Log[10]~"average count"))
  
  # detected in >= 5 cells, ave.counts >=5 but not too high
  genes.to.keep <- num.cells > 50 & ave.counts >= 10^-2  & ave.counts <10^1  
  summary(genes.to.keep)
  sce <- sce[genes.to.keep, ]
  rownames(sce) = make.names(rownames(sce), unique = TRUE)
  subref = as.Seurat(sce, counts = "counts", data = "logcounts")
  rm(sce)
  Idents(subref) = as.factor(subref$celltypes)
  
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

subref <- RunUMAP(subref, dims = 1:30, n.neighbors = 30, min.dist = 0.1, reduction = 'pca', reduction.key = 'umap_')
DimPlot(subref, group.by = "celltype", label = TRUE, reduction = 'umap') 

saveRDS(subref, file = paste0(RdataDir, 'seuratObject_snRNAseq_subset_for_NicheNet_', species, '.rds')) 

# reload the processed seurat object
subref = readRDS(file = paste0(RdataDir, 'seuratObject_snRNAseq_subset_for_NicheNet_', species, '.rds')) 
#subref = subset(x = subref, downsample = 1000) # downsample the CM and EC for the sake of speed
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
celltype_noBZ = setdiff(celltypes_all, celltypes_BZ)
rm(celltypes_all)

receiver_cells = 'CM_BZ'

# define the niches
niches = list(
  "BZ_niche" = list(
    "sender" = setdiff(celltypes_BZ, receiver_cells),
    "receiver" = receiver_cells),
  "noBZ_niche" = list(
    "sender" = setdiff(celltype_noBZ, gsub('_BZ',  '', receiver_cells)),
    "receiver" = gsub('_BZ',  '', receiver_cells))
) 

# ## receiver
# receiver = "Proliferating_CM"
# expressed_genes_receiver = get_expressed_genes(receiver, subref, pct = 0.10)
# background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
# 
# ## sender
# sender_celltypes = c('Mono_Macrophages', 'Proliferating_CM', 'Neutrophil', 'Injury_specific_EC')
# # lapply to get the expressed genes of every sender cell type separately here
# list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, subref, 0.10) 
# expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

##########################################
# step 2. Calculate differential expression between the niches 
# for unknown reasons, this DE step takes really long time to finish
# issue sovled because find markers for all genes, should run only for receptors and ligands
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

#save(DE_sender, DE_receiver, 
#     file = paste0(RdataDir, 'seuratObject_snRNAseq_subset_for_NicheNet_DE_sender_receiver_', 
#                                 species, '.Rdata')) 

DE_sender = DE_sender %>% 
  mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), 
                             ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
DE_receiver = DE_receiver %>% 
  mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), 
                             ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))

expression_pct = 0.10 # expected percentage of cells expressing the genes 
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
include_spatial_info_sender = FALSE 
# if not spatial info to include: put this to false # user adaptation required on own dataset
include_spatial_info_receiver = FALSE 
# if spatial info to include: put this to true # user adaptation required on own dataset

# user adaptation required on own dataset
spatial_info = tibble(celltype_region_oi = "CAF_High", 
                      celltype_other_region = "myofibroblast_High", 
                      niche =  "pEMT_High_niche", celltype_type = "sender") 
specificity_score_spatial = "lfc"

# this is how this should be defined if you don't have spatial info
# mock spatial info

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
  # # add a neutral spatial score for all sender celltypes (for none of them, spatial is relevant in this case)
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
  
} else {
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
lfc_cutoff = 0.15 # recommended for 10x as min_lfc cutoff. 
specificity_score_targets = "min_lfc"

# here DE all genes, not only for ligand and receptors as before
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

top_n_target = 500

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

ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, 
                                                          ligand_target_matrix = ligand_target_matrix, 
                                                          top_n_target = top_n_target)


##########################################
# step 5. Calculate (scaled) expression of ligands, receptors and targets 
# across cell types of interest (log expression values and expression fractions) 
##########################################
features_oi = union(lr_network$ligand, lr_network$receptor) %>% 
  union(ligand_activities_targets$target) %>% setdiff(NA)

# save dotplot 
dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()), 
                                           features = features_oi, assay = assay_oi))

dotplot_data  = dotplot$data
colnames(dotplot_data) = c('expression', 'fraction', 'gene', 'celltype', 'expression_scaled')

exprs_tbl = dotplot_data %>% as_tibble()
exprs_tbl = exprs_tbl %>% 
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
prioritizing_weights = c("scaled_ligand_score" = 2,
                         "scaled_ligand_expression_scaled" = 1,
                         "ligand_fraction" = 1,
                         "scaled_ligand_score_spatial" = 2, 
                         "scaled_receptor_score" = 0.5,
                         "scaled_receptor_expression_scaled" = 0.5,
                         "receptor_fraction" = 1, 
                         "ligand_scaled_receptor_expression_fraction" = 1,
                         "scaled_receptor_score_spatial" = 0,
                         "scaled_activity" = 0,
                         "scaled_activity_normalized" = 2,
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

ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  select(niche, sender, receiver, ligand, prioritization_score) %>% 
  group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% 
  ungroup() %>% distinct() %>% 
  inner_join(top_ligand_niche_df) %>% 
  filter(niche == top_niche) %>% 
  group_by(niche) %>% 
  top_n(50, prioritization_score) %>% 
  ungroup() # get the top50 ligands per niche

receiver_oi = "CM_BZ" 

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

ggsave(paste0(outDir, '/Ligand_receptors_LFC.pdf'), width = 10, height = 12)


exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(receiver_oi, 
                                                                    prioritized_tbl_oi,  
                                                prioritization_tables$prioritization_tbl_ligand_receptor,  
                                                prioritization_tables$prioritization_tbl_ligand_target, 
                                                output$exprs_tbl_ligand,  
                                                output$exprs_tbl_target, 
                                                lfc_cutoff, ligand_target_matrix, 
                                                plot_legend = FALSE, 
                                                heights = NULL, widths = NULL)
exprs_activity_target_plot$combined_plot

filtered_ligands = ligand_prioritized_tbl_oi %>% 
  filter(receiver == receiver_oi) %>% 
  top_n(20, prioritization_score) %>% 
  pull(ligand) %>% unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(ligand %in% filtered_ligands) %>% 
  select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
  distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% 
  group_by(ligand) %>% filter(receiver == receiver_oi) %>% 
  top_n(2, prioritization_score) %>% ungroup() 

exprs_activity_target_plot = make_ligand_activity_target_exprs_plot(receiver_oi, 
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

ggsave(paste0(outDir, '/Combined_plots_top20_ligand.pdf'), width = 45, height = 12)

## Circos plot of prioritized ligand-receptor pairs
filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(15, prioritization_score) %>% pull(ligand) %>% unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 

colors_sender = brewer.pal(n = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% length(), name = 'Spectral') %>% magrittr::set_names(prioritized_tbl_oi$sender %>% unique() %>% sort())
colors_receiver = c("lavender")  %>% magrittr::set_names(prioritized_tbl_oi$receiver %>% unique() %>% sort())

circos_output = make_circos_lr(prioritized_tbl_oi, colors_sender, colors_receiver)


