##########################################################################
##########################################################################
# Project: main script for cross species analysis 
# Script purpose: 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Feb 20 15:26:14 2024
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
library(liana)
source('functions_scRNAseq.R')
source('functions_Visium.R')
source('functions_cccInference.R')

version.analysis = '_20240220'

resDir = paste0("../results/cross_species", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

mem_used()

########################################################
########################################################
# Section I: prepare the scRNA-seq for each species
# 
########################################################
########################################################
##########################################
## axolotl snRNA-seq data with single batch 
##########################################
refs_file = paste0('/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_subtypes_final_20221117.rds')

refs = readRDS(file = refs_file)
table(refs$subtypes)
length(table(refs$subtypes))

refs$subtypes = droplevels(refs$subtypes) 
length(table(refs$subtypes)) 

refs$condition = gsub('_scRNA', '', refs$condition)

## prepare the celltype to use and also specify the time-specific subtypes
refs$celltype_toUse = as.character(refs$subtypes)
length(table(refs$celltype_toUse))

refs$celltype_toUse = gsub('Mo/Macs', 'Mo.Macs', refs$celltype_toUse)
refs$celltype_toUse = gsub("[(]", '', refs$celltype_toUse)
refs$celltype_toUse = gsub("[)]", '', refs$celltype_toUse)
refs$celltype_toUse = gsub('CM_ven_Robo2', 'CM_Robo2', refs$celltype_toUse)
refs$celltype_toUse = gsub('CM_ven_Cav3_1', 'CM_Cav3.1', refs$celltype_toUse)

table(refs$celltype_toUse)
length(table(refs$celltype_toUse))

refs$subtypes = refs$celltype_toUse # clean the special symbols
refs$subtypes = as.factor(refs$subtypes) 
table(refs$subtypes)
length(table(refs$subtypes))

refs$celltypes = as.character(refs$celltype_toUse)

refs$celltypes[grep('CM_|CMs_|_CM|_CM_', refs$subtypes)] = 'CM'
refs$celltypes[grep('EC_|_EC', refs$subtypes)] = 'EC'
refs$celltypes[grep('FB_', refs$subtypes)] = 'FB'
refs$celltypes[grep('B_cells', refs$subtypes)] = 'Bcell'
refs$celltypes[grep('T_cells', refs$subtypes)] = 'Tcell'

refs$celltypes[grep('Macrophages|_MF', refs$subtypes)] = 'Macrophages'
refs$celltypes[grep('Megakeryocytes', refs$subtypes)] = 'Megakeryocytes'
refs$celltypes[grep('RBC', refs$subtypes)] = 'RBC'
refs$celltypes[grep('Neu_', refs$subtypes)] = 'Neutrophil'
refs$celltypes[grep('Mo.Macs', refs$subtypes)] = 'Mo.Macs'

table(refs$celltypes)
length(table(refs$celltypes))

ax = refs 
rm(refs)

saveRDS(ax, file = paste0(RdataDir, 'ax_scRNAseq_saved_1batch.rds'))

rm(ax)

##########################################
## neonatal scRNA-seq data after batch correction
##########################################
nm = readRDS(file = paste0('../data/data_examples/ref_scRNAseq_neonatalMice_clean.v1.2.rds'))

head(rownames(nm))

#DefaultAssay(cms) = 'integrated'
p1 = DimPlot(nm, group.by = 'celltype', label = TRUE, repel = TRUE)
p2 = DimPlot(nm, group.by = 'condition',label = TRUE, repel = TRUE)

p1 + p2

ggsave(paste0(resDir, '/scRNAseq_neonatalMice.pdf'), 
       width = 16, height = 6)


DimPlot(nm, group.by = 'dataset', label = TRUE, repel = TRUE)

ggsave(paste0(resDir, '/scRNAseq_dataSource_neonatalMice.pdf'), 
       width = 8, height = 6)

saveRDS(nm, file = paste0(RdataDir, 'nm_scRNAseq_8batchCorrected.rds'))

rm(nm)

##########################################
# import the neonatal mouse data before batch correction
##########################################
aa = readRDS(file = paste0('../results/scRNAseq_neonatalMouse_20231108/Rdata/', 
                           'Seurat.obj_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_regress.nUMI.rds'))

aa$species = 'nm'
aa$dataset = factor(aa$dataset, levels = c('Wang2020', 'Cui2020'))
aa$batch = paste0(aa$dataset, '_', aa$condition)
aa$batch = factor(aa$batch)

aa$time[which(aa$time == 'D1')] = 'd1'
aa$time[which(aa$time == 'D3')] = 'd3'

aa$timepoints = aa$time

aa$celltype = aa$BroadID
aa$subtype = aa$FineID

aa = subset(aa, cells = colnames(aa)[which(aa$celltype != 'glial')])

p1 = DimPlot(aa, group.by = 'subtype', label = TRUE, repel = TRUE, raster=FALSE)
p2 = DimPlot(aa, group.by = 'batch', label = TRUE, repel = TRUE) 
p3 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE, raster=FALSE)
p4 = DimPlot(aa, group.by = 'dataset', label = TRUE, repel = TRUE, raster=FALSE)
(p1 + p4) /(p3 + p2)

ggsave(paste0(resDir, 'neonatalMice_CM.Cui2020_noCM.Wang2020_P1_overview.pdf'), 
       width = 16, height = 12)

aa = DietSeurat(aa, counts = TRUE, data = TRUE, scale.data = FALSE, assays = 'RNA')
aa = NormalizeData(aa) %>% FindVariableFeatures(nfeatures = 3000) %>%  
  ScaleData() %>% 
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

p1 = DimPlot(aa, group.by = 'subtype', label = TRUE, repel = TRUE, raster=FALSE)
p2 = DimPlot(aa, group.by = 'batch', label = TRUE, repel = TRUE) 
p3 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE, raster=FALSE)
p4 = DimPlot(aa, group.by = 'dataset', label = TRUE, repel = TRUE, raster=FALSE)
(p1 + p4) /(p3 + p2)

ggsave(paste0(resDir, 'neonatalMice_CM.Cui2020_noCM.Wang2020_P1_overview_noBatchCorrection.pdf'), 
       width = 16, height = 12)

saveRDS(aa, file = paste0(resDir, 'nm_scRNAseq_8batches_noCorrection.rds'))

rm(aa)

##########################################
## adult mice after batch correction 
##########################################
mm = readRDS(file = paste0('../data/data_examples/ref_scRNAseq_adultMice_clean.v1.rds'))
head(rownames(mm))

mm$condition = factor(mm$condition, levels = c('d0', 'd1', 'd3', 'd5', 'd7', 'd14', 'd28'))

p1 = DimPlot(mm, reduction = 'umap', group.by = 'celltype', raster = T,shuffle= T, pt.size = 2, 
             label = TRUE, repel = TRUE)

p2 = DimPlot(mm, reduction = 'umap', group.by = 'condition',raster = T,shuffle= T, pt.size = 2, 
       label = TRUE, repel = TRUE)
#p2 = FeaturePlot(mm, features = 'Axl') +  
#  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))

p1 + p2

#outDir = paste0(resDir, '/noCM.Forte2020_CM.Ren2020_integrated_LRanalysis/')
#if(!dir.exists(outDir)) dir.create(outDir)

ggsave(paste0(resDir, '/scRNAseq_adultMice_Forte2020_Ren2020_Integration.pdf'), 
       width = 16, height = 6)

DimPlot(mm, group.by = 'dataset', label = TRUE, repel = TRUE)

ggsave(paste0(resDir, '/scRNAseq_dataSource_adultMice.pdf'), 
       width = 8, height = 6)

mm$batch = mm$dataset

#saveRDS(ax, file = paste0(RdataDir, 'ax_scRNAseq.rds'))
saveRDS(mm, file = paste0(RdataDir, 'mm_scRNAseq_2batches_withCorrection.rds'))

rm(mm)

# mm = readRDS(file = paste0(RdataDir, 'mm_scRNAseq_2batches_withCorrection.rds'))

##########################################
# adult mice without batch correction 
##########################################
bb = readRDS(file = paste0('../results/Rdata/',
                           'Forte2020_logNormalize_allgenes_majorCellTypes_subtypes.rds'))
cms = readRDS(file =  paste0('../results/Rdata/', 
                             'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization_umap_subtypes.rds'))

bb$dataset = 'Forte2020'
cms$dataset = 'Ren2020'

features.common = intersect(rownames(bb), rownames(cms))
length(union(rownames(bb), rownames(cms)))

mm = merge(bb, y = cms, add.cell.ids = c("Forte2020", "Ren2020"), project = "adultHeart")

rm(list=c('bb', 'cms'))

mm$species = 'mm'
jj = which(is.na(mm$timepoints))
mm$timepoints[grep('0w_', mm$sample)] = 'd0'
mm$timepoints[grep('2w_', mm$sample)] = 'd14'
mm$condition = mm$timepoints
mm$time = mm$timepoints

mm = NormalizeData(mm)
mm = FindVariableFeatures(mm, selection.method = "vst", nfeatures = 3000)
mm <- ScaleData(mm, verbose = FALSE)
mm <- RunPCA(mm, npcs = 30, verbose = FALSE)

ElbowPlot(mm, ndims = 30)

mm <- FindNeighbors(mm, reduction = "pca", dims = 1:20)
#mm <- FindClusters(mm, resolution = 0.2)
mm <- RunUMAP(mm, reduction = "pca", dims = 1:30, n.neighbors = 50, min.dist = 0.05) 

p1 = DimPlot(mm, group.by = 'subtype', label = TRUE, repel = TRUE, raster=FALSE)
p2 = DimPlot(mm, group.by = 'condition', label = TRUE, repel = TRUE) 
p3 = DimPlot(mm, group.by = 'celltype', label = TRUE, repel = TRUE, raster=FALSE)
p4 = DimPlot(mm, group.by = 'dataset', label = TRUE, repel = TRUE, raster=FALSE)
(p1 + p4) /(p3 + p2)

ggsave(paste0(resDir, 'adultMice_overview_noBatchCorrection.pdf'), 
       width = 16, height = 12)

mm$batch = mm$dataset

saveRDS(mm, file = paste0(resDir, 'mm_scRNAseq_2batch_noCorrection.rds'))

rm(mm)




########################################################
########################################################
# Section II: cross-species cell types/states alignment 
# SAMAP 
########################################################
########################################################


## define the one-on-one ortholog between axololt and mice
an_orthologs = data.frame(ref = rownames(ax), query = rownames(ax))
rownames(an_orthologs) = an_orthologs$ref
an_orthologs$query = sapply(an_orthologs$query, 
                            function(x){firstup(unlist(strsplit(as.character(x), '-'))[1])})

jj = which(!is.na(match(an_orthologs$query, rownames(aa))))
an_orthologs = an_orthologs[jj, ]

counts = table(an_orthologs$query)
gg_uniq = names(counts)[which(counts == 1)]
jj2 = which(!is.na(match(an_orthologs$query, gg_uniq)))

an_orthologs = an_orthologs[jj2, ]

ax = subset(ax, features = an_orthologs$ref)

aa = subset(aa, features = an_orthologs$query)


counts = ax@assays$RNA@counts
metadata = ax@meta.data
counts = counts[match(an_orthologs$ref, rownames(counts)), ]
rownames(counts) = an_orthologs$query

new_ax <- CreateSeuratObject(counts=counts, assay = 'RNA', meta.data = metadata)

rm(list = c('counts', 'metadata'))

aa = merge(aa, y = new_ax, add.cell.ids = c("m", "ax"), project = "heartReg")

rm(new_ax)
#aa = DietSeurat(aa, counts = TRUE, data = TRUE, scale.data = FALSE, assays = 'RNA')
#jj = which(is.na(aa$batch))
#aa$batch[jj] = aa$dataset[]






########################################################
########################################################
# Section III: cross-species cell types/states neighborhood comparison and differential analysis
# 
########################################################
########################################################


mm = readRDS(file = paste0(RdataDir, 'mm_scRNAseq.rds'))


########################################################
########################################################
# Section III : Differential LR analysis in ST 
# 
########################################################
########################################################
outDir = paste0(resDir, '/LR_differentialAnalysis/')
if(!dir.exists(outDir)) dir.create(outDir)

##########################################
# test the LIANA + tensor_cell2cell
# original code : https://saezlab.github.io/liana/articles/liana_cc2tensor.html
##########################################
library(tidyverse, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)
library(reticulate, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(liana, quietly = TRUE)
library(ExperimentHub, quietly = TRUE)

eh <- ExperimentHub()
# Get Data
(sce <- eh[["EH2259"]])


# basic feature filtering
sce <- sce[rowSums(counts(sce) >= 1) >= 5, ]

# basic outlier filtering
qc <- scater::perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- scater::isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]

# Remove doublets
sce <- sce[, sce$multiplets=="singlet"]

# Set rownames to symbols
rownames(sce) <- rowData(sce)$SYMBOL

# log-transform
sce <- scuttle::logNormCounts(sce)

# Create a label unique for every sample
sce$context <- paste(sce$stim, sce$ind, sep="|")

# Plot
sce %>%
  get_abundance_summary(sample_col = "context",
                        idents_col = "cell", 
                        min_cells = 10, # min cells per sample
                        min_samples = 3, # min samples
                        min_prop = 0.2 # min prop of samples
  ) %>%
  plot_abundance_summary()

# filter non abundant celltypes
sce <- liana::filter_nonabundant_celltypes(sce,
                                           sample_col = "context",
                                           idents_col = "cell")

# Run LIANA by sample
sce <- liana_bysample(sce = sce,
                      sample_col = "context",
                      idents_col = "cell",
                      method = "sca", # we use SingleCellSignalR's score alone
                      expr_prop = 0, # expression proportion threshold
                      inplace=TRUE, # saves inplace to sce
                      return_all = FALSE # whether to return non-expressed interactions 
)
