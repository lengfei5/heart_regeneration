##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: prepare the scRNA-seq data for each species 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Nov 11 10:55:32 2024
##########################################################################
##########################################################################
rm(list = ls())
# required libraries
library(data.table)
require(Seurat)
library(SeuratObject)
library(SeuratDisk)
require(sctransform)
require(ggplot2)
library(dplyr)
library(patchwork)
require(tictoc)
library(DropletUtils)
library(edgeR)
library(future)
library(pryr) # monitor the memory usage

source('functions_scRNAseq.R')
source('functions_Visium.R')

options(future.globals.maxSize = 160000 * 1024^2)
mem_used()

version.analysis = '_20240220'
resDir = paste0("../results/cross_species", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

outDir = paste0(resDir, '/samap_out/')
system(paste0('mkdir -p ', outDir))


########################################################
########################################################
# Section I : mouse data
# 
########################################################
########################################################

##########################################
# import the neonatal mouse data before batch correction
# because the samap requires count tables
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

ggsave(paste0(outDir, 'neonatalMice_CM.Cui2020_noCM.Wang2020_P1_overview.pdf'), 
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

ggsave(paste0(outDir, 'neonatalMice_CM.Cui2020_noCM.Wang2020_P1_overview_noBatchCorrection.pdf'), 
       width = 16, height = 12)

saveRDS(aa, file = paste0(outDir, 'nm_scRNAseq_no.batchCorrection.rds'))

##########################################
# import adult mice data and merge with neonatal mouse 
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

ggsave(paste0(outDir, 'adultMice_overview_noBatchCorrection.pdf'), 
       width = 16, height = 12)

length(intersect(rownames(mm), rownames(aa)))
length(union(rownames(mm), rownames(aa)))

saveRDS(aa, file = paste0(outDir, 'mm_scRNAseq_no.batchCorrection.rds'))


aa = merge(aa, y = mm, add.cell.ids = c("nm", "mm"), project = "heartReg")

aa$model = aa$species
aa$species = 'mouse'

jj = which(is.na(aa$batch))
aa$batch[jj] = aa$dataset[jj]

aa = NormalizeData(aa) %>% FindVariableFeatures(nfeatures = 3000) %>%  
  ScaleData() %>% 
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

p1 = DimPlot(aa, group.by = 'subtype', label = TRUE, repel = TRUE, raster=FALSE)
p2 = DimPlot(aa, group.by = 'batch', label = TRUE, repel = TRUE) 
p3 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE, raster=FALSE)
p4 = DimPlot(aa, group.by = 'dataset', label = TRUE, repel = TRUE, raster=FALSE)
(p1 + p4) /(p3 + p2)

ggsave(paste0(outDir, 'neonatalMice_adultMice_noBatchCorrection.pdf'), 
       width = 20, height = 16)


##########################################
# transcript-to-gene mapping
##########################################
aa = readRDS(file = paste0(outDir, 'nm_mm_scRNAseq_merged_v1.rds'))

ggs = rtracklayer::import(paste0('/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/',
                                 'Mus_musculus.GRCm38.87.gtf'))

ggs = ggs[which(ggs$type == 'transcript')]
ggs = as.data.frame(ggs)
ggs = ggs[, c(17, 10, 12)]

cds = read.delim(file = paste0('/groups/tanaka/People/current/jiwang/src/samap_directory/example_data/',
                               'transcriptomes/mm_transcript_gene.txt'), 
                 sep = ' ', header = FALSE)
cds = cds[, c(1, 4, 7)]
cds$V1 = gsub('[>]', '', cds$V1)
cds$V4 = gsub('gene:', '', cds$V4)
cds$V7 = gsub('gene_symbol:', '', cds$V7)

mm = match(rownames(aa), cds$V7)
length(which(is.na(mm)))

missed = (rownames(aa)[which(is.na(mm))])

colnames(cds) = c('transcriptID', 'geneID', 'geneSymbol')
cds = cds[, c(1, 3, 2)]

write.table(cds, file = paste0(outDir, 'mice_transcript_gene.txt'), sep = '\t', row.names = FALSE, 
            col.names = TRUE, quote = FALSE)

mm = match(rownames(aa), cds$geneSymbol)
length(which(is.na(mm)))
features = rownames(aa)[which(!is.na(mm))]

aa = subset(aa, features = features)

saveRDS(aa, file = paste0(outDir, 'nm_mm_scRNAseq_merged.rds'))

saveFile = paste0(outDir, 'nm_mm_scRNAseq_merged.h5Seurat')
SaveH5Seurat(aa, filename = saveFile, overwrite = TRUE)
Convert(saveFile, dest = "h5ad", overwrite = TRUE)


########################################################
########################################################
# Section II: axolotl data  
# 
########################################################
########################################################
ax = readRDS(file = paste0(RdataDir, 'ax_scRNAseq.rds'))

ax$species = 'ax'
ax$model = 'ax'
ax$dataset = 'ax_Bassat2024'
ax$batch = ax$dataset

# aa = readRDS(file = paste0(RdataDir, 'nm_mm_scRNAseq_merged_v1.rds'))
# 
# aa = DietSeurat(aa, counts = TRUE, data = TRUE, scale.data = FALSE, assays = 'RNA')
# table(aa$batch)
# 
# aa$dataset[which(aa$dataset == 'Cui2020')] = 'nm_Cui2020'
# aa$dataset[which(aa$dataset == 'Wang2020')] = 'nm_Wang2020'
# aa$dataset[which(aa$dataset == 'Forte2020')] = 'mm_Forte2020'
# aa$dataset[which(aa$dataset == 'Ren2020')] = 'mm_Ren2020'
# 
# aa$batch[which(aa$dataset == 'mm_Forte2020')] = 'mm_Forte2020'
# aa$batch[which(aa$dataset == 'mm_Ren2020')] = 'mm_Ren2020'

ax$subtype = ax$celltypes

ax$celltype = ax$celltypes

ax$celltype[grep('FB_', ax$celltype)] = 'FB'
ax$celltype[grep('B_', ax$celltype)] = 'B'
ax$celltype[grep('CM_', ax$celltype)] = 'CM'
ax$celltype[grep('EC_', ax$celltype)] = 'EC'
ax$celltype[grep('Mo.Macs_', ax$celltype)] = 'Mo.Macs'

ax = subset(ax, cells = colnames(ax)[which(ax$celltype != 'Neuronal' & ax$celltype != 'RBC' &
                                             ax$celltype != 'Proliferating_RBC')])

ax = subset(ax, cells = colnames(ax)[which(ax$celltype != 'Neu_IL1R1')])



##########################################
# transcript-to-gene mapping for axolotl
##########################################
ax = readRDS(file = paste0(outDir, 'ax_scRNAseq.rds'))

ggs = read.delim(paste0('/groups/tanaka/People/current/jiwang/Genomes/axolotl/Transcriptomics/',
                                 'AmexT_v47_transcripts_genes_t2g.txt'), sep = '\t', header = FALSE)

cds = read.delim(file = paste0('/groups/tanaka/People/current/jiwang/src/samap_directory/example_data/',
                               'transcriptomes/ax_transcript.txt'), 
                 sep = '\t', header = FALSE)
#cds = cds[, c(1, 4, 7)]
cds$V1 = gsub('[>]', '', cds$V1)
#cds$V4 = gsub('gene:', '', cds$V4)
#cds$V7 = gsub('gene_symbol:', '', cds$V7)

cds$geneID = ggs$V2[match(cds$V1, ggs$V1)]
rm(ggs)

ggs = data.frame(gene = rownames(ax))
ggs$geneID = sapply(ggs$gene, get_geneID)

cds$gene = ggs$gene[match(cds$geneID, ggs$geneID)] 

cds = cds[which(!is.na(cds$gene)), ]

mm = match(rownames(ax), cds$gene)
length(which(is.na(mm)))

missed = (rownames(aa)[which(is.na(mm))])

colnames(cds) = c('transcriptID', 'geneID', 'geneSymbol')
cds = cds[, c(1, 3, 2)]

write.table(cds, file = paste0(outDir, 'ax_transcript_gene.txt'), sep = '\t', row.names = FALSE, 
            col.names = TRUE, quote = FALSE)

mm = match(rownames(ax), cds$geneSymbol)
length(which(is.na(mm)))
#features = rownames(aa)[which(!is.na(mm))]
#aa = subset(aa, features = features)

#saveRDS(ax, file = paste0(outDir, 'ax_scRNAseq.rds'))

library(SeuratDisk)
VariableFeatures(ax)
#Idents(aa) = droplevels(aa$condition)

saveFile = paste0(outDir, 'ax_scRNAseq.h5Seurat')
SaveH5Seurat(aa, filename = saveFile, overwrite = TRUE)
Convert(saveFile, dest = "h5ad", overwrite = TRUE)



########################################################
########################################################
# Section III:
# 
########################################################
########################################################
##########################################
# import the zebrafish data 
##########################################
dr = readRDS(file = paste0('../published_dataset/zebrafish/Hu_Junker_2022/zebrafish_heart_processing/Rdata',
                           '/all_heart_annotation_normalized.pca.umap_filtered.duplex.dead.rds'))

##########################################
# correct cell labels for scANVI
##########################################
aa = readRDS(paste0(RdataDir, 'nm_mm_ax_scRNAseq_merged_v2.rds'))

aa$celltype[which(aa$celltype == 'T_cells')] = 'T'
aa$celltype[which(aa$celltype == 'prolife.Mphage')] = 'MP'
aa$celltype[which(aa$celltype == 'MHCII.Mphage')] = 'MP'

aa$celltype[which(is.na(aa$celltype))] = 'Unknown'

aa$celltype[which(aa$subtype == 'CM_Cav3.1')] = 'Unknown'

library(SeuratDisk)
VariableFeatures(aa)
#Idents(aa) = droplevels(aa$condition)

saveFile = paste0(outDir, 'nm_mm_ax_scRNAseq_merged_celllabels_v2.h5Seurat')
SaveH5Seurat(aa, filename = saveFile, overwrite = TRUE)
Convert(saveFile, dest = "h5ad", overwrite = TRUE)


########################################################
########################################################
# Section IV: human scRNA-seq data
# 
########################################################
########################################################

