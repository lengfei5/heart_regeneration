##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: analyze the bulk RNA-seq from Plasmidsaurus
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Oct  9 09:53:27 2026
##########################################################################
##########################################################################
rm(list = ls())

require(openxlsx)
require(ggplot2)
require(DESeq2)
require(dplyr)
require(gridExtra)

RNA.functions = '/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)
#source('Functions_rnaseq.R')

version.Data = 'bulk_dnAXL'
version.analysis = paste0("_", version.Data, "_20261009")

## Directories to save results 
dataDir = "/groups/tanaka/People/current/jiwang/projects/heart_regeneration/bulkRNAseq_dnAXL_KDTXBX"
#design.file = paste0(dataDir, 'sampleInfos_parsed.txt')

resDir = paste0(dataDir, "/results/")
#tabDir =  paste0(resDir, "/tables/")

RdataDir = paste0(resDir, "Rdata/")
tableDir = paste0(resDir, 'tables4plots/')

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tableDir)){dir.create(tableDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}


##########################################
# import the tables from featurecount
##########################################
xlist = list.files(path=paste0(dataDir, '/featurecounts_Q10'),
                   pattern = "*_featureCounts.txt$", full.names = TRUE) ## list of data set to merge

#ggs = read.table(file = xlist[1], header = TRUE, sep = '\t')

all = cat.countTable(xlist, countsfrom = 'featureCounts')

colnames(all)[2:ncol(all)] = gsub('.sort_umiDedup_featureCounts.txt', '', colnames(all)[2:ncol(all)])
#all = all[, c(1:3)]

all = data.frame(all, stringsAsFactors = FALSE)

#colnames(all)[2:3] = c('187242_nonUMI', '187242_UMI')
#all = data.frame(all, ggs[match(all$gene, ggs$Geneid), ], stringsAsFactors = FALSE)

annot_all = readRDS(file = paste0('/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                                  'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))
annot = readRDS(paste0('/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))

#mm2 = match(all$gene, annot_all$geneID)
#all = data.frame(all, annot_all$chr_gene[mm2], annot_all$start_gene[mm2], annot_all$end_gene[mm2], stringsAsFactors = FALSE)
mm = match(all$gene, annot$geneID)

ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
all$gene[!is.na(mm)] = ggs[!is.na(mm)]

##########################################
# make DESeq object, QCs
##########################################
#all$total = apply(as.matrix(all[, c(2:3)]), 1, sum)
#all = all[order(-all$total), ]
raw = as.matrix(all[, -1])
rownames(raw) = all$gene

design = data.frame(sample = colnames(raw), condition = c(rep('Ctrl', 3), rep('dnAXL', 3)))
                    
#design$batch = paste0(design$request, '_', design$protocol)
rownames(design) = colnames(raw)

dds <- DESeqDataSetFromMatrix(raw, DataFrame(design), design = ~ condition)

saveRDS(dds, file = paste0(RdataDir, 'dds_dnAXL_plasmidsaurus.rds'))

ss = rowSums(counts(dds))

hist(log10(ss), breaks = 100);abline(v = log10(20), lwd = 2.0, col = 'red')
cat(length(which(ss>10)), ' gene selected \n')

dds = dds[which(ss>10), ]

dds = estimateSizeFactors(dds)

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = c('condition'), returnData = FALSE)
print(pca)


saveRDS(dds, file = paste0(RdataDir, 'dds_normalized_dnAXL_plasmidsaurus.rds'))

dds = readRDS(file = paste0(RdataDir, 'dds_normalized_dnAXL_plasmidsaurus.rds'))

dds$condition <- relevel(dds$condition, ref = "Ctrl")

dds <- estimateDispersions(dds, fitType = 'parametric')
plotDispEsts(dds, ymin = 10^-3); abline(h = 0.1, col = 'blue', lwd = 2.0)

dds = nbinomWaldTest(dds, betaPrior = TRUE)
resultsNames(dds)

res = results(dds, contrast=c("condition", 'dnAXL', 'Ctrl'), alpha = 0.1)
colnames(res) = paste0(colnames(res), "_dnAXL.vs.Ctrl")
res = data.frame(res[, c(2, 5, 6)])

cpm = fpm(dds)

res = data.frame(cpm, res, stringsAsFactors = FALSE)

write.csv2(res,
          file = paste0(tableDir, 'dnAXL_plasmidsaurus_cpm_DEgenes.csv'), 
          row.names = TRUE)


##########################################
# compare those two samples with CMs from scRNA-seq data 
##########################################
library(Seurat)
#library(decoupleR)
library(tictoc)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(data.table)

source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')

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


Idents(refs) = refs$celltype_toUse

res = read.csv2(file = paste0(tableDir, 'dnAXL_plasmidsaurus_cpm_DEgenes.csv'), row.names = c(1))
cpm = res[, c(1:6)]

ggs_id2 = sapply(rownames(cpm), 
                 function(x) {test = unlist(strsplit(as.character(x), '_')); return(test[length(test)])})
rownames(cpm) = ggs_id2


markers = FindMarkers(refs, ident.1 = 'CM_Prol_IS', ident.2 = "CM_Robo2",
                      logfc.threshold = 0.25,
                      test.use = "wilcox",
                      min.pct = 0.1,
                      only.pos = TRUE)

ggs_id = get_geneID(rownames(markers))
mm = match(ggs_id, rownames(cpm))

markers = markers[!is.na(mm), ]


ntop = 500
ggs = rownames(markers)[1:ntop]


cm_pseudo = AverageExpression(refs, assays = "RNA", features = ggs,
                              group.by = c('celltype_toUse'), 
                              layer = 'data')
cm_pseudo = data.frame(cm_pseudo)
cm_pseudo = cm_pseudo[, grep('CM_Prol_IS|CM_Robo2', colnames(cm_pseudo))]

rownames(cm_pseudo) = get_geneID(rownames(cm_pseudo))

mm = match(rownames(cm_pseudo), rownames(cpm))
cpm_2 = cpm[mm, ]


cors = cor(x = cpm_2, y = cm_pseudo, method = 'spearman')

#res = data.frame(cors = cors, condition = c(rep('Ctrl', 3), rep('dnAXL', 3)))

