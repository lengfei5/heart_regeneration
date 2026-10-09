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

