##########################################################################
##########################################################################
# Project: Elad heart regeneration 
# Script purpose: analyze pilot bulk RNA-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Apr 20 11:13:02 2022
##########################################################################
##########################################################################
rm(list = ls())

RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)
source('Functions_rnaseq.R')
require(openxlsx)
require(ggplot2)
require(DESeq2)
require(dplyr)
require(gridExtra)


version.Data = 'bulk_smartseq3';
version.analysis = paste0("_", version.Data, "_20220423")

## Directories to save results 
dataDir = "/Volumes/groups/tanaka/People/current/jiwang/projects/heart_regeneration/R13267_smartseq3_pilot"
#design.file = paste0(dataDir, 'sampleInfos_parsed.txt')

resDir = paste0(dataDir, "/results/", version.Data)
#tabDir =  paste0(resDir, "/tables/")

RdataDir = paste0(resDir, "/Rdata/")
tableDir = paste0(resDir, 'tables4plots/')

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tableDir)){dir.create(tableDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

xlist = list.files(path=paste0(dataDir, '/featurecounts_Q30'),
                   pattern = "*_featureCounts.txt$", full.names = TRUE) ## list of data set to merge

#ggs = read.table(file = xlist[1], header = TRUE, sep = '\t')

all = cat.countTable(xlist, countsfrom = 'featureCounts')

all = all[, c(1:3)]

all = data.frame(all, stringsAsFactors = FALSE)

colnames(all)[2:3] = c('187242_nonUMI', '187242_UMI')

#all = data.frame(all, ggs[match(all$gene, ggs$Geneid), ], stringsAsFactors = FALSE)

annot_all = readRDS(file = paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                                  'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))
annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))

#mm2 = match(all$gene, annot_all$geneID)
#all = data.frame(all, annot_all$chr_gene[mm2], annot_all$start_gene[mm2], annot_all$end_gene[mm2], stringsAsFactors = FALSE)
mm = match(all$gene, annot$geneID)

ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
all$gene[!is.na(mm)] = ggs[!is.na(mm)]

#all = data.frame(all, stringsAsFactors = FALSE)
all$total = apply(as.matrix(all[, c(2:3)]), 1, sum)

all = all[order(-all$total), ]

#all$total = apply(as.matrix(all[, c(2:3)]), 1, sum)


all$cpm_nonUMI = all$`187242_nonUMI`/sum(all$`187242_nonUMI`)*10^6 
all$cpm.total = all$total/sum(all$total) * 10^6



colnames(all)[4:6] = c('chr', 'start', 'end')

aa = all[order(-all$X187242_nonUMI), ]
aa$start = log2(as.numeric(as.character(aa$start)))
aa$end = log2(as.numeric(as.character(aa$end)))

write.table(all, file = paste0(resDir, '/raw_readCounts_cpm_strand.specific.txt'), sep = '\t', 
            quote = FALSE, col.names = TRUE, row.names = FALSE)

##########################################
# compare the two different counting 
##########################################
aa = read.table(file = paste0(resDir, '/raw_counts_cpm_umi_nonumi_strand.specific.txt'), header = TRUE, sep = '\t')
bb = read.table(file = paste0(resDir, '/raw_counts_cpm_umi_nonumi.txt'), header = TRUE, sep = '\t')
aa = aa[order(-aa$X187242_nonUMI), ]

bb = bb[match(aa$gene, bb$gene),]
ss1 = apply(aa[, c(2:3)], 1, sum)
ss2 = apply(bb[, c(2,3)], 1, sum)

jj = which(ss1>0 | ss2>0)
aa = aa[jj, ]
bb = bb[jj, ]

plot(bb$X187242_nonUMI, aa$X187242_nonUMI, log ='xy');
abline(0, 1, lwd =2.0, col = 'red')

plot(bb$X187242_UMI, aa$X187242_UMI, log ='xy');
abline(0, 1, lwd =2.0, col = 'red')

plot(bb$X187242_UMI, bb$X187242_nonUMI, log ='xy');
abline(0, 1, lwd =2.0, col = 'red')


