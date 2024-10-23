##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: process and analyze the mouse heart data from Torres group
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Jul 10 14:13:28 2024
##########################################################################
##########################################################################
rm(list = ls())

RNA.functions = '/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)

# setup for data import and sequencing QCs
version.analysis = 'Rxxxx_mouse_bulkRNAseq_polyA_2020710'

resDir = paste0("../results/", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../Rxxxx_bulkRNAseq_mouse/featurecounts_Q30/'

Collect.QCs.stat = FALSE

# sps = readRDS(file = paste0('../data/annotations/curated_signaling.pathways_gene.list_v3.rds'))
# sps = unique(sps$gene)
# 
# tfs = readRDS(file = paste0('../data/annotations/curated_human_TFs_Lambert.rds'))
# tfs = unique(tfs$`HGNC symbol`)
# 
# tfs = toupper(unique(c(tfs,sps)))
# 
# source('functions_utility.R')

########################################################
########################################################
# Section I : data processing and sequencing quality controls
# 
########################################################
########################################################

##################################################
## Import count table
##################################################
source(RNA.functions)
source(RNA.QC.functions)

xlist = list.files(path=dataDir,
                   pattern = "*featureCounts.txt$", full.names = TRUE) ## list of data set to merge

all = cat.countTable(xlist, countsfrom = 'featureCounts')

samples = colnames(all)
samples = samples[-1]
design = data.frame(t(sapply(samples, function(x){unlist(strsplit(as.character(x), '_'))[1:3]})))
rownames(design) = NULL

colnames(design)[c(1,2,3)] = c('SampleID', 'condition',  'treat')
design$SampleID = paste0(design$SampleID, '_', design$condition, '_', design$treat)

counts = process.countTable(all=all, design = design[, c(1, 3, 2)], merge.technicalRep.sameID = FALSE)

rownames(counts) = counts$gene
#counts = as.matrix(counts[, -1])

save(design, counts, file=paste0(RdataDir, 'design_rawCounts', version.analysis, '.Rdata'))

##########################################
# gene names converted from ensID to gene symbol 
##########################################
load(file=paste0(RdataDir, 'design_rawCounts', version.analysis, '.Rdata'))

colnames(counts)[-1] = design$SampleID

## convert gene names
annot = read.delim(paste0('/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/', 
                          'ens_BioMart_GRCm38.p6.txt'), sep = '\t', header = TRUE)

# annot = annot[which(annot$Gene.type == 'protein_coding'), ]
annot = annot[!is.na(match(annot$Chromosome.scaffold.name, as.character(c(1:19, 'X', 'Y')))), ]

mm = match(counts$gene, annot$Gene.stable.ID)
cat(length(which(!is.na(mm))), ' genes in the count table \n')

counts = counts[!is.na(mm), ]

rownames(counts) = counts$gene

counts$gene = annot$Gene.name[match(counts$gene, annot$Gene.stable.ID)]

counts = counts[!is.na(counts$gene), ]
gg.uniq = unique(counts$gene)
counts = counts[match(gg.uniq, counts$gene), ]
rownames(counts) = counts$gene
counts = counts[, -1]

save(design, counts, file=paste0(RdataDir, 'design_rawCounts_geneSymbols_', version.analysis, '.Rdata'))

##########################################
# QCs of replicates and conditions
##########################################
load(file=paste0(RdataDir, 'design_rawCounts_geneSymbols_', version.analysis, '.Rdata'))

design$condition = paste0(design[, c(2)], "_", design[,3])



QC.for.cpm = FALSE
if(QC.for.cpm){
  
  source(RNA.functions)
  source(RNA.QC.functions)
  
  raw = as.matrix(counts)
  
  write.csv2(raw, file = paste0(resDir, '/countTable_geneSymbol.csv'), row.names = TRUE, 
             col.names = TRUE, quote = FALSE)
  
  #kk = which(design$SampleID != '161040' & design$condition != 'N2B27')
  #raw = raw[, -kk]
  
  ss = apply(as.matrix(raw), 1, sum)
  raw = raw[which(ss > 0), ]
  
  pdfname = paste0(resDir, "/Data_qulity_assessment", version.analysis, ".pdf")
  pdf(pdfname, width = 20, height = 16)
  
  Check.RNAseq.Quality(read.count=raw, design.matrix = design[ , c(1, 3)], 
                       lowlyExpressed.readCount.threshold=20)
  
  dev.off()
  
  
}

########################################################
########################################################
# Section II : normalization and pairwise comparison with DESeq2 
# 
########################################################
########################################################
require(DESeq2)
require(ggplot2)
library(ggrepel)
require(gridExtra)
library(dplyr)
library(ggrepel)
library(patchwork)
require(pheatmap)
library(org.Mm.eg.db)
library(enrichplot)
library(clusterProfiler)
library(stringr)


load(file = paste0(RdataDir, 'design_rawCounts', version.analysis, '.Rdata'))
design$condition = paste0(design[, c(2)], "_", design[,3], "_", design[,4], '_', design[,5])

dds <- DESeqDataSetFromMatrix(counts, DataFrame(design), design = ~ condition)

ss = rowSums(counts(dds))

hist(log10(ss), breaks = 100, main = 'log10(sum of reads for each gene)')

cutoff.gene = 20
cat(length(which(ss > cutoff.gene)), 'genes selected \n')

dds <- dds[ss > cutoff.gene, ]

# normalization and dimensionality reduction
dds = estimateSizeFactors(dds)
fpm = fpm(dds, robust = TRUE)

#save(dds, design.matrix,  file = paste0(RdataDir, 'TM3_dds_normalized.Rdata'))
#save(fpm, design, file = paste0(tfDir, '/RNAseq_fpm_fitered.cutoff.', cutoff.gene, '.Rdata'))
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = c('timepoint', 'Wash', 'Culture', 'RA'), returnData = FALSE)
print(pca)

pca2save = as.data.frame(plotPCA(vsd, intgroup = c('timepoint', 'Wash', 'Culture', 'RA'), 
                                 returnData = TRUE, ntop = 3000))
pca2save$time = as.factor(pca2save$timepoint)
pca2save$cc = paste0(pca2save$Wash, '_', pca2save$Culture, '_', pca2save$RA)

ggplot(data=pca2save, aes(PC1, PC2, label = name, color= cc, shape = time))  + 
  geom_point(size=3) + 
  scale_shape_manual(values=1:nlevels(pca2save$time)) +
  #geom_text(hjust = 0.5, nudge_y = 0.5, size=3.0) + 
  geom_text_repel(data= pca2save, size = 3) 

ggsave(paste0(resDir, '/PCA_quantseq_allSamples.pdf'),  width=20, height = 14)

save(dds, design, file = paste0(RdataDir, 'design_dds_', version.analysis, '.Rdata'))