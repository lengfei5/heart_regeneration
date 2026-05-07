##########################################################################
##########################################################################
# Project: heart regeneration
# Script purpose: process the newly mapped axolotl snRNA-seq data with UKYAmexF1
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon May 4 15:12:47 2026
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_R13591_intron.exon_UKYAmexF1_20260504'

resDir = paste0("../results/snRNAseq_axolotl", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R13591_axolotl_multiome/Amex_snRNA_mapping_UKYAmexF1'

source('functions_scRNAseq.R')
source('functions_Visium.R')
species = 'axloltl_scRNAseq'

require(Seurat)
#require(sctransform)3
library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
require(tictoc)
options(future.globals.maxSize = 80000 * 1024^2)

mem_used()


########################################################
########################################################
# Section I: import the scRNAseq data by kalisto
# check QCs 
########################################################
########################################################
source('functions_scRNAseq.R')

design = data.frame(sampleID = seq(197249, 197253), 
                    condition = c(paste0('Amex_scRNA_d', c(0, 1, 4, 7, 14))), stringsAsFactors = FALSE)

for(n in 1:nrow(design))
{
  # n = 1
  cat('-----------', design$condition[n], '-------------\n')
  
  # load nf output and process
  topdir = paste0(dataDir, '/', design$condition[n], '/', design$condition[n], '/')
  
  aa = make_SeuratObj_scRNAseq(topdir = topdir,
                               saveDir = paste0(resDir, '/', design$condition[n], '_', design$sampleID[n], '/'), 
                               changeGeneName.axolotl = FALSE, 
                               defaultDrops.only = TRUE)
  
  aa$condition = design$condition[n]
  
  if(n == 1) {
    scn = aa
  }else{
    scn = merge(scn, aa)
  }
  
}

rm(aa)

save(design, scn, 
     file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, version.analysis, '.Rdata'))


##########################################
# General QCs and gene filtering
# 
##########################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, version.analysis, '.Rdata'))

# scn = aa

## filter genes here 
aa = CreateSeuratObject(counts = scn@assays$RNA@counts[, scn$iscell_dd],
                        meta.data = scn@meta.data[scn$iscell_dd, ], 
                        assay = 'RNA',
                        min.cells = 20, 
                        min.features = 50)

rm(scn)

# Cell QC metrics: percentage of Mt, nb of counts, nb of genes 
# get MT% (genes from UKY_Amex annotation)
mtgenes =  read.table(file = paste0('/groups/tanaka/People/current/jiwang/Genomes/axolotl_new/',
                                    'UKY_AmexF1_1/annotations/RefSeq/',
                                    'UKY_AmexF1_mitochondrion_genes.txt'), 
                    header = TRUE, sep = '\t')

mtgenes = unique(mtgenes$gene)

mtgenes = mtgenes[!is.na(match(mtgenes, rownames(aa)))]

xx = PercentageFeatureSet(aa, col.name = "percent.mt", assay = "RNA", features = mtgenes)
aa[['percent.mt']] = xx$percent.mt
rm(xx)

Idents(aa) = aa$condition


pdfname = paste0(resDir, '/QCs_nFeatures_nCounts_percentMT',  version.analysis,  '_v2.pdf')
pdf(pdfname, width=12, height = 8)

levels = design$condition

table(aa$condition) %>%
  as.data.frame() %>%
  as_tibble() %>%
  mutate(condition = factor(Var1, levels=levels)) %>%
  #mutate(cellNbs = integer(Freq))
  ggplot(aes(x=condition, y=Freq, fill = condition)) +
  geom_bar(stat="identity", width=0.5) +
  theme_classic() + 
  scale_fill_brewer(palette="Dark2")+
  labs( x = '', y = 'detected cell # from cellRanger barcodes' )  +
  theme(axis.text.x = element_text(angle = 0, size = 10))

VlnPlot(aa, features = 'nFeature_RNA', y.max = 10000) +
  geom_hline(yintercept = c(500, 2500, 3000))

VlnPlot(aa, features = 'nCount_RNA', y.max = 50000)
VlnPlot(aa, features = 'percent.mt', y.max = 100)

FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "percent.mt")

dev.off()


## second time cell filtering 
aa <- subset(aa, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 60)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_QCs_cellFiltered.rds')) 