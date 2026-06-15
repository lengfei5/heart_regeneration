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

VlnPlot(aa, features = 'nFeature_RNA', y.max = 12000) +
  geom_hline(yintercept = c(1000, 2500, 3000))

VlnPlot(aa, features = 'nCount_RNA', y.max = 50000)
VlnPlot(aa, features = 'percent.mt', y.max = 20)

FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "percent.mt")

dev.off()


## second time cell filtering 
aa <- subset(aa, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 & percent.mt < 20)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_QCs_cellFiltered.rds')) 


##########################################
# import the cell annotation from the previous analysis
# keep only cells in the final version
##########################################
refs_file = paste0('/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_subtypes_final_20221117.rds')

refs = readRDS(file = refs_file)
table(refs$subtypes)
length(table(refs$subtypes))

refs$subtypes = droplevels(refs$subtypes) 
length(table(refs$subtypes)) 
## only 41 subtype annotations from Elad, with additional annotation "doublet" with 0 cell

## prepare the celltype to use and also specify the time-specific subtypes
refs$celltype_toUse = as.character(refs$subtypes)
length(table(refs$celltype_toUse))

refs$condition = gsub('_scRNA', '', refs$condition)
refs$celltype_toUse = gsub('Mo/Macs', 'Mo.Macs', refs$celltype_toUse)
refs$celltype_toUse = gsub("[(]", '', refs$celltype_toUse)
refs$celltype_toUse = gsub("[)]", '', refs$celltype_toUse)

table(refs$celltype_toUse)
length(table(refs$celltype_toUse))

aa$condition = gsub('scRNA_','', aa$condition)
aa$cellid = colnames(aa)
aa$cellid = sapply(aa$cellid, function(x){unlist(strsplit(as.character(x), '-'))[1]})
aa$cellid = paste0(aa$condition, '_', aa$cellid)

refs$cellid = colnames(refs)
refs$cellid = sapply(refs$cellid, function(x){unlist(strsplit(as.character(x), '-'))[1]})
refs$cellid = paste0(refs$condition, '_', refs$cellid)

mm = match(aa$cellid, refs$cellid)

aa$subtypes = refs$celltype_toUse[mm]

# define coarse clusters
aa$celltypes = as.character(aa$subtypes)

aa$celltypes[grep('CM_|CMs_|_CM|_CM_', aa$subtypes)] = 'CM'
aa$celltypes[grep('EC_|_EC', aa$subtypes)] = 'EC'
aa$celltypes[grep('FB_', aa$subtypes)] = 'FB'
aa$celltypes[grep('B_cells', aa$subtypes)] = 'Bcell'
aa$celltypes[grep('T_cells', aa$subtypes)] = 'Tcell'

aa$celltypes[grep('Macrophages|_MF|Mo.Macs_', aa$subtypes)] = 'Macrophages'
aa$celltypes[grep('Megakeryocytes', aa$subtypes)] = 'Megakeryocytes'
aa$celltypes[grep('RBC', aa$subtypes)] = 'RBC'

aa$celltypes[grep('Neu_', aa$subtypes)] = 'Neu'

aa = subset(aa, cells = colnames(aa)[which(!is.na(aa$subtypes))])

VlnPlot(aa, features = 'percent.mt', y.max = 10)

aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000)
all.genes <- rownames(aa)

aa <- ScaleData(aa, features = all.genes)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)
ElbowPlot(aa, ndims = 50)

aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

p1 = DimPlot(aa, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()
p2 = DimPlot(aa, label = TRUE, group.by = 'celltypes',  repel = TRUE) + NoLegend()

p1 + p2

ggsave(filename = paste0(resDir, '/umap_Elad_doubletRM_cleaned_manualAnnot.pdf'), width = 20, height = 8)

saveRDS(aa, file = paste0(RdataDir, 'aa_annotated_no_doublets_Elad.rds'))
