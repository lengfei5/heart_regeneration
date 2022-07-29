##########################################################################
##########################################################################
# Project: heart regeneration
# Script purpose: analyze the scRNA-seq from scMultiome
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Jul 20 15:45:18 2022
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_R13591_intron.exon.20220729'

resDir = paste0("../results/sc_multiome", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R13591_axolotl_multiome'

source('functions_scRNAseq.R')
source('functions_Visium.R')
require(Seurat)
#require(sctransform)
library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
require(tictoc)
mem_used()

species = 'axloltl_scRNAseq'

########################################################
########################################################
# Section I: import the scRNAseq data by kalisto
# check QCs 
########################################################
########################################################
design = data.frame(sampleID = seq(197249, 197253), 
                    condition = c(paste0('Amex_scRNA_d', c(0, 1, 4, 7, 14))), stringsAsFactors = FALSE)

source('functions_scRNAseq.R')

for(n in 1:nrow(design))
{
  # n = 1
  cat('-----------', design$condition[n], '-------------\n')
  
  # load nf output and process
  topdir = paste0(dataDir, '/', design$condition[n], '/', design$condition[n], '/')
  
  aa = make_SeuratObj_scRNAseq(topdir = topdir,
                             saveDir = paste0(resDir, '/', design$condition[n], '_', design$sampleID[n], '/'), 
                             changeGeneName.axolotl = TRUE, 
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
# get MT% (genes curated from NCBI chrMT genes)
mtgenes = c("COX1", "COX2", "COX3", "ATP6", "ND1", "ND5", "CYTB", "ND2", "ND4", "ATP8", "MT-CO1", "COI", "LOC9829747")
mtgenes = c(mtgenes, paste0("MT", mtgenes), paste0("MT-", mtgenes))

ggs = sapply(rownames(aa), function(x) unlist(strsplit(as.character(x), '-'))[1])
mtgenes = rownames(aa)[!is.na(match(ggs, mtgenes))]

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

########################################################
########################################################
# Section : normalization and quick clustering
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_QCs_cellFiltered.rds')) 

Normalize_with_sctransform = FALSE

if(!Normalize_with_sctransform){
  
  aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # tic()
  # aa = Normalize_with_scran(aa)
  # toc()
  
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
  all.genes <- rownames(aa)
  aa <- ScaleData(aa, features = all.genes)
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)
  
  # saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_lognormamlized_pca.rds')) 
  
}else{
  
  tic()
  test <- SCTransform(aa, ncells = 3000, assay = "RNA", verbose = FALSE, 
                      variable.features.n = 3000, return.only.var.genes = TRUE, vst.flavor = "v2")
  toc()
  
  saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_SCTnormamlized.Rdata'))
  
  aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
}

# aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_lognormamlized_pca.rds')) 

ElbowPlot(aa, ndims = 30)

aa <- FindNeighbors(aa, dims = 1:30)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)

aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.3)
DimPlot(aa, label = TRUE, repel = TRUE) + ggtitle("Unsupervised clustering")

ggsave(filename = paste0(resDir, '/first_test_umap_v1.pdf'), width = 8, height = 6)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_lognormamlized_pca_umap_v2.rds'))

########################################################
########################################################
# Section :  cell type annotation  
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_lognormamlized_pca_umap_v2.rds'))

DimPlot(aa, label = TRUE, repel = TRUE) + ggtitle("scNuc (multiome)")

features = rownames(aa)[grep('LY6', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('PTPRC', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))
ggsave(filename = paste0(resDir, '/FeaturePlot_CD45.pdf'), width = 8, height = 6)

features = rownames(aa)[grep('VIM|COL1A2|FSTL1|POSTN', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('MYH6|ACTN2|NPPA|TNNT2|GATA4', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('CD68|CD8A|CD74|CSF1R|ITGAM', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('MKI67|CCNB2|PCNA-|CDK1-', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

#DimPlot(aa, reduction = "umap", label = TRUE, repel = TRUE, group.by = "nCount_RNA") +
#  NoLegend()

### cluster markers
# markers = FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# saveRDS(markers, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_markers_v2.rds')) 
markers = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_markers_v2.rds')) 

markers %>%
  filter(!str_detect(gene, '^(AMEX|LOC)')) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

saveRDS(top10, file = paste0(RdataDir, 'top10_markerGenes_coarseCluster.rds'))

xx = subset(aa, downsample = 500)

DoHeatmap(xx, features = top10$gene) + NoLegend()
ggsave(filename = paste0(resDir, '/first_test_clusterMarkers_v2.pdf'), width = 10, height = 30)

