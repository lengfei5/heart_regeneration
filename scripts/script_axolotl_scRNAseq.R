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

version.analysis = '_R13591_20220720'

resDir = paste0("../results/sc_multiome", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R13591_axolotl_multiome'

source('functions_scRNAseq.R')
source('functions_Visium.R')
library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
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

check.QC.each.condition = TRUE

for(n in 1:nrow(design))
{
  # n = 1
  cat('-----------', design$condition[n], '-------------\n')
  
  # load nf output and process
  topdir = paste0(dataDir, '/', design$condition[n], '/')
  
  aa = make_SeuratObj_scRNAseq(topdir = topdir,
                             saveDir = paste0(resDir, '/', design$condition[n], '_', design$sampleID[n], '/'), 
                             keyname = design$condition[n], 
                             changeGeneName.axolotl = FALSE, 
                             QC.umi = TRUE)
  
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
aa = CreateSeuratObject(counts = scn@assays$RNA@counts[, aa$iscell_dd],
                   meta.data = scn@meta.data[aa$iscell_dd, ], 
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
# mtgenes = mtgenes[mtgenes %in% g[,1]]
# srat = PercentageFeatureSet(srat, col.name = "percent.mt", assay = "Spatial",
#                             features = mtgenes)
xx = PercentageFeatureSet(aa, col.name = "percent.mt", assay = "RNA", features = mtgenes)
aa[['percent.mt']] = xx$percent.mt

rm(xx)

Idents(aa) = aa$condition

pdfname = paste0(resDir, '/QCs_gene_marker_check_', design$condition[n], version.analysis,  '.pdf')
pdf(pdfname, width=16, height = 8)

VlnPlot(aa, features = 'nFeature_RNA', y.max = 5000)
VlnPlot(aa, features = 'nCount_RNA', y.max = 10000)
VlnPlot(aa, features = 'percent.mt', y.max = 100)

FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "percent.mt")

##########################################
## filter cells again
##########################################
aa <- subset(aa, subset = nFeature_RNA > 100 & nFeature_RNA < 3000 & percent.mt < 60)

Normalize_with_sctransform = FALSE

########################################################
########################################################
# Section : normalization and quick clustering
# 
########################################################
########################################################

if(Normalize_with_sctransform){
  aa <- SCTransform(aa, assay = "RNA", verbose = FALSE, 
                    variable.features.n = 3000, return.only.var.genes = FALSE)
  aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
  
}else{
  # aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
  
  aa = Normalize_with_scran(aa)
  
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)
  all.genes <- rownames(aa)
  aa <- ScaleData(aa, features = all.genes)
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)
   
}

ElbowPlot(aa, ndims = 30)

aa <- FindNeighbors(aa, dims = 1:10)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)

aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 30, min.dist = 0.01)
DimPlot(aa, label = TRUE, repel = TRUE) + ggtitle("Unsupervised clustering")

ggsave(filename = paste0(resDir, '/first_test_umap.pdf'), width = 8, height = 6)


saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_normamlized_clustered_umap.Rdata'))

features = rownames(aa)[grep('VIM|COL1A2|FSTL1|POSTN', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('MYH6|ACTN2|NPPA|TNNT2|GATA4', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('CD68|CD8A|CD74|CSF1R|ITGAM', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('MKI67|CCNB2|PCNA-|CDK1-', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

dev.off()

FeaturePlot(aa, )

DimPlot(aa, reduction = "umap", label = TRUE, repel = TRUE, group.by = "nCount_RNA") +
  NoLegend()

markers = FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

markers %>%
  filter(!str_detect(gene, '^(AMEX|LOC)')) %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC) -> top10


DoHeatmap(aa, features = top10$gene) + NoLegend()

ggsave(filename = paste0(resDir, '/first_test_clusterMarkers.pdf'), width = 6, height = 12)
