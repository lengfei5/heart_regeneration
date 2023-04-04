##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: scATAC-seq analysis: DA analysis and motif analysis
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Mar 21 10:27:06 2023
##########################################################################
##########################################################################

rm(list = ls())
version.analysis = '_R13591_atac_reseq_20221115'

resDir = paste0("../results/sc_multiome", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R14353_ax_snATAC_reseq'

source('functions_scATAC.R')
source('functions_scRNAseq.R')
source('functions_Visium.R')
#source('functions_scRNAseq.R')

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(patchwork)
require(SeuratObject)
library(data.table)

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
require(tictoc)
library(future)
options(future.globals.maxSize = 80 * 1024^3)
set.seed(1234)
mem_used()

##########################################
# import annotation and metadata
##########################################
library(BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M)

library(ballgown)
gtf_axolotl = paste0("/groups/tanaka/People/current/jiwang/scripts/axolotl_multiome/r_package/", 
                     "AmexT_v47.FULL_corr_chr_cut.gtf")

granges_axolotl = ballgown::gffReadGR(gtf_axolotl)
# adding a gene biotype, as that's necessary for TSS metaprofile
granges_axolotl$gene_biotype = "protein_coding"

# with need to add the "proper" gene name
# basedir = "/links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/"
# gene_name_match = read.table(paste0(basedir, "AmexT_v47.FULL_t2g_note.txt"), sep = "\t")[,2:3]
# gene_name_match = gene_name_match[!duplicated(gene_name_match$V2), ]
# rownames(gene_name_match) = gene_name_match$V2
# newgenenames = gene_name_match[granges_axolotl$gene_id,2]
# granges_axolotl$gene_name = newgenenames

species = 'axloltl_scATAC'

#design = data.frame(sampleID = seq(197254, 197258), 
#                    condition = c(paste0('Amex_scATAC_d', c(0, 1, 4, 7, 14))), stringsAsFactors = FALSE)
#design$timepoint = gsub('Amex_scATAC_', '', design$condition) 

##########################################
# import data  
##########################################
Processing.scATAC = FALSE
if(Processing.scATAC){
  srat_cr = readRDS(file = paste0(RdataDir, 
                                  'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                  'scATAC.merged.peaks.cr.',
                                  '584K.annot_38280cells.rds'))
  # normalize ATAC and UMAP
  DefaultAssay(srat_cr) <- "ATAC"
  
  Filter.cells.with.scATAC = FALSE
  if(Filter.cells.with.scATAC){
    # quick filtering 
    srat_cr <- subset(
      x = srat_cr,
      subset = nCount_ATAC < 100000 &
        nCount_RNA < 25000 &
        nCount_ATAC > 200 &
        nCount_RNA > 1000 &
        nucleosome_signal < 6 &
        TSS.enrichment > 1
    )
  }
  
  Idents(srat_cr) = as.factor(srat_cr$condition)
  
  srat_cr <- RunTFIDF(srat_cr)
  srat_cr = FindTopFeatures(srat_cr, min.cutoff = 'q5')
  srat_cr <- RunSVD(srat_cr)
  
  DepthCor(srat_cr, n = 30)
  
  cordat = DepthCor(srat_cr, reduction = "lsi", n = 30)$data
  dims_use = cordat$Component[abs(cordat$counts)<0.3]
  
  dims_use = c(2:30)
  print(dims_use)
  
  srat_cr = FindNeighbors(object = srat_cr, reduction = 'lsi', dims = dims_use, 
                          force.recalc = T, graph.name = "thegraph")
  srat_cr = FindClusters(object = srat_cr, verbose = FALSE, algorithm = 3, 
                         graph.name = "thegraph", resolution = 1)
  
  srat_cr <- RunUMAP(object = srat_cr, reduction = 'lsi', dims = dims_use, n.neighbors = 30, min.dist = 0.1, 
                     reduction.name = "umap_lsi")
  
  DimPlot(object = srat_cr, label = TRUE, reduction = 'umap_lsi', group.by = 'celltypes') + NoLegend()
  
  p1 = DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap', group.by = 'subtypes') + 
    NoLegend() + ggtitle('snRNA-seq')
  p2 = DimPlot(object = srat_cr, label = TRUE, repel = TRUE, group.by = 'subtypes') + 
    NoLegend() + ggtitle('scATAC-seq')
  
  p1 + p2
  
  ggsave(filename = paste0(resDir, '/multiome_snRNA_umap_snATAC_lsi.umap.pdf'), height = 8, width = 20)
  
  
  saveRDS(srat_cr, file = paste0(RdataDir, 
                                 'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                 'scATAC.merged.peaks.cr_filtered_lsi.umap_',
                                 '584K.annot_38280cells.rds'))
    
}

########################################################
########################################################
# Section : Find differentially accessible peaks and motif enrichment analysis  
# 
########################################################
########################################################
library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(chromVAR)
#library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M)
library(patchwork)
set.seed(1234)

srat_cr = readRDS(file = paste0(RdataDir, 
              'seuratObj_multiome_snRNA.annotated.normalized.umap_',
              'scATAC.merged.peaks.cr_filtered_lsi.umap_',
              '584K.annot_38280cells.rds'))

DefaultAssay(srat_cr) <- 'ATAC'

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', 
              all_versions = FALSE)
)

# add motif information
srat_cr <- AddMotifs(
  object = srat_cr,
  genome = BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M,
  pfm = pfm
  
)

library(BiocParallel)
#register(MulticoreParam(32)) # Use multiple cores
register(SerialParam())

srat_cr <- RunChromVAR(
  object = srat_cr,
  genome = BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M
)

saveRDS(srat_cr, file = paste0(RdataDir, 'atac_seuratObject_motifClass_chromVAR_v2.1.rds'))

Plot.motif.activity.ChromVar = FALSE
if(Plot.motif.activity.ChromVar){
  
  srat_cr = readRDS(file = paste0(RdataDir, 'atac_seuratObject_motifClass_chromVAR_v2.1.rds'))
    
  DefaultAssay(srat_cr) <- 'chromvar'
  
  ss = colSums(srat_cr@assays$chromvar@data)
  
  srat_cr = subset(srat_cr, cells = colnames(srat_cr)[which(!is.na(ss))])
  
  # look at the activity of Mef2c
  
  p2 <- FeaturePlot(
    object = srat_cr,
    #features = "MA0019.1",
    features = "MA0497.1",
    min.cutoff = 'q5',
    max.cutoff = 'q95',
    pt.size = 0.1
  )
  
  p1 + p2
  
  
  
  
}







