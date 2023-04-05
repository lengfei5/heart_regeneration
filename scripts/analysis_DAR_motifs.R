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

###  import annotation and metadata
library(BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M)
library(ballgown)
gtf_axolotl = paste0("/groups/tanaka/People/current/jiwang/scripts/axolotl_multiome/r_package/", 
                     "AmexT_v47.FULL_corr_chr_cut.gtf")

granges_axolotl = ballgown::gffReadGR(gtf_axolotl)

# adding a gene biotype, as that's necessary for TSS metaprofile
granges_axolotl$gene_biotype = "protein_coding"

species = 'axloltl_scATAC'

########################################################
########################################################
# Section :  import data  
# process the scATAC-seq data: add UMAP from topic model 
# # update the celltypes and subtypes

########################################################
########################################################
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
  
  
  ##########################################
  # add topic umap  
  ##########################################
  ## seurat object after filtering
  srat_cr = readRDS(file = paste0(RdataDir, 
                                  'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                  'scATAC.merged.peaks.cr_filtered_lsi.umap_',
                                  '584K.annot_38280cells.rds'))
  
  seurat_obj = readRDS(file = paste0(RdataDir, 
                                     'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                     'scATAC.merged.peaks.cr_584K.features_38247cells_200topics.rds')) 
  
  umap.topics = seurat_obj@reductions$umap@cell.embeddings
  
  colnames(umap.topics) = paste0('umaptopic_', 1:ncol(umap.topics))
  dimensions = ncol(umap.topics)
  
  srat_cr = subset(srat_cr, cells = rownames(umap.topics))
  
  umap.topics = umap.topics[match(colnames(srat_cr), rownames(umap.topics)), ]
  srat_cr[['umap_topics']] = Seurat::CreateDimReducObject(embeddings=umap.topics,
                                                          key='UMAPTOPICS_',
                                                          assay='ATAC')
  rm(seurat_obj)
  
  saveRDS(srat_cr, file = paste0(RdataDir, 
                                 'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                 'scATAC.merged.peaks.cr_filtered_umap.lsi',
                                 '584K.features_38247cells_umap.topics.rds'))
  
  
  ##########################################
  # update the celltype, subtypes and umap from Elad 
  ##########################################
  srat_cr = readRDS(file = paste0(RdataDir, 
                                  'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                  'scATAC.merged.peaks.cr_filtered_umap.lsi',
                                  '584K.features_38247cells_umap.topics.rds'))
  
  aa = readRDS(file = 
                 '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets_2022_10_17.rds')
  
  Idents(aa) = aa$subtypes
  DimPlot(aa, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()
  
  cells = intersect(colnames(aa), colnames(srat_cr))
  srat_cr = subset(srat_cr, cells = cells)
  
  # replace the subtypes 
  srat_cr$subtypes = aa$subtypes[match(colnames(srat_cr), colnames(aa))]
  DimPlot(srat_cr, reduction = 'umap',  
          label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()
  
  # keep the same umap 
  umap.mat = aa@reductions$umap@cell.embeddings
  umap.mat = umap.mat[match(colnames(srat_cr), rownames(umap.mat)), ]
  #colnames(umap.topics) = paste0('umaptopic_', 1:ncol(umap.topics))
  #dimensions = ncol(umap.topics)
  srat_cr[['umap']] = Seurat::CreateDimReducObject(embeddings=umap.mat,
                                                          key='UMAP_',
                                                          assay='RNA')
  
  DefaultAssay(srat_cr) <- 'RNA'
  
  DimPlot(srat_cr, reduction = 'umap',  
          label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()
  
  saveRDS(srat_cr, file = paste0(RdataDir, 
                                 'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                 'scATAC.merged.peaks.cr_filtered_umap.lsi',
                                 '584K.features_37680cells_umap.topics_updated.umap.subtypes.rds'))
  
        
}

########################################################
########################################################
# Section II : Differential Accessible Regions (DAR analysis)
# 
########################################################
########################################################
srat_cr = readRDS(file = paste0(RdataDir, 
                                'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                'scATAC.merged.peaks.cr_filtered_umap.lsi',
                                '584K.features_37680cells_umap.topics_updated.umap.subtypes.rds'))



########################################################
########################################################
# Section III: Motif enrichment analysis  
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
                                'scATAC.merged.peaks.cr_filtered_umap.lsi',
                                '584K.features_37680cells_umap.topics_updated.umap.subtypes.rds'))

DefaultAssay(srat_cr) <- 'RNA'
DimPlot(srat_cr, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()

DefaultAssay(srat_cr) <- 'ATAC'
DimPlot(srat_cr, label = TRUE, reduction = 'umap_topics',
        group.by = 'subtypes',  repel = TRUE) + NoLegend()


# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', 
              all_versions = FALSE)
)


Save_motif_to_tfs = FALSE
if(Save_motif_to_tfs){
  motif_tfs = data.frame(motif = names(pfm), tf = name(pfm), stringsAsFactors = FALSE)
  motif_tfs$tf = toupper(motif_tfs$tf)
  motif_tfs$tf = gsub('::', '_', motif_tfs$tf)
  motif_tfs$tf = gsub("\\s*\\([^\\)]+\\)","",as.character(motif_tfs$tf))
  
  saveRDS(motif_tfs, file = paste0(RdataDir, 
                                   'motif_to_tfs_pfm_JASPAR2020_CORE_vertebrate.rds'))
  
}


Run_ChromVar = TRUE
if(Run_ChromVar){
  
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
  
  saveRDS(srat_cr, file = paste0(RdataDir, 'atac_seuratObject_motifClass_chromVAR_v3.rds'))
  
}


Plot.motif.activity.ChromVar = FALSE
if(Plot.motif.activity.ChromVar){
  
  aa = readRDS(file = paste0(RdataDir, 'atac_seuratObject_motifClass_chromVAR_v2.1.rds'))
    
  DefaultAssay(aa) <- 'chromvar'
  
  ss = colSums(aa@assays$chromvar@data)
  
  aa = subset(aa, cells = colnames(aa)[which(!is.na(ss))])
  
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

