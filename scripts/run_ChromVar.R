##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: run ChromVar for scATAC-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Apr  5 11:22:56 2023
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

#DefaultAssay(srat_cr) <- 'RNA'
#DimPlot(srat_cr, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()

DefaultAssay(srat_cr) <- 'ATAC'
#DimPlot(srat_cr, label = TRUE, reduction = 'umap_topics',
#        group.by = 'subtypes',  repel = TRUE) + NoLegend()


# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', 
              all_versions = FALSE)
)


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