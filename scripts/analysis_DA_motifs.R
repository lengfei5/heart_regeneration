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

design = data.frame(sampleID = seq(197254, 197258), 
                    condition = c(paste0('Amex_scATAC_d', c(0, 1, 4, 7, 14))), stringsAsFactors = FALSE)
design$timepoint = gsub('Amex_scATAC_', '', design$condition) 



##########################################
# import data  
##########################################
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


########################################################
########################################################
# Section : Find differentially accessible peaks and motif enrichment analysis  
# 
########################################################
########################################################
DefaultAssay(seurat.cistopic) <- 'peaks'

seurat.cistopic = compute.motif.enrichment(seurat.cistopic)

DefaultAssay(seurat.cistopic) <- 'chromvar'

# tbx-38, med-2, elt-3, pal-1
motifs2check = rownames(seurat.cistopic)[grep('tbx-35|tbx-38|MA0542.1|MA0924.1 pal-1|MA0546.1|end-1|MA0547.1 skn-1|hnd-1|unc-130|mom-2|med-2', 
                                              rownames(seurat.cistopic))]
FeaturePlot(
  object = seurat.cistopic,
  features = motifs2check,
  pt.size = 0.1,
  max.cutoff = 'q90',
  min.cutoff = 'q10',
  ncol = 3
)

saveRDS(seurat.cistopic, file =  paste0(RdataDir, 'atac_LDA_seurat_object_motifClassChromVar.rds'))



