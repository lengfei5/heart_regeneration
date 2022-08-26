##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: analyze scATAC-seq data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Aug 25 14:24:06 2022
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_R13591_atac.20220825'

resDir = paste0("../results/sc_multiome", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')


if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R13591_axolotl_multiome'

source('functions_scATAC.R')
source('functions_scRNAseq.R')
source('functions_Visium.R')


library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(patchwork)

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
require(tictoc)
options(future.globals.maxSize = 64000 * 1024^2)
set.seed(1234)
mem_used()

library(data.table)


##########################################
# axoltol genome and annotation (fragmented version) 
##########################################

# this (bespoke) package hosts the axolotl genome (comments from Tomas)
# package is very large, cannot be installed in home directory and should be in a fast disk
#install.packages("/groups/tanaka/People/current/jiwang/scripts/axolotl_multiome
# /r_package/BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M_1.0.0.tar.gz",
#                 repos = NULL, type = "source")
library(BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M)

library(ballgown)
gtf_axolotl = "/groups/tanaka/People/current/jiwang/scripts/axolotl_multiome/r_package/AmexT_v47.FULL_corr_chr_cut.gtf"
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

########################################################
########################################################
# Section I: 1) check QCs for individual sample
# 2) recall peaks for eahc sample and save it
# 3) peak consensus and recount 
########################################################
########################################################
design = data.frame(sampleID = seq(197254, 197258), 
                    condition = c(paste0('Amex_scATAC_d', c(0, 1, 4, 7, 14))), stringsAsFactors = FALSE)
design$timepoint = gsub('Amex_scATAC_', '', design$condition) 

source('functions_scRNAseq.R')

for(n in 1:nrow(design))
{
  # n = 1
  cat('-----------', design$condition[n], '-------------\n')
  
  # load nf output and process
  topdir = paste0(dataDir, '/multiome_', design$timepoint[n], '/outs')
  
  counts <- Read10X_h5(filename = paste0(topdir, "/filtered_feature_bc_matrix.h5"))
  fragpath <- paste0(topdir, "/atac_fragments.tsv.gz")
  metadata <- read.csv(
    file = paste0(topdir, '/per_barcode_metrics.csv'),
    header = TRUE,
    row.names = 1
  )
  
  
  chrom_assay <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = granges_axolotl,
    min.cells = 10,
    min.features = 50
  )
  
  bb <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    
  )
  rm(chrom_assay)
  
  bb$condition = design$condition[n]
  mm = match(colnames(bb), metadata$gex_barcode)
  bb$atac_barcode = metadata$atac_barcode[mm]
  bb$total_fragments = metadata$atac_fragments[mm]
  bb$peak_region_fragments = metadata$atac_peak_region_fragments[mm]
  bb$pct_reads_in_peaks <- bb$peak_region_fragments / bb$total_fragments * 100
  
  DefaultAssay(bb) <- "ATAC"
  
  bb <- NucleosomeSignal(bb)
  bb <- TSSEnrichment(bb, fast = FALSE)
  
  bb$high.tss <- ifelse(bb$TSS.enrichment > 5, 'High', 'Low')
  TSSPlot(bb, group.by = 'high.tss') + NoLegend()
  
  Idents(bb) = bb$condition
  VlnPlot(bb, 
          features = c("nCount_ATAC", "total_fragments", "peak_region_fragments"), 
          ncol = 3) 
  
  geom_hline(yintercept = c(500, 1000))
  
  VlnPlot(bb, features = 'TSS.enrichment', y.max = 6)
  
  VlnPlot(
    object = bb,
    features = c("TSS.enrichment", "nucleosome_signal"),
    ncol = 2,
    pt.size = 0
  )
  
}


rm(aa)

save(design, scn, 
     file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, version.analysis, '.Rdata'))


