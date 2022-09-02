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
# Section I: check QCs for individual sample from the cellranger-arc output 

########################################################
########################################################
design = data.frame(sampleID = seq(197254, 197258), 
                    condition = c(paste0('Amex_scATAC_d', c(0, 1, 4, 7, 14))), stringsAsFactors = FALSE)
design$timepoint = gsub('Amex_scATAC_', '', design$condition) 

source('functions_scRNAseq.R')
srat = list()

for(n in 2:nrow(design))
{
  # n = 5
  cat('----------- : ', n, ':',  design$condition[n], '-------------\n')
  
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
    min.features = 20
  )
  
  bb <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    
  )
  rm(chrom_assay)
  
  bb$condition = design$condition[n]
  
  mm = match(colnames(bb), metadata$gex_barcode)
  metadata = metadata[, c(1,2, 3, 21:30)]
  
  bb = AddMetaData(bb, metadata = metadata)
  
  bb$pct_reads_in_peaks <- bb$atac_peak_region_fragments / bb$atac_fragments * 100
  bb$pct_usable_fragments = bb$atac_fragments/bb$atac_raw_reads
  
  DefaultAssay(bb) <- "ATAC"
  
  bb <- NucleosomeSignal(bb)
  bb <- TSSEnrichment(bb, fast = FALSE)
  
  Idents(bb) = bb$condition
  bb$high.tss <- ifelse(bb$TSS.enrichment > 5, 'High', 'Low')
  
  pdfname = paste0(resDir, '/scATAC_QCs_cellRangerOutput_', design$condition[n],  version.analysis,  '.pdf')
  pdf(pdfname, width=12, height = 8)
  
  TSSPlot(bb, group.by = 'high.tss') + NoLegend()
  
  VlnPlot(bb, features = "nCount_ATAC", ncol = 1, y.max = 5000) +
  geom_hline(yintercept = c(500, 1000))
  
  VlnPlot(bb, 
          features = c("nCount_ATAC", "atac_fragments", "pct_reads_in_peaks"), 
          ncol = 3)
  
  cat(median(bb$nCount_ATAC), '--', median(bb$atac_fragments), ' -- ', median(bb$pct_reads_in_peaks), ' \n')
  
  VlnPlot(
    object = bb,
    features = c("TSS.enrichment", "nucleosome_signal"),
    ncol = 2,
    pt.size = 0
  )
  
  # quick filtering 
  bb <- subset(
    x = bb,
    subset = nCount_ATAC > 200 &
      nCount_ATAC < 20000 &
      #pct_reads_in_peaks > 15 &
      #blacklist_ratio < 0.05 &
      nucleosome_signal < 4 &
      TSS.enrichment > 1
  )
  
  bb = RunTFIDF(bb)
  bb = FindTopFeatures(bb, min.cutoff = 5)
  bb = RunSVD(bb)
  
  #DepthCor(bb, n = 30)
  #cordat = DepthCor(bb, n = 30)$data
  dims_use = c(2:30)
  print(dims_use)
  
  bb = FindNeighbors(object = bb, reduction = 'lsi', dims = dims_use, 
                                   force.recalc = T, graph.name = "thegraph")
  bb = FindClusters(object = bb, verbose = FALSE, algorithm = 3, 
                                  graph.name = "thegraph", resolution = 1)
  bb <- RunUMAP(object = bb, reduction = 'lsi', dims = dims_use)
  DimPlot(object = bb, label = TRUE) + NoLegend()
  
  dev.off()
  
  srat[[n]] = bb
  rm(bb)
  
}


########################################################
########################################################
# Section : recall peaks using macs2 
# 1) recall peaks for eahc sample and save it
# 2) peak consensus and recount 
########################################################
########################################################
srat_macs = list()
for(n in 1:nrow(design))
{
  table(bb$thegraph_res.1)
  
  macs = CallPeaks(object = bb, 
                    group.by = "thegraph_res.1",
                    #macs2.path = "/software/2020/software/macs2/2.2.5-foss-2018b-python-3.6.6/bin/macs2",
                    macs2.path = "/groups/tanaka/People/current/jiwang/local/anaconda3/envs/macs3/bin/macs3",
                    effective.genome.size = 2.0e+10, verbose = TRUE)
  
  saveRDS(macs, file = paste0(RdataDir, 'macs2_peaks_', design$condition[n], '.rds'))
  
  combined_macs = macs
  
  peakwidths = width(combined_macs)
  combined_macs = combined_macs[peakwidths > 12]
  
  # there is a problem with coordinates starting at 0 for some reason...
  combined_macs = restrict(combined_macs, start = 1)
  combined_macs
  
  
  feat = FeatureMatrix(fragments = frags_l[[s]], features = combined_macs, 
                       cells = colnames(counts_l[[s]]$Peaks))
  ass = CreateChromatinAssay(feat, fragments = frags_l[[s]], sep = c(":", "-"), min.cells = 0,
                             min.features = 100, annotation = granges_axolotl)
  srat_macs_l[[s]] = CreateSeuratObject(ass, assay = "ATAC")
  srat_macs_l[[s]]$dataset = s
  srat_macs_l[[s]]$animal = strsplit(s, "_")[[1]][1]
  
  
  
  
}


rm(aa)

save(design, scn, 
     file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, version.analysis, '.Rdata'))


