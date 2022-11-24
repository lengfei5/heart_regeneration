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

version.analysis = '_R13591_atac_reseq_20221115'

resDir = paste0("../results/sc_multiome", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')


if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R14353_ax_snATAC_reseq'

source('functions_scATAC.R')
source('functions_scRNAseq.R')
source('functions_Visium.R')


library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(patchwork)
require(SeuratObject)

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

design = data.frame(sampleID = seq(197254, 197258), 
                    condition = c(paste0('Amex_scATAC_d', c(0, 1, 4, 7, 14))), stringsAsFactors = FALSE)
design$timepoint = gsub('Amex_scATAC_', '', design$condition) 

source('functions_scRNAseq.R')
library(future)


########################################################
########################################################
# Section I: merge all cellranger peaks as peak consensus
# and quantify the count tables 
# 
########################################################
########################################################
## import scRNA seq data as reference to select cells
scRNA_file = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_subtypes_final_20221117.rds'
refs = readRDS(file = scRNA_file)
table(refs$subtypes)

# merge the peak first 
for(n in 1:nrow(design))
{
  # n = 1
  cat('----------- : ', n, ':',  design$condition[n], '-------------\n')
  topdir = paste0(dataDir, '/multiome_', design$timepoint[n], '/outs')
  
  p = read.table(paste0(topdir, "/atac_peaks.bed"), 
                     col.names = c("chr", "start", "end"))
  p = makeGRangesFromDataFrame(p)
  cat('---', length(p), ' peaks \n')
  
  if(n == 1){
    peaks = p
  }else{
    peaks = union(peaks, p)
  }
}

length(peaks)
combined.peaks = peaks
rm(peaks)


##########################################
# combine peaks filtering
##########################################
peakwidths = width(combined.peaks)
combined.peaks = combined.peaks[peakwidths > 50]
# there is a problem with coordinates starting at 0 for some reason...
combined.peaks = restrict(combined.peaks, start = 1)
cat(length(combined.peaks), ' combined peaks \n')

srat_cr = list()

for(n in 1:nrow(design))
#for(n in 2:nrow(design))
{
  # n = 1
  cat('----------- : ', n, ':',  design$condition[n], '-------------\n')
  
  # load nf output and process
  topdir = paste0(dataDir, '/multiome_', design$timepoint[n], '/outs')
  #counts <- Read10X_h5(filename = paste0(topdir, "/filtered_feature_bc_matrix.h5"))
  counts <- Read10X_h5(filename = paste0(topdir, "/raw_feature_bc_matrix.h5"))
  fragpath <- paste0(topdir, "/atac_fragments.tsv.gz")
  
  bc_rna = colnames(refs)[which(refs$condition == paste0("Amex_scRNA_", design$timepoint[n]))]
  cells_rna = sapply(bc_rna, function(x) {test = unlist(strsplit(as.character(x), '-'))[1]; paste0(test, '-1')})
  cells_rna = as.character(cells_rna)
  
  cells_peak = colnames(counts$Peaks)
  cat(length(cells_peak), ' cells from atac \n')
  cat(length(cells_rna), ' cell from rna \n')
  sum(!is.na(match(cells_rna, cells_peak)))
  
  #frags_l = CreateFragmentObject(path = fragpath, cells = colnames(counts$Peaks))
  frags_l = CreateFragmentObject(path = fragpath, cells = cells_rna)
  
  # slow step takes 16 mins without parall computation
  tic()
  feat = FeatureMatrix(fragments = frags_l, features = combined.peaks, cells = cells_rna)
  saveRDS(feat, file = paste0(RdataDir, 'snATAC_FeatureMatrix_', design$condition[n], '.rds'))
  toc()
  
  # Do not change the atac-seq barcode names, because fragment file is still connected with the cell bc;
  # so it is easier to change the cell barcodes of RNA assay
  # colnames(feat) = bc_rna 
  
  chrom_assay <- CreateChromatinAssay(
    counts = feat,
    sep = c(":", "-"),
    fragments = frags_l,
    annotation = granges_axolotl,
    min.cells = 0,
    min.features = 0 
  )
  
  bb <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    
  )
  
  rm(chrom_assay)
  
  metadata <- read.csv(
    file = paste0(topdir, '/per_barcode_metrics.csv'),
    header = TRUE,
    row.names = 1
  )
  
  bb$cell_bc = colnames(bb)
  bb$condition = gsub('_scATAC', '', design$condition[n])
  
  mm = match(bb$cell_bc, metadata$gex_barcode)
  metadata = metadata[mm, c(1,2, 3, 21:30)]
  #rownames(metadata) = colnames(bb)
  
  bb = AddMetaData(bb, metadata = metadata)
  
  bb$pct_reads_in_peaks <- bb$atac_peak_region_fragments / bb$atac_fragments * 100
  bb$pct_usable_fragments = bb$atac_fragments/bb$atac_raw_reads
  
  DefaultAssay(bb) <- "ATAC"
  
  bb = NucleosomeSignal(bb)
  bb = TSSEnrichment(bb, fast = FALSE)
  
  srat_cr[[n]] = bb
  
}

saveRDS(srat_cr, file = (paste0(RdataDir, 'seuratObj_scATAC_beforeMerged.peaks.cellranger.584K_v1.rds')))

srat_cr = readRDS(file = paste0(RdataDir, 'seuratObj_scATAC_beforeMerged.peaks.cellranger.584K_v1.rds'))
srat_reduced = Reduce(merge, srat_cr)

saveRDS(srat_reduced, file = (paste0(RdataDir, 'seuratObj_scATAC_merged.peaks.cellranger.584K_v1.rds')))

##########################################
# filtering and normalization
##########################################
srat_cr = readRDS(file = paste0(RdataDir, 'seuratObj_scATAC_merged.peaks.cellranger.584K_v1.rds'))
design$condition = gsub('_scATAC', '', design$condition)

levels = design$condition
srat_cr$condition = factor(srat_cr$condition, levels = levels)
Idents(srat_cr) = srat_cr$condition

# import the snRNA-seq data to merge multiome
cat('snRNA-seq --', basename(scRNA_file), '\n')
refs = readRDS(file = scRNA_file)

table(refs$subtypes)

Merge_scATAC_snRNA = FALSE
if(Merge_scATAC_snRNA){
  
  xx <- RenameCells(
    refs,
    new.names = colnames(x = srat_cr[["ATAC"]])
  )
  
  srat_cr[['RNA']] = xx[['RNA']]
  
  # refs[['ATAC']] = srat_cr[['ATAC']]
  # refs = AddMetaData(refs, metadata = srat_cr@meta.data)
  # saveRDS(refs, file = (paste0(RdataDir, 'seuratObj_snRNA_annotated_scATAC_merged.peaks.cellranger.441K_v0.rds')))
  # refs = readRDS(file = (paste0(RdataDir, 'seuratObj_snRNA_annotated_scATAC_merged.peaks.cellranger.441K_v0.rds')))
  
  metadata = refs@meta.data
  metadata = metadata[, grep('DF.classifications|pANN_', colnames(metadata), invert = TRUE)]
  metadata = metadata[, -1]
  colnames(metadata)[-c(1:2)] = paste0(colnames(metadata)[-c(1,2)], '_RNA')
  
  srat_cr = AddMetaData(srat_cr, metadata = metadata)
  saveRDS(srat_cr, file = paste0(RdataDir, 
                                 'seuratObj_multiome_snRNA.annotated_scATAC.merged.peaks.cellranger.441K_v2.rds'))
  
  #DefaultAssay(refs)  = 'ATAC'
  #refs = NucleosomeSignal(refs)
}

srat_cr <- NucleosomeSignal(srat_cr)
srat_cr <- TSSEnrichment(srat_cr, fast = FALSE)

srat_cr$high.tss <- ifelse(srat_cr$TSS.enrichment > 5, 'High', 'Low')

saveRDS(srat_cr, file = (paste0(RdataDir, 'seuratObj_scATAC_merged.peaks.cellranger_QCs.rds')))


TSSPlot(srat_cr, group.by = 'high.tss') + NoLegend()

VlnPlot(srat_cr, features = "nCount_ATAC", ncol = 1, y.max = 10000, group.by = 'condition', pt.size = 0., log = TRUE) +
  geom_hline(yintercept = c(1000, 1500, 2000))

ggsave(filename = paste0(resDir, '/QCs_nCount_ATAC_cellRangerPeaks.pdf'), height =8, width = 12)

VlnPlot(srat_cr, features = c("atac_fragments"), y.max = 10^6, group.by = 'condition', pt.size = 0, log = TRUE) +
  geom_hline(yintercept = c(10000, 20000))
ggsave(filename = paste0(resDir, '/QCs_atac_fragments_cellRangerPeaks_20k.pdf'), height =8, width = 12)

VlnPlot(srat_cr, features = c("atac_raw_reads"), y.max = 10^6, group.by = 'condition', pt.size = 0, log = TRUE) +
  geom_hline(yintercept = c(10000, 50000))
ggsave(filename = paste0(resDir, '/QCs_atac_rawReads_cellRangerPeaks_50k.pdf'), height =8, width = 12)


VlnPlot(srat_cr, features = c("pct_reads_in_peaks"))
ggsave(filename = paste0(resDir, '/QCs_pct_readsWithinPeaks_cellRangerPeaks.pdf'), height =8, width = 12 )


VlnPlot(object = srat_cr, features = c("TSS.enrichment"), pt.size = 0, y.max = 10) +
  geom_hline(yintercept = c(1, 2))
ggsave(filename = paste0(resDir, '/QCs_TSS.enrichment_cellRangerPeaks.pdf'), height =8, width = 12 )

VlnPlot(object = srat_cr, features = c("nucleosome_signal"), pt.size = 0)
ggsave(filename = paste0(resDir, '/QCs_nucleosome_signal_cellRangerPeaks.pdf'), height =8, width = 12 )

# quick filtering 
srat_cr <- subset(
  x = srat_cr,
  subset = nCount_ATAC > 100 &
    nCount_ATAC < 20000 &
    #pct_reads_in_peaks > 15 &
    #blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 1
)

srat_cr = RunTFIDF(srat_cr)
srat_cr = FindTopFeatures(srat_cr, min.cutoff = 'q5')
srat_cr = RunSVD(srat_cr)

DepthCor(srat_cr, n = 30)

#cordat = DepthCor(srat_cr, n = 30)$data
cordat = DepthCor(srat_cr, reduction = "lsi", n = 30)$data
dims_use = cordat$Component[abs(cordat$counts)<0.3]

dims_use = c(2:30)
print(dims_use)

srat_cr = FindNeighbors(object = srat_cr, reduction = 'lsi', dims = dims_use, 
                     force.recalc = T, graph.name = "thegraph")
srat_cr = FindClusters(object = srat_cr, verbose = FALSE, algorithm = 3, 
                    graph.name = "thegraph", resolution = 1)

srat_cr <- RunUMAP(object = srat_cr, reduction = 'lsi', dims = dims_use, n.neighbors = 30, min.dist = 0.1)

DimPlot(object = srat_cr, label = TRUE, repel = TRUE) + NoLegend()
ggsave(filename = paste0(resDir, '/cellRangerPeaks_umap_v1.pdf'), height =8, width = 12 )

DimPlot(object = srat_cr, label = TRUE, repel = TRUE, split.by = 'condition') + NoLegend()

ggsave(filename = paste0(resDir, '/cellRangerPeaks_umap_perCondition_v1.pdf'), height =8, width = 30 )


