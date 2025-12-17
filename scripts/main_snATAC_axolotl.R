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
# axoltol genome and annotation (fragmented version) 
##########################################
# this (bespoke) package hosts the axolotl genome (comments from Tomas)
# package is very large, cannot be installed in home directory and should be in a fast disk
#install.packages("/groups/tanaka/People/current/jiwang/scripts/axolotl_multiome
# /r_package/BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M_1.0.0.tar.gz",
#                 repos = NULL, type = "source")
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


########################################################
########################################################
# Section I: merge all cellranger peaks as peak consensus
# and quantify the count tables 
# 
########################################################
########################################################
## import scRNA seq data as reference to select cells
scRNA_file = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets_20221004_2.rds'
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

########################################################
########################################################
# Section II : combine snRNA-seq and scATAC-seq
# - continue scATAC-seq analysis
########################################################
########################################################
srat_cr = readRDS(file = paste0(RdataDir, 'seuratObj_scATAC_merged.peaks.cellranger.584K_v1.rds'))
design$condition = gsub('_scATAC', '', design$condition)

levels = design$condition
srat_cr$condition = factor(srat_cr$condition, levels = levels)
Idents(srat_cr) = srat_cr$condition

##########################################
# merge snRNA and scATAC data 
##########################################
Merge_scATAC_snRNA = FALSE
if(Merge_scATAC_snRNA){
  
  # import the snRNA-seq data to merge multiome
  cat('snRNA-seq --', basename(scRNA_file), '\n')
  refs = readRDS(file = scRNA_file)
  
  table(refs$subtypes)
  
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
  saveRDS(srat_cr, 
          file = paste0(RdataDir, 
                        'seuratObj_multiome_snRNA.annotated_scATAC.merged.peaks.cellranger.584K_38280cells.rds'))
  
  rm(xx, refs)
  ##########################################
  # modify the gene names in the ATAC-seq annotation to have the same as RNA assay
  ##########################################
  Run.annotation.motification = FALSE
  if(Run.annotation.motification){
    annotation = srat_cr@assays$ATAC@annotation
    
    DefaultAssay(srat_cr) <- "RNA"
    ggs = rownames(srat_cr)
    
    get_geneID = function(srat_cr)
    {
      return(sapply(srat_cr, function(x) {test = unlist(strsplit(as.character(x), '-')); return(test[length(test)])}))
      
    }
    
    ######
    geneids = get_geneID(ggs)
    ggs = data.frame(ggs, geneids)
    
    mm = match(annotation$gene_id, ggs$geneids)
    ii = which(!is.na(mm))
    jj = mm[ii]
    annotation$gene_name[ii] = ggs$ggs[jj]
    
    # for(n in 1:length(geneids))
    # {
    #   cat(n, '\n')
    #   annotation$gene_name[which(annotation$gene_id == geneids[n])] = ggs[n]
    #   
    # }
    
    # tx_id required for the annotation plot in CoveragePlot
    # https://github.com/stuart-lab/signac/issues/1159
    annotation$tx_id = annotation$transcript_id
    
    # save the modified annotation
    #saveRDS(annotation , "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq/Annotation_atac.rds")
    
  }
  
  srat_cr@assays$ATAC@annotation <- annotation
  
  saveRDS(srat_cr, 
          file = paste0(RdataDir, 
                        'seuratObj_multiome_snRNA.annotated_scATAC.merged.peaks.cr.584K.annot_38280cells.rds'))
  
  
}

########################################################
########################################################
# Section III : analyze the scATAC and scRNA
########################################################
########################################################
srat_cr = readRDS(file = paste0(RdataDir, 
                               'seuratObj_multiome_snRNA.annotated_scATAC.merged.peaks.cr.',
                               '584K.annot_38280cells.rds'))

#srat_cr <- NucleosomeSignal(srat_cr)
#srat_cr <- TSSEnrichment(srat_cr, fast = FALSE)
srat_cr$high.tss <- ifelse(srat_cr$TSS.enrichment > 5, 'High', 'Low')
#saveRDS(srat_cr, file = (paste0(RdataDir, 'seuratObj_scATAC_merged.peaks.cellranger_QCs.rds')))


make.QC.plots = FALSE
if(make.QC.plots){
  #TSSPlot(srat_cr, group.by = 'high.tss') + NoLegend()
  
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
  
}

##########################################
# reprocess snRNÅ-seq data
##########################################
DefaultAssay(srat_cr) <- "RNA"
srat_cr$subtypes = srat_cr$subtypes_RNA
srat_cr$celltypes = as.character(srat_cr$subtypes)

srat_cr$celltypes[grep('CM_|CMs_|_CM|_CM_', srat_cr$subtypes)] = 'CM'
srat_cr$celltypes[grep('EC_|_EC', srat_cr$subtypes)] = 'EC'
srat_cr$celltypes[grep('FB_', srat_cr$subtypes)] = 'FB'
srat_cr$celltypes[grep('B_cells', srat_cr$subtypes)] = 'Bcell'

srat_cr$celltypes[grep('Macrophages|_MF', srat_cr$subtypes)] = 'Macrophages'
srat_cr$celltypes[grep('Megakeryocytes', srat_cr$subtypes)] = 'Megakeryocytes'
srat_cr$celltypes[grep('RBC', srat_cr$subtypes)] = 'RBC'

## need to rerun normalizaiton, PCA and UMAP
DefaultAssay(srat_cr) = 'RNA'
# renormalize the RNA data
srat_cr <- NormalizeData(srat_cr, normalization.method = "LogNormalize", scale.factor = 10000)

use.SCT.normalization = FALSE
if(use.SCT.normalization){
  library(tictoc)
  tic()
  srat_cr <- SCTransform(srat_cr, ncells = 3000, assay = "RNA", verbose = TRUE, 
                         variable.features.n = 5000, 
                         return.only.var.genes = TRUE, vst.flavor = "v2")
  toc()
  
  srat_cr <- RunPCA(srat_cr, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(srat_cr, ndims = 30)
  
  srat_cr <- RunUMAP(srat_cr, dims = 1:40, n.neighbors = 50, min.dist = 0.3,  reduction.name = "umap_sct")
  DimPlot(srat_cr, label = TRUE, repel = TRUE, group.by = 'celltypes', reduction = 'umap_sct') + NoLegend()
  
  DimPlot(srat_cr, label = TRUE, repel = TRUE, group.by = 'subtypes', reduction = 'umap_sct') + NoLegend()
  
  saveRDS(srat_cr, paste0("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq/",
                          "seuratObj_multiome_snRNA.annotated_scATAC.merged.peaks.cellranger.441K_v4_testByJK_", 
                          "logNormal_SCT_umap_lsiUmap_sctUmap.rds"))
}

srat_cr <- FindVariableFeatures(srat_cr, selection.method = "vst", nfeatures = 8000)
srat_cr <- ScaleData(srat_cr)
srat_cr <- RunPCA(srat_cr, features = VariableFeatures(object = srat_cr), verbose = FALSE)

ElbowPlot(srat_cr, ndims = 50)
srat_cr <- RunUMAP(srat_cr, dims = 1:30, n.neighbors = 30, min.dist = 0.1, 
                   reduction.name = "umap")

DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap') + NoLegend()

DimPlot(srat_cr, label = TRUE, group.by = 'celltypes', repel = TRUE, reduction = 'umap') + NoLegend()

saveRDS(srat_cr, file = paste0(RdataDir, 
              'seuratObj_multiome_snRNA.annotated.normalized.umap_',
              'scATAC.merged.peaks.cr.',
              '584K.annot_38280cells.rds'))

##########################################
# process and normalize the ATAC-seq data 
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

Test_cisTopic = FALSE
if(Test_cisTopic){
  
}


srat_cr <- RunTFIDF(srat_cr)
srat_cr = FindTopFeatures(srat_cr, min.cutoff = 'q5')
srat_cr <- RunSVD(srat_cr)

DepthCor(srat_cr, n = 30)

#cordat = DepthCor(srat_cr, reduction = "lsi", n = 30)$data
#dims_use = cordat$Component[abs(cordat$counts)<0.3]

srat_cr = FindNeighbors(object = srat_cr, reduction = 'lsi', dims = dims_use, 
                        force.recalc = T, graph.name = "thegraph")
srat_cr = FindClusters(object = srat_cr, verbose = FALSE, algorithm = 3, 
                       graph.name = "thegraph", resolution = 1)

dims_use = c(2:30)
print(dims_use)

srat_cr <- RunUMAP(object = srat_cr, reduction = 'lsi', dims = 2:30, n.neighbors = 30, min.dist = 0.1, 
                   reduction.name = "umap_lsi")


DimPlot(object = srat_cr, label = TRUE, reduction = 'umap_lsi') + NoLegend()

p1 = DimPlot(srat_cr, label = TRUE, repel = TRUE, reduction = 'umap', group.by = 'subtypes') + 
  NoLegend() + ggtitle('snRNA-seq')
p2 = DimPlot(object = srat_cr, label = TRUE, repel = TRUE, group.by = 'subtypes') + 
  NoLegend() + ggtitle('scATAC-seq')

p1 + p2
ggsave(filename = paste0(resDir, '/multiome_snRNA_snATAC_filtered.pdf'), height = 8, width = 20)

#DimPlot(object = srat_cr, label = TRUE, repel = TRUE, split.by = 'condition') + NoLegend()
#ggsave(filename = paste0(resDir, '/cellRangerPeaks_umap_perCondition_v1.pdf'), height =8, width = 30 )

FeaturePlot(srat_cr, features = c('nCount_ATAC', 'nucleosome_signal'), reduction = 'umap_lsi')


##########################################
# link peaks and coveragePlots
# 
##########################################
DefaultAssay(srat_cr) = 'ATAC'
srat_cr <- RegionStats(srat_cr, genome = BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M)

# srat_cr <- LinkPeaks(
#   object = srat_cr,
#   peak.assay = "ATAC",
#   expression.assay = "RNA",
#   
#   genes.use = c("CD68-AMEX60DD012740")
#   
# )

features = rownames(srat_cr@assays$RNA)[grep('ITGAM', rownames(srat_cr@assays$RNA))]

DefaultAssay(srat_cr) <- "RNA"
FeaturePlot(srat_cr, features = features[1], order = TRUE, cols = c('gray', 'red'))

features = features[1]

## SCT normalization doesn't seem to be as good as lognormal normalizaiton, so use the 'RNA' rather 'SCT'
DefaultAssay(srat_cr) <- "ATAC"
srat_cr <- LinkPeaks(
  object = srat_cr,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  genes.use = features
)

srat_cr$subtypes = srat_cr$subtypes_RNA
Idents(srat_cr) = as.factor(srat_cr$subtypes)
idents.plot <- c("Proliferating_CM",'CM_Atria', 'CM_OFT', "Mono_Macrophages", "B_cells", "EC")

p1 <- CoveragePlot(
  object = srat_cr,
  region = features,
  features = features,
  expression.assay = "RNA",
  idents = idents.plot ,
  extend.upstream = 500,
  extend.downstream = 10000
)


patchwork::wrap_plots(p1)

ggsave(filename = paste0(resDir, '/LinkPeaks_test.example.ITGAM.pdf'), height =8, width = 12)

