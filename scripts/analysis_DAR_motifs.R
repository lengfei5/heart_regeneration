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
source('functions_scATAC.R')
library(ArchR)


library(JASPAR2020)
library(TFBSTools)
library(chromVAR)

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
require(tictoc)
library(future)
require(pheatmap)
require(RColorBrewer)
options(future.globals.maxSize = 80 * 1024^3)

mem_used()

###  import annotation and metadata
set.seed(1234)
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

DefaultAssay(srat_cr) <- 'RNA'
DimPlot(srat_cr, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()

DefaultAssay(srat_cr) <- 'ATAC'
DimPlot(srat_cr, label = TRUE, reduction = 'umap',
        group.by = 'celltypes',  repel = TRUE) + NoLegend()


# identify the DARs using the celltypes 
DefaultAssay(srat_cr) <- 'ATAC'

srat_cr = subset(srat_cr, cells = colnames(srat_cr)[which(srat_cr$celltypes != 'Neuronal')])
Idents(srat_cr) = srat_cr$celltypes

da_peaks <- FindAllMarkers(
  object = srat_cr,
  #ident.1 = cc,
  #ident.2 = ,
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'atac_peak_region_fragments'
)

# saveRDS(da_peaks, file = paste0(RdataDir, 'DAR_major.celltype_v2.rds'))
da_peaks = readRDS(file = paste0(RdataDir, 'DAR_major.celltype_v2.rds'))

cat(length(unique(da_peaks$gene)), 'unique peaks found \n')

#da_peaks = da_peaks[which(da_peaks$cluster != 'Neuronal'), ]
#subs = subset(srat_cr, cells = colnames(srat_cr)[which(srat_cr$celltypes != 'Neuronal')])

source('functions_scATAC.R')
library(ArchR)

peak.mat = aggregate_peak_signals_by_groups(srat_cr, group_by = 'celltypes', assay = 'ATAC') 
peak.mat = peak.mat[which(!is.na(match(rownames(peak.mat), da_peaks$gene))), ]

pheatmap(peak.mat, 
         cluster_rows=TRUE,
         #cutree_rows = 6,
         show_rownames=FALSE, 
         fontsize_row = 4,
         #color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(8), 
         color = ArchR::paletteContinuous(set = "solarExtra", n = 100),
         show_colnames = TRUE,
         scale = 'row',
         cluster_cols=TRUE, 
         #annotation_col=df,
         #gaps_col = gaps_col,
         legend = TRUE,
         treeheight_row = 15,
         annotation_legend = FALSE, 
         #annotation_colors = annot_colors,
         #clustering_callback = callback,
         #breaks = seq(-2, 2, length.out = 8),
         clustering_method = 'complete', 
         #cutree_rows = ,
         #breaks = seq(-range, range, length.out = 20),
         #gaps_row =  c(22, 79),
         legend_labels = FALSE,
         width = 4, height = 12, 
         filename = paste0(resDir, '/heatmap_DAR_v2.pdf'))

##########################################
# motif activities for major cell types
##########################################
source('functions_scATAC.R')
#motif_tf = readRDS(paste0(RdataDir, 'motif_to_tfs_pfm_JASPAR2020_CORE_vertebrate.rds'))
#motif_tf$name = paste0(motif_tf$tf, '_', motif_tf$motif)
motif_tf = readRDS(file = paste0(RdataDir, 'motif_to_tfs_pfm_JASPAR2020_CORE_vertebrate_v1.rds'))

aa = readRDS(file = paste0(RdataDir, 'atac_seuratObject_motifClass_chromVAR_v3.rds'))
aa = subset(aa, cells = colnames(srat_cr))

## motif enrichment analysis by groups
DefaultAssay(aa) <- 'ATAC'
Idents(aa) = aa$celltypes

groups = unique(aa$celltypes)

motif.mat = matrix(NA, ncol = length(groups), nrow = nrow(motif_tf))
colnames(motif.mat) = groups
rownames(motif.mat) = motif_tf$motif

for(n in 1:length(groups))
{
  # n = 1
  cat(n, '--', groups[n], '\n')
  da_peaks <- FindMarkers(
    object = aa,
    ident.1 = groups[n],
    ident.2 = NULL,
    only.pos = TRUE,
    test.use = 'LR',
    min.pct = 0.05,
    latent.vars = 'nCount_ATAC'
  )
  
  # get top differentially accessible peaks
  top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
  
  # test enrichment
  enriched.motifs <- FindMotifs(
    object = aa,
    features = top.da.peak
  )
  
  motif.mat[,n] = enriched.motifs$pvalue[match(rownames(motif.mat), enriched.motifs$motif)]
  
}

saveRDS(motif.mat, file = paste0(RdataDir, 'enriched_motif_pvalues_v1.rds'))

#motif.mat = aggregate_chromVar_scores_by_groups(aa, group_by = 'celltypes', assay = 'chromvar') 
#motif.mat = peak.mat[which(!is.na(match(rownames(peak.mat), da_peaks$gene))), ]
motif.mat = readRDS(file = paste0(RdataDir, 'enriched_motif_pvalues_v1.rds'))

motif.mat = -log10(motif.mat)

rownames(motif.mat) = motif_tf$name[match(rownames(motif.mat), motif_tf$motif)]

maxs = apply(motif.mat, 1, max)
o1 = order(-maxs)

motif.mat = motif.mat[o1, ]

#cat(length(which(maxs>10)), ' enrichment motifs \n')


ntop.per.group = 10
kk = c()
for(n in 1:ncol(motif.mat))
{
  # n = 1
  test = motif.mat[order(-motif.mat[,n]), ]
  test = test[1:ntop.per.group, ]
  kk = c(kk, match(rownames(test), rownames(motif.mat)))  
  
}
kk = unique(kk)

max.limit = 40
res = motif.mat[kk, ]
res[which(res>max.limit)] = max.limit

pheatmap(res, 
         cluster_rows=TRUE,
         #cutree_rows = 6,
         show_rownames=TRUE, 
         fontsize_row = 8,
         #color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(8), 
         #color = ArchR::paletteContinuous(set = "solarExtra", n = 100),
         color = ArchR::paletteContinuous(set = "comet", n = 100),
         show_colnames = TRUE,
         #scale = 'row',
         cluster_cols=TRUE, 
         #annotation_col=df,
         #gaps_col = gaps_col,
         legend = TRUE,
         treeheight_row = 15,
         annotation_legend = FALSE, 
         #annotation_colors = annot_colors,
         #clustering_callback = callback,
         #breaks = seq(-2, 2, length.out = 8),
         clustering_method = 'complete', 
         #cutree_rows = ,
         #breaks = seq(-range, range, length.out = 20),
         #gaps_row =  c(22, 79),
         legend_labels = FALSE,
         width = 5, height = 12, 
         filename = paste0(resDir, '/heatmap_enrichmentMotifs_v1.pdf'))


##########################################
# ChromVar analysis
##########################################
motif_tf = readRDS(paste0(RdataDir, 'motif_to_tfs_pfm_JASPAR2020_CORE_vertebrate_v1.rds'))

aa = readRDS(file = paste0(RdataDir, 'atac_seuratObject_motifClass_chromVAR_v3.rds'))
aa = subset(aa, cells = colnames(srat_cr))

DefaultAssay(aa) <- 'chromvar'

ss = colSums(aa@assays$chromvar@data)
length(which(is.na(ss)))

data = aa@assays$chromvar@data
data[which(is.na(data))] = 0

aa@assays$chromvar@data = data

#aa = subset(aa, cells = colnames(aa)[which(!is.na(ss))])

# look at the activity of individual motif

tf_sel = 'FEV'
feature = motif_tf$motif[which(motif_tf$tf == tf_sel)]

cat(tf_sel, ' -- ', feature, '\n')

FeaturePlot(
  object = aa,
  reduction = 'umap',
  #features = "MA0019.1",
  features = feature,
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1
) + ggtitle(paste0(tf_sel, ' -- ', feature))


DimPlot(aa, group.by = 'subtypes')


########################################################
########################################################
# Section : trajectory analysis to serach for the CM drivers 
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 'atac_seuratObject_motifClass_chromVAR_v3.rds'))
aa = subset(aa, cells = colnames(srat_cr))

##########################################
# test different methods for CM trajectory analysis 
##########################################
Test_trajectory_DM_EPgraph = FALSE
if(Test_trajectory_DM_EPgraph){
  aa = subset(aa, cells = colnames(aa)[which(aa$celltypes == 'CM')])
  aa$subtypes = droplevels(aa$subtypes)
  
  ref =  readRDS(file = paste0("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/", 
                              "CM_subset_for_velocity.rds"))
  
  mm = match(colnames(aa), colnames(ref))
  cells_CM = colnames(aa)[which(!is.na(mm))]
  mm = mm[which(!is.na(mm))]
  aa = subset(aa, cells = cells_CM)
  aa$subtypes = ref$subtypes[mm]
  
  CMmyLevels <- c("CM_ven_(Robo2)", "CM_ven_(Cav3_1)", "CM_IS","CM_Prol_IS", "CM_Prol_1", "CM_Prol_3")
  aa$subtypes = factor(aa$subtypes, levels = CMmyLevels)
  
  cols = c("#4CC9F0", "#49A2F0",
           "#941F56", "#C61010",  
           "#4361EE",  "#3F37C9")
  
  
  DimPlot(aa, dims = c(1, 2), label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE,
          cols = cols
  )
  
  #aa = subset(aa, cells = colnames(aa)[which(aa$subtypes == 'CM_IS'|aa$subtypes == "CM_Prol_IS")])
  
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
  # 
  aa <- ScaleData(aa)
  # 
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), weight.by.var = TRUE, verbose = FALSE)
  # 
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 20, min.dist = 0.3)
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)
  aa <- ScaleData(aa)
  
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), weight.by.var = TRUE, verbose = FALSE)
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.3)
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  
}





##########################################
# DAR of CM subtypes 
##########################################
Idents(srat_cr) = srat_cr$celltypes
sub.obj = subset(srat_cr, idents = 'CM')

Idents(sub.obj) = sub.obj$subtypes

da_peaks <- FindAllMarkers(
  object = sub.obj,
  #ident.1 = cc,
  #ident.2 = ,
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'atac_peak_region_fragments'
)

# saveRDS(da_peaks, file = paste0(RdataDir, 'DAR_subtypes_CM_v1.rds'))
da_peaks = readRDS(file = paste0(RdataDir, 'DAR_subtypes_CM_v1.rds'))

cat(length(unique(da_peaks$gene)), 'unique peaks found \n')

#da_peaks = da_peaks[which(da_peaks$cluster != 'Neuronal'), ]
#subs = subset(srat_cr, cells = colnames(srat_cr)[which(srat_cr$celltypes != 'Neuronal')])

source('functions_scATAC.R')
library(ArchR)

peak.mat = aggregate_peak_signals_by_groups(sub.obj, group_by = 'subtypes', assay = 'ATAC') 
peak.mat = peak.mat[which(!is.na(match(rownames(peak.mat), da_peaks$gene))), ]

pheatmap(peak.mat, 
         cluster_rows=TRUE,
         #cutree_rows = 6,
         show_rownames=FALSE, 
         fontsize_row = 4,
         #color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(8), 
         color = ArchR::paletteContinuous(set = "solarExtra", n = 100),
         show_colnames = TRUE,
         scale = 'row',
         cluster_cols=TRUE, 
         #annotation_col=df,
         #gaps_col = gaps_col,
         legend = TRUE,
         treeheight_row = 15,
         annotation_legend = FALSE, 
         #annotation_colors = annot_colors,
         #clustering_callback = callback,
         #breaks = seq(-2, 2, length.out = 8),
         clustering_method = 'complete', 
         #cutree_rows = ,
         #breaks = seq(-range, range, length.out = 20),
         #gaps_row =  c(22, 79),
         legend_labels = FALSE,
         width = 4, height = 12, 
         filename = paste0(resDir, '/heatmap_DAR_CM_subtypesv2.pdf'))




###########################################################################
###########################################################################
### not used from here 
###########################################################################
###########################################################################

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

