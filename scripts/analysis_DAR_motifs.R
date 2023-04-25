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
                 paste0('/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/',
                        'aa_subtypes_final_20221117.rds'))
  
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

##########################################
# integrate motif and chromvar in the same object
# update the FINAL version of subtype annotation by Elad
##########################################
srat_cr = readRDS(file = paste0(RdataDir, 
                                'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                'scATAC.merged.peaks.cr_filtered_umap.lsi',
                                '584K.features_37680cells_umap.topics_updated.umap.subtypes.rds'))

DefaultAssay(srat_cr) <- 'RNA'
DimPlot(srat_cr, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()

DefaultAssay(srat_cr) <- 'ATAC'
DimPlot(srat_cr, label = TRUE, reduction = 'umap',
        group.by = 'celltypes',  repel = TRUE) + NoLegend()

motif_tf = readRDS(file = paste0(RdataDir, 'motif_to_tfs_pfm_JASPAR2020_CORE_vertebrate_v1.rds'))

chromvar = readRDS(file = paste0(RdataDir, 'atac_seuratObject_motifClass_chromVAR_v3.rds'))

aa = srat_cr
aa[['chromvar']] = chromvar[['chromvar']]
motifs = chromvar@assays$ATAC@motifs 
aa@assays$ATAC@motifs = motifs

DefaultAssay(aa) <- 'chromvar'

ss = colSums(aa@assays$chromvar@data)
length(which(is.na(ss)))

data = aa@assays$chromvar@data
data[which(is.na(data))] = 0

aa@assays$chromvar@data = data

saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                          'scATAC.merged.peaks.cr_filtered_umap.lsi',
                          '_motifs_chromvar.rds'))

rm(srat_cr)        
rm(chromvar)
rm(motifs)
rm(data)

aa = readRDS(file = paste0(RdataDir, 
                           'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                           'scATAC.merged.peaks.cr_filtered_umap.lsi',
                           '_motifs_chromvar.rds'))

refs = readRDS(file = paste0('/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/',
                             'aa_subtypes_final_20221117.rds'))

cells_sels = intersect(colnames(aa), colnames(refs))

aa = subset(aa, cells = cells_sels)

aa$subtypes = refs$subtypes[match(colnames(aa), colnames(refs))]

DimPlot(aa, group.by = 'subtypes', label = TRUE, repel = TRUE) + NoLegend()

saveRDS(aa, file = paste0(RdataDir, 
                          'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                          'scATAC.merged.peaks.cr_filtered_umap.lsi',
                          '_motifs_chromvar_subtypes.final.rds'))


### re-annoate the cell types
srat_cr = readRDS(file = paste0(RdataDir, 
                          'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                          'scATAC.merged.peaks.cr_filtered_umap.lsi',
                          '_motifs_chromvar_subtypes.final.rds'))

refs = srat_cr

refs$celltypes = as.character(refs$subtypes)

refs$celltypes[grep('CM_|CMs_|_CM|_CM_', refs$subtypes)] = 'CM'
refs$celltypes[grep('EC_|_EC', refs$subtypes)] = 'EC'
refs$celltypes[grep('FB_', refs$subtypes)] = 'FB'
refs$celltypes[grep('B_cells', refs$subtypes)] = 'B'
refs$celltypes[grep('T_cells', refs$subtypes)] = 'T'

refs$celltypes[grep('Macrophages|_MF', refs$subtypes)] = 'Macrophages'
refs$celltypes[grep('Megakeryocytes', refs$subtypes)] = 'Megakeryocytes'
refs$celltypes[grep('RBC', refs$subtypes)] = 'RBC'

refs$celltypes[grep('Mo/Macs', refs$subtypes)] = 'Mo.Macs'
refs$celltypes[grep('Neu_', refs$subtypes)] = 'Neu'

saveRDS(refs, file = paste0(RdataDir, 
                            'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                            'scATAC.merged.peaks.cr_filtered_umap.lsi',
                            '584K.features_37680cells_umap.topics_updated.umap.subtypes_celltypes.rds'))

#srat_cr = refs
#rm(refs)



########################################################
########################################################
# Section II : Differential Accessible Regions (DAR analysis)
# 
########################################################
########################################################
srat_cr = readRDS(file = paste0(RdataDir, 
                                'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                'scATAC.merged.peaks.cr_filtered_umap.lsi',
                                '584K.features_37680cells_umap.topics_updated.umap.subtypes_celltypes.rds'))

DefaultAssay(srat_cr) <- 'RNA'
DimPlot(srat_cr, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()

DimPlot(srat_cr, reduction = 'umap_topics',
        label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()

DefaultAssay(srat_cr) <- 'ATAC'
DimPlot(srat_cr, label = TRUE, reduction = 'umap',
        group.by = 'celltypes',  repel = TRUE) + NoLegend()

DimPlot(srat_cr, label = TRUE, reduction = 'umap_topics',
        group.by = 'celltypes',  repel = TRUE) + NoLegend()


# identify the DARs using the celltypes 
DefaultAssay(srat_cr) <- 'ATAC'

srat_cr = subset(srat_cr, cells = colnames(srat_cr)[which(srat_cr$celltypes != 'Neuronal')])
Idents(srat_cr) = srat_cr$celltypes


##########################################
# DA analysis either using FindAllMarkers or FindMarkers 
##########################################
for(pct_cutoff in  c(0.04, 0.03))
{
  pct_cutoff = 0.04
  cat('pct cutoff --- ', pct_cutoff, '---\n')
  da_peaks <- FindAllMarkers(
    object = srat_cr,
    #ident.1 = cc,
    #ident.2 = ,
    min.pct = pct_cutoff,
    test.use = 'LR',
    latent.vars = 'atac_peak_region_fragments'
    
  )
  saveRDS(da_peaks, file = paste0(RdataDir, 'DAR_major.celltype_FindAllMarkers_v4_pct.',
                                  pct_cutoff, '.rds'))
  
  ## too much time consuming 
  # clusters = unique(srat_cr$celltypes)
  # da_peaks = c()
  # for(n in 1:length(clusters))
  # {
  #   # n = 6; m = 7
  #   for(m in setdiff(1:length(clusters), n))
  #   {
  #     cat(clusters[n], ' vs. ', clusters[m], '\n')
  #     
  #     da_test <- FindMarkers(
  #       object = srat_cr,
  #       ident.1 = clusters[n],
  #       ident.2 = clusters[m],
  #       min.pct = pct_cutoff,
  #       test.use = 'LR', 
  #       latent.vars = 'atac_peak_region_fragments'
  #     )
  #     
  #     da_test$cluster = clusters[n]
  #     da_test$gene = rownames(da_test)
  #     da_peaks = rbind(da_peaks, da_test)
  #     
  #   }
    
  #}
  #saveRDS(da_peaks, file = paste0(RdataDir, 'DAR_major.celltype_FindMarkers_subtypes.one.vs.one_v4_pct.',
  #                                pct_cutoff, '.rds'))
  
}

#da_xx = readRDS(file = paste0(RdataDir, 'DAR_major.elltype_v3.rds'))
#da_peaks = readRDS(file = paste0(RdataDir, 'DAR_major.celltype_v2.rds'))

da_peaks = readRDS(file = paste0(RdataDir, 'DAR_major.celltype_FindAllMarkers_v4_pct.0.04.rds'))
cat(length(unique(da_peaks$gene)), 'unique peaks found \n')

#da_peaks = da_peaks[which(da_peaks$cluster != 'Neuronal'), ]
#subs = subset(srat_cr, cells = colnames(srat_cr)[which(srat_cr$celltypes != 'Neuronal')])

source('functions_scATAC.R')
library(ArchR)

peak.mat = aggregate_peak_signals_by_groups(srat_cr, group_by = 'celltypes', assay = 'ATAC') 
peak.mat = peak.mat[which(!is.na(match(rownames(peak.mat), da_peaks$gene))), ]
#groups = c("CM", "EC", "FB", "Mo.Macs", 'Neu', "Megakeryocytes", 'B', 'T', 'RBC')

# scale the matrix 
mat = binarySort_matrix(mat = peak.mat, limits = c(-2, 2), cutOff = 1, clusterCols = TRUE)
  
pheatmap(mat, 
         cluster_rows=FALSE,
         #cutree_rows = 6,
         show_rownames=TRUE, 
         fontsize_row = 10,
         #color = colorRampPalette(rev(brewer.pal(n = 8, name ="RdBu")))(8), 
         color = ArchR::paletteContinuous(set = "solarExtra", n = 100),
         show_colnames = FALSE,
         scale = 'none',
         cluster_cols=FALSE, 
         #annotation_col=df,
         #gaps_col = gaps_col,
         legend = TRUE,
         #treeheight_row = 15,
         annotation_legend = FALSE, 
         #annotation_colors = annot_colors,
         #clustering_callback = callback,
         #breaks = seq(-2, 2, length.out = 8),
         #clustering_method = 'complete', 
         #cutree_rows = ,
         #breaks = seq(-range, range, length.out = 20),
         #gaps_row =  c(22, 79),
         legend_labels = FALSE,
         width = 10, height = 4, 
         filename = paste0(resDir, '/heatmap_DAR_majoy.celltypes_v4.pdf'))

##########################################
# motif activities for major cell types
##########################################
source('functions_scATAC.R')
#motif_tf = readRDS(paste0(RdataDir, 'motif_to_tfs_pfm_JASPAR2020_CORE_vertebrate.rds'))
#motif_tf$name = paste0(motif_tf$tf, '_', motif_tf$motif)
motif_tf = readRDS(file = paste0(RdataDir, 'motif_to_tfs_pfm_JASPAR2020_CORE_vertebrate_v1.rds'))

#aa = readRDS(file = paste0(RdataDir, 'atac_seuratObject_motifClass_chromVAR_v3.rds'))
#aa = subset(aa, cells = colnames(srat_cr))

aa = srat_cr

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

saveRDS(motif.mat, file = paste0(RdataDir, 'enriched_motif_pvalues_v2.rds'))

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
# Section II: subsets of CM 
#  
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 
                           'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                           'scATAC.merged.peaks.cr_filtered_umap.lsi',
                           '_motifs_chromvar_subtypes.final.rds'))

# identify the DARs using the celltypes 
DefaultAssay(aa) <- 'ATAC'

aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != 'Neuronal')])
Idents(aa) = aa$celltypes


PLOT_chromVar_example = FALSE
if(PLOT_chromVar_example){
  DefaultAssay(aa) <- 'chromvar'
  motif_tf[grep('SMAD', motif_tf$tf), ]
  feature = motif_tf$motif[grep('SMAD', motif_tf$tf)]
  FeaturePlot(aa, features = feature, min.cutoff = 'q5', max.cutoff = 'q95')
  
}

##########################################
# the CM, machropage, FB
##########################################
for(celltype_sel in c('CM', 'EC', 'FB'))
{
  celltype_sel = 'CM'
  
  sub_obj =subset(aa, cells = colnames(aa)[which(aa$celltypes == celltype_sel)])
  sub_obj$subtypes = droplevels(sub_obj$subtypes)
  
  if(celltype_sel == 'CM'){
    ref =  readRDS(file = paste0("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/", 
                                 "CM_subset_for_velocity.rds"))
    
    mm = match(colnames(sub_obj), colnames(ref))
    cells_CM = colnames(sub_obj)[which(!is.na(mm))]
    mm = mm[which(!is.na(mm))]
    sub_obj = subset(sub_obj, cells = cells_CM)
    sub_obj$subtypes = ref$subtypes[mm]
    
    CMmyLevels <- c("CM_ven_(Robo2)", "CM_ven_(Cav3_1)", "CM_IS","CM_Prol_IS", "CM_Prol_1", "CM_Prol_3")
    sub_obj$subtypes = factor(sub_obj$subtypes, levels = CMmyLevels)
    
  }
  
  
  DefaultAssay(sub_obj) <- 'RNA'
  p1 = DimPlot(sub_obj, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()
  
  DefaultAssay(sub_obj) <- 'ATAC'
  p2 = DimPlot(sub_obj, label = TRUE, reduction = 'umap_topics',
          group.by = 'subtypes',  repel = TRUE) + NoLegend()
  p1 + p2
  
  # identify the DARs using the celltypes 
  DefaultAssay(sub_obj) <- 'ATAC'
  
  #sub_obj = subset(sub_obj, cells = colnames(sub_obj)[which(sub_obj$celltypes != 'Neuronal')])
  Idents(sub_obj) = sub_obj$subtypes
  
  da_peaks <- FindAllMarkers(
    object = sub_obj,
    #ident.1 = cc,
    #ident.2 = ,
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'atac_peak_region_fragments'
  )
  
  saveRDS(da_peaks, file = paste0(RdataDir, 'DAR_celltype_v2.rds'))
  da_peaks = readRDS(file = paste0(RdataDir, 'DAR_subtypes_CMsubet_v1.rds'))
  
  cat(length(unique(da_peaks$gene)), 'unique peaks found \n')
  
  #da_peaks = da_peaks[which(da_peaks$cluster != 'Neuronal'), ]
  #subs = subset(sub_obj, cells = colnames(sub_obj)[which(sub_obj$celltypes != 'Neuronal')])
  
  source('functions_scATAC.R')
  library(ArchR)
  
  peak.mat = aggregate_peak_signals_by_groups(sub_obj, group_by = 'subtypes', assay = 'ATAC') 
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
           filename = paste0(resDir, '/heatmap_DAR_CMsubset.pdf'))
  
  ##########################################
  # motif activities for major cell types
  ##########################################
  source('functions_scATAC.R')
  #motif_tf = readRDS(paste0(RdataDir, 'motif_to_tfs_pfm_JASPAR2020_CORE_vertebrate.rds'))
  #motif_tf$name = paste0(motif_tf$tf, '_', motif_tf$motif)
  motif_tf = readRDS(file = paste0(RdataDir, 'motif_to_tfs_pfm_JASPAR2020_CORE_vertebrate_v1.rds'))
  
  chromvar = readRDS(file = paste0(RdataDir, 'atac_seuratObject_motifClass_chromVAR_v3.rds'))
  chromvar = subset(chromvar, cells = colnames(sub_obj))
  
  chromvar$subtypes = sub_obj$subtypes[match(colnames(chromvar), colnames(sub_obj))]
  
  sub_obj = chromvar
  ## motif enrichment analysis by groups
  DefaultAssay(sub_obj) <- 'chromvar'
  Idents(sub_obj) = sub_obj$subtypes
  
  FeaturePlot(sub_obj, features = 'MA0002.2')
  
  groups = unique(as.character(sub_obj$subtypes))
  
  motif.mat = matrix(NA, ncol = length(groups), nrow = nrow(motif_tf))
  colnames(motif.mat) = groups
  rownames(motif.mat) = motif_tf$motif
  
  for(n in 1:length(groups))
  {
    # n = 1
    cat(n, '--', groups[n], '\n')
    da_peaks <- FindMarkers(
      object = sub_obj,
      ident.1 = groups[n],
      ident.2 = NULL,
      only.pos = TRUE,
      test.use = 'LR',
      min.pct = 0.05,
      latent.vars = 'nCount_ATAC'
    )
    
    # get top differentially accessible peaks
    top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.05, ])
    
    # test enrichment
    enriched.motifs <- FindMotifs(
      object = sub_obj,
      features = top.da.peak
    )
    
    motif.mat[,n] = enriched.motifs$pvalue[match(rownames(motif.mat), enriched.motifs$motif)]
    
  }
  
  saveRDS(motif.mat, file = paste0(RdataDir, 'enriched_motif_pvalues_CMsubset_v1.rds'))
  
  motif.mat = readRDS(file = paste0(RdataDir, 'enriched_motif_pvalues_CMsubset_v1.rds'))
  
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
  
  res = motif.mat[kk, ]
  
  max.limit = 40
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
           filename = paste0(resDir, '/heatmap_enrichmentMotifs_CMsubset_v1.pdf'))
  
  
}


