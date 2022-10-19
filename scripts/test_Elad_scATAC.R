##########################################################################
##########################################################################
# Project: Heart regeneration
# Purpose: snATAC-seq analysis for cardiac regeneration project
##########################################################################
##########################################################################

resDir = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/res_Elad'
#setwd("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq")

library(Signac)
library(Seurat)
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
library(BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M)
library(GenomicRanges)
library(future)
library(ballgown)

set.seed(1234)

#install.packages("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq/BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M_1.0.0.tar.gz", repos=NULL, type="source")

aa <- readRDS("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq/seuratObj_multiome_snRNA.annotated_scATAC.merged.peaks.cellranger.441K_v2.rds")
aa$subtypes_RNA -> Idents(aa)

##########################################
# modify the gene names in the ATAC-seq annotation to have the same as RNA assay
##########################################
# library(ballgown)
# gtf_axolotl = paste0("/groups/tanaka/People/current/jiwang/scripts/axolotl_multiome/r_package/", 
#                      "AmexT_v47.FULL_corr_chr_cut.gtf")
# 
# granges_axolotl = ballgown::gffReadGR(gtf_axolotl)
# # adding a gene biotype, as that's necessary for TSS metaprofile
# granges_axolotl$gene_biotype = "protein_coding"

annotation = aa@assays$ATAC@annotation

DefaultAssay(aa) <- "RNA"
ggs = rownames(aa)
geneids = get_geneID(ggs)

for(n in 1:length(geneids))
{
  cat(n, '\n') 
  annotation$gene_name[which(annotation$gene_id == geneids[n])] = ggs[n]
}



# tx_id required for the annotation plot in CoveragePlot
# https://github.com/stuart-lab/signac/issues/1159
annotation$tx_id = annotation$transcript_id 

saveRDS(annotation, file = paste0(RdataDir, 'modified_Amex47_annotation_multiome.rds'))

aa@assays$ATAC@annotation = annotation



## continue with analysis
VlnPlot(
  object = aa,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

aa <- subset(
  x = aa,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 200 &
    nCount_RNA > 1000 &
    nucleosome_signal < 6 &
    TSS.enrichment > 1
)


DefaultAssay(aa) <- "RNA"
aa <- SCTransform(aa)
aa <- RunPCA(aa)
aa <- FindNeighbors(aa, dims = 1:30)


aa <- FindClusters(aa, resolution = 0.5)

aa <- RunUMAP(aa, dims = 1:30)
DimPlot(aa, group.by = "subtypes_RNA")


DefaultAssay(aa) <- "ATAC"
aa <- FindTopFeatures(aa, min.cutoff = 5)
aa <- RunTFIDF(aa)
aa <- RunSVD(aa)

aa <- RegionStats(aa, genome = BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M)

aa <- LinkPeaks(
  object = aa,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  
  genes.use = c("PTPRC-AMEX60DD018338")
)

##########################################
# try to modify the ATAC annotation to have the same gene names in both RNA and ATAC
# it worked 
##########################################
DefaultAssay(aa) <- "ATAC"
aa <- LinkPeaks(
  object = aa,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  
  genes.use = c("PTPRC-AMEX60DD018338")
)


idents.plot <- c("Proliferating_CM", "Mono_Macrophages")

CoveragePlot(
  object = aa,
  region = "PTPRC-AMEX60DD018338",
  features = "PTPRC-AMEX60DD018338",
  expression.assay = "RNA",
  idents = idents.plot,
  extend.upstream = 500,
  extend.downstream = 10000,
  annotation = 'gene'
)

##########################################
# test Elad's subclustering of FB 
##########################################
DefaultAssay(aa) = 'RNA'

FB_subset <- subset(aa,  subtypes_RNA %in% c("FB_1", "FB_2", "FB_3"))

FB_subset <- FindVariableFeatures(FB_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(FB_subset)
FB_subset <- ScaleData(FB_subset, features = all.genes)

FB_subset <- RunPCA(FB_subset, features = VariableFeatures(object = FB_subset))


FB_subset <- FindNeighbors(FB_subset, dims = 1:10)


FB_subset <- FindClusters(FB_subset, resolution = 0.1)


FB_subset <- RunUMAP(FB_subset, dims = 1:10)


FB.cluster.ids <- c("FB_1","FB_2","Doublet","FB_3","FB_4")
names(FB.cluster.ids) <- levels(FB_subset)
FB_subset <- RenameIdents(FB_subset, FB.cluster.ids)


FB_subset$subtypes = Idents(FB_subset)

cell.sels = colnames(FB_subset)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes_RNA = as.character(aa$subtypes_RNA)
aa$subtypes_RNA[mm] <- as.character(FB_subset$subtypes)
aa$subtypes_RNA = as.factor(aa$subtypes_RNA)


########################################################
########################################################
# Section : check the scATAC-seq with the code of Elad
# 
########################################################
########################################################
resDir = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/res_Elad'
#setwd("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq")

library(Signac)
library(Seurat)
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
library(BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M)
library(GenomicRanges)
library(future)
library(ballgown)

set.seed(1234)

#install.packages("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq/BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M_1.0.0.tar.gz", repos=NULL, type="source")

aa <- readRDS("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq/seuratObj_multiome_snRNA.annotated_scATAC.merged.peaks.cellranger.441K_v3.rds")

DefaultAssay(aa) = 'RNA'
aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000)
all.genes <- rownames(aa)
aa <- ScaleData(aa, features = all.genes)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)

ElbowPlot(aa, ndims = 50)
aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)
DimPlot(aa, label = TRUE, repel = TRUE) + NoLegend()

# bb <- readRDS("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets_20221004_2.rds")
Idents(aa) = aa$subtypes_RNA
#aa$subtypes_RNA -> Idents(aa)
FB_subset <- subset(aa,  subtypes_RNA %in% c("FB_1", "FB_2", "FB_3"))

FB_subset <- FindVariableFeatures(FB_subset, selection.method = "vst", nfeatures = 3000)

FB_subset <- ScaleData(FB_subset)

FB_subset <- RunPCA(FB_subset, features = VariableFeatures(object = FB_subset))
ElbowPlot(FB_subset, ndims = 50)

FB_subset <- FindNeighbors(FB_subset, dims = 1:20)

FB_subset <- FindClusters(FB_subset, resolution = 0.3)

FB_subset <- RunUMAP(FB_subset, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

DimPlot(FB_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(FB_subset, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "subtypes_RNA")
DimPlot(FB_subset, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "condition")
DimPlot(FB_subset, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "seurat_clusters")
DimPlot(FB_subset, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "newId")

features = rownames(FB_subset)[grep('PECAM', rownames(EC_subset))]
FeaturePlot(FB_subset, features = features, order = TRUE, cols = c('gray', 'red'))

VlnPlot(FB_subset, features = features)

VlnPlot(aa, features = features) + NoLegend()

DimPlot(aa, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(aa, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "subtypes")
#### 2 = doublets


FB.cluster.ids <- c("FB_1","FB_2","Doublet","FB_3","FB_4")
names(FB.cluster.ids) <- levels(FB_subset)
FB_subset <- RenameIdents(FB_subset, FB.cluster.ids)


FB_subset$subtypes_RNA = Idents(FB_subset)


cell.sels = colnames(FB_subset)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes_RNA = as.character(aa$subtypes_RNA)
aa$subtypes_RNA[mm] <- as.character(FB_subset$subtypes_RNA)
aa$subtypes_RNA = as.factor(aa$subtypes_RNA)

aa$subtypes_RNA <- Idents(aa)

#saveRDS(aa , "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq/seuratObj_multiome_snRNA.annotated_scATAC.merged.peaks.cellranger.441K_v3.rds")


##########################################
# modify the gene names in the ATAC-seq annotation to have the same as RNA assay
##########################################
# library(ballgown)
# gtf_axolotl = paste0("/groups/tanaka/People/current/jiwang/scripts/axolotl_multiome/r_package/",
#                      "AmexT_v47.FULL_corr_chr_cut.gtf")
#
# granges_axolotl = ballgown::gffReadGR(gtf_axolotl)
# # adding a gene biotype, as that's necessary for TSS metaprofile
# granges_axolotl$gene_biotype = "protein_coding"

annotation = aa@assays$ATAC@annotation

DefaultAssay(aa) <- "RNA"
ggs = rownames(aa)

########
#get_geneName = function(aa)
#{
#  return(sapply(aa, function(x) {test = unlist(strsplit(as.character(x), '-')); test = test[-length(test)];
#  paste0(test, collapse = '-')}))
#}

get_geneID = function(aa)
{
  return(sapply(aa, function(x) {test = unlist(strsplit(as.character(x), '-')); return(test[length(test)])}))
  
}

########
# geneids = get_geneID(ggs)
# 
# for(n in 1:length(geneids))
# {
#   cat(n, '\n')
#   annotation$gene_name[which(annotation$gene_id == geneids[n])] = ggs[n]
# }

# tx_id required for the annotation plot in CoveragePlot
# https://github.com/stuart-lab/signac/issues/1159
#annotation$tx_id = annotation$transcript_id
#aa@assays$ATAC@annotation = annotation

## instead of rerun the annotation modification, I load the one I run before
annotation = readRDS(paste0('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/results/', 
                            'sc_multiome_R13591_intron.exon.20220729/Rdata/modified_Amex47_annotation_multiome.rds'))
#annotation = readRDS("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq/Annotation_atac.rds")
#saveRDS(annotation , "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq/Annotation_atac.rds")
aa@assays$ATAC@annotation <- annotation

## continue with analysis
VlnPlot(
  object = aa,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

aa <- subset(
  x = aa,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 200 &
    nCount_RNA > 1000 &
    nucleosome_signal < 6 &
    TSS.enrichment > 1
)


DefaultAssay(aa) <- "RNA"
aa$subtypes = aa$subtypes_RNA
aa$celltypes = as.character(aa$subtypes)

aa$celltypes[grep('CM_|CMs_|_CM|_CM_', aa$subtypes)] = 'CM'
aa$celltypes[grep('EC_|_EC', aa$subtypes)] = 'EC'
aa$celltypes[grep('FB_', aa$subtypes)] = 'FB'
aa$celltypes[grep('B_cells', aa$subtypes)] = 'Bcell'

aa$celltypes[grep('Macrophages|_MF', aa$subtypes)] = 'Macrophages'
aa$celltypes[grep('Megakeryocytes', aa$subtypes)] = 'Megakeryocytes'
aa$celltypes[grep('RBC', aa$subtypes)] = 'RBC'


# saveRDS(aa, paste0("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq/",
#                     "seuratObj_multiome_snRNA.annotated_scATAC.merged.peaks.cellranger.441K_v3_testByJK.rds"))


DefaultAssay(aa) = 'RNA'
# renormalize the RNA data
aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 8000)
aa <- ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)

ElbowPlot(aa, ndims = 50)
aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1, 
              reduction.name = "umap")

DimPlot(aa, label = TRUE, repel = TRUE, reduction = 'umap') + NoLegend()

DimPlot(aa, label = TRUE, group.by = 'celltypes', repel = TRUE, reduction = 'umap') + NoLegend()

DimPlot(aa, label = TRUE, group.by = 'celltypes', repel = TRUE, reduction = 'umap_lsi') + NoLegend()

# normalize ATAC and UMAP
DefaultAssay(aa) <- "ATAC"
aa <- FindTopFeatures(aa, min.cutoff = 5)
aa <- RunTFIDF(aa)
aa <- RunSVD(aa)

DepthCor(aa)

aa <- RunUMAP(object = aa, reduction = 'lsi', dims = 2:50, n.neighbors = 30, min.dist = 0.1, 
              reduction.name = "umap_lsi")

DimPlot(object = aa, label = TRUE, group.by = 'celltypes', reduction = 'umap_lsi') + NoLegend()


DefaultAssay(aa) = 'RNA'

library(tictoc)
tic()
aa <- SCTransform(aa, ncells = 3000, assay = "RNA", verbose = TRUE, 
                    variable.features.n = 5000, 
                  return.only.var.genes = TRUE, vst.flavor = "v2")
toc()

aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 30)

aa <- RunUMAP(aa, dims = 1:40, n.neighbors = 50, min.dist = 0.3,  reduction.name = "umap_sct")
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', reduction = 'umap_sct') + NoLegend()

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'subtypes', reduction = 'umap_sct') + NoLegend()

saveRDS(aa, paste0("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq/",
                   "seuratObj_multiome_snRNA.annotated_scATAC.merged.peaks.cellranger.441K_v4_testByJK_", 
                   "logNormal_SCT_umap_lsiUmap_sctUmap.rds"))

##########################################
# link peaks and coveragePlots
#
##########################################
DefaultAssay(aa) = 'ATAC'
aa <- RegionStats(aa, genome = BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M)

# aa <- LinkPeaks(
#   object = aa,
#   peak.assay = "ATAC",
#   expression.assay = "RNA",
#   
#   genes.use = c("CD68-AMEX60DD012740")
#   
# )

features = rownames(aa@assays$RNA)[grep('ITGAM', rownames(aa@assays$RNA))]

DefaultAssay(aa) <- "RNA"
FeaturePlot(aa, features = features[1], order = TRUE, cols = c('gray', 'red'))

features = features[1]

## SCT normalization doesn't seem to be as good as lognormal normalizaiton, so use the 'RNA' rather 'SCT'
DefaultAssay(aa) <- "ATAC"
aa <- LinkPeaks(
  object = aa,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  genes.use = features
)

#bb$subtypes -> Idents(aa)
idents.plot <- c("Proliferating_CM",'CM_Atria', 'CM_OFT', "Mono_Macrophages", "B_cells", "EC")

p1 <- CoveragePlot(
  object = aa,
  region = features,
  features = features,
  expression.assay = "RNA",
  idents = idents.plot ,
  extend.upstream = 500,
  extend.downstream = 10000
)


patchwork::wrap_plots(p1)

