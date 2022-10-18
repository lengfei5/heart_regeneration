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

DimPlot(aa, reduction = '')

bb <- readRDS("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets_20221004_2.rds")

aa$subtypes_RNA -> Idents(aa)

FB_subset <- subset(aa,  subtypes_RNA %in% c("FB_1", "FB_2", "FB_3"))

FB_subset <- FindVariableFeatures(FB_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(FB_subset)
FB_subset <- ScaleData(FB_subset, features = all.genes)

FB_subset <- RunPCA(FB_subset, features = VariableFeatures(object = FB_subset))


FB_subset <- FindNeighbors(FB_subset, dims = 1:10)


FB_subset <- FindClusters(FB_subset, resolution = 0.1)


FB_subset <- RunUMAP(FB_subset, dims = 1:10)



DimPlot(FB_subset, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(FB_subset, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "subtypes")
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

saveRDS(aa , "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq/seuratObj_multiome_snRNA.annotated_scATAC.merged.peaks.cellranger.441K_v3.rds")


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
geneids = get_geneID(ggs)

for(n in 1:length(geneids))
{
  cat(n, '\n')
  annotation$gene_name[which(annotation$gene_id == geneids[n])] = ggs[n]
}

# tx_id required for the annotation plot in CoveragePlot
# https://github.com/stuart-lab/signac/issues/1159
#annotation$tx_id = annotation$transcript_id
#aa@assays$ATAC@annotation = annotation

annotation = readRDS("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq/Annotation_atac.rds")
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
#aa <- SCTransform(aa)
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

#saveRDS(aa , "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/snATACseq/seuratObj_multiome_snRNA.annotated_scATAC.merged.peaks.cellranger.441K_v3.rds")

aa <- LinkPeaks(
  object = aa,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  
  genes.use = c("CD68-AMEX60DD012740")
)


##########################################
# try to modify the ATAC annotation to have the same gene names in both RNA and ATAC
# it worked
##########################################
features = rownames(bb)[grep('TOP2A', rownames(bb))]
FeaturePlot(bb, features = features, order = TRUE, cols = c('gray', 'red'))

DefaultAssay(aa) <- "ATAC"
aa <- LinkPeaks(
  object = aa,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  
  genes.use = features
)

#bb$subtypes -> Idents(aa)

idents.plot <- c("Proliferating_CM", "Mono_Macrophages", "B_cells", "EC")


p1 <- CoveragePlot(
  object = aa,
  region = features,
  features = features,
  expression.assay = "SCT",
  idents = idents.plot ,
  extend.upstream = 500,
  extend.downstream = 10000
)

patchwork::wrap_plots(p1)

