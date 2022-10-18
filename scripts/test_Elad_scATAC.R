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

