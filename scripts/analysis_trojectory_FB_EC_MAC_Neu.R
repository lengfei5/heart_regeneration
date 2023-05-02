##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Apr 27 10:55:49 2023
##########################################################################
##########################################################################
rm(list = ls())
version.analysis = '_R13591_atac_reseq_20221115'

resDir = paste0("../results/sc_multiome", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

outDir = paste0('../results/sc_multiome_R13591_intron.exon.20220729', 
                '/FB_EC_Mac_Neu_trajectory_test/')
system(paste0('mkdir -p ', outDir))

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


##########################################
# import the scRNA and scATAC data
##########################################
aa = readRDS(file = paste0(RdataDir,
                           'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                           'scATAC.merged.peaks.cr_filtered_umap.lsi',
                           '584K.features_37680cells_umap.topics_updated.umap.subtypes_celltypes.rds'))

aa$time = gsub('Amex_', '', aa$condition)
aa$cell.ids = sapply(colnames(aa), function(x) unlist(strsplit(as.character(x), '-'))[1]) 
aa$cell.ids = paste0(aa$cell.ids, '_', aa$time)


# identify the DARs using the celltypes 
DefaultAssay(aa) <- 'ATAC'

aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != 'Neuronal')])
Idents(aa) = aa$celltypes

motif_tf = readRDS(file = paste0(RdataDir, 'motif_to_tfs_pfm_JASPAR2020_CORE_vertebrate_v1.rds'))
chromvar = readRDS(file = paste0(RdataDir, 'atac_seuratObject_motifClass_chromVAR_v3.rds'))
DefaultAssay(chromvar) <- 'chromvar'
ss = colSums(chromvar@assays$chromvar@data)
length(which(is.na(ss)))
data = chromvar@assays$chromvar@data
data[which(is.na(data))] = 0
chromvar@assays$chromvar@data = data

########################################################
########################################################
# Section I : EC 
# double check the subtype annotations 
# Prepare trajectory analysis with PAGA and RNA velocity
# run Diffusion map and 3D visualization
########################################################
########################################################
celltype_sel = 'EC'

sub_obj = subset(aa, cells = colnames(aa)[which(aa$celltypes == celltype_sel)])
sub_obj$subtypes = droplevels(sub_obj$subtypes)

## Endothecial cells are many and look pretty clear for the trajectory
ref =  readRDS(file = paste0("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/", 
                             "EC_subset_final_20221117.rds"))

mm = match(colnames(sub_obj), colnames(ref))
cat(length(which(is.na(mm))), ' cells without updated subtypes \n')
cat(length(which(is.na(match(colnames(ref), colnames(sub_obj))))), 
    ' cells filtered in the scATAC analysis\n')

cells_keep = colnames(sub_obj)[which(!is.na(mm))]
mm = mm[which(!is.na(mm))]

sub_obj = subset(sub_obj, cells = cells_keep)
sub_obj$subtypes = ref$subtypes[mm]

sub_obj$subtypes = droplevels(sub_obj$subtypes)

ECmyLevels <- c("EC", "EC_(NOS3)","EC_(WNT4)","EC_(LHX6)",  
                "EC_(CEMIP)","EC_Prol", "EC_IS_(LOX)","EC_IS_(IARS1)","EC_IS_Prol")

#factor(Idents(sub_obj, levels= ECmyLevels))
Idents(sub_obj) <- factor(sub_obj$subtypes, levels= ECmyLevels)

DimPlot(ref, group.by = 'subtypes')

umap.embedding = ref@reductions$umap@cell.embeddings
umap.embedding = umap.embedding[match(colnames(sub_obj), rownames(umap.embedding)), ]

sub_obj[['umap']] = Seurat::CreateDimReducObject(embeddings=umap.embedding,
                                                 key='UMAP_',
                                                 assay='RNA')
rm(umap.embedding)
DefaultAssay(sub_obj) <- 'RNA'

subtype_version = '_4subtypes'

refine_subpopulation_for_trajectory = TRUE
if(refine_subpopulation_for_trajectory){
  
  subtypes_sel = c("EC", "EC_IS_(IARS1)", "EC_IS_(LOX)", "EC_IS_Prol")
  sub_obj = subset(sub_obj,  cells = colnames(sub_obj)[which(!is.na(match(sub_obj$subtypes, subtypes_sel)))])
  
  sub_obj$subtypes = droplevels(sub_obj$subtypes)
  
  ## redo the clustering in case needed in the downstream analysis
  sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000)
  sub_obj <- ScaleData(sub_obj)
  sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), weight.by.var = TRUE, 
                    verbose = FALSE)
  ElbowPlot(sub_obj, ndims = 50)
  
  sub_obj <- RunUMAP(sub_obj, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
  
  DimPlot(sub_obj, group.by = 'subtypes', label = TRUE, repel = TRUE)
  
}else{
  ## redo the clustering in case needed in the downstream analysis
  sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000)
  sub_obj <- ScaleData(sub_obj)
  sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), weight.by.var = TRUE, 
                    verbose = FALSE)
  ElbowPlot(sub_obj, ndims = 50)
  
}

sub_obj <- FindNeighbors(sub_obj, dims = 1:20)
sub_obj <- FindClusters(sub_obj, verbose = FALSE, algorithm = 3, resolution = 0.3)

p1 = DimPlot(sub_obj, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()
p2 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)

p1 + p2

sub_obj$time = gsub('d', '', sub_obj$time)
sub_obj$clusters = sub_obj$seurat_clusters

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_subset_', celltype_sel, subtype_version,
                         '_reclustered.pdf'), 
       height = 5, width = 14)

DefaultAssay(sub_obj) <- 'RNA'
p1 = DimPlot(sub_obj, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()

DefaultAssay(sub_obj) <- 'ATAC'
p2 = DimPlot(sub_obj, label = TRUE, reduction = 'umap_topics',
             group.by = 'subtypes',  repel = TRUE) + NoLegend()
p1 + p2

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_subset_', celltype_sel,  subtype_version, '.pdf'), 
       height = 6, width = 14)

DefaultAssay(sub_obj) <- 'RNA'
DimPlot(sub_obj, label = TRUE, group.by = 'subtypes', split.by = 'condition', repel = TRUE) + NoLegend()

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_subset_', celltype_sel, subtype_version, 
                         '_bytimePoint.pdf'), 
       height = 5, width = 20)


##########################################
# preapre the spliced and unspliced matrix   
##########################################
DefaultAssay(sub_obj) = 'RNA'

saveRDS(sub_obj, file = paste0(outDir, 'seuratObj_multiome_snRNA_scATAC_',  celltype_sel, subtype_version, 
                               '.rds'))

source('utility_velocity.R')
# process_spliced_unspliced_kallisto(aa)
mnt = preapre_dataFile_for_RNAvelocity_PAGA(seuratObj = sub_obj)

Idents(mnt) = mnt$condition
mnt = subset(mnt, downsample = 2000)

saveFile = paste0('RNAmatrix_umap_kalisto.velocity_spliced_unspliced_',
                  'EC_subtypes_4subtypes_timepoints.all_downsample.10k.h5Seurat')

SaveH5Seurat(mnt, filename = paste0(outDir, saveFile), 
             overwrite = TRUE)
Convert(paste0(outDir, saveFile), 
        dest = "h5ad", overwrite = TRUE)


##########################################
# 3D DM to visualize the trajectory
##########################################
library(slingshot, quietly = FALSE)
library(destiny, quietly = TRUE)
library(mclust, quietly = TRUE)
library(scater)
library(SingleCellExperiment)
library(scran)
library(RColorBrewer)

DefaultAssay(sub_obj) = 'RNA'
sub_obj_diet = DietSeurat(sub_obj, 
                          counts = TRUE, 
                          data = TRUE,
                          scale.data = TRUE,
                          features = rownames(sub_obj), 
                          assays = c('RNA'), 
                          dimreducs = c('umap'), 
                          graphs = NULL, 
                          misc = TRUE)
sce = as.SingleCellExperiment(sub_obj_diet)

rm(sub_obj_diet)

dec <- modelGeneVar(sce)

#plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
#curve(metadata(dec)$trend(x), col="blue", add=TRUE)

nb_features = 3000; 
top.hvgs <- getTopHVGs(dec, n=nb_features)

sce <- runPCA(sce, subset_row=top.hvgs, ncomponents = 100)
# reducedDimNames(sce)
ll.pca = reducedDim(sce, 'PCA')[, c(1:50)]

n_neighbors = 100; n_eigs = 50; sigma = 'global';

tic()
dm <- DiffusionMap(ll.pca, sigma = sigma, k = n_neighbors, n_eigs = n_eigs, distance = 'euclidean')
toc()

cells = names(dm$DC1)
metadata = sub_obj@meta.data
dcs = data.frame(dm@eigenvectors, stringsAsFactors = FALSE)

dcs = dcs[match(rownames(metadata), cells), ]

dcs = as.matrix(dcs)

sub_obj[["DC"]] <- CreateDimReducObject(embeddings = as.matrix(dcs), key = "DC_", 
                                        assay = DefaultAssay(sub_obj))

rm(metadata)

#sub_obj = RunUMAP(sub_obj, reduction = "DC", dims = 1:30, n.neighbors = 100, min.dist = 0.3, 
      # metric = 'euclidean',
#             reduction.name = "dc_umap")
#DimPlot(sub_obj, reduction = 'dc_umap', label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE)

## try to make 3d for DC1, DC2 and DC3
## example from https://plotly-r.com/d-charts.html
library(plotly)
dcs = as.data.frame(dcs)
dcs$subtypes = sub_obj$subtypes[match(rownames(dcs), colnames(sub_obj))]
dcs$condition = sub_obj$condition[match(rownames(dcs), colnames(sub_obj))]


plot_ly(data.frame(dcs), x = ~DC1, y = ~DC2, z = ~DC3, size = 3) %>%
  add_markers(color = ~ subtypes)

plot_ly(data.frame(dcs), x = ~DC1, y = ~DC2, z = ~DC4, size = 3) %>%
  add_markers(color = ~ subtypes)

plot_ly(data.frame(dcs), x = ~DC1, y = ~DC3, z = ~DC4, size = 3) %>%
  add_markers(color = ~ subtypes)


p1 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
             group.by = 'subtypes')
p2 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 3), label = TRUE, repel = TRUE,
             group.by = 'subtypes')

p3 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 4), label = TRUE, repel = TRUE,
             group.by = 'subtypes')

p4 = DimPlot(sub_obj, reduction = 'DC', dims = c(3, 4), label = TRUE, repel = TRUE,
             group.by = 'subtypes')

(p1 / p2)|(p3 / p4)

ggsave(filename = paste0(outDir, 'trajectory_test_', celltype_sel, subtype_version, '.pdf'), 
       width = 12, height = 16)

p1 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
             group.by = 'condition')
p2 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 3), label = TRUE, repel = TRUE,
             group.by = 'condition')
p1/p2

DimPlot(sub_obj, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
        group.by = 'subtypes', split.by = 'condition')
ggsave(filename = paste0(outDir, 'trajectory_test_', celltype_sel, subtype_version, '_DC_split.by.time.pdf'), 
       width = 6, height = 18)

save(sub_obj, dcs, file = paste0(outDir, 'trajectory_test_DM_', 
                                 celltype_sel, subtype_version, 
                                 '.Rdata'))


########################################################
########################################################
# Section II : FB 
# double check the subtype annotations 
# Prepare trajectory analysis with PAGA and RNA velocity
# run Diffusion map and 3D visualization
########################################################
########################################################
celltype_sel = 'FB'

sub_obj = subset(aa, cells = colnames(aa)[which(aa$celltypes == celltype_sel)])
sub_obj$subtypes = droplevels(sub_obj$subtypes)

DimPlot(sub_obj, group.by = 'subtypes', label = TRUE, repel = FALSE)

DefaultAssay(sub_obj) <- 'RNA'

subtype_version = '_subtypes.3subtypes_timepoints.all'

refine_subpopulation_for_trajectory = FALSE
if(refine_subpopulation_for_trajectory){
  
  subtypes_sel = c("FB_(PKD1)", "FB_IS_(TFPI2)", "FB_IS_(TNC)")
  sub_obj = subset(sub_obj,  cells = colnames(sub_obj)[which(!is.na(match(sub_obj$subtypes, subtypes_sel)))])
  
  sub_obj$subtypes = droplevels(sub_obj$subtypes)
  
}

## redo the clustering in case needed in the downstream analysis
sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000)
sub_obj <- ScaleData(sub_obj)
sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), weight.by.var = TRUE, 
                  verbose = FALSE)
ElbowPlot(sub_obj, ndims = 50)

sub_obj <- RunUMAP(sub_obj, dims = 1:20, n.neighbors = 30, min.dist = 0.3)

DimPlot(sub_obj, group.by = 'subtypes', label = TRUE, repel = TRUE)

sub_obj <- FindNeighbors(sub_obj, dims = 1:20)
sub_obj <- FindClusters(sub_obj, verbose = FALSE, algorithm = 3, resolution = 0.3)

p1 = DimPlot(sub_obj, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()
p2 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)

p1 + p2

sub_obj$time = gsub('d', '', sub_obj$time)
sub_obj$clusters = sub_obj$seurat_clusters

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_subset_', celltype_sel, subtype_version,
                         '_reclustered.pdf'), 
       height = 5, width = 14)

DefaultAssay(sub_obj) <- 'RNA'
p1 = DimPlot(sub_obj, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()

DefaultAssay(sub_obj) <- 'ATAC'
p2 = DimPlot(sub_obj, label = TRUE, reduction = 'umap_topics',
             group.by = 'subtypes',  repel = TRUE) + NoLegend()
p1 + p2

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_subset_', celltype_sel,  subtype_version, '.pdf'), 
       height = 6, width = 14)

DefaultAssay(sub_obj) <- 'RNA'
DimPlot(sub_obj, label = TRUE, group.by = 'subtypes', split.by = 'condition', repel = TRUE) + NoLegend()

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_subset_', celltype_sel, subtype_version, 
                         '_bytimePoint.pdf'), 
       height = 5, width = 20)


##########################################
# preapre the spliced and unspliced matrix   
##########################################
DefaultAssay(sub_obj) = 'RNA'

saveRDS(sub_obj, file = paste0(outDir, 'seuratObj_multiome_snRNA_scATAC_',  celltype_sel, subtype_version, 
                               '.rds'))

source('utility_velocity.R')
# process_spliced_unspliced_kallisto(aa)
mnt = preapre_dataFile_for_RNAvelocity_PAGA(seuratObj = sub_obj)

Idents(mnt) = mnt$condition
#mnt = subset(mnt, downsample = 2000)

saveFile = paste0('RNAmatrix_umap_kalisto.velocity_spliced_unspliced_',
                  celltype_sel, subtype_version, '.h5Seurat')

SaveH5Seurat(mnt, filename = paste0(outDir, saveFile), 
             overwrite = TRUE)
Convert(paste0(outDir, saveFile), 
        dest = "h5ad", overwrite = TRUE)


##########################################
# 3D DM to visualize the trajectory
##########################################
library(slingshot, quietly = FALSE)
library(destiny, quietly = TRUE)
library(mclust, quietly = TRUE)
library(scater)
library(SingleCellExperiment)
library(scran)
library(RColorBrewer)

DefaultAssay(sub_obj) = 'RNA'
sub_obj_diet = DietSeurat(sub_obj, 
                          counts = TRUE, 
                          data = TRUE,
                          scale.data = TRUE,
                          features = rownames(sub_obj), 
                          assays = c('RNA'), 
                          dimreducs = c('umap'), 
                          graphs = NULL, 
                          misc = TRUE)
sce = as.SingleCellExperiment(sub_obj_diet)

rm(sub_obj_diet)

dec <- modelGeneVar(sce)

#plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
#curve(metadata(dec)$trend(x), col="blue", add=TRUE)

nb_features = 3000; 
top.hvgs <- getTopHVGs(dec, n=nb_features)

sce <- runPCA(sce, subset_row=top.hvgs, ncomponents = 100)
# reducedDimNames(sce)
ll.pca = reducedDim(sce, 'PCA')[, c(1:50)]

n_neighbors = 100; n_eigs = 50; sigma = 'global';

tic()
dm <- DiffusionMap(ll.pca, sigma = sigma, k = n_neighbors, n_eigs = n_eigs, distance = 'euclidean')
toc()

cells = names(dm$DC1)
metadata = sub_obj@meta.data
dcs = data.frame(dm@eigenvectors, stringsAsFactors = FALSE)

dcs = dcs[match(rownames(metadata), cells), ]

dcs = as.matrix(dcs)

sub_obj[["DC"]] <- CreateDimReducObject(embeddings = as.matrix(dcs), key = "DC_", 
                                        assay = DefaultAssay(sub_obj))

rm(metadata)

#sub_obj = RunUMAP(sub_obj, reduction = "DC", dims = 1:30, n.neighbors = 100, min.dist = 0.3, 
# metric = 'euclidean',
#             reduction.name = "dc_umap")
#DimPlot(sub_obj, reduction = 'dc_umap', label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE)

## try to make 3d for DC1, DC2 and DC3
## example from https://plotly-r.com/d-charts.html
library(plotly)
dcs = as.data.frame(dcs)
dcs$subtypes = sub_obj$subtypes[match(rownames(dcs), colnames(sub_obj))]
dcs$condition = sub_obj$condition[match(rownames(dcs), colnames(sub_obj))]


plot_ly(data.frame(dcs), x = ~DC1, y = ~DC2, z = ~DC3, size = 3) %>%
  add_markers(color = ~ subtypes)

plot_ly(data.frame(dcs), x = ~DC1, y = ~DC2, z = ~DC4, size = 3) %>%
  add_markers(color = ~ subtypes)

plot_ly(data.frame(dcs), x = ~DC1, y = ~DC3, z = ~DC4, size = 3) %>%
  add_markers(color = ~ subtypes)


p1 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
             group.by = 'subtypes')
p2 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 3), label = TRUE, repel = TRUE,
             group.by = 'subtypes')

p3 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 4), label = TRUE, repel = TRUE,
             group.by = 'subtypes')

p4 = DimPlot(sub_obj, reduction = 'DC', dims = c(3, 4), label = TRUE, repel = TRUE,
             group.by = 'subtypes')

(p1 / p2)|(p3 / p4)

ggsave(filename = paste0(outDir, 'trajectory_test_', celltype_sel, subtype_version, '.pdf'), 
       width = 12, height = 16)

p1 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
             group.by = 'condition')
p2 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 3), label = TRUE, repel = TRUE,
             group.by = 'condition')
p1/p2

DimPlot(sub_obj, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
        group.by = 'subtypes', split.by = 'condition')
ggsave(filename = paste0(outDir, 'trajectory_test_', celltype_sel, subtype_version, '_DC_split.by.time.pdf'), 
       width = 6, height = 18)

save(sub_obj, dcs, file = paste0(outDir, 'trajectory_test_DM_', 
                                 celltype_sel, subtype_version, 
                                 '.Rdata'))


########################################################
########################################################
# Section II : Macrophage and neutrophil
# double check the subtype annotations 
# Prepare trajectory analysis with PAGA and RNA velocity
# run Diffusion map and 3D visualization
########################################################
########################################################
celltype_sel = 'myeloid'

sub_obj = subset(aa, cells = colnames(aa)[which(aa$celltypes == 'Mo.Macs'|
                                                  aa$celltypes == 'Neu')])
sub_obj$subtypes = droplevels(sub_obj$subtypes)

DimPlot(sub_obj, group.by = 'subtypes', label = TRUE, repel = FALSE)

DefaultAssay(sub_obj) <- 'RNA'

subtype_version = '_subtypes.all_timepoints.all'

refine_subpopulation_for_trajectory = FALSE
if(refine_subpopulation_for_trajectory){
  
  subtypes_sel = c("FB_(PKD1)", "FB_IS_(TFPI2)", "FB_IS_(TNC)")
  sub_obj = subset(sub_obj,  cells = colnames(sub_obj)[which(!is.na(match(sub_obj$subtypes, subtypes_sel)))])
  
  sub_obj$subtypes = droplevels(sub_obj$subtypes)
  
}

## redo the clustering in case needed in the downstream analysis
sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000)
sub_obj <- ScaleData(sub_obj)
sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), weight.by.var = TRUE, 
                  verbose = FALSE)
ElbowPlot(sub_obj, ndims = 50)

sub_obj <- RunUMAP(sub_obj, dims = 1:20, n.neighbors = 30, min.dist = 0.1)

DimPlot(sub_obj, group.by = 'subtypes', label = TRUE, repel = TRUE)

sub_obj <- FindNeighbors(sub_obj, dims = 1:20)
sub_obj <- FindClusters(sub_obj, verbose = FALSE, algorithm = 3, resolution = 0.3)

p1 = DimPlot(sub_obj, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()
p2 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)

p1 + p2

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_subset_', celltype_sel, subtype_version,
                         '_reclustered.pdf'), 
       height = 5, width = 14)


sub_obj$time = gsub('d', '', sub_obj$time)
sub_obj$clusters = sub_obj$seurat_clusters

DefaultAssay(sub_obj) <- 'RNA'
p1 = DimPlot(sub_obj, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()

DefaultAssay(sub_obj) <- 'ATAC'
p2 = DimPlot(sub_obj, label = TRUE, reduction = 'umap_topics',
             group.by = 'subtypes',  repel = TRUE) + NoLegend()
p1 + p2

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_subset_', celltype_sel,  subtype_version, '.pdf'), 
       height = 6, width = 14)

DefaultAssay(sub_obj) <- 'RNA'
DimPlot(sub_obj, label = TRUE, group.by = 'subtypes', split.by = 'condition', repel = TRUE) + NoLegend()

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_subset_', celltype_sel, subtype_version, 
                         '_bytimePoint.pdf'), 
       height = 5, width = 20)


##########################################
# preapre the spliced and unspliced matrix   
##########################################
DefaultAssay(sub_obj) = 'RNA'

saveRDS(sub_obj, file = paste0(outDir, 'seuratObj_multiome_snRNA_scATAC_',  celltype_sel, subtype_version, 
                               '.rds'))

source('utility_velocity.R')
# process_spliced_unspliced_kallisto(aa)
mnt = preapre_dataFile_for_RNAvelocity_PAGA(seuratObj = sub_obj)

Idents(mnt) = mnt$condition
#mnt = subset(mnt, downsample = 2000)

saveFile = paste0('RNAmatrix_umap_kalisto.velocity_spliced_unspliced_',
                  celltype_sel, subtype_version, '.h5Seurat')

SaveH5Seurat(mnt, filename = paste0(outDir, saveFile), 
             overwrite = TRUE)
Convert(paste0(outDir, saveFile), 
        dest = "h5ad", overwrite = TRUE)


##########################################
# 3D DM to visualize the trajectory
##########################################
library(slingshot, quietly = FALSE)
library(destiny, quietly = TRUE)
library(mclust, quietly = TRUE)
library(scater)
library(SingleCellExperiment)
library(scran)
library(RColorBrewer)

DefaultAssay(sub_obj) = 'RNA'
sub_obj_diet = DietSeurat(sub_obj, 
                          counts = TRUE, 
                          data = TRUE,
                          scale.data = TRUE,
                          features = rownames(sub_obj), 
                          assays = c('RNA'), 
                          dimreducs = c('umap'), 
                          graphs = NULL, 
                          misc = TRUE)
sce = as.SingleCellExperiment(sub_obj_diet)

rm(sub_obj_diet)

dec <- modelGeneVar(sce)

#plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
#curve(metadata(dec)$trend(x), col="blue", add=TRUE)

nb_features = 3000; 
top.hvgs <- getTopHVGs(dec, n=nb_features)

sce <- runPCA(sce, subset_row=top.hvgs, ncomponents = 100)
# reducedDimNames(sce)
ll.pca = reducedDim(sce, 'PCA')[, c(1:50)]

n_neighbors = 100; n_eigs = 50; sigma = 'global';

tic()
dm <- DiffusionMap(ll.pca, sigma = sigma, k = n_neighbors, n_eigs = n_eigs, distance = 'euclidean')
toc()

cells = names(dm$DC1)
metadata = sub_obj@meta.data
dcs = data.frame(dm@eigenvectors, stringsAsFactors = FALSE)

dcs = dcs[match(rownames(metadata), cells), ]

dcs = as.matrix(dcs)

sub_obj[["DC"]] <- CreateDimReducObject(embeddings = as.matrix(dcs), key = "DC_", 
                                        assay = DefaultAssay(sub_obj))

rm(metadata)

#sub_obj = RunUMAP(sub_obj, reduction = "DC", dims = 1:30, n.neighbors = 100, min.dist = 0.3, 
# metric = 'euclidean',
#             reduction.name = "dc_umap")
#DimPlot(sub_obj, reduction = 'dc_umap', label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE)

## try to make 3d for DC1, DC2 and DC3
## example from https://plotly-r.com/d-charts.html
library(plotly)
dcs = as.data.frame(dcs)
dcs$subtypes = sub_obj$subtypes[match(rownames(dcs), colnames(sub_obj))]
dcs$condition = sub_obj$condition[match(rownames(dcs), colnames(sub_obj))]


plot_ly(data.frame(dcs), x = ~DC1, y = ~DC2, z = ~DC3, size = 3) %>%
  add_markers(color = ~ subtypes)

plot_ly(data.frame(dcs), x = ~DC1, y = ~DC2, z = ~DC4, size = 3) %>%
  add_markers(color = ~ subtypes)

plot_ly(data.frame(dcs), x = ~DC1, y = ~DC3, z = ~DC4, size = 3) %>%
  add_markers(color = ~ subtypes)


p1 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
             group.by = 'subtypes')
p2 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 3), label = TRUE, repel = TRUE,
             group.by = 'subtypes')

p3 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 4), label = TRUE, repel = TRUE,
             group.by = 'subtypes')

p4 = DimPlot(sub_obj, reduction = 'DC', dims = c(3, 4), label = TRUE, repel = TRUE,
             group.by = 'subtypes')

(p1 / p2)|(p3 / p4)

ggsave(filename = paste0(outDir, 'trajectory_test_', celltype_sel, subtype_version, '.pdf'), 
       width = 12, height = 16)

p1 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
             group.by = 'condition')
p2 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 3), label = TRUE, repel = TRUE,
             group.by = 'condition')
p1/p2

DimPlot(sub_obj, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
        group.by = 'subtypes', split.by = 'condition')
ggsave(filename = paste0(outDir, 'trajectory_test_', celltype_sel, subtype_version, '_DC_split.by.time.pdf'), 
       width = 6, height = 18)

save(sub_obj, dcs, file = paste0(outDir, 'trajectory_test_DM_', 
                                 celltype_sel, subtype_version, 
                                 '.Rdata'))


