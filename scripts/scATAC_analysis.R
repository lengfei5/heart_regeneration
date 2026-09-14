##########################################################################
##########################################################################
# Project: Aleks' single-cell MS project 
# Script purpose: analyze scATAC-seq with cisTopic
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Mar  2 14:43:04 2020
##########################################################################
##########################################################################
source.my.script <- function(name.of.function){
  tryCatch(path <- rstudioapi::getSourceEditorContext()$path,
           error = function(e){
             install.packages("rstudioapi")
             path <-  rstudioapi::getSourceEditorContext()$path})
  source.path <- sub(basename(path), "", path)
  source(paste0(source.path,name.of.function))
}

## set up the paths for the data and results
tryCatch(path <- rstudioapi::getSourceEditorContext()$path, 
         error = function(e){
           install.packages("rstudioapi")
           path <-  rstudioapi::getSourceEditorContext()$path})
source.path <- sub(basename(path), "", path)

user <- "results_jiwang/"
setwd(paste0("../", user))

version.DATA = 'scATAC_earlyEmbryo'
version.analysis =  paste0(version.DATA, '_20200302')
resDir = paste0("results/", version.analysis)
tabDir = paste0("results/", version.analysis, "/tables/")
RdataDir = paste0("results/", version.analysis, "/Rdata/")
if(!dir.exists("results/")){dir.create("results/")}
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

source.my.script('scATAC_functions.R')
#source.my.script('scATAC_cluster_annotation_functions.R')

########################################################
########################################################
# Section : peak annotatoin to have overview and QC for scATAC-seq peaks
# 
########################################################
########################################################
source.my.script('scATAC_functions.R')
filtered_mtx_dir = paste0("../output_cellranger.ce11_scATACpro/filtered_matrix_peaks_barcodes")

run_scATAC_peak_annotation(paste0(filtered_mtx_dir, '/peaks.bed'))

########################################################
########################################################
# Section : transformation comparison: lsi, lsi_log and LDA
# 
########################################################
########################################################
source.my.script('scATAC_functions.R')
tenx.bmat = load_tenx_atac(paste0(filtered_mtx_dir, '/matrix.mtx'), 
                                 paste0(filtered_mtx_dir, '/peaks.bed'), 
                                 paste0(filtered_mtx_dir, '/barcodes.tsv'))

# select features showing in > 50 cells
ss = Matrix::rowSums(tenx.bmat)
sum(ss>=50)

tenx.bmat = filter_features(tenx.bmat, cells=50)

# Binarize the matrix for consistency
tenx.bmat@x[tenx.bmat@x > 1] = 1

##########################################
# compare lsi and lsi_log
##########################################
tenx.seurat.lsi = lsi_workflow(tenx.bmat, dims=2:100, log_scale_tf=FALSE, reduction='pca', resolution=0.8)
tenx.seurat.lsi_log = lsi_workflow(tenx.bmat, dims=2:100, log_scale_tf=TRUE, reduction='pca.l2', resolution=0.8)

plot_clustering_comparison(tenx.seurat.lsi,
                           tenx.seurat.lsi_log,
                           reduction='umap',
                           description1='LSI',
                           description2='LSI logTF',
                           cluster_column1='peaks_snn_res.0.8',
                           cluster_column2='peaks_snn_res.0.8')

saveRDS(tenx.seurat.lsi, file =  paste0(RdataDir, 'atac_lsi_seurat_object.rds'))
saveRDS(tenx.seurat.lsi_log, file =  paste0(RdataDir, 'atac_lsi.log_seurat_object.rds'))

##########################################
# compare LAD from cistopic with lsi-log
##########################################
source.my.script('scATAC_functions.R')

Run.cistopic_workflow = TRUE
if(Run.cistopic_workflow){
  tic('run model for cisTopic')
  tenx.cistopic = cistopic_workflow(tenx.bmat, topic = seq(25, 120, by=5))
  saveRDS(tenx.cistopic, file =  paste0(RdataDir, 
                                        'atac_cisTopic_runWrapLDAModel_localRunning.rds'))
  toc()
  
}else{
  #tenx.cistopic = readRDS(file = paste0(RdataDir,
  #                                      'atac_cisTopic_runWrapLDAModel_localRunning.rds'))
  tenx.cistopic = readRDS(file =  paste0(RdataDir, 
                                         'atac_cisTopic_runWrapLDAModel_localRunning.rds'))
  
}

cisTopicObject = tenx.cistopic
par(mfrow=c(3,3))
cisTopicObject <- selectModel(cisTopicObject, type='maximum')
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
cisTopicObject <- selectModel(cisTopicObject, type='derivative')

seurat.cistopic = seurat_workflow_on_cistopic(tenx.cistopic, resolution=0.8, method='Z-score')

plot_clustering_comparison(tenx.seurat.lsi,
                            seurat.cistopic,
                            reduction='umap',
                            description1='LSI logTF',
                            description2='cisTopic',
                            cluster_column1='peaks_snn_res.1.5',
                            cluster_column2='peaks_snn_res.0.8')

##########################################
# select transformation method 
##########################################
resolution = 0.8
nb.pcs = 100; n.neighbors = 30; min.dist = 0.3;

tenx.seurat.lsi = FindNeighbors(object = tenx.seurat.lsi, 
                            reduction='pca', nn.eps=0.25, dims=2:nb.pcs)
tenx.seurat.lsi = FindClusters(object = tenx.seurat.lsi, 
                           n.start=20, resolution=resolution,
                           algorithm = 3)
tenx.seurat.lsi <- RunUMAP(object = tenx.seurat.lsi, 
                               reduction = 'pca', 
                               dims = 2:nb.pcs, n.neighbors = n.neighbors, 
                               min.dist = min.dist)
p1 = DimPlot(object = tenx.seurat.lsi_log, label = TRUE, reduction = 'umap') + 
  NoLegend()

tenx.seurat.lsi_log = FindNeighbors(object = tenx.seurat.lsi_log, 
                            reduction='pca.l2', nn.eps=0.25, dims=2:nb.pcs)
tenx.seurat.lsi_log = FindClusters(object = tenx.seurat.lsi_log, 
                           n.start=20, resolution=resolution,
                           algorithm = 3)

nb.pcs = 50; n.neighbors = 30; min.dist = 0.3;
tenx.seurat.lsi_log <- RunUMAP(object = tenx.seurat.lsi_log, 
                               reduction = 'pca.l2', 
                               dims = 2:nb.pcs, n.neighbors = n.neighbors, 
                               min.dist = min.dist)
p2 = DimPlot(object = tenx.seurat.lsi_log, label = TRUE, reduction = 'umap') + 
  NoLegend()

nb.pcs = ncol(seurat.cistopic[['pca']])
seurat.cistopic = FindNeighbors(object = seurat.cistopic,
                            reduction='pca', nn.eps=0.25, dims=1:nb.pcs)
seurat.cistopic = FindClusters(object = seurat.cistopic, 
                           n.start=20, resolution=0.8,
                           algorithm = 3)

n.neighbors = 30; min.dist = 0.3;
seurat.cistopic <- RunUMAP(object = seurat.cistopic, 
                           reduction = 'pca', 
                           dims = 1:nb.pcs, 
                           n.neighbors = n.neighbors, 
                           min.dist = min.dist)

Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
p1 = DimPlot(seurat.cistopic, label = TRUE, pt.size = 0.5, label.size = 8) + 
  NoLegend()

Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.1
p2 = DimPlot(seurat.cistopic, label = TRUE, pt.size = 0.5, label.size = 8) + 
  NoLegend()

p1 + p2
# nb.pcs = 35; n.neighbors = 30; min.dist = 0.25;
# seurat.cistopic <- RunUMAP(object = seurat.cistopic, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
# DimPlot(object = seurat.cistopic, label = TRUE, reduction = 'umap') + NoLegend()
Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
saveRDS(seurat.cistopic, file =  paste0(RdataDir, 'atac_LDA_seurat_object.rds'))


########################################################
########################################################
# Section : cluster annotations using cell-type specific features (motifs) or cis-regulatory elements 
# or using scRNA-seq data
########################################################
########################################################
#seurat.cistopic = load( file =  paste0(RdataDir, 'atac_LDA_seurat_object.rds'))
seurat.cistopic = readRDS(file =  paste0(RdataDir, 'atac_LDA_seurat_object.rds'))

DimPlot(seurat.cistopic, label = TRUE, pt.size = 0.4, label.size = 8) + NoLegend()

nb.pcs = ncol(seurat.cistopic[['pca']]); 
n.neighbors = 30; min.dist = 0.3;
seurat.cistopic <- RunUMAP(object = seurat.cistopic, reduction = 'pca', 
                           dims = 1:nb.pcs, 
                           n.neighbors = n.neighbors, min.dist = min.dist)
DimPlot(seurat.cistopic, label = TRUE, pt.size = 0.5, label.size = 8) + 
  NoLegend()

##########################################
# check the feature distribution for each clusters 
##########################################

VlnPlot(seurat.cistopic, features = 'nCount_peaks', log = TRUE, pt.size = 0.5) + 
  NoLegend()

VlnPlot(seurat.cistopic, features = 'nFeature_peaks', log = TRUE, pt.size = 0.3) + 
  NoLegend()


##########################################
# generate bigwigs for each cluster for visualization
##########################################
barcode.cluster = data.frame(Barcode = colnames(seurat.cistopic),
                             Cluster = seurat.cistopic$peaks_snn_res.0.8, 
                             stringsAsFactors = FALSE)
write.table(barcode.cluster, file = paste0(filtered_mtx_dir, '/barcodes_and_clusters.txt'), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')


##########################################
# Manually annotate cistopic-based clusters
# mainly based on two papers:
# Hashimshony et al. 2015 and Warner et al. 2019
##########################################
seurat.cistopic = readRDS(file =  paste0(RdataDir, 'atac_LDA_seurat_object.rds'))

DefaultAssay(seurat.cistopic) <- 'peaks'
Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
new.cluster.ids <- c("P0", "1.ABaxx", "2.MS", "3.AB", "4.ABa/p/MS/E", "5.E",
                     '6.E', '7.AB/EMS','8.ABp', '9.ABaxx', '10.MS', 
                     '11.ABp', '12.AB', '13.MS', '14.MS', '15.ABp',
                     '16.MS', '17.ABp', '18.AB','19.C/Ea', '20.ABa/p/MS/E', '21.P4')
names(new.cluster.ids) <- levels(seurat.cistopic)
seurat.cistopic <- RenameIdents(seurat.cistopic, new.cluster.ids)

DimPlot(seurat.cistopic, reduction = "umap", label = TRUE, pt.size = 2,  label.size = 10, repel = TRUE) + NoLegend()

########################################################
########################################################
# Section : Find differentially accessible peaks and motif enrichment analysis  
# 
########################################################
########################################################
DefaultAssay(seurat.cistopic) <- 'peaks'

seurat.cistopic = compute.motif.enrichment(seurat.cistopic)

DefaultAssay(seurat.cistopic) <- 'chromvar'

# tbx-38, med-2, elt-3, pal-1
motifs2check = rownames(seurat.cistopic)[grep('tbx-35|tbx-38|MA0542.1|MA0924.1 pal-1|MA0546.1|end-1|MA0547.1 skn-1|hnd-1|unc-130|mom-2|med-2', 
                                              rownames(seurat.cistopic))]
FeaturePlot(
  object = seurat.cistopic,
  features = motifs2check,
  pt.size = 0.1,
  max.cutoff = 'q90',
  min.cutoff = 'q10',
  ncol = 3
)

saveRDS(seurat.cistopic, file =  paste0(RdataDir, 'atac_LDA_seurat_object_motifClassChromVar.rds'))

##########################################
# compute gene activity score to annotate clusters
##########################################
Compute.gene.activity.scores = FALSE
if(Compute.gene.activity.scores){
  peaks.chrM = grep('chrM', rownames(seurat.cistopic))
  if(length(peaks.chrM)>0) seurat.cistopic = seurat.cistopic[-peaks.chrM, ]
  
  source.my.script('scATAC_functions.R')
  
  #fragment.file = '/Volumes/groups/cochella/jiwang/Projects/Aleks/R8898_scATAC/cellranger_atac_wbcel235/outs/fragments.tsv.gz'
  fragment.file = '../output_cellranger.ce11_scATACpro/cellranger_atac_ce11/fragments.tsv.gz'
  seurat.cistopic = compute.gene.acitivity.scores(seurat.cistopic, 
                                                  method = 'seurat',
                                                  fragment.file = fragment.file,
                                                  regions = 'promoter.geneBody',
                                                  size.promoter.upstream = 2000,
                                                  saveDir = RdataDir
                                                  )
  
  # normalize the gene activity matrix
  seurat.cistopic <- NormalizeData(
    object = seurat.cistopic,
    assay = 'RNA',
    normalization.method = 'CLR',
    scale.factor = median(seurat.cistopic$nCount_RNA)
  )

}

