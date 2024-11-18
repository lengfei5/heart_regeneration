##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: process and analyze the neonatal mice scRNA-seq/snRNA-seq
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Nov  8 13:29:36 2023
##########################################################################
##########################################################################
rm(list = ls())

# required libraries
library(data.table)
require(Seurat)
library(SeuratObject)
require(sctransform)
require(ggplot2)
library(dplyr)
library(patchwork)
require(tictoc)
library(DropletUtils)
library(edgeR)
library(future)
library(pryr) # monitor the memory usage

source('functions_scRNAseq.R')
source('functions_Visium.R')

options(future.globals.maxSize = 160000 * 1024^2)
mem_used()

version.analysis = '_20231108'

resDir = paste0("../results/scRNAseq_neonatalMouse", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

Normalization = 'lognormal' # ('lognormal or SCT')
dataDir = "/groups/tanaka/Collaborations/Jingkui-Elad/Mouse_data_shoval/"


########################################################
########################################################
# Section I : First process the neonatal cardiomyocyte from 
# Cui et al., Dev Cell 2020
# still used here, because the CM dataset from Shoval was too stringently filtered 
########################################################
########################################################
Process.Cardiomyocyte.Cui.et.al.2020 = FALSE
if(Process.Cardiomyocyte.Cui.et.al.2020){
  dataDir = '../published_dataset/neonatal_mice/Cui_2020_DevCell/GSE130699_RAW'
  
  outDir = paste0(resDir, '/QCs')
  if(!dir.exists(outDir)) dir.create(outDir)
  
  matrix_list = list.files(path =  dataDir, 
                           pattern = '*_matrix.mtx.gz',
                           full.names = TRUE)
  design = data.frame(sampleID = basename(matrix_list),
                      condition = basename(matrix_list))
  
  design$condition = sapply(design$condition, function(x){unlist(strsplit(x, '_'))[2]})
  
  
  
  # import data from cellranger output
  for(n in 1:nrow(design))
  {
    # n = 1
    cat(n, ' : ', design$condition[n], '\n')
    
    exp = Matrix::readMM(paste0(dataDir, '/', design$sampleID[n])) #read matrix
    bc = read.csv(paste0(dataDir, "/", gsub('matrix.mtx.gz', 'barcodes.tsv.gz', design$sampleID[n])), 
                  header = F, stringsAsFactors = F)
    g = read.csv(paste0(dataDir, "/", gsub('matrix.mtx.gz', 'genes.tsv.gz', design$sampleID[n])), 
                 header = F, stringsAsFactors = F, sep = '\t')
    
    ## make unique gene names
    g$name = g$V2
    gg.counts = table(g$V2)
    gg.dup = names(gg.counts)[which(gg.counts>1)]
    index.dup = which(!is.na(match(g$V2, gg.dup)))
    g$name[index.dup] = paste0(g$V2[index.dup], '_', g$V1[index.dup])
    
    colnames(exp) = bc$V1
    rownames(exp) = g$name
    
    count.data = exp
    rm(exp);
    
    cat('get empty drops with UMI rank \n')
    
    # get emptyDrops and default cutoff cell estimates
    iscell_dd = defaultDrops(count.data, expected = 5000) # default cell estimate, similar to 10x cellranger
    sum(iscell_dd, na.rm=TRUE)
    
    ## not used the emptyDrops too slow 
    # eout = emptyDrops(count.data, lower = 200)
    # eout$FDR[is.na(eout$FDR)] = 1
    # iscell_ed = eout$FDR<=0.01
    # sum(iscell_ed, na.rm=TRUE)
    
    meta = data.frame(row.names = colnames(count.data), condition = design$condition[n],
                      iscell_dd = iscell_dd)
    
    # plot rankings for number of UMI
    br.out <- barcodeRanks(count.data)
    
    pdf(paste0(outDir, "/UMIrank_emptyDrop_", design$condition[n], "_", design$sampleID[n],  ".pdf"), 
        height = 6, width =10, useDingbats = FALSE)
    
    plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
    
    o <- order(br.out$rank)
    lines(br.out$rank[o], br.out$fitted[o], col="red")
    abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
    abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
    abline(v = sum(iscell_dd), col = 'darkgreen', lwd = 2.0)
    abline(v = c(3000, 5000, 8000, 10000, 12000), col = 'gray')
    text(x = c(3000, 5000, 8000, 10000, 12000), y =10000, labels = c(3000, 5000, 8000, 10000, 12000), 
         col = 'red')
    legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
           legend=c("knee", "inflection"))
    
    dev.off()
    
    # use defaultDrop to select cells.
    aa = CreateSeuratObject(counts = count.data[, iscell_dd],
                            meta.data = meta[iscell_dd, ], 
                            min.cells = 20, min.features = 100)
    aa$cell.id = paste0(colnames(aa), '_', design$condition[n], '_', design$sampleID[n])
    
    if(n == 1) {
      mnt = aa
    }else{
      mnt = merge(mnt, aa)
    }
  }
  
  
  mnt[["percent.mt"]] <- PercentageFeatureSet(mnt, pattern = "^mt-")
  
  save(design, mnt, 
      file = paste0(RdataDir, 'metadata_counts_Cui.et.al.2020.Rdata'))
  
  
  ##########################################
  # make Seurat object with metadata and counts  
  ##########################################
  load(file = paste0(RdataDir, 'metadata_counts_Cui.et.al.2020.Rdata'))
  aa = mnt
  rm(mnt)
  
  plot1 <- FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  
  Idents(aa) = factor(aa$condition)
  
  p1 = VlnPlot(aa, features = 'nFeature_RNA', y.max = 5000) +
    geom_hline(yintercept = c(200, 500, 1000)) + NoLegend()
  p2 = VlnPlot(aa, features = 'nCount_RNA', y.max = 50000) + NoLegend()
  p3 = VlnPlot(aa, features = 'percent.mt', y.max = 100) + NoLegend()
  
  p1 | p2 | p3
  
  aa = subset(aa, subset = nFeature_RNA > 200  & nCount_RNA < 20000 &  percent.mt < 25)
  
  if(SCT.normalization){
    aa <- SCTransform(aa, assay = "RNA", verbose = FALSE, variable.features.n = 3000, return.only.var.genes = FALSE, 
                      min_cells=5) 
    
  }else{
    aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
    aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
    
    plot1 <- VariableFeaturePlot(aa)
    
    top10 <- head(VariableFeatures(aa), 10)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    plot1 + plot2
    
    all.genes <- rownames(aa)
    aa <- ScaleData(aa, features = all.genes)
    
  }
  
  aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 30)
  
  aa <- FindNeighbors(aa, dims = 1:20)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
  
  if(Normalization == 'SCT'){
    aa = RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
    #saveRDS(aa, file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_SCT_umap.rds'))
  }else{
    
    aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 50, min.dist = 0.3)
    
    DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
    
    saveRDS(aa, file =  paste0(RdataDir, 'Seurat.obj_neonatalMice_Normalization_umap_CM_Cui2020.rds'))
    
  }
  
}

########################################################
########################################################
# Section II : process the data shared by Shoval for Wang et al., 2020 Cell Reports and 
# Cui et al., Dev Cell 2020 processed myself
########################################################
########################################################
double_check_two_seuratObjects_Wang2020 = FALSE
if(double_check_two_seuratObjects_Wang2020){
  xx = readRDS(file = paste0(dataDir, '/ZW_regeneration_scRNAseq.rds'))
  aa = readRDS(file = paste0(dataDir, '/Wang_DevCell_2020_scRNAseq.rds'))
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'orig.ident', raster=FALSE)
  p2 = DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'orig.ident', raster=FALSE)
  
  p1 + p2
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'FineID', raster=FALSE)
  p2 = DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'FineID', raster=FALSE)
  
  p1 + p2
  
  identical(aa, xx) # the same R object
  
  saveRDS(aa, file = paste0(RdataDir, 'Wang_DevCell_2020_scRNAseq.rds'))
  
}

##########################################
# futher processing of CM data 
##########################################
Further_cleaning_CM = FALSE
if(Further_cleaning_CM){
  aa = readRDS(file =  paste0(RdataDir, 'Seurat.obj_neonatalMice_Normalization_umap_CM_Cui2020.rds'))
  
  aa$condition = gsub('P11Sham', 'P8_Sham_d3', aa$condition)
  aa$condition = gsub('P1MID1', 'P1_MI_d1', aa$condition)
  aa$condition = gsub('P1MID3', 'P1_MI_d3', aa$condition)
  aa$condition = gsub('P2Sham', 'P1_Sham_d1', aa$condition)
  aa$condition = gsub('P4Sham', 'P1_Sham_d3', aa$condition)
  aa$condition = gsub('P8MID1', 'P8_MI_d1', aa$condition)
  aa$condition = gsub('P8MID3', 'P8_MI_d3', aa$condition)
  aa$condition = gsub('P9Sham', 'P8_Sham_d1', aa$condition)
  
  Idents(aa) = factor(aa$condition)
  p1 = VlnPlot(aa, features = 'nFeature_RNA', y.max = 5000) +
    geom_hline(yintercept = c(200, 500, 1000)) + NoLegend()
  p2 = VlnPlot(aa, features = 'nCount_RNA', y.max = 50000) + NoLegend()
  p3 = VlnPlot(aa, features = 'percent.mt', y.max = 100) + NoLegend()
  
  p1 | p2 | p3
  
  ## try to remove the doublets
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  sels = c()
  cc = unique(aa$condition)
  for(n in 1:length(cc))
  {
    jj = which(aa$condition == cc[n])
    sels = c(sels, jj[which(aa$nFeature_RNA[jj] < quantile(aa$nFeature_RNA[jj], 0.99) & 
                              aa$nFeature_RNA[jj] > quantile(aa$nFeature_RNA[jj], 0.01) & 
                              aa$nFeature_RNA[jj] > 500 & 
                           aa$nCount_RNA[jj] < quantile(aa$nCount_RNA[jj], 0.99) &
                             aa$percent.mt[jj] < 25)]) 
  }
  
  xx = subset(aa, cells = colnames(aa)[sels])
  p1 = VlnPlot(xx, features = 'nFeature_RNA', y.max = 5000) +
    geom_hline(yintercept = c(200, 500, 1000)) + NoLegend()
  p2 = VlnPlot(xx, features = 'nCount_RNA', y.max = 10000) + NoLegend()
  p3 = VlnPlot(xx, features = 'percent.mt', y.max = 30) + NoLegend()
  p1 | p2 | p3
  
  ggsave(filename = paste0(RdataDir, 'CM_Cui2020_futher_cleaning.pdf'), 
         width = 16, height = 6)
  
  #aa = subset(aa, subset = nFeature_RNA > 200  & nCount_RNA < 20000 &  percent.mt < 25)
  aa = xx
  rm(xx)
  
  saveRDS(aa, file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM_Cui2020_furtherCleaning.rds'))
  
  aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
  
  plot1 <- VariableFeaturePlot(aa)
  
  top10 <- head(VariableFeatures(aa), 10)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2
  
  all.genes <- rownames(aa)
  aa <- ScaleData(aa, features = all.genes)
  
  aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 30)
  
  aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 50, min.dist = 0.3)
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  ## doublet identification with Doublets
  Removal_Doublets = FALSE
  library(DoubletFinder)
  require(Seurat)
  
  aa = readRDS(file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM_Cui2020_furtherCleaning.rds'))
  
  source('functions_scRNAseq.R')
  xx = detect_doubletCell_scRNAseq(aa)
  aa = xx
  rm(xx)
  
  saveRDS(aa, file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM_Cui2020_furtherCleaning_identDoublets.rds'))
  
  aa = readRDS(file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM_Cui2020_furtherCleaning_identDoublets.rds'))
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'DF_out', raster=FALSE)
  
  p1 + p2  
  
  aa = subset(aa, cells = colnames(aa)[which(aa$DF_out == 'Singlet')])
  
  Idents(aa) = factor(aa$condition)
  p1 = VlnPlot(aa, features = 'nFeature_RNA', y.max = 5000) +
    geom_hline(yintercept = c(200, 500, 1000)) + NoLegend()
  p2 = VlnPlot(aa, features = 'nCount_RNA', y.max = 20000) + NoLegend()
  p3 = VlnPlot(aa, features = 'percent.mt', y.max = 100) + NoLegend()
  
  p1 | p2 | p3
  
  aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000)
  
  plot1 <- VariableFeaturePlot(aa)
  
  top10 <- head(VariableFeatures(aa), 10)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2
  
  all.genes <- rownames(aa)
  aa <- ScaleData(aa, features = all.genes)
  
  aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 30)
  
  aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 50, min.dist = 0.3)
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  saveRDS(aa, file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM_Cui2020_furtherCleaning_rmDoublets.rds'))
  
  
}

Compare_shoval.version_vs.my.version = FALSE
if(Compare_shoval.version_vs.my.version)
{
  aa = readRDS(file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM_Cui2020_furtherCleaning_rmDoublets.rds'))
  
  source('functions_dataIntegration.R')
  ## test batch correction as in the original paper
  ref.combined = IntegrateData_Seurat_CCA(aa, group.by = 'condition')
  aa = ref.combined
  
  aa <- RunUMAP(aa, reduction = "pca", dims = 1:30, n.neighbors = 30, 
                min.dist = 0.3)
  DimPlot(aa, reduction = "umap", group.by = group.by, label = TRUE,
          repel = TRUE, raster=FALSE)
  
  aa <- FindNeighbors(aa, reduction = "pca", dims = 1:30)
  aa <- FindClusters(aa, resolution = 1.0)
  p1 = DimPlot(aa, reduction = "umap", label = TRUE, repel = TRUE, raster=FALSE)
  
  ## use those two markers to clean CMs from contamination
  FeaturePlot(aa, features = c('Myh6', 'Tnnt2'))
  p2 = VlnPlot(aa, features = c('Myh6', 'Tnnt2'))
  
  p1 / p2
  
  
  xx = readRDS(paste0(dataDir, 'CMs.rds'))
  DimPlot(xx, reduction = "umap", label = TRUE, repel = TRUE, raster=FALSE)
  
  xx <- NormalizeData(xx, normalization.method = "LogNormalize", scale.factor = 10000)
  xx <- FindVariableFeatures(xx, selection.method = "vst", nfeatures = 2000)
  xx <- ScaleData(xx)
  
  xx <- RunPCA(xx, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(xx, ndims = 30)
  
  xx <- RunUMAP(xx, dims = 1:30, n.neighbors = 30, min.dist = 0.3)
  
  DimPlot(xx, label = TRUE, repel = TRUE, raster=FALSE)
  
  xx$cell.ids = colnames(xx)
  xx$cell.ids = sapply(xx$cell.ids, function(x){unlist(strsplit(as.character(x), '_'))[1]})
  xx$cell.ids = paste0(xx$cell.ids, '_', xx$orig.ident)
  
  aa$cell.id = gsub('_filtered_gene_bc_matrices_matrix.mtx.gz', '', aa$cell.id)
  aa$cell.id = sapply(aa$cell.id, 
                      function(x){ 
                        x = unlist(strsplit(as.character(x), '_'));
                        x = x[c(1, length(x))];
                        x = paste0(x, collapse = '_')
                        x = gsub('-1', '', x)
                        x
                      })
  
  mm = match(xx$cell.ids, aa$cell.id)
  kk = which(!is.na(mm))
  jj = mm[kk]
  
  xx = subset(xx, cells = colnames(xx)[kk])
  yy = subset(aa, cells = colnames(aa)[jj])
  
  xx <- NormalizeData(xx, normalization.method = "LogNormalize", scale.factor = 10000)
  xx <- FindVariableFeatures(xx, selection.method = "vst", nfeatures = 2000)
  xx <- ScaleData(xx)
  xx <- RunPCA(xx, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(xx, ndims = 30)
  xx <- RunUMAP(xx, dims = 1:30, n.neighbors = 30, min.dist = 0.3)
  DimPlot(xx, label = TRUE, repel = TRUE, raster=FALSE)
  
  DefaultAssay(aa) = 'RNA'
  xx = subset(aa, cells = colnames(aa)[jj])
  xx <- NormalizeData(xx, normalization.method = "LogNormalize", scale.factor = 10000)
  xx <- FindVariableFeatures(xx, selection.method = "vst", nfeatures = 2000)
  xx <- ScaleData(xx)
  xx <- RunPCA(xx, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(xx, ndims = 30)
  xx <- RunUMAP(xx, dims = 1:30, n.neighbors = 30, min.dist = 0.3)
  DimPlot(xx, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
  
}

## the conclusion is that they are the same data; but the cell filtering by Shoval was more stringent
## the gene features are more than my version
## decide to use Shoval's version due to the manual CM cleaning based on Myh6 and Tnnt2 
Process.CM.Cui.et.al.2020_Shoval.version = FALSE
if(Process.CM.Cui.et.al.2020_Shoval.version){
  aa = readRDS(paste0(dataDir, 'CMs.rds'))
  
  aa[["pca_old"]] = aa[["pca"]]
  aa[["umap_old"]] = aa[["umap"]]
  
  aa$condition = aa$orig.ident
  aa$condition = gsub('P11Sham', 'P8_Sham_d3', aa$condition)
  aa$condition = gsub('P1MID1', 'P1_MI_d1', aa$condition)
  aa$condition = gsub('P1MID3', 'P1_MI_d3', aa$condition)
  aa$condition = gsub('P2Sham', 'P1_Sham_d1', aa$condition)
  aa$condition = gsub('P4Sham', 'P1_Sham_d3', aa$condition)
  aa$condition = gsub('P8MID1', 'P8_MI_d1', aa$condition)
  aa$condition = gsub('P8MID3', 'P8_MI_d3', aa$condition)
  aa$condition = gsub('P9Sham', 'P8_Sham_d1', aa$condition)
  
  Idents(aa) = factor(aa$condition)
  p1 = VlnPlot(aa, features = 'nFeature_RNA', y.max = 5000) +
    geom_hline(yintercept = c(200, 500, 1000)) + NoLegend()
  p2 = VlnPlot(aa, features = 'nCount_RNA', y.max = 50000) + NoLegend()
  p3 = VlnPlot(aa, features = 'percent.mt', y.max = 100) + NoLegend()
  p1 | p2 | p3
  
  
  DimPlot(aa, reduction = "umap", label = TRUE, repel = TRUE, raster=FALSE)
  
  aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000)
  aa <- ScaleData(aa)
  
  aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 30)
  
  aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.3)
  
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)
  
  ggsave(filename = paste0(resDir, '/CM_Cui2020_futher_beforeBatchCorrection.pdf'), 
         width = 10, height = 6)
  
  ## run batch correction
  source('functions_dataIntegration.R')
  ref.combined = IntegrateData_Seurat_CCA(aa, group.by = 'condition')
  aa = ref.combined
  
  aa <- RunUMAP(aa, reduction = "pca", dims = 1:30, n.neighbors = 30, 
                min.dist = 0.3)
  DimPlot(aa, reduction = "umap", group.by = group.by, label = TRUE,
          repel = TRUE, raster=FALSE)
  
  saveRDS(aa, file = paste0(RdataDir, 
                            'Seurat.obj_neonatalMice_CM_Cui2020_shoval.version_integrated.rds'))
  
  
}

########################################################
########################################################
# Section III : batch correction and data integration of CMs and noCMs 
#  
########################################################
########################################################

##########################################
# test the batch correction for CM across time points and treatment 
##########################################
process_CM_P1 = FALSE
if(process_CM_P1){
  
  aa = readRDS(file = paste0(RdataDir, 
                             'Seurat.obj_neonatalMice_CM_Cui2020_shoval.version_integrated.rds'))
  
  DefaultAssay(aa) = 'RNA'
  aa$age = sapply(aa$condition, function(x){unlist(strsplit(x, '_'))[1]})
  aa$injury = sapply(aa$condition, function(x){unlist(strsplit(x, '_'))[2]})
  aa$time = sapply(aa$condition, function(x){unlist(strsplit(x, '_'))[3]})
  
  aa = subset(aa, cells = colnames(aa)[which(aa$age == 'P1')])
  
  DimPlot(aa)
  
  aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000)
  aa <- ScaleData(aa)
  
  aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 30)
  
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.3)
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)
  
  ggsave(filename = paste0(resDir, '/CM_Cui2020_P1_noDataIntegration.pdf'), 
         width = 10, height = 6)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
  aa <- FindNeighbors(aa, reduction = "pca", dims = 1:10)
  aa <- FindClusters(aa, resolution = 0.4)
  
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)
  
  p1 + p2
  
  aa$celltypes = paste0('CM', aa$seurat_clusters)
  
  saveRDS(aa, file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM_Cui2020_shoval.version_P1_clusterd.rds'))
  
  ##########################################
  # test batch correction with different integration methods
  # ideally 
  ##########################################
  Test_DataIntegration = FALSE
  if(Test_DataIntegration){
    require(dplyr)
    require(tibble)
    require(ggplot2)
    library(viridis)
    library(hrbrthemes)
    
    aa = readRDS(file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM_Cui2020_shoval.version_P1_clusterd.rds'))
    aa$condition = factor(aa$condition, levels = c('P1_Sham_d3', 'P1_Sham_d1', 'P1_MI_d1', 'P1_MI_d3'))
    
    outDir = paste0(resDir, '/CM_Cui2020_P1/')
    if(!dir.exists(outDir)) dir.create(outDir)
    
    DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
    DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition', split.by = 'condition')
    
    group.by = 'condition'
    feature_sels = c('Nppa', 'Nppb', 'Ankrd1', 'Foxp1', 'Myh7', 'Des', 'Uchl1')
    integration_methods = c('noDataIntegration', 'Seurat_CCA', 'fastMNN') 
    
    for(method in integration_methods)
    {
      #method = "fastMNN"
      cat('------------------------\n')
      cat('method -- ', method, '\n')
      cat('------------------------\n')
      source('functions_dataIntegration.R')
      
      ## no data integration 
      if(method == 'noDataIntegration'){ref.combined = aa}
      
      if(method == 'Seurat_CCA'){
        ref.combined = IntegrateData_Seurat_CCA(aa, 
                                                group.by = group.by, 
                                                merge.order = matrix(c(-2, -3, 1, -1, -4, 2), ncol = 2),
                                                redo.normalization.scaling = FALSE,
                                                correct.all = TRUE)
      }
      
      if(method == 'Seurat_RPCA'){
        ref.combined = IntegrateData_Seurat_RPCA(aa, group.by = group.by, 
                                                 redo.normalization.scaling = FALSE,
                                                 correct.all = TRUE)
      }
      
      if(method == 'fastMNN'){
        ref.combined = IntegrateData_runFastMNN(aa, group.by = group.by,
                                                merge.order = c('P1_Sham_d3', 'P1_Sham_d1',
                                                                'P1_MI_d1', 'P1_MI_d3'),
                                                correct.all = TRUE)
      }
      
      p1 = DimPlot(ref.combined, group.by = 'celltypes', label = TRUE, repel = TRUE, raster=FALSE) + 
        ggtitle(method)
      p2 = DimPlot(ref.combined, group.by = 'condition', label = TRUE, repel = TRUE) +
        ggtitle(method)
      p1 + p2
      ggsave(filename = paste0(outDir, '/CM_Cui2020_P1_merged_dataIntegration_', method, '.pdf'), 
             width = 16, height = 6)
      
      FeaturePlot(ref.combined, features = feature_sels)
      ggsave(filename = paste0(outDir, '/CM_Cui2020_P1_merged_borderMarkers_dataIntegration_', method, '.pdf'), 
             width = 16, height = 12)
      
      if(method == 'fastMNN'){
        ref.combined <- FindNeighbors(ref.combined, reduction = "mnn", dims = 1:10)
        ref.combined <- FindClusters(ref.combined, resolution = 0.4)
      }else{
        ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:10)
        ref.combined <- FindClusters(ref.combined, resolution = 0.4)
      }
      
      p1 = DimPlot(ref.combined, group.by = 'seurat_clusters', label = TRUE, repel = TRUE)
      p2 = DimPlot(ref.combined, group.by = 'seurat_clusters', split.by = group.by, label = TRUE, 
                   repel = TRUE)
      p1 / p2
      ggsave(filename = paste0(outDir, '/CM_Cui2020_P1_merged_clusters.by.dataset', 
                               '_dataIntegration_', method, '.pdf'), 
             width = 16, height = 12)
      
      fraction = table(ref.combined$seurat_clusters, ref.combined$condition)
      for(n in 1:ncol(fraction)) fraction[,n] = fraction[,n]/sum(fraction[,n])
      as.data.frame(fraction) %>% 
        ggplot(aes(fill=Var1, y=Freq, x=Var2)) + 
        geom_bar(position="fill", stat="identity")
        #scale_fill_viridis(discrete = T)
      
      ggsave(filename = paste0(outDir, '/CM_Cui2020_P1_merged_cluster_fractions', 
                               '_dataIntegration_', method, '.pdf'), 
             width = 12, height = 6)
      
      if(method != 'noDataIntegration'){
        saveRDS(ref.combined, 
                file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM_Cui2020_P1_dataIntegration_', 
                              method, '.rds'))
      }
    }
    
  }
}

##########################################
# batch correction test for noCM  
##########################################
process_and_batchCorrection_noCM = FALSE
if(process_and_batchCorrection_noCM){
  
  aa = readRDS(file = paste0(RdataDir, 'Wang_DevCell_2020_scRNAseq.rds'))
  aa$treatment = aa$condition
  aa$condition = paste0(aa$stage, '_', aa$condition, '_', aa$time)
  
  DimPlot(aa, group.by = 'FineID', label = TRUE, repel = TRUE, raster=FALSE)
  DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE, raster=FALSE)
  
  ## subset cells only from P1
  aa = subset(aa, cells = colnames(aa)[which(aa$stage == 'P1')])
  
  aa = DietSeurat(aa, counts = TRUE, data = TRUE, scale.data = TRUE, assays = 'RNA')
  aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000)
  aa <- ScaleData(aa, vars.to.regress = c('nUMI', 'percent.mito'), features = rownames(aa))
  
  aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 30)
  
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.3)
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)
  
  ggsave(filename = paste0(resDir, '/Wang2020_P1_noDataIntegration.pdf'), 
         width = 10, height = 6)
  
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
  aa <- FindNeighbors(aa, reduction = "pca", dims = 1:10)
  aa <- FindClusters(aa, resolution = 0.4)
  
  saveRDS(aa, file = paste0(RdataDir, 'Wang_DevCell_2020_scRNAseq_P1_regressedout.rds'))
  
  ##########################################
  # test batch correction within the dataset
  ##########################################
  Test_DataIntegration = FALSE
  if(Test_DataIntegration){
    require(dplyr)
    require(tibble)
    require(ggplot2)
    library(viridis)
    library(hrbrthemes)
    
    aa = readRDS(file = paste0(RdataDir, 'Wang_DevCell_2020_scRNAseq_P1_regressedout.rds'))
    aa$condition = factor(aa$condition, levels = c('P1_Sham_D1', 'P1_Sham_D3',  'P1_MI_D1', 'P1_MI_D3'))
    
    outDir = paste0(resDir, '/noCM_Wang2020_P1/')
    if(!dir.exists(outDir)) dir.create(outDir)
    
    DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
    DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'FineID', split.by = 'condition')
    
    group.by = 'condition'
    feature_sels = c('Nppa', 'Nppb', 'Ankrd1', 'Foxp1', 'Myh7', 'Des', 'Uchl1')
    integration_methods = c('noDataIntegration', 'Seurat_CCA', 'fastMNN') 
    
    for(method in integration_methods)
    {
      # method = "Seurat_CCA"
      
      cat('------------------------\n')
      cat('method -- ', method, '\n')
      cat('------------------------\n')
      source('functions_dataIntegration.R')
      
      ## no data integration 
      if(method == 'noDataIntegration'){ref.combined = aa}
      
      if(method == 'Seurat_CCA'){
        ref.combined = IntegrateData_Seurat_CCA(aa, 
                                                group.by = group.by, 
                                                merge.order = matrix(c(-1, -3, 1, -2, -4, 2), ncol = 2),
                                                redo.normalization.scaling = FALSE,
                                                correct.all = TRUE)
      }
      
      if(method == 'Seurat_RPCA'){
        ref.combined = IntegrateData_Seurat_RPCA(aa, group.by = group.by, 
                                                 redo.normalization.scaling = FALSE,
                                                 correct.all = TRUE)
      }
      
      if(method == 'fastMNN'){
        ref.combined = IntegrateData_runFastMNN(aa, group.by = group.by,
                                                merge.order = c('P1_Sham_D3', 'P1_Sham_D1',
                                                                'P1_MI_D1', 'P1_MI_D3'),
                                                correct.all = TRUE)
      }
      
      
      p1 = DimPlot(ref.combined, group.by = 'FineID', label = TRUE, repel = TRUE, raster=FALSE) + 
        ggtitle(method)
      p2 = DimPlot(ref.combined, group.by = 'condition', label = TRUE, repel = TRUE) +
        ggtitle(method)
      p1 + p2
      ggsave(filename = paste0(outDir, '/noCM_Wang2020_P1_merged_dataIntegration_', method, '.pdf'), 
             width = 16, height = 6)
      
      #FeaturePlot(ref.combined, features = feature_sels)
      #ggsave(filename = paste0(outDir, '/CM_Cui2020_P1_merged_borderMarkers_dataIntegration_', method, '.pdf'), 
      #       width = 16, height = 12)
      
      if(method == 'fastMNN'){
        ref.combined <- FindNeighbors(ref.combined, reduction = "mnn", dims = 1:10)
        ref.combined <- FindClusters(ref.combined, resolution = 1.0)
      }else{
        ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:20)
        ref.combined <- FindClusters(ref.combined, resolution = 1.0)
      }
      
      p1 = DimPlot(ref.combined, group.by = 'seurat_clusters', label = TRUE, repel = TRUE)
      p2 = DimPlot(ref.combined, group.by = 'seurat_clusters', split.by = group.by, label = TRUE, 
                   repel = TRUE)
      p1 / p2
      ggsave(filename = paste0(outDir, '/noCM_Wang2020_P1_merged_clusters.by.dataset', 
                               '_dataIntegration_', method, '.pdf'), 
             width = 16, height = 12)
      
      fraction = table(ref.combined$seurat_clusters, ref.combined$condition)
      for(n in 1:ncol(fraction)) fraction[,n] = fraction[,n]/sum(fraction[,n])
      as.data.frame(fraction) %>% 
        ggplot(aes(fill=Var1, y=Freq, x=Var2)) + 
        geom_bar(position="fill", stat="identity")
      #scale_fill_viridis(discrete = T)
      
      ggsave(filename = paste0(outDir, '/noCM_Wang2020_P1_merged_cluster_fractions', 
                               '_dataIntegration_', method, '.pdf'), 
             width = 12, height = 6)
      
      if(method != 'noDataIntegration'){
        saveRDS(ref.combined, 
                file = paste0(RdataDir, 'Seurat.obj_neonatalMice_noCM_Wang2020_P1_dataIntegration_', 
                              method, '.rds'))
      }
    }
    
  }
}


########################################################
########################################################
## Section IV : merge CMs and noCMs
## also test how to correct the batches and data integration 
## 1) batch corrected separately for CM and noCM and merged 
## 2) merge CM and noCM first and correct the batch in data integration
########################################################
########################################################
cms = readRDS(file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM_Cui2020_shoval.version_P1_clusterd.rds'))
aa = readRDS(file = paste0(RdataDir, 'Wang_DevCell_2020_scRNAseq_P1_regressedout.rds'))


##########################################
# merge CMs and noCMs 
# also test if correct the batch effect on each dataset
##########################################

cms[['integrated']] <- NULL

cms$dataset = 'Cui2020'
aa$dataset = 'Wang2020'

cms$BroadID = 'CM'
cms$FineID = cms$celltypes
cms$nUMI = cms$nCount_RNA
aa$celltypes = aa$BroadID
cms$percent.mito = cms$percent.mt

features.common = intersect(rownames(aa), rownames(cms))

aa = subset(aa, features = features.common)
cms = subset(cms, features = features.common)

aa = merge(aa, y = cms, add.cell.ids = c("Wang2020", "Cui2020"), project = "neonatalMice")
aa$condition = gsub('_D', '_d', aa$condition)

aa$condition = factor(aa$condition, levels = c('P1_Sham_d1', 'P1_Sham_d3',  'P1_MI_d1', 'P1_MI_d3'))

aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
aa <- ScaleData(aa, vars.to.regress = c('nCount_RNA'), features = features.common)

aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
ElbowPlot(aa, ndims = 30)

aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.3)
DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)

saveRDS(aa, file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_regress.nUMI_v2.rds'))

##########################################
# test Data integration after data merging 
##########################################
Test_DataIntegration = FALSE
if(Test_DataIntegration){
  
  aa = readRDS(file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_regress.nUMI.rds'))
  
  aa$dataset = factor(aa$dataset, levels = c('Wang2020', 'Cui2020'))
  aa$batch = paste0(aa$dataset, '_', aa$condition)
  aa$batch = factor(aa$batch)
  
  outDir = paste0(resDir, '/CM_Cui2020_noCM_Wang2020_P1_merged/')
  if(!dir.exists(outDir)) dir.create(outDir)
  
  integration_methods = c('noDataIntegration',  
                          'fastMNN',
                          'Seurat_RPCA')
  
  for(method in integration_methods)
  {
    
    method = "Seurat_RPCA"
    
    source('functions_dataIntegration.R')
    
    ## no data integration 
    if(method == 'noDataIntegration'){ref.combined = aa}
    
    if(method == 'Seurat_RPCA'){
      ref.combined = IntegrateData_Seurat_RPCA(aa, 
                                               group.by = 'batch', 
                                               nfeatures = 3000,
                                               merge.order = matrix(c(-2, -3, 1, -8, -5, 5, 3,
                                                                      -1, -4, 2, -7, -6, 4, 6), ncol = 2),
                                               redo.normalization.scaling = TRUE,
                                               correct.all = TRUE)
      
      DefaultAssay(ref.combined) = 'integrated'
      
    }
    
    if(method == 'fastMNN'){
      # 
      source('functions_dataIntegration.R')
      ref.combined = IntegrateData_runFastMNN(aa, group.by = 'batch', 
                                              nfeatures = 3000,
                                              merge.order = list(list(4,3,1,2), list(8,7,5,6)),
                                              correct.all = TRUE)
      DefaultAssay(ref.combined) = 'mnn.reconstructed'
      
    }
    
    p1 = DimPlot(ref.combined, group.by = 'FineID', label = TRUE, repel = TRUE, raster=FALSE) + 
      ggtitle(method)
    p2 = DimPlot(ref.combined, group.by = 'condition', label = TRUE, repel = TRUE) +
      ggtitle(method)
    p3 = DimPlot(ref.combined, group.by = 'BroadID', label = TRUE, repel = TRUE, raster=FALSE) + 
      ggtitle(method)
    
    p1 /(p3 + p2)
    
    ggsave(filename = paste0(outDir, '/CM_Cui2020_noCM_Wang2020_P1_merged_dataIntegration_', method, '.pdf'), 
           width = 16, height = 12)
    
    DimPlot(ref.combined, group.by = 'FineID', split.by = 'dataset', label = TRUE, repel = TRUE)
    ggsave(filename = paste0(outDir, '/CM_Cui2020_noCM_Wang2020_P1_merged_subtypes.by.dataset', 
                             '_dataIntegration_', method, '.pdf'), 
           width = 16, height = 6)
    
    DimPlot(ref.combined, group.by = 'BroadID', split.by = 'dataset', label = TRUE, repel = TRUE)
    ggsave(filename = paste0(outDir, '/CM_Cui2020_noCM_Wang2020_P1_merged_celltypes.by.dataset', 
                             '_dataIntegration_', method, '.pdf'), 
           width = 16, height = 6)
    
    DimPlot(ref.combined, group.by = 'FineID', split.by = 'condition', label = TRUE, repel = TRUE, ncol = 2)
    ggsave(filename = paste0(outDir, '/CM_Cui2020_noCM_Wang2020_P1_merged_subtypes.by.conditions', 
                             '_dataIntegration_', method, '.pdf'),  width = 16, height = 12)
    
    DimPlot(ref.combined, group.by = 'BroadID', split.by = 'condition', label = TRUE, repel = TRUE, ncol = 2)
    ggsave(filename = paste0(outDir, '/CM_Cui2020_noCM_Wang2020_P1_merged_celltypes.by.conditions', 
                             '_dataIntegration_', method, '.pdf'),  width = 16, height = 12)
    
    if(method != 'noDataIntegration'){
      saveRDS(ref.combined, 
              file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_dataIntegration_', 
                            method, '_allgeneCorrected.rds'))
      
    }
    
  }
  
}


########################################################
########################################################
# Section V : refine the subtypes after fastMNN integration of CM and noCM
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 
                           'Seurat.obj_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_dataIntegration_',
                           'fastMNN_allgeneCorrected.rds'))


## overview of fastMNN integerated data
aa$time[which(aa$time == 'D1')] = 'd1'
aa$time[which(aa$time == 'D3')] = 'd3'

aa$timepoints = aa$time
aa$condition = gsub('P1_', '', aa$condition)

aa$celltype = aa$BroadID
aa$subtype = aa$FineID
  
p1 = DimPlot(aa, reduction = 'umap', group.by = 'dataset', raster = T, pt.size = 2, shuffle= T,label = T)
p2 = DimPlot(aa, reduction = 'umap', group.by = 'celltype', raster = T, shuffle= T,label = T)
p3 = DimPlot(aa, reduction = 'umap', group.by = 'subtype', raster = T, shuffle= T,label = T)
p4 = DimPlot(aa, reduction = 'umap', group.by = 'timepoints', raster = T, shuffle= T,label = T)

(p1 + p2)/(p3 + p4)

ggsave(filename = paste0(resDir, '/UMAP_scRNAseq_refrence_dataset_timepoints_celltypes_v1.pdf'), 
       width = 32, height = 16)


aa$subtype_old = aa$subtype
aa$subtype = NA
aa$strigentFiltered = NA

s.genes <- firstup(cc.genes$s.genes)
g2m.genes <- firstup(cc.genes$g2m.genes)
aa <- CellCycleScoring(aa, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

aa = subset(aa,  cells = colnames(aa)[which(aa$celltype != 'glial')])
aa$celltypes <- NULL
jj = which(aa$dataset == 'Wang2020')
aa$percent.mt[jj] = aa$percent.mito[jj] * 100 
aa$percent.mito = NULL

DefaultAssay(aa) = 'mnn.reconstructed'
#aa@assays$mnn.reconstructed@counts = aa@assays$RNA@counts

##########################################
# double check and refine subtypes CM
##########################################
refine.subtypes_adultMice_CM = FALSE
if(refine.subtypes_adultMice_CM){
  
  mcells = 'CM'
  outDir = paste0(resDir, '/', mcells)
  if(!dir.exists(outDir)) dir.create(outDir)
  
  ax = subset(aa, cells = colnames(aa)[which(aa$celltype == mcells)])
  table(ax$celltype)
  table(ax$subtype_old)
  
  DefaultAssay(ax) = 'mnn.reconstructed'
  
  ax$condition = factor(ax$condition, levels = c("Sham_d1", "MI_d1", "Sham_d3", "MI_d3"))
  
  ax <- FindVariableFeatures(ax, selection.method = "vst", nfeatures = 2000)
  ax <- ScaleData(ax)
  
  ax <- RunPCA(ax, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(ax, ndims = 30)
  
  # UMAP to visualize subtypes
  ax <- RunUMAP(ax, dims = 1:10, n.neighbors = 20, min.dist = 0.1)
  
  DimPlot(ax, reduction = 'umap', group.by = 'subtype_old')
  
  ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, '_subcelltypes.pdf'), 
         width = 10, height = 8)
  
  FeaturePlot(ax, features = c("Myh6", 'Nppa', 'Agrn'))
  ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, '_MarkerGenes.pdf'), 
         width = 10, height = 8)
  
  ax <- FindNeighbors(ax, dims = 1:10)
  
  ax <- FindClusters(ax, verbose = FALSE, algorithm = 3, resolution = 0.3)
  p1 = DimPlot(ax, reduction = 'umap', group.by = 'seurat_clusters')
  p0 = DimPlot(ax, reduction = 'umap', group.by = 'subtype')
  
  p0 + p1
  
  VlnPlot(ax, features = 'Nppa')
  DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4,group.by = "condition")
  
  ## check the quality of CM
  pdf(paste0(outDir, "/QCs_doubleCheck.pdf"), 
      width=10, height=8)
  DimPlot(object = ax, reduction = 'umap', raster = T, pt.size = 4, shuffle= T,label = T)
  
  DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4,group.by = "condition")
  DimPlot(object = ax, reduction = 'umap',raster = T, shuffle= T, pt.size = 4,group.by = "dataset")
  #DimPlot(object = xx, reduction = 'umap',raster = T, shuffle= T, pt.size = 3,group.by = "time",
  #        cols = c("grey50",rev(viridis_pal(option = "plasma")(length(unique(xx@meta.data$time))-1))))
  FeaturePlot(ax, raster = T,pt.size = 2, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
              order = T)
  
  VlnPlot(ax, features = 'nFeature_RNA')
  VlnPlot(ax, features = 'nCount_RNA')
  VlnPlot(ax, features = 'percent.mt')
  
  dev.off()
  
  Test_regressOut_factors = FALSE
  if(Test_regressOut_factors){
    
    p1 = VlnPlot(ax, features = 'nFeature_RNA', group.by = 'condition') +
      geom_hline(yintercept= c(500, 800, 1000))
    p2 = VlnPlot(ax, features = 'nCount_RNA',  group.by = 'condition')
    p3 = VlnPlot(ax, features = 'percent.mt',  group.by = 'condition') + 
      geom_hline(yintercept= c(50, 60))
    
    p1 + p2 + p3
    
    ax = subset(ax, cells = colnames(ax)[which(ax$nFeature_RNA > 500 & ax$seurat_clusters != '5')])
    
    ax <- FindVariableFeatures(ax, selection.method = "vst", nfeatures = 2000)
    ax <- ScaleData(ax)
    
    ax <- RunPCA(ax, verbose = FALSE, weight.by.var = TRUE)
    ElbowPlot(ax, ndims = 30)
    
    
    # UMAP to visualize subtypes
    ax <- RunUMAP(ax, dims = 1:10, n.neighbors = 20, min.dist = 0.1)
    
    DimPlot(ax, reduction = 'umap', group.by = 'subtype_old')
    
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, '_subcelltypes.pdf'), 
           width = 10, height = 8)
    
    FeaturePlot(ax, features = c("Myh6", 'Nppa', 'Agrn'))
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, '_MarkerGenes.pdf'), 
           width = 10, height = 8)
    
    ax <- FindNeighbors(ax, dims = 1:10)
    ax <- FindClusters(ax, verbose = FALSE, algorithm = 3, resolution = 0.3)
    
    p1 = DimPlot(ax, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE, repel = TRUE)
    p0 = DimPlot(ax, reduction = 'umap', group.by = 'subtype_old', label = TRUE, repel = TRUE)
    p2 = DimPlot(ax, reduction = 'umap', group.by = 'Phase')
    p0 + p1 + p2
    
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                  '_afterFiltering.highMt.pdf'), 
           width = 12, height = 6)
    
    
    VlnPlot(ax, features = 'Nppa')
    
    DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4, group.by = "seurat_clusters",
            split.by = 'condition', label = TRUE, repel = TRUE)
    
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                  '_afterFiltering.highMt_clusterPerCondition.pdf'), 
           width = 14, height = 8)
    
  }
  
  ## re-annotate
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(0))))] = 'CM1'
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(1))))] = 'CM2_injurySpec_proliferating'
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(2))))] = 'CM3_injurySpec'
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(3))))] = 'CM4_injurySpec_d1'
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(4))))] = 'CM5_proliferating'
  
  DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4, group.by = "subtype",
          label = TRUE, repel = TRUE)
  
  ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                '_afterFiltering_reAnnot.pdf'), 
         width = 12, height = 8)
  
  # save the annotation and processing 
  mm = match(colnames(ax), colnames(aa))
  cat(length(which(is.na(mm))), 'cells missing \n')
  aa$subtype[mm] = ax$subtype
  aa$strigentFiltered[mm] = 'keep'
  
  saveRDS(ax, file = paste0(RdataDir, 
                            'Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_', mcells, 'refinedCelltypes_20240130.rds'))
  
  saveRDS(aa, file = paste0(RdataDir, 
                            'Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_refineSubtypes_', mcells, '_20240130.rds'))
  
}

##########################################
# EC 
##########################################
DimPlot(aa, reduction = 'umap', group.by = 'celltype', label = TRUE, repel = TRUE)

refine.subtypes_adultMice_EC = FALSE
if(refine.subtypes_adultMice_EC){
  
  mcells = 'EC'
  outDir = paste0(resDir, '/', mcells)
  if(!dir.exists(outDir)) dir.create(outDir)
  
  ax = subset(aa, cells = colnames(aa)[which(aa$celltype == mcells)])
  table(ax$celltype)
  table(ax$subtype_old)
  
  ax$condition = factor(ax$condition, levels = c("Sham_d1", "MI_d1", "Sham_d3", "MI_d3"))
  
  DimPlot(ax, group.by = 'subtype_old')
  
  ax <- FindVariableFeatures(ax, selection.method = "vst", nfeatures = 2000)
  ax <- ScaleData(ax)
  
  ax <- RunPCA(ax, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(ax, ndims = 30)
  
  # UMAP to visualize subtypes
  ax <- RunUMAP(ax, dims = 1:10, n.neighbors = 20, min.dist = 0.1)
  
  DimPlot(ax, reduction = 'umap', group.by = 'subtype_old')
  
  ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, '_subcelltypes.pdf'), 
         width = 10, height = 8)
  
  #FeaturePlot(ax, features = c("Myh6", 'Nppa', 'Agrn'))
  #ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, '_MarkerGenes.pdf'), 
  #       width = 10, height = 8)
  
  ax <- FindNeighbors(ax, dims = 1:10)
  
  ax <- FindClusters(ax, verbose = FALSE, algorithm = 3, resolution = 0.4)
  p0 = DimPlot(ax, reduction = 'umap', group.by = 'subtype_old')
  p1 = DimPlot(ax, reduction = 'umap', group.by = 'seurat_clusters')
  
  p0 + p1
  
  ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, '_redo.clustering.pdf'), 
         width = 16, height = 8)
  
  #VlnPlot(ax, features = 'Nppa')
  DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4, group.by = "condition")
  
  ## check the quality of CM
  pdf(paste0(outDir, "/QCs_doubleCheck.pdf"), 
      width=10, height=8)
  
  DimPlot(ax, reduction = 'umap', group.by = 'subtype_old')
  DimPlot(object = ax, reduction = 'umap', raster = T, pt.size = 4, shuffle= T,label = T)
  #DimPlot(ax, reduction = 'umap', group.by = 'Phase')
  DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4,group.by = "condition")
  DimPlot(object = ax, reduction = 'umap',raster = T, shuffle= T, pt.size = 4,group.by = "Phase")
  #DimPlot(object = xx, reduction = 'umap',raster = T, shuffle= T, pt.size = 3,group.by = "time",
  #        cols = c("grey50",rev(viridis_pal(option = "plasma")(length(unique(xx@meta.data$time))-1))))
  FeaturePlot(ax, raster = T,pt.size = 2, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
              order = T)
  
  VlnPlot(ax, features = 'nFeature_RNA')
  VlnPlot(ax, features = 'nCount_RNA')
  VlnPlot(ax, features = 'percent.mt')
  
  dev.off()
  
  Test_filtering.lowQualityCells_regressOut.factors = FALSE
  if(Test_filtering.lowQualityCells_regressOut.factors){
    
    groupby = 'seurat_clusters'
    p1 = VlnPlot(ax, features = 'nFeature_RNA', group.by = groupby) +
      geom_hline(yintercept= c(500, 800, 1000))
    p2 = VlnPlot(ax, features = 'nCount_RNA',  group.by = groupby) +
      geom_hline(yintercept= c(1000, 2000))
    p3 = VlnPlot(ax, features = 'percent.mt',  group.by = groupby) + 
      geom_hline(yintercept= c(50, 60))
    
    p1 + p2 + p3
    
    # stringent filtering for EC 
    ax = subset(ax, cells = colnames(ax)[which(ax$nCount_RNA >1000 & ax$nFeature_RNA > 500)])
    
    ax <- FindVariableFeatures(ax, selection.method = "vst", nfeatures = 2000)
    ax <- ScaleData(ax)
    
    ax <- RunPCA(ax, verbose = FALSE, weight.by.var = TRUE)
    ElbowPlot(ax, ndims = 30)
    
    # UMAP to visualize subtypes
    ax <- RunUMAP(ax, dims = 1:20, n.neighbors = 20, min.dist = 0.1)
    
    DimPlot(ax, reduction = 'umap', group.by = 'subtype_old')
    
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, '_afterFiltering_subcelltypes.pdf'), 
           width = 10, height = 8)
    
    FeaturePlot(ax, features = c("Myh6", 'Nppa', 'Agrn'))
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, '_MarkerGenes.pdf'), 
           width = 10, height = 8)
    
    ax <- FindNeighbors(ax, dims = 1:20)
    ax <- FindClusters(ax, verbose = FALSE, algorithm = 3, resolution = 0.3)
    
    p1 = DimPlot(ax, reduction = 'umap', group.by = 'seurat_clusters',  label = TRUE, repel = TRUE)
    p0 = DimPlot(ax, reduction = 'umap', group.by = 'subtype_old',  label = TRUE, repel = TRUE)
    p2 = DimPlot(ax, reduction = 'umap', group.by = 'Phase')
    p0 + p1 + p2
    
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                  '_afterFiltering_subcelltypes_reclustering.pdf'), 
           width = 20, height = 8)
    
    
    DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4, group.by = "seurat_clusters",
            split.by = 'condition', label = TRUE, repel = TRUE, ncol = 2)
    
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                  '_afterFiltering_clusterPerCondition.pdf'), 
           width = 16, height = 12)
    
  }
  
  ## won't change the previous annotations
  ax$subtype = ax$subtype_old
  
  DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4, group.by = "subtype",
          label = TRUE, repel = TRUE)
  
  ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                '_afterFiltering_reAnnot.pdf'), 
         width = 12, height = 8)
  
  # save the annotation and processing 
  mm = match(colnames(ax), colnames(aa))
  cat(length(which(is.na(mm))), 'cells missing \n')
  aa$subtype[mm] = ax$subtype
  aa$strigentFiltered[mm] = 'keep'
  
  saveRDS(ax, file = paste0(RdataDir, 
                            'Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_', mcells, 'refinedCelltypes_20240130.rds'))
  
  saveRDS(aa, file = paste0(RdataDir, 
                            'Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_refineSubtypes_', mcells, '_20240130.rds'))
  
  
}


##########################################
# FB
##########################################
DimPlot(aa, reduction = 'umap', group.by = 'celltype', label = TRUE, repel = TRUE)

refine.subtypes_adultMice_FB = FALSE
if(refine.subtypes_adultMice_FB){
  
  mcells = 'FB'
  outDir = paste0(resDir, '/', mcells)
  if(!dir.exists(outDir)) dir.create(outDir)
  
  ax = subset(aa, cells = colnames(aa)[which(aa$celltype == mcells)])
  table(ax$celltype)
  table(ax$subtype_old)
  
  ax$condition = factor(ax$condition, levels = c("Sham_d1", "MI_d1", "Sham_d3", "MI_d3"))
  
  DimPlot(ax, group.by = 'subtype_old')
  
  ax <- FindVariableFeatures(ax, selection.method = "vst", nfeatures = 2000)
  ax <- ScaleData(ax)
  
  ax <- RunPCA(ax, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(ax, ndims = 30)
  
  # UMAP to visualize subtypes
  ax <- RunUMAP(ax, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
  
  DimPlot(ax, reduction = 'umap', group.by = 'subtype_old')
  
  ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, '_subcelltype_old.pdf'), 
         width = 10, height = 8)
  
  #FeaturePlot(ax, features = c("Myh6", 'Nppa', 'Agrn'))
  #ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, '_MarkerGenes.pdf'), 
  #       width = 10, height = 8)
  
  ax <- FindNeighbors(ax, dims = 1:20)
  
  ax <- FindClusters(ax, verbose = FALSE, algorithm = 3, resolution = 0.3)
  p0 = DimPlot(ax, reduction = 'umap', group.by = 'subtype_old', label = TRUE, repel = FALSE)
  p1 = DimPlot(ax, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE, repel = FALSE)
  
  p0 + p1
  
  ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, '_redo.clustering.pdf'), 
         width = 16, height = 8)
  
  ###### Filtering the cluster 8, 9
  Quick_Filtering = FALSE
  if(Quick_Filtering){
    ax = subset(ax, idents = c("5"), invert = TRUE)
    
    ax <- FindVariableFeatures(ax, selection.method = "vst", nfeatures = 3000)
    ax <- ScaleData(ax)
    ax <- RunPCA(ax, verbose = FALSE, weight.by.var = TRUE)
    ElbowPlot(ax, ndims = 30)
    
    ax <- RunUMAP(ax, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
    
    DimPlot(ax, reduction = 'umap', group.by = 'subtype_old')
    
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                  '_subcelltype_old_1stFiltering.pdf'), 
           width = 10, height = 8)
    
    ax <- FindNeighbors(ax, dims = 1:20)
    
    ax <- FindClusters(ax, verbose = FALSE, algorithm = 3, resolution = 0.4)
    p0 = DimPlot(ax, reduction = 'umap', group.by = 'subtype_old', label = TRUE, repel = FALSE)
    p1 = DimPlot(ax, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE, repel = FALSE)
    
    p0 + p1
    
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                  '_redo.clustering_quickFiltering.pdf'), 
           width = 16, height = 8)
    
  }
  
  DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4, group.by = "condition")
  
  ## check the quality of CM
  pdf(paste0(outDir, "/QCs_doubleCheck.pdf"), 
      width=10, height=8)
  
  DimPlot(ax, reduction = 'umap', group.by = 'subtype_old', label = T, repel = TRUE)
  DimPlot(object = ax, reduction = 'umap', raster = T, pt.size = 4, shuffle= T,label = T)
  DimPlot(ax, reduction = 'umap', group.by = 'Phase')
  DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4,group.by = "condition")
  #DimPlot(object = ax, reduction = 'umap',raster = T, shuffle= T, pt.size = 4,group.by = "group")
  #DimPlot(object = xx, reduction = 'umap',raster = T, shuffle= T, pt.size = 3,group.by = "time",
  #        cols = c("grey50",rev(viridis_pal(option = "plasma")(length(unique(xx@meta.data$time))-1))))
  FeaturePlot(ax, raster = T,pt.size = 2, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
              order = T)
  
  VlnPlot(ax, features = 'nFeature_RNA')
  VlnPlot(ax, features = 'nCount_RNA')
  VlnPlot(ax, features = 'percent.mt')
  
  dev.off()
  
  Test_filtering.lowQualityCells_regressOut.factors = FALSE
  if(Test_filtering.lowQualityCells_regressOut.factors){
    
    groupby = 'seurat_clusters'
    p1 = VlnPlot(ax, features = 'nFeature_RNA', group.by = groupby) +
      geom_hline(yintercept= c(500, 800, 1000))
    p2 = VlnPlot(ax, features = 'nCount_RNA',  group.by = groupby) +
      geom_hline(yintercept= c(1000, 2000))
    p3 = VlnPlot(ax, features = 'percent.mt',  group.by = groupby) + 
      geom_hline(yintercept= c(50, 60))
    
    p1 + p2 + p3
    
    # stringent filtering for EC and discard cluster 1
    ax = subset(ax, cells = colnames(ax)[which(ax$nCount_RNA >1000 & ax$nFeature_RNA > 500
                                                 )])
    
    ax <- FindVariableFeatures(ax, selection.method = "vst", nfeatures = 2000)
    ax <- ScaleData(ax)
    
    ax <- RunPCA(ax, verbose = FALSE, weight.by.var = TRUE)
    ElbowPlot(ax, ndims = 30)
    
    # UMAP to visualize subtypes
    ax <- RunUMAP(ax, dims = 1:20, n.neighbors = 20, min.dist = 0.1)
    
    DimPlot(ax, reduction = 'umap', group.by = 'subtype_old')
    
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                  '_afterFiltering_subcelltypes.pdf'), 
           width = 10, height = 8)
    
    #FeaturePlot(ax, features = c("Myh6", 'Nppa', 'Agrn'))
    #ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, '_MarkerGenes.pdf'), 
    #       width = 10, height = 8)
    
    ax <- FindNeighbors(ax, dims = 1:20)
    ax <- FindClusters(ax, verbose = FALSE, algorithm = 3, resolution = 1.0)
    
    p1 = DimPlot(ax, reduction = 'umap', group.by = 'seurat_clusters',  label = T, repel = TRUE)
    p0 = DimPlot(ax, reduction = 'umap', group.by = 'subtype_old',  label = T, repel = TRUE)
    p2 = DimPlot(ax, reduction = 'umap', group.by = 'Phase')
    p0 + p1 + p2
    
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                  '_afterFiltering_subcelltypes_reclustering.pdf'), 
           width = 20, height = 8)
    
    
    DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4, group.by = "seurat_clusters",
            split.by = 'condition', label = TRUE, repel = TRUE, ncol = 2)
    
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                  '_afterFiltering_clusterPerCondition.pdf'), 
           width = 16, height = 12)
    
  }
  
  p1 = DimPlot(ax, reduction = 'umap', group.by = 'seurat_clusters',  label = T, repel = TRUE)
  p0 = DimPlot(ax, reduction = 'umap', group.by = 'subtype_old',  label = T, repel = TRUE)
  p0 + p1
  
  ## re-annotate
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(0))))] = 'FB4'
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(1))))] = 'FB3'
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(2, 6))))] = 'FB2'
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(3, 4))))] = 'FB1'
  
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(5))))] = 'FB5_pro'
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(7))))] = 'FB6_injurySpec'
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(8))))] = 'FB7_injurySpec_pro'
  
  DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4, group.by = "subtype",
          label = TRUE, repel = TRUE)
  
  ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                '_afterFiltering_reAnnot.pdf'), 
         width = 12, height = 8)
  
  # save the annotation and processing 
  mm = match(colnames(ax), colnames(aa))
  cat(length(which(is.na(mm))), 'cells missing \n')
  aa$subtype[mm] = ax$subtype
  aa$strigentFiltered[mm] = 'keep'
  
  saveRDS(ax, file = paste0(RdataDir, 
                            'Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_', mcells, 'refinedCelltypes_20240202.rds'))
  
  saveRDS(aa, file = paste0(RdataDir, 
                            'Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_refineSubtypes_', mcells, '_20240202.rds'))
  
}


##########################################
# immune cells, here mainly for macrophage
##########################################
DimPlot(aa, reduction = 'umap', group.by = 'celltype', label = TRUE, repel = TRUE)

mm = which(aa$subtype_old == 'SMC' | aa$subtype_old == 'Pericyte'|
             aa$subtype_old == 'Gra' | aa$subtype_old == 'EPI'|
             aa$subtype_old == 'T cells' | aa$subtype_old == 'B cells'| aa$subtype_old == 'DC-like')

aa$subtype[mm] = aa$subtype_old[mm]
aa$strigentFiltered[mm] = 'keep'

refine.subtypes_adultMice_immuneCells = FALSE
if(refine.subtypes_adultMice_immuneCells){
  
  mcells = 'ImmuneCells'
  
  outDir = paste0(resDir, '/', mcells)
  if(!dir.exists(outDir)) dir.create(outDir)
  
  ax = subset(aa, cells = colnames(aa)[which(aa$subtype_old == 'Macrophage' | aa$subtype_old == 'Monocyte')])
  
  table(ax$celltype)
  table(ax$subtype_old)
  
  ax$condition = factor(ax$condition, levels = c("Sham_d1", "MI_d1", "Sham_d3", "MI_d3"))
  
  DimPlot(ax, group.by = 'subtype_old', label = TRUE, repel = TRUE, pt.size = 3)
  
  ax <- FindVariableFeatures(ax, selection.method = "vst", nfeatures = 2000)
  ax <- ScaleData(ax)
  
  ax <- RunPCA(ax, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(ax, ndims = 30)
  
  # UMAP to visualize subtypes
  ax <- RunUMAP(ax, dims = 1:10, n.neighbors = 30, min.dist = 0.1)
  
  DimPlot(ax, reduction = 'umap', group.by = 'subtype_old',label = TRUE, repel = TRUE, pt.size = 1)
  
  ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                '_subcelltype_old.pdf'), width = 10, height = 8)
  
  
  ax <- FindNeighbors(ax, dims = 1:20)
  ax <- FindClusters(ax, verbose = FALSE, algorithm = 3, resolution = 0.3)
  
  p0 = DimPlot(ax, reduction = 'umap', group.by = 'subtype_old', label = TRUE, repel = FALSE)
  p1 = DimPlot(ax, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE, repel = FALSE)
  
  p0 + p1
  
  ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, '_redo.clustering.pdf'), 
         width = 16, height = 8)
  
  DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4, group.by = "condition")
  
  ## check the quality of CM
  pdf(paste0(outDir, "/QCs_doubleCheck.pdf"), 
      width=10, height=8)
  
  DimPlot(ax, reduction = 'umap', group.by = 'subtype_old', label = T, repel = TRUE)
  DimPlot(object = ax, reduction = 'umap', raster = T, pt.size = 4, shuffle= T,label = T)
  DimPlot(ax, reduction = 'umap', group.by = 'Phase')
  DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4,group.by = "condition")
  #DimPlot(object = ax, reduction = 'umap',raster = T, shuffle= T, pt.size = 4,group.by = "group")
  #DimPlot(object = xx, reduction = 'umap',raster = T, shuffle= T, pt.size = 3,group.by = "time",
  #        cols = c("grey50",rev(viridis_pal(option = "plasma")(length(unique(xx@meta.data$time))-1))))
  FeaturePlot(ax, raster = T,pt.size = 2, features = c("nFeature_RNA","nCount_RNA","percent.mt"),
              order = T)
  
  VlnPlot(ax, features = 'nFeature_RNA')
  VlnPlot(ax, features = 'nCount_RNA')
  VlnPlot(ax, features = 'percent.mt')
  
  dev.off()
  
  Test_filtering.lowQualityCells_regressOut.factors = FALSE
  if(Test_filtering.lowQualityCells_regressOut.factors){
    
    groupby = 'seurat_clusters'
    p1 = VlnPlot(ax, features = 'nFeature_RNA', group.by = groupby) +
      geom_hline(yintercept= c(500, 800, 1000))
    p2 = VlnPlot(ax, features = 'nCount_RNA',  group.by = groupby) +
      geom_hline(yintercept= c(1000, 2000))
    p3 = VlnPlot(ax, features = 'percent.mt',  group.by = groupby) + 
      geom_hline(yintercept= c(50, 60))
    
    p1 + p2 + p3
    
    # stringent filtering for EC and discard cluster 1
    ax = subset(ax, cells = colnames(ax)[which(ax$nCount_RNA >1000 & ax$nFeature_RNA > 500)])
    
    ax <- FindVariableFeatures(ax, selection.method = "vst", nfeatures = 2000)
    ax <- ScaleData(ax)
    
    ax <- RunPCA(ax, verbose = FALSE, weight.by.var = TRUE)
    ElbowPlot(ax, ndims = 30)
    
    # UMAP to visualize subtypes
    ax <- RunUMAP(ax, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
    
    DimPlot(ax, reduction = 'umap', group.by = 'subtype_old', label = TRUE, repel = TRUE)
    
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                  '_afterFiltering_subcelltypes.pdf'), 
           width = 10, height = 8)
    
    
    ax <- FindNeighbors(ax, dims = 1:20)
    ax <- FindClusters(ax, verbose = FALSE, algorithm = 3, resolution = 0.3)
    
    p1 = DimPlot(ax, reduction = 'umap', group.by = 'seurat_clusters',  label = T, repel = TRUE)
    p0 = DimPlot(ax, reduction = 'umap', group.by = 'subtype_old',  label = T, repel = TRUE)
    p2 = DimPlot(ax, reduction = 'umap', group.by = 'Phase')
    p0 + p1 + p2
    
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                  '_afterFiltering_subcelltypes_reclustering.pdf'), 
           width = 20, height = 8)
    
    
    DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4, group.by = "seurat_clusters",
            split.by = 'condition', label = TRUE, repel = TRUE, ncol = 2)
    
    ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                  '_afterFiltering_clusterPerCondition.pdf'), 
           width = 16, height = 12)
    
  }
  
  p1 = DimPlot(ax, reduction = 'umap', group.by = 'seurat_clusters',  label = T, repel = TRUE)
  p0 = DimPlot(ax, reduction = 'umap', group.by = 'subtype_old',  label = T, repel = TRUE)
  p0 + p1
  
  ## re-annotate
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(2))))] = 'Monocyte'
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(0))))] = "Macrophage_injurySpec"
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(1))))] = 'Macrophage_d3'
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(3))))] = 'Macrophage_pro'
  ax$subtype[which(!is.na(match(ax$seurat_clusters, c(4))))] = 'Macrophage_injurySpec_d3'
  
  DimPlot(object = ax, reduction = 'umap', raster = T,shuffle= T, pt.size = 4, group.by = "subtype",
          label = TRUE, repel = TRUE)
  
  ggsave(paste0(outDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_', mcells, 
                '_afterFiltering_reAnnot.pdf'), 
         width = 12, height = 8)
  
  ax = subset(ax, cells = colnames(ax)[which(!is.na(ax$subtype))])
  
  # save the annotation and processing 
  mm = match(colnames(ax), colnames(aa))
  cat(length(which(is.na(mm))), 'cells missing \n')
  aa$subtype[mm] = ax$subtype
  aa$strigentFiltered[mm] = 'keep'
  
  mm = which(aa$subtype_old == 'Pericyte')
  aa$subtype[mm] = aa$subtype_old[mm]
  aa$strigentFiltered[mm] = 'keep'
  
  saveRDS(ax, file = paste0(RdataDir, 
                            'Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_', mcells, 
                            'refinedCelltypes_20240202.rds'))
  
  saveRDS(aa, file = paste0(RdataDir, 
                            'Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_refineSubtypes_', mcells, 
                            '_20240202.rds'))
  
}


##########################################
# update the reference 
# replaced the corrected expression matrix with Seurat_RPCA
# keep mnn dimension reduction and corrected expression matrix from Seurat_RPCA, same as adult mice 
##########################################
refs = readRDS(file = paste0(RdataDir, 
                             'Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_refineSubtypes_', 
                             "ImmuneCells", '_20240202.rds'))

DimPlot(refs, reduction = 'umap', group.by = 'celltype')
DimPlot(refs, reduction = 'umap', group.by = 'subtype')

refs = subset(refs, cells = colnames(refs)[which(!is.na(refs$subtype))])

DefaultAssay(refs) = 'mnn.reconstructed'

xx = readRDS(file = paste0('../results/scRNAseq_neonatalMouse_20231108/Rdata/', 
                           'Seurat.obj_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_dataIntegration_',
                           'Seurat_RPCA', '_allgeneCorrected.rds'))

DefaultAssay(xx) <- "integrated"

xx = subset(xx, cells = colnames(refs))
refs[['integrated']] = xx[['integrated']]

#refs@assays$mnn.reconstructed@counts = refs@assays$RNA@counts
# refs[['RNA']] = NULL ## save mnn.reconstructed 

#refs <- FindVariableFeatures(refs, selection.method = "vst", nfeatures = 5000)
#refs <- ScaleData(refs, verbose = FALSE)
#refs <- RunPCA(refs, npcs = 50, verbose = FALSE)

#ElbowPlot(refs, ndims = 30)

#refs <- FindNeighbors(refs, reduction = "pca", dims = 1:20)
#refs <- FindClusters(refs, resolution = 0.2)

#refs <- RunUMAP(refs, reduction = "pca", dims = 1:20, n.neighbors = 30, min.dist = 0.1)

DimPlot(refs, reduction = 'umap', group.by = 'celltype', raster = T,shuffle= T, pt.size = 2, 
        label = TRUE, repel = TRUE)

DimPlot(refs, reduction = 'umap', group.by = 'subtype',raster = T,shuffle= T, pt.size = 2, 
        label = TRUE, repel = TRUE)

saveRDS(refs, file = paste0('../data/data_examples/ref_scRNAseq_neonatalMice_clean.v1.2.rds'))


##########################################
# check Axl gene
##########################################
Check_Axl_inCM =FALSE
if(Check_Axl_inCM){
  aa = readRDS(file = paste0('../data/data_examples/ref_scRNAseq_neonatalMice_clean.v1.2.rds'))
  
  #DefaultAssay(cms) = 'integrated'
  p1 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)
  p2 = DimPlot(aa, group.by = 'condition',label = TRUE, repel = TRUE)
  
  p1 / p2
  p3 = FeaturePlot(aa, features = c('Axl', 'Gas6', 'Nppa', 'Myh6'), max.cutoff = 'q99')
  
  (p1 + p2) / p3
  
  ggsave(paste0(resDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_UMAP_',
                '_afterFiltering_clusterPerCondition_Gas6_Axl.pdf'), 
         width = 16, height = 16)
  
  
  p1 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)
  #p2 = DimPlot(aa, group.by = 'condition',label = TRUE, repel = TRUE)
  #p1 / p2
  #p3 = FeaturePlot(aa, features = c('Axl', 'Gas6', 'Nppa', 'Myh6'), max.cutoff = 'q99')
  #(p1 + p2) / p3
  p2 = FeaturePlot(aa, features = 'Axl') +  
    scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))
  #p2 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE)
  p1 + p2
  
  ggsave(paste0(resDir, '/Ref_neonatalMice_CM.Cui2020_noCM.Wang2020_P1_integration_UMAP_',
                '_afterFiltering_cl_Axl.pdf'), 
         width = 14, height = 6)
  
  
}
