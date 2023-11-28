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

version.analysis = '_20231108'

resDir = paste0("../results/scRNAseq_neonatalMouse", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

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
options(future.globals.maxSize = 160000 * 1024^2)
library(pryr) # monitor the memory usage

mem_used()

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
  aa = readRDS(file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM_Cui2020_shoval.version_P1_clusterd.rds'))
  
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
  DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition', split.by = 'condition')
  
  aa$condition = factor(aa$condition, levels = c('P1_Sham_d3', 'P1_Sham_d1', 'P1_MI_d1', 'P1_MI_d3'))
  
  Test_DataIntegration = FALSE
  if(Test_DataIntegration){
    
    source('functions_dataIntegration.R')
    ref.combined = IntegrateData_Seurat_CCA(aa, 
                                            group.by = 'condition', 
                                            merge.order = matrix(c(-2, -3, 1, -1, -4, 2), ncol = 2),
                                            correct.all = TRUE)
    DimPlot(ref.combined, group.by = 'condition') + ggtitle('Seurat_CCA')
    
    ggsave(filename = paste0(resDir, '/CM_Cui2020_P1_dataIntegration_CCA.pdf'), 
           width = 10, height = 6)
    
    DimPlot(ref.combined, group.by = 'condition', split.by = 'condition')
    ggsave(filename = paste0(resDir, '/CM_Cui2020_P1_dataIntegration_CCA_perCondition.pdf'), 
           width = 10, height = 6)
    
    
    source('functions_dataIntegration.R')
    ref.combined = IntegrateData_Seurat_RPCA(aa, 
                                            group.by = 'condition', 
                                            #merge.order = matrix(c(-1, -3, 1, -2, -4, 2), ncol = 2),
                                            correct.all = TRUE)
    DimPlot(ref.combined, group.by = 'condition') + ggtitle('Seurat_RPCA')
    
    ggsave(filename = paste0(resDir, '/CM_Cui2020_P1_dataIntegration_RPCA.pdf'), 
           width = 10, height = 6)
    
    DimPlot(ref.combined, group.by = 'condition', split.by = 'condition')
    ggsave(filename = paste0(resDir, '/CM_Cui2020_P1_dataIntegration_RPCA_perCondition.pdf'), 
           width = 10, height = 6)
    
    source('functions_dataIntegration.R')
    refs.merged = IntegrateData_runHarmony(aa, group.by = 'condition', max.iter.harmony = 10)
    
    DimPlot(refs.merged, group.by = 'condition')
    
    DimPlot(refs.merged, group.by = 'condition', split.by = 'condition') + ggtitle('Harmony')
    ggsave(filename = paste0(resDir, '/CM_Cui2020_P1_dataIntegration_Harmony.pdf'), 
           width = 10, height = 6)
    
    DimPlot(refs.merged, group.by = 'condition', split.by = 'condition')
    
    p1 = DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE, group.by = 'condition')
    p2 = DimPlot(refs.merged, group.by = 'condition') + ggtitle('Harmony')
    
    p1 /p2
    
    ggsave(filename = paste0(resDir, '/CM_Cui2020_P1_dataIntegration_noIntegratoin.vs.Harmony.pdf'), 
           width = 10, height = 12)
    
    
    ### fastMNN is the most satisfying of integration 
    source('functions_dataIntegration.R')
    refs.merged = IntegrateData_runFastMNN(aa, group.by = 'condition',
                                           merge.order = c('P1_Sham_d3', 'P1_Sham_d1', 'P1_MI_d1', 'P1_MI_d3'), 
                                           correct.all = TRUE)
    
    DimPlot(refs.merged, group.by = 'condition') + ggtitle('fastMNN')
    ggsave(filename = paste0(resDir, '/CM_Cui2020_P1_dataIntegration_fastMNN.pdf'), 
           width = 10, height = 6)
    
    DimPlot(refs.merged, group.by = 'condition', split.by = 'condition')
    
    p1 = DimPlot(aa, label = TRUE, repel = TRUE, raster=FALSE)
    p2 = DimPlot(refs.merged) + ggtitle('fastMNN')
    
    p1 /p2
    
    ggsave(filename = paste0(resDir, '/CM_Cui2020_P1_dataIntegration_noIntegration.vs.fastMNN.pdf'), 
           width = 10, height = 12)
    
    refs.merged <- FindNeighbors(refs.merged, reduction = "mnn", dims = 1:20)
    refs.merged <- FindClusters(refs.merged, resolution = 0.4)
    
    DimPlot(refs.merged, split.by = 'condition') + ggtitle('fastMNN')
    ggsave(filename = paste0(resDir, '/CM_Cui2020_P1_dataIntegration_fastMNN_clusters.pdf'), 
           width = 16, height = 6)
    
    saveRDS(refs.merged, 
            file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM_Cui2020_shoval.version_P1_clusterd_fastMNN.rds'))
    
    FeaturePlot(refs.merged, features = c('Nppa', 'Ankrd1'))
    ggsave(filename = paste0(resDir, '/CM_Cui2020_P1_dataIntegration_fastMNN_Nppa.ANkrd1.pdf'), 
           width = 10, height = 6)
    
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
  aa = readRDS(file = paste0(RdataDir, 'Wang_DevCell_2020_scRNAseq_P1_regressedout.rds'))
  
  aa$condition = factor(aa$condition, levels = c('P1_Sham_D1', 'P1_Sham_D3',  'P1_MI_D1', 'P1_MI_D3'))
  
  Test_DataIntegration = FALSE
  if(Test_DataIntegration){
    ## no data integration 
    p1 = DimPlot(aa, group.by = 'FineID', label = TRUE, repel = TRUE, raster=FALSE)
    p2 = DimPlot(aa, group.by = 'condition', label = TRUE, repel = TRUE)
    
    p1 + p2
    
    ggsave(filename = paste0(resDir, '/Wang2020_P1_celltypes_conditions_noDataIntegration.pdf'), 
           width = 14, height = 6)
    
    DimPlot(aa, group.by = 'FineID', split.by = 'condition')
    ggsave(filename = paste0(resDir, '/noCM_Wang2020_P1_celltypes_noDataIntegration_perCondition.pdf'), 
           width = 16, height = 6)
    
    
    ## Seurat CCA
    source('functions_dataIntegration.R')
    ref.combined = IntegrateData_Seurat_CCA(aa, 
                                            group.by = 'condition', 
                                            merge.order = matrix(c(-1, -3, 1, -2, -4, 2), ncol = 2),
                                            correct.all = TRUE)
    
    p1 = DimPlot(ref.combined, group.by = 'FineID', label = TRUE, repel = TRUE, raster=FALSE) + 
      ggtitle('Seurat_CCA')
    p2 = DimPlot(ref.combined, group.by = 'condition', label = TRUE, repel = TRUE) +
      ggtitle('Seurat_CCA')
    p1 + p2
    
    ggsave(filename = paste0(resDir, '/noCM_Wang2020_P1_dataIntegration_CCA.pdf'), 
           width = 16, height = 6)
    
    DimPlot(ref.combined, group.by = 'FineID', split.by = 'condition')
    ggsave(filename = paste0(resDir, '/noCM_Wang2020_P1_dataIntegration_CCA_perCondition.pdf'), 
           width = 16, height = 6)
    
    saveRDS(ref.combined, 
            file = paste0(RdataDir, 'Seurat.obj_neonatalMice_noCM_Wang2020_P1_SeuratCCA.rds'))
    
    #ref.combined = readRDS(file = paste0(RdataDir, 'Seurat.obj_neonatalMice_noCM_Wang2020_P1_SeuratCCA.rds'))
    
    ## Seurat RPCA
    source('functions_dataIntegration.R')
    ref.combined = IntegrateData_Seurat_RPCA(aa, 
                                             group.by = 'condition', 
                                             #merge.order = matrix(c(-1, -3, 1, -2, -4, 2), ncol = 2),
                                             correct.all = TRUE)
    
    p1 = DimPlot(ref.combined, group.by = 'FineID', label = TRUE, repel = TRUE, raster=FALSE) + 
      ggtitle('Seurat_RPCA')
    p2 = DimPlot(ref.combined, group.by = 'condition', label = TRUE, repel = TRUE) + 
      ggtitle('Seurat_RPCA')
    
    p1 + p2
    
    ggsave(filename = paste0(resDir, '/noCM_Wang2020_P1_dataIntegration_RPCA.pdf'), 
           width = 16, height = 6)
    
    DimPlot(ref.combined, group.by = 'FineID', split.by = 'condition')
    ggsave(filename = paste0(resDir, '/noCM_Wang2020_P1_dataIntegration_RPCA_perCondition.pdf'), 
           width = 16, height = 6)
    
    saveRDS(ref.combined,
            file = paste0(RdataDir, 'Seurat.obj_neonatalMice_noCM_Wang2020_P1_SeuratRPCA.rds'))
    
    ## Harmony
    source('functions_dataIntegration.R')
    refs.merged = IntegrateData_runHarmony(aa, group.by = 'condition', max.iter.harmony = 10)
    
    p1 = DimPlot(refs.merged, group.by = 'FineID', label = TRUE, repel = TRUE, raster=FALSE) + 
      ggtitle('Harmony')
    p2 = DimPlot(refs.merged, group.by = 'condition', label = TRUE, repel = TRUE) + 
      ggtitle('Harmony')
    
    p1 + p2
    
    ggsave(filename = paste0(resDir, '/noCM_Wang2020_P1_dataIntegration_Harmony.pdf'), 
           width = 16, height = 6)
    
    DimPlot(ref.combined, group.by = 'FineID', split.by = 'condition')
    ggsave(filename = paste0(resDir, '/noCM_Wang2020_P1_dataIntegration_Harmony_perCondition.pdf'), 
           width = 16, height = 6)
    
    saveRDS(ref.combined,
            file = paste0(RdataDir, 'Seurat.obj_neonatalMice_noCM_Wang2020_P1_Harmony.rds'))
    
    
    ### fastMNN is the most satisfying of integration 
    source('functions_dataIntegration.R')
    refs.merged = IntegrateData_runFastMNN(aa, group.by = 'condition',
                                           merge.order = c('P1_Sham_D1', 'P1_Sham_D3',  'P1_MI_D1', 'P1_MI_D3'), 
                                           correct.all = TRUE)
    
    p1 = DimPlot(refs.merged, group.by = 'FineID', label = TRUE, repel = TRUE, raster=FALSE) + 
      ggtitle('fastMNN')
    p2 = DimPlot(refs.merged, group.by = 'condition', label = TRUE, repel = TRUE) + 
      ggtitle('fastMNN')
    
    p1 + p2
    
    ggsave(filename = paste0(resDir, '/noCM_Wang2020_P1_dataIntegration_fastMNN.pdf'), 
           width = 16, height = 6)
    
    DimPlot(ref.combined, group.by = 'FineID', split.by = 'condition')
    ggsave(filename = paste0(resDir, '/noCM_Wang2020_P1_dataIntegration_fastMNN_perCondition.pdf'), 
           width = 16, height = 6)
    
    saveRDS(ref.combined,
            file = paste0(RdataDir, 'Seurat.obj_neonatalMice_noCM_Wang2020_P1_fastMNN.rds'))
    
    
    FeaturePlot(refs.merged, features = c('Nppa', 'Ankrd1'))
    ggsave(filename = paste0(resDir, '/noCM_Wang2020_P1_dataIntegration_fastMNN_Nppa.ANkrd1.pdf'), 
           width = 10, height = 6)
    
  }
  
  
  
  
  
}






