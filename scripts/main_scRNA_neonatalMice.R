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
# Section II : process the data shared by Shoval
# Wang et al., 2020 Cell Reports and 
# Cui et al., Dev Cell 2020
########################################################
########################################################
dataDir = "/groups/tanaka/Collaborations/Jingkui-Elad/Mouse_data_shoval/"

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
  
  
  saveRDS(aa, file = paste0(RdataDir, 'Seurat.obj_neonatalMice_CM_Cui2020_furtherCleaning_rmDoublets.rds'))

  
}





xx = readRDS(file = paste0(dataDir, '/seurObj_CM_norm_PCA_UMAP_jan21.rds')) # with 18818 cells
#xx = readRDS(file = paste0(dataDir, '/CMs.rds')) # with 16934 cells
identical(aa, xx)

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
p2 = DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'orig.ident', raster=FALSE)

p1 + p2

mm = match(colnames(xx), colnames(aa))

p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
p2 = DimPlot(xx, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)

p1 + p2


##########################################
# Import the scRNA-seq from shoval
# original data from 
##########################################
process_otherCelltypes_Wang.et.al.2020 = FALSE
if(process_otherCelltypes_Wang.et.al.2020){
  
 
  
  aa[["percent.mt"]] <- PercentageFeatureSet(aa, pattern = "^mt-")
  #plot1 <- FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "percent.mt")
  #plot2 <- FeatureScatter(aa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #plot1 + plot2
  
  DefaultAssay(aa) = 'RNA'
  Idents(aa) = factor(aa$condition)
  
  p1 = VlnPlot(aa, features = 'nFeature_RNA', y.max = 10000) +
    geom_hline(yintercept = c(200, 500, 1000)) + NoLegend()
  p2 = VlnPlot(aa, features = 'nCount_RNA', y.max = 100000) + NoLegend()
  #p3 = VlnPlot(aa, features = 'percent.mt', y.max = 100) + NoLegend()
  
  p1 | p2 | p3
  
  aa = subset(aa, subset = nFeature_RNA > 200  & nCount_RNA < 25000 &  percent.mt < 25)
  
  
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
    
    DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'FineID', raster=FALSE)
    
    saveRDS(aa, file =  paste0(RdataDir, 'Seurat.obj_neonatalMice_Normalization_umap_CM_Cui2020.rds'))
    
  }
  
}


########################################################
########################################################
# Section II : 
# 
########################################################
########################################################


