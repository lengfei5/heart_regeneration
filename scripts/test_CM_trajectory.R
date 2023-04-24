##########################################################################
##########################################################################
# Project: heart regeration
# Script purpose: test trajectory
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Apr 18 13:07:41 2023
##########################################################################
##########################################################################
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

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
require(tictoc)
library(future)
options(future.globals.maxSize = 80 * 1024^3)
set.seed(1234)
mem_used()

outDir = paste0(resDir, '/CM_trajectory_test/')
system(paste0('mkdir -p ', outDir))

##########################################
# prepare the CM data
##########################################
aa =  readRDS(file = paste0("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/", 
                            "CM_subset_for_velocity.rds"))
aa$time = gsub('Amex_scRNA_', '', aa$condition)
aa$cell.ids = sapply(colnames(aa), function(x) unlist(strsplit(as.character(x), '-'))[1]) 
aa$cell.ids = paste0(aa$cell.ids, '_', aa$time)


levels = c("Amex_scRNA_d0", "Amex_scRNA_d1",
           "Amex_scRNA_d4", "Amex_scRNA_d7", 
           "Amex_scRNA_d14")
Idents(aa) = factor(aa$condition, levels = levels)

#DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

CMmyLevels <- c(  "CM_ven_(Robo2)", "CM_ven_(Cav3_1)", "CM_IS","CM_Prol_IS", "CM_Prol_1", "CM_Prol_3")
aa$subtypes = factor(aa$subtypes, levels = CMmyLevels)

cols = c("#4CC9F0", "#49A2F0",
         "#941F56", "#C61010",  
         "#4361EE",  "#3F37C9")

DimPlot(aa, dims = c(1, 2), label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE,
        cols = cols
)



xx = aa
xx@reductions$umap@cell.embeddings = -xx@reductions$umap@cell.embeddings
DimPlot(xx, dims = c(1, 2), label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE,
        cols = cols
)

aa = xx
rm(xx)

ggsave(filename = paste0(outDir, 'UMAP_CMsubsets_for.trajectory.test.pdf'), width = 10, height = 8)

saveRDS(aa, file = paste0(RdataDir, '"CM_subset_for_trajectory_analysis.rds"'))
########################################################
########################################################
# Section : test DM first 
# 
########################################################
########################################################
library(slingshot, quietly = FALSE)
library(destiny, quietly = TRUE)
library(mclust, quietly = TRUE)
library(scater)
library(SingleCellExperiment)
library(scran)
library(RColorBrewer)

aa = readRDS(file = paste0(RdataDir, '"CM_subset_for_trajectory_analysis.rds"'))

Remove_CM.prol.1_CM.prol.3 = TRUE
if(Remove_CM.prol.1_CM.prol.3){
  aa = subset(aa, cells = colnames(aa)[which(aa$subtypes != 'CM_Prol_1' & aa$subtypes != 'CM_Prol_3'
                                             & aa$condition != 'Amex_scRNA_d14')]) 
  
  #aa = subset(aa, cells = colnames(aa)[which(aa$subtypes != 'CM_Prol_1' & aa$subtypes != 'CM_Prol_3')]) 
  
  aa$subtypes = droplevels(aa$subtypes)
}


Test_Batch_correction_timepoint = FALSE
if(Test_Batch_correction_timepoint){
  aa.list <- SplitObject(aa, split.by = "condition")
  
  aa.list <- lapply(X = aa.list, FUN = function(x) {
    x<- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  })
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = aa.list)
  
  Use.softer.rpca = TRUE
  if(Use.softer.rpca){
    aa.list <- lapply(X = aa.list, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features, verbose = FALSE)
    })
    
    aa.anchors <- FindIntegrationAnchors(object.list = aa.list, anchor.features = features, 
                                         reduction = "rpca")
  }else{
    aa.anchors <- FindIntegrationAnchors(object.list = aa.list, anchor.features = features)
  }
  
  aa.combined <- IntegrateData(anchorset = aa.anchors)
  
  DefaultAssay(aa.combined) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  aa.combined <- ScaleData(aa.combined, verbose = FALSE)
  aa.combined <- RunPCA(aa.combined, npcs = 50, verbose = FALSE, weight.by.var = FALSE)
  
  aa.combined$condition = factor(aa.combined$condition, levels = levels(aa$condition))
  
  p1 = DimPlot(aa.combined, label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE) 
  p2 = DimPlot(aa.combined, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  p1 + p2
  
  aa.combined <- FindNeighbors(aa.combined, reduction = "pca", dims = 1:30)
  aa.combined <- FindClusters(aa.combined, resolution = 0.5)
  
  p3 = DimPlot(aa.combined, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  aa.combined$cluster = aa.combined$seurat_clusters
  
  aa.combined <- RunUMAP(aa.combined, reduction = "pca", dims = 1:30, n.neighbors = 100, min.dist = 0.3)
  #DimPlot(aa.combined, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  #p1 + p2+ p3
  DimPlot(aa.combined, label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE)
  
  saveRDS(aa.combined, file = paste0(RdataDir, 'CMsubset_batch_corrected.rds'))
  
  aa = readRDS(file = paste0(RdataDir, 'CMsubset_batch_corrected.rds'))
  
}

sce = as.SingleCellExperiment(aa)
#rm(bb)

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
metadata = aa@meta.data
dcs = data.frame(dm@eigenvectors, stringsAsFactors = FALSE)

dcs = dcs[match(rownames(metadata), cells), ]

dcs = as.matrix(dcs)

aa[["DC"]] <- CreateDimReducObject(embeddings = as.matrix(dcs), key = "DC_", assay = DefaultAssay(aa))

rm(metadata)

#aa = RunUMAP(aa, reduction = "DC", dims = 1:30, n.neighbors = 100, min.dist = 0.3, metric = 'euclidean',
#             reduction.name = "dc_umap")
#DimPlot(aa, reduction = 'dc_umap', label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE)

## try to make 3d for DC1, DC2 and DC3
## example from https://plotly-r.com/d-charts.html
library(plotly)
dcs = as.data.frame(dcs)
dcs$subtypes = aa$subtypes[match(rownames(dcs), colnames(aa))]
dcs$condition = aa$condition[match(rownames(dcs), colnames(aa))]

plot_ly(data.frame(dcs), x = ~DC1, y = ~DC2, z = ~DC4, size = 3) %>%
  add_markers(color = ~ subtypes)

plot_ly(data.frame(dcs), x = ~DC2, y = ~DC3, z = ~DC4, size = 3) %>%
  add_markers(color = ~ subtypes)

plot_ly(data.frame(dcs), x = ~DC1, y = ~DC2, z = ~DC3, size = 3) %>%
  add_markers(color = ~ condition)

p1 = DimPlot(aa, reduction = 'DC', dims = c(1, 2), cols = cols)
p2 = DimPlot(aa, reduction = 'DC', dims = c(1, 3), cols = cols)
p3 = DimPlot(aa, reduction = 'DC', dims = c(1, 4), cols = cols)
p1 / p2 / p3

DimPlot(aa, reduction = 'DC', dims = c(2, 3), label = TRUE, repel = TRUE,
        group.by = 'subtypes')

DimPlot(aa, reduction = 'DC', dims = c(1, 3), label = TRUE, repel = TRUE,
        group.by = 'subtypes')


ggsave(filename = paste0(outDir, 'CMsubsets_DM__no.day14.2D_batchCorrected.pdf'), width = 10, height = 8)

p1 = DimPlot(aa, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
        group.by = 'subtypes', cols = cols)
p2 = DimPlot(aa, reduction = 'DC', dims = c(1, 3), label = TRUE, repel = TRUE,
             group.by = 'subtypes')

DimPlot(aa, reduction = 'DC', dims = c(2, 3), label = TRUE, repel = TRUE,
        group.by = 'subtypes')

p1 / p2

saveRDS(aa, file = paste0(RdataDir, 'CMsubset_batch_corrected_DM.rds'))


##########################################
# trajectory using elastic principle graph
# example code from https://github.com/Albluca/ElPiGraph.R/blob/master/guides/base.md
##########################################
aa = readRDS(file = paste0(RdataDir, 'CMsubset_batch_corrected_DM.rds'))

dcs = as.data.frame(aa@reductions$DC@cell.embeddings)
dcs$subtypes = aa$subtypes[match(rownames(dcs), colnames(aa))]
dcs$condition = aa$condition[match(rownames(dcs), colnames(aa))]

library(plotly)
plot_ly(data.frame(dcs), x = ~DC_1, y = ~DC_2, z = ~DC_4, size = 3) %>%
  add_markers(color = ~ subtypes)

plot_ly(data.frame(dcs), x = ~DC_2, y = ~DC_3, z = ~DC_4, size = 3) %>%
  add_markers(color = ~ subtypes)

plot_ly(data.frame(dcs), x = ~DC_1, y = ~DC_2, z = ~DC_3, size = 3) %>%
  add_markers(color = ~ condition)


library("ElPiGraph.R")

#CurveEPG <- computeElasticPrincipalCurve(X = curve_data, NumNodes = 10)
TreeEPG <- computeElasticPrincipalTree(X = as.matrix(dcs[, c(1:5)]), 
                                       NumNodes = 60, 
                                       Lambda = .03, Mu = .01, 
                                       Do_PCA = FALSE,
                                       ShowTimer = FALSE, ReduceDimension = NULL)

graphEPG = computeElasticPrincipalGraph(Data = as.matrix(dcs[, c(1:5)]), 
                                       NumNodes = 60, 
                                       Lambda = .03, 
                                       Mu = .01, 
                                       Do_PCA = FALSE,
                                       ShowTimer = FALSE, ReduceDimension = NULL)

if(trajectory_pseudotime_method == 'slingshot'){ #
  rd =  aa[['DC']]@cell.embeddings[, c(1, 2)]
  cl = aa$subtypes
  
  # reclustering 
  set.seed(2022)
  res_kmean = kmeans(rd, centers = 20)
  cl2 <- res_kmean$cluster
  
  centers <- res_kmean$centers[res_kmean$cluster, ] 
  distances <- sqrt(rowSums((rd - centers)^2))
  
  hist(distances, breaks = 100)
  abline(v = 0.005, col = 'red')
  outliers <- which(distances>0.005)
  cat(length(outliers), ' outliers found \n')
  
  cols = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[res_kmean$cluster]
  
  plot(rd, pch=19, col=cols, cex=0.2)
  points(res_kmean$centers, col='darkred', pch=15, cex=1)
  points(rd[outliers, ], pch="+", col = 'darkgray', cex=0.8)
  
  colData(sce)$kmeans <- cl2
  # 
  aa$dc_clusters = cl2
  #aa$dc_clusters[which(aa$celltypes == '14'| aa$celltypes == '11')] = NA
  DimPlot(aa, reduction = 'DC', dims = c(1, 2), group.by = 'dc_clusters', label = TRUE, repel = TRUE)
  # ggsave(filename = paste0(outDir, 'DMcomponents_clustering_for_trajectory_slingshot_outliers.pdf'), 
  #        width = 16, height = 8)
  
  table(aa$condition, aa$dc_clusters)
  
  ### manually change the cluster labels
  cl2 = aa$dc_clusters
  cl2[which(cl2== 12 | cl2 == 8 | cl2==13)] = NA # day4_RA
  #cl2[which(cl2 == 10 | cl2 == 2)] = 10 # day4_noRA 
  colData(sce)$kmeans <- cl2
  
  aa$dc_clusters = cl2
  #aa$dc_clusters[!is.na(aa$dc_clusters_outliers)] = NA
  DimPlot(aa, reduction = 'DC', dims = c(1, 2), group.by = 'dc_clusters', label = TRUE, repel = TRUE)
  
  ggsave(filename = paste0(outDir, 
                           'DMcomponents_clustering.manaulCorrection_for_trajectory_slingshot_outliers.pdf'), 
         width = 16, height = 8)
  

  table(aa$condition, aa$dc_clusters)
  
  
  #kk_clean = which(!is.na(aa$dc_clusters))
  #rd3 = rd2[-outliers, ]
  #cl3 = cl2[-outliers]
  
  rd =  aa[['DC']]@cell.embeddings[, c(1,2)]
  cl2 = aa$dc_clusters
  
  index_keep = which(!is.na(cl2))
  rd = rd[index_keep, ]
  cl2 = cl2[index_keep]
  
  # reclustering 
  set.seed(2022)
  res_kmean = kmeans(rd, centers = 6)
  cl2 <- res_kmean$cluster
  cols = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[res_kmean$cluster]
  
  plot(rd, pch=19, col=cols, cex=0.2)
  points(res_kmean$centers, col='darkred', pch=15, cex=1)
  
  #cl2 = aa$dc_clusters
  
  #colData(sce)$kmeans <- cl2
  aa$dc_clusters = NA
  aa$dc_clusters = cl2[match(names(cl2), colnames(aa))]
  DimPlot(aa, reduction = 'DC', dims = c(1, 2), group.by = 'dc_clusters', label = TRUE, repel = TRUE)
  
  lin1 <- getLineages(rd, cl2, start.clus = '4', end.clus = c('1') )
  
  #cl.uniq = unique(cl)
  cols = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[cl2]
  #names(cols) = cl.uniq
  
  pdf(paste0(outDir, "Slingshot_getLineage_withStarting.cluster.pdf"),
      height = 6, width =10, useDingbats = FALSE)
  
  plot(rd, col = cols, asp = 1, pch = 1, cex = 0.5)
  lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')
  
  dev.off()
  
  crv1 <- getCurves(lin1, shrink = 1, stretch = 2,  smoother = "smooth.spline", thresh = 0.1, maxit = 50, 
                    shrink.method = 'density')
  # crv1 <- getCurves(lin1, shrink = TRUE, stretch = 2,  smoother = "loess",
  #                   thresh = 0.001, maxit = 50, 
  #                   shrink.method = 'density')
  
  pdf(paste0(outDir, "Slingshot_getCurve_withLineage.pdf"),
      height = 6, width =10, useDingbats = FALSE)
  
  
  #crv1
  plot(rd, col = cols, asp = 1, pch = 16, cex = 0.6)
  lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')
  
  dev.off()
  
  
  scrv = slingCurves(crv1)
  pst = as.data.frame(slingPseudotime(crv1))
  
  pst$Lineage1 = pseudotime.scaling(pst$Lineage1)
  pst$weight.lineage1 = scrv$Lineage1$w
  pst$weight.lineage2 = scrv$Lineage2$w
  
  aa$dpt = NA
  aa$dpt = pst$Lineage1[match(colnames(aa), rownames(pst))]
  
  DimPlot(aa, reduction = 'DC', dims = c(1, 2), group.by = 'subtypes', label = TRUE, repel = TRUE)
  FeaturePlot(aa, features = 'dpt', reduction = 'DC', dims = c(1:2), coord.fixed = TRUE) +
    xlim(c(-0.06, 0.02)) +ylim(c(-0.025,0.04))
  
  ggsave(filename = paste0(outDir, 
                           'CMsubset_pseudtotime_dpt.pdf'),  width = 10, height = 6)
  
  save(aa, sce, lin1, crv1, 
       file = paste0(outDir, 
                     'CMsubset_noDay14_DM_slingshot_lineage_pseudottime.Rdata'))
  
  
  
}

##########################################
# not used here 
##########################################
Trajectory_pseudotime_principlecurve = FALSE
if(Trajectory_pseudotime_principlecurve){
  library(destiny)
  library(princurve)
  
  pseudotime.scaling = function(X) {
    return((X - min(X))/diff(range(X)))
  }
  
  # save the pseudotime estimation 
  #pdfname = paste0(resDir, "/pseudotime_estimation_v1.pdf")
  #pdf(pdfname, width=12, height = 10)
  #par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  aa$pseudot_noRA = NA
  aa$pseudot_RA = NA 
  
  DimPlot(aa, reduction = 'DC', group.by = 'subtypes', label = TRUE, repel = TRUE)
  table(aa$condition, aa$dc_clusters)
  
  
  dcs_all = data.frame(aa[['DC']]@cell.embeddings[, c(1,2)])
  
  ## no RA conditions
  #cluster_sels = c(12, 4, 8, 7,9,10)
  mm = which(!is.na(aa$dc_clusters))
  dcs = dcs_all[mm, ]
  
  cols = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[aa$dc_clusters[mm]]
  plot(dcs, col = cols, cex = 0.2)
  dcs[, 1] = - dcs[,1]
  
  princurve = principal_curve(x = as.matrix(dcs[, c(1:2)]), 
                              #start = start,
                              smoother = 'smooth_spline', stretch = 2)
  
  cols = c(brewer.pal(9,"Set1"), brewer.pal(8,"Set2"))[aa$dc_clusters]
  plot(dcs[, c(1:2)], cex = 0.1)
  #points(start, cex = 0.4, col = 'green')
  lines(princurve$s[order(princurve$lambda),], lty=1,lwd=4,col="red",type = "l")
  
  
}

