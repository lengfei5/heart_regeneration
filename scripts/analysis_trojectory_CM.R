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
Prepare_CMsubset = FALSE
if(Prepare_CMsubset){
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
  
}

########################################################
########################################################
# Section II: test DM first 
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

DimPlot(aa, group.by = 'subtypes', label = TRUE, repel = TRUE)


aa$pst = NA

Remove_CM.prol.1_CM.prol.3 = FALSE
if(Remove_CM.prol.1_CM.prol.3){
  
  Idents(aa) = as.factor(aa$subtypes)
  sub_obj = subset(aa, idents = c('CM_ven_(Cav3_1)', 'CM_IS', 'CM_Prol_IS'))
  #aa = subset(aa, cells = colnames(aa)[which(aa$subtypes != 'CM_Prol_1' & aa$subtypes != 'CM_Prol_3')]) 
  
  #aa = subset(aa, cells = colnames(aa)[which(aa$subtypes != 'CM_Prol_1' & aa$subtypes != 'CM_Prol_3')]) 
  sub_obj$subtypes = droplevels(sub_obj$subtypes)
  
  table(sub_obj$subtypes)
  
}

sce = as.SingleCellExperiment(sub_obj)
#rm(bb)

dec <- modelGeneVar(sce)

#plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
#curve(metadata(dec)$trend(x), col="blue", add=TRUE)

nb_features = 1000; 
top.hvgs <- getTopHVGs(dec, n=nb_features)

sce <- runPCA(sce, subset_row=top.hvgs, ncomponents = 100)
# reducedDimNames(sce)
ll.pca = reducedDim(sce, 'PCA')[, c(1:50)]

n_neighbors = 50; n_eigs = 50; sigma = 'local';

tic()

dm <- DiffusionMap(ll.pca, sigma = sigma, 
                   #k = n_neighbors, 
                   n_eigs = n_eigs, distance = 'euclidean')

toc()

cells = names(dm$DC1)

metadata = sub_obj@meta.data
dcs = data.frame(dm@eigenvectors, stringsAsFactors = FALSE)

dcs = dcs[match(rownames(metadata), cells), ]

dcs = as.matrix(dcs)

sub_obj[["DC"]] <- CreateDimReducObject(embeddings = as.matrix(dcs), key = "DC_", assay = DefaultAssay(sub_obj))

rm(metadata)

DimPlot(sub_obj, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
        group.by = 'subtypes')

DimPlot(sub_obj, reduction = 'DC', dims = c(1, 3), label = TRUE, repel = TRUE,
        group.by = 'subtypes')

DimPlot(sub_obj, reduction = 'DC', dims = c(2, 3), label = TRUE, repel = TRUE,
        group.by = 'subtypes')
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
  add_markers(color = ~ subtypes)

p1 = DimPlot(aa, reduction = 'DC', dims = c(1, 2), cols = cols)
p2 = DimPlot(aa, reduction = 'DC', dims = c(1, 3), cols = cols)
p3 = DimPlot(aa, reduction = 'DC', dims = c(1, 4), cols = cols)
p1 / p2 / p3

DimPlot(sub_obj, reduction = 'DC', dims = c(2, 3), label = TRUE, repel = TRUE,
        group.by = 'subtypes')

DimPlot(sub_obj, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
        group.by = 'subtypes')


p1 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
        group.by = 'subtypes')
p2 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 3), label = TRUE, repel = TRUE,
             group.by = 'subtypes')

p1 / p2

ggsave(filename = paste0(outDir, 'trajectory_test_DM_CM_subtypes.no.prol1.prol3_timepoints.all.pdf'), 
       width = 8, height = 16)

DimPlot(aa, reduction = 'DC', dims = c(2, 3), label = TRUE, repel = TRUE,
        group.by = 'subtypes')

save(aa, dcs, file = paste0(RdataDir, 
                            'trajectory_test_DM_CM_subtypes.no.prol1.prol3_timepoints.all.Rdata'))


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
# trajectory using elastic principle graph
# example code from https://github.com/Albluca/ElPiGraph.R/blob/master/guides/base.md
# not used here
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

##########################################
# test batch correction across time points
##########################################
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


########################################################
########################################################
# Section II: pseudotime with scanpy for injury-specific and injury-nospecific
# 
########################################################
########################################################
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
# import the multiome data
##########################################
aa = readRDS(file = paste0('../results/sc_multiome_R13591_atac_reseq_20221115/Rdata/',
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

motif_tf = readRDS(file = paste0('../results/sc_multiome_R13591_atac_reseq_20221115/Rdata/', 
                                 'motif_to_tfs_pfm_JASPAR2020_CORE_vertebrate_v1.rds'))
chromvar = readRDS(file = paste0('../results/sc_multiome_R13591_atac_reseq_20221115/Rdata/', 
                                 'atac_seuratObject_motifClass_chromVAR_v3.rds'))
DefaultAssay(chromvar) <- 'chromvar'
ss = colSums(chromvar@assays$chromvar@data)
length(which(is.na(ss)))
data = chromvar@assays$chromvar@data
data[which(is.na(data))] = 0
chromvar@assays$chromvar@data = data

##########################################
# subset CM sels 
##########################################
celltype_sel = 'CM'

sub_obj = subset(aa, cells = colnames(aa)[which(aa$celltypes == celltype_sel)])
sub_obj$subtypes = droplevels(sub_obj$subtypes)

## make sure the CM subtypes are the latest version
ref =  readRDS(file = paste0("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/", 
                             "CM_subset_for_velocity.rds"))

mm = match(colnames(sub_obj), colnames(ref))

cat(length(which(is.na(mm))), ' cells without updated subtypes \n')
cat(length(which(is.na(match(colnames(ref), colnames(sub_obj))))), 
    ' cells filtered in the scATAC analysis\n')

cells_keep = colnames(sub_obj)[which(!is.na(mm))]
mm = mm[which(!is.na(mm))]

sub_obj = subset(sub_obj, cells = cells_keep)
sub_obj$subtypes = ref$subtypes[mm]

sub_obj$subtypes = droplevels(sub_obj$subtypes)

Idents(sub_obj) <- factor(sub_obj$subtypes)

DimPlot(ref, group.by = 'subtypes')

umap.embedding = ref@reductions$umap@cell.embeddings
umap.embedding = umap.embedding[match(colnames(sub_obj), rownames(umap.embedding)), ]

sub_obj[['umap']] = Seurat::CreateDimReducObject(embeddings=umap.embedding,
                                                 key='UMAP_',
                                                 assay='RNA')
rm(umap.embedding)
DefaultAssay(sub_obj) <- 'RNA'

subtype_version = '_injurySubtypes'

sel_subtypes = c('CM_IS', 'CM_ven_(Cav3_1)', 'CM_Prol_IS')
# sub_obj = subset(sub_obj, cells = colnames(sub_obj)[which(!is.na(match(sub_obj$subtypes, sel_subtypes)))])

sub_obj = subset(sub_obj, cells = colnames(sub_obj)[which(!is.na(match(sub_obj$subtypes, sel_subtypes)))])
sub_obj$subtypes = droplevels(sub_obj$subtypes)

table(sub_obj$subtypes)

refine_subpopulation_for_trajectory = FALSE
if(refine_subpopulation_for_trajectory){
  
  #subtypes_sel = c("EC", "EC_IS_(IARS1)", "EC_IS_(LOX)", "EC_IS_Prol")
  subtypes_sel = setdiff(unique(sub_obj$subtypes), c("EC_IS_(IARS1)", "EC_IS_(LOX)", "EC_IS_Prol"))
  
  sub_obj = subset(sub_obj,  cells = colnames(sub_obj)[which(!is.na(match(sub_obj$subtypes, subtypes_sel)))])
  
  sub_obj$subtypes = droplevels(sub_obj$subtypes)
  
  ## redo the clustering in case needed in the downstream analysis
  sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000)
  sub_obj <- ScaleData(sub_obj)
  sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), weight.by.var = TRUE, 
                    verbose = FALSE)
  ElbowPlot(sub_obj, ndims = 50)
  
  sub_obj <- RunUMAP(sub_obj, dims = 1:10, n.neighbors = 30, min.dist = 0.1)
  
  DimPlot(sub_obj, group.by = 'subtypes', label = TRUE, repel = TRUE)
  
}else{
  ## redo the clustering in case needed in the downstream analysis
  sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 1000)
  sub_obj <- ScaleData(sub_obj)
  sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), weight.by.var = TRUE, 
                    verbose = FALSE)
  ElbowPlot(sub_obj, ndims = 50)
  
  sub_obj <- RunUMAP(sub_obj, dims = 1:10, n.neighbors = 30, min.dist = 0.3)
  DimPlot(sub_obj, group.by = 'subtypes', label = TRUE, repel = TRUE)
  
}

sub_obj <- FindNeighbors(sub_obj, dims = 1:10)
sub_obj <- FindClusters(sub_obj, verbose = FALSE, algorithm = 3, resolution = 0.3)

p1 = DimPlot(sub_obj, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()
p2 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)

p1 + p2

sub_obj$time = gsub('d', '', sub_obj$time)
sub_obj$clusters = sub_obj$seurat_clusters

ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_subset_', celltype_sel, subtype_version,
                         '_reclustered.pdf'), 
       height = 5, width = 14)

# DefaultAssay(sub_obj) <- 'RNA'
# p1 = DimPlot(sub_obj, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()
# 
# DefaultAssay(sub_obj) <- 'ATAC'
# p2 = DimPlot(sub_obj, label = TRUE, reduction = 'umap_topics',
#              group.by = 'subtypes',  repel = TRUE) + NoLegend()
# p1 + p2
# 
# ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_subset_', celltype_sel,  subtype_version, '.pdf'), 
#        height = 6, width = 14)
# 
# DefaultAssay(sub_obj) <- 'RNA'
# DimPlot(sub_obj, label = TRUE, group.by = 'subtypes', split.by = 'condition', repel = TRUE) + NoLegend()
# 
# ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_subset_', celltype_sel, subtype_version, 
#                          '_bytimePoint.pdf'), 
#        height = 5, width = 20)


##########################################
# preapre the spliced and unspliced matrix   
##########################################
DefaultAssay(sub_obj) = 'RNA'

#saveRDS(sub_obj, file = paste0(outDir, 'seuratObj_multiome_snRNA_scATAC_',  celltype_sel, 
#                               subtype_version, '.rds'))

source('utility_velocity.R')
mnt = preapre_dataFile_for_RNAvelocity_PAGA(seuratObj = sub_obj)

Idents(mnt) = mnt$condition
table(mnt$condition)

#mnt = subset(mnt, downsample = 2000)

saveFile = paste0('RNAmatrix_umap_kalisto.velocity_spliced_unspliced_',
                  'CM_subtypes', subtype_version, '_timepoints.all_downsample.h5Seurat')

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




