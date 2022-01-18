##########################################################################
##########################################################################
# Project: Elad's heart regeneration 
# Script purpose: the common functions used for adult, neonadal mice and axolotl
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Oct  1 10:37:45 2021
##########################################################################
##########################################################################
require(Seurat)
require(SeuratObject)
require(ggplot2)
require(tibble)
require(dplyr)
require(tictoc)
library(patchwork)

firstup <- function(x) {
  x = tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


########################################################
########################################################
# Section : process the nf output for axolotl 
# original codes from nf pipeline from Ashley and Tomas
########################################################
########################################################
make_SeuratObj_visium = function(topdir = './', saveDir = './results')
{
  #topdir = "./" # source dir
  library(Seurat)
  library(DropletUtils)
  library(edgeR)
  
  if(!dir.exists(saveDir)) dir.create(saveDir)
  
  ##########################################
  # missed function in DropletUtils
  ##########################################
  
  estimateAmbience <- function(m, lower=100, by.rank=NULL, round=TRUE, good.turing=TRUE) {
    m <- .rounded_to_integer(m, round)
    totals <- .intColSums(m)
    lower <- .get_lower(totals, lower, by.rank=by.rank)
    
    if (good.turing) {
      a <- .compute_ambient_stats(m, totals, lower=lower)
      output <- numeric(nrow(m))
      output[!a$discard] <- a$ambient.prop
      names(output) <- rownames(m)
    } else {
      ambient <- totals <= lower
      output <- rowSums(m[,ambient,drop=FALSE])
    }
    
    output
  }
  
  .rounded_to_integer <- function(m, round=TRUE) {
    if (round) {
      cs <- colSums(m)
      rs <- rowSums(m)
      if (!isTRUE(all.equal(cs, round(cs))) ||
          !isTRUE(all.equal(rs, round(rs))))
      {
        m <- round(m)
      }
    }
    m
  }
  
  .intColSums <- function(m) {
    # Enforcing discreteness mainly for emptyDrops()'s Monte Carlo step.
    as.integer(round(colSums(m)))
  }
  
  #' @importFrom Matrix rowSums
  .compute_ambient_stats <- function(m, totals, lower) {
    # This doesn't invalidate 'totals', by definition.
    discard <- rowSums(m) == 0
    if (any(discard)) {
      m <- m[!discard,,drop=FALSE]
    }
    
    # Computing the average profile from the ambient cells.
    ambient <- totals <= lower # lower => "T" in the text.
    ambient.m <- m[,ambient,drop=FALSE]
    ambient.prof <- rowSums(ambient.m)
    
    if (sum(ambient.prof)==0) {
      stop("no counts available to estimate the ambient profile")
    }
    ambient.prop <- .safe_good_turing(ambient.prof)
    
    list(
      m=m, # this MUST have the same number of columns as input.
      discard=discard,
      ambient=ambient,
      ambient.m=ambient.m,
      ambient.prop=ambient.prop
    )
  }
  
  .get_lower <- function(totals, lower, by.rank) {
    if (is.null(by.rank)) {
      lower
    } else if (by.rank >= length(totals)) {
      stop("not have enough columns for supplied 'by.rank'")
    } else {
      totals[order(totals, decreasing=TRUE)[by.rank+1]]
    }
  }
  
  #' @importFrom edgeR goodTuringProportions
  .safe_good_turing <- function(ambient.prof) {
    ambient.prob <- goodTuringProportions(ambient.prof)
    
    still.zero <- ambient.prob<=0
    if (any(still.zero)) {
      pseudo.prob <- 1/sum(ambient.prof)
      ambient.prob[still.zero] <- pseudo.prob/sum(still.zero)
      ambient.prob[!still.zero] <- ambient.prob[!still.zero] * (1 - pseudo.prob)
    }
    
    ambient.prob
  }
  
  
  # import read in from kallisto output
  exp = Matrix::readMM(paste0(topdir, "genecounts.mtx")) #read matrix
  bc = read.csv(paste0(topdir, "genecounts.barcodes.txt"), header = F, stringsAsFactors = F)
  g = read.csv(paste0(topdir, "genecounts.genes.txt"), header = F, stringsAsFactors = F)
  dimnames(exp) = list(paste0(bc$V1,"-1"),g$V1) # number added because of seurat format for barcodes
  count.data = Matrix::t(exp)
  
  # plot rankings for number of UMI
  br.out <- barcodeRanks(count.data)
  
  pdf("UMIrank.pdf", height = 5, width = 5, useDingbats = F)
  
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
  
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
         legend=c("knee", "inflection"))
  
  dev.off()
  
  # UMI duplication
  umi = read.table(paste0(topdir, "umicount.txt"), sep = "\t", header = F, stringsAsFactors = F)
  
  sumUMI = c()
  sumi = sum(umi$V4)
  
  for(i in 0:250){
    sumUMI = c(sumUMI, sum(umi$V4[umi$V4>i])/sumi)
  }
  
  pdf("UMIduplication.pdf", height = 3.5, width = 7, useDingbats = F)
  
  par(mfrow = c(1,2))
  
  plot(sumUMI, ylim = c(0,1), pch = 20, col = "grey30", ylab = "% of total reads",
       xlab = "More than xx UMI")
  diffUMI = sumUMI[-length(sumUMI)] - sumUMI[-1]
  plot(diffUMI, ylim = c(0,0.2), pch = 20, col = "grey30", ylab = "Change in % of total reads",
       xlab = "More than xx UMI")
  
  dev.off()
  
  # create Seurat object and add image information
  srat <- CreateSeuratObject(counts = count.data, assay = "Spatial") # create object
  
  image <- Read10X_Image(image.dir = paste0(topdir, "mock/outs/spatial/"),
                         filter.matrix = FALSE) # read in the images
  image <- image[Cells(x = srat)] # filter image by the spots
  DefaultAssay(object = image) <- "Spatial" # set default assay
  srat[["slice1"]] <- image # slice name might be changed
  
  # estimate ambient proportion
  amb_prop = estimateAmbience(count.data)[rownames(srat@assays$Spatial@meta.features)]
  
  srat@assays$Spatial@meta.features = data.frame(row.names = rownames(srat@assays$Spatial@meta.features),
                                                 "ambient_prop" = amb_prop)
  
  # get MT% (genes curated from NCBI chrMT genes)
  # mtgenes = c("COX1", "COX2", "COX3", "ATP6", "ND1", "ND5", "CYTB", "ND2", "ND4",
  #             "ATP8", "MT-CO1", "COI", "LOC9829747")
  # mtgenes = c(mtgenes, paste0("MT", mtgenes), paste0("MT-", mtgenes))
  # mtgenes = mtgenes[mtgenes %in% g[,1]]
  # srat = PercentageFeatureSet(srat, col.name = "percent.mt", assay = "Spatial",
  #                             features = mtgenes)
  
  #saveRDS(srat, file = "srat.RDS")
  
  return(srat)
   
}

########################################################
########################################################
# Section : cell type deconvolution analysis 
# 
########################################################
########################################################
##########################################
# prepare reference data and double check the main cell types and subtypes 
##########################################
Double.check.adult.cardiomyocyte.major.celltypes.subtypes = function(aa)
{
  p1 = DimPlot(aa, reduction = 'umap', group.by = 'seurat_clusters')
  p2 = DimPlot(aa, reduction = "umap", group.by = c("CellType"))
  p1 + p2 
  
  p3 = DimPlot(aa, reduction = 'umap', group.by = 'condition')
  p4 = DimPlot(aa, reduction = 'umap', group.by = 'sample')
  (p1 + p2)/(p3 + p4)
  ggsave(paste0(resDir, '/Umap_Ren2020_week0.week2_newClusters_vs_cellType.original_conditions_samples_seuratNorm.pdf'), 
           width = 12, height = 8)
  
  # double check the cardiomyocyte markers from Elad 
  gg.examples =  c('Nppa', 'Nppb', 'Myh6', 'Tnnc1', 'Tnni3', 
                   'Tnnt2', 'Actn2', 'Gata4', 'Nkx2-5', 'Gja1', 
                   'Myl2', 'Tpm1', 'Ryr2', 'Atp2a2', 'Acta1')
  
  p5 = FeaturePlot(aa, reduction = 'umap', features = gg.examples)
  
  p2 / p5 + ggsave(paste0(resDir, '/Umap_cellType.original_FeatuerPlot_cardiomyocyoteMarkers_seuratNorm.pdf'), 
                   width = 14, height = 20)
  
  p6 = VlnPlot(aa, features = gg.examples, group.by = 'CellType')
  
  p6 + ggsave(paste0(resDir, '/Umap_VlnPlot_cardiomyocyoteMarkers_seuratNorm.pdf'), 
              width = 14, height = 10)
  
  p1 + p2 + ggsave(paste0(resDir, '/Umap_newClusters_vs_cellType.original_seuratNorm.pdf'), 
                   width = 12, height = 8)
  
  #DimPlot(aa, reduction = "umap", group.by = c("SubCluster"))
  
  FeaturePlot(aa, reduction = 'umap', features = c('Csf1r', 'Cd163'))
  FeaturePlot(aa, reduction = 'umap', features = c('S100a3', 'S100a9'))
  
  ##########################################
  # double check the subtypes of CM 
  ##########################################
  aa$celltype = aa$CellType
  aa$subtype = aa$SubCluster
  
  mcells = 'CM'
  ax = subset(aa, cells = colnames(aa)[which(aa$CellType == mcells)])
  table(ax$celltype)
  table(ax$subtype)
  
  ax <- FindVariableFeatures(ax, selection.method = "vst", nfeatures = 5000)
  ax <- ScaleData(ax, features = rownames(ax))
  
  ax <- RunPCA(ax, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(ax, ndims = 30)
  
  # UMAP to visualize subtypes
  ax <- RunUMAP(ax, dims = 1:10, n.neighbors = 20, min.dist = 0.05, n_threads = 6)
  
  DimPlot(ax, reduction = 'umap', group.by = 'subtype') + 
    ggtitle(paste0(mcells, '-', ' cells UMAP (', Normalization, ' nfeature = 5000, ndim=10, neighbors=20, mdist=0.05)'))
  
  ggsave(paste0(resDir, '/Ref_Ren2020_UMAP_', mcells, '_subcelltypes.pdf'), 
         width = 10, height = 8)
  
  ax <- FindNeighbors(ax, dims = 1:10)
  
  ax <- FindClusters(ax, verbose = FALSE, algorithm = 3, resolution = 0.5)
  p1 = DimPlot(ax, reduction = 'umap', group.by = 'seurat_clusters')
  p0 = DimPlot(ax, reduction = 'umap', group.by = 'subtype')
  
  p0 + p1
  
  ax$subtype = paste0('CM', ax$seurat_clusters)
  Idents(ax) = ax$subtype
  
  FeaturePlot(ax, features = c('Sln', 'Myl4', 'Myl7', 'Nppa', 'Myl1', # Atrial CM markers
                               'Fabp3', 'Pln', 'Myl2', 'Myl3', 'Pin', 'Mb' # Ventricular CMs
                               ))
  
 
  VlnPlot(ax, features = c('Sln', 'Myl4', 'Myl7', 'Nppa', 'Myl1', # Atrial CM markers
                           'Fabp3', 'Pln', 'Myl2', 'Myl3', 'Pin', 'Mb' # Ventricular CMs
  ))
  
  ggsave(paste0(resDir, '/VlnPLot_Atrial_Ventricular_markerGenes_', mcells, '_subtypes.pdf'), width = 12, height = 10)
  
  ax.markers <- FindAllMarkers(ax, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.3)
  # saveRDS(ax.markers, file = paste0(RdataDir, 'Forte2020_logNormalize_allgenes_majorCellTypes_markerGenes.rds'))
  
  ax.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top10
  
  DoHeatmap(ax, features = top10$gene)
  
  ggsave(paste0(resDir, '/heatmap_markerGenes_', mcells, '_subtypes.pdf'), width = 12, height = 26)
  
  aa$subtype[match(colnames(ax), colnames(aa))] = ax$subtype
  saveRDS(aa, file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization_umap_subtypes.rds'))
  
  rm(aa)
  rm(ax)
  
}

Double.check.adult.non.cardiomyocyte.major.celltypes.subtypes = function(aa)
{
  # modify the cell annotations
  aa$celltype = as.character(aa$my_annot)
  
  # major cell type FB
  aa$celltype[which(aa$celltype == '0 - Fibro-I')] = 'FB1'
  aa$celltype[which(aa$celltype == '9 - Fibro-II')] = 'FB2'
  aa$celltype[which(aa$celltype == '11 - Fibro-III')] = 'FB3'
  aa$celltype[which(aa$celltype == '15 - Fibro-IV')] = 'FB4'
  aa$celltype[which(aa$celltype == '3 - MyoF')] = 'myoFB'
  
  # major cell type EC
  aa$celltype[which(aa$celltype == '2 - EC-I')] = 'EC1'
  aa$celltype[which(aa$celltype == '18 - EC-II')] = 'EC2'
  aa$celltype[which(aa$celltype == '19 - EC-III')] = 'EC3'
  aa$celltype[which(aa$celltype == '17 - Lymph-EC')] = 'Lymph.EC'
  
  # major cell type smooth muscle cells
  aa$celltype[which(aa$celltype == '16 - SMC')] = 'SMC'
  
  # immune cells
  aa$celltype[which(aa$celltype == '6 - B cell')] = 'B'
  aa$celltype[which(aa$celltype == '4 - GRN')] = 'GN'
  aa$celltype[which(aa$celltype == '12 - NK/T')] = 'NK.T'
  aa$celltype[which(aa$celltype == '13 - Monon/DC')] = 'MCT.DC'
  aa$celltype[which(aa$celltype == '5 - Chil3 Mono')] = 'MCT.Chil3'
  
  aa$celltype[which(aa$celltype == '1 - Trem2 Macs')] = 'Mphage.Trem2'
  aa$celltype[which(aa$celltype == '10 - Arg1 Macs')] = 'Mphage.Argl'
  aa$celltype[which(aa$celltype == '7 - MHCII Macs')] = 'MHCII.Mphage'
  aa$celltype[which(aa$celltype == '14 - Proliferating Macs')] = 'prolife.Mphage'
  
  # merge major cell types
  aa$subtype = aa$celltype
  
  aa$celltype[grep('FB', aa$subtype)]  = 'FB'
  aa$celltype[grep('EC1|EC2|EC3|Lymph.EC', aa$subtype)]  = 'EC'
  aa$celltype[grep('SMC', aa$subtype)] = 'SMC'
  
  aa$celltype[grep('GN|MCT|Mphage|NK', aa$subtype)]  = 'immune'
  aa$celltype[which(aa$subtype == 'B')] = 'immune'
  
  p0 = DimPlot(aa, reduction = 'umap_0.05', group.by = 'celltype') + ggtitle('Shoval UMAP')
  
  p1 = DimPlot(aa, reduction = 'umap', group.by = 'celltype') + ggtitle(paste0(Normalization, ' - Elad umap'))
  
  p0 + p1 
  ggsave(paste0(resDir, '/Ref_Forte2020_majorCelltypes_Shoval.umap_vs_Elad.umap.', Normalization, '.pdf'), 
         width = 18, height = 8)
  
  
  for(mcells in c('all', 'FB', 'EC', 'immune'))
  {
    # mcells = 'immune.others'
    
    if(mcells == 'all'){
      ax = aa
      
      Idents(ax) = ax$celltype
      
    }else{
      ax = subset(aa, cells = colnames(aa)[which(aa$celltype == mcells & aa$subtype != 'B' & aa$subtype != 'GN' &
                                                   aa$subtype != 'NK.T')] )
      table(ax$celltype)
      table(ax$subtype)
      
      ax <- FindVariableFeatures(ax, selection.method = "vst", nfeatures = 3000)
      ax <- ScaleData(ax, features = rownames(ax))
      
      ax <- RunPCA(ax, verbose = FALSE, weight.by.var = TRUE)
      ElbowPlot(ax, ndims = 30)
      
      # UMAP to visualize subtypes
      ax <- RunUMAP(ax, dims = 1:30, n.neighbors = 30, min.dist = 0.05, n_threads = 6)
      
      DimPlot(ax, reduction = 'umap', group.by = 'subtype') + 
      ggtitle(paste0(mcells, '-', ' cells UMAP (', Normalization, ' nfeature = 3000, ndim=30, neighbors=30, mdist=0.05)'))
      
      ggsave(paste0(resDir, '/Ref_Forte2020_UMAP_', mcells, '_excluding B GN NK.T_subcelltypes.pdf'), 
             width = 10, height = 8)
      
      ax <- FindNeighbors(ax, dims = 1:20)
      
      ax <- FindClusters(ax, verbose = FALSE, algorithm = 3, resolution = 0.4)
      p1 = DimPlot(ax, reduction = 'umap', group.by = 'seurat_clusters')
      p0 = DimPlot(ax, reduction = 'umap', group.by = 'subtype')
      
      p0 + p1
      # marker genes and heatmaps
      Idents(ax) = ax$subtype
      
      gg.examples = c('Col1a2', 'Vim', 'Fstl1', 'Ddr2', 'Acta2', 'Postn', 'Tcf21', 'Pdgfra', 'Col3a1', 'Col1a1', 
                      'Gsn', 'Fbln2', 'Sparc', 'Mmp2', 'Msln', 'Rspo1', 'Lum', 'Col8a1')
      
      gg.examples = c('Tek', 'Pecam1', 'Emcn', 'Cdh5', 'Kdr', 'Vwf', 'Fabp4', 'Tie1', 'Flt1', 'Epas1', 'Ednrb', 'Gpihbp1', 'Npr3')
      VlnPlot(ax, features = gg.examples)
      ggsave(paste0(resDir, '/VlnPlot_markerGenes_', mcells, '_subtypes.pdf'), width = 20, height = 16)
      
      
    }
    
    ax.markers <- FindAllMarkers(ax, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
    # saveRDS(ax.markers, file = paste0(RdataDir, 'Forte2020_logNormalize_allgenes_majorCellTypes_markerGenes.rds'))
    
    ax.markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC) -> top10
    
    DoHeatmap(ax, features = top10$gene)
    
    ggsave(paste0(resDir, '/heatmap_markerGenes_', mcells, '_subtypes.pdf'), width = 12, height = 20)
    
  }
  
  ##########################################
  # add distinct immune cells as major cell types  
  ##########################################
  jj = which(aa$subtype == 'B')
  aa$celltype[jj] = aa$subtype[jj] 
  
  jj = which(aa$subtype == 'GN')
  aa$celltype[jj] = aa$subtype[jj] 
  
  jj = which(aa$subtype == 'NK.T')
  aa$celltype[jj] = aa$subtype[jj] 
  
  jj = which(aa$subtype == 'prolife.Mphage')
  aa$celltype[jj] = aa$subtype[jj] 
  
  jj = which(aa$subtype == 'MHCII.Mphage')
  aa$celltype[jj] = aa$subtype[jj] 
  
  jj = which(aa$subtype == 'MCT.DC')
  aa$celltype[jj] = 'immune.others'
  
  jj = which(aa$celltype == 'immune')
  aa$celltype[jj] = 'immune.others'
  
  p0 = DimPlot(aa, reduction = 'umap', group.by = 'celltype')
  p1 = DimPlot(aa, reduction = 'umap', group.by = 'subtype')
  
  p0 + p1
  ggsave(paste0(resDir, '/heatmap_markerGenes_', mcells, '_majoyCelltypes_subtypes.pdf'), width = 12, height = 20)
  
  rm(ax)
  
  xx = DietSeurat(aa, counts = TRUE, data = TRUE, scale.data = FALSE, assays = 'RNA', dimreducs = 'umap')
  rm(aa)
  saveRDS(xx, file = paste0(RdataDir, 'Forte2020_logNormalize_allgenes_majorCellTypes_subtypes.rds'))
   
}


##########################################
# two test functions 
# - test to combine the adult heart datasets without integration
# - test to convert the batch corrected expression to counts (not good idea)
##########################################
Combine.adult.mice.heart.without.integration = function(refs)
{
  # test refs without integration
  Test_refs_withoutIntegration = FALSE
  if(Test_refs_){
    DefaultAssay(refs) <- "RNA"
    refs = FindVariableFeatures(refs, selection.method = "vst", nfeatures = 3000)
    
    # Run the standard workflow for visualization and clustering
    refs <- ScaleData(refs, verbose = FALSE)
    refs <- RunPCA(refs, npcs = 30, verbose = FALSE)
    
    ElbowPlot(refs, ndims = 30)
    
    refs <- FindNeighbors(refs, reduction = "pca", dims = 1:20)
    refs <- FindClusters(refs, resolution = 0.5)
    
    refs <- RunUMAP(refs, reduction = "pca", dims = 1:30, n.neighbors = 50, min.dist = 0.05) 
    
    # Visualization
    p1 <- DimPlot(refs, reduction = "umap", group.by = "dataset")
    p2 <- DimPlot(refs, reduction = "umap", group.by = "annot.ref", label = TRUE,
                  repel = TRUE)
    p1 + p2 + ggsave(paste0(resDir, '/Forte2020_Ren2020_noCorrection_', Normalization, '.pdf'), 
                     width = 24, height = 10)
    
  }
  
}

##########################################
# try to reversely calculated batch-corrected UMI counts using corrected gene expression matrix from Seurat 
##########################################
Convert.batch.corrected.expression.matrix.to.UMIcount = function(refs){
  refs = readRDS(file = paste0(RdataDir, 
                               'SeuratObj_adultMiceHeart_refCombine_Forte2020.nonCM_Ren2020CM_cleanAnnot_logNormalize_v1.rds'))
  
  #jj = which(metadata$dataset == 'Ren2020')
  #aa = readRDS(file = paste0(RdataDir, 'Forte2020_logNormalize_allgenes.rds'))
  #cms = readRDS(file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization_umap.rds'))
  p1 <- DimPlot(refs, reduction = "umap", group.by = "dataset")
  p2 <- DimPlot(refs, reduction = "umap", group.by = "celltype", label = TRUE,
                repel = TRUE)
  p1 + p2
  
  VlnPlot(refs, features = c("nCount_RNA"), group.by = 'dataset')
  
  metadata = refs@meta.data   
  Ec = refs@assays$integrated@data
  
  counts = refs@assays$RNA@counts
  counts = counts[match(rownames(Ec), rownames(counts)), ]
  #E = refs@assays$RNA@data
  
  ss = colSums(counts)
  ccts = expm1(Ec)/10000
  ccts = t(t(ccts)*ss)
  ccts = round(ccts)
  
  metadata$nCount_RNA = ss
  
  aa <- CreateSeuratObject(counts = ccts, project = "adult", min.cells = 50, min.features = 500)
  aa = AddMetaData(aa, metadata, col.name = NULL) 
  
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
  
  aa <- FindNeighbors(aa, dims = 1:10)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
  
  aa = RunUMAP(aa, reduction = "pca", dims = 1:30, n.neighbors = 30, min.dist = 0.05) 
  
  # Visualization
  p1 <- DimPlot(aa, reduction = "umap", group.by = "dataset")
  p2 <- DimPlot(aa, reduction = "umap", group.by = "celltype", label = TRUE,
                repel = TRUE)
  p1 + p2 + ggsave(paste0(resDir, '/refCombined_correctUMIcounts_Forte2020_Ren2020__', Normalization, '.pdf'), 
                   width = 24, height = 10)
  
  refs = aa
  #aa = readRDS(file = paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization.rds'))
  
  saveRDS(refs, file = paste0(RdataDir, 
                              'SeuratObj_adultMiceHeart_refCombine_Forte2020.nonCM_Ren2020CM_cleanAnnot_correctedUMIcounts_v1.rds'))
  
}

##########################################
# Run cell type deconvolution with RCTD for all slices
# inputs are seurat object st and prepared refs
# The deconvolution will be done for each condition of st with a for loop
# For each condition, major cell type is first considered and then subtype in references further considered
# 
##########################################
Run.celltype.deconvolution.RCTD = function(st, refs, Normalization = 'lognormalize', 
                                           save.RCTD.in.seuratObj = FALSE, 
                                           plot.RCTD.summary = TRUE, 
                                           PLOT.scatterpie = TRUE,
                                           run.RCTD.subtype = FALSE)
{
  require(RCTD)
  require(Matrix)
  # slice = 'adult.day4'; Normalization = 'lognormalize';
  #refs = readRDS(file = paste0(RdataDir, 
  #                             'SeuratObj_adultMiceHeart_refCombine_Forte2020.nonCM_Ren2020CM_cleanAnnot_logNormalize_v1.rds'))
  
  # import cardiomyocyte reference 
  # aa = readRDS(file = paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization.rds'))
  # 
  # if(Normalization == 'SCT'){
  #   aa = RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
  #   #saveRDS(aa, file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_SCT_umap.rds'))
  # }else{
  #   aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 20, min.dist = 0.05)
  #   #saveRDS(aa, file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization_umap.rds'))
  # }
  # 
  # aa = aa[ ,!is.na(aa$CellType)]
  # Idents(aa) = aa$CellType
  
  # aa.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
  # gene used for cell type signature
  #genes.used =  intersect(aa.markers$gene[!is.na(match(aa.markers$gene, rownames(aa)))], 
  #                        aa.markers$gene[!is.na(match(aa.markers$gene, rownames(stx)))])
  
  ##########################################
  # prepare references for RTCD
  # original code from https://raw.githack.com/dmcable/RCTD/dev/vignettes/spatial-transcriptomics.html
  # first for major cell types 
  # secondly for subtypes
  ##########################################
  E.corrected <- refs@assays$integrated@data
  E.corrected = expm1(E.corrected) # linear scale of corrected gene expression
  
  ### Create the Reference object for major cell types
  cell_types <- refs@meta.data$celltype
  names(cell_types) <- colnames(refs) # create cell_types named list
  cell_types <- as.factor(cell_types) # convert to factor data type
  reference <- Reference(E.corrected, cell_types, nUMI = NULL, require_int = FALSE)
  
  ## Examine reference object (optional)
  print(dim(reference@counts)) #observe Digital Gene Expression matrix
  #> [1] 384 475
  table(reference@cell_types) #number of occurences for each cell type
  
  ### create reference objects for subtypes
  if(run.RCTD.subtype){
    subtypes <- refs@meta.data$subtype
    names(subtypes) <- colnames(refs) # create cell_types named list
    subtypes <- as.factor(subtypes) # convert to factor data type
    reference_subtype <- Reference(E.corrected, subtypes, nUMI = NULL, require_int = FALSE)
    table(reference_subtype@cell_types)
    
  }
  
  rm(E.corrected)
  
  ##########################################
  # loop over all conditions of st
  ##########################################
  cat('visium conditions :\n')
  print(table(st$condition))
  cc = names(table(st$condition))
  
  for(n in 1:length(cc))
  #for(n in 1:2)
  {
    # n = 2
    slice = cc[n]
    stx = st[, which(st$condition == slice)]
    
    resultsdir <- paste0(resDir, '/RCTD/', slice)
    system(paste0('mkdir -p ', resultsdir))
    
    ##########################################
    # prepare ST data for RTCD
    # original code from https://raw.githack.com/dmcable/RCTD/dev/vignettes/spatial-transcriptomics.html
    ##########################################
    counts <- stx@assays$Spatial@counts
    #counts = counts[!is.na(match(rownames(counts), genes.used)), ]
    coords <- eval(parse(text = paste0('stx@images$',  slice, '@coordinates')))
    coords = coords[, c(4, 5)]
    
    nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
    
    ### Create SpatialRNA object
    puck <- SpatialRNA(coords, counts, nUMI)
    
    ## Examine SpatialRNA object (optional)
    print(dim(puck@counts)) # observe Digital Gene Expression matrix
    hist(log(puck@nUMI,2)) # histogram of log_2 nUMI
    
    print(head(puck@coords)) # start of coordinate data.frame
    barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 
    
    # This list can be restricted if you want to crop the puck e.g. 
    # puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
    # on the plot:
    plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0, round(quantile(puck@nUMI,0.9))), 
                         title ='plot of nUMI') 
    ggsave(paste0(resultsdir, '/RCTD_nUMI_plot_', slice, '.pdf'), width = 12, height = 10)
    
    myRCTD <- create.RCTD(puck, reference, max_cores = 16, gene_cutoff = 0.000125, fc_cutoff = 0.5, 
                          gene_cutoff_reg = 2e-04,  fc_cutoff_reg = 0.75,
                          UMI_min = 100, UMI_max = 2e+07, 
                          CELL_MIN_INSTANCE = 25)
    
    tic()
    myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
    
    saveRDS(myRCTD, file = paste0(RdataDir, 'RCTD_results_refCombined_Forte2020.Ren2020_', slice, '.rds'))
    
    toc()
    
    # myRCTD = readRDS(file = paste0(RdataDir, 'RCTD_results_refCombined_Forte2020.Ren2020_', slice, '.rds'))
    
    results <- myRCTD@results
    
    # normalize the cell type proportions to sum to 1.
    norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
    cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
    
    spatialRNA <- myRCTD@spatialRNA
    
    # make the plots 
    if(plot.RCTD.summary){
      # Plots the confident weights for each cell type as in full_mode (saved as 'results/cell_type_weights_unthreshold.pdf')
      plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
      
      # Plots all weights for each cell type as in full_mode. (saved as 'results/cell_type_weights.pdf')
      plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
      # Plots the weights for each cell type as in doublet_mode. (saved as 'results/cell_type_weights_doublets.pdf')
      
      plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, results$results_df) 
      # Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 'results/cell_type_occur.pdf')
      plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
      
      # makes a map of all cell types, (saved as 
      # 'results/all_cell_types.pdf')
      plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir)
      
    }
    
    ##########################################
    # visualization using scatterpie plot 
    # (original code from 
    # https://github.com/MarcElosua/SPOTlight_deconvolution_analysis/blob/master/analysis/pancreas_PDAC/3b-PDAC_piecharts_immune.Rmd)
    ##########################################
    if(PLOT.scatterpie){
      require(scatterpie)
      require(cowplot)
      library(RColorBrewer)
      
      # set the color vectors for all cell types
      # n <- 60
      # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      # col_vector  <-  unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      # 
      # set color vector
      getPalette <- colorRampPalette(brewer.pal(8, "Dark2"))
      # use a panel of colors from https://gotellilab.github.io/GotelliLabMeetingHacks/NickGotelli/ColorPalettes.html
      tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77",
                         "#661100", "#CC6677", "#882255", "#AA4499")
      cell_types_plt <- sort(unique(cell_type_names))
      #col_vector = getPalette(length(cell_types_plt))
      col_vector = tol10qualitative
      
      col_ct <- col_vector[seq_len(length(cell_types_plt))]
      names(col_ct) = cell_types_plt
      
      #plot(1:length(col_ct), 1:length(col_ct), col = getPalette(length(col_ct)))
      
      ## Preprocess coordinates 
      spatial_coord <-  spatialRNA@coords %>%
        tibble::rownames_to_column("ID")
      
      #colnames(spatial_coord)[2:3] = c('x', 'y')
      spatial_coord = data.frame(spatial_coord, norm_weights[match(spatial_coord$ID, rownames(norm_weights)), ])
      
      ggplot() + geom_scatterpie(aes(x=x, y=y), data=spatial_coord,
                                 cols=cell_type_names, 
                                 color = NA,
                                 alpha = 1, 
                                 pie_scale = 0.4) +
        ggplot2::scale_fill_manual(values = col_ct) + # try to change the color of cell types
        #ggplot2::coord_fixed(ratio = 1) +
        coord_flip() + 
        #ggplot2::scale_y_reverse() +
        ggplot2::scale_x_reverse() + # flip first and reverse x to match seurat Spatial plots
        cowplot::theme_half_open(11, rel_small = 1) +
        ggplot2::theme_void() + 
        
        ggplot2::labs(title = "Spatial scatterpie") +
        ggplot2::theme(
          # plot.background = element_rect(fill = "#FFFFFF"),
          # panel.background = element_blank(),
          # plot.margin = margin(20, 20, 20, 20),
          plot.title = ggplot2::element_text(hjust = 0.5, size = 20)) +
        ggplot2::guides(fill = guide_legend(ncol = 1))
      
      ggsave(paste0(resultsdir, '/RCTD_scatterpie_', slice, '.pdf'), width = 22, height = 16)
      
    }
    
    ##########################################
    # save and assign the cell type into Seurat object
    # default FALSE
    ##########################################
    if(save.RCTD.in.seuratObj){
      cts = results$results_df
      cts$cellType = NA
      cts$cellType[which(cts$spot_class == 'singlet')] = as.character(cts$first_type[which(cts$spot_class == 'singlet')])
      cts$cellType[which(cts$spot_class == 'doublet_certain')] = 
        paste0(cts$first_type[which(cts$spot_class == 'doublet_certain')], 
               '_', cts$second_type[which(cts$spot_class == 'doublet_certain')])
      
      stx$celltype = cts$cellType[match(colnames(stx), rownames(cts))]
      
      xx = SpatialDimPlot(stx, group.by = 'celltype', images = slice, stroke = 0, interactive = FALSE)
      
    }
    
    ##########################################
    # finer annotation by running RCTD with subtypes 
    ##########################################
    if(run.RCTD.subtype){
      cat('run RCTD with subtypes \n')
      cat('to construct ...\n')
      
    }
    
  }
  
}


########################################################
########################################################
# Section : functions for SpatialDE analysis
# test SPARK first (original code from https://xzhoulab.github.io/SPARK/02_SPARK_Example/)
########################################################
########################################################
Find.SpatialDE = function(st, use.sparkX = TRUE)
{
  library('SPARK')
  require(parallel)
  
  cat('visium conditions :\n')
  print(table(st$condition))
  cc = names(table(st$condition))
  
  for(n in 1:length(cc))
  {
    # n = 3
    slice = cc[n]
    stx = st[, which(st$condition == slice)]
    
    resultsdir <- paste0(resDir, '/SpatialDE_SPARK/')
    system(paste0('mkdir -p ', resultsdir))
    
    rawcount = stx@assays$Spatial@counts
    coords <- eval(parse(text = paste0('stx@images$',  slice, '@coordinates')))
    coords = coords[, c(4, 5)]
    
    ##########################################
    # Analyze the data with SPARK-X
    ##########################################
    if(use.sparkX){
      info <- cbind.data.frame(x=as.numeric(coords[,1]),
                               y=as.numeric(coords[, 2]))
      rownames(info) <- colnames(rawcount)
      location        <- as.matrix(info)
      
      tic()
      
      sparkX <- sparkx(rawcount,location,numCores=8,option="mixture")
      
      toc()
      
      head(sparkX$res_mtest)
      
      mtest = sparkX$res_mtest
      mtest = mtest[order(mtest$adjustedPval), ]
      
      write.csv(mtest, file = paste0(resultsdir, 'SpatialDE_sparkX_', slice, '.csv'), row.names = TRUE)
      
    }else{
      info <- cbind.data.frame(x=as.numeric(coords[,1]),
                               y=as.numeric(coords[, 2]),
                               total_counts=apply(rawcount,2,sum))
      
      rownames(info) <- colnames(rawcount)
      
      ## filter genes and cells/spots and 
      spark <- CreateSPARKObject(counts=rawcount, 
                                 location=info[,1:2],
                                 percentage = 0.1, 
                                 min_total_counts = 10)
      
      ## total counts for each cell/spot
      spark@lib_size <- apply(spark@counts, 2, sum)
      
      ## Estimating Parameter Under Null
      tic()
      spark <- spark.vc(spark, 
                        covariates = NULL, 
                        lib_size = spark@lib_size, 
                        num_core = 8,
                        verbose = F)
      
      toc()
      
      ## Calculating pval
      tic()
      spark <- spark.test(spark, 
                          check_positive = T, 
                          verbose = F)
      toc()
      
      ## Output the final results, i.e., combined p-values, adjusted p-values, etc.
      head(spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")])
      
    }
    
  }
 
}



