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
make_SeuratObj_visium = function(topdir = './', saveDir = './results', changeGeneName = TRUE, keyname = 'slice1', 
                                 QC.umi = TRUE)
{
  #topdir = "./" # source dir
  library(Seurat)
  library(DropletUtils)
  library(edgeR)
  
  if(!dir.exists(saveDir)) dir.create(saveDir)
  
  # import read in from kallisto output
  cat('import reads from kallisto output\n')
  exp = Matrix::readMM(paste0(topdir, "genecounts.mtx")) # UMI count matrix
  
  bc = read.csv(paste0(topdir, "genecounts.barcodes.txt"), header = F, stringsAsFactors = F)
  g = read.csv(paste0(topdir, "genecounts.genes.txt"), header = F, stringsAsFactors = F)
  
  if(changeGeneName){
    cat('change gene names \n')
    # change the gene names before making Seurat object
    annot = readRDS(paste0('/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
          'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse_manual_v1.rds'))
    
        
    mm = match(g$V1, annot$geneID)
    ggs = paste0(annot$gene.symbol.toUse[mm], '|',  annot$geneID[mm])
    g$V1[!is.na(mm)] = ggs[!is.na(mm)]
    
  }
  
  dimnames(exp) = list(paste0(bc$V1,"-1"),g$V1) # number added because of seurat format for barcodes
  count.data = Matrix::t(exp)
  
  
  if(QC.umi){
    # plot rankings for number of UMI
    cat('plot UMI rank \n')
    br.out <- barcodeRanks(count.data)
    
    pdf(paste0(saveDir, "UMIrank.pdf"), height = 8, width = 8, useDingbats = FALSE)
    
    plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
    
    o <- order(br.out$rank)
    lines(br.out$rank[o], br.out$fitted[o], col="red")
    abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
    abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
    legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
           legend=c("knee", "inflection"))
    
    dev.off()
    
    # UMI duplication 
    cat('umi duplication check \n')
    umi = read.table(paste0(topdir, "umicount.txt"), sep = "\t", header = F, stringsAsFactors = F)
    colnames(umi) = c('cell.bc', 'umi', 'kallisto.seqIndex', 'counts')
    
    sumUMI = c()
    sumi = sum(umi$counts)
    
    for(i in 0:250){
      sumUMI = c(sumUMI, sum(umi$counts[umi$counts>i])/sumi)
    }
    
    counts.umi = c()
    for(i in 1:100){
      counts.umi = c(counts.umi, length(which(umi$counts == i))/nrow(umi))
    }
    
    
    pdf(paste0(saveDir, "UMIduplication.pdf"), height = 6, width = 12, useDingbats = F)
    
    par(mfrow = c(1,1))
    plot(1:100, counts.umi, ylim = c(0, 1), col = 'gray30', 
         ylab = '% of umi', xlab = 'duplication nb of umi', pch = 20)
    
    
    par(mfrow = c(1,2))
    
    plot(sumUMI, ylim = c(0,1), pch = 20, col = "grey30", ylab = "% of total reads",
         xlab = "More than xx UMI")
    
    diffUMI = sumUMI[-length(sumUMI)] - sumUMI[-1]
    plot(diffUMI, ylim = c(0,0.2), pch = 20, col = "grey30", ylab = "Change in % of total reads",
         xlab = "More than xx UMI")
    
    dev.off()
    
  }
  
  #xx = as.data.frame(count.data)
  
  # create Seurat object and add image information
  cat('creat seurat object with counts and image \n')
  srat <- CreateSeuratObject(counts = count.data, assay = "Spatial") # create object
  
  image <- Read10X_Image(image.dir = paste0(topdir, "mock/outs/spatial/"),
                         filter.matrix = FALSE) # read in the images
  image <- image[Cells(x = srat)] # filter image by the spots
  DefaultAssay(object = image) <- "Spatial" # set default assay
  srat[[keyname]] <- image # slice name might be changed
  
  # estimate ambient proportion
  cat('estimate ambience prop \n')
  amb_prop = DropletUtils::estimateAmbience(count.data)[rownames(srat@assays$Spatial@meta.features)]
  
  srat@assays$Spatial@meta.features = data.frame(row.names = rownames(srat@assays$Spatial@meta.features),
                                                 "ambient_prop" = amb_prop)
  
  # keep only the spots identified for tissue
  cat('keep only spots that have been deteremed to be over tissue \n')
  tissue = read.csv(paste0(topdir, 'mock/outs/spatial/tissue_positions_list.csv'), header =FALSE)
  tissue = tissue[which(tissue$V2>0), ]
  mm = match(tissue$V1, colnames(srat))
  
  srat = subset(srat, cells = tissue$V1[!is.na(mm)])
  cat(ncol(srat), ' spots determined to be tissue  \n')
  
  srat = srat[rowSums(srat@assays$Spatial@counts) > 0, ]
  cat(nrow(srat), ' features with >0 UMI \n')
  
  saveRDS(srat, file = paste0(saveDir, "srat.RDS"))
  
  return(srat)
   
}


########################################################
########################################################
# Section : test SPATA to manual select border zone and remote zones
# 
########################################################
########################################################
manual_selection_spots_image_Spata = function(aa, slice = 'Amex_d7')
{
  #dyn.load("/software/f2021/software/proj/7.2.1-gcccore-10.2.0/lib/libproj.so")
  #dyn.load("/software/f2021/software/gdal/3.2.1-foss-2020b/lib/libgdal.so") 
  #library(rgdal)
  require(SPATA2)
  # slice = cc[n]
  
  aa <- ScaleData(aa, features = rownames(aa), assay = 'Spatial')
  spata_obj = transformSeuratToSpata(aa, sample_name = unique(aa$condition), method = 'spatial', assay_name = 'Spatial', 
                                     assay_slot = 'scale.data', 
                                     image_name = slice, coords_from = 'umap')
  
  #setActiveExpressionMatrix(spata_obj, 'data')
  
  spata_obj <- createSegmentation(object = spata_obj)
  
  plotSegmentation(object = spata_obj, pt_size = 1.9) +
  ggplot2::scale_y_reverse()
  
  #coord_flip() + 
    #ggplot2::scale_y_reverse() +
  #  ggplot2::scale_x_reverse()  # flip first and reverse x to match seurat Spatial plots
  
  #getFeatureVariables(spata_obj, features = "segmentation", return = "data.frame")
  
  aa$segmentation = 'others'
  aa$segmentation[match(getSegmentDf(spata_obj, segment_names = c('border_zone'))$barcodes, colnames(aa))] = 'border_zone'
  aa$segmentation[match(getSegmentDf(spata_obj, segment_names = c('remote_zone1'))$barcodes, colnames(aa))] = 'remote_zone1'
  aa$segmentation[match(getSegmentDf(spata_obj, segment_names = c('remote_zone2'))$barcodes, colnames(aa))] = 'remote_zone2'
  
  return(aa)
  
}

########################################################
########################################################
# Section : run bayesSpace to detect border zone
# original code from 
# https://edward130603.github.io/BayesSpace/articles/ji_SCC.html
########################################################
########################################################
run_bayesSpace = function(st, 
                          outDir = paste0(resDir, '/bayesSpace/'))
{
  ## aa is a seurat objet with one slice / one image
  require(SingleCellExperiment)
  library(BayesSpace)
  library(ggplot2)
  library(patchwork)
  library(scran)
  library(scuttle)
  require(RColorBrewer)
  
  system(paste0('mkdir -p ', outDir))
  
  cat('-- check visium conditions -- \n')
  print(table(st$condition))
  cc = names(table(st$condition))
  
  for(n in 1:length(cc))
  {
    resultsdir <- paste0(RCTD_out, '/', slice)
    ## aa = st
    # DefaultAssay(aa) = 'Spatial'
    slice = names(table(aa$condition))
    scc <- as.SingleCellExperiment(aa, assay = 'Spatial')
    coords <- eval(parse(text = paste0('aa@images$',  slice, '@coordinates')))
    
    ## flip first x y and reverse x to match seurat Spatial plots
    scc$row = coords$row
    scc$col = coords$col
    scc$imagerow = coords$imagerow
    scc$imagecol = coords$imagecol
    
    clusters <- quickCluster(scc)
    scc <- computeSumFactors(scc, clusters=clusters)
    summary(sizeFactors(scc))
    
    scc <- logNormCounts(scc)
    
    set.seed(101)
    dec <- scran::modelGeneVar(scc)
    top <- scran::getTopHVGs(dec, n = 2000)
    
    set.seed(102)
    scc <- scater::runPCA(scc, subset_row=top)
    
    ## Add BayesSpace metadata
    scc <- spatialPreprocess(scc, platform="Visium", skip.PCA=TRUE)
    
    scc <- qTune(scc, qs=seq(4, 15))
    qPlot(scc)
    
    ggsave(filename =  paste0(resDir, "/BayesSpace_SpatialCluster_cluster.nb.selection_", species, '_', slice, ".pdf"), 
           width = 12, height = 8)
    
    # sptial clustering 
    q <- 8  # Number of clusters
    d <- 15  # Number of PCs
    
    library(mclust) ## Here we run mclust externally so the random seeding is consistent with ## original analyses
    Y <- reducedDim(scc, "PCA")[, seq_len(d)]
    set.seed(101)
    init <- Mclust(Y, q, "EEE", verbose=FALSE)$classification
    
    ## Run BayesSpace clustering
    set.seed(100)
    scc <- spatialCluster(scc, q=q, d=d, platform='Visium', init=init,
                          nrep=10000, gamma=3)
    
    
    # scc = readRDS(file = paste0(RdataDir, "/BayesSpace_SpatialSlustered_", species, '_', slice,  "_with.", q, "clusters.rds"))
    p1 = SpatialDimPlot(aa, group.by = 'spatial_domain_manual', images = slice)
    
    palette <- RColorBrewer::brewer.pal(q, "Paired")
    #palette =  c(brewer.pal(name="Dark2", n = (q-8)), brewer.pal(name="Paired", n = 8))
    
    p2 = clusterPlot(scc, palette=palette, size=0.1) +
      labs(title= paste0("Spot-level clustering : ",  slice)) 
    # + 
    #coord_flip() + 
    #ggplot2::scale_y_reverse() +
    #ggplot2::scale_x_reverse()  # flip first and reverse x to match seurat Spatial plots
    
    p1 + p2  
    
    ggsave(filename =  paste0(resDir, "/BayesSpace_SpatialSlustered_", species, '_', slice, "manual_vs_bayeSpace_with.", 
                              q, ".clusters.pdf"), width = 16, height = 8)
    
    saveRDS(scc, file = paste0(RdataDir, "/BayesSpace_SpatialSlustered_", species, '_', slice,  "_with.", q, "clusters.rds"))
    
    aa$spatial_domain_bayeSpace = NA
    aa$spatial_domain_bayeSpace[match(colnames(scc), colnames(aa))] = scc$spatial.cluster
    
    SpatialDimPlot(aa, group.by = 'spatial_domain_bayeSpace', images = slice)
    
    saveRDS(aa, file = paste0(paste0(RdataDir, "/SeuratObj_spatialDomain_BayesSpace_", 
                                     species, '_', slice,  "_with.", q, "clusters.rds")))
    
    Find.top.markers.for.spatial.clusters = FALSE
    if(Find.top.markers.for.spatial.clusters){
      # top markers of spatial clusters 
      library(dplyr)
      
      ## Convert SCE to Seurat object and use BayesSpace cluster as identifier
      sobj <- Seurat::CreateSeuratObject(counts=logcounts(scc),
                                         assay='Spatial',
                                         meta.data=as.data.frame(colData(scc)))
      sobj <- Seurat::SetIdent(sobj, value = "spatial.cluster")
      
      
      ## Scale data
      sobj = Seurat::ScaleData(sobj)
      #sobj@assays$Spatial@scale.data <-
      #  sobj@assays$Spatial@data %>% as.matrix %>% t %>% scale %>% t
      
      ## Select top n markers from each cluster (by log fold change)
      top_markers <- Seurat::FindAllMarkers(sobj, assay='Spatial', slot='data',
                                            group.by='spatial.cluster',
                                            only.pos=TRUE) 
      
      top_markers %>% group_by(cluster) %>%
        top_n(10, avg_log2FC) -> top10
      
      ## Plot expression of markers
      Seurat::DoHeatmap(sobj, features = top10$gene, 
                        group.by = "spatial.cluster", group.colors=palette, 
                        angle=0, size=4, label = FALSE, raster=FALSE) + 
        guides(col = FALSE)
      
      SpatialFeaturePlot(aa, features = "CTSS-AMEX60DD007394" )
      
    }
    
    Run.bayesSpace.enhanced.clustering = FALSE # this takes > 1 hour at least long time
    if(Run.bayesSpace.enhanced.clustering){
      ## Run BayesSpace enhanced clustering
      set.seed(100)
      scc.enhanced <- spatialEnhance(scc, q=q, d=d, platform="Visium",
                                     nrep=20000, gamma=3, verbose=TRUE,
                                     jitter_scale=5.5, jitter_prior=0.3,
                                     save.chain=FALSE)
      
      
      # We compared the two clusterings using clusterPlot()
      enhanced.plot <- clusterPlot(scc.enhanced, palette=palette, size=0.05) +
        labs(title= paste0("Enhanced clustering :", slice)) +
        coord_flip() + 
        ggplot2::scale_x_reverse()  # flip first and reverse x to match seurat Spatial plots
      
      spot.plot + enhanced.plot
      
    }
    
  }
  
  
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
      
      gg.examples = c('Tek', 'Pecam1', 'Emcn', 'Cdh5', 'Kdr', 'Vwf', 'Fabp4', 'Tie1', 
                      'Flt1', 'Epas1', 'Ednrb', 'Gpihbp1', 'Npr3')
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
  
  ggsave(filename = paste0(resDir, '/AdultMice_scRNAref_overView.pdf'),  width = 20, height = 8)
  
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
Run.celltype.deconvolution.RCTD = function(st, # spatial transcriptome seurat object 
                                           refs, # scRNA-seq reference with annotation names celltype_toUse
                                           RCTD_out = '../results/RCTD_out', 
                                           require_int_SpatialRNA = FALSE,
                                           max_cores = 16,
                                           Normalization = 'lognormalize', 
                                           save.RCTD.in.seuratObj = FALSE, 
                                           plot.RCTD.summary = TRUE, 
                                           PLOT.scatterpie = TRUE,
                                           run.RCTD.subtype = FALSE)
{
  library(spacexr)
  require(Matrix)
  
  ##########################################
  # prepare references for RTCD
  # original code from https://raw.githack.com/dmcable/RCTD/dev/vignettes/spatial-transcriptomics.html
  # first for major cell types 
  # secondly for subtypes
  ##########################################
  cat('-- prepare reference --\n')
  
  #E.corrected <- refs@assays$integrated@data
  #E.corrected = expm1(E.corrected) # linear scale of corrected gene expression
  E_refs = GetAssayData(object = refs, slot = "data")
  E_refs = expm1(E_refs) # linear scale of corrected gene expression
  
  ### Create the Reference object for major cell types
  cell_types <- refs@meta.data$celltype_toUse
  names(cell_types) <- colnames(refs) # create cell_types named list
  cell_types <- as.factor(cell_types) # convert to factor data type
  reference <- Reference(E_refs, cell_types, nUMI = NULL, require_int = FALSE)
  
  rm(E_refs, cell_types)
  
  ## Examine reference object (optional)
  print(dim(reference@counts)) #observe Digital Gene Expression matrix
  #> [1] 384 475
  table(reference@cell_types) #number of occurences for each cell type
  
  ### create reference objects for subtypes
  # if(run.RCTD.subtype){
  #   subtypes <- refs@meta.data$subtype
  #   names(subtypes) <- colnames(refs) # create cell_types named list
  #   subtypes <- as.factor(subtypes) # convert to factor data type
  #   reference_subtype <- Reference(E.corrected, subtypes, nUMI = NULL, require_int = FALSE)
  #   table(reference_subtype@cell_types)
  # }
  # 
  # #rm(E.corrected)
  
  ##########################################
  # loop over all conditions of st for now
  ##########################################
  cat('-- check visium conditions -- \n')
  print(table(st$condition))
  cc = names(table(st$condition))
  
  for(n in 1:length(cc))
  #for(n in c(1, 2, 4))
  {
    # n = 4
    cat('slice -- ', cc[n], '\n')
    slice = cc[n]
    stx = st[, which(st$condition == slice)]
    
    resultsdir <- paste0(RCTD_out, '/', slice)
    system(paste0('mkdir -p ', resultsdir))
    
    ##########################################
    # prepare ST data for RTCD
    # original code from https://raw.githack.com/dmcable/RCTD/dev/vignettes/spatial-transcriptomics.html
    ##########################################
    counts = GetAssayData(object = st, slot = "counts")
    #counts <- stx@assays$Spatial@counts
    
    coords <- eval(parse(text = paste0('stx@images$',  slice, '@coordinates')))
    coords = coords[, c(4, 5)]
    
    nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
    
    ### Create SpatialRNA object, require_int = FALSE, because of kallisto output
    puck <- SpatialRNA(coords, as.matrix(counts), nUMI, require_int = require_int_SpatialRNA)
    
    ## Examine SpatialRNA object (optional)
    print(dim(puck@counts)) # observe Digital Gene Expression matrix
    hist(log(puck@nUMI, 2)) # histogram of log_2 nUMI
    
    print(head(puck@coords)) # start of coordinate data.frame
    barcodes <- colnames(puck@counts) # pixels to be used (a list of barcode names). 
    
    # This list can be restricted if you want to crop the puck e.g. 
    # puck <- restrict_puck(puck, barcodes) provides a basic plot of the nUMI of each pixel
    # on the plot:
    plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0, round(quantile(puck@nUMI,0.9))), 
                         title ='plot of nUMI') 
    ggsave(paste0(resultsdir, '/RCTD_nUMI_plot_', slice, '.pdf'), width = 12, height = 10)
    
    # make RCTD object
    myRCTD <- create.RCTD(puck, reference, max_cores = max_cores, 
                          gene_cutoff = 0.000125, fc_cutoff = 0.5, 
                          gene_cutoff_reg = 2e-04,  
                          fc_cutoff_reg = 0.75,
                          UMI_min = 20, UMI_max = 2e+07, 
                          CELL_MIN_INSTANCE = 30)
    
    doubleCheck.markerGenes = FALSE
    if(doubleCheck.markerGenes){
      require(pheatmap)
      yy = myRCTD@cell_type_info$info[[1]]
      ggs = myRCTD@internal_vars$gene_list_reg
      yy = yy[!is.na(match(rownames(yy), ggs)), ]
      #yy[which(yy == 0)] = 10^-7
      
      #yy[which(is.na(yy))] 
      yy = log10(yy + 10^-7)
      
      pheatmap(yy, cluster_cols =  FALSE, scale = 'row', show_rownames = FALSE)
      ggsave(filename = paste0(resultsdir, '/MarkerGenes_mainCelltypes_used_inRCTD.pdf'), width = 12, height = 10)
        
    }
    
    # run RCTD main function
    tic()
    myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
    saveRDS(myRCTD, file = paste0(resultsdir, '/RCTD_out_doubletMode_', slice, '.rds'))
    toc()
    
    ##########################################
    # check the result
    ##########################################
    myRCTD = readRDS(file = paste0(resultsdir, '/RCTD_out_doubletMode_', slice, '.rds'))
    results <- myRCTD@results
    
    # normalize the cell type proportions to sum to 1.
    norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
    cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
    
    spatialRNA <- myRCTD@spatialRNA
    
    # make the plots 
    if(plot.RCTD.summary){
      
      # Plots the confident weights for each cell type as in full_mode 
      # (saved as 'results/cell_type_weights_unthreshold.pdf')
      plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
      
      # Plots all weights for each cell type as in full_mode. (saved as 'resultrs/cell_type_weights.pdf')
      plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
      
      # Plots the weights for each cell type as in doublet_mode. 
      # (saved as 'results/cell_type_weights_doublets.pdf')
      plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, results$results_df) 
      
      # Plots the number of confident pixels of each cell type in 'full_mode'. 
      # (saved as 'results/cell_type_occur.pdf')
      plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
      
      # makes a map of all cell types, (saved as 
      # 'results/all_cell_types.pdf')
      plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir)
      
      # doublets
      #obtain a dataframe of only doublets
      doublets <- results$results_df[results$results_df$spot_class == "doublet_certain",] 
      # Plots all doublets in space (saved as 
      # 'results/all_doublets.pdf')
      plot_doublets(spatialRNA, doublets, resultsdir, cell_type_names) 
      
      # Plots all doublets in space for each cell type (saved as 
      # 'results/all_doublets_type.pdf')
      plot_doublets_type(spatialRNA, doublets, resultsdir, cell_type_names) 
      # a table of frequency of doublet pairs 
      doub_occur <- table(doublets$second_type, doublets$first_type) 
      # Plots a stacked bar plot of doublet ocurrences (saved as 
      # 'results/doublet_stacked_bar.pdf')
      
      plot_doub_occur_stack(doub_occur, resultsdir, cell_type_names) 
      
    }
    
    ##########################################
    # visualization using scatterpie plot 
    # (original code from 
    # https://github.com/MarcElosua/SPOTlight_deconvolution_analysis/blob/master/
    # analysis/pancreas_PDAC/3b-PDAC_piecharts_immune.Rmd)
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
      getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
      
      # use a panel of colors from https://gotellilab.github.io/GotelliLabMeetingHacks/NickGotelli/ColorPalettes.html
      #tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77",
      #                   "#661100", "#CC6677", "#882255", "#AA4499")
      
      ## Preprocess coordinates 
      spatial_coord <-  spatialRNA@coords %>%
        tibble::rownames_to_column("ID")
      
      weights = norm_weights[match(spatial_coord$ID, rownames(norm_weights)), ]
      celltype_keep = colnames(weights)
      
      Process_celltype_weight_singlet = FALSE
      if(Process_celltype_weight_singlet){
        ## process the cell type weights
        dfs = data.frame(results$results_df, results$weights_doublet)
        colnames(dfs)[10:11] = c('first_type_weights', 'second_type_weights')
        
        #weights = weights[, match(cell_types_plt, colnames(weights))]
        
        celltype_keep = c()
        for(j in 1:nrow(weights))
        {
          # j = 1
          cat(j, '\n')
          if(dfs$spot_class[j] == 'singlet'){
            weights[j, which(colnames(weights) == dfs$first_type[j])] = 1.0
            weights[j, which(colnames(weights) != dfs$first_type[j])] = 0.0
            celltype_keep = c(celltype_keep, as.character(dfs$first_type[j]))
          }
          
          if(dfs$spot_class[j] == 'doublet_certain'){
            weights[j, which(colnames(weights) == dfs$first_type[j])] = dfs$first_type_weights[j]
            weights[j, which(colnames(weights) == dfs$second_type[j])] = dfs$second_type_weights[j]
            weights[j, which(colnames(weights) != dfs$first_type[j] & colnames(weights) != dfs$second_type[j])]  = 0.0
            celltype_keep = c(celltype_keep, c(as.character(dfs$first_type[j]), as.character(dfs$second_type[j])))
          }
        }
        
      }
      
      #ss = colSums(weights)
      cell_types_plt = unique(celltype_keep)
      #ss = ss[match(cell_types_plt, names(ss))]
      #cell_types_counts = table(results$results_df$first_type)
      #cell_types_counts = cell_types_counts[which(cell_types_counts>0)]
      #cell_types_counts = cell_types_counts[order(-cell_types_counts)]
      #cell_types_plt <- cell_types_plt[order(-ss)]
      col_vector = getPalette(length(cell_types_plt))
      #names(col_vector) = cell_types_plt
      #col_ct = readRDS(file = paste0(RdataDir, 'main_celltype_colors_4RCTD.rds'))
      col_ct <- col_vector[seq_len(length(cell_types_plt))]
      #names(col_ct) = cell_types_plt
      #saveRDS(col_ct, file = paste0(RdataDir, 'main_celltype_colors_4RCTD.rds'))
      
      plot(1:length(col_ct), 1:length(col_ct), col = col_ct)
      text(1:length(col_ct), 1:length(col_ct), names(col_ct), offset = 2, adj = 0.5,  col = col_ct)
      
      spatial_coord = data.frame(spatial_coord, weights)
      spatial_coord$x = as.numeric(spatial_coord$x)
      spatial_coord$y = as.numeric(spatial_coord$y)
      
      
      ggplot() + geom_scatterpie(aes(x=x, y=y), data=spatial_coord,
                                 cols=colnames(spatial_coord)[-c(1:3)], 
                                 #color = NA,
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
      
      ggsave(paste0(resultsdir, '/RCTD_scatterpie_', slice, '_v3.pdf'), width = 22, height = 16)
      
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
# Section : Analyze the cell type proximity
# original code from https://github.com/madhavmantri/chicken_heart/blob/master/scripts/anchor_integeration.R
# and the paper (https://www.nature.com/articles/s41467-021-21892-z#data-availability)
# 
########################################################
########################################################
run_neighborhood_analysis = function(aa)
{
  # import the manual spatial domain or by bayesSpace
  slice = "adult.day4"
  q = 8
  
  aa =  readRDS(file = paste0(paste0(RdataDir, "/SeuratObj_spatialDomain_BayesSpace_", 
                            species, '_', slice,  "_with.", q, "clusters.rds")))
  
  SpatialDimPlot(aa, group.by = 'spatial_domain_bayeSpace', images = slice, label = TRUE)
  
  # import the cell type deconvolution result
  myRCTD = readRDS(file = paste0(RdataDir, 'RCTD_results_refCombined_Forte2020.Ren2020_doubletMode_rmSMC_', slice, '.rds'))
  results <- myRCTD@results
  
  # normalize the cell type proportions to sum to 1.
  norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
  
  # use a threshold to binarize weights
  weights = norm_weights >= 0.05 
  
  index_ctl = match(colnames(aa)[which(aa$spatial_domain_bayeSpace == '4')], rownames(weights))
  index_border = match(colnames(aa)[which(aa$spatial_domain_bayeSpace == '7')], rownames(weights))
  
  ctl = weights[index_ctl, -5]
  border = weights[index_border, -5]
  
  ss1 = apply(ctl, 2, sum)
  ss2 = apply(border, 2, sum)
  
  #bgs = bgs/sum(bgs)
  #border = border[, which(bgs>0)]
  #bgs = bgs[which(bgs>0)]
  cell_type_names <- colnames(ctl)
  
  coc = matrix(0, nrow = length(cell_type_names), ncol = length(cell_type_names))
  colnames(coc) = cell_type_names
  rownames(coc) = cell_type_names
  
  for(m in 1:ncol(coc))
  {
    # m = 2
    jj = which(border[,m] == TRUE)
    jj0 = which(ctl[,m] == TRUE)
    
    if(length(jj) == 1){
      xx = border[jj, ]
    }
    if(length(jj) > 1){
      xx = apply(border[jj, ], 2, sum)
    }
    
    if(length(jj0) == 1){
      xx0 = ctl[jj0, ]
    }
    if(length(jj0) > 1){
      xx0 = apply(ctl[jj0, ], 2, sum)
    }
    
    #xx = xx/length(jj)
    xx = xx/nrow(border)
    xx0 = (xx0+0.1)/nrow(ctl)
    #xx = xx/bgs[m]
    coc[m,] = xx/xx0
    
  }
  
  coc = log10(coc + 10^-2)
  library(corrplot)
  require(RColorBrewer)
  corrplot(coc, method = 'color', is.corr = FALSE, hclust.method = c("ward.D2"),
           col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(10))
  
  ggsave(filename =  paste0(resDir, "/celltype_cooccurence_score_border_vs_", species, '_', slice, ".pdf"), 
         width = 12, height = 8)
  
  
}
##########################################
# The original code from somewhere for the cell type proximity network analysis 
# paper https://www.nature.com/articles/s41467-021-21892-z
# https://github.com/madhavmantri/chicken_heart/blob/master/scripts/anchor_integeration.R
##########################################
analyze.celltype.proximity.network = function()
{
  # setRepositories(ind = 1:2)
  # options(unzip = "internal")
  # devtools::install_github("satijalab/seurat", ref = "spatial")
  # library(Seurat, lib.loc = "/home/mm2937/miniconda3/r-seurat-3.1.1-r35h0357c0b_0/lib/R/library/")
  library(Seurat)
  packageVersion("Seurat")
  
  library(Matrix); library(stringr)
  library(readr); library(here)
  library(fitdistrplus); library(dplyr)
  library(SeuratData); library(ggplot2)
  library(cowplot); library(reticulate)
  library(pals); library(monocle)
  setwd("/workdir/mm2937/chicken/") 
  
  # load(here("robjs", "chicken_normalised_scanorama3.Robj"))
  # load("robjs/chicken_visium.4.Robj")
  
  DefaultAssay(chicken.integrated) <- "RNA"
  DefaultAssay(chicken_visium) <- "Spatial"
  
  # Find gene anchors between scRNAseq and spatila RNAseq datasets
  anchors <- FindTransferAnchors(reference = chicken.integrated, query = chicken_visium, reduction = "cca")
  
  # Transfer labels from scRNAseq to spatial using the anchor based approach
  chicken.integrated$cellname <- colnames(chicken.integrated)
  table(chicken.integrated$celltypes.0.5)
  predictions.assay <- TransferData(anchorset = anchors, refdata = chicken.integrated$celltypes.0.5, prediction.assay = TRUE, 
                                    weight.reduction = chicken_visium[["pca"]])
  # save(anchors, predictions.assay, file = "robjs/anchors_and_prediction_assay.RData")
  # load("robjs/anchors_and_prediction_assay.RData")
  
  # Adding cell type predictions to original seurat object
  chicken_visium[["predictions"]] <- predictions.assay
  dim(GetAssayData(chicken_visium, assay = "predictions"))
  
  # Adding cell type predictions in meta data as well
  chicken_visium <- AddMetaData(chicken_visium, metadata = as.data.frame(t(GetAssayData(chicken_visium, assay = "predictions"))))
  head(chicken_visium@meta.data)
  
  # Define cell type with maximum prediction score as spot type 
  prediction.scores <- as.data.frame(t(GetAssayData(chicken_visium, assay = "predictions")))
  prediction.scores$max <- NULL
  sum(is.na(prediction.scores))
  prediction.scores$celltype_prediction <- NA
  dim(prediction.scores)
  for(i in 1:nrow(prediction.scores)){
    prediction.scores$celltype_prediction[i] <- 
      colnames(prediction.scores)[prediction.scores[i,1:15] == max(prediction.scores[i,1:15])]
  }
  
  table(prediction.scores$celltype_prediction)
  chicken_visium$celltype_prediction <- prediction.scores$celltype_prediction
  
  # save(chicken_visium, file="robjs/chicken_visium.4.prediction.1.Robj")
  load("robjs/chicken_visium.4.prediction.1.Robj")
  
  #####################  This sections calclulates the celltype pair neighborhood maps  ############################
  # Example shown for D10 (Run this 4 times for individual stages
  prediction.scores <- as.data.frame(t(GetAssayData(chicken_visium, assay = "predictions")))
  prediction.scores$max <- NULL
  dim(prediction.scores)
  prediction.scores.1 <- prediction.scores[colnames(chicken_visium)[chicken_visium$orig.ident == "D10"],]
  dim(prediction.scores.1)
  interaction_matrix = matrix(0, ncol = length(unique(chicken_visium$celltype_prediction)), 
                              nrow = length(unique(chicken_visium$celltype_prediction)))
  rownames(interaction_matrix) <- unique(chicken_visium$celltype_prediction)
  colnames(interaction_matrix) <- unique(chicken_visium$celltype_prediction)
  for(i in 1:nrow(prediction.scores.1)){
    temp <- colnames(sort(prediction.scores.1[i,prediction.scores.1[i,] > 0], decreasing = T))
    if(length(temp) == 2){
      interaction_matrix[temp[1], temp[2]] <- interaction_matrix[temp[1], temp[2]] + 1
    } else if(length(temp) == 3){
      interaction_matrix[temp[1], temp[2]] <- interaction_matrix[temp[1], temp[2]] + 1
      interaction_matrix[temp[2], temp[3]] <- interaction_matrix[temp[2], temp[3]] + 1
      interaction_matrix[temp[1], temp[3]] <- interaction_matrix[temp[1], temp[3]] + 1
    } else if(length(temp) >= 4){
      interaction_matrix[temp[1], temp[2]] <- interaction_matrix[temp[1], temp[2]] + 1
      interaction_matrix[temp[2], temp[3]] <- interaction_matrix[temp[2], temp[3]] + 1
      interaction_matrix[temp[3], temp[4]] <- interaction_matrix[temp[3], temp[4]] + 1
      interaction_matrix[temp[1], temp[3]] <- interaction_matrix[temp[1], temp[3]] + 1
      interaction_matrix[temp[1], temp[4]] <- interaction_matrix[temp[1], temp[4]] + 1
      interaction_matrix[temp[2], temp[4]] <- interaction_matrix[temp[2], temp[4]] + 1
    }
  }
  
  interaction_matrix <- interaction_matrix + t(interaction_matrix)
  colnames(interaction_matrix)
  temp <- colnames(interaction_matrix)[!colnames(interaction_matrix) %in% c("Erythrocytes", 
                                                                            "Macrophages", 
                                                                            "Mitochondria enriched cardiomyocytes")]
  interaction_matrix <- interaction_matrix[temp, temp]
  
  library(pals)
  color_pelette <- rev(as.vector(kelly()[3:(2+length(levels(chicken.integrated$celltypes.0.5)))]))
  names(color_pelette) <- levels(chicken.integrated$celltypes.0.5)
  
  # write.csv(interaction_matrix, file = "interactions-D10.csv")
  # interaction_matrix <- read.csv("interactions-D10.csv", row.names = 1)
  interaction_matrix[lower.tri(interaction_matrix)] <- 0
  
  library(circlize)
  color_used <- color_pelette[colnames(interaction_matrix)]
  row_col <- color_used
  row_col[names(row_col) != "TMSB4X high cells"] <- "#cecece"
  
  col <- matrix(rep(color_used, each = ncol(interaction_matrix), T), nrow = nrow(interaction_matrix), 
                ncol = ncol(interaction_matrix))
  rownames(col) <- rownames(interaction_matrix)
  colnames(col) <- colnames(interaction_matrix)
  
  chordDiagram(interaction_matrix, grid.col = color_used, col = col, annotationTrack = "grid")
  
  col[rownames(col)!= "TMSB4X high cells", colnames(col) != "TMSB4X high cells"] <- "#cecece"
  chordDiagram(interaction_matrix, grid.col = color_used, col = col, annotationTrack = "grid")
  
  # col[rownames(col)!= "TMSB4X high cells", colnames(col) != "TMSB4X high cells"] <- "#cecece"
  
  pdf(file="Fig1e.4.pdf",
      width=1.5, height=1.5, paper="special", bg="transparent",
      fonts="Helvetica", colormodel = "rgb")
  chordDiagram(interaction_matrix, grid.col = color_used, col = col, annotationTrack = "grid")
  dev.off()
  
  #####################  This section saves cell ids for visium samples  ############################
  
  table(chicken_visium$orig.ident)
  colnames(chicken_visium)
  Images(chicken_visium)
  sample_cell_id_map <- data.frame(sample = chicken_visium$orig.ident, 
                                   cell_id = str_split_fixed(colnames(chicken_visium), "_", 2)[,2])
  head(sample_cell_id_map)
  # save(sample_cell_id_map, file="robjs/sample_cell_id_map.Robj")
  load("robjs/sample_cell_id_map.Robj")
  
  
  #####################  This section calculates the cell spot similarity map ############################
  # Transfer cellnames from scRNAseq to spatial using the anchor based approach to get a cell spot similairy map
  chicken.integrated$cellname <- colnames(chicken.integrated)
  predictions.assay <- TransferData(anchorset = anchors, refdata = chicken.integrated$cellname, prediction.assay = TRUE, 
                                    weight.reduction = chicken_visium[["pca"]])
  
  # Adding cellname predictions to original seurat object
  chicken_visium[["predictions_cells"]] <- predictions.assay
  dim(GetAssayData(chicken_visium, assay = "predictions_cells"))
  
  # save(chicken_visium, file="robjs/chicken_visium.4.prediction.1.Robj")
  load("robjs/chicken_visium.4.prediction.1.Robj")
  
  
  #####################  This section uses the cell spot similarity map and
  ##################### defines spot type in 2 differet ways (optional: not used in manuscript) 
  ############################
  
  # Spot type defined by cell type with maxium 
  prediction.scores <- as.data.frame(t(GetAssayData(chicken_visium, assay = "predictions_cells")))
  # prediction.scores$max <- NULL
  sum(is.na(prediction.scores))
  prediction.scores$cellprediction_max <- NA
  dim(prediction.scores)
  for(i in 1:nrow(prediction.scores)){
    prediction.scores$cellprediction_max[i] <- colnames(prediction.scores)[prediction.scores[i,1:22315] == 
                                                                             prediction.scores$max[i]]
  }
  prediction.scores$cellprediction_max
  sum(is.na(prediction.scores$cellprediction_max))
  temp <- str_replace_all(prediction.scores$cellprediction_max, pattern = "V-", replacement = "V_") 
  temp <- str_replace_all(temp, pattern = "D4-", replacement = "D4_")
  sum(!temp %in% rownames(chicken.integrated@meta.data))
  prediction.scores$celltype_prediction_max <- chicken.integrated$celltypes.0.5[temp]
  dim(chicken.integrated)
  
  prediction.scores$celltype_prediction_mode <- NA
  dim(prediction.scores)
  for(i in 1:nrow(prediction.scores)){
    temp <- table(chicken.integrated$celltypes.0.5[prediction.scores[i,1:22315] > 0.0])
    prediction.scores$celltype_prediction_mode[i] <- names(temp)[temp == max(temp)]
  }
  table(prediction.scores$celltype_prediction_mode)
  
  sum(prediction.scores$celltype_prediction_max == chicken_visium$celltype_prediction) 
  sum(prediction.scores$celltype_prediction_mode == chicken_visium$celltype_prediction)
  
  chicken_visium$celltype_prediction_max <- prediction.scores$celltype_prediction_max
  chicken_visium$celltype_prediction_mode <- prediction.scores$celltype_prediction_mode
  
  # save(chicken_visium, file="robjs/chicken_visium.4.prediction.1.Robj")
  load("robjs/chicken_visium.4.prediction.1.Robj")
  
  
}

########################################################
########################################################
# Section : functions for SpatialDE analysis
# test SPARK first (original code from https://xzhoulab.github.io/SPARK/02_SPARK_Example/)
########################################################
########################################################
Find.SpatialDE = function(aa, use.method = 'sparkX')
{
  
  # use.method = c('spark', 'sparkX', 'moransi')
  cat('visium conditions :\n')
  print(table(aa$condition))
  slice = names(table(aa$condition))
  
  stx = aa
  
  ##########################################
  # Analyze the data with SPARK-X or monansi in Seurat
  ##########################################
  if(use.method == 'monansi'){
    DefaultAssay(stx) <- "SCT"
    stx <- FindSpatiallyVariableFeatures(stx, 
                                         assay = "SCT", 
                                         slot = "scale.data", 
                                         features = VariableFeatures(stx)[1:1000],
                                         selection.method = "moransi", 
                                         x.cuts = 100, y.cuts = 100)
    
  }else{
    require(SPARK)
    require(parallel)
    
    if(use.method == 'sparkX'){
      
      rawcount = stx@assays$Spatial@counts
      coords <- eval(parse(text = paste0('stx@images$',  slice, '@coordinates')))
      coords = coords[, c(4, 5)]
      
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
      
      #resultsdir <- paste0(resDir, '/SpatialDE_SPARKX/')
      #system(paste0('mkdir -p ', resultsdir))
      #write.csv(mtest, file = paste0(resultsdir, 'SpatialDE_sparkX_', slice, '.csv'), row.names = TRUE)
      
    }
    
    if(use.method == 'spark'){
      
      resultsdir <- paste0(resDir, '/SpatialDE_SPARK/')
      system(paste0('mkdir -p ', resultsdir))
      
      rawcount = stx@assays$Spatial@counts
      coords <- eval(parse(text = paste0('stx@images$',  slice, '@coordinates')))
      coords = coords[, c(4, 5)]
      
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
  
  return(mtest)
 
}

########################################################
########################################################
# Section : functions of cell-cell interaction, ligand-receptor analysis
# 
########################################################
########################################################
run_LIANA = function() # original code from https://saezlab.github.io/liana/articles/liana_tutorial.html
{
  require(tidyverse)
  require(magrittr)
  require(liana)
  require(scran)
  require(scater)
  require(DelayedMatrixStats)
  require(genefilter)
  require(matrixStats)
  require(MatrixGenerics)
  require(sparseMatrixStats)
  
  # Resource currently included in OmniPathR (and hence `liana`) include:
  show_resources()
  
  # Resource currently included in OmniPathR (and hence `liana`) include:
  show_methods()
  
  # liana_path <- system.file(package = "liana")
  # testdata <- readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
  # 
  # testdata %>% glimpse()
  
  sels =c(which(refs$celltypes == 'CM')[1:1000], 
          which(refs$celltypes == 'FB')[1:1000], 
          which(refs$celltypes == 'Macrophages')[1:1000], 
          which(refs$celltypes == 'Neutrophil'))
  
  subref = subset(refs, cells = colnames(refs)[sels])
  
  #rownames(subref) = toupper(rownames(subref))
  
  Idents(subref) = subref$celltype
  
  # Run liana
  # liana_test <- liana_wrap(testdata, method = 'cellphonedb', resource = 'CellPhoneDB')
  source('functions_scRNAseq.R')
  sce <- as.SingleCellExperiment(subref)
  colLabels(sce) = as.factor(sce$celltypes)
  rownames(sce) = toupper(get_geneName(rownames(sce)))
  
  raw = counts(sce)
  
  expr = logcounts(sce)
  #logcounts(sce) = exp(expr)
  ss1 = apply(raw, 1, sum)
  ss2 = apply(expr, 1, sum)
  cc1 = apply(raw, 2, sum)
  cc2 = apply(expr, 2, sum)
  
  sce = sce[which(ss1>10 & ss2> 10), which(cc1>10 & cc2 > 0)]
  
  ## !! not working 
  liana_test <- liana_wrap(sce,  
                           resource = c("Consensus", 'CellPhoneDB',  "CellChatDB",  "CellTalkDB"), 
                           assay.type = "logcounts")
  
  #> Warning in .filter_sce(sce): 3465 genes and/or 0 cells were removed as they had
  #> no counts!
  
  # Liana returns a list of results, each element of which corresponds to a method
  liana_test %>% glimpse
  
  # We can aggregate these results into a tibble with consensus ranks
  liana_test <- liana_test %>%
    liana_aggregate(resource = 'Consensus')
  
  liana_test %>%
    liana_dotplot(source_groups = c("FB"),
                  target_groups = c("CM", "Macrophages", "Neutrophil"),
                  ntop = 20)
  ggsave(filename = paste0(resDir, '/liana_LR_prediction_FB_to_CM_.pdf'), width = 12, height = 8)
  
  
  liana_test %>% 
    filter(source =="CM") %>%
    #top_n(25,  desc(aggregate_rank)) %>%
    liana_dotplot(source_groups = c("CM"),
                  target_groups = c("FB", "Macrophages", "Neutrophil"),
                  ntop = 20)
  ggsave(filename = paste0(resDir, '/liana_LR_prediction_CM_to_FB_prolife.Mphage.pdf'), width = 12, height = 8)
  
  # liana_test %>% 
  #   filter(source =="prolife.Mphage") %>%
  #   top_n(25,  desc(aggregate_rank)) %>%
  #   liana_dotplot(source_groups = c("prolife.Mphage"),
  #                 target_groups = c("FB", "CM"))
  # ggsave(filename = paste0(resDir, '/liana_LR_prediction_prolife.Mphage_to_FB_CM.pdf'), width = 12, height = 8)
  # 
  
  liana_trunc <- liana_test %>%
    # only keep interactions concordant between methods
    filter(aggregate_rank <= 0.01) # this can be FDR-corr if n is too high
  
  heat_freq(liana_trunc)
  
  p <- chord_freq(liana_trunc,
                  source_groups = c("CM", "FB", "Macrophages", "Neutrophil"),
                  target_groups = c("CM", "FB", "Macrophages", "Neutrophil"))
  
  # run any method of choice
  # Load Sce testdata
  # sce <- readRDS(file.path(liana_path , "testdata", "input", "testsce.rds"))
  # 
  # # RUN CPDB alone
  # cpdb_test <- liana_wrap(sce,
  #                         method = 'cellphonedb',
  #                         resource = c('CellPhoneDB'),
  #                         permutation.params = list(nperms=100,
  #                                                   parallelize=FALSE,
  #                                                   workers=4))
  # 
  # # Plot toy results
  # cpdb_test %>%
  #   # filter(pvalue <= 0.05) %>% # only keep interactions with p-val <= 0.05
  #   # invert size (low p-value/high specificity = larger dot size)
  #   # + add a small value to avoid Infinity for 0s
  #   mutate(pvalue = -log10(pvalue + 1e-10)) %>% 
  #   liana_dotplot(source_groups = c("c"),
  #                 target_groups = c("c", "a", "b"),
  #                 specificity = "pvalue",
  #                 magnitude = "lr.mean",
  #                 show_complex = TRUE)
  # 
  # # Run liana re-implementations with the CellPhoneDB resource
  # complex_test <- liana_wrap(testdata,
  #                            method = c('natmi', 'sca', 'logfc'),
  #                            resource = c('CellPhoneDB'))
  # #> Warning in .filter_sce(sce): 3465 genes and/or 0 cells were removed as they had
  # #> no counts!
  # 
  # complex_test %>% liana_aggregate()
  # 
  # 
  
}

run_NicheNet = function() # original code from https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md
{
  library(nichenetr)
  library(Seurat) # please update to Seurat V4
  library(tidyverse)
  
  ## import NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks
  dataPath_nichenet = '../data/NicheNet/'
  ligand_target_matrix = readRDS(paste0(dataPath_nichenet,  "ligand_target_matrix.rds"))
  ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
    
  lr_network = readRDS(paste0(dataPath_nichenet, "lr_network.rds"))
  head(lr_network)
    
  weighted_networks = readRDS(paste0(dataPath_nichenet,  "weighted_networks.rds"))
  head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
  
  head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
  
  ## Read in the expression data of interacting cells
  seuratObj = readRDS(paste0(dataPath_nichenet,  "seuratObj.rds"))
  seuratObj@meta.data %>% head()
  
  sels =c(#which(refs$celltype == 'CM')[1:1000], 
    which(refs$celltype == 'FB')[1:1000], 
    which(refs$celltype == 'prolife.Mphage'))
  subref = subset(refs, cells = colnames(refs)[sels])
  Idents(subref) = subref$celltype
  subref@meta.data$celltype %>% table()
  
  # note that the number of cells of some cell types is very low and should preferably be higher for a real application
  seuratObj@meta.data$celltype %>% table() 
  
  DimPlot(seuratObj, reduction = "tsne")
  
  seuratObj@meta.data$aggregate %>% table()
  
  DimPlot(seuratObj, reduction = "tsne", group.by = "aggregate")
  
  ## Perform the NicheNet analysis
  # indicated cell types should be cell class identities
  # check via: 
  # seuratObj %>% Idents() %>% table()
  nichenet_output = nichenet_seuratobj_aggregate(
    seurat_obj = seuratObj, 
    receiver = "CD8 T", 
    condition_colname = "aggregate", 
    condition_oi = "LCMV", 
    condition_reference = "SS", 
    sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"),
    
    ligand_target_matrix = ligand_target_matrix, 
    lr_network = lr_network, 
    weighted_networks = weighted_networks, 
    organism = "mouse")
  
  ## Interpret the NicheNet analysis output
  nichenet_output$ligand_activities
  
  nichenet_output$top_ligands
  nichenet_output$ligand_expression_dotplot
  
  nichenet_output$ligand_differential_expression_heatmap
  
  nichenet_output$ligand_target_heatmap
  
  nichenet_output$ligand_target_heatmap + 
    scale_fill_gradient2(low = "whitesmoke",  high = "royalblue", breaks = c(0,0.0045,0.009)) + 
    xlab("anti-LCMV response genes in CD8 T cells") + 
    ylab("Prioritized immmune cell ligands")
  
  nichenet_output$ligand_activity_target_heatmap
  
}

# original code from https://saezlab.github.io/liana/articles/liana_nichenet.html
run_liana_nitchenet = function()
{
  
  library(tidyverse)
  library(liana)
  library(nichenetr)
  library(Seurat)
  library(ggrepel)
  library(cowplot)
  options(timeout=180) # required for downloading single-cell expression on slow connection
  
  # single-cell expression matrix described in Puram et al. 2017
  hnscc_expression <-  readRDS(url("https://zenodo.org/record/3260758/files/hnscc_expression.rds"))
  expression <- hnscc_expression$expression
  sample_info <- hnscc_expression$sample_info
  colnames(sample_info) <- make.names(colnames(sample_info))
  
  # filter samples based on vignette's information and add cell type
  tumors_remove <-  c("HN10", "HN", "HN12", "HN13", "HN24", "HN7", "HN8", "HN23")
  sample_info <- sample_info %>%
    subset( !(tumor %in% tumors_remove) & Lymph.node == 0) %>%
    # fix some cell type identity names
    mutate(cell_type = ifelse(classified..as.cancer.cell == 1, "Tumor", non.cancer.cell.type)) %>%
    subset(cell_type %in% c("Tumor", "CAF"))
  
  # cell ID as rownames
  rownames(sample_info) <- sample_info$cell
  
  # subset expression to selected cells
  expression <- expression[sample_info$cell, ]
  
  # model weights
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
  
  # gene set of interest
  geneset_oi <- read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), 
                         col_types = cols(), col_names = "gene") %>%
    pull(gene) %>%
    .[. %in% rownames(ligand_target_matrix)]
  
  # create seurat object
  seurat_object <- Seurat::CreateAssayObject(counts = expm1(t(expression))) %>%
    Seurat::CreateSeuratObject(., meta.data = sample_info) %>%
    Seurat::NormalizeData()
  
  # set cell identity to cell type
  Idents(seurat_object) <- seurat_object@meta.data$cell_type
  
  liana_results <- liana_wrap(seurat_object) %>%
    liana_aggregate()
  
  
  # filter results to cell types of interest
  caf_tumor_results <- liana_results %>%
    subset(source == "CAF" & target == "Tumor")
  
  # filter results to top N interactions
  n <- 50
  top_n_caf_tumor <- caf_tumor_results %>%
    arrange(aggregate_rank) %>%
    slice_head(n = n) %>%
    mutate(id = fct_inorder(paste0(ligand, " -> ", receptor)))
  
  # visualize median rank
  top_n_caf_tumor %>%
    ggplot(aes(y = aggregate_rank, x = id)) +
    geom_bar(stat = "identity") +
    xlab("Interaction") + ylab("LIANA's aggregate rank") +
    theme_cowplot() +
    theme(axis.text.x = element_text(size = 8, angle = 60, hjust = 1, vjust = 1))
  
  
  # get ligands and filter to those included in NicheNet's ligand-target matrix
  ligands <- unique(top_n_caf_tumor$ligand)
  ligands <- ligands[ligands %in% colnames(ligand_target_matrix)]
  ligands
  
  background_genes <- expression[sample_info$cell[sample_info$cell_type == "Tumor"], ] %>%
    apply(2,function(x){10*(2**x - 1)}) %>%
    apply(2,function(x){log2(mean(x) + 1)}) %>%
    .[. >= 4] %>%
    names()
  
  nichenet_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_genes,
    ligand_target_matrix = ligand_target_matrix, potential_ligands = ligands
  )
  
  # prepare data for visualization
  vis_liana_nichenet <- top_n_caf_tumor %>%
    inner_join(nichenet_activities, by = c("ligand" = "test_ligand")) %>%
    arrange(pearson) %>%
    mutate(ligand = fct_inorder(ligand))
  
  # prepare NicheNet figure
  nichenet_scores_plot <- vis_liana_nichenet %>%
    group_by(ligand) %>%
    summarize(pearson = mean(pearson)) %>%
    ggplot(aes(y = ligand, x = pearson)) +
    geom_bar(stat = "identity") +
    ggtitle("NicheNet") +
    xlab("Pearson's score") +
    theme_cowplot() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_line(color = "white"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
  
  # prepare LIANA figure
  liana_receptor_heatmap <- vis_liana_nichenet %>%
    ggplot(aes(y = ligand, x = receptor, fill = aggregate_rank)) +
    geom_tile() +
    theme_cowplot() +
    ggtitle("LIANA") +
    ylab("Ligand") + xlab("Receptor") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_line(colour = "gray", linetype = 2),
          legend.position = "left")
  
  # combine plots
  plot_grid(liana_receptor_heatmap, nichenet_scores_plot,
            align = "h", nrow = 1, rel_widths = c(0.8,0.3))
  
  
}
