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
  
  #saveRDS(aa, file = paste0(RdataDir, 'Forte2020_logNormalize_allgenes_majorCellTypes_subtypes.rds'))
  
  aa = readRDS(file = paste0(RdataDir, 'Forte2020_logNormalize_allgenes_majorCellTypes_subtypes.rds'))
  
  p0 = DimPlot(aa, reduction = 'umap_0.05', group.by = 'celltype') + ggtitle('Shoval UMAP')
  
  p1 = DimPlot(aa, reduction = 'umap', group.by = 'celltype') + ggtitle(paste0(Normalization, ' - Elad umap'))
  
  p0 + p1 
  ggsave(paste0(resDir, '/Ref_Forte2020_majorCelltypes_Shoval.umap_vs_Elad.umap.', Normalization, '.pdf'), 
         width = 18, height = 8)
  
  
  for(mcells in c('all', 'FB', 'EC', 'immune'))
  {
    # mcells = 'immune'
    
    if(mcells == 'all'){
      ax = aa
      
      Idents(ax) = ax$celltype
      
    }else{
      ax = subset(aa, cells = colnames(aa)[which(aa$celltype == mcells & aa$subtype != 'B' & aa$subtype != 'GN')] )
      table(ax$celltype)
      table(ax$subtype)
      
      ax <- FindVariableFeatures(ax, selection.method = "vst", nfeatures = 3000)
      ax <- ScaleData(ax, features = rownames(ax))
      
      ax <- RunPCA(ax, verbose = FALSE, weight.by.var = TRUE)
      ElbowPlot(ax, ndims = 30)
      
      # UMAP to visualize subtypes
      ax <- RunUMAP(ax, dims = 1:20, n.neighbors = 30, min.dist = 0.05, n_threads = 6)
      
      DimPlot(ax, reduction = 'umap', group.by = 'subtype') + 
      ggtitle(paste0(mcells, '-', ' cells UMAP (', Normalization, ' nfeature = 3000, ndim=30, neighbors=30, mdist=0.05)'))
      
      ggsave(paste0(resDir, '/Ref_Forte2020_UMAP_', mcells, '_excluding B and GN_subcelltypes.pdf'), 
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

}

##########################################
# run cell type deconvolution with RCTD
##########################################
Run.celltype.deconvolution.RCTD = function(stx = std1, slice = 'adult.day1', Normalization = 'lognormalize')
{
  # slice = 'adult.day4'; Normalization = 'lognormalize';
  refs = readRDS(file = paste0(RdataDir, 
                               'SeuratObj_adultMiceHeart_refCombine_Forte2020.nonCM_Ren2020CM_cleanAnnot_logNormalize_v1.rds'))
  
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
  # prepare reference for RTCD
  # original code from https://raw.githack.com/dmcable/RCTD/dev/vignettes/spatial-transcriptomics.html
  # mainly for coarse CM cells
  ##########################################
  E.corrected <- refs@assays$integrated@data
  E.corrected = expm1(E.corrected) # linear scale of corrected gene expression
  
  #counts = counts[!is.na(match(rownames(counts), genes.used)), ]
  #rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
  cell_types <- refs@meta.data$celltype
  names(cell_types) <- colnames(refs) # create cell_types named list
  cell_types <- as.factor(cell_types) # convert to factor data type
  
  #nUMI <- meta_data$nCount_RNA 
  #names(nUMI) <- meta_data$CellID # create nUMI named list
  
  ### Create the Reference object
  reference <- Reference(E.corrected, cell_types, nUMI = NULL, require_int = FALSE)
    
  ## Examine reference object (optional)
  print(dim(reference@counts)) #observe Digital Gene Expression matrix
  #> [1] 384 475
  table(reference@cell_types) #number of occurences for each cell type
  
  ##########################################
  # prepare ST data for RTCD
  # original code from https://raw.githack.com/dmcable/RCTD/dev/vignettes/spatial-transcriptomics.html
  ##########################################
  counts <- stx@assays$Spatial@counts
  #counts = counts[!is.na(match(rownames(counts), genes.used)), ]
  coords <- eval(parse(text = paste0('stx@images$',  slice, '@coordinates')))
  coords = coords[, c(4, 5)]
  #rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
  #rownames(coords) <- coords$barcodes; 
  #coords$barcodes <- NULL # Move barcodes to rownames
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
  
  myRCTD <- create.RCTD(puck, reference, max_cores = 16, gene_cutoff = 0.000125, fc_cutoff = 0.5, 
                        gene_cutoff_reg = 2e-04,  fc_cutoff_reg = 0.75,
                        UMI_min = 100, UMI_max = 2e+07, 
                        CELL_MIN_INSTANCE = 25)
  
  
  tic()
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
  saveRDS(myRCTD, file = paste0(RdataDir, 'RCTD_results_refCombined_Forte2020.Ren2020_', slice, '.rds'))
  toc()
  
  results <- myRCTD@results
  
  # normalize the cell type proportions to sum to 1.
  norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
  cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
  
  spatialRNA <- myRCTD@spatialRNA
  resultsdir <- resDir
  
  # make the plots 
  # Plots the confident weights for each cell type as in full_mode (saved as 
  # 'results/cell_type_weights_unthreshold.pdf')
  plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
  
  # Plots all weights for each cell type as in full_mode. (saved as 
  # 'results/cell_type_weights.pdf')
  plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
  # Plots the weights for each cell type as in doublet_mode. (saved as 
  # 'results/cell_type_weights_doublets.pdf')
  
  plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet, 
                       results$results_df) 
  # Plots the number of confident pixels of each cell type in 'full_mode'. (saved as 
  # 'results/cell_type_occur.pdf')
  plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
  
  # makes a map of all cell types, (saved as 
  # 'results/all_cell_types.pdf')
  plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir) 
  
  ##########################################
  # save and assign the cell type 
  ##########################################
  cts = results$results_df
  cts$cellType = NA
  cts$cellType[which(cts$spot_class == 'singlet')] = as.character(cts$first_type[which(cts$spot_class == 'singlet')])
  cts$cellType[which(cts$spot_class == 'doublet_certain')] = 
    paste0(cts$first_type[which(cts$spot_class == 'doublet_certain')], 
           '_', cts$second_type[which(cts$spot_class == 'doublet_certain')])
  
  stx$celltype = cts$cellType[match(colnames(stx), rownames(cts))]
  
  SpatialDimPlot(stx, group.by = 'celltype', images = slice, stroke = 0, interactive = FALSE)
  
  # remove the CM annotated spots from Visium and prepare for the second-round RCTD
  
  
  ##########################################
  # finer annotation for non-cardiomyocyte using reference Forte2020
  ##########################################
  
  
  
  
}


