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
Run.celltype.deconvolution.RCTD = function(stx = std1, slice = 'adult.day1', Normalization = 'lognormalize')
{
  # slice = 'adult.day4'; Normalization = 'lognormalize'
  
  # import cardiomyocyte reference 
  aa = readRDS(file = paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization.rds'))
  
  if(Normalization == 'SCT'){
    aa = RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
    #saveRDS(aa, file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_SCT_umap.rds'))
  }else{
    aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 20, min.dist = 0.05)
    #saveRDS(aa, file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization_umap.rds'))
  }
  
  aa = aa[ ,!is.na(aa$CellType)]
  Idents(aa) = aa$CellType
  
  aa.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.5)
  
  # gene used for cell type signature
  genes.used =  intersect(aa.markers$gene[!is.na(match(aa.markers$gene, rownames(aa)))], 
                          aa.markers$gene[!is.na(match(aa.markers$gene, rownames(stx)))])
  
  ##########################################
  # prepare reference for RTCD
  # original code from https://raw.githack.com/dmcable/RCTD/dev/vignettes/spatial-transcriptomics.html
  # mainly for coarse CM cells
  ##########################################
  counts <- aa@assays$RNA@counts
  counts = counts[!is.na(match(rownames(counts), genes.used)), ]
  #rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
  meta_data <- aa@meta.data
  cell_types <- meta_data$CellType
  names(cell_types) <- meta_data$CellID # create cell_types named list
  cell_types <- as.factor(cell_types) # convert to factor data type
  nUMI <- meta_data$nCount_RNA 
  names(nUMI) <- meta_data$CellID # create nUMI named list
  
  ### Create the Reference object
  reference <- Reference(counts, cell_types, nUMI)
    
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
  
  myRCTD <- create.RCTD(puck, reference, max_cores = 8, gene_cutoff = 10^-10, fc_cutoff = 0.01, gene_cutoff_reg = 10^-10, 
                        CELL_MIN_INSTANCE = 25)
  
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
  
  saveRDS(myRCTD, file = paste0(RdataDir, 'RCTD_results_ref_Ren2020_', slice, '.rds'))
  
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


