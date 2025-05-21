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
make_SeuratObj_visium = function(topdir = './', saveDir = './results', changeGeneName = TRUE, 
                                 keyname = 'slice1', 
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
# Section : cell type deconvolution analysis 
# 
########################################################
########################################################
assess_sutypes_similarity_for_RCTD = function(st, refs, 
                                              require_int_SpatialRNA = FALSE,
                                              RCTD_out)
{
  library(spacexr)
  require(Matrix)
  
  ## remove the subtypes in Atria
  sutypes_atria = c('EC_WNT4', 'EC_CEMIP', 'EC_LHX6', 'FB_VWA2', 'FB_TNXB', 'CM_PM_HCN4', 'CM_OFT',
                    'CM_Atria', 'CM_Atria_Tagln', 'Neuronal', 'Proliferating_RBC'
                    )
  
  ## remove the proliferating cells
  #csc[which(rownames(csc) == 'EC_IS_Prol'), ] = NA
  sutypes_atria = c('EC_WNT4', 'EC_CEMIP', 'EC_LHX6', 'FB_VWA2', 'FB_TNXB', 'CM_PM_HCN4', 'CM_OFT',
                    'CM_Atria', 'CM_Atria_Tagln', 'Neuronal', 'Proliferating_RBC',
                    'Proliferating_Megakeryocytes', 'Mo.Macs_Prol', 'EC_Prol', 'CM_Prol_1',
                    'CM_Prol_2', 'CM_Prol_3', 'FB_Prol', 'B_cells_Prol',
                    'Neu_IL1R1', 'EC_IS_Prol'
  )
  
  
  mm = match(refs$celltype_toUse, sutypes_atria)
  refs = subset(refs, cells = colnames(refs)[which(is.na(mm))])
  
  cat('-- prepare global reference for all conditions --\n')
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
  
  
  cat('-- check visium conditions -- \n')
  st$condition = droplevels(factor(st$condition))
  print(table(st$condition))
  cc = names(table(st$condition))
  
  n = 1
  slice = cc[n]
  stx = st[, which(st$condition == slice)]
  
  ##########################################
  # prepare ST data for RTCD
  # original code from https://raw.githack.com/dmcable/RCTD/dev/vignettes/spatial-transcriptomics.html
  ##########################################
  counts = GetAssayData(object = stx, slot = "counts")
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
  
  # make RCTD object
  myRCTD <- create.RCTD(puck, reference, max_cores = 1, 
                        gene_cutoff = 0.000125, fc_cutoff = 0.5, 
                        gene_cutoff_reg = 2e-04,  
                        fc_cutoff_reg = 0.75,
                        UMI_min = 20, UMI_max = 2e+07, 
                        CELL_MIN_INSTANCE = 30)
  
  markers = myRCTD@cell_type_info$info[[1]]
  ggs = myRCTD@internal_vars$gene_list_reg
  
  mm = match(ggs, rownames(markers))
  
  markers = markers[mm, ]
  sampleDists <- dist(t(markers))
  
  library("RColorBrewer")
  library(pheatmap)
  sampleDistMatrix <- as.matrix(sampleDists, method = "euclidean")
  
  #rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")) )(255)
  
  pdfname = paste0(RCTD_out, '/Similarity_subtypes_atria.proliferatingExcluded_v2.pdf')
  pdf(pdfname, width=10, height = 7)
  
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  
  dev.off()
  
  
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
                                           mapping_ref = 'condition',
                                           condition.specific.ref = FALSE,
                                           condition.specific_celltypes = NULL,
                                           RCTD_mode = 'doublet',
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
  if(!condition.specific.ref){
    cat('-- prepare global reference for all conditions --\n')
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
    
  }
  
  ##########################################
  # loop over all conditions of st for now
  ##########################################
  cat('-- check visium conditions -- \n')
  st$condition = droplevels(factor(st$condition))
  print(table(st$condition))
  cc = names(table(st$condition))
  
  for(n in 1:length(cc))
  #for(n in c(1, 2, 4))
  {
    # n = 1
    cat('slice -- ', cc[n], '\n')
    
    ## prepare condition-specific refs
    if(condition.specific.ref | !is.null(condition.specific_celltypes)){
      cat('-- prepare refs for ', cc[n], '\n')
      
      if(is.null(condition.specific_celltypes)){
        cat('subset refs based on the condition -- ', cc[n], '\n')
        
        if(mapping_ref == 'condition'){
          refs_sub = subset(refs, condition == cc[n])
        }else{
          tt = paste0(unlist(strsplit(as.character(cc[n]), '_'))[1:2], collapse = '_')
          refs_sub = subset(refs, condition == tt)
        }
        
        celltypes_nbs = table(refs_sub$celltype_toUse)
        celltypes_nbs = celltypes_nbs[which(celltypes_nbs>30)]
        refs_sub = subset(refs_sub, cells = colnames(refs_sub)[!is.na(match(refs_sub$celltype_toUse, 
                                                                            names(celltypes_nbs)))])
      }else{
        cat('subset refs based on the  -- condition.specific_celltypes ', cc[n], '\n')
        
        if(mapping_ref == 'condition'){
          ii_col = which(colnames(condition.specific_celltypes) == cc[n])
        }else{
          tt = paste0(unlist(strsplit(as.character(cc[n]), '_'))[1:2], collapse = '_')
          ii_col = which(colnames(condition.specific_celltypes) == tt)
        }
        
        if(length(ii_col) != 1) stop(paste0('Error -- no celltypes found for ', cc[n]))
        
        celltypes_sel = 
          rownames(condition.specific_celltypes)[which(!is.na(condition.specific_celltypes[, ii_col]))]
        
        refs_sub = subset(refs, cells = colnames(refs)[!is.na(match(refs$celltype_toUse, 
                                                                            celltypes_sel))])
        
      }
      
      
      E_refs = GetAssayData(object = refs_sub, slot = "data")
      E_refs = expm1(E_refs) # linear scale of corrected gene expression
      
      ### Create the Reference object for major cell types
      cell_types <- refs_sub@meta.data$celltype_toUse
      names(cell_types) <- colnames(refs_sub) # create cell_types named list
      cell_types <- as.factor(cell_types) # convert to factor data type
      reference <- Reference(E_refs, cell_types, nUMI = NULL, require_int = FALSE)
      
      rm(E_refs, cell_types, refs_sub)
      
      ## Examine reference object (optional)
      print(dim(reference@counts)) #observe Digital Gene Expression matrix
      #> [1] 384 475
      table(reference@cell_types) #number of occurences for each cell type
      
      
    }
    
    slice = cc[n]
    stx = st[, which(st$condition == slice)]
    
    resultsdir <- paste0(RCTD_out, '/', slice)
    system(paste0('mkdir -p ', resultsdir))
    
    ##########################################
    # prepare ST data for RTCD
    # original code from https://raw.githack.com/dmcable/RCTD/dev/vignettes/spatial-transcriptomics.html
    ##########################################
    counts = GetAssayData(object = stx, slot = "counts")
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
    cat('----RCTD running mode : ', RCTD_mode, ' ----\n')
    myRCTD <- run.RCTD(myRCTD, doublet_mode = RCTD_mode)
    saveRDS(myRCTD, file = paste0(resultsdir, '/RCTD_out_doubletMode_', slice, '.rds'))
    toc()
    
    ##########################################
    # check the result
    ##########################################
    # myRCTD = readRDS(file = paste0(resultsdir, '/RCTD_out_doubletMode_', slice, '.rds'))
    myRCTD = readRDS(file = paste0(resultsdir, '/RCTD_out_doubletMode_', slice, '.rds'))
    results <- myRCTD@results
    
    cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
    spatialRNA <- myRCTD@spatialRNA
    
    # extract the weights
    if(RCTD_mode == 'doublet'){
      # normalize the cell type proportions to sum to 1
      norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
      
    }else{
      #cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
      norm_weights = matrix(0, ncol = length(cell_type_names), nrow = length(results))
      
      colnames(norm_weights) = cell_type_names
      for(jj in 1:nrow(norm_weights))
      {
        wts = results[[jj]]
        ii_w = match(wts$cell_type_list, colnames(norm_weights))
        norm_weights[jj, ii_w] = wts$sub_weights
      }
      
      rownames(norm_weights) = colnames(stx)
      
    }
    
    
    # make the plots 
    if(plot.RCTD.summary & RCTD_mode == 'doublet'){
      
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

plot.RCTD.results = function(st, 
                             species = 'axolotl',
                             RCTD_out = '../results/RCTD_out',
                             RCTD_mode = 'multi',
                             plot.RCTD.summary = FALSE,
                             PLOT.scatterpie = TRUE)
{
  library(spacexr)
  require(Matrix)
  require(scatterpie)
  library(tidyr)
  library(dplyr)
  require(cowplot)
  library(RColorBrewer)
  
  cat('-- check visium conditions -- \n')
  print(table(st$condition))
  nb_cells = table(st$condition)
  cc = names(nb_cells[which(nb_cells>0)])
  
  #RCTD_out = paste0(resDir, '/RCTD_coarse_out_v1')
  #RCTD_out = paste0(resDir, '/RCTD_subtype_out')
  #RCTD_out = paste0(resDir, '/RCTD_subtype_out_v3.5')
  #RCTD_out = paste0(resDir, '/RCTD_subtype_out_v4_FBsubtypes')
  
  cat('-- RCTD output folder : \n -- ', RCTD_out, '\n')
  
  for(n in 1:length(cc))
  #for(n in c(1, 2, 4))
  {
    # n = 3
    cat('slice -- ', cc[n], '\n')
    slice = cc[n]
    #stx = st[, which(st$condition == slice)]
    stx = subset(st, condition == slice)
    #stx@images = stx@images[[n]]
    
    resultsdir <- paste0(RCTD_out, '/', slice)
    system(paste0('mkdir -p ', resultsdir))
    
    ## you may change this to a more accessible directory on your computer.
    #resultsdir <- paste0(RCTD_out, '/', slice, '_Plots') 
    #dir.create(resultsdir)
    
    ##########################################
    # check the result
    # it seems that the full mode works better than doublet mode
    ##########################################
    file = paste0(resultsdir, '/RCTD_out_doubletMode_', slice, '.rds')
    
    if(file.exists(file)){
      myRCTD = readRDS(file = file)
    }else{
      slice_name = paste0(unlist(strsplit(as.character(slice), '_'))[1:2], collapse = '_')
      file = paste0(resultsdir, '/', slice_name,  '/RCTD_out_doubletMode_', 
                    slice_name, '.rds')
      
      if(file.exists(file)){
        myRCTD = readRDS(file = file)
      }else{
        cat('-- RCTD output for', slice, ' not found -- \n')
        next()
      }
    }
   
    results <- myRCTD@results
    
    cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
    spatialRNA <- myRCTD@spatialRNA
    
    # extract the weights
    if(RCTD_mode == 'doublet'){
      # normalize the cell type proportions to sum to 1
      norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
      
    }else{
      #cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
      norm_weights = matrix(0, ncol = length(cell_type_names), nrow = length(results))
      
      colnames(norm_weights) = cell_type_names
      for(jj in 1:nrow(norm_weights))
      {
        wts = results[[jj]]
        ii_w = match(wts$cell_type_list, colnames(norm_weights))
        norm_weights[jj, ii_w] = wts$sub_weights
      }
      
      rownames(norm_weights) = colnames(stx)
      
    }
    
    # normalize the cell type proportions to sum to 1.
    #norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
    #norm_weights = normalize_weights(results$weights) 
    #cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
    #spatialRNA <- myRCTD@spatialRNA
    
    stx$condition = as.factor(stx$condition)
    stx$condition = droplevels(stx$condition)
    DefaultAssay(stx) = 'SCT'
    
    if(species == 'axolotl'){
      SpatialFeaturePlot(stx, images =  slice, #slot = 'data',
                         features = "NPPA-AMEX60DD051098") +
        ggtitle('NPPA-AMEX60DD051098')
    }
    if(species == 'mouse_adult'|species == 'mouse_neonatal'){
      SpatialFeaturePlot(stx, images =  slice, #slot = 'data',
                         features = "Nppa") +
        ggtitle('Nppa')
    }
   
    #coord_flip() + 
    #ggplot2::scale_y_reverse() # rotate the coordinates to have the same as RCTD 
    ggsave(filename =  paste0(RCTD_out, "/Spatial_patterning_", slice, "_NPPA.pdf"), 
           width = 12, height = 8)
    
    # make the plots 
    if(plot.RCTD.summary & RCTD_mode == 'doublet'){
      # stx$condition = droplevels(stx$condition)
      # DefaultAssay(stx) = 'SCT'
      # 
      # SpatialFeaturePlot(stx, images =  slice, slot = 'data',
      #                    features = "LOC115076416-AMEX60DD051098", image.alpha = 0.5) +
      #   ggtitle('NPPA-AMEX60DD051098') +
      #   coord_flip() + 
      #   ggplot2::scale_y_reverse() # rotate the coordinates to have the same as RCTD 
      # 
      # ggsave(filename =  paste0(resultsdir, "/Spatial_patterning_", slice, "_NPPA.pdf"), width = 12, height = 8)
      
      # Plots the confident weights for each cell type as in full_mode 
      # (saved as 'results/cell_type_weights_unthreshold.pdf')
      plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
      
      # Plots all weights for each cell type as in full_mode. (saved as 'resultrs/cell_type_weights.pdf')
      plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights) 
      
      # Plots the weights for each cell type as in doublet_mode: doublet model
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
      # set color vector
      getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
      
      ## Preprocess coordinates 
      spatial_coord <-  spatialRNA@coords %>%
        tibble::rownames_to_column("ID")
      
      weights = norm_weights[match(spatial_coord$ID, rownames(norm_weights)), ]
      
      # set the color vectors for all cell types
      # n <- 60
      # qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
      # col_vector  <-  unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      # 
      # use a panel of colors from https://gotellilab.github.io/GotelliLabMeetingHacks/NickGotelli/ColorPalettes.html
      #tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77",
      #                   "#661100", "#CC6677", "#882255", "#AA4499")
      
      Process_celltype_weight_singlet = FALSE 
      if(Process_celltype_weight_singlet){ ### if use the doublet mode 
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
      }else{ # clean the weight: force to be 0 if lower than certain threshold, e.g. 0.05 or 0.1
        cutoff_weight = 0.01
        weights_processed = as.matrix(weights)
        weights_processed[which(weights_processed < cutoff_weight)] = 0
        weights_processed = normalize_weights(weights_processed)
        
        weights = weights_processed;
        rm(weights_processed)
        weights = weights[, apply(weights, 2, sum)>0]
        
      }
      
      cell_types_plt = unique(colnames(weights))
      col_vector = getPalette(length(cell_types_plt))
      
      #names(col_vector) = cell_types_plt
      #col_ct = readRDS(file = paste0(RdataDir, 'main_celltype_colors_4RCTD.rds'))
      col_ct <- col_vector[seq_len(length(cell_types_plt))]
      #names(col_ct) = cell_types_plt
      #saveRDS(col_ct, file = paste0(RdataDir, 'main_celltype_colors_4RCTD.rds'))
      
      #plot(1:length(col_ct), 1:length(col_ct), col = col_ct)
      #text(1:length(col_ct), 1:length(col_ct), names(col_ct), offset = 2, adj = 0.5,  col = col_ct)
      
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
      
      ggsave(paste0(RCTD_out, '/RCTD_scatterpie_', slice, '.pdf'), width = 16, height = 10)
      
      
      pdfname = paste0(RCTD_out, '/Spatial_distribution_each_celltype_', slice, '.pdf')
      pdf(pdfname, width=16, height = 10)
      
      for(m in 1:length(cell_types_plt)){
        # m = 8
        sub.spatial = spatial_coord[,c(1:3, (3+m))]
        colnames(sub.spatial)[ncol(sub.spatial)] = 'celltype'
        sub.spatial$celltype = as.numeric(as.character(sub.spatial$celltype))
        sub.spatial = data.frame(sub.spatial)
        
        p1 = ggplot(data = sub.spatial, aes(x = x, y = y)) +
          geom_point(aes(size = celltype), color = col_ct[m]) +  
          #scale_size_continuous(range = c(0, 0.7)) +
          #geom_raster() +
          ggplot2::labs(title = cell_types_plt[m]) +
          coord_flip() + 
          #ggplot2::scale_y_reverse() +
          ggplot2::scale_x_reverse()+
          ggplot2::theme_void() + 
          ggplot2::theme(
            # plot.background = element_rect(fill = "#FFFFFF"),
            # panel.background = element_blank(),
            # plot.margin = margin(20, 20, 20, 20),
            plot.title = ggplot2::element_text(hjust = 0.5, size = 20)) +
          ggplot2::guides(fill = guide_legend(ncol = 1)) + 
          scale_size(limits = c(0,1), breaks = seq(0, 1, by = 0.2))
        
        plot(p1)
        
        make_special_plot = FALSE
        if(make_special_plot){
          p1
          
          ggsave(paste0(resDir, '/RCTD_res_', slice, '_', cell_types_plt[m],  '.pdf'),
                 width = 10, height = 8)
          
        }
        
      }
      
      dev.off()
      
    }
    
  }
  
}

Run_imputation_snRNAseq_visium = function(stx, refs, 
                                          RCTD_out = paste0('../results/visium_axolotl_R12830_resequenced_20220308/',
                                                            'RCTD_subtype_out_42subtypes_ref.time.specific_v4.3'), 
                                          slice, normalized_weights = TRUE)
{
  library(spacexr)
  require(Matrix)
  require(scatterpie)
  require(cowplot)
  library(RColorBrewer)
  ## preapre the reference 
  cat('-- prepare global reference --\n')
  
  E_refs = GetAssayData(object = refs, slot = "data")
  E_refs = expm1(E_refs) # linear scale of corrected gene expression
  
  ### summarize the expression for each subtypes
  cell_types_ref <- unique(refs@meta.data$celltypes)
  
  X = matrix(NA, ncol = length(cell_types_ref), nrow = nrow(E_refs))
  colnames(X) = cell_types_ref
  rownames(X) = rownames(E_refs)
  
  for(n in 1:ncol(X))
  {
    cat(n, '--', colnames(X)[n], '\n')
    X[,n] = rowMeans(E_refs[,which(refs@meta.data$celltypes == colnames(X)[n])])
  }
  
  cell_types_ref = gsub('CM_Robo2', 'CM_ven_Robo2', cell_types_ref)
  cell_types_ref = gsub('CM_Cav3.1', 'CM_ven_Cav3_1', cell_types_ref)
  colnames(X) = cell_types_ref
  
  ## import the inferred cell type proportions
  #RCTD_out = paste0('../results/visium_axolotl_R12830_resequenced_20220308/',
  #                  'RCTD_subtype_out_42subtypes_ref.time.specific_v4.3')
  
  resultsdir <- paste0(RCTD_out)
 
  myRCTD = readRDS(file = paste0(resultsdir, '/RCTD_out_doubletMode_', slice, '.rds'))
  results <- myRCTD@results
  
  # normalize the cell type proportions to sum to 1.
  #norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
  norm_weights = normalize_weights(results$weights) 
  cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
  weights = results$weights
  
  mm = match(cell_type_names, cell_types_ref)
  
  X1 = X[,mm]
  if(normalized_weights){
    y = as.matrix(X1) %*% t(as.matrix(norm_weights))
  }else{
    y = as.matrix(X1) %*% t(as.matrix(weights))
  }
  
  # y = log2(y+1)
  
  cell.shared = intersect(colnames(stx), colnames(y))
  if(length(cell.shared) < ncol(stx)){
    stx = subset(stx, cells = cell.shared)
  }
  y = y[, match(cell.shared, colnames(y))]
  
  stx[['imputated']] = CreateAssayObject(data = y)
  
  return(stx)
  
}
########################################################
########################################################
# Section : define boarder zone, injury zone and remote zone
# 1) manually with SPATA2
########################################################
########################################################
manual_selection_spots_image_Spata2 = function(st,
                                               outDir = paste0(resDir, '/manual_borderZones_spata2/')
                                               )
{
  #dyn.load("/software/f2021/software/proj/7.2.1-gcccore-10.2.0/lib/libproj.so")
  #dyn.load("/software/f2021/software/gdal/3.2.1-foss-2020b/lib/libgdal.so") 
  # outDir = '/groups/tanaka/Collaborations/Jingkui-Elad/visium_axolotl_reseq/spata2_manual_regions/'
  # st = readRDS(file = paste0(outDir, 'axolotl_visium_newRep.rds'))
  
  # outDir = paste0(resDir, '/manual_Ventricle_spata2/')
  require(SPATA2)
  
  # slice = cc[n]
  cat('-- check visium conditions -- \n')
  print(table(st$condition))
  cc = names(table(st$condition))
  system(paste0('mkdir -p ', outDir))
  
  for(n in 1:length(cc))
  {
    # n = 4
    cat('slice -- ', cc[n], '\n')
    slice = cc[n]
    aa = st[, which(st$condition == slice)]
    DefaultAssay(aa) = 'Spatial'
    #aa@images = aa@images$Amex_d4
    
    resultsdir <- paste0(outDir, '', slice)
    system(paste0('mkdir -p ', resultsdir))
    
    # reconstruct Seurat object with associated image; 
    # because aa still have multiple images, doesnot work to transform to spata 
    srat <- CreateSeuratObject(counts = GetAssayData(aa, slot = 'counts'), assay = "Spatial", 
                               meta.data = aa@meta.data) 
    #image <- image[Cells(x = srat)] # filter image by the spots
    #DefaultAssay(object = image) <- "Spatial" # set default assay
    srat[[slice]] <- eval(parse(text = paste0('aa@images$',  slice)))
    
    aa = srat; rm(srat)
    aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
    aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 1000)
    
    aa <- ScaleData(aa, features = rownames(aa), assay = 'Spatial')
    aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)
    aa = RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
    # spata_obj = asSPATA2(aa, 
    #                      sample_name = slice, 
    #                      platform = "Undefined",
    #                      assay_name = 'Spatial', 
    #                      img_name = slice,
    #                      img_scale_fct = "lowres", 
    #                      assay_modality = 'gene')
    spata_obj = transformSeuratToSpata(aa, 
                                       sample_name = slice, method = 'spatial',
                           assay_name = 'Spatial',
                           assay_slot = 'scale.data',
                           image_name = slice)
    
    #setActiveExpressionMatrix(spata_obj, 'data')
    spata_obj <- createSegmentation(object = spata_obj)
    
    
    
    plotSegmentation(object = spata_obj, pt_size = 1.9) +
      ggplot2::scale_y_reverse()
    
    saveRDS(spata_obj, file = paste0(outDir, '/segemented_spata2_', slice, '.rds'))
    
    if(n == 1) { st$seg2 = NA }
    
    st$seg2[which(st$condition == slice)] = NA
    st$seg2[match(getSegmentDf(spata_obj, segment_names = c('Ventricle'))$barcodes, 
                  colnames(st))] = 'Ventricle'
    SpatialDimPlot(st, group.by = 'seg2', ncol = 4)
    
    #spata_obj = readRDS(file = paste0(outDir, '/segemented_spata2_', slice, '.rds'))
    #st$seg2[match(getSegmentDf(spata_obj, segment_names = c('Ventricle'))$barcodes, 
    #              colnames(st))] = 'Ventricle'
    #aa$segmentation[match(getSegmentDf(spata_obj, segment_names = c('border_zone'))$barcodes, 
    #                      colnames(aa))] = 'border_zone'
    #aa$segmentation[match(getSegmentDf(spata_obj, segment_names = c('remote_zone1'))$barcodes, 
    #                      colnames(aa))] = 'remote_zone1'
    #aa$segmentation[match(getSegmentDf(spata_obj, segment_names = c('remote_zone2'))$barcodes, 
    #                      colnames(aa))] = 'remote_zone2'
    
    saveRDS(st, file = paste0(RdataDir, 'seuratObject_allVisiusmst_',
                              'filtered.spots_time_conditions_manualSegmentation_ventricleRegions', 
                              version.analysis, '.rds'))
    
  }
}


########################################################
# 2) run bayesSpace to detect border zone
# original code from https://bioconductor.org/packages/release/bioc/vignettes/BayesSpace/inst/doc/BayesSpace.html
########################################################
run_bayesSpace = function(st, 
                          outDir = paste0(resDir, '/bayesSpace/'),
                          Find.top.markers.for.spatial.clusters = FALSE,
                          Run.bayesSpace.enhanced.clustering = FALSE)
{
  ## aa is a seurat objet with one slice / one image
  require(SingleCellExperiment)
  library(BayesSpace)
  library(ggplot2)
  library(patchwork)
  library(scran)
  library(scuttle)
  require(RColorBrewer)
  library(mclust) ## Here we run mclust externally so the random seeding is consistent with ## original analyses
  
  system(paste0('mkdir -p ', outDir))
  
  cat('-- check visium conditions -- \n')
  print(table(st$condition))
  cc = names(table(st$condition))
  
  for(n in 1:length(cc))
  {
    # n = 1
    
    cat('slice -- ', cc[n], '\n')
    slice = cc[n]
    stx = st[, which(st$condition == slice)]
    DefaultAssay(stx) = 'Spatial'
    
    #slice = names(table(stx$condition))
    scc <- as.SingleCellExperiment(stx, assay = 'Spatial')
    coords <- eval(parse(text = paste0('stx@images$',  slice, '@coordinates')))
    
    resultsdir <- paste0(outDir, '', slice)
    system(paste0('mkdir -p ', resultsdir))
    
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
    set.seed(102)
    scc <- spatialPreprocess(scc, platform="Visium", skip.PCA=TRUE)
    
    cat('---- select nb of clusters ----\n ')
    scc <- qTune(scc, qs=seq(4, 15))
    qPlot(scc)
    
    ggsave(filename =  paste0(resultsdir, "/SpatialCluster_nb.clusters.selection_", slice, ".pdf"), 
           width = 12, height = 8)
    
    
    # sptial clustering 
    d <- 15  # Number of PCs recommended by the paper
    
    for(q in c(5:15)){ # enumerate the number of clusters
      
      # q <- 13  # Number of clusters
      cat('---- nb of clusters : ', q, '----\n')
      #Y <- reducedDim(scc, "PCA")[, seq_len(d)]
      #set.seed(101)
      #init <- Mclust(Y, q, "EEE", verbose=FALSE)$classification
      
      ## Run BayesSpace clustering
      set.seed(100)
      #scc <- spatialCluster(scc, q=q, d=d, platform='Visium', init=init,
      #                      nrep=10000, gamma=3)
      scc = spatialCluster(scc, q=q, 
                           use.dimred = "PCA",
                           d=d, 
                           platform="Visium",
                           init.method="mclust", 
                           model="t", 
                           gamma=3,
                           nrep = 20000,
                           burn.in = 1000,
                           save.chain=FALSE)
      
      # scc = readRDS(file = paste0(RdataDir, "/BayesSpace_SpatialSlustered_", species, '_', slice,  "_with.", q, "clusters.rds"))
      # p1 = SpatialDimPlot(aa, group.by = 'spatial_domain_manual', images = slice)
      
      palette =  c(brewer.pal(name="Paired", n = 12), brewer.pal(name="Set3", n = 12))[1:q]
     
      p2 = clusterPlot(scc, palette=palette, size=0.1) +
        labs(title= paste0("spatial domains : ", q))  +
        coord_flip() + 
        ggplot2::scale_x_reverse() 
      #ggplot2::scale_x_reverse()  # flip first and reverse x to match seurat Spatial plots
      
      p2  
      
      ggsave(filename =  paste0(resultsdir, "/BayesSpace_SpatialSlustered_", slice, "_cluster_", 
                                q, ".pdf"), width = 16, height = 8)
      saveRDS(scc, file = paste0(resultsdir, "/BayesSpace_SpatialClustered_", slice,  "_with_clusters_", q, ".rds"))
      
      #aa$spatial_domain_bayeSpace = NA
      #aa$spatial_domain_bayeSpace[match(colnames(scc), colnames(aa))] = scc$spatial.cluster
      
      #SpatialDimPlot(aa, group.by = 'spatial_domain_bayeSpace', images = slice)
      
      #saveRDS(aa, file = paste0(paste0(RdataDir, "/SeuratObj_spatialDomain_BayesSpace_", 
      #                                 species, '_', slice,  "_with.", q, "clusters.rds")))
      
    }
    
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
    
    # this takes > 1 hour at least long time
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
# Section : Analyze the cell type proximity
# original code from https://github.com/madhavmantri/chicken_heart/blob/master/scripts/anchor_integeration.R
# and the paper (https://www.nature.com/articles/s41467-021-21892-z#data-availability)
# 
########################################################
########################################################
run_neighborhood_analysis = function(st, 
                              outDir = '../results/neighborhood_test/',
                              RCTD_out = '../results/visium_axolotl_R12830_resequenced_20220308/RCTD_subtype_out_v2',
                              background = 'remote'
                              )
{
  library(corrplot)
  require(RColorBrewer)
  require(spacexr)
  library(ggplot2)
  library(SPOTlight)
  library(SingleCellExperiment)
  library(SpatialExperiment)
  library(scater)
  library(scran)
  
  system(paste0('mkdir -p ', outDir))
  
  cat('-- check visium conditions -- \n')
  st$condition = droplevels(as.factor(st$condition))
  print(table(st$condition))
  cc = names(table(st$condition))
  
  Idents(st) = st$condition
  
  for(n in 1:length(cc))
  {
    # n = 2
    cat('slice -- ', cc[n], '\n')
    slice = cc[n]
    stx = subset(st, condition == slice)
    DefaultAssay(stx) = 'Spatial'
    
    stx$segmentation = as.factor(stx$segmentation)
    Idents(stx) = stx$segmentation
    SpatialDimPlot(stx, images = slice)
    
    print(table(stx$segmentation))
    
    # import the cell type deconvolution result
    myRCTD = readRDS(file = paste0(RCTD_out, '/', slice,  '/RCTD_out_doubletMode_', slice, '.rds'))
    results <- myRCTD@results
    
    norm_weights = normalize_weights(results$weights) 
    #cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
    
    SPOTlight::plotCorrelationMatrix(as.matrix(norm_weights))
    ggsave(paste0(outDir, '/CorrelationMatrix_allRegions', slice, '.pdf'), width = 10, height = 10)
    
    #plotInteractions(as.matrix(norm_weights), which = "heatmap", metric = "jaccard")
    #plotInteractions(as.matrix(norm_weights), which = "network")
    
    
    # use a threshold to binarize weights
    weights = norm_weights >= 0.1 
    index_border = match(colnames(stx)[which(stx$segmentation == 'BZ')], rownames(weights))
    
    SPOTlight::plotCorrelationMatrix(as.matrix(norm_weights[index_border, ]))
    ggsave(paste0(outDir, '/CorrelationMatrix_BZ_', slice, '.pdf'), width = 10, height = 10)
    
    #plotInteractions(as.matrix(norm_weights[index_border, ]), which = "network")
    #ggsave(paste0(outDir, '/InteractionNetwork_BZ_', slice, '.pdf'), width = 14, height = 14)
    
    if(background == 'remote'){
      index_bg = match(colnames(stx)[grep('Remote|remote', stx$segmentation)], rownames(weights))
    }
    if(background == 'remote.others'){
      index_bg = match(colnames(stx)[unique(c(grep('Remote|remote', stx$segmentation), 
                                              which(is.na(stx$segmentation))))], rownames(weights))
    }
    
    if(background == 'whole.slice'){
      index_bg = c(1:nrow(weights))
    }
    
    cat(length(index_border), ' spots in BZ\n')
    cat(length(index_bg), ' spots in background \n')
    
    border = weights[index_border, ]
    bg = weights[index_bg, ]
    
    # constrcut the co-localization matrix
    cell_type_names <- colnames(weights)
    coc = matrix(NA, nrow = length(cell_type_names), ncol = length(cell_type_names))
    colnames(coc) = cell_type_names
    rownames(coc) = cell_type_names
    
    #ss1 = apply(bg, 2, sum)
    #ss2 = apply(border, 2, sum)
    nb.spots.border = nrow(border)
    nb.spots.bg = nrow(bg)
    expect_border = apply(border, 2, sum)
    expect_bg = apply(bg, 2, sum)
    
    for(m in 1:ncol(coc))
    {
      # m = 1
      cat(colnames(coc)[m], '\n')
      
      # first compute the border
      jj = which(border[,m] == TRUE)
      if(length(jj) == 1){
        cond_border = border[jj, ]
      }else{
        cond_border = apply(border[jj,], 2, sum)
      }
      
      # compute bg
      jj0 = which(bg[,m] == TRUE)
      if(length(jj0) == 1){
        cond_bg = bg[jj0, ]
      }else{
        cond_bg = apply(bg[jj0, ], 2, sum)
      }
      
      #score_border = cond_border/(expect_border)
      #score_bg = cond_bg /expect_bg
      #kk_score = which(!is.na(score_border) & !is.na(score_bg))
      #coc[, m] = cond_border/(expect_border + 0.5) * expect_bg/(cond_bg + 0.5)
      #coc[kk_score, m] = score_border[kk_score]/(score_bg[kk_score] + 0.1)
      
      for(k in 1:nrow(coc)){
        # k = 2
        dat <- data.frame(
          "withHits" = c(cond_border[k], cond_bg[k]),
          "noHits" =  c( expect_border[k] - cond_border[k], expect_bg[k] - cond_bg[k]),
          row.names = c("forground", "background"),
          stringsAsFactors = FALSE
        )
        
        coc[k, m] = fisher.test(x = dat, 
                                alternative = "greater")$p.value
        
      }
      
      
    }
    
    #ss1 = apply(coc, 1, sum)
    #ss2 = apply(coc, 2, sum)
    #sels = which(ss1>0 & ss2>0)
    #coc = coc[sels, sels]
    
    coc = -log10(coc)
    
    pdfname = paste0(outDir, "/cooccurence_score_BZ_subtypes_", slice, ".pdf")
    pdf(pdfname, width=12, height = 10)
    
    corrplot(coc, method = 'color', is.corr = FALSE, 
             #col.lim = c(-2.5, 2.5),
             hclust.method = c("complete"),
             col = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(6))
    
    dev.off()
    
  }
  
}


########################################################
########################################################
# test colocalization analysis wiht mistyR 
# original code from : 
# https://github.com/lengfei5/visium_heart/blob/master/st_snRNAseq/05_colocalization/colocalize_misty.R
########################################################
########################################################
run_misty_colocalization_analysis = function(st, 
                                             outDir,
                                             RCTD_out,
                                             condSpec_celltypes = NULL,
                                             segmentation_annots = c('all', 'BZ', 'RZ', 'Intact')
                                             )
{
  library(tidyverse)
  library(Seurat)
  library(mistyR)
  library(future)
  require(spacexr)
  # data manipulation
  library(Matrix)
  library(tibble)
  library(dplyr)
  library(purrr)
  
  # normalization
  library(sctransform)
  
  # resource
  #library(progeny)
  
  # setup parallel execution
  options(future.globals.maxSize = 80000 * 1024^2)
  # plan(multisession)
  source("utility_misty.R")
  #future::plan(future::multisession)
  
  # outDir = paste0(resDir, '/colocalization_misty')
  system(paste0('mkdir -p ', outDir))
  
  # Pipeline definition:
  run_colocalization <- function(slide, 
                                 assay, 
                                 useful_features, 
                                 slide_id, 
                                 cv.folds = 10,
                                 misty_out_alias = paste0(outDir, "/cm_")) 
  {
    # Define assay of each view ---------------
    view_assays <- list("main" = assay,
                        "juxta" = assay,
                        "para" = assay)
    
    # Define features of each view ------------
    view_features <- list("main" = useful_features, 
                          "juxta" = useful_features,
                          "para" = useful_features)
    
    # Define spatial context of each view -----
    view_types <- list("main" = "intra", 
                       "juxta" = "juxta",
                       "para" = "para")
    
    # Define additional parameters (l in case of paraview,
    # n of neighbors in case of juxta) --------
    view_params <- list("main" = NULL, 
                        "juxta" = 5,
                        "para" = 15)
    
    misty_out <- paste0(misty_out_alias, 
                        slide_id, "_", assay)
    
    run_misty_seurat(visium.slide = slide,
                     view.assays = view_assays,
                     view.features = view_features,
                     view.types = view_types,
                     view.params = view_params,
                     slide_id = slide_id,
                     cv.folds = cv.folds,
                     spot.ids = NULL,
                     out.alias = misty_out)
    
    return(misty_out)
    
  }
  
  # Main ------------------------------------------------------------------------
  cat('-- check visium conditions -- \n')
  st$condition = droplevels(as.factor(st$condition))
  print(table(st$condition))
  cc = names(table(st$condition))
  cat('ST conditions : \n')
  print(cc)
  tt = as.character(sapply(cc, function(x){unlist(strsplit(as.character(x), '_'))[2]}))
  cat('ST time :\n')
  print(tt)
  
  Idents(st) = st$condition
  
  for(n in 1:length(cc))
  {
    # n = 1
    cat(n, '  slice -- ', cc[n], '\n')
    slice = cc[n]
    stx = subset(st, condition == slice)
    DefaultAssay(stx) = 'Spatial'
    
    # import the cell type deconvolution result
    file_RCTD = paste0(RCTD_out,  '/RCTD_out_doubletMode_', slice, '.rds')
    if(!file.exists(file_RCTD)){
      file_RCTD = paste0(RCTD_out, '/', slice, '/RCTD_out_doubletMode_', slice, '.rds')
      if(!file.exists(file_RCTD)){
        cat(file_RCTD, ' not found \n');
      }
    }
    
    myRCTD = readRDS(file = file_RCTD)
    results <- myRCTD@results
    
    weights = results$weights
    
    #colnames(weights) = gsub('CM_ven_Robo2', 'CM_Robo2', colnames(weights))
    #colnames(weights) = gsub('CM_ven_Cav3_1', 'CM_Cav3.1', colnames(weights))
    
    norm_weights = normalize_weights(weights)
    weights = as.matrix(t(weights))
    norm_weights = as.matrix(t(norm_weights))
    
    if(!is.null(condSpec_celltypes)){
      cat('-- condition-specific cell types were specified --\n')
      
      celltypes_sels = condSpec_celltypes[[which(names(condSpec_celltypes) == tt[n])]]
      mm = match(celltypes_sels, rownames(weights))
      #celltypes_sels[which(is.na(mm))]
      
      if(length(which(is.na(mm))) == 0){
        weights = weights[mm, ]
        norm_weights = norm_weights[mm, ]
        cat('--', nrow(weights), ' celltypes to use -- \n')
      }else{
        stop(paste0(celltypes_sels[which(is.na(mm))], collapse = ' '), ' not found \n')
      }
    }else{
      cat('-- NO condition-specific cell types were specified --\n')
      
      # celltypes_sels = unique() 
      # mm = match(celltypes_sels, rownames(weights))
      # #celltypes_sels[which(is.na(mm))]
      # 
      # if(length(which(is.na(mm))) == 0){
      #   weights = weights[mm, ]
      #   norm_weights = norm_weights[mm, ]
      #   cat('--', nrow(weights), ' celltypes to use -- \n')
      # }else{
      #   stop(paste0(celltypes_sels[which(is.na(mm))], collapse = ' '), ' not found \n')
      # }
      
    }
    
    
    # filter the subtypes without matching
    mm = match(colnames(stx), colnames(norm_weights))
    stx = subset(stx, cells = colnames(stx)[which(!is.na(mm))])
    norm_weights = norm_weights[, mm[which(!is.na(mm))]]
    weights = weights[ ,mm[which(!is.na(mm))]]
    
    stx[["RCTD_propos"]] = CreateAssayObject(data = norm_weights)
    stx[["RCTD_density"]] = CreateAssayObject(data = weights)
    
    for(assay_label in c("RCTD_propos", "RCTD_density"))
    {
      for(segment in segmentation_annots)
      {
        # assay_label = 'RCTD_density'; segment = 'all'
        cat(' --- selected assay ', assay_label, '; selected region ', segment, '---\n')
        
        slide_id = slice
        assay <- assay_label
        DefaultAssay(stx) <- assay
        
        if(segment == 'all'){
          stx_sel = stx
          cv.folds = 10
          
        }else{
          #jj = grep(segment, stx$segmentation)
          jj = which(stx$segmentation == segment)
          if(length(jj)>0){
            stx_sel = subset(stx, cells = colnames(stx)[jj])
            cv.folds = 5
            
          }else{
            cat('no spots annoated in segment: ', segment, '\n')
            stx_sel = NULL
          }
          
        }
        
        if(!is.null(stx_sel)){
          xx = GetAssayData(stx_sel, slot = 'data', assay = assay)
          sds = apply(xx, 1, sd)
          useful_features <- rownames(stx_sel)[which(sds>0.001)]
          
          mout <- run_colocalization(slide = stx_sel,
                                     useful_features = useful_features,
                                     slide_id = slide_id,
                                     assay = assay,
                                     cv.folds = cv.folds,
                                     misty_out_alias = paste0(outDir, 'out_misty/'))
          
          misty_res_slide <- collect_results(mout)
          
          plot_folder <- paste0(outDir, "Plots_", assay_label)
          system(paste0("mkdir -p ", plot_folder))
          
          pdf(file = paste0(plot_folder, "/", slide_id, "_", segment,  "_summary_plots.pdf"), 
              height = 10, width = 16)
          
          mistyR::plot_improvement_stats(misty_res_slide)
          mistyR::plot_view_contributions(misty_res_slide)
          
          mistyR::plot_interaction_heatmap(misty_res_slide, "intra", cutoff = 0)
          mistyR::plot_interaction_communities(misty_res_slide, "intra", cutoff = 0.5)
          
          mistyR::plot_interaction_heatmap(misty_res_slide, "juxta_5", cutoff = 0)
          mistyR::plot_interaction_communities(misty_res_slide, "juxta_5", cutoff = 0.5)
          
          mistyR::plot_interaction_heatmap(misty_res_slide, "para_15", cutoff = 0)
          mistyR::plot_interaction_communities(misty_res_slide, "para_15", cutoff = 0.5)
          
          dev.off()
          
          transform_tibble_matrix = function(imp)
          {
            # imp = intra_imp 
            predictors = unique(imp$Predictor)
            targets = unique(imp$Target)
            
            res = matrix(NA, ncol = length(predictors), nrow = length(targets))
            colnames(res) = predictors[order(predictors, decreasing = FALSE)]
            rownames(res) = targets[order(targets, decreasing = TRUE)]
            
            for(nn in 1:nrow(imp)){
              # nn = 1
              index_col = which(colnames(res) == imp$Predictor[nn])
              index_row = which(rownames(res) == imp$Target[nn])
              res[index_row, index_col] = imp$Importance[nn]
              
            }
            return(res)
            
          }
          
          importances = as.data.frame(misty_res_slide$importances)
          importances = importances[, -1]
          importances$Importance[which(importances$Importance<0)] = 0
          
          
          intra_imp = importances[which(importances$view == 'intra'), ]
          intra_imp = transform_tibble_matrix(intra_imp)
          
          write.csv2(intra_imp, file = paste0(plot_folder, "/", slide_id, "_", segment,  
                                              "_summary_table_intra.csv"), quote = FALSE, 
                     row.names = TRUE)
          
          juxta_imp = importances[which(importances$view == 'juxta_5'), ]
          juxta_imp = transform_tibble_matrix(juxta_imp)
          
          write.csv2(juxta_imp, file = paste0(plot_folder, "/", slide_id, "_", segment,  
                                              "_summary_table_juxta5.csv"), quote = FALSE, 
                     row.names = TRUE)
          
          para_imp = importances[which(importances$view == 'para_15'), ]
          para_imp = transform_tibble_matrix(para_imp)
          
          write.csv2(para_imp, file = paste0(plot_folder, "/", slide_id, "_", segment,  
                                             "_summary_table_para15.csv"), quote = FALSE, 
                     row.names = TRUE)
          
        }
       
      }
     
    }
  
  }
}


##########################################
# compare the bz vs remote regions and also the control samples at day 0
##########################################
get_comparable_matrix = function(xx, yy)
{
  # xx = bz; yy = ctl1;
  res = matrix(NA, nrow = nrow(xx), ncol = ncol(xx))
  
  jj1 = match(colnames(xx), colnames(yy))
  jj2 = match(rownames(xx), rownames(yy))
  
  if(length(which(is.na(jj1))) == 0 & length(which(is.na(jj2))) == 0){
    res = as.matrix(yy[jj2, jj1])
    
  }else{
    for(m in 1:length(jj1))
    {
      if(!is.na(jj1[m])){
        res[,m] = yy[jj2, jj1[m]]
      }
    }
  }
  
  colnames(res) = colnames(xx)
  rownames(res) = rownames(xx)
  
  
  res[which(is.na(res))] = 0
  return(res)
  
}

##summarize_colocalization with parameters bz_rz_intra, bz_rz_juxta, bz_rz_para,
##bz_c_intra, bz_c_juxta, bz_c_para
summarize_misty_mode = function(x1, x2, x3, method = c('median', 'max', 'mean', 'sum'))
{
  col_names = colnames(x1)
  row_names = rownames(x1)
  
  res = matrix(0, ncol = ncol(x1), nrow = nrow(x1))
  x1 = as.matrix(x1)
  x1[which(is.na(x1))] = 0
  x2 = as.matrix(x2)
  x2[which(is.na(x2))] = 0
  x3 = as.matrix(x3)
  x3[which(is.na(x3))] = 0
  
  for(i in 1:ncol(x1))
  {
    for(j in 1:nrow(x1))
    {
      if(col_names[i] != row_names[j]){
        #if(x1[j, i] > 0 | x2[j, i] > 0 | x3[j, i] > 0){
        if(method == 'max'){res[j, i] = max(c(x1[j, i], x2[j, i], x3[j, i]));}  
        if(method == 'sum'){res[j, i] = sum(c(x1[j, i], x2[j, i], x3[j, i]));}
        if(method == 'median') {res[j, i] = median(c(x1[j, i], x2[j, i], x3[j, i]));}
        if(method == 'mean') {res[j, i] = mean(c(x1[j, i], x2[j, i], x3[j, i]));}
        
      }
    }
  }
  
  colnames(res) = col_names
  rownames(res) = row_names
  
  return(res)
  
}

convert_tibble_for_save = function(sample = 'Amex_d0_294946', region = 'Intact', xx)
{
  xx = data.frame(sample = rep(sample, nrow(xx)), 
                  region = region,
                  source = rownames(xx), 
                  xx, stringsAsFactors = FALSE)
  
  xx = pivot_longer(xx, cols = colnames(xx)[-c(1:3)],
                    names_to = "target", 
                    values_to = "importance")
  
  return(xx)
}

transform_zscore = function(x)
{
  # x = ctl_intra
  x = (x + t(x))/2.0
  
  xx = as.numeric(as.matrix(x))
  xx = (xx - mean(xx, na.rm = TRUE))/sd(xx, na.rm = TRUE)
  xx = matrix(xx, nrow = nrow(x), ncol = ncol(x), byrow = FALSE)
  colnames(xx) = colnames(x)
  rownames(xx) = rownames(x)
  
  return(xx)
  
}


##########################################
# summarize the misty analysis and plot the network with graph
##########################################
summarize_cell_neighborhood_misty = function(st, 
                                      outDir,
                                      time = c('d4', 'd7'),
                                      misty_mode = c('density'),
                                      summary_method = 'median',
                                      cutoff = 0.2,
                                      resolution = 0.7,
                                      segmentation_annots = c('all', 'BZ', 'RZ', 'Intact'), 
                                      controls = c('RZ', 'Intact')
                                      )
{
  # time = c('d1', 'd4', 'd7', 'd14'); segmentation_annots = c('all', 'BZ', 'RZ', 'Intact'); 
  # misty_mode = c('density'); summary_method = 'median';
  cat(' -- significance test of celltype proximity for Misty output -- \n')
  cat(paste0('--- misty mode : ', misty_mode, '\n'))
  library("plyr")
  library("reshape2")
  library("ggplot2")
  library(igraph)
  library(ggraph)
  library(CellChat)
  library(tidyr)
  
  testDir = paste0(outDir, '/signficant_neighborhood_', misty_mode, '/')
  if(!dir.exists(testDir)){
    system(paste0('mkdir -p ', testDir))
  }
  
  cat('-- check visium conditions -- \n')
  st$condition = droplevels(as.factor(st$condition))
  print(table(st$condition))
  cc = names(table(st$condition))
  cat('ST conditions : \n')
  print(cc)
  
  tt = as.character(sapply(cc, function(x){unlist(strsplit(as.character(x), '_'))[2]}))
  cat('ST time :\n')
  print(tt)
  
  Idents(st) = st$condition
  
  ## process first the intact samples
  ctl_intra = read.csv2(file = paste0(outDir, 'Plots_RCTD_', misty_mode, '/',
                                      'Amex_d0_294946_Intact_summary_table_intra', '.csv'), row.names = c(1))
  ctl_intra = ctl_intra[c(nrow(ctl_intra):1), ]
  ctl_intra = transform_zscore(ctl_intra)
  
  ctl_juxta = read.csv2(file = paste0(outDir, 'Plots_RCTD_', misty_mode, '/',
                                      'Amex_d0_294946_Intact_summary_table_juxta5.csv'), row.names = c(1))
  ctl_juxta = get_comparable_matrix(ctl_intra, ctl_juxta)
  ctl_juxta = transform_zscore(ctl_juxta)
  
  ctl_para = read.csv2(file = paste0(outDir, 'Plots_RCTD_', misty_mode, '/',
                                     'Amex_d0_294946_Intact_summary_table_para15.csv'), row.names = c(1))
  ctl_para = get_comparable_matrix(ctl_intra, ctl_para)
  ctl_para = transform_zscore(ctl_para)
  
  ctl1 = summarize_misty_mode(ctl_intra, ctl_juxta, ctl_para, method = summary_method)
  rm(list = c('ctl_intra', 'ctl_juxta', 'ctl_para'))
  
  prox_ctl = convert_tibble_for_save(sample = 'd0_294946', region = 'Intact', ctl1)
  
  write.table(ctl1, 
              file = paste0(testDir, 'cell_cell_colocalization_summary_', 'Amex_d0_294946', '.txt'), 
              quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
  
  ## plot the heatmap of cell-cell colocalization
  pdfname = paste0(testDir, 'Plot_', 'Amex_d0_294946', '.pdf')
  pdf(pdfname, width=14, height = 10)
  
  for(cutoff in seq(0, max(ctl1, na.rm = TRUE)*0.6, by=0.1))
  {
    # cutoff = 0.
    set2.blue <- "#8DA0CB"
    
    bz_allx = ctl1
    bz_allx[bz_allx < cutoff] <- 0
    plot.data = reshape2::melt(bz_allx, value.name = 'Importance')
    colnames(plot.data)[1:2] = c('Target', 'Predictor')
    
    results.plot <- ggplot2::ggplot(
      plot.data,
      ggplot2::aes(
        x = Predictor,
        y = reorder(Target, desc(Target))
      )
    ) +
      ggplot2::geom_tile(ggplot2::aes(fill = Importance)) + 
      ggplot2::scale_fill_gradient2(
        low = "white",
        mid = "white",
        high = set2.blue,
        midpoint = 0.2
      ) +
      ggplot2::theme_classic() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::labs(y = 'Target') +
      ggplot2::coord_equal() +
      ggplot2::ggtitle(paste0('colocalization heatmap with cutoff : ', cutoff))
    
    print(results.plot)
    
    . <- NULL
    A = ctl1
    A[A < cutoff | is.na(A)] <- 0
    #A = A[c(nrow(A):1), ]
    
    #A[is.na(A)] <- 0
    G <- igraph::graph.adjacency(A, mode = "plus", weighted = TRUE) %>%
      igraph::set.vertex.attribute("name", value = names(igraph::V(.))) %>%
      igraph::delete.vertices(which(igraph::degree(.) == 0))
    
    ## try make similar plot but the cellchat circle style
    ## original code from cellchat function netVisual_circle
    cols_C = scPalette(nrow(A))
    names(cols_C) = colnames(A)
    my_netVisual_circle(A,
                        color.use = cols_C,
                        #vertex.weight = groupSize, 
                        weight.scale = TRUE, 
                        label.edge= FALSE, 
                        title.name = paste0("Cell-Cell colocalization with cutoff : ", cutoff))
    
    
  }
  
  dev.off()
  
  
  
  ctl2 = read.csv2(file = paste0(outDir, 'Plots_RCTD_', misty_mode, '/',
                                 'Amex_d0_294949_Intact_summary_table_intra', '.csv'), 
                   row.names = c(1))
  
  ctl_intra = read.csv2(file = paste0(outDir, 'Plots_RCTD_', misty_mode, '/',
                                      'Amex_d0_294949_Intact_summary_table_intra', '.csv'), 
                        row.names = c(1))
  ctl_intra = ctl_intra[c(nrow(ctl_intra):1), ]
  ctl_intra = transform_zscore(ctl_intra)
  
  ctl_juxta = read.csv2(file = paste0(outDir, 'Plots_RCTD_', misty_mode, '/',
                                      'Amex_d0_294949_Intact_summary_table_juxta5.csv'), 
                        row.names = c(1))
  ctl_juxta = get_comparable_matrix(ctl_intra, ctl_juxta)
  ctl_juxta = transform_zscore(ctl_juxta)
  
  ctl_para = read.csv2(file = paste0(outDir, 'Plots_RCTD_', misty_mode, '/',
                                     'Amex_d0_294949_Intact_summary_table_para15.csv'), 
                       row.names = c(1))
  ctl_para = get_comparable_matrix(ctl_intra, ctl_para)
  ctl_para = transform_zscore(ctl_para)
  
  ctl2 = summarize_misty_mode(ctl_intra, ctl_juxta, ctl_para, method = summary_method)
  rm(list = c('ctl_intra', 'ctl_juxta', 'ctl_para'))
  
  prox_ctl = bind_rows(prox_ctl, 
                       convert_tibble_for_save(sample = 'd0_294949', region = 'Intact', ctl2))
  
  
  write.table(ctl2, 
              file = paste0(testDir, 'cell_cell_colocalization_summary_', 'Amex_d0_294949', '.txt'), 
              quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
  
  ## plot the heatmap of cell-cell colocalization
  pdfname = paste0(testDir, 'Plot_', 'Amex_d0_294949', '.pdf')
  pdf(pdfname, width=14, height = 10)
  
  for(cutoff in seq(0, max(ctl1, na.rm = TRUE)*0.6, by=0.1))
  {
    # cutoff = 0.
    set2.blue <- "#8DA0CB"
    
    bz_allx = ctl2
    bz_allx[bz_allx < cutoff] <- 0
    plot.data = reshape2::melt(bz_allx, value.name = 'Importance')
    colnames(plot.data)[1:2] = c('Target', 'Predictor')
    
    results.plot <- ggplot2::ggplot(
      plot.data,
      ggplot2::aes(
        x = Predictor,
        y = reorder(Target, desc(Target))
      )
    ) +
      ggplot2::geom_tile(ggplot2::aes(fill = Importance)) + 
      ggplot2::scale_fill_gradient2(
        low = "white",
        mid = "white",
        high = set2.blue,
        midpoint = 0.2
      ) +
      ggplot2::theme_classic() + 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::labs(y = 'Target') +
      ggplot2::coord_equal() +
      ggplot2::ggtitle(paste0('colocalization heatmap with cutoff : ', cutoff))
    
    print(results.plot)
    
    . <- NULL
    A = ctl2
    A[A < cutoff | is.na(A)] <- 0
    #A = A[c(nrow(A):1), ]
    
    #A[is.na(A)] <- 0
    G <- igraph::graph.adjacency(A, mode = "plus", weighted = TRUE) %>%
      igraph::set.vertex.attribute("name", value = names(igraph::V(.))) %>%
      igraph::delete.vertices(which(igraph::degree(.) == 0))
    
    ## try make similar plot but the cellchat circle style
    ## original code from cellchat function netVisual_circle
    cols_C = scPalette(nrow(A))
    names(cols_C) = colnames(A)
    my_netVisual_circle(A,
                        color.use = cols_C,
                        #vertex.weight = groupSize, 
                        weight.scale = TRUE, 
                        label.edge= FALSE, 
                        title.name = paste0("Cell-Cell colocalization with cutoff : ", cutoff))
    
    
  }
  
  dev.off()
  
  
  ## process the time points
  for(n in 1:length(time))
  {
    # n = 2
    kk = which(tt == time[n])
    
    for(k in kk)
    {
      # k = kk[1]
      cat('---- condition :', cc[k], '---- \n')
      
      intra = read.csv2(file = paste0(outDir, 'Plots_RCTD_', misty_mode, '/',
                                   cc[k], '_BZ_summary_table_intra', '.csv'), row.names = c(1))
      intra = intra[c(nrow(intra):1), ]
      intra = transform_zscore(intra)
      
      juxta = read.csv2(file = paste0(outDir, 'Plots_RCTD_',  misty_mode, '/',
                                   cc[k], '_BZ_summary_table_juxta5.csv'), row.names = c(1))
      juxta = get_comparable_matrix(intra, juxta)
      juxta = transform_zscore(juxta)
      
      para = read.csv2(file = paste0(outDir, 'Plots_RCTD_', misty_mode, '/',
                                   cc[k], '_BZ_summary_table_para15.csv'), row.names = c(1))
      para = get_comparable_matrix(intra, para)
      para = transform_zscore(para)
      
      bz = summarize_misty_mode(intra, juxta, para, method = summary_method)
      
      rm(list = c('intra', 'juxta', 'para'))
      
      prox_all = bind_rows(prox_ctl, 
                           convert_tibble_for_save(sample = cc[k], 
                                                   region = 'BZ', 
                                                   bz))
      
      intra = read.csv2(file = paste0(outDir, 'Plots_RCTD_', misty_mode, '/',
                                      cc[k], '_RZ_summary_table_intra', '.csv'), row.names = c(1))
      intra = intra[c(nrow(intra):1), ]
      intra = transform_zscore(intra)
      
      juxta = read.csv2(file = paste0(outDir, 'Plots_RCTD_',  misty_mode, '/',
                                      cc[k], '_RZ_summary_table_juxta5.csv'), row.names = c(1))
      juxta = get_comparable_matrix(intra, juxta)
      juxta = transform_zscore(juxta)
      
      para = read.csv2(file = paste0(outDir, 'Plots_RCTD_', misty_mode, '/',
                                     cc[k], '_RZ_summary_table_para15.csv'), row.names = c(1))
      para = get_comparable_matrix(intra, para)
      para = transform_zscore(para)
      
      rz = summarize_misty_mode(intra, juxta, para, method = summary_method)
      
      rm(list = c('intra', 'juxta', 'para'))
      
      prox_all = bind_rows(prox_all, 
                           convert_tibble_for_save(sample = cc[k], 
                                                   region = 'RZ', 
                                                   rz))
      rz = get_comparable_matrix(bz, rz)
      
      bz_all = bz - rz
      #bz_all[bz_all<0] = 0
      
      write.table(bz_all, 
                  file = paste0(testDir, 'cell_cell_colocalization_summary_', cc[k], '.txt'), 
                  quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
      
      # bz_all = pairs
      ## plot the heatmap of cell-cell colocalization
      pdfname = paste0(testDir, 'Plot_', cc[k], '.pdf')
      pdf(pdfname, width=14, height = 10)
      
      for(cutoff in seq(0, max(bz_all, na.rm = TRUE)*0.6, by=0.1))
      {
        # cutoff = 0.
        set2.blue <- "#8DA0CB"
        
        bz_allx = bz_all
        bz_allx[bz_allx < cutoff] <- 0
        plot.data = reshape2::melt(bz_allx, value.name = 'Importance')
        colnames(plot.data)[1:2] = c('Target', 'Predictor')
        
        results.plot <- ggplot2::ggplot(
          plot.data,
          ggplot2::aes(
            x = Predictor,
            y = reorder(Target, desc(Target))
          )
        ) +
          ggplot2::geom_tile(ggplot2::aes(fill = Importance)) + 
          ggplot2::scale_fill_gradient2(
            low = "white",
            mid = "white",
            high = set2.blue,
            midpoint = 0.2
          ) +
          ggplot2::theme_classic() + 
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
          ggplot2::labs(y = 'Target') +
          ggplot2::coord_equal() +
          ggplot2::ggtitle(paste0('colocalization heatmap with cutoff : ', cutoff))
        
        print(results.plot)
        
        . <- NULL
        A = bz_all
        A[A < cutoff | is.na(A)] <- 0
        #A = A[c(nrow(A):1), ]
        
        #A[is.na(A)] <- 0
        G <- igraph::graph.adjacency(A, mode = "plus", weighted = TRUE) %>%
          igraph::set.vertex.attribute("name", value = names(igraph::V(.))) %>%
          igraph::delete.vertices(which(igraph::degree(.) == 0))
        
        ## try make similar plot but the cellchat circle style
        ## original code from cellchat function netVisual_circle
        cols_C = scPalette(nrow(A))
        names(cols_C) = colnames(A)
        my_netVisual_circle(A,
                            color.use = cols_C,
                            #vertex.weight = groupSize, 
                            weight.scale = TRUE, 
                            label.edge= FALSE, 
                            title.name = paste0("Cell-Cell colocalization with cutoff : ", cutoff))
        
        
      }
      
      dev.off()
      
    }
    
    ## barpot to compare the replicates, BZ, RZ and also Intact
    prox_all = prox_all %>% mutate(segment = paste0(region, '_', sample)) 
    
    pdfname = paste0(testDir, 'Plot_neighborhood_rep_regions_', time[n], '.pdf')
    pdf(pdfname, width=16, height = 8)
    
    if(length(kk) == 2){
      
      bz1 = read.table(file = paste0(testDir, 'cell_cell_colocalization_summary_', cc[kk[1]], '.txt'), 
                  row.names = c(1),  header = TRUE, sep = '\t')
      
      bz2 = read.table(file = paste0(testDir, 'cell_cell_colocalization_summary_', cc[kk[2]], '.txt'), 
                       row.names = c(1),  header = TRUE, sep = '\t')
      
      bz2 = get_comparable_matrix(bz1, bz2)
      
      bz1 = as.numeric(as.matrix(bz1))
      #bz1[which(bz1<0)] = 0
      bz2 = as.numeric(as.matrix(bz2))
      #bz2[which(bz2<0)] = 0
      plot(bz1, bz2, xlab = cc[kk[1]], ylab = cc[kk[2]]);
      abline(0, 1, col = 'red', lwd =1.0)
      cat(time[n],  'correlation between two reps -- ', cor(bz1, bz2), '\n')
      
    }
    
    
    for(s in unique(prox_all$source))
    {
      # s = "CM.Prol.IS"
      #cat(s, '\n')
      p1 = prox_all %>% filter(source == s) %>% #filter(target == "T.cells")
        ggplot(aes(x=target, y=importance, fill=segment)) +
        geom_bar(stat="identity", position= position_dodge(), width = 0.7) +
        scale_x_discrete(drop = FALSE) + 
        scale_fill_brewer(palette="Paired", direction = -1) +
        theme_minimal() + 
        #ylim(-1, 2) +
        theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 1)) + 
        ggtitle(paste0("cell types colocalizing with ", s))
        
      plot(p1)
    }
    
    
    dev.off()
    
              
    
    
      
  }
  
}

my_netVisual_circle <-function(net, color.use = NULL,title.name = NULL, 
                               sources.use = NULL, targets.use = NULL, idents.use = NULL, 
                               remove.isolate = FALSE, top = 1,
                               weight.scale = FALSE, 
                               vertex.weight = 20, vertex.weight.max = NULL, 
                               vertex.size.max = NULL, vertex.label.cex=1,vertex.label.color= "black",
                               edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6, 
                               label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                               edge.curved=0.2, shape='circle', layout=in_circle(), 
                               margin=0.2, vertex.size = NULL,
                               arrow.width=1,arrow.size = 0.2,
                               text.x = 0, text.y = 1.5)
{
  # net = A; color.use = cols_C; weight.scale = TRUE; label.edge= FALSE
  if (!is.null(vertex.size)) {
    warning("'vertex.size' is deprecated. Use `vertex.weight`")
  }
  if (is.null(vertex.size.max)) {
    if (length(unique(vertex.weight)) == 1) {
      vertex.size.max <- 5
    } else {
      vertex.size.max <- 15
    }
  }
  options(warn = -1)
  thresh <- stats::quantile(net, probs = 1-top)
  net[net < thresh] <- 0
  
  if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use)) ) {
    if (is.null(rownames(net))) {
      stop("The input weighted matrix should have rownames!")
    }
    cells.level <- rownames(net)
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- cells.level[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- cells.level[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    if (!is.null(idents.use)) {
      if (is.numeric(idents.use)) {
        idents.use <- cells.level[idents.use]
      }
      df.net <- filter(df.net, (source %in% idents.use) | (target %in% idents.use))
    }
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  
  if (is.null(color.use)) {
    color.use = scPalette(nrow(net))
    names(color.use) <- rownames(net)
  } else {
    if (is.null(names(color.use))) {
      stop("The input `color.use` should be a named vector! \n")
    }
  }
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx.isolate <- intersect(idx1, idx2)
    if (length(idx.isolate) > 0) {
      net <- net[-idx.isolate, ]
      net <- net[, -idx.isolate]
      color.use = color.use[-idx.isolate]
      if (length(unique(vertex.weight)) > 1) {
        vertex.weight <- vertex.weight[-idx.isolate]
      }
    }
  }
  
  g <- graph_from_adjacency_matrix(net, mode = "lower", weighted = TRUE)
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5
  
  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }
  
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  igraph::E(g)$loop.angle <- rep(0, length(igraph::E(g)))
  
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
  if (!is.null(title.name)) {
    text(text.x,text.y,title.name, cex = 1.1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}


## this is a network plot function from Giotto
## (https://github.com/RubD/Giotto/blob/master/R/spatial_interaction_visuals.R)
cellProximityNetwork = function(gobject,
                                CPscore,
                                remove_self_edges = FALSE,
                                self_loop_strength = 0.1,
                                color_depletion = 'lightgreen',
                                color_enrichment = 'red',
                                rescale_edge_weights = TRUE,
                                edge_weight_range_depletion = c(0.1, 1),
                                edge_weight_range_enrichment = c(1, 5),
                                layout = c('Fruchterman', 'DrL', 'Kamada-Kawai'),
                                only_show_enrichment_edges = F,
                                edge_width_range = c(0.1, 2),
                                node_size = 4,
                                node_text_size = 6,
                                show_plot = NA,
                                return_plot = NA,
                                save_plot = NA,
                                save_param =  list(),
                                default_save_name = 'cellProximityNetwork') {
  
  # extract scores
  
  # data.table variables
  cell_1 = cell_2 = unified_int = color = size = name = NULL
  
  CPscores = CPscore[['enrichm_res']]
  CPscores[, cell_1 := strsplit(as.character(unified_int), split = '--')[[1]][1], by = 1:nrow(CPscores)]
  CPscores[, cell_2 := strsplit(as.character(unified_int), split = '--')[[1]][2], by = 1:nrow(CPscores)]
  
  # create igraph with enrichm as weight edges
  igd = igraph::graph_from_data_frame(d = CPscores[,c('cell_1', 'cell_2', 'enrichm')], directed = F)
  
  if(remove_self_edges == TRUE) {
    igd = igraph::simplify(graph = igd, remove.loops = TRUE, remove.multiple = FALSE)
  }
  
  edges_sizes = igraph::get.edge.attribute(igd, 'enrichm')
  post_edges_sizes = edges_sizes[edges_sizes > 0]
  neg_edges_sizes = edges_sizes[edges_sizes <= 0]
  
  # rescale if wanted
  if(rescale_edge_weights == TRUE) {
    pos_edges_sizes_resc = scales::rescale(x = post_edges_sizes, to = edge_weight_range_enrichment)
    neg_edges_sizes_resc = scales::rescale(x = neg_edges_sizes, to = edge_weight_range_depletion)
    edges_sizes_resc = c(pos_edges_sizes_resc, neg_edges_sizes_resc)
  } else {
    edges_sizes_resc = c(post_edges_sizes, neg_edges_sizes)
  }
  
  # colors
  edges_colors = ifelse(edges_sizes > 0, 'enriched', 'depleted')
  
  
  # create coordinates for layout
  if(class(layout) %in% c('data.frame', 'data.table')) {
    if(ncol(layout) < 2) {
      stop('custom layout needs to have at least 2 columns')
    }
    
    if(nrow(layout) != length(igraph::E(igd))) {
      stop('rows of custom layout need to be the same as number of edges')
    }
    
  } else {
    layout = match.arg(arg = layout, choices = c('Fruchterman', 'DrL', 'Kamada-Kawai'))
  }
  
  
  
  
  #iplot = igraph::plot.igraph(igd, edge.color = edges_colors, edge.width = edges_sizes_resc, layout = coords)
  
  igd = igraph::set.edge.attribute(graph = igd, index = igraph::E(igd), name = 'color', value = edges_colors)
  igd = igraph::set.edge.attribute(graph = igd, index = igraph::E(igd), name = 'size', value = as.numeric(edges_sizes_resc))
  
  ## only show attractive edges
  if(only_show_enrichment_edges == TRUE) {
    colors = igraph::get.edge.attribute(igd, name = 'color')
    subvertices_ids = which(colors == 'enriched')
    igd = igraph::subgraph.edges(graph = igd, eids = subvertices_ids)
    
    # get new rescale vector (in case vector id is lost)
    edges_sizes_resc = igraph::E(igd)$size
  }
  
  ## get coordinates layouts
  if(layout == 'Fruchterman') {
    coords = igraph::layout_with_fr(graph = igd, weights = edges_sizes_resc)
  } else if(layout == 'DrL') {
    coords = igraph::layout_with_drl(graph = igd, weights = edges_sizes_resc)
  } else if(layout == 'Kamada-Kawai') {
    coords = igraph::layout_with_kk(graph = igd, weights = edges_sizes_resc)
  } else {
    stop('\n Currently no other layouts have been implemented \n')
  }
  
  
  #longDT = as.data.table(igraph::as_long_data_frame(igd))
  #return(longDT)
  #return(list(igd, coords))
  
  ## create plot
  gpl = ggraph::ggraph(graph = igd, layout = coords)
  gpl = gpl + ggraph::geom_edge_link(ggplot2::aes(color = factor(color), edge_width = size, edge_alpha = size), show.legend = F)
  gpl = gpl + ggraph::geom_edge_loop(ggplot2::aes(color = factor(color), edge_width = size, edge_alpha = size, strength = self_loop_strength), show.legend = F)
  gpl = gpl + ggraph::scale_edge_color_manual(values = c('enriched' = color_enrichment, 'depleted' = color_depletion))
  gpl = gpl + ggraph::scale_edge_width(range = edge_width_range)
  gpl = gpl + ggraph::scale_edge_alpha(range = c(0.1,1))
  gpl = gpl + ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = node_text_size)
  gpl = gpl + ggraph::geom_node_point(size = node_size)
  gpl = gpl + ggplot2::theme_bw() + ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                                   panel.border = ggplot2::element_blank(),
                                                   axis.title = ggplot2::element_blank(),
                                                   axis.text = ggplot2::element_blank(),
                                                   axis.ticks = ggplot2::element_blank())
  
  
  # print, return and save parameters
  show_plot = ifelse(is.na(show_plot), readGiottoInstructions(gobject, param = 'show_plot'), show_plot)
  save_plot = ifelse(is.na(save_plot), readGiottoInstructions(gobject, param = 'save_plot'), save_plot)
  return_plot = ifelse(is.na(return_plot), readGiottoInstructions(gobject, param = 'return_plot'), return_plot)
  
  ## print plot
  if(show_plot == TRUE) {
    print(gpl)
  }
  
  ## save plot
  if(save_plot == TRUE) {
    do.call('all_plots_save_function', c(list(gobject = gobject, plot_object = gpl, default_save_name = default_save_name), save_param))
  }
  
  ## return plot
  if(return_plot == TRUE) {
    return(gpl)
  }
  
  
}



##########################################
# test Giotto for neighborhood analysis  
##########################################
analyze.celltype.proximity.network_Giotto = function()
{
  library(Giotto)
  installGiottoEnvironment()
  
  
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
  