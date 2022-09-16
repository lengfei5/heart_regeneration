##########################################################################
##########################################################################
# Project:
# Script purpose: analyze the RCTD output results
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Sep 16 11:03:42 2022
##########################################################################
##########################################################################
library(spacexr)
require(Matrix)

cat('-- check visium conditions -- \n')
print(table(st$condition))
cc = names(table(st$condition))

for(n in 1:length(cc))
#for(n in c(1, 2, 4))
{
  # n = 2
  cat('slice -- ', cc[n], '\n')
  slice = cc[n]
  #stx = st[, which(st$condition == slice)]
  stx = subset(st, condition == slice)
  stx@images = stx@images[[n]]
  
  resultsdir <- paste0(RCTD_out, '/', slice)
  system(paste0('mkdir -p ', resultsdir))
  #resultsdir <- paste0(RCTD_out, '/', slice, '_Plots') ## you may change this to a more accessible directory on your computer.
  #dir.create(resultsdir)
  
  ##########################################
  # check the result
  ##########################################
  myRCTD = readRDS(file = paste0(resultsdir, '/RCTD_out_doubletMode_', slice, '.rds'))
  results <- myRCTD@results
  
  # normalize the cell type proportions to sum to 1.
  #norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
  norm_weights = normalize_weights(results$weights) 
  cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
  
  spatialRNA <- myRCTD@spatialRNA
  
  # make the plots 
  if(plot.RCTD.summary){
    stx$condition = droplevels(stx$condition)
    DefaultAssay(stx) = 'SCT'
    
    SpatialFeaturePlot(stx, images =  slice, slot = 'data',
                       features = "LOC115076416-AMEX60DD051098", image.alpha = 0.5) +
      ggtitle('NPPA-AMEX60DD051098') +
      coord_flip() + 
      ggplot2::scale_y_reverse() # rotate the coordinates to have the same as RCTD 
    
    ggsave(filename =  paste0(resultsdir, "/Spatial_patterning_", slice, "_NPPA.pdf"), width = 12, height = 8)
      
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
  
}