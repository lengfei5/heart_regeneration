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
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered', species, '.Rdata'))
st$condition = factor(st$condition, levels = design$condition)

plot.cell2location.summary = FALSE
PLOT.scatterpie = TRUE


cell2location_file = '../results/visium_axolotl_R12830_resequenced_20220308/cell2location_coarse_out/predictions_cell2loc_v1.csv'
resultsdir = paste0(resDir, '/cell2location_coarse_out/summary')

#cell2location_file = '../results/visium_axolotl_R12830_resequenced_20220308/cell2location_subtype_out/predictions_cell2loc.csv'
#resultsdir = paste0(resDir, '/cell2location_subtype_out/summary')
system(paste0('mkdir -p ', resultsdir))

res = read.csv(cell2location_file)
mm = match(colnames(st), res$CellID)
length(mm); length(which(!is.na(mm)))

res = res[mm, ]
rownames(res) = res$CellID
res = res[, c(17:ncol(res))]
res = normalize_weights(res)

st = AddMetaData(st, metadata = res)

cat('-- start to check cell decomposition result condition by conditions -- \n')

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
  #stx@images = stx@images[[n]]
  
  ##########################################
  # check the result from cell2location
  ##########################################
  results = res[which(st$condition == cc[n]), ]
  
  # normalize the cell type proportions to sum to 1.
  norm_weights = normalize_weights(results) 
  cell_type_names <- colnames(norm_weights)
  
  # spatialRNA <- myRCTD@spatialRNA
  
  # make the plots 
  stx$condition = droplevels(stx$condition)
  DefaultAssay(stx) = 'SCT'
  
  SpatialFeaturePlot(stx, images =  slice, #slot = 'data',
                     features = "LOC115076416-AMEX60DD051098") +
    ggtitle('NPPA-AMEX60DD051098')
    #coord_flip() + 
    #ggplot2::scale_y_reverse() # rotate the coordinates to have the same as RCTD 
  
  ggsave(filename =  paste0(resultsdir, "/Spatial_patterning_", slice, "_NPPA.pdf"), width = 12, height = 8)
  
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
    # use a panel of colors from https://gotellilab.github.io/GotelliLabMeetingHacks/NickGotelli/ColorPalettes.html
    #tol10qualitative=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77",
    #                   "#661100", "#CC6677", "#882255", "#AA4499")
    # set color vector
    getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
    
    ## Preprocess coordinates 
    spatial_coord <-  GetTissueCoordinates(stx, image = slice) %>%
      tibble::rownames_to_column("ID")
    
    weights = norm_weights[match(spatial_coord$ID, rownames(norm_weights)), ]
    
    # clean the weight: force to be 0 if lower than certain threshold, e.g. 0.05 or 0.1
    cutoff_weight = 0.05
    weights_processed = as.matrix(weights)
    weights_processed[which(weights_processed < cutoff_weight)] = 0
    weights_processed = normalize_weights(weights_processed)
    
    weights = weights_processed;
    rm(weights_processed)
    weights = weights[, apply(weights, 2, sum)>0]
    
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
    #spatial_coord$x = as.numeric(spatial_coord$x)
    #spatial_coord$y = as.numeric(spatial_coord$y)
    
    ggplot() + geom_scatterpie(aes(x=imagerow, y=imagecol), data=spatial_coord,
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
    
    ggsave(paste0(resultsdir, '/cell2location_scatterpie_', slice, '_v0.pdf'), width = 16, height = 10)
    
    
    pdfname = paste0(resultsdir, '/Spatial_distribution_each_celltype_', slice, '.pdf')
    pdf(pdfname, width=16, height = 10)
    
    for(m in 1:length(cell_types_plt)){
      # m = 5
      sub.spatial = spatial_coord[,c(1:3, (3+m))]
      colnames(sub.spatial)[ncol(sub.spatial)] = 'celltype'
      sub.spatial$celltype = as.numeric(as.character(sub.spatial$celltype))
      sub.spatial = data.frame(sub.spatial)
      
      p1 = ggplot(data = sub.spatial, aes(x=imagerow, y=imagecol, size = celltype)) +
        geom_point(color = col_ct[m]) +  
        #scale_size_continuous(name = "celltype") +
        #scale_size(range = c(0, 0.8)) + 
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
        ggplot2::guides(fill = guide_legend(ncol = 1))
      
      plot(p1)
      
    }
    
    dev.off()
    
    
  }
  
}