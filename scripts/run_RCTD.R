##########################################
##########################################
# an example of running RCTD to deconvolve cell composition in Visium data 
# the main input data: visium data and reference (scRNA-seq/snRNA-seq data), both saved as Seurat object here
# the RCTD paper Cable et al. 2022 (https://www.nature.com/articles/s41587-021-00830-w)
# the RCTD vignette https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html 
##########################################
##########################################
rm(list = ls())

library(spacexr) 
require(Matrix)
require(RColorBrewer)
require(scatterpie)
require(cowplot)

## load the visium and reference data
st = readRDS(file = paste0('data_examples/st_visium.rds'))
refs = readRDS(file = paste0('data_examples/ref_scRNAseq.rds'))

# select the cell labels to use in the deconvolution
refs$celltype_toUse = refs$celltype

table(refs$celltype_toUse)
length(table(refs$celltype_toUse))

## preapre the paramters for RCTD subtypes
DefaultAssay(refs) = 'integrated'
DefaultAssay(st) = 'Spatial'
require_int_SpatialRNA = FALSE # important parameter, depending on the count data is integer or not in ref
max_cores = 16

resultsdir <- 'test'
system(paste0('mkdir -p ', resultsdir))

##########################################
# ## prepare the reference
##########################################
cat('-- prepare global reference --\n')

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

table(reference@cell_types) #number of occurences for each cell type

##########################################
# prepare ST data for RTCD
##########################################
counts = GetAssayData(object = st, slot = "counts")

# select the image coordinates
coords <- eval(parse(text = paste0('st@images$',  'adult.day4', '@coordinates'))) 
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
myRCTD <- create.RCTD(puck, reference, max_cores = max_cores, 
                      gene_cutoff = 0.000125, fc_cutoff = 0.5, 
                      gene_cutoff_reg = 2e-04,  
                      fc_cutoff_reg = 0.75,
                      UMI_min = 20, UMI_max = 2e+07, 
                      CELL_MIN_INSTANCE = 30)


# run RCTD main function; this could take some hours, 
# depending how many spots in the visium data, nb of celltypes and nb of cores

tic()
myRCTD <- run.RCTD(myRCTD, doublet_mode = "doublet")
saveRDS(myRCTD, file = paste0(resultsdir, '/RCTD_out_doubletMode.rds'))
toc()


##########################################
# RCTD summary plots
##########################################
myRCTD = readRDS(file = paste0(resultsdir, '/RCTD_out_doubletMode.rds'))
results <- myRCTD@results

# normalize the cell type proportions to sum to 1.
norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names

spatialRNA <- myRCTD@spatialRNA

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

##########################################
# scatterbar plot
##########################################
myRCTD = readRDS(file = paste0(resultsdir, '/RCTD_out_doubletMode.rds'))
results <- myRCTD@results

slice = "adult.day4"

# normalize the cell type proportions to sum to 1.
norm_weights = normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names

spatialRNA <- myRCTD@spatialRNA


# set color vector
getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

## Preprocess coordinates 
spatial_coord <-  spatialRNA@coords %>%
  tibble::rownames_to_column("ID")

weights = norm_weights[match(spatial_coord$ID, rownames(norm_weights)), ]

cell_types_plt = unique(colnames(weights))
col_vector = getPalette(length(cell_types_plt))
col_ct <- col_vector[seq_len(length(cell_types_plt))]

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

