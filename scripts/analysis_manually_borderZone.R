##########################################################################
##########################################################################
# Project:
# Script purpose: manually define the injury, border zones and remote regions
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Oct  4 15:03:48 2022
##########################################################################
##########################################################################
rm(list =ls())

require(ggplot2)
require(Seurat)
require(SPATA2) # installation https://themilolab.github.io/SPATA2/articles/spata-v2-installation.html

load(file = paste0('/groups/tanaka/Collaborations/Jingkui-Elad/visium_axolotl_reseq/',
                  'seuratObject_design_variableGenes_umap.clusteredaxolotl.Rdata'))
st$condition = factor(st$condition, levels = design$condition)

outDir = paste0('/groups/tanaka/Collaborations/Jingkui-Elad/visium_axolotl_reseq/', 
                'spata2_manual_regions/')

cat('-- check visium conditions -- \n')
print(table(st$condition))
cc = names(table(st$condition))
system(paste0('mkdir -p ', outDir))

for(n in 1:length(cc))
{
  # n = 1
  cat('slice -- ', cc[n], '\n')
  slice = cc[n]
  aa = st[, which(st$condition == slice)]
  DefaultAssay(aa) = 'Spatial'
  
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
  
  spata_obj = transformSeuratToSpata(aa, sample_name = slice, method = 'spatial', 
                                     assay_name = 'Spatial', 
                                     assay_slot = 'scale.data', 
                                     image_name = slice)
  
  #setActiveExpressionMatrix(spata_obj, 'data')
  spata_obj <- createSegmentation(object = spata_obj)
  
  plotSegmentation(object = spata_obj, pt_size = 1.9) +
    ggplot2::scale_y_reverse()
  
  # save the segment in the object of each slice
  aa$segmentation = NA
  segments = as.character(levels(getSegmentNames(spata_obj)))
  segments = segments[which(segments != 'none')]
  
  for(seg in segments){
    cat(seg, '\n')
    aa$segmentation[match(getSegmentDf(spata_obj, segment_names = seg)$barcodes, colnames(aa))] = seg
  }
  
  saveRDS(aa, file = paste0(outDir, 'visium_manual_segmentation_', slice, '.rds'))
  
}

