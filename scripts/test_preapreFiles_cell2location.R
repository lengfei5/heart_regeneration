##########################################################################
##########################################################################
# Project:
# Script purpose: prepare files for cell2location test
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Oct  3 11:15:27 2022
##########################################################################
##########################################################################

##########################################
# test cell2location: prepare loom files
##########################################
require(SeuratDisk)
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered', species, '.Rdata'))
st$condition = factor(st$condition, levels = design$condition)
st_list = Seurat::SplitObject(st, split.by = 'condition')

for(n in 1:length(st_list)){
  # n = 1
  stx = st_list[[n]]
  stx$condition = droplevels(stx$condition)
  cc = unique(stx$condition)
  DefaultAssay(stx) = 'Spatial'
  
  #stx = DietSeurat(stx, counts = TRUE, data = TRUE, scale.data = FALSE, assays = 'Spatial', graphs = NULL, images = NULL)
  stx = CreateSeuratObject(counts = GetAssayData(stx, slot = "counts"), assay = "Visium", 
                           meta.data = stx@meta.data) # create object
  
  stx <- as.loom(stx, filename = paste0("../data/visium_", cc, ".loom"), verbose = FALSE, overwrite = TRUE)
  stx
  
  stx$close_all()
  
}

# or save all visium data into one loop file
DefaultAssay(st) = 'Spatial'
st_all = CreateSeuratObject(counts = GetAssayData(st, slot = "counts"), assay = "Visium", 
                            meta.data = st@meta.data) # create object

st_all <- as.loom(st_all, filename = paste0("../data/visium_Amex_all.loom"), verbose = FALSE, overwrite = TRUE)

st_all$close_all()
