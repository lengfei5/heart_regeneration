##########################################################################
##########################################################################
# Project: Heart regeneration 
# Script purpose: analyze axolot visium data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Jan 18 10:16:57 2022
##########################################################################
##########################################################################

# setup for data import and sequencing QCs
rm(list = ls())

species = 'axolotl'

version.analysis = '_R12830_2020118'

resDir = paste0("../results/visium_axolotl", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R12830_visium/nf_out'

source('functions_Visium.R')
library(pryr) # monitor the memory usage
require(ggplot2)
mem_used()

########################################################
########################################################
# Section I: import the processed visium data by nf from Tomas and 
# first processing and QCs 
########################################################
########################################################
design = data.frame(sampleID = seq(183623, 183626), 
                    condition = c(paste0('Amex_d', c(1, 4, 7, 14))), stringsAsFactors = FALSE)
varibleGenes = c()

for(n in 1:nrow(design))
{
  # n = 1
  cat('-----------', design$condition[n], '-------------\n')
  # load nf output and process
  source('functions_Visium.R')
  aa = make_SeuratObj_visium(topdir = paste0(dataDir, '/', design$condition[n], '_', design$sampleID[n], '/'), 
                             saveDir = paste0(resDir, '/', design$condition[n], '_', design$sampleID[n], '/'))
  
  aa$condition = design$condition[n]
  #aa <- SCTransform(aa, assay = "Spatial",  method = "glmGamPoi", verbose = FALSE)
  
  aa <- SCTransform(aa, assay = "Spatial", verbose = FALSE, variable.features.n = 3000, return.only.var.genes = FALSE)
  
  test.clustering.each.condtiion = FALSE
  if(test.clustering.each.condtiion){
    aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
    ElbowPlot(aa)
    
    aa <- FindNeighbors(aa, dims = 1:10)
    aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)
    aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 30, min.dist = 0.05)
    
    DimPlot(aa, reduction = "umap", group.by = c("ident"))
    
  }
  
  varibleGenes = unique(c(varibleGenes, VariableFeatures(aa)))
  cat(design$condition[n], ' : ',  ncol(aa), ' spot found \n')
  
  # merge slices from different time points and 
  if(n == 1) {
    st = aa
  }else{
    st = merge(st, aa)
  }
  
  remove(aa)
  
}


save(design, varibleGenes, st, file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, '.Rdata'))

##########################################
# cell and gene filtering
##########################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, '.Rdata'))


##########################################
# gene and cell filtering (to add later)
##########################################
Filtering.cells.genes = FALSE
if(Filtering.cells.genes){
  #st[["percent.mt"]] <- PercentageFeatureSet(st, pattern = "^Mt-")
  # Visualize QC metrics as a violin plot
  VlnPlot(st, features = c("nCount_Spatial", "nFeature_Spatial"), ncol = 2)
  
  Idents(st) = st$condition
  FeatureScatter(st, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
  
  #plot1 <- VlnPlot(st, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  #plot2 <- SpatialFeaturePlot(aa, features = "nCount_Spatial") + theme(legend.position = "right")
  #wrap_plots(plot1, plot2)
  
}

#st = SCTransform(st, assay = "Spatial", verbose = FALSE)

DefaultAssay(st) <- "SCT"
VariableFeatures(st) <- varibleGenes

st <- RunPCA(st, verbose = FALSE)
ElbowPlot(st)

st <- FindNeighbors(st, dims = 1:20)
st <- FindClusters(st, verbose = FALSE, resolution = 0.5)

st <- RunUMAP(st, dims = 1:20, n.neighbors = 30, min.dist = 0.05)

DimPlot(st, reduction = "umap", group.by = c("ident", "condition")) 

ggsave(filename = paste0(resDir, '/UMAP_all.timepoints_', species, '.pdf'), width = 16, height = 8)

FeaturePlot(st, features = c("Myh6", 'Nppa'))

SpatialFeaturePlot(st, features = 'Myh6', image.alpha = 0.5)

save(design, varibleGenes, st, file = paste0(RdataDir, 'seuratObject_design_variableGenes_adultMice_umap.clustered.Rdata'))

QCs.with.marker.genes = FALSE
if(QCs.with.marker.genes){
  # import marker genes
  library(openxlsx)
  #aa = read.xlsx('../data/Markers_updated_v2.xlsx', sheet = 1, colNames = TRUE)
  aa = read.csv('../data/Tzahor_geneList.csv', header = TRUE)
  
  aa = aa[-c(1), c(1:2)]
  
  markers = c()
  for(n in 1:ncol(aa))
  {
    markers = c(markers, aa[!is.na(aa[,n]), n])
  }
  markers = markers[which(markers != '')]
  markers = firstup(tolower(unique(markers)))
  markers = gsub(' ', '', markers)
  
  markers[is.na(match(markers, rownames(st)))]
  
  markers = markers[!is.na(match(markers, rownames(st)))]
  
  pdfname = paste0(resDir, "/check_detected_celltypes_using_AdditionalMarkerGenes_", species, '_', colnames(aa)[1],   ".pdf")
  pdf(pdfname, width = 16, height = 8)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  for(n in 1:length(markers))
  {
    cat(markers[n], '\n')
    p3 = SpatialFeaturePlot(st, features = markers[n], image.alpha = 0.5)
    
    plot(p3)
  }
  
  dev.off()
  
}


