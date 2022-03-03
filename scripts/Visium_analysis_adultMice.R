##########################################################################
##########################################################################
# Project: heart regeneration project
# Script purpose: First script to analyze the processed Visium data by spaceranger of 10x
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Aug 27 11:39:22 2021
##########################################################################
##########################################################################
# setup for data import and sequencing QCs
rm(list = ls())

version.analysis = '_R11934_20210827'

resDir = paste0("../results/visium_adultMice", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R11934_visium'

source('functions_Visium.R')
library(pryr) # monitor the memory usage
require(ggplot2)
mem_used()

species = 'mouse_adult'

########################################################
########################################################
# Section I: import the processed visium data by spaceranger
# first processing and QCs 
########################################################
########################################################
#design = data.frame(seq(166904, 166911), c(paste0("adult.day", c(14, 7, 4, 1)), 
#                                      paste0('neonatal.day', c(1, 4, 7, 14))), stringsAsFactors = FALSE)
#design = data.frame(seq(166907, 166904), c(paste0("adult.day", c(1, 4, 7, 14))), stringsAsFactors = FALSE)

design = data.frame(seq(166908, 166911), c(paste0('neonatal.day', c(1, 4, 7, 14))), stringsAsFactors = FALSE)

colnames(design) = c('sampleID', 'condition')

varibleGenes = c()
for(n in 1:nrow(design))
{
  # n = 1
  
  # load output from spaceranger
  aa = Seurat::Load10X_Spatial(
    data.dir = paste0(dataDir, '/output_', design$sampleID[n], '/', design$sampleID[n],  '/outs'),
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice =  design$condition[n],
    filter.matrix = TRUE,
    to.upper = FALSE
  )
  
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

# save(design, varibleGenes, st, file = paste0(RdataDir, 'seuratObject_design_variableGenes_mouse_adult.Rdata'))
#save(design, varibleGenes, st, file = paste0(RdataDir, 'seuratObject_design_variableGenes_mouse_neonadal.Rdata'))

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

SpatialFeaturePlot(st, features = 'Agrn', image.alpha = 0.6)

save(design, varibleGenes, st, 
     file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, '_umap.clustered.Rdata'))

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

########################################################
########################################################
# Section II : cell type deconvolution for each time point 
# 
########################################################
########################################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, '_umap.clustered.Rdata'))
source('functions_Visium.R')

##########################################
# spatial domain searching and potential define remote regions and border zone 
##########################################
obj.list <- SplitObject(st, split.by = "condition")
# select day4
aa = obj.list[[2]]
aa$sampleID = design$sampleID[which(design$condition == names(table(aa$condition)))]

# import manually defined spatial domain by Elad
sdomain = read.csv('/groups/tanaka/Collaborations/Jingkui-Elad/Mouse_Visium_annotations/Anno_166906.csv')
sdomain = sdomain[which(sdomain$Anno_1 != ''), ]
aa$spatial_domain_manual = NA

cells = gsub('_2_1', '',  colnames(aa))
aa$spatial_domain_manual[match(sdomain$Barcode, cells)] = sdomain$Anno_1

aa = run_bayesSpace(aa)

##########################################
# cell type deconvolution
##########################################
refs = readRDS(file = paste0('../results/Rdata/', 
                             'Seurat.obj_adultMiceHeart_Forte2020.nonCM_Ren2020CM_refCombined_cleanAnnot_logNormalize_v3.rds'))

refs$celltype[which(refs$celltype == 'immune.others')] = 'Mphage.MCT'
refs = subset(refs, cells = colnames(refs)[which(refs$celltype != 'SMC')])

st = Run.celltype.deconvolution.RCTD(st, refs)

##########################################
# cell proximity analysis 
##########################################
run_cell_proximity_analysis(aa)


########################################################
########################################################
# Section IV: spatial organization of cell types and genes  
# 
########################################################
########################################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, '_umap.clustered.Rdata'))

source('functions_Visium.R')
st = Find.SpatialDE(st)


########################################################
########################################################
# Section : cell-cell signaling pathway analysis
# 
########################################################
########################################################


