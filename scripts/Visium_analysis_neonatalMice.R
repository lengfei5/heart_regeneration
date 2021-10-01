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
version.analysis = '_R11934_20210827_neonatal'

resDir = paste0("../results/visium_mouse", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R11934_visium'

source('functions_Visium.R')

########################################################
########################################################
# Section : import the processed visium data by spaceranger
# first start with neonatal mice samples
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
  cat(design$condition[n], ' -- ', design$sampleID[n], ' :\n')
  cat(ncol(aa), ' spots ', nrow(aa), 'genes detected \n')
  
  aa$condition = design$condition[n]
  
  
  ##########################################
  # gene and cell filtering (to add later)
  ##########################################
  Filtering.cells.genes = FALSE
  if(Filtering.cells.genes){
    
    # remove some spots
    cd = aa@images$neonatal.day1@coordinates
    if(n ==1 ){
      cropped = which(cd$imagerow > 6000 & cd$imagecol > 6000)
      aa = aa[, cropped]
      cat(ncol(aa), ' spots ', nrow(aa), 'genes left after spot removal\n')
      #aa = subset(aa, neonatal.day1_imagerow > 8000 & neonatal.day1_imagecol > 8000, invert = FALSE)
    }
    
    
    aa[['percent.mt']] = PercentageFeatureSet(aa, pattern = "^Mt", assay = 'Spatial')
    
    VlnPlot(aa, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), ncol = 3)
    
    plot1 <- SpatialFeaturePlot(aa, features = "nCount_Spatial") + theme(legend.position = "bottom")
    plot2 <- SpatialFeaturePlot(aa, features = "nFeature_Spatial") + theme(legend.position = "bottom")
    plot3 <- SpatialFeaturePlot(aa, features = "percent.mt") + theme(legend.position = "bottom")
    wrap_plots(plot1, plot2, plot3)
    
    
    #st[["percent.mt"]] <- PercentageFeatureSet(st, pattern = "^Mt-")
    # Visualize QC metrics as a violin plot
    #VlnPlot(st, features = c("nCount_Spatial", "nFeature_Spatial"), ncol = 2)
    
    Idents(st) = st$condition
    FeatureScatter(st, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
    
    #plot1 <- VlnPlot(st, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
    #plot2 <- SpatialFeaturePlot(aa, features = "nCount_Spatial") + theme(legend.position = "right")
    #wrap_plots(plot1, plot2)
    
  }
  
  ##########################################
  # normalization 
  ##########################################
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
save(design, varibleGenes, st, file = paste0(RdataDir, 'seuratObject_design_variableGenes_mouse_neonadal.Rdata'))

##########################################
# cell and gene filtering
##########################################
species = 'mouse_adult'
species = 'mouse_neonadal'

#load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_mouse_adult.Rdata'))
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, '.Rdata'))



#st = SCTransform(st, assay = "Spatial", verbose = FALSE)

DefaultAssay(st) <- "SCT"
VariableFeatures(st) <- varibleGenes

st <- RunPCA(st, verbose = FALSE)
ElbowPlot(st)
st <- FindNeighbors(st, dims = 1:20)
st <- FindClusters(st, verbose = FALSE)

st <- RunUMAP(st, dims = 1:10, n.neighbors = 20, min.dist = 0.1)

DimPlot(st, reduction = "umap", group.by = c("ident", "condition"))


# import marker genes
library(openxlsx)
#aa = read.xlsx('../data/Markers_updated_v2.xlsx', sheet = 1, colNames = TRUE)
aa = read.csv('../data/Tzahor_geneList.csv', header = TRUE)

aa = aa[-c(1), c(1:2)]
#aa = aa[-c(1), -c(1:3)]

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

#markers = unique(c(markers, rownames(st)[grep('Ly6g', rownames(st))]))

#markers = c('Timp1', 'Timp2', 'Timp3', 'Timp4')

pdfname = paste0(resDir, "/check_detected_celltypes_using_AdditionalMarkerGenes_", species, '_', colnames(aa)[1],   ".pdf")
pdf(pdfname, width = 16, height = 8)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

# p1 = DimPlot(st, reduction = "umap", group.by = c("ident", "condition"))
# plot(p1)
# p2 = SpatialDimPlot(st)
# plot(p2)

for(n in 1:length(markers))
{
  cat(markers[n], '\n')
  p3 = SpatialFeaturePlot(st, features = markers[n], image.alpha = 0.5)
  
  plot(p3)
}

dev.off()


