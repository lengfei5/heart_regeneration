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

resDir = paste0("../results/visium_neonatalMice", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R11934_visium'

source('functions_Visium.R')

########################################################
########################################################
# Section I : import the processed visium data by spaceranger
# first start with neonatal mice samples
########################################################
########################################################
design = data.frame(seq(166908, 166911), c(paste0('neonatal.day', c(1, 4, 7, 14))), stringsAsFactors = FALSE)
colnames(design) = c('sampleID', 'condition')
design.adult = data.frame(seq(166907, 166904), c(paste0("adult.day", c(1, 4, 7, 14))), stringsAsFactors = FALSE)
colnames(design.adult) = colnames(design)

design = data.frame(rbind(design, design.adult))
rm(design.adult)
design$species = NA
design$species[1:4] = 'neonatal'
design$species[5:8] = 'adult'

for(n in 1:nrow(design))
{
  n = 8
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
  pdfname = paste0(resDir, '/QCs_gene_marker_check_', design$condition[n], '.pdf')
  pdf(pdfname, width=16, height = 8)
  
  #cd = aa@images$neonatal.day1@coordinates
  cd = eval(parse(text = paste0('aa@images$', design$condition[n], '@coordinates')))
  plot(cd$imagerow, cd$imagecol)
  
  Filtering.cells.genes = TRUE
  if(Filtering.cells.genes){
    
    # remove some spots that don't need to be considered 
    if(n ==1 ){
      plot(cd$imagerow, cd$imagecol)
      abline(v = 6000, col = 'red')
      abline(h = 6000, col = 'red')
      
      cropped = which(cd$imagerow > 6000 & cd$imagecol > 6000)
      aa = aa[, cropped]
      cat(ncol(aa), ' spots ', nrow(aa), 'genes left after spot removal\n')
      #aa = subset(aa, neonatal.day1_imagerow > 8000 & neonatal.day1_imagecol > 8000, invert = FALSE)
    }
    
    if(n == 2){
      plot(cd$imagerow, cd$imagecol)
      abline(v = 6000, col = 'red')
      abline(h = 15000, col = 'red')
      
      cropped = which(cd$imagecol < 15000)
      aa = aa[, cropped]
      cat(ncol(aa), ' spots ', nrow(aa), 'genes left after spot removal\n')
      
    }
    
    if(n == 3){
      plot(cd$imagerow, cd$imagecol)
      abline(v = 6000, col = 'red')
      abline(h = 15000, col = 'red')
      
      cropped = which(cd$imagerow > 6000)
      aa = aa[, cropped]
      cat(ncol(aa), ' spots ', nrow(aa), 'genes left after spot removal\n')
      
    }
    
    if(n == 4){
      plot(cd$imagerow, cd$imagecol)
      abline(v = 14500, col = 'red')
      abline(h = 15000, col = 'red')
      
      cropped = which(cd$imagerow <14500)
      aa = aa[, cropped]
      cat(ncol(aa), ' spots ', nrow(aa), 'genes left after spot removal\n')
      
    }
    
    # Cell QC metrics: percentage of Mt, nb of counts, nb of genes 
    aa[['percent.mt']] = PercentageFeatureSet(aa, pattern = "^Mt", assay = 'Spatial')
    
    Idents(aa) = design$condition[n]
    
    p1 = VlnPlot(aa, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), ncol = 3, pt.size = 1.0)
    plot(p1)
    
    p1 = FeatureScatter(aa, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
    p2 = FeatureScatter(aa, feature1 = "nCount_Spatial", feature2 = "percent.mt")
    p3 = FeatureScatter(aa, feature1 = "nFeature_Spatial", feature2 = "percent.mt")
    plot(wrap_plots(p1, p2, p3))
    
    plot1 <- SpatialFeaturePlot(aa, features = "nCount_Spatial") + theme(legend.position = "bottom")
    plot2 <- SpatialFeaturePlot(aa, features = "nFeature_Spatial") + theme(legend.position = "bottom")
    plot3 <- SpatialFeaturePlot(aa, features = "percent.mt") + theme(legend.position = "bottom")
    plot(wrap_plots(plot1, plot2, plot3))
    
    if(design$species[n] == 'neonatal'){ pct.mt.cutoff = 0.6; }
    if(design$species[n] == 'adult') { pct.mt.cutoff = 1 }
    
    aa = aa[, which(aa$percent.mt < pct.mt.cutoff)]
    
    cat(ncol(aa), ' spots left after spot filtering \n')
    
  }
  
  ##########################################
  # normalization 
  ##########################################
  #aa <- SCTransform(aa, assay = "Spatial",  method = "glmGamPoi", verbose = FALSE)
  # min_cells = 5 only use genes that have been detected in at least this many cells 
  # (determine the gene nb of output)
  aa <- SCTransform(aa, assay = "Spatial", verbose = FALSE, variable.features.n = 3000, 
                    return.only.var.genes = FALSE, 
                    min_cells=5) 
  
  aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa)
  
  aa <- FindNeighbors(aa, dims = 1:10)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
  aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 30, min.dist = 0.1)
  
  DimPlot(aa, reduction = "umap", group.by = c("ident"))
  
  # check cell cycle scores and phase assignments
  s.genes <- firstup(cc.genes.updated.2019$s.genes)
  g2m.genes <- firstup(cc.genes.updated.2019$g2m.genes)
  aa <- CellCycleScoring(aa, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
  
  #head(aa[[]])
  
  p1 = SpatialDimPlot(aa, group.by = "Phase", label.size = 7)
  p2 = DimPlot(aa, reduction = "umap", group.by = c("Phase"))
  
  plot(wrap_plots(p1, p2))
  
  examples  = c('Mki67', 'Cdk1', 'Cdk4', 'Birc5',  
                'Myh6', 'Tnnc1', 'Nppa', 'Col1a2', 'Vim', 'Emcn', 'Cd68', 'Csf1r')
  
  for(m in 1:length(examples)){
    # m = 1
    gg = examples[m]
    if(length(which(rownames(aa) == gg)) > 0){
      p1 = FeaturePlot(aa, reduction = 'umap', features = gg)
      p2 = SpatialFeaturePlot(aa, features = gg, image.alpha = 0.5)
      
      plot(wrap_plots(p2, p1))
      
    }
    
  }
  
  dev.off()
  saveRDS(aa, file = paste0(RdataDir, 'seuratObject_design_st_mouse_', design$condition[n], '.rds'))
  
}

#save(design, varibleGenes, st, file = paste0(RdataDir, 'seuratObject_design_variableGenes_mouse_adult.Rdata'))
#save(design, varibleGenes, st, file = paste0(RdataDir, 'seuratObject_design_variableGenes_mouse_neonadal.Rdata'))

##########################################
# cell and gene filtering
##########################################
species = 'mouse_neonadal'

load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, '.Rdata'))
# saveRDS(st, file = paste0(RdataDir, 'seuratObject_visium_UMIcounts_Image_neonatal.mouse_2share.rds'))

DefaultAssay(st) <- "SCT"
VariableFeatures(st) <- varibleGenes

st <- RunPCA(st, verbose = FALSE)
ElbowPlot(st)
st <- FindNeighbors(st, dims = 1:20)
st <- FindClusters(st, verbose = FALSE)

st <- RunUMAP(st, dims = 1:10, n.neighbors = 20, min.dist = 0.1)

DimPlot(st, reduction = "umap", group.by = c("ident", "condition"))

SpatialFeaturePlot(st, features = 'Cd24a', image.alpha = 0.5)
SpatialFeaturePlot(st, features = 'Agrn', image.alpha = 0.5)


# import marker genes
library(openxlsx)
#aa = read.xlsx('../data/Markers_updated_v2.xlsx', sheet = 1, colNames = TRUE)
#aa = read.csv('../data/Tzahor_geneList.csv', header = TRUE)
#aa = read.xlsx('../data/Neonate_visium_gene_Lingling.xlsx', sheet = 1, colNames = TRUE)
aa = read.xlsx('../data/Neonate_visium_GeneList_Lingling_20230226.xlsx', sheet = 1, colNames = TRUE)

#aa = aa[-c(1), c(1:2)]
#aa = aa[-c(1), -c(1:3)]

markers = aa$Gene.Name
markers = markers[!is.na(markers)]
#for(n in 1:ncol(aa))
#{
#  markers = c(markers, aa[!is.na(aa[,n]), 1])
#}
markers = markers[which(markers != '')]
markers = firstup(tolower(unique(markers)))
markers = gsub(' ', '', markers)

markers[is.na(match(markers, rownames(st)))]

markers = markers[!is.na(match(markers, rownames(st)))]
#markers = unique(c(markers, rownames(st)[grep('Ly6g', rownames(st))]))
#markers = c('Timp1', 'Timp2', 'Timp3', 'Timp4')

#markers = c('F13a1', 'Apoe', 'Igfbp4', 'C1qa', 'Vldlr', 'Lrp5', 'Lrp6', 'Itga9', 'Cspg4', 'Cr1l')
#markers = c("Ntn1", "Mfap5", "Ubb", "Spon2", "Sparc", "Comp", "Ncanb", "Cthrc1", "Gpr42", "Ffr2")

mm = match(markers, rownames(st))
markers[is.na(mm)]

pdfname = paste0(resDir, "/Neondal_visium_4Lingling_v20230226.pdf")
pdf(pdfname, width = 8, height = 6)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

for(n in 1:length(markers))
#for(n in 1:10)
{
  # n = 1
  cat(n, '--', markers[n])
  if(length(which(rownames(st) == markers[n]))==1){
    cat('\n')
    p3 = SpatialFeaturePlot(st, features = markers[n], image.alpha = 0., 
                            images = "neonatal.day4",
                            pt.size.factor = 2.5,
                            slot = "data",
                            #min.cutoff = 'q1', 
                            max.cutoff = 'q98')
    #+
    #ggplot2::scale_fill_continuous(limits = c(0.0,1.0),breaks = c(0.0, 0.5, 1.0))
    plot(p3)
    
  }else{
    cat(' NOT FOUND\n')
  }
}

dev.off()

## test two features, there is no support for blend option
# SpatialFeaturePlot(object = st,
#             features = c("Egr1", "Col1a1"),
#             image.alpha = 0., 
#             images = "neonatal.day4",
#             pt.size.factor = 2.5,
#             slot = "data",
#             #min.cutoff = 'q1', 
#             max.cutoff = 'q98',
#             blend = TRUE, 
#             blend.threshold = 1, 
#             order = T)

