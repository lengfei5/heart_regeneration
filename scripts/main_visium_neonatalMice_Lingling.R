##########################################################################
##########################################################################
# Project: heart regeneration project
# Script purpose: analyze the processed Visium data by spaceranger of 10x
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Aug 27 11:39:22 2021
##########################################################################
##########################################################################
rm(list = ls())

library(pryr) # monitor the memory usage
require(ggplot2)

library(Seurat) # update to Seurat V5

library(tidyverse)
library(circlize)
library(RColorBrewer)
#require(scran)
#require(scater)
#library(nichenetr)
source('functions_scRNAseq.R')
source('functions_Visium.R')
source('functions_cccInference.R')

version.analysis = '_R11934_20210827_neonatal'
resDir = paste0("../results/visium_neonatalMice", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R11934_visium_mice'

species = 'mouse_neonadal'

design = data.frame(seq(166908, 166911), 
                    c(paste0('neonatal.day', c(1, 4, 7, 14))), 
                    stringsAsFactors = FALSE)

colnames(design) = c('sampleID', 'condition')

########################################################
########################################################
# Section I : import the processed visium data by spaceranger
# first start with neonatal mice samples
########################################################
########################################################
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
  pdfname = paste0(resDir, '/QCs_gene_marker_check_', design$condition[n], '.pdf')
  pdf(pdfname, width=16, height = 8)
  
  #cd = aa@images$neonatal.day1@coordinates
  cd = eval(parse(text = paste0('aa@images$', design$condition[n], '@coordinates')))
  plot(cd$imagerow, cd$imagecol, main = 'before manual crop')
  
  Manually_crop_images = FALSE
  if(Manually_crop_images){
    
    cat('manually crop images \n')
    
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
    
  }
  
  cd = eval(parse(text = paste0('aa@images$', design$condition[n], '@coordinates')))
  plot(cd$imagerow, cd$imagecol, main = 'after manual crop')
  
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
  
  p1 = VlnPlot(aa, features = c("nCount_Spatial"), pt.size = 1.0) + NoLegend() +
    geom_hline(yintercept = 500)
  
  p2 = VlnPlot(aa, features = c("nFeature_Spatial"), pt.size = 1.0) + NoLegend() +
    geom_hline(yintercept = 300)
  
  plot(wrap_plots(p1, p2))
  
  sels = which(aa$nCount_Spatial > 500 & aa$nFeature_Spatial > 300)
  aa= aa[, sels]
  cat(ncol(aa), ' spots after cell filtering \n')
  
  ##########################################
  # normalization 
  ##########################################
  # min_cells = 5 only use genes that have been detected in at least this many cells 
  aa <- SCTransform(aa, assay = "Spatial", verbose = FALSE, variable.features.n = 3000, 
                    return.only.var.genes = FALSE, 
                    min_cells = 5)
  
  aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa)
  
  aa <- FindNeighbors(aa, dims = 1:10)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
  aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 30, min.dist = 0.1)
  
  DimPlot(aa, reduction = "umap", group.by = c("ident"))
  
  saveRDS(aa, file = paste0(RdataDir, 'seuratObject_design_st_mouse_', design$condition[n], '_nocropping.rds'))
  
  varibleGenes = unique(c(varibleGenes, VariableFeatures(aa)))
  
  # merge slices from different time points and 
  if(n == 1) {
    st = aa
  }else{
    st = merge(st, aa)
  }
  
  remove(aa)
  
  dev.off()
  
}

save(st, design, varibleGenes, file = paste0(RdataDir, 
                                             'seuratObject_cell.gene.filtering_design_variableGenes_nocropping_', 
                                             species, '.Rdata'))

##########################################
# crop the images  
##########################################
Crop_images_test = FALSE
if(Crop_images_test){
  
  load(file = paste0(RdataDir, 
                     'seuratObject_cell.gene.filtering_design_variableGenes_nocropping_', 
                     species, '.Rdata'))
  
  cd = st@images$neonatal.day14@coordinates
  plot(cd$imagerow, cd$imagecol, main = 'before manual crop')
  
  Idents(st) = factor(st$condition, levels = c("neonatal.day1", "neonatal.day4", 
                                               "neonatal.day7", 
                                               "neonatal.day14"))
  
  DefaultAssay(st) <- "SCT"
  
  SpatialFeaturePlot(st, features = 'Tnni3', 
                     crop = FALSE, 
                     image.alpha = 0.1, 
                     #images = "neonatal.day4",
                     pt.size.factor = 1,
                     slot = "data",
                     keep.scale = 'feature',
                     min.cutoff = 'q1', 
                     max.cutoff = 'q98'
  )
  
  a1 = readRDS(paste0(RdataDir, 'seuratObject_design_st_mouse_neonatal.day1_nocropping.rds'))
  #ImageFeaturePlot(a1, features = 'Tnni3')
  
  library(ggplot2)
  library(patchwork)
  library(semla)
  library(SeuratData)
  library(STutility)
  
  a1 = GetStaffli(a1)
  a1<- CropImages(a1, crop.geometry.list = list("1" = "500x500+500+500"))
  
  brain <- LoadData("stxBrain", type = "anterior1")
  # Make Seurat object compatible with semla
  brain_semla <- UpdateSeuratForSemla(brain)
  
  
  # Make Seurat object compatible with semla
  xx <- UpdateSeuratForSemla(a1, image_type = 'tissue_lowres')
  
  cd = a1@images$neonatal.day1@coordinates
  plot(cd$imagerow, cd$imagecol, main = 'before manual crop')
  
  cropped = which(cd$imagerow > 6000 & cd$imagerow < 16000 & cd$imagecol >6000 & cd$imagecol <16000)
  a1 = a1[, cropped]
  
  cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
  cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
  cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)
  
  coord <- GetTissueCoordinates(object = a1@images$neonatal.day1)
  # calculate the aspect ratio of rows to columns
  myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
  # force the image into the right aspect ratio
  #SpatialDimPlot(a1, crop = TRUE, pt.size.factor = 5) + theme(aspect.ratio = myratio)
  
  SpatialFeaturePlot(a1, features = 'Tnni3', 
                     crop = FALSE, 
                     image.alpha = 0.1, 
                     #images = "neonatal.day4",
                     pt.size.factor = 1,
                     slot = "data",
                     keep.scale = 'feature',
                     min.cutoff = 'q1', 
                     max.cutoff = 'q98'
  ) + theme(aspect.ratio = 0.5)
  
  a2 = readRDS(paste0(RdataDir, 'seuratObject_design_st_mouse_neonatal.day14_nocropping.rds'))
  cd = a2@images$neonatal.day1@coordinates
  plot(cd$imagerow, cd$imagecol, main = 'before manual crop')
  
  coord <- GetTissueCoordinates(object = a2@images$neonatal.day1)
  # calculate the aspect ratio of rows to columns
  myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
  # force the image into the right aspect ratio
  SpatialDimPlot(a2, crop = TRUE) + theme(aspect.ratio = myratio)
  
  
  cropped = which(cd$imagerow > 6000 & cd$imagerow < 16000 & cd$imagecol >6000 & cd$imagecol <16000)
  a2 = a2[, cropped]
  
  SpatialFeaturePlot(a2, features = 'Tnni3', 
                     crop = TRUE, 
                     image.alpha = 0.1, 
                     #images = "neonatal.day4",
                     pt.size.factor = 1,
                     slot = "data",
                     keep.scale = 'feature',
                     min.cutoff = 'q1', 
                     max.cutoff = 'q98'
  )
  
  #xx = Crop(a1[['neonatal.day14']], x = c(1750, 3000), y = c(3750, 5250), coords = "tissue")
 
  xx = merge(a1, a2)
  
  SpatialFeaturePlot(xx, features = 'Tnni3', 
                     crop = TRUE, 
                     image.alpha = 0.1, 
                     #images = "neonatal.day4",
                     pt.size.factor = 1,
                     slot = "data",
                     keep.scale = 'feature',
                     min.cutoff = 'q1', 
                     max.cutoff = 'q98'
  )
  
}


########################################################
########################################################
# Section II: make plots for selected marker genes
# 
########################################################
########################################################
load(file = paste0(RdataDir, 'seuratObject_cell.gene.filtering_design_variableGenes_', 
                   species, '.Rdata'))

#load(file = paste0(RdataDir, 'seuratObject_no.cell.gene.filtering_design_variableGenes_mouse_neonadal.Rdata'))
Idents(st) = factor(st$condition, levels = c("neonatal.day1", "neonatal.day4", 
                                             "neonatal.day7", 
                                             "neonatal.day14"))

DefaultAssay(st) <- "SCT"

SpatialFeaturePlot(st, features = 'Tnni3', 
                   crop = FALSE, 
                   image.alpha = 0.1, 
                   #images = "neonatal.day4",
                   pt.size.factor = 1.0,
                   slot = "data",
                   keep.scale = 'feature',
                   min.cutoff = 'q1', 
                   max.cutoff = 'q98'
)


VlnPlot(st, features = c("nCount_Spatial"), pt.size = 1.0) + NoLegend() +
  geom_hline(yintercept = 500)

VlnPlot(st, features = c("nFeature_Spatial"), pt.size = 1.0) + NoLegend() +
  geom_hline(yintercept = 300)

DefaultAssay(st) <- "SCT"

SpatialFeaturePlot(st, features = 'Cd24a', image.alpha = 0.5)
SpatialFeaturePlot(st, features = 'Agrn', image.alpha = 0.5)


##########################################
# make plots for Lingling with day4
##########################################
# import marker genes
#library(openxlsx)
#aa = read.xlsx('../data/Neonate_visium_GeneList_Lingling_20230226.xlsx', sheet = 1, colNames = TRUE)

markers = c('Vim', 'Vcan', 'Tnni3', 'Tnc', 'Thbs1', 'Ptgis', 'Pdha1', 
            'Nppa', 'Myh6', 'Mpc2', 'Lum', 'Lgmn', 'Lgals3', 'Lamp1', 
            'Itgb5', 'Itgb1', 'Hexb', 'Hexa', 'Gm2a', 'Fn1', 'Egr1', 'Dlat', 
            'Csrp2', 'Col1a1', 'Cdkn1a', 'C1qa', 'C1qb', 'Bgn', 'Anxa2')


pdfname = paste0(resDir, "/Neondal_visium_4Lingling.manuscript_v2.pdf")
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

##########################################
# make plots for Lingling across all time points
##########################################
# import marker genes
library(openxlsx)
markers = read.xlsx('../data/2024_03_12_geneList_forVisium_Lingling_Neonates.xlsx', 
                    sheet = 1, colNames = TRUE)

markers = markers$Gene.list
mm = match(markers, rownames(st))

# markers = c('Vim', 'Vcan', 'Tnni3', 'Tnc', 'Thbs1', 'Ptgis', 'Pdha1', 
#             'Nppa', 'Myh6', 'Mpc2', 'Lum', 'Lgmn', 'Lgals3', 'Lamp1', 
#             'Itgb5', 'Itgb1', 'Hexb', 'Hexa', 'Gm2a', 'Fn1', 'Egr1', 'Dlat', 
#             'Csrp2', 'Col1a1', 'Cdkn1a', 'C1qa', 'C1qb', 'Bgn', 'Anxa2')

SCT.normalizaiton.allSlice = FALSE
if(SCT.normalizaiton.allSlice){
  st <- SCTransform(st, assay = "Spatial", verbose = FALSE, variable.features.n = 3000, 
                    return.only.var.genes = FALSE,
                    min_cells = 5)
  
}


pdfname = paste0(resDir, "/Neondal_visium_4Lingling.manuscript_20240312_normalization.all.slice_nocropping.pdf")
pdf(pdfname, width = 25, height = 6)
par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)

for(n in 1:length(markers))
#for(n in 1:10)
{
  # n = 1
  cat(n, '--', markers[n])
  if(length(which(rownames(st) == markers[n]))==1){
    cat('\n')
    p3 = SpatialFeaturePlot(st, features = markers[n], image.alpha = 0.1,
                            crop = FALSE,
                            #images = "neonatal.day4",
                            pt.size.factor = 1.0,
                            slot = "data",
                            keep.scale = 'feature',
                            min.cutoff = 'q1', 
                            max.cutoff = 'q98'
                            ) 
    
    plot(p3)
    
  }else{
    cat(' NOT FOUND\n')
  }
}

dev.off()


##########################################
# extract gene expression in the injuried regions 
##########################################
Import.manual.spatial.domains = TRUE
if(Import.manual.spatial.domains){
  #require(SPATA2) # installation https://themilolab.github.io/SPATA2/articles/spata-v2-installation.html
  
  st$segmentation = NA
  
  #sdomain = read.csv('/groups/tanaka/Collaborations/Jingkui-Elad/Mouse_Visium_annotations/Anno_166906.csv')
  inputDir = '/groups/tanaka/Collaborations/Jingkui-Elad/Mouse_Visium_annotations/'
  cc = design$condition
  
  #manual_selection_spots_image_Spata
  for(n in 1:length(cc))
  {
    # n = 1
    slice = cc[n]
    cat('slice -- ', slice, '\n')
    
    sdomain = read.csv(file = paste0(inputDir, 'Anno_', design$sampleID[which(design$condition == slice)], 
                                     '.csv'))
    
    colnames(sdomain) = c('Barcode', 'Anno')
    sdomain = sdomain[which(sdomain$Anno != ''), ]
    table(sdomain$Anno)
    
    kk = which(st$condition == slice)
    cells_slice = sapply(colnames(st)[kk], function(x){unlist(strsplit(as.character(x), '_'))[1]})
    
    mm = match(sdomain$Barcode, cells_slice)
    head(mm)
    cat(length(which(is.na(mm))), ' cells with no mapping out of ', length(mm),  ' \n ')
    
    st$segmentation[kk[mm[which(!is.na(mm))]]] = sdomain$Anno[which(!is.na(mm))]
    
    rm(sdomain)
    
  }
  
  st$segmentation = gsub('infarct', 'Infarct', st$segmentation)
  st$segmentation = gsub('Inj', 'Injury', st$segmentation)
  st$segmentation = gsub('Injuryury', 'Injury', st$segmentation)
  st$segmentation = gsub('Infarct', 'Injury', st$segmentation)
  st$segmentation = gsub('remote', 'Remote_2', st$segmentation)
  st$segmentation = gsub('Remote_3_septum', 'Remote_3', st$segmentation)
  
  st$segmentation = as.factor(st$segmentation)
  Idents(st) = as.factor(st$segmentation)
  SpatialDimPlot(st, image.alpha = 0.3, pt.size.factor = 1.0, label.size = 10, crop = FALSE)
  
  ggsave(paste0(resDir, '/Manual_segmentation_spata2_Elad_v2.pdf'), width = 16, height = 6)
  
  st$segmentation[which(st$segmentation != 'Injury')] = NA
  st$segmentation = droplevels(st$segmentation)
  Idents(st) = as.factor(st$segmentation)
  SpatialDimPlot(st, image.alpha = 0.3, pt.size.factor = 1.0, label.size = 10, crop = FALSE)
  ggsave(paste0(resDir, '/Manual_InjuriedZone_Elad_noCropping.pdf'), width = 16, height = 6)
  
  #save(st, design, file = paste0(RdataDir,
  #                               'seuratObject_mouse_adult_cell.gene.filtered_umap.clustered_manualSegmentation', 
  #                               species, '.Rdata'))
  
  
}

load(file = paste0(RdataDir,
                  'seuratObject_mouse_adult_cell.gene.filtered_umap.clustered_manualSegmentation', 
                   species, '.Rdata'))

mm = match(markers, rownames(st))
xx = st@assays$SCT@data
xx = xx[mm, ]
xx = t(xx)
xx = data.frame(xx)

xx = data.frame(condition = st$condition,  segmentation = st$segmentation,  xx, stringsAsFactors = FALSE)
xx = xx[which(xx$segmentation == 'Injury'), ]

write.csv2(xx, file = paste0(resDir, '/markGenes_expression_injuriedZones_acrossTimePoint.csv'),
           quote = FALSE)

