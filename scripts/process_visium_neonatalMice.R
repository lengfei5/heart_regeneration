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
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
library(circlize)
library(RColorBrewer)
require(scran)
require(scater)
source('functions_scRNAseq.R')
source('functions_Visium.R')
source('functions_cccInference.R')

dataPath_nichenet = '../data/NicheNet/'

version.analysis = '_R11934_20210827_neonatal'

resDir = paste0("../results/visium_neonatalMice", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R11934_visium_mice'

species = 'mouse_neonadal'


########################################################
########################################################
# Section I : import the processed visium data by spaceranger
# first start with neonatal mice samples
########################################################
########################################################
design = data.frame(seq(166908, 166911), 
                    c(paste0('neonatal.day', c(1, 4, 7, 14))), 
                    stringsAsFactors = FALSE)

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
  pdfname = paste0(resDir, '/QCs_gene_marker_check_', design$condition[n], '.pdf')
  pdf(pdfname, width=16, height = 8)
  
  #cd = aa@images$neonatal.day1@coordinates
  cd = eval(parse(text = paste0('aa@images$', design$condition[n], '@coordinates')))
  plot(cd$imagerow, cd$imagecol, main = 'before manual crop')
  
  Manually_crop_images = TRUE
  if(Manually_crop_images){
    
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
  
  saveRDS(aa, file = paste0(RdataDir, 'seuratObject_design_st_mouse_', design$condition[n], '.rds'))
  
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

save(st, design, varibleGenes, file = paste0(RdataDir, 'seuratObject_cell.gene.filtering_design_variableGenes_', 
                                             species, '.Rdata'))

##########################################
# filtering cells and genes 
##########################################
load(file = paste0(RdataDir, 'seuratObject_cell.gene.filtering_design_variableGenes_', 
                   species, '.Rdata'))
table(st$condition)

Idents(st) = factor(st$condition, 
                    levels = c("neonatal.day1", "neonatal.day4", "neonatal.day7", "neonatal.day14"))

p1 = VlnPlot(st, features = c("nCount_Spatial"), pt.size = 1.0) + NoLegend() +
  geom_hline(yintercept = 500)

p2 = VlnPlot(st, features = c("nFeature_Spatial"), pt.size = 1.0) + NoLegend() +
  geom_hline(yintercept = 300)

p1 / p2

DefaultAssay(st) <- "SCT"

SpatialFeaturePlot(st, features = 'Cd24a', image.alpha = 0.5)
SpatialFeaturePlot(st, features = 'Agrn', image.alpha = 0.5)


##########################################
# make plots for Lingling 
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

########################################################
########################################################
# Section II: filtering in the same way as adult mice 
# and cell type deconvolution
########################################################
########################################################
##########################################
# cell filtering and gene filtering (ribosome and mt genes)
##########################################
Filtering.cells.genes = FALSE
if(Filtering.cells.genes){
  load(file = paste0(RdataDir, 'seuratObject_cell.gene.filtering_design_variableGenes_', 
                     species, '.Rdata'))
  table(st$condition)
  
  Idents(st) = factor(st$condition, 
                      levels = c("neonatal.day1", "neonatal.day4", "neonatal.day7", "neonatal.day14"))
  
  st[["percent.mt"]] <- PercentageFeatureSet(st, pattern = "^mt-")
  
  # Visualize QC metrics as a violin plot
  VlnPlot(st, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), ncol = 3)
  
  cells =  colnames(st)[which(st$nCount_Spatial > 500 & st$nFeature_Spatial > 300)]
  ggs = rownames(st)[grep('^Rp[sl]|^mt-', rownames(st), invert = TRUE)]
  
  st = subset(st, cells = cells, features = ggs)
  
  cat(ncol(st), ' spots after cell filtering \n')
  cat(nrow(st), ' genes left after MT and Ribo filtering \n')
  
  ## redo the SCT nomralization 
  st <- SCTransform(st, assay = "Spatial", verbose = FALSE, variable.features.n = 3000, 
                    return.only.var.genes = FALSE, 
                    min_cells = 5
                    ) 
  
  st <- RunPCA(st, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(st)
  
  st <- FindNeighbors(st, dims = 1:10)
  st <- FindClusters(st, verbose = FALSE, algorithm = 3, resolution = 0.5)
  st <- RunUMAP(st, dims = 1:10, n.neighbors = 30, min.dist = 0.1)
  
  DimPlot(st, reduction = "umap", group.by = c("ident"))
  DimPlot(st, reduction = "umap", group.by = c("condition"))
  
  # check cell cycle scores and phase assignments
  s.genes <- firstup(cc.genes.updated.2019$s.genes)
  g2m.genes <- firstup(cc.genes.updated.2019$g2m.genes)
  st <- CellCycleScoring(st, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
  
  
  saveRDS(st, file = paste0(RdataDir,
                            'seuratObject_neonatalMouse_cell.gene.filtered_umap.clustered.rds'))
  
}


##########################################
# cell type deconvolution for subtypes
##########################################
source('functions_Visium.R')

st = readRDS(file = paste0(RdataDir,
                           'seuratObject_neonatalMouse_cell.gene.filtered_umap.clustered.rds'))



refs = readRDS(file = paste0('../data/data_examples/ref_scRNAseq_neonatalMice_clean.v1.2.rds'))

refs$subtype = gsub(' ','.', refs$subtype)
refs$subtype = gsub('_','.', refs$subtype)
refs$subtype = gsub('-','.', refs$subtype)

p1 = DimPlot(refs, reduction = 'umap', group.by = 'dataset',  label = TRUE, repel = TRUE)
p2 = DimPlot(refs, reduction = 'umap', group.by = 'celltype',  label = TRUE, repel = TRUE)
p3 = DimPlot(refs, reduction = 'umap', group.by = 'subtype',  label = TRUE, repel = TRUE)
p4 = DimPlot(refs, reduction = 'umap', group.by = 'timepoints',  label = TRUE, repel = TRUE)

(p1 + p2)/(p3 + p4)

ggsave(filename = paste0(resDir, '/UMAP_scRNAseq_refrence_dataset_timepoints_celltypes.pdf'), 
       width = 24, height = 18)

Use_fineGrained_subtypes = FALSE
if(Use_fineGrained_subtypes){
  
  refs$celltype_toUse = as.character(refs$subtype)
  
  length(table(refs$celltype_toUse))
  table(refs$celltype_toUse)
  
  DimPlot(refs, reduction = 'umap', group.by = 'celltype_toUse')
  
  ## prepare the celltype to use and also specify the time-specific subtypes
  table(refs$condition)
  
  table(refs$celltype_toUse)
  length(table(refs$celltype_toUse))
  st$condition = factor(st$condition)
  table(st$condition)
  
  ## preapre the paramters for RCTD subtypes
  DefaultAssay(refs) = 'integrated'
  DefaultAssay(st) = 'Spatial'
  require_int_SpatialRNA = FALSE
  
  condition.specific.ref = FALSE
  
  outDir = paste0(resDir, '/celltype_deconvolution')
  RCTD_out = paste0(outDir, '/RCTD_', length(table(refs$celltype_toUse)), 'Subtype_ref_v0.1')
  
  max_cores = 16
  
  Idents(st) = st$condition
  st = subset(st, idents = c('neonatal.day4', 'neonatal.day7')); 
  st$condition = droplevels(st$condition)
  
  source('functions_Visium.R')
  Run.celltype.deconvolution.RCTD(st, refs, 
                                  condition.specific.ref = condition.specific.ref,
                                  #condition.specific_celltypes = condition.specific_celltypes,
                                  require_int_SpatialRNA = require_int_SpatialRNA,
                                  max_cores = max_cores,
                                  RCTD_out = RCTD_out,
                                  plot.RCTD.summary = FALSE, 
                                  PLOT.scatterpie = FALSE
                                  
  )
  
  source('functions_Visium.R')
  plot.RCTD.results(st = st, 
                    RCTD_out = RCTD_out,
                    species = species,
                    plot.RCTD.summary = FALSE)
  
  
}


