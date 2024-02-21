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
dataPath_nichenet = '../data/NicheNet/'

version.analysis = '_R11934_20240131'

resDir = paste0("../results/visium_adultMice", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

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
design = data.frame(seq(166907, 166904), 
                    c(paste0("adult.day", c(1, 4, 7, 14))), 
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
# saveRDS(st, file = paste0(RdataDir, 'seuratObject_visium_UMIcounts_Image_', species, '_2share.rds'))

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

## reload the seurat object
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, '_umap.clustered.Rdata'))

make_plots_4Lingling = FALSE
if(make_plots_4Lingling){
  
  DefaultAssay(st) <- "SCT"
  
  # import marker genes
  library(openxlsx)
  #aa = read.xlsx('../data/Markers_updated_v2.xlsx', sheet = 1, colNames = TRUE)
  #aa = read.csv('../data/Tzahor_geneList.csv', header = TRUE)
  aa = read.xlsx('../data/Neonate_visium_gene_Lingling.xlsx', sheet = 1, colNames = TRUE)
  #aa = aa[-c(1), c(1:2)]
  #aa = aa[-c(1), -c(1:3)]
  
  markers = aa$Gene.names
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
  
  pdfname = paste0(resDir, "/adultMice_visium_4Lingling_v2.pdf")
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
                              images = "adult.day4",
                              pt.size.factor = 1.6,
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
  
}

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
  
  markers = c('F13a1', 'Apoe', 'Igfbp4', 'C1qa', 'Vldlr', 'Lrp5', 'Lrp6', 'Itga9', 'Cspg4', 'Cr1l')
  
  markers = c("Ntn1", "Mfap5", "Ubb", "Spon2", "Sparc", "Comp", "Ncanb", "Cthrc1", "Gpr42", "Ffr2")
  mm = match(markers, rownames(st))
  markers[is.na(mm)]
  
  pdfname = paste0(resDir, "/Visium_AdditionalMarkerGenes_Bern_", species, ".pdf")
  pdf(pdfname, width = 16, height = 8)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  for(n in 1:length(markers))
  {
    if(length(which(rownames(st) == markers[n]))==1){
      cat(markers[n], '\n')
      p3 = SpatialFeaturePlot(st, features = markers[n], image.alpha = 0.5)
      plot(p3)
    }else{
      cat(markers[n], ' NOT FOUND\n')
    }
  }
  
  dev.off()
  
}

##########################################
# cell filtering and gene filtering (ribosome and mt genes)
##########################################
Filtering.cells.genes = FALSE
if(Filtering.cells.genes){
  
  load(file = paste0("../results/visium_adultMice_R11934_20210827/Rdata/", 
                     'seuratObject_design_variableGenes_mouse_adult_umap.clustered.Rdata'))
  
  Idents(st) = factor(st$condition, levels = c('adult.day1', 'adult.day4', 'adult.day7', 'adult.day14'))
  
  st[["percent.mt"]] <- PercentageFeatureSet(st, pattern = "^mt-")
  
  # Visualize QC metrics as a violin plot
  VlnPlot(st, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), ncol = 3)
  
  
  FeatureScatter(st, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
  
  #plot1 <- VlnPlot(st, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  #plot2 <- SpatialFeaturePlot(aa, features = "nCount_Spatial") + theme(legend.position = "right")
  #wrap_plots(plot1, plot2)
  
  cells =  colnames(st)[which(st$nCount_Spatial > 500 & st$nFeature_Spatial > 300)]
  ggs = rownames(st)[grep('^Rp[sl]|^mt-', rownames(st), invert = TRUE)]
  
  st = subset(st, cells = cells, features = ggs)
  
  cat(ncol(st), ' spots after cell filtering \n')
  cat(nrow(st), ' genes left after MT and Ribo filtering \n')
  
  ## redo the SCT nomralization 
  st <- SCTransform(st, assay = "Spatial", verbose = FALSE, variable.features.n = 3000, 
                    return.only.var.genes = FALSE, 
                    min_cells=5 ) 
  
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
  
  
  saveRDS(st, file = paste0("../results/visium_adultMice_R11934_20210827/Rdata/",
                            'seuratObject_mouse_adult_cell.gene.filtered_umap.clustered.rds'))
  
}


########################################################
########################################################
# Section II : cell type deconvolution for each time point 
# 
########################################################
########################################################
source('functions_Visium.R')

st = readRDS(file = paste0("../results/visium_adultMice_R11934_20210827/Rdata/",
                           'seuratObject_mouse_adult_cell.gene.filtered_umap.clustered.rds'))

refs = readRDS(file = paste0('../data/data_examples/ref_scRNAseq_adultMice_clean.v1.rds'))

#jj = which(refs$dataset == 'Ren2020')
#refs$timepoints[jj] = refs$condition[jj]
#refs$condition = refs$timepoints

p1 = DimPlot(refs, reduction = 'umap', group.by = 'dataset',  label = TRUE, repel = TRUE)
p2 = DimPlot(refs, reduction = 'umap', group.by = 'celltype',  label = TRUE, repel = TRUE)
p3 = DimPlot(refs, reduction = 'umap', group.by = 'subtype',  label = TRUE, repel = TRUE)
p4 = DimPlot(refs, reduction = 'umap', group.by = 'timepoints',  label = TRUE, repel = TRUE)

(p1 + p2)/(p3 + p4)

ggsave(filename = paste0(resDir, '/UMAP_scRNAseq_refrence_dataset_timepoints_celltypes.pdf'), 
       width = 24, height = 12)

##########################################
# cell type deconvolution for cell types
##########################################
Use_coarse_celltypes = FALSE
if(Use_coarse_celltypes){
  refs$celltype_toUse = as.character(refs$celltype)
  length(table(refs$celltype_toUse))
  table(refs$celltype_toUse)
  DimPlot(refs, reduction = 'umap', group.by = 'celltype')
  
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
  RCTD_out = paste0(outDir, '/RCTD_9celltypes_ref_v0.1')
  max_cores = 16
  # st = subset(st, condition == 'Amex_d4')
  
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
  
  #saveRDS(condition.specific_celltypes, paste0(RdataDir, 'RCTD_refs_condition_specificity.rds'))
  #saveRDS(refs, file = paste0(RdataDir, 'RCTD_refs_subtypes_final_20221117.rds'))
  
  source('functions_Visium.R')
  plot.RCTD.results(st = st, 
                    RCTD_out = RCTD_out,
                    species = species,
                    plot.RCTD.summary = FALSE)
  
}

##########################################
# cell type deconvolution for subtypes
##########################################
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
  RCTD_out = paste0(outDir, '/RCTD_', length(table(refs$celltype_toUse)), 'Subtype_ref_v1.0')
  max_cores = 16
  
  # st = subset(st, condition == 'adult.day7'); st$condition = droplevels(st$condition)
  
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

########################################################
########################################################
# Section III: spatial organization of cell types and genes  
# 
########################################################
########################################################
Spatial_variableGenes_analysis = FALSE
if(Spatial_variableGenes_analysis){
  load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, 
                     '_umap.clustered.Rdata'))
  
  source('functions_Visium.R')
  st = Find.SpatialDE(st)
  
}

########################################################
########################################################
# Section IV : cell-to-cell communication analysis
# 
# 1) step: define spatial domain 
# 2) step: neighborhood enrichment analysis
# 3) step: ligand-receptor analysis
########################################################
########################################################
source('functions_Visium.R')
st = readRDS(file = paste0("../results/visium_adultMice_R11934_20210827/Rdata/",
                           'seuratObject_mouse_adult_cell.gene.filtered_umap.clustered.rds'))

cc = c('adult.day1', 'adult.day4', 'adult.day7', 'adult.day14')
st$condition = factor(st$condition, levels = cc)

cat('visium conditions :\n')
print(table(st$condition))


VlnPlot(st, features = 'nFeature_Spatial', group.by = 'condition') +
  geom_hline(yintercept = c(300, 500, 1000, 2000))

ggsave(paste0(resDir, '/QCs_nFeatures_ST_perCondition_afterFiltering.pdf'), width = 12, height = 8)

VlnPlot(st, features = 'nFeature_SCT', group.by = 'condition') +
  geom_hline(yintercept = c(500, 1000, 2000))

ggsave(paste0(resDir, '/QCs_nFeatures_ST_perCondition_afterFiltering.pdf'), width = 12, height = 8)


##########################################
# step 1) Spatial domain searching and potential define remote regions and border zone
# here using computational methods to define regions of interest or cell niches
##########################################
## import manually defined spatial domains
Import.manual.spatial.domains = TRUE
if(Import.manual.spatial.domains){
  require(SPATA2) # installation https://themilolab.github.io/SPATA2/articles/spata-v2-installation.html
  
  st$segmentation = NA
  #sdomain = read.csv('/groups/tanaka/Collaborations/Jingkui-Elad/Mouse_Visium_annotations/Anno_166906.csv')
  inputDir = '/groups/tanaka/Collaborations/Jingkui-Elad/Mouse_Visium_annotations/'
  
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
  
  st$segmentation = gsub('Mixed Inj_BZ', 'Mixed_Inj_BZ', st$segmentation)
  st$segmentation = gsub('Proximal BZ', 'Proximal_BZ', st$segmentation)
  st$segmentation = gsub('Remote_3_septum', 'Remote_3', st$segmentation)
  
  st$segmentation = as.factor(st$segmentation)
  Idents(st) = as.factor(st$segmentation)
  SpatialDimPlot(st)
  
  ggsave(paste0(resDir, '/Manual_segmentation_spata2_Elad_v2.pdf'), width = 16, height = 6)
  
  
  save(st, design, file = paste0(RdataDir,
                                 'seuratObject_mouse_adult_cell.gene.filtered_umap.clustered_manualSegmentation', 
                                 species, '.Rdata'))
  
  
}else{ ### run bayesSpace to systematic spatial domain searching
  source('functions_Visium.R')
  run_bayesSpace(st, outDir = paste0(resDir, '/bayesSpace/'))
  
}

##########################################
# Step 2): cell neighborhood analysis
##########################################
##########################################
Run_Neighborhood_Enrichment_Analysis = FALSE
if(Run_Neighborhood_Enrichment_Analysis){
  
  # load ST data with additional region annotations
  load(file = paste0(RdataDir,
                     'seuratObject_mouse_adult_cell.gene.filtered_umap.clustered_manualSegmentation', 
                     species, '.Rdata'))
  
  st$segmentation = gsub('Distal_BZ', 'BZ', st$segmentation)
  st$segmentation = gsub('Proximal_BZ', 'BZ', st$segmentation)
  st$segmentation = gsub('Remote_1|Remote_2|Remote_3|Remote_3_septum', 'Remote', st$segmentation)
  st$segmentation = as.factor(st$segmentation)
  
  table(st$segmentation, st$condition)
  
  Idents(st) = as.factor(st$segmentation)
  SpatialDimPlot(st)
  
  ## snRNA-seq reference  
  refs = readRDS(file = paste0('../data/data_examples/ref_scRNAseq_adultMice_clean.v1.rds'))
  
  jj = which(refs$dataset == 'Ren2020')
  refs$timepoints[jj] = refs$condition[jj]
  refs$condition = refs$timepoints
  refs$subtypes = refs$subtype
  
  table(refs$subtypes)
  length(table(refs$subtypes))
  
  refs$subtypes = as.factor(refs$subtypes) 
  
  RCTD_out = paste0(resDir, '/celltype_deconvolution/', 
                    'RCTD_32Subtype_ref_v1.2')
  
  out_misty = paste0(resDir, '/neighborhood_test/Run_misty_1.0/')
  
  levels(refs$subtypes)
  
  source('functions_Visium.R')
  run_misty_colocalization_analysis(st, 
                                    outDir = out_misty,
                                    RCTD_out = RCTD_out,
                                    condSpec_celltypes = NULL
  )
  
  
}

########################################################
# step 3) ligand-receptor prediction analysis with 
# LIANA and NicheNet
# time-specifc and space-specific niches for nichenet
########################################################
# set parameter for ligand-receptor analysis
misty_Dir = paste0(out_misty, 'Plots_RCTD_density/')
misty_cutoff = 1

outDir_liana = paste0(resDir, '/LR_pairs/LIANA_v0.5_mistyCutoff_', misty_cutoff)

if(!dir.exists(outDir_liana)) system(paste0('mkdir -p ', outDir_liana))

# run LIANA day by day
timepoint_specific = TRUE

times_slice = levels(st$condition)
#times_slice = c('d1', 'd7', 'd14')

refs$celltypes = gsub('_', '.', refs$subtype)

subtypes = unique(refs$celltypes)

for(n in 1:length(times_slice))
{
  # n = 1
  source('functions_cccInference.R')
  time = times_slice[n]
  cat(' run LIANA for time -- ', time, '\n')
  
  outDir = paste(outDir_liana, '/', time, collapse = '')
  outDir = gsub(' ', '', outDir)
  
  ## select the interacting subtype pairs  
  intra = read.csv2(paste0(misty_Dir, time, '_all_summary_table_intra.csv'), row.names = 1)
  juxta = read.csv2(paste0(misty_Dir, time, '_all_summary_table_juxta5.csv'), row.names = 1)
  
  intra = intra > misty_cutoff
  juxta = juxta > misty_cutoff
  
  intra[which(is.na(intra))] = FALSE
  juxta[which(is.na(juxta))] = FALSE
  
  pairs = intra + juxta > 0
  #pairs = pairs[which(rownames(pairs) != 'RBC'), which(colnames(pairs) != 'RBC')]
  ss_row = apply(pairs, 1, sum)
  ss_col = apply(pairs, 2, sum)
  
  pairs = pairs[ ,which(ss_col >= 1)] # at least interacting with 1 receivers
  pairs = pairs[which(ss_row >= 3), ] # at least have 3 senders 
  
  #colnames(pairs) = gsub("Cav3_1", "Cav3.1", gsub('Mo_Macs', 'Mo.Macs', gsub('[.]','_', colnames(pairs))))
  #rownames(pairs) = gsub("Cav3_1", "Cav3.1", gsub('Mo_Macs', 'Mo.Macs', gsub('[.]','_', rownames(pairs))))
  
  cat(match(colnames(pairs), subtypes), '\n')
  cat(match(rownames(pairs), subtypes), '\n')
  
  celltypes_BZ_timeSpecific = vector("list", nrow(pairs))
  for(m in 1:nrow(pairs))
  {
    # m = 16
    celltypes_BZ_timeSpecific[[m]] = colnames(pairs)[which(pairs[m, ] == TRUE)]
    #x = pairs[which(rownames(pairs) == 'CM.Prol.IS'), ]
    names(celltypes_BZ_timeSpecific)[m] = rownames(pairs)[m]
    
  }
  
  run_LIANA(refs, 
            timepoint_specific = TRUE,
            include_autocrine = TRUE,
            celltypes_timeSpecific = celltypes_BZ_timeSpecific,
            outDir = outDir,
            species = 'mm'
  )
  
  #res = aggregate_output_LIANA(paste(outDir, time))
  
}


## double check the ligand and receptor expression distribution
#FeaturePlot(refs, features = rownames(refs)[grep('EGFC|VIPR2', rownames(refs))])
