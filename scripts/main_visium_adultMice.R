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
design = data.frame(seq(166907, 166904), c(paste0("adult.day", c(1, 4, 7, 14))), stringsAsFactors = FALSE)

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

########################################################
########################################################
# Section II : cell type deconvolution for each time point 
# 
########################################################
########################################################
source('functions_Visium.R')

load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, 
                   '_umap.clustered.Rdata'))
refs = readRDS(file = paste0('../data/data_examples/ref_scRNAseq.rds'))

jj = which(refs$dataset == 'Ren2020')
refs$timepoints[jj] = refs$condition[jj]
refs$condition = refs$timepoints

p1 = DimPlot(refs, reduction = 'umap', group.by = 'dataset')
p2 = DimPlot(refs, reduction = 'umap', group.by = 'celltype')
p3 = DimPlot(refs, reduction = 'umap', group.by = 'subtype')
p4 = DimPlot(refs, reduction = 'umap', group.by = 'timepoints')

(p1 + p2)/(p3 + p4)

ggsave(filename = paste0(resDir, '/UMAP_scRNAseq_refrence_dataset_timepoints_celltypes.pdf'), 
       width = 32, height = 16)

##########################################
# cell type deconvolution for cell types
##########################################
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


##########################################
# cell type deconvolution for subtypes
##########################################
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
max_cores = 32

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


########################################################
########################################################
# Section III : proximicity enrichment analysis
# 
########################################################
########################################################
##########################################
# spatial domain searching and potential define remote regions and border zone 
##########################################
Run_BayesSpace_for_spatialDomain = FALSE
if(Run_BayesSpace_for_spatialDomain){
  obj.list <- SplitObject(st, split.by = "condition")
  # select day4
  aa = obj.list[[2]]
  
  saveRDS(aa, file = paste0('data_examples/st_visium.rds'))
  
  aa$sampleID = design$sampleID[which(design$condition == names(table(aa$condition)))]
  
  # import manually defined spatial domain by Elad
  sdomain = read.csv('/groups/tanaka/Collaborations/Jingkui-Elad/Mouse_Visium_annotations/Anno_166906.csv')
  sdomain = sdomain[which(sdomain$Anno_1 != ''), ]
  aa$spatial_domain_manual = NA
  
  cells = gsub('_2_1', '',  colnames(aa))
  aa$spatial_domain_manual[match(sdomain$Barcode, cells)] = sdomain$Anno_1
  
  aa = run_bayesSpace(aa)
  
}

##########################################
# cell proximity analysis 
##########################################
run_cell_proximity_analysis(aa)

########################################################
########################################################
# Section IV: ligand-receptor-target prediction 
# LIANA and NICHNET  
########################################################
########################################################
source('functions_Visium.R')
run_LIANA()

run_NicheNet()


########################################################
########################################################
# Section V: spatial organization of cell types and genes  
# 
########################################################
########################################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, '_umap.clustered.Rdata'))

source('functions_Visium.R')
st = Find.SpatialDE(st)

