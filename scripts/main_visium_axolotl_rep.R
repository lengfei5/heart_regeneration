##########################################################################
##########################################################################
# Project: Heart regeneration 
# Script purpose: analyze axolot visium data
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Jan 18 10:16:57 2022
##########################################################################
##########################################################################
rm(list = ls())

species = 'axolotl'
version.analysis = '_R17246_R12830_allVisium_20240905'
dataDir = '../R17246_visium_axolotl/nf_out'

resDir = paste0("../results/visium_axolotl", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

source('functions_Visium.R')
library(pryr) # monitor the memory usage
require(ggplot2)
options(future.globals.maxSize = 120000 * 1024^2)

mem_used()

########################################################
########################################################
# Section I: import the processed visium data by nf from Tomas and 
# first processing and QCs 
########################################################
########################################################
design = data.frame(sampleID = c(seq(294946, 294949),seq(183623, 183626)), 
                    time = c(paste0('Amex_d', c(0, 4, 7, 0, 1, 4, 7, 14))), 
                    stringsAsFactors = FALSE)

design$condition = paste0(design$time, '_', design$sampleID)

varibleGenes = c()
check.QC.each.condition = TRUE

for(n in 1:nrow(design))
{
  # n = 1
  cat('-----------', design$condition[n], '--', design$sampleID[n], '-------------\n')
  # load nf output and process
  source('functions_Visium.R')
  aa = make_SeuratObj_visium(topdir = paste0(dataDir, '/', design$condition[n],  '/'), 
                             saveDir = paste0(resDir, '/', 
                                              design$condition[n], '/'), 
                             keyname = design$condition[n], 
                             QC.umi = TRUE)
  
  aa$condition = design$condition[n]
  aa$sampleid = design$sampleID[n]
  #aa <- SCTransform(aa, assay = "Spatial",  method = "glmGamPoi", verbose = FALSE)
  
  aa = subset(aa, subset = nCount_Spatial > 10) # 10 umi from the umi rank
  aa <- SCTransform(aa, assay = "Spatial", verbose = FALSE, 
                    variable.features.n = 3000, return.only.var.genes = FALSE)
  
  if(check.QC.each.condition){
    
    pdfname = paste0(resDir, '/QCs_gene_marker_check_', design$condition[n], version.analysis,  '.pdf')
    pdf(pdfname, width=16, height = 8)
    
    # Cell QC metrics: percentage of Mt, nb of counts, nb of genes 
    # get MT% (genes curated from NCBI chrMT genes)
    mtgenes = c("COX1", "COX2", "COX3", "ATP6", "ND1", "ND5", "CYTB", "ND2", 
                "ND4", "ATP8", "MT-CO1", "COI", "LOC9829747")
    
    mtgenes = c(mtgenes, paste0("MT", mtgenes), paste0("MT-", mtgenes))
    
    ggs = sapply(rownames(aa), function(x) unlist(strsplit(as.character(x), '-'))[1])
    mtgenes = rownames(aa)[!is.na(match(ggs, mtgenes))]
    # mtgenes = mtgenes[mtgenes %in% g[,1]]
    # srat = PercentageFeatureSet(srat, col.name = "percent.mt", assay = "Spatial",
    #                             features = mtgenes)
    xx = PercentageFeatureSet(aa, col.name = "percent.mt", assay = "SCT", features = mtgenes)
    aa[['percent.mt']] = xx$percent.mt
    
    Idents(aa) = design$condition[n]
    
    p1 = VlnPlot(aa, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), 
                 ncol = 3, pt.size = 1.0)
    plot(p1)
    
    p1 = FeatureScatter(aa, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
    p2 = FeatureScatter(aa, feature1 = "nCount_Spatial", feature2 = "percent.mt")
    p3 = FeatureScatter(aa, feature1 = "nFeature_Spatial", feature2 = "percent.mt")
    plot(wrap_plots(p1, p2, p3))
    
    plot1 <- SpatialFeaturePlot(aa, features = "nCount_Spatial") + theme(legend.position = "bottom")
    plot2 <- SpatialFeaturePlot(aa, features = "nFeature_Spatial") + theme(legend.position = "bottom")
    plot3 <- SpatialFeaturePlot(aa, features = "percent.mt") + theme(legend.position = "bottom")
    plot(wrap_plots(plot1, plot2, plot3))
    
    #if(design$species[n] == 'neonatal'){ pct.mt.cutoff = 0.6; }
    #if(design$species[n] == 'adult') { pct.mt.cutoff = 1 }
    #aa = aa[, which(aa$percent.mt < pct.mt.cutoff)]
    #cat(ncol(aa), ' spots left after spot filtering \n')
    
    aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
    ElbowPlot(aa)
    
    aa <- FindNeighbors(aa, dims = 1:10)
    aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.7)
    aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 30, min.dist = 0.05)
    
    DimPlot(aa, reduction = "umap", group.by = c("ident"))
    
    features = rownames(aa)[grep('MYH6|NPPA|CLU-AMEX60DD032706', rownames(aa))]
    #FeaturePlot(aa, features = features)
    
    SpatialFeaturePlot(aa, features = features[2])
    
    dev.off()
    
    
    saveRDS(aa, file = paste0(RdataDir, 'seuratObject_design_st_', design$condition[n], 
                               version.analysis, '.rds'))
    
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

save(design, varibleGenes, st, 
     file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, version.analysis, 
                   '.Rdata'))

########################################################
########################################################
# Section II: QC and marker gene checking
# 
########################################################
########################################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, 
                   version.analysis, '.Rdata'))

st$condition = factor(st$condition, levels = c('Amex_d0_294946', 'Amex_d0_294949',
                                                  'Amex_d1_183623',
                                                  'Amex_d4_294947', 'Amex_d4_183624',
                                                  'Amex_d7_294948', 'Amex_d7_183625', 
                                                  'Amex_d14_183626'))

Idents(st) = st$condition

SpatialFeaturePlot(st, features = 'nCount_Spatial', ncol = 4, image.alpha = 0.1)
SpatialFeaturePlot(st, features = 'nFeature_Spatial', ncol = 4)

##########################################
# gene and cell filtering (to add later)
##########################################
Filtering.cells.genes = FALSE
if(Filtering.cells.genes){
  #st[["percent.mt"]] <- PercentageFeatureSet(st, pattern = "^Mt-")
  
  # Visualize QC metrics as a violin plot
  #Idents(st) = st$condition
  
  VlnPlot(st, features = c("nCount_Spatial"), log = FALSE) + 
    geom_hline(yintercept = c(300, 500)) + 
    scale_y_continuous(limits = c(0, 5000))
  
  ggsave(paste0(resDir, '/QCs_reseq_nUMI.counts_linear.pdf'), width = 12, height = 8)
  
  VlnPlot(st, features = c("nFeature_Spatial")) + 
    geom_hline( yintercept = c(100, 200, 500)) + 
    scale_y_continuous(limits = c(0, 2000))
  ggsave(paste0(resDir, '/QCs_reseq_nFeatures_Spatial.pdf'), width = 12, height = 8)
  
  FeatureScatter(st, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") + 
  geom_hline(yintercept = c(100, 200)) + 
  geom_vline(xintercept = c(100, 300, 500))
  
  FeatureScatter(st, feature1 = "nCount_Spatial", feature2 = "nCount_SCT")
  
  VlnPlot(st, features = c("nCount_SCT"))  
    #geom_hline( yintercept = c(2500, 5000)) + 
    #scale_y_continuous(limits = c(0, 30000))
  ggsave(paste0(resDir, '/QCs_reseq_nCounts_SCT.pdf'), width = 12, height = 8)
  
  VlnPlot(st, features = c("nFeature_SCT"))  
  #geom_hline( yintercept = c(2500, 5000)) + 
  #scale_y_continuous(limits = c(0, 30000))
  ggsave(paste0(resDir, '/QCs_reseq_nCounts_SCT.pdf'), width = 12, height = 8)
  
  #plot1 <- VlnPlot(st, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  #plot2 <- SpatialFeaturePlot(aa, features = "nCount_Spatial") + theme(legend.position = "right")
  #wrap_plots(plot1, plot2)
  
  features = rownames(st)[grep('MYH6|NPPA|CLU-AMEX60DD032706', rownames(st))]
  #FeaturePlot(aa, features = features)
  
  SpatialFeaturePlot(st, features = features[2], ncol = 4)
  ggsave(paste0(resDir, '/QCs_Visium_featurePlots_NPPA.pdf'), width = 12, height = 8)
  
  SpatialFeaturePlot(st, features = features[4], ncol = 4)
  ggsave(paste0(resDir, '/QCs_Visium_featurePlots_MYH6.pdf'), width = 12, height = 8)
  
  
  
  p1 = VlnPlot(st, features = c("nCount_Spatial"), pt.size = 0.7) + NoLegend() +
    geom_hline(yintercept = 200)
  
  p2 = VlnPlot(st, features = c("nFeature_Spatial"), pt.size = 0.7) + NoLegend() +
    geom_hline(yintercept = 100)
  
  plot(wrap_plots(p1, p2))
  
  
  cat(ncol(st), ' spots before cell filtering \n')
  sels = which(st$nCount_Spatial > 300 & st$nFeature_Spatial > 200)
  st= st[, sels]
  cat(ncol(st), ' spots after cell filtering \n')
  
  p1 = SpatialFeaturePlot(st, features = features[2])
  
  p2 = SpatialFeaturePlot(st, features = features[4])
  
  p1 / p2
 
  saveRDS(st, file = paste0(RdataDir, 'seuratObject_st_filtered.spots', version.analysis, '.rds'))
  
}


st = readRDS(file = paste0(RdataDir, 'seuratObject_st_filtered.spots', version.analysis, '.rds'))

##########################################
# test SCT normalization
##########################################
#st = SCTransform(st, assay = "Spatial", verbose = FALSE)
DefaultAssay(st) <- "SCT"
VariableFeatures(st) <- varibleGenes

st <- RunPCA(st, verbose = FALSE)
ElbowPlot(st)

st <- FindNeighbors(st, dims = 1:10)
st <- FindClusters(st, verbose = FALSE, resolution = 1.0)

st <- RunUMAP(st, dims = 1:20, n.neighbors = 20, min.dist = 0.05)
SpatialDimPlot(st, image.alpha = 0.5)

DimPlot(st, reduction = "umap", group.by = c("ident", "condition")) 

ggsave(filename = paste0(resDir, '/UMAP_all.timepoints_', species, '.pdf'), width = 16, height = 8)

#features = rownames(st)[grep('MYH6|NPPA|AXL|CLU-AMEX60DD032706', rownames(st))]
#FeaturePlot(st, features = features)
#SpatialFeaturePlot(st, features = features[1])

save(design, varibleGenes, st, 
     file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered', 
                   species, version.analysis,'.Rdata'))

##########################################
# QCs with marker genes
##########################################
st = readRDS(file = paste0(RdataDir, 'seuratObject_st_filtered.spots', version.analysis, '.rds'))

SpatialFeaturePlot(st, features = 'nCount_Spatial', image.alpha = 0.5, ncol = 4)

ggsave(paste0(resDir, '/QCs_allVisium_SpatialFeatures_nCount.pdf'), width = 12, height = 8)

SpatialFeaturePlot(st, features = 'nFeature_Spatial', image.alpha = 0.5, ncol = 4)
ggsave(paste0(resDir, '/QCs_allVisium_SpatialFeatures_nFeature.pdf'), width = 12, height = 8)

VlnPlot(st, features = c("nCount_Spatial"), log = FALSE) + 
  geom_hline(yintercept = c(300, 500)) + 
  scale_y_continuous(limits = c(0, 20000))

ggsave(paste0(resDir, '/QCs_allVisium_VlinPlot_nCount.pdf'), width = 12, height = 8)

VlnPlot(st, features = c("nFeature_Spatial"), log = FALSE) + 
  geom_hline(yintercept = c(500, 1000)) + 
  scale_y_continuous(limits = c(0, 3000))
ggsave(paste0(resDir, '/QCs_allVisium_VlinPlot_nFeature.pdf'), width = 12, height = 8)

VlnPlot(st, features = c("percent.mt"), log = FALSE) + 
  geom_hline(yintercept = c(30, 50)) + 
  scale_y_continuous(limits = c(0, 70))
ggsave(paste0(resDir, '/QCs_allVisium_VlinPlot_pctMT.pdf'), width = 12, height = 8)


QCs.with.marker.genes = FALSE
if(QCs.with.marker.genes){
  # import marker genes
  library(openxlsx)
  aa = read.xlsx('../data/Markers_updated_v2.xlsx', sheet = 1, colNames = TRUE)
  #aa = read.csv('../data/Tzahor_geneList.csv', header = TRUE)
  
  #aa = aa[-c(1), c(1:2)]
  
  markers = c()
  for(n in 1:ncol(aa))
  {
    markers = c(markers, aa[!is.na(aa[,n]), n])
  }
  markers = unique(c(markers, c('hmgb2', 'top2a', 'AGRN')))
  
  markers = unique(c(markers, 
                     c("Ntn1", "Mfap5", "Ubb", "Spon2", "Sparc", "Comp", 
                       "Ncanb", "Cthrc1", "Gpr42", "Ffr2")))
                   
  markers = markers[which(markers != '')]
  markers = firstup(tolower(unique(markers)))
  markers = gsub(' ', '', markers)
  
  xx = c()
  for(n in 1:length(markers))
  {
    xx = c(xx, rownames(st)[grep(paste0('^', toupper(markers[n])), rownames(st))])
  }
  
  xx = unique(xx)
  markers = xx
  
  #xx = rownames(st)[grep('AGRN', rownames(st))]
  #SpatialFeaturePlot(st, features = xx, image.alpha = 0.5)
  
  #markers[is.na(match(markers, rownames(st)))]
  #markers = markers[!is.na(match(markers, rownames(st)))]
  xx = rownames(st)[grep('CLU', rownames(st))]
  p0 = SpatialFeaturePlot(st, features = features[4], image.alpha = 0.2)
  p1 = SpatialDimPlot(st, image.alpha = 0.5)
  plot(wrap_plots(p0, p1, nrow = 2))
  
  pdfname = paste0(resDir, "/check_detected_celltypes_using_AdditionalMarkerGenes_Bern_", 
                   species, ".pdf")
  pdf(pdfname, width = 16, height = 8)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  for(n in 1:length(markers))
  {
    cat(markers[n], '\n')
    p3 = SpatialFeaturePlot(st, features = markers[n], image.alpha = 0.5, ncol = 4)
    
    plot(p3)
    
  }
  
  dev.off()
  
}


########################################################
########################################################
# Section III: start the analysis 
# analyze the axolotl visium per condition
# To identify the boarder zone
########################################################
########################################################
Data_Exploration = FALSE
if(Data_Exploration){
  ##########################################
  # loop over all conditions of st
  ##########################################
  #st = readRDS(file = paste0(RdataDir, 'seuratObject_st_filtered.spots', version.analysis, '.rds'))
  st = readRDS(file = paste0(RdataDir, 'seuratObject_st_filtered.spots_time_conditions', 
                             version.analysis, '.rds'))
  #load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered', species, '.Rdata'))
  #st = Seurat::SplitObject(st, split.by = 'condition')
  
  #st$condition = factor(st$condition, levels = design$condition)
  
  DefaultAssay(st) = 'Spatial'
  st <- NormalizeData(st, normalization.method = "LogNormalize", scale.factor = 10000)
  st <- FindVariableFeatures(st, selection.method = "vst", nfeatures = 3000)
  all.genes <- rownames(st)
  st <- ScaleData(st, features = all.genes)
  
  st <- RunPCA(st, verbose = FALSE, features = VariableFeatures(object = st), weight.by.var = TRUE)
  ElbowPlot(st, ndims = 30)
  
  st <- RunUMAP(st, dims = 1:20, n.neighbors = 30, min.dist = 0.05)
  
  DimPlot(st, group.by = 'condition', label = TRUE, repel = TRUE)
  
  cat('visium conditions :\n')
  print(table(design$condition))
  cc = design$condition
  
  use.SCTransform = TRUE
  
  for(n in 1:length(cc))
    #for(n in 1:2)
  {
    # n = 4
    aa = readRDS(file = paste0(RdataDir, 'seuratObject_design_st_', cc[n],  '.rds'))
    
    if(use.SCTransform){
      DefaultAssay(aa) <- "SCT"
    }else{
      
      DefaultAssay(aa) = 'Spatial'
      aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
      aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000)
      all.genes <- rownames(aa)
      aa <- ScaleData(aa, features = all.genes)
      
    }
    
    aa <- RunPCA(aa, verbose = FALSE, features = VariableFeatures(object = aa), weight.by.var = TRUE)
    ElbowPlot(aa)
    
    DimPlot(aa, reduction = 'pca')
    
    par(mfrow = c(2,2))
    plot(aa@reductions$pca@cell.embeddings[,1], aa$nCount_Spatial, xlab = 'PC1')
    plot(aa@reductions$pca@cell.embeddings[,2], aa$nCount_Spatial, xlab = 'PC2')
    plot(aa@reductions$pca@cell.embeddings[,3], aa$nCount_Spatial, xlab = 'PC3')
    plot(aa@reductions$pca@cell.embeddings[,4], aa$nCount_Spatial, xlab = 'PC4')
    
    aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.05)
    
    aa <- FindNeighbors(aa, dims = 1:20)
    aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 1.5)
    
    p1 = DimPlot(aa, reduction = "umap", group.by = c("ident"), label = TRUE, label.size = 8)
    p2 <- SpatialDimPlot(aa, label = TRUE, label.size = 5)
    p1 + p2
    
    
    ggsave(filename =  paste0(resDir, "/Spatial_patterning_", species, '_', cc[n], ".pdf"), 
           width = 12, height = 8)
    
    features = rownames(st)[grep('MYH6|NPPA|CLU-AMEX60DD032706', rownames(st))]
    FeaturePlot(st, features = features)
    
    SpatialFeaturePlot(st, features = features[2])
    
    ##########################################
    # to check if there is border zone, manually select border zone and remote zones from images 
    ##########################################
    source('functions_Visium.R')
    aa = manual_selection_spots_image_Spata(aa, slice = cc[n])
    
    SpatialDimPlot(aa, group.by = 'segmentation', label = TRUE, label.size = 5)
    
    ggsave(filename =  paste0(resDir, "/manual_segmentation_", species, '_', cc[n], ".pdf"), 
           width = 12, height = 8)
    
    Idents(aa) = aa$segmentation
    
    border_markers1 = FindMarkers(aa, ident.1 = 'border_zone', ident.2 = 'remote_zone1', only.pos = FALSE, 
                                  min.pct = 0.2, logfc.threshold = 0.25)
    
    border_markers2 = FindMarkers(aa, ident.1 = 'border_zone', ident.2 = 'remote_zone2', only.pos = FALSE, 
                                  min.pct = 0.2, logfc.threshold = 0.25)
    
    SpatialFeaturePlot(aa, features = rownames(border_markers1)[c(1:4)])
    
    SpatialFeaturePlot(aa, features = rownames(border_markers2)[c(1:4)])
    
    
    ##########################################
    # test BayesSpace to predict spatial domain
    ##########################################
    source('functions_Visium.R')
    
    ##########################################
    # # test SpatialDE
    ##########################################
    source('functions_Visium.R')
    ggs = Find.SpatialDE(aa, use.method = 'sparkX')
    
    markers = rownames(ggs)[which(ggs$adjustedPval<10^-6)]
    markers = markers[c(1:300)]
    
    pdfname = paste0(resDir, "/Spatial_patterningGenes_", species, '_', cc[n], ".pdf")
    pdf(pdfname, width = 20, height = 8)
    par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
    
    for(n in 1:length(markers))
    {
      cat(n, '--', markers[n], '\n')
      p3 = SpatialFeaturePlot(st, features = markers[n], alpha = c(0.1, 1))
      
      plot(p3)
      
    }
    
    dev.off()
    
    
    ##########################################
    #  # test different clustering methods
    ##########################################
    source('functions_Visium.R')
    aa = findClusters_SC3(aa)
    
    ggsave(filename = paste0(resDir, '/Visium_Clustering_SptialDimPLot', cc[n], '.pdf'), 
           width = 16, height = 8)
    
    # find marker of those clusters
    aa.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.3)
    
    aa.markers %>%
      group_by(cluster) %>%
      top_n(n = 15, wt = avg_log2FC) -> top10
    
    DoHeatmap(aa, features = top10$gene)
    
    ggsave(paste0(resDir, '/heatmap_markerGenes_', cc[n], '.pdf'), width = 12, height = 20)
    
    
  }
  
}

########################################################
########################################################
# Section IV : cell type deconvolution for each time point with RCTD (not cell2location)
# using the snRNA-seq annotated by Elad
########################################################
########################################################
st = readRDS(file = paste0(RdataDir, 'seuratObject_st_filtered.spots_time_conditions', 
                               version.analysis, '.rds'))
#load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered', species, '.Rdata'))
#st = readRDS(file = paste0(RdataDir, 'seuratObject_st_filtered.spots', version.analysis, '.rds'))
#st$time = design$time[match(st$sampleid, design$sampleID)]

#saveRDS(st, file = paste0(RdataDir, 'seuratObject_st_filtered.spots_time_conditions', 
#                          version.analysis, '.rds'))

#st$condition = factor(st$condition, levels = design$condition)

table(st$condition)

# refined subtypes by Elad
refs_file = paste0('/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_subtypes_final_20221117.rds')

refs = readRDS(file = refs_file)
table(refs$subtypes)
length(table(refs$subtypes))

refs$subtypes = droplevels(refs$subtypes) 
length(table(refs$subtypes)) 
## only 41 subtype annotations from Elad, with additional annotation "doublet" with 0 cell

## prepare the celltype to use and also specify the time-specific subtypes
refs$celltype_toUse = as.character(refs$subtypes)
length(table(refs$celltype_toUse))

refs$condition = gsub('_scRNA', '', refs$condition)
refs$celltype_toUse = gsub('Mo/Macs', 'Mo.Macs', refs$celltype_toUse)
refs$celltype_toUse = gsub("[(]", '', refs$celltype_toUse)
refs$celltype_toUse = gsub("[)]", '', refs$celltype_toUse)

table(refs$celltype_toUse)
length(table(refs$celltype_toUse))



Import_manually_timeSpecific_celltypes = FALSE
if(Import_manually_timeSpecific_celltypes){
  condition.specific_celltypes = openxlsx::read.xlsx('../data/subtypes_timespecific_v2.xlsx', 
                                                     rowNames = TRUE)
  colnames(condition.specific_celltypes) = gsub('day', 'd', 
                                                paste0('Amex_', colnames(condition.specific_celltypes)))
  
  names = gsub('MO/Macs', 'Mo.Macs', rownames(condition.specific_celltypes))
  names = gsub("[(]", '', names)
  names = gsub("[)]", '', names)
  names = gsub('_PROL', '_Prol', names)
  names = gsub('_prol', '_Prol', names)
  names = gsub("CM_ROBO2", "CM_ven_Cav3_1", names)
  names = gsub("CM_CAV3.1", "CM_ven_Robo2", names)
  names = gsub("CM_IS_Prol", "CM_Prol_IS", names)
  names = gsub("CM_Prol1", "CM_Prol_1", names)
  names = gsub("CM_Prol2", "CM_Prol_2", names)
  names = gsub("CM_Prol3", "CM_Prol_3", names)
  names = gsub("CM_atria$", "CM_Atria", names)
  names = gsub("CM_atria_tagln", "CM_Atria_Tagln", names)
  names = gsub("RBC_Prol", "Proliferating_RBC", names)
  names = gsub("Megakeryocytes_Prol", "Proliferating_Megakeryocytes", names)
  
  rownames(condition.specific_celltypes) = names
  
}else{
  cutoff = 0.25 # fraciton of subtypes in the time point
  csc = table(refs$celltype_toUse, refs$condition)
  #colnames(csc) = gsub('Amex_', '', colnames(csc))
  csc = csc[, c(1,2, 4, 5, 3)]
  for(n in 1:ncol(csc))
  {
    csc[,n] = csc[,n]/sum(csc[,n])*100
  }
  
  csc[which(csc<=0.25)] = NA
  #csc[which(!is.na(csc))] = 
  
  ## further removing the proliferating cells
  csc[which(rownames(csc) == 'Proliferating_RBC'), ] = NA
  csc[which(rownames(csc) == 'Proliferating_Megakeryocytes'), ] = NA
  csc[which(rownames(csc) == 'Neuronal'), ] = NA
  csc[which(rownames(csc) == 'Mo.Macs_Prol'), ] = NA
  
  csc[which(rownames(csc) == 'EC_Prol'), ] = NA
  #csc[which(rownames(csc) == 'EC_IS_Prol'), ] = NA
  
  csc[which(rownames(csc) == 'CM_Prol_3'), ] = NA
  csc[which(rownames(csc) == 'CM_Prol_2'), ] = NA
  csc[which(rownames(csc) == 'CM_Prol_1'), ] = NA
  
  csc[which(rownames(csc) == 'B_cells_Prol'), ] = NA
  
  condition.specific_celltypes = csc
  
  
}


celltypes = unique(refs$celltype_toUse)
mm = match(rownames(condition.specific_celltypes), celltypes)
cat(length(which(is.na(mm))), ' cell types missing in the condition.specific_celltypes \n')


# st = subset(st, condition == 'Amex_d4')
saveRDS(st, file = paste(resDir, 'RCTD_st_axolotl_allVisium.rds'))

saveRDS(refs, file = paste0(resDir, '/RCTD_refs_subtypes_final_20221117.rds'))

saveRDS(condition.specific_celltypes, 
        file = paste0(resDir, '/RCTD_refs_condition_specificity_v3.rds'))


source('functions_Visium.R')
## preapre the paramters for RCTD subtypes
DefaultAssay(refs) = 'RNA'
DefaultAssay(st) = 'Spatial'
require_int_SpatialRNA = FALSE

condition.specific.ref = TRUE

RCTD_out = paste0(resDir, '/RCTD_subtype_out_41subtypes_ref.time.specific_v2.0')
max_cores = 16

cc = names(table(st$condition))


source('functions_Visium.R')
Run.celltype.deconvolution.RCTD(st = st,
                                refs, 
                                mapping_ref = 'time',
                                condition.specific.ref = condition.specific.ref,
                                condition.specific_celltypes = condition.specific_celltypes,
                                require_int_SpatialRNA = require_int_SpatialRNA,
                                max_cores = max_cores,
                                RCTD_out = RCTD_out,
                                plot.RCTD.summary = FALSE, 
                                PLOT.scatterpie = FALSE
                                
)


source('functions_Visium.R')
plot.RCTD.results(st = st,
                  RCTD_out = RCTD_out,
                  plot.RCTD.summary = FALSE)


## prepare the parameters for RCTD coarse cell types
Run.RCTD.coarse.celltypes = FALSE
if(Run.RCTD.coarse.celltypes){
  # define coarse clusters
  refs$celltypes = as.character(refs$subtypes)
  
  refs$celltypes[grep('CM_|CMs_|_CM|_CM_', refs$subtypes)] = 'CM'
  refs$celltypes[grep('EC_|_EC', refs$subtypes)] = 'EC'
  refs$celltypes[grep('FB_', refs$subtypes)] = 'FB'
  refs$celltypes[grep('B_cells', refs$subtypes)] = 'Bcell'
  
  refs$celltypes[grep('Macrophages|_MF', refs$subtypes)] = 'Macrophages'
  refs$celltypes[grep('Megakeryocytes', refs$subtypes)] = 'Megakeryocytes'
  refs$celltypes[grep('RBC', refs$subtypes)] = 'RBC'
  
  refs$celltype_toUse = refs$celltypes
  DefaultAssay(refs) = 'RNA'
  DefaultAssay(st) = 'Spatial'
  require_int_SpatialRNA = FALSE
  RCTD_out = paste0(resDir, '/RCTD_coarse_out_v1')
  
  max_cores = 16
  
  source('functions_Visium.R')
  
  Run.celltype.deconvolution.RCTD(st, refs, 
                                  require_int_SpatialRNA = require_int_SpatialRNA,
                                  max_cores = max_cores,
                                  RCTD_out = RCTD_out
  )
  
}

########################################################
########################################################
# Section V: manually define borzder zones, injury zones

########################################################
########################################################
source('functions_Visium.R')
st = readRDS(file = paste0(RdataDir, 'seuratObject_st_filtered.spots_time_conditions', 
                           version.analysis, '.rds'))


Idents(st) = st$condition
cat('visium conditions :\n')
print(table(st$condition))
cc = names(table(st$condition))


##########################################
# step 1) Spatial domain searching and potential define remote regions and border zone
# here using computational methods to define regions of interest or cell niches
##########################################
## import manually defined spatial domains
Import.manual.spatial.domains = FALSE
if(Import.manual.spatial.domains){
  require(SPATA2) # installation https://themilolab.github.io/SPATA2/articles/spata-v2-installation.html
  
  st$segmentation = NA
  
  #inputDir = '/groups/tanaka/Collaborations/Jingkui-Elad/visium_axolotl_reseq/spata2_manual_regions/'
  #saveRDS(aa, file = paste0(inputDir, 'axolotl_visium_newRep.rds'))
  inputDir = paste0(resDir, '/manual_borderZones_spata2/')
  
  #manual_selection_spots_image_Spata
  for(n in 1:length(cc))
  {
    # n = 1
    cat('slice -- ', cc[n], '\n')
    slice = cc[n]
    
    file = paste0(inputDir, 'segemented_spata2_', slice, '.rds')
    
    if(file.exists(file)){
      aa = readRDS(file = file)
      plotSegmentation(object = aa, pt_size = 1.9) +
        ggplot2::scale_y_reverse()
      
      segs = getSegmentNames(aa)
      for(s in segs)
      {
        cat('-- segmentation : ', s, '--\n')
        st$segmentation[match(getSegmentDf(aa, segment_names = c(s))$barcodes, 
                              colnames(st))] = s
      }
      
      #Idents(aa) = as.factor(aa$segmentation)
      
      #SpatialDimPlot(aa)
      #mm = match(colnames(aa), colnames(st))
      #jj = which(!is.na(mm))
      #st$segmentation[mm[jj]] = aa$segmentation[jj]
      
    }else{
      cat('no manual segmentation found for -- ', slice, '\n' )
    }
        
  }
  
  st$segmentation = as.factor(st$segmentation)
  Idents(st) = as.factor(st$segmentation)
  SpatialDimPlot(st, ncol = 4)
  
  # add the segmentation from the old replicates
  st$cell.id = colnames(st)
  st$cell.id = sapply(st$cell.id, function(x) unlist(strsplit(as.character(x), '[-]'))[1])
  st$cell.id = paste0(st$condition, '_', st$cell.id)
  
  aa = st;
  rm(st)
  
  load(file = paste0('../results/Rdata/', 
                     'seuratObject_design_variableGenes_umap.clustered_manualSegmentation', 
                     'axolotl', '.Rdata'))
  rm(design)
  st$condition = as.character(st$condition)
  st$condition[which(st$condition == 'Amex_d1')] = 'Amex_d1_183623'
  st$condition[which(st$condition == 'Amex_d4')] = 'Amex_d4_183624'
  st$condition[which(st$condition == 'Amex_d7')] = 'Amex_d7_183625'
  st$condition[which(st$condition == 'Amex_d14')] = 'Amex_d14_183626'
  st$cell.id = colnames(st)
  st$cell.id = sapply(st$cell.id, function(x) unlist(strsplit(as.character(x), '[-]'))[1])
  st$cell.id = paste0(st$condition, '_', st$cell.id)
  
  mm = match(st$cell.id, aa$cell.id)
  
  aa$segmentation = as.character(aa$segmentation)
  st$segmentation = as.character(st$segmentation)
  aa$segmentation[mm[which(!is.na(mm))]] = st$segmentation[which(!is.na(mm))]
  
  st = aa;
  st$segmentation = as.factor(st$segmentation)
  Idents(st) = as.factor(st$segmentation)
  SpatialDimPlot(st, ncol = 4)
  
  rm(aa)
  
  st$segmentation[which(st$segmentation == 'border_zone')] = 'BZ'
  st$segmentation[which(st$segmentation == 'Border_zone')] = 'BZ'
  
  st$segmentation[which(st$segmentation == 'injury_zone')] = 'Injury'
  st$segmentation[which(st$segmentation == 'InjUry_zone')] = 'Injury'
  
  st$segmentation[which(st$segmentation == 'Remote1')] = 'RZ1'
  st$segmentation[which(st$segmentation == 'Remote2')] = 'RZ2'
  
  st$segmentation = droplevels(st$segmentation)
  
  Idents(st) = as.factor(st$segmentation)
  SpatialDimPlot(st, ncol = 4)
  
  ggsave(paste0(resDir, '/Manual_segmentation_spata2_Elad.pdf'), width = 16, height = 10)
  
  saveRDS(st, file = paste0(RdataDir, 'seuratObject_allVisiusmst_',
                            'filtered.spots_time_conditions_manualSegmentation', 
                            version.analysis, '.rds'))
  
  
  
}

##########################################
# compare the injury zone, border zone, remote and intact   
##########################################
st = readRDS(file = paste0(RdataDir, 'seuratObject_allVisiusmst_',
                           'filtered.spots_time_conditions_manualSegmentation', 
                           version.analysis, '.rds'))

st = subset(st, cells = which(!is.na(st$segmentation)))

#st <- SCTransform(st, assay = "Spatial", verbose = FALSE, variable.features.n = 3000)
# DefaultAssay(st) <- "SCT"
st = NormalizeData(st, normalization.method = "LogNormalize", scale.factor = 10000)
st = FindVariableFeatures(st, selection.method = "vst", nfeatures = 2000)
st = ScaleData(st,  
               vars.to.regress = c('nCount_Spatial'), 
               assay = 'Spatial')

st <- RunPCA(st, verbose = FALSE)
ElbowPlot(st)

st <- FindNeighbors(st, dims = 1:20)
st <- FindClusters(st, verbose = FALSE, resolution = 1.0)

st <- RunUMAP(st, dims = 1:20, n.neighbors = 30, min.dist = 0.1)

st$segmentation = as.character(st$segmentation)
st$segmentation[which(st$segmentation == 'RZ1')] = 'RZ'
st$segmentation[which(st$segmentation == 'RZ2')] = 'RZ'
st$segmentation[which(st$segmentation == 'Intact1')] = 'Intact'
st$segmentation[which(st$segmentation == 'Intact2')] = 'Intact'

st$segmentation = as.factor(st$segmentation)

p1 = DimPlot(st, reduction = "umap", group.by = c("condition"), label = TRUE, repel = TRUE) 
p2 = DimPlot(st, reduction = "umap", group.by = c("segmentation"), label = TRUE, repel = TRUE) 

p1 / p2

ggsave(paste0(resDir, '/Compare_diffRegions_intact_injuried_time_logNormal_regressed.nCounts.pdf'), 
       width = 10, height = 12)



st$dataset = 'batch1'
st$dataset[grep('_2949', st$condition)] = 'batch2'

source('functions_dataIntegration.R')
st_bc = IntegrateData_Seurat_CCA(st, group.by = 'dataset', k.weight = 100)

st_bc = IntegrateData_runFastMNN(st, group.by = 'dataset', assays = 'Spatial',
                                 nfeatures = 1000,
                                 correct.all = FALSE)

st_bc <- RunUMAP(st_bc, reduction = "mnn", dims = 1:20, 
                       n.neighbors = 30, min.dist = 0.1)

st_bc$condition = factor(st_bc$condition, levels = c('Amex_d0_294946', 'Amex_d0_294949',
                                                     'Amex_d1_183623', 'Amex_d4_294947', 
                                                     'Amex_d4_183624', 'Amex_d7_294948', 
                                                     'Amex_d7_183625', 'Amex_d14_183626'))

p1 = DimPlot(st_bc, reduction = "umap", group.by = c("condition"), label = TRUE, repel = TRUE) 
p2 = DimPlot(st_bc, reduction = "umap", group.by = c("segmentation"), label = TRUE, repel = TRUE) 

p1 / p2

ggsave(paste0(resDir, '/Compare_diffRegions_intact_injuried_time_logNormal_seurat.RunfastMNN.pdf'), 
       width = 10, height = 12)
#st$time = gsub('Amex_', '', st$time)
#st$cc = paste0(st$time, '_', st$segmentation)
#p3 = DimPlot(st, reduction = "umap", group.by = c("cc"), label = TRUE, repel = TRUE) 

#p1 / p2 /p3
SpatialPlot(st, group.by = 'segmentation', image.alpha = 0.5, ncol = 4)

ggsave(filename = paste0(resDir, '/UMAP_all.timepoints_', species, '.pdf'), width = 16, height = 8)


