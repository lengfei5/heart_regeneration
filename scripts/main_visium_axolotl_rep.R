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
#st = readRDS(file = paste0(RdataDir, 'seuratObject_st_filtered.spots_time_conditions', 
#                               version.analysis, '.rds'))
st = readRDS(file = paste0(RdataDir, 'seuratObject_allVisiusmst_',
              'filtered.spots_time_conditions_manualSegmentation_ventricleRegions', 
              version.analysis, '.rds'))

Idents(st) = st$seg_ventricle

## subset only the ventricle regions 
st = subset(st, cells = colnames(st)[which(!is.na(st$seg_ventricle))])

SpatialPlot(st, group.by = 'seg_ventricle', ncol = 4)

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

##########################################
# double check the subtype similarity used by RCTD 
##########################################
Assess_subtypes_similarity_for_RCTD = FALSE
if(Assess_subtypes_similarity_for_RCTD){
  
  ## preapre the paramters for RCTD subtypes
  DefaultAssay(refs) = 'RNA'
  DefaultAssay(st) = 'Spatial'
  require_int_SpatialRNA = FALSE
  
  RCTD_out = paste0(resDir, '/RCTD_out', 
                    '/subtype_similarity/')
  system(paste0('mkdir -p ', RCTD_out))
  
  source('functions_Visium.R')
  assess_sutypes_similarity_for_RCTD(st, refs, 
                                  mapping_ref = 'time',
                                  require_int_SpatialRNA = require_int_SpatialRNA,
                                  RCTD_out = RCTD_out)
  
  
  
}

##########################################
# manually set time-specifici subtypes used for RCTD
##########################################
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
  
  csc[which(rownames(csc) == 'Neu_IL1R1'), ] = NA
  
  ## further removing the proliferating cells
  csc[which(rownames(csc) == 'Proliferating_RBC'), ] = NA
  csc[which(rownames(csc) == 'Proliferating_Megakeryocytes'), ] = NA
  csc[which(rownames(csc) == 'Neuronal'), ] = NA
  csc[which(rownames(csc) == 'Mo.Macs_Prol'), ] = NA
  csc[which(rownames(csc) == 'B_cells_Prol'), ] = NA
  
  csc[which(rownames(csc) == 'EC_Prol'), ] = NA
  csc[which(rownames(csc) == 'EC_IS_Prol'), ] = NA
  
  csc[which(rownames(csc) == 'CM_Prol_3'), ] = NA
  csc[which(rownames(csc) == 'CM_Prol_2'), ] = NA
  csc[which(rownames(csc) == 'CM_Prol_1'), ] = NA
  csc[which(rownames(csc) == 'CM_Prol_1'), ] = NA
  
  csc[which(rownames(csc) == 'FB_Prol'), ] = NA
  
  ## remove the subtypes in Atria
  csc[which(rownames(csc) == 'EC_WNT4'), ] = NA
  csc[which(rownames(csc) == 'EC_CEMIP'), ] = NA
  csc[which(rownames(csc) == 'EC_LHX6'), ] = NA
  csc[which(rownames(csc) == 'FB_VWA2'), ] = NA
  csc[which(rownames(csc) == 'FB_TNXB'), ] = NA
  csc[which(rownames(csc) == 'CM_PM_HCN4'), ] = NA
  csc[which(rownames(csc) == 'CM_OFT'), ] = NA
  csc[which(rownames(csc) == 'CM_Atria'), ] = NA
  csc[which(rownames(csc) == 'CM_Atria_Tagln'), ] = NA
  
  condition.specific_celltypes = csc
  
}


celltypes = unique(refs$celltype_toUse)
mm = match(rownames(condition.specific_celltypes), celltypes)
cat(length(which(is.na(mm))), ' cell types missing in the condition.specific_celltypes \n')

source('functions_Visium.R')
## preapre the paramters for RCTD subtypes
DefaultAssay(refs) = 'RNA'
DefaultAssay(st) = 'Spatial'
require_int_SpatialRNA = FALSE
condition.specific.ref = TRUE
RCTD_mode = 'doublet'

RCTD_out = paste0(resDir, '/RCTD_out', 
                  '/RCTD_subtype_out_41subtypes_ref.time.specific_v3.7_ventricleRegion')
system(paste0('mkdir -p ', RCTD_out))


# st = subset(st, condition == 'Amex_d4')
saveRDS(st, file = paste0(RCTD_out, '/RCTD_st_axolotl_allVisium.rds'))

saveRDS(refs, file = paste0(RCTD_out, '/RCTD_refs_subtypes_final_20221117.rds'))

saveRDS(condition.specific_celltypes, 
        file = paste0(RCTD_out, '/RCTD_refs_condition_specificity.rds'))



max_cores = 16
cc = names(table(st$condition))

source('functions_Visium.R')
Run.celltype.deconvolution.RCTD(st = st,
                                refs, 
                                mapping_ref = 'time',
                                condition.specific.ref = condition.specific.ref,
                                condition.specific_celltypes = condition.specific_celltypes,
                                require_int_SpatialRNA = require_int_SpatialRNA,
                                RCTD_mode = RCTD_mode,
                                max_cores = max_cores,
                                RCTD_out = RCTD_out,
                                plot.RCTD.summary = FALSE, 
                                PLOT.scatterpie = FALSE
)


### plot the RCTD outputs
source('functions_Visium.R')

## only ventricle 
RCTD_out = paste0(resDir, '/RCTD_out', 
                  '/RCTD_subtype_out_41subtypes_ref.time.specific_v3.7_ventricleRegion') 

## whole heart
RCTD_out = paste0(resDir, '/RCTD_out', 
                  '/RCTD_allVisium_subtype_out_41subtypes_ref.time.specific_v3.0') 


RCTD_mode = 'doublet'

plot.RCTD.results(st = st,
                  RCTD_out = RCTD_out,
                  RCTD_mode = RCTD_mode,
                  plot.RCTD.summary = FALSE,
                  celltypeProp_cutoff2show = 0.05
                  )


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
      cat('no manual segmentation found for -- ', slice, '\n')
      
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
  
  
  ## add the ventricle segmentation
  st = readRDS(file = paste0(RdataDir, 'seuratObject_allVisiusmst_',
                             'filtered.spots_time_conditions_manualSegmentation', 
                             version.analysis, '.rds'))
  xx = readRDS(file = paste0(RdataDir, 'seuratObject_allVisiusmst_',
                     'filtered.spots_time_conditions_manualSegmentation_ventricleRegions', 
                     version.analysis, '.rds'))
 
  st$seg_ventricle = NA
  st$seg_ventricle = xx$seg2[match(colnames(st), colnames(xx))]
  
  ## all previous manual segmentation were included in ventricle annotation
  jj = which(!is.na(st$segmentation) & is.na(st$seg_ventricle))
  
  saveRDS(st, file = paste0(RdataDir, 'seuratObject_allVisiusmst_',
                            'filtered.spots_time_conditions_manualSegmentation_ventricleRegions', 
                            version.analysis, '.rds'))
  
  
}

##########################################
# Compare replicates from two batches 
# the injury zone, border zone, remote and intact   
##########################################
st = readRDS(file = paste0(RdataDir, 'seuratObject_allVisiusmst_',
                           'filtered.spots_time_conditions_manualSegmentation_ventricleRegions', 
                           version.analysis, '.rds'))


#st = readRDS(file = paste0(RdataDir, 'seuratObject_allVisiusmst_',
#                           'filtered.spots_time_conditions_manualSegmentation', 
#                           version.analysis, '.rds'))

Idents(st) = as.factor(st$segmentation)
SpatialDimPlot(st, ncol = 4)

st$segmentation = as.character(st$segmentation)
st$segmentation[which(st$segmentation == 'RZ1')] = 'RZ'
st$segmentation[which(st$segmentation == 'RZ2')] = 'RZ'
st$segmentation[which(st$segmentation == 'Intact1')] = 'Intact'
st$segmentation[which(st$segmentation == 'Intact2')] = 'Intact'

st$segmentation = as.factor(st$segmentation)

Idents(st) = as.factor(st$segmentation)
SpatialDimPlot(st, ncol = 4)


outDir = paste0(resDir, '/Replicates_comparison_with_segmentation/')
system(paste0('mkdir -p ', outDir))

Test_Data_Integration = FALSE
if(Test_Data_Integration){
  
  Batch_correction_only_segmentation = TRUE
  if(Batch_correction_only_segmentation){
    st = subset(st, cells = which(!is.na(st$segmentation)))
    st = subset(st, cells = which(st$condition != 'Amex_d1_183623' & st$condition != 'Amex_d14_183626'))
    st$condition = droplevels(st$condition)
  }
  
  #st <- SCTransform(st, assay = "Spatial", verbose = FALSE, variable.features.n = 3000)
  # DefaultAssay(st) <- "SCT"
  st = NormalizeData(st, normalization.method = "LogNormalize", scale.factor = 10000)
  st = FindVariableFeatures(st, selection.method = "vst", nfeatures = 5000)
  st = ScaleData(st,  
                 vars.to.regress = c('nCount_Spatial', 'percent.mt'), 
                 assay = 'Spatial')
  
  st <- RunPCA(st, verbose = FALSE)
  ElbowPlot(st, ndims = 30)
  
  st <- FindNeighbors(st, dims = 1:20)
  st <- FindClusters(st, verbose = FALSE, resolution = 1.0)
  
  st <- RunUMAP(st, dims = 1:20, n.neighbors = 30, min.dist = 0.1)
  
  p1 = DimPlot(st, reduction = "umap", group.by = c("condition"), label = TRUE, repel = TRUE) 
  p2 = DimPlot(st, reduction = "umap", group.by = c("segmentation"), label = TRUE, repel = TRUE) 
  
  p1 / p2
  
  ggsave(paste0(outDir, '/Compare_diffRegions_intact_injuried_time_logNormal_regressed.nCounts_reps.pdf'), 
         width = 10, height = 12)
  
  
  st$dataset = 'batch1'
  st$dataset[grep('_2949', st$condition)] = 'batch2'
  
  p0 = DimPlot(st, reduction = "umap", group.by = c("condition"), label = TRUE, repel = TRUE) 
  p1 = DimPlot(st, reduction = "umap", group.by = c("time"), label = TRUE, repel = TRUE) 
  p2 = DimPlot(st, reduction = "umap", group.by = c("dataset"), label = TRUE, repel = TRUE) 
  p3 = DimPlot(st, reduction = "umap", group.by = c("segmentation"), label = TRUE, repel = TRUE) 
  
  (p0 + p2) / (p1 + p3)
  ggsave(paste0(outDir, '/Compare_Replicates_diffRegions_logNormal_regressed.nCounts_reps.pdf'), 
         width = 14, height = 10)
  
  
  source('functions_dataIntegration.R')
  #st_bc = IntegrateData_Seurat_RPCA(st, group.by = 'dataset', k.weight = 50)
  
  st_bc = IntegrateData_runFastMNN(st, 
                                   group.by = 'dataset', 
                                   assays = 'Spatial',
                                   nfeatures = 2000,
                                   correct.all = FALSE)
  
  st_bc <- RunUMAP(st_bc, reduction = "mnn", dims = 1:20, 
                   n.neighbors = 30, min.dist = 0.3)
  
  st_bc$condition = factor(st_bc$condition, levels = c('Amex_d0_294946', 'Amex_d0_294949',
                                                       'Amex_d4_294947', 'Amex_d4_183624', 
                                                       'Amex_d7_294948', 'Amex_d7_183625'))
  
  p1 = DimPlot(st_bc, reduction = "umap", group.by = c("condition"), label = TRUE, repel = TRUE) 
  p2 = DimPlot(st_bc, reduction = "umap", group.by = c("segmentation"), label = TRUE, repel = TRUE) 
  p3 = DimPlot(st_bc, reduction = "umap", group.by = c("time"), label = TRUE, repel = TRUE) 
  p4 = DimPlot(st_bc, reduction = "umap", group.by = c("dataset"), label = TRUE, repel = TRUE) 
  
  p1 / p2
  
  ggsave(paste0(outDir, '/Compare_diffRegions_intact_injuried_time_logNormal_seurat.RunfastMNN_rep.pdf'), 
         width = 10, height = 8)
  
  
  (p2 + p4) / (p3 + p1) 
  ggsave(paste0(outDir, '/Compare_Replicates_diffRegions_logNormal_regressed.nCounts_fastMNN_reps.pdf'), 
         width = 14, height = 10)
  
  
  
  #st$time = gsub('Amex_', '', st$time)
  #st$cc = paste0(st$time, '_', st$segmentation)
  #p3 = DimPlot(st, reduction = "umap", group.by = c("cc"), label = TRUE, repel = TRUE) 
  
  #p1 / p2 /p3
  SpatialPlot(st, group.by = 'segmentation', image.alpha = 0.5, ncol = 4)
  
  ggsave(filename = paste0(resDir, '/UMAP_all.timepoints_', species, '.pdf'), width = 16, height = 8)
  
}


Test_Spatial_alignment = FALSE
if(Test_Spatial_alignment){
  
  outDir = paste0(resDir, '/Replicates_comparison_with_segmentation/')
  system(paste0('mkdir -p ', outDir))
  
  SpatialDimPlot(st, ncol = 4)
  
  SpatialDimPlot(st, group.by = 'seg_ventricle',  ncol = 4)
  mm = which(is.na(st$segmentation) & !is.na(st$seg_ventricle))
  st$segmentation = as.character(st$segmentation)
  st$segmentation[mm] = 'ventricle'
  
  SpatialDimPlot(st, group.by =  'segmentation', ncol = 4)
  
  mm = which(is.na(st$segmentation))
  st$segmentation[mm] = 'others'
  
  SpatialDimPlot(st, group.by =  'segmentation', ncol = 4)
  
  saveRDS(st, file = paste0(outDir, 'seuratObject_allVisiusmst_',
                             'filtered.spots_time_conditions_manualSegmentation_for_scSLAT.rds'))
  
  st = readRDS(file = paste0(outDir, 'seuratObject_allVisiusmst_',
                             'filtered.spots_time_conditions_manualSegmentation_for_scSLAT.rds'))
  
  
  SpatialDimPlot(st, group.by =  'segmentation', ncol = 4)
  
  
  
  library(SeuratDisk)
  
  cc = names(table(st$condition))
  cat('ST conditions : \n')
  print(cc)
  tt = as.character(sapply(cc, function(x){unlist(strsplit(as.character(x), '_'))[2]}))
  cat('ST time :\n')
  print(tt)
  
  Manually_modifying_segmentAnnot = FALSE
  
  for(n in c(4, 5))
  {
    # n = 5
    cat(n, '  slice -- ', cc[n], '\n')
    slice = cc[n]
    Idents(st) = factor(st$condition)
    stx = subset(st, condition == slice)
    DefaultAssay(stx) = 'Spatial'
    
    if(Manually_modifying_segmentAnnot){
      if(slice == 'Amex_d4_294947'){
        p1 = SpatialDimPlot(stx, group.by =  'segmentation', images = slice)
        
        features = rownames(stx)[grep('MYH7-AMEX60DD009525|NPPA-AMEX60DD051098', rownames(stx))]
                
        p2 = SpatialFeaturePlot(st, features = features, images = slice)
        
        p1 + p2
        # now remove additional cells, use SpatialDimPlots to visualize what to remove
        stx$image_row = stx@images$Amex_d4_294947@coordinates$imagerow
        stx$image_col = stx@images$Amex_d4_294947@coordinates$imagecol
        
        SpatialDimPlot(stx, 
                       cells.highlight = WhichCells(stx, expression = segmentation == "ventricle" &
                                                      `MYH6;MYH7-AMEX60DD009525` < 30 & 
                                                      (image_row > 9000 & image_row < 12000) & 
                                                      image_col < 11000),
                       images = slice)
                                                    
        cell_sels = WhichCells(stx, expression = segmentation == "ventricle" &
                                   `MYH6;MYH7-AMEX60DD009525` < 30 & 
                                   (image_row > 9000 & image_row < 12000) & 
                                   image_col < 11000)
        
        stx$segmentation[match(cell_sels, colnames(stx))] = 'Injury'  
        
        SpatialDimPlot(stx, 
                       cells.highlight = WhichCells(stx, expression = segmentation == "others" &
                                                      `MYH6;MYH7-AMEX60DD009525` < 30 & 
                                                      (image_row > 9000 & image_row < 12000) & 
                                                      image_col < 11000),
                       images = slice)
        
        cell_sels = WhichCells(stx, expression = segmentation == "others" &
                     `MYH6;MYH7-AMEX60DD009525` < 30 & 
                     (image_row > 9000 & image_row < 12000) & 
                     image_col < 11000)
        stx$segmentation[match(cell_sels, colnames(stx))] = 'Injury'  
        
        p1 = SpatialDimPlot(stx, group.by =  'segmentation', images = slice)
        
        SpatialDimPlot(stx, 
                       cells.highlight = WhichCells(stx, expression = segmentation == "others" &
                                                      (image_row < 9000 ) & 
                                                      image_col < 12500),
                       images = slice)
        cell_sels = WhichCells(stx, expression = segmentation == "others" &
                                 (image_row < 9000 ) & 
                                 image_col < 12500)
        
        stx$segmentation[match(cell_sels, colnames(stx))] = 'ventricle'  
        
        SpatialDimPlot(stx, 
                       cells.highlight = WhichCells(stx, expression = segmentation == "others" &
                                                      image_col > 15500),
                       images = slice)
        
        cell_sels = WhichCells(stx, expression = segmentation == "others" &
                                 image_col > 15500)
        
        stx$segmentation[match(cell_sels, colnames(stx))] = 'OFT'
        
        SpatialDimPlot(stx, 
                       cells.highlight = WhichCells(stx, expression = segmentation == "others" &
                                                      image_col < 15500),
                       images = slice)
        
        cell_sels = WhichCells(stx, expression = segmentation == "others" &
                                 image_col < 15500)
        
        stx$segmentation[match(cell_sels, colnames(stx))] = 'Atria'
        
        
        saveRDS(stx, file = paste0(outDir, 'ST_visium_manaulSegementaiton_Atria_OFT_', slice, '.rds'))
        
        
        stx = readRDS(file = paste0(outDir, 'ST_visium_manaulSegementaiton_Atria_OFT_', slice, '.rds'))
        SpatialDimPlot(stx, group.by =  'segmentation', images = slice)
        ggsave(paste0(outDir, '/ST_visium_manaulSegementaiton_Atria_OFT_', slice, '.pdf'),
               width = 12, height = 8)
        
        cell_sels = WhichCells(stx, expression = segmentation == "Atria")
        stx$segmentation[match(cell_sels, colnames(stx))] = 'OFT'
        
        cell_sels = WhichCells(stx, expression = segmentation == "others")
        stx$segmentation[match(cell_sels, colnames(stx))] = 'Atria'
        
        cell_sels = WhichCells(stx, expression = segmentation == "RZ")
        stx$segmentation[match(cell_sels, colnames(stx))] = 'ventricle'
        
        SpatialDimPlot(stx, group.by =  'segmentation', images = slice)
        ggsave(paste0(outDir, '/ST_visium_manaulSegementaiton_Atria_OFT_', slice, '.pdf'),
               width = 12, height = 8)
        
        saveRDS(stx, file = paste0(outDir, 'ST_visium_manaulSegementaiton_Atria_OFT_', slice, '_v2.rds'))
        
        
      }
      
      if(slice == 'Amex_d4_183624'){
        p1 = SpatialDimPlot(stx, group.by =  'segmentation', images = slice)
        
        features = rownames(stx)[grep('MYH7-AMEX60DD009525|NPPA-AMEX60DD051098', rownames(stx))]
        
        p2 = SpatialFeaturePlot(st, features = features, images = slice)
        
        p1 + p2
        
        # now remove additional cells, use SpatialDimPlots to visualize what to remove
        stx$image_row = stx@images$Amex_d4_183624@coordinates$imagerow
        stx$image_col = stx@images$Amex_d4_183624@coordinates$imagecol
        
        SpatialDimPlot(stx, 
                       cells.highlight = WhichCells(stx, expression = segmentation == "others" &
                                                      (image_row < 9000) & 
                                                      image_col > 11000),
                       images = slice)
        
        cell_sels = WhichCells(stx, expression = segmentation == "others" &
                                 (image_row < 9000) & 
                                 image_col > 11000)
        
        stx$segmentation[match(cell_sels, colnames(stx))] = 'ventricle'  
        
        SpatialDimPlot(stx, group.by =  'segmentation', images = slice)
        
        
        SpatialDimPlot(stx, 
                       cells.highlight = WhichCells(stx, expression = segmentation == "others" &
                                                      (image_row < 11000) & 
                                                      image_col < 11000),
                       images = slice)
        
        cell_sels = WhichCells(stx, expression = segmentation == "others" &
                                 (image_row < 11000) & 
                                 image_col < 11000)
        
        stx$segmentation[match(cell_sels, colnames(stx))] = 'OFT'  
        
        p1 = SpatialDimPlot(stx, group.by =  'segmentation', images = slice)
        
        SpatialDimPlot(stx, 
                       cells.highlight = WhichCells(stx, expression = segmentation == "OFT" &
                                                      (image_row > 10000 ) & 
                                                      image_col > 9500),
                       images = slice)
        
        cell_sels = WhichCells(stx, expression = segmentation == "OFT" &
                                 (image_row > 10000 ) & 
                                 image_col > 9500)
        
        stx$segmentation[match(cell_sels, colnames(stx))] = 'others'  
        
        SpatialDimPlot(stx, 
                       cells.highlight = WhichCells(stx, expression = segmentation == "others"),
                       images = slice)
        
        cell_sels =  WhichCells(stx, expression = segmentation == "others")
        
        stx$segmentation[match(cell_sels, colnames(stx))] = 'Atria'
        
        SpatialDimPlot(stx, 
                       cells.highlight = WhichCells(stx, expression = segmentation == "others" &
                                                      image_col < 15500),
                       images = slice)
        
        cell_sels = WhichCells(stx, expression = segmentation == "others" &
                                 image_col < 15500)
        
        stx$segmentation[match(cell_sels, colnames(stx))] = 'Atria'
        
        SpatialDimPlot(stx, group.by =  'segmentation', images = slice)
        
        ggsave(paste0(outDir, '/ST_visium_manaulSegementaiton_Atria_OFT_', slice, '.pdf'),
               width = 12, height = 8)
        
        
        saveRDS(stx, file = paste0(outDir, 'ST_visium_manaulSegementaiton_Atria_OFT_', slice, '.rds'))
        
        stx = readRDS(file = paste0(outDir, 'ST_visium_manaulSegementaiton_Atria_OFT_', slice, '.rds'))
        SpatialDimPlot(stx, group.by =  'segmentation', images = slice)
                
                
        cell_sels = WhichCells(stx, expression = segmentation == "RZ")
        stx$segmentation[match(cell_sels, colnames(stx))] = 'ventricle'
        
        SpatialDimPlot(stx, group.by =  'segmentation', images = slice)
        ggsave(paste0(outDir, '/ST_visium_manaulSegementaiton_Atria_OFT_', slice, '.pdf'),
               width = 12, height = 8)
        
        saveRDS(stx, file = paste0(outDir, 'ST_visium_manaulSegementaiton_Atria_OFT_', slice, '_v2.rds'))
        
        
      }
      
    }
    
    coordinates = eval(parse(text = paste0("data.frame(stx@images$", 
                                           slice, 
                                           "@coordinates)[, c(2, 3)]")))
    
    write.csv(coordinates, file = paste0(outDir, 'image_coordinates_', slice, '.csv'), 
              row.names = TRUE, quote = FALSE)
    
    mnt = stx
    VariableFeatures(mnt) = NULL
    
    DefaultAssay(mnt) = 'Spatial'
    
    mnt = DietSeurat(mnt, 
                     counts = TRUE, 
                     data = TRUE,
                     scale.data = FALSE,
                     features = rownames(mnt), 
                     assays = c('Spatial'), 
                     dimreducs = NULL,
                     graphs = NULL, 
                     misc = TRUE
    )
    
    VariableFeatures(mnt)
    
    Idents(mnt) = droplevels(mnt$condition)
    
    saveFile = paste0('ST_visium_', slice, '.h5Seurat')
    
    SaveH5Seurat(mnt, filename = paste0(outDir, saveFile), overwrite = TRUE)
    Convert(paste0(outDir, saveFile), dest = "h5ad", overwrite = TRUE)
    
  }
  
}

########################################################
########################################################
# Section VI: cell neighborhood analysis
# 
########################################################
########################################################
st = readRDS(file = paste0(RdataDir, 'seuratObject_allVisiusmst_',
                           'filtered.spots_time_conditions_manualSegmentation_ventricleRegions', 
                           version.analysis, '.rds'))

Idents(st) = st$seg_ventricle
st = subset(st, cells = colnames(st)[which(!is.na(st$seg_ventricle))])

SpatialPlot(st, group.by = 'seg_ventricle', ncol = 4)

table(st$segmentation, st$condition)

SpatialDimPlot(st, ncol = 4)

## snRNA-seq reference  
refs = readRDS(file = paste0(resDir, '/RCTD_refs_subtypes_final_20221117.rds'))
refs$subtypes = refs$celltype_toUse # clean the special symbols
refs$celltypes = refs$celltype_toUse


table(refs$subtypes)
length(table(refs$subtypes))

refs$subtypes = as.factor(refs$subtypes) 

refs$celltypes = gsub('CM_ven_Robo2', 'CM_Robo2', refs$celltypes)
refs$celltypes = gsub('CM_ven_Cav3_1', 'CM_Cav3.1', refs$celltypes)

st$segmentation = as.character(st$segmentation)
st$segmentation[which(st$segmentation == 'RZ1')] = 'RZ'
st$segmentation[which(st$segmentation == 'RZ2')] = 'RZ'
st$segmentation[which(st$segmentation == 'Intact1')] = 'Intact'
st$segmentation[which(st$segmentation == 'Intact2')] = 'Intact'

st$segmentation = as.factor(st$segmentation)
SpatialDimPlot(st, group.by = 'segmentation', ncol = 4)
table(st$segmentation, st$condition)

RCTD_out = paste0(resDir, '/RCTD_out', 
                  '/RCTD_subtype_out_41subtypes_ref.time.specific_v3.7_ventricleRegion')
outDir = paste0(resDir, '/neighborhood_test/',  basename(RCTD_out), '/')

levels(refs$subtypes)


Run_Neighborhood_Enrichment_Analysis = FALSE
if(Run_Neighborhood_Enrichment_Analysis){
  
  # # condition-specific subtypes selected
  # #condSpec_celltypes = readxl::read_xlsx("../data/neighbourhood_analysis_list_short.xlsx")
  # #condSpec_celltypes = as.data.frame(condSpec_celltypes)
  # condSpec_celltypes = list(d1 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", 'EC_IS_IARS1',
  #                                  "FB_TNXB",'FB_IS_TFPI2',
  #                                  'Mo.Macs_SNX22', "Neu_DYSF",
  #                                  "CM_Cav3.1", "CM_Robo2", 'CM_IS',
  #                                  "Megakeryocytes","RBC"),
  #                           
  #                           d4 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", 'EC_IS_IARS1',
  #                                  "EC_IS_LOX",
  #                                  "FB_PKD1", "FB_TNXB",
  #                                  "Mo.Macs_resident",  "Mo.Macs_FAXDC2", 'Mo.Macs_SNX22', 'Neu_DYSF',
  #                                  "CM_Cav3.1", "CM_Robo2", 'CM_IS',  'CM_Prol_IS',
  #                                  "Megakeryocytes", 'RBC'),
  #                           d7 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", "EC_IS_LOX",
  #                                  "FB_PKD1", "FB_TNXB", "FB_VWA2",
  #                                  "Mo.Macs_resident", "Mo.Macs_FAXDC2", 'Neu_DYSF',
  #                                  "CM_Robo2", "CM_Cav3.1",  'CM_IS',  'CM_Prol_IS',
  #                                  "Megakeryocytes", 'RBC'),
  #                           d14 = c('EC', "EC_CEMIP", "EC_LHX6", 'EC_NOS3', "EC_WNT4", "EC_IS_LOX",
  #                                   "FB_PKD1", "FB_TNXB", "FB_VWA2",
  #                                   "Mo.Macs_resident", 'Neu_DYSF',
  #                                   "CM_Robo2", "CM_Cav3.1",  'CM_IS',
  #                                   "Megakeryocytes", 'RBC')
  # )
  # 
  # celltypes_interest = c()
  # for(n in 1:length(condSpec_celltypes))
  # {
  #   celltypes_interest = unique(c(celltypes_interest, condSpec_celltypes[[n]]))
  # }
  # 
  # condSpec_celltypes$d0 = celltypes_interest
  
  
  source('functions_Visium.R')
  run_misty_colocalization_analysis(st,
                                    outDir = outDir,
                                    RCTD_out = RCTD_out,
                                    condSpec_celltypes = NULL,
                                    segmentation_annots = c('all', 'BZ', 'RZ', 'Intact')
  )
  
  
  ## cell-cell co-localization
  source('functions_Visium.R')
  summarize_cell_neighborhood_misty(st,
                             outDir = outDir, 
                             time = c('d1', 'd4', 'd7', 'd14'),
                             misty_mode = c('density'),
                             summary_method = 'median'
                             #segmentation_annots = c('all', 'BZ', 'RZ', 'Intact'),
                             #controls = c('RZ', 'Intact')
                             )
  
  
}

########################################################
########################################################
# Section VII : ligand-receptor prediction analysis with 
# LIANA and NicheNet
# time-specifc and space-specific niches for nichenet
########################################################
########################################################
st = readRDS(file = paste0(RdataDir, 'seuratObject_allVisiusmst_',
                           'filtered.spots_time_conditions_manualSegmentation_ventricleRegions', 
                           version.analysis, '.rds'))

Idents(st) = st$seg_ventricle
st = subset(st, cells = colnames(st)[which(!is.na(st$seg_ventricle))])

SpatialPlot(st, group.by = 'seg_ventricle', ncol = 4)


table(st$segmentation, st$condition)

SpatialDimPlot(st, ncol = 4)

## snRNA-seq reference  
refs = readRDS(file = paste0(resDir, '/RCTD_refs_subtypes_final_20221117.rds'))
refs$subtypes = refs$celltype_toUse # clean the special symbols
refs$celltypes = refs$celltype_toUse

table(refs$subtypes)
length(table(refs$subtypes))

refs$subtypes = as.factor(refs$subtypes) 

refs$celltypes = gsub('CM_ven_Robo2', 'CM_Robo2', refs$celltypes)
refs$celltypes = gsub('CM_ven_Cav3_1', 'CM_Cav3.1', refs$celltypes)


##########################################
# Part 1) manually specific sub-populations to compare
# sender cells, receiver cells
# BZ-specific and Remote-specific populations (Nichenet specific)
##########################################
# define a list of cell type for each time point, either manual defined or from neighborhood enrichment analysis
#celltypes = c('EC', 'EC_NOS3', 'EC_IS_IARS1', 'FB_IS_TFPI2', 'Mo.Macs_SNX22', 'Neu_IL1R1', 
#              'CM_IS', "RBC")
version_testing_short = FALSE
if(version_testing){
  
  timepoint_specific = TRUE
  celltypes_BZ_timeSpecific = list(day1 = c('EC', 'EC_NOS3', 'EC_IS_IARS1', 'FB_IS_TFPI2', 
                                            'Mo.Macs_SNX22',
                                            'Neu_IL1R1',
                                            'CM_IS', "RBC"),
                                   day4 = c('EC_IS_LOX', 'EC_IS_Prol', 'Mo.Macs_SNX22', 
                                            'Neu_DYSF', 
                                            'CM_IS',
                                            'CM_Prol_IS', 'RBC'),
                                   day7 = c('EC_IS_LOX', 'EC_IS_Prol', 'Mo.Macs_FAXDC2', 'Neu_DYSF', 'Neu_IL1R1',
                                            'CM_IS',
                                            'CM_Prol_IS', 'RBC'),
                                   day14 = c('EC_IS_LOX', 'EC_IS_Prol', 'FB_PKD1', 'Neu_DYSF', 'CM_IS', 'Megakeryocytes',
                                             'RBC')
  )
  celltypes_RZ_timeSpecific = list(day1 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2'),
                                   day4 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2'),
                                   day7 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2'),
                                   day14 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2')
  )
  receivers_BZ_timeSpecific = list(day1 = c("CM_IS"),
                                   day4 = c('CM_Prol_IS', "CM_IS"),
                                   day7 = c('CM_Prol_IS', "CM_IS"),
                                   day14 = c('CM_Prol_IS', "CM_IS")
  )
  
  receivers_RZ_timeSpecific = list(day1 = c("CM_Robo2"),
                                   day4 = c("CM_Robo2"),
                                   day7 = c("CM_Robo2"),
                                   day14 = c("CM_Robo2")
  )
  
}


version_testing_long = FALSE
if(version_testing_long){
  timepoint_specific = TRUE
  
  celltypes_BZ_timeSpecific = list(day1 = c('EC', 'EC_NOS3', "FB_TNXB", 
                                            "CM_Cav3.1", "CM_Robo2", "CM_IS", 
                                            "Megakeryocytes", "RBC"),
                                   day4 = c('EC', 'EC_IS_LOX', "FB_PKD1", 
                                            "Mo.Macs_FAXDC2", 'Mo.Macs_SNX22', 
                                            'Neu_DYSF', "CM_Robo2", 'CM_Prol_IS', "CM_IS",
                                            "Megakeryocytes" ,'RBC'),
                                   day7 = c('EC_NOS3', 'EC_IS_LOX', "FB_PKD1","FB_TNXB", 
                                            "Mo.Macs_FAXDC2", 'Neu_DYSF',
                                            "CM_Robo2", 'CM_Prol_IS', "CM_IS",
                                            "Megakeryocytes" ,'RBC'),
                                   day14 = c('EC_NOS3', "FB_PKD1", "FB_TNXB", 
                                             "Mo.Macs_resident", 'Neu_DYSF', "CM_IS", 'CM_Prol_IS',
                                             "CM_Robo2", 'RBC')
  )
  
  celltypes_RZ_timeSpecific = list(day1 = c('EC', 'EC_NOS3', "FB_TNXB", 'Mo.Macs_SNX22',
                                            "CM_Cav3.1", "CM_Robo2", 'RBC'),
                                   day4 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2', 'CM_Prol_IS', "RBC"),
                                   day7 = c('EC', 'EC_NOS3', 'FB_PKD1', "FB_TNXB", 'CM_Robo2'),
                                   day14 = c('EC', 'EC_NOS3', 'FB_PKD1', 'CM_Robo2', "RBC")
  )
  
  
  receivers_BZ_timeSpecific = list(day1 = c("CM_IS"),
                                   day4 = c('CM_Prol_IS', "CM_IS"),
                                   day7 = c('CM_Prol_IS', "CM_IS"),
                                   day14 = c('CM_Prol_IS', "CM_IS")
  )
  
  receivers_RZ_timeSpecific = list(day1 = c("CM_Robo2"),
                                   day4 = c("CM_Robo2"),
                                   day7 = c("CM_Robo2"),
                                   day14 = c("CM_Robo2")
  )
  
}

##########################################
# Part 2) # run LR analysis for all pairs 
##########################################
# set parameter for ligand-receptor analysis
out_misty =  paste0(resDir, "/neighborhood_test/RCTD_subtype_out_41subtypes_",
                    "ref.time.specific_v3.7_ventricleRegion/signficant_neighborhood_density/")

outDir_version = paste0(resDir, '/Ligand_Receptor_analysis/LIANA_allPairs_ventricleRegion_v3.7/')
if(!dir.exists(outDir_version)) system(paste0('mkdir -p ', outDir_version))

refs$celltypes = refs$celltype_toUse
subtypes = unique(refs$celltypes)

Select_specificPairs = TRUE

#times_slice = c('d1', 'd4', 'd7', 'd14')
times_slice = c('d0_294946', 'd4_294947', 'd7_183625', "d14_183626")

for(n in 1:length(times_slice))
{
  # n = 4
  source('functions_cccInference.R')
  
  time = times_slice[n]
  cat(' run LIANA for time -- ', time, '\n')
  
  outDir = paste(outDir_version, time, collapse = '')
  outDir = gsub(' ', '', outDir)
  
  ## select the interacting subtype pairs  
  pairs = read.table(file = paste0(out_misty, 'cell_cell_colocalization_summary_Amex_', time, '.txt'),
                     sep = '\t', row.names = c(1), header = TRUE)
  
  
  if(Select_specificPairs){
    
    if(time == 'd0_294946'){
      subtypes_sel = c("CM.ven.Robo2", "CM.ven.Cav3.1", "RBC",  "Mo.Macs.resident")
    }
    
    if(time == 'd4_294947') {
      subtypes_sel = c("CM.Prol.IS", "CM.ven.Robo2", 'FB.PKD1', 'Mo.Macs.FAXDC2',  "EC.IS.LOX")
    }
    
    if(time == 'd7_183625'){
      subtypes_sel = c("CM.Prol.IS", "CM.ven.Robo2", 'FB.IS.TNC', 'Mo.Macs.SNX22',  'B.cells.FOXO1',
                        "EC.NOS3")
    }
    
    if(time == 'd14_183626'){
      subtypes_sel = c("CM.ven.Robo2", "EC", "EC.NOS3", 'FB.IS.TNC', "Mo.Macs.resident")
    }
    
    
    ii1 = which(!is.na(match(colnames(pairs), subtypes_sel)))
    jj1 = which(!is.na(match(rownames(pairs), subtypes_sel)))
    
    pairs = pairs[jj1, ii1] 
    pairs[pairs < 0 ] = 0
    
    #pairs[which(rownames(pairs) == 'CM.Prol.IS'), which(rownames(pairs) == 'EC.IS.LOX')] = 1
    pairs = pairs > 0.1
    
    ss_col = apply(pairs, 2, sum)
    ss_row = apply(pairs, 1, sum)
    
    pairs = pairs[which(ss_row>=1), which(ss_col >= 1)] # at least interacting with 1 receivers
    
  }
  # }else{
  #   pairs[is.na(pairs)] = 10
  #   pairs = pairs > 1.6
  #   
  #   ss_row = apply(pairs, 1, sum)
  #   ss_col = apply(pairs, 2, sum)
  #   
  #   pairs = pairs[ ,which(ss_col >= 1)] # at least interacting with 1 receivers
  #   pairs = pairs[which(ss_row >= 1), ] # at least have 3 senders 
  #   
  # }
  
  subtypes_names = gsub('Mo_Macs', 'Mo.Macs', gsub('[.]','_', colnames(pairs)))
  cat(match(subtypes_names, subtypes), '\n')
  
  colnames(pairs) = subtypes_names
  rownames(pairs) = subtypes_names
  
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
            outDir = outDir
  )
  
  source("functions_cccInference.R")
  res = aggregate_output_LIANA(liana_out = outDir)
  
  #require(cellcall)
  library(SeuratData)
  library(Connectome)
  library(cowplot)
  
  res = res[order(-res$sca.LRscore), ]
  
  which(res$ligand == 'GAS6' & res$receptor == 'AXL')
  res[which(res$ligand == 'GAS6' & res$receptor == 'AXL'), ]
  
  colnames(res)[1:2] = c('source', 'target')
  res$weight_norm = res$sca.LRscore
  res$pair = paste0(res$ligand, ' - ', res$receptor)
  res$vector = paste0(res$source, ' - ', res$target)
  res$edge = paste0(res$source, ' - ', res$ligand, ' - ', res$receptor, ' - ', res$target)
  res$source.ligand = paste0(res$source, ' - ', res$ligand)
  res$receptor.target = paste0(res$receptor, ' - ', res$target)
  
  write.table(res, 
              file = paste0(outDir, '/LR_interactions_allPairs_LIANA.txt'), 
              quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')
  
  ## plot the heatmap of cell-cell colocalization
  source('functions_cccInference.R')
  
  pdfname = paste0(outDir, '/LR_interactions_allPairs_LIANA_tops.pdf')
  pdf(pdfname, width=12, height = 8)
  
  for(ntop in c(100, 200, 300))
  {
    # ntop = 300
    test = res[c(1:ntop), ]
    test = test[-which(test$ligand == 'SPON1'), ] ## for unknow reason this ligand making problem for plots
    
    cells.of.interest = unique(c(test$source, test$target))
    cell_color = randomcoloR::distinctColorPalette(length(cells.of.interest))
    names(cell_color) <- cells.of.interest
    
    my_CircosPlot(test, 
                  weight.attribute = 'weight_norm',
                  cols.use = cell_color,
                  sources.include = cells.of.interest,
                  targets.include = cells.of.interest,
                  lab.cex = 0.5,
                  title = paste('LR scores top :', ntop))
    
  }
  
  dev.off()
  
  
}


## double check the ligand and receptor expression distribution
#FeaturePlot(refs, features = rownames(refs)[grep('EGFC|VIPR2', rownames(refs))])

########################################################
##########################################
# Part 3) # run NicheNet analysis for all pairs 
##########################################
##########################################
# diff Nichenet for ligand-receptor analysis
# original code from https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md
########################################################
### all pairs of subtypes
version_testing_all.subtype.pairs = FALSE
if(version_testing_all.subtype.pairs){
  
  # define cell subtype pairs in the border zone 
  celltypes_BZ_timeSpecific = list(day1 = list(CM_IS = c("CM_Cav3.1", "CM_Robo2", "EC_WNT4", 
                                                         'Mo.Macs_SNX22','Neu_IL1R1', "RBC")
                                               #CM_Prol_IS = c("CM_Cav3.1", "CM_Robo2", "EC_WNT4", 
                                               #           'Mo.Macs_SNX22','Neu_IL1R1', "RBC")
  ),
  
  day4 = list(CM_Prol_IS = c('CM_IS',"CM_Robo2", 'EC', 'FB_TNXB',  
                             'Mo.Macs_SNX22', "Mo.Macs_FAXDC2")
              
  ),
  
  day7 = list(CM_Prol_IS = c("CM_Robo2",'EC_IS_LOX', "EC_NOS3",
                             "FB_PKD1", "FB_TNXB",
                             'Mo.Macs_FAXDC2','Neu_DYSF','RBC')),
  
  day14 = list(CM_IS = c("CM_Prol_3", "CM_Robo2", 'EC_IS_LOX', 
                         'FB_PKD1', "FB_TNXB", "Mo.Macs_resident", 'RBC'))
  
  )
  
  # define cell subtype pairs in the remote zone as control for NicheNet
  celltypes_RZ_timeSpecific = list(day1 = list(CM_Robo2 = c("CM_Cav3.1", "EC_IS_IARS1", "EC_WNT4")),
                                   day4 = list(CM_Robo2 = c('CM_Robo2',"FB_TNXB", "Megakeryocytes",
                                                            'Mo.Macs_resident', 'RBC')),
                                   day7 = list(CM_Robo2 = c("EC_NOS3", "EC_WNT4", "Mo.Macs_resident", "RBC")),
                                   day14 = list(CM_Robo2 = c("CM_Cav3.1", "EC_IS_IARS1", "FB_TNXB", 
                                                             "Mo.Macs_resident"))
  )
  
}



outDir_version = paste0(resDir, '/Ligand_Receptor_analysis/DiffNicheNet_v5.1_allpairs_intraOnly')

for(n in 1:length(celltypes_BZ_timeSpecific))
{
  # n = 2
  #source('functions_cccInference.R')
  time = names(celltypes_BZ_timeSpecific)[n]
  cat(' run DiffNicheNet for time -- ', time, '\n')
  outDir = paste(outDir_version, '/', time, collapse = '')
  outDir = gsub(' ', '', outDir)
  
  system(paste0('mkdir -p ', outDir))
  
  source('functions_cccInference_backup.R')
  
  run_Diff_NicheNet(refs = refs, 
                    timepoint_specific = TRUE,
                    include_autocrine = TRUE,
                    celltypes_BZ_specificDay = celltypes_BZ_timeSpecific[[n]],
                    celltypes_RZ_specificDay = celltypes_RZ_timeSpecific[[n]],
                    outDir = outDir
  )
  
  # extract_tables_from_res_Diff_NicheNet(outDir)
  
}

########################################################
########################################################
# Section VIII: additional analysis and plots
# 
########################################################
########################################################

##########################################
# ## quick construction of GAS6-AXL to targets 
##########################################
library(nichenetr)
library(tidyverse)

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
ligand_tf_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_tf_matrix.rds"))

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
sig_network = readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
gr_network = readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))


ligands_all = "GAS6" # this can be a list of multiple ligands if required
targets_all = unique(table_targets$target[which(table_targets$ligand == ligands_all)])

targets_all = c("EGFR", "ETS1","LEF1", 'MYH9', 'SMAD3', 'STAT1')


active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, 
                                                     ligands_all = ligands_all, 
                                                     targets_all = targets_all, 
                                                     weighted_networks = weighted_networks)

# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and 
# gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% 
  mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% 
  mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, 
                                                  ligands_all = ligands_all, 
                                                  targets_all = targets_all, 
                                                  sig_color = "indianred", 
                                                  gr_color = "steelblue")

# To render the graph: uncomment following line of code
DiagrammeR::render_graph(graph_min_max, layout = "kk")


data_source_network = infer_supporting_datasources(signaling_graph_list = active_signaling_network, 
                                                   lr_network = lr_network, 
                                                   sig_network = sig_network, 
                                                   gr_network = gr_network)
head(data_source_network) 



