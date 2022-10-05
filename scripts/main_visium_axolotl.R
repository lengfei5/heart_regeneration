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
version.analysis = '_R12830_resequenced_20220308'
dataDir = '../R12830_visium_reseqenced/nf_out'
resDir = paste0("../results/visium_axolotl", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

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
check.QC.each.condition = TRUE

for(n in 1:nrow(design))
{
  # n = 3
  cat('-----------', design$condition[n], '-------------\n')
  # load nf output and process
  source('functions_Visium.R')
  aa = make_SeuratObj_visium(topdir = paste0(dataDir, '/', design$condition[n], '_', design$sampleID[n], '/'), 
                             saveDir = paste0(resDir, '/', design$condition[n], '_', design$sampleID[n], '/'), 
                             keyname = design$condition[n], 
                             QC.umi = TRUE)
  
  aa$condition = design$condition[n]
  
  #aa <- SCTransform(aa, assay = "Spatial",  method = "glmGamPoi", verbose = FALSE)
  
  aa = subset(aa, subset = nCount_Spatial > 10) # 10 umi from the umi rank
  aa <- SCTransform(aa, assay = "Spatial", verbose = FALSE, variable.features.n = 3000, return.only.var.genes = FALSE)
  
  if(check.QC.each.condition){
    
    pdfname = paste0(resDir, '/QCs_gene_marker_check_', design$condition[n], version.analysis,  '.pdf')
    pdf(pdfname, width=16, height = 8)
    
    # Cell QC metrics: percentage of Mt, nb of counts, nb of genes 
    # get MT% (genes curated from NCBI chrMT genes)
    mtgenes = c("COX1", "COX2", "COX3", "ATP6", "ND1", "ND5", "CYTB", "ND2", "ND4", "ATP8", "MT-CO1", "COI", "LOC9829747")
    mtgenes = c(mtgenes, paste0("MT", mtgenes), paste0("MT-", mtgenes))
    
    ggs = sapply(rownames(aa), function(x) unlist(strsplit(as.character(x), '-'))[1])
    mtgenes = rownames(aa)[!is.na(match(ggs, mtgenes))]
    # mtgenes = mtgenes[mtgenes %in% g[,1]]
    # srat = PercentageFeatureSet(srat, col.name = "percent.mt", assay = "Spatial",
    #                             features = mtgenes)
    xx = PercentageFeatureSet(aa, col.name = "percent.mt", assay = "SCT", features = mtgenes)
    aa[['percent.mt']] = xx$percent.mt
    
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
    
    
    saveRDS(aa, file = paste0(RdataDir, 'seuratObject_design_st_', design$condition[n], version.analysis, '.rds'))
    
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
     file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, version.analysis, '.Rdata'))

########################################################
########################################################
# Section : QC and marker gene checking
# 
########################################################
########################################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, version.analysis, '.Rdata'))

##########################################
# gene and cell filtering (to add later)
##########################################
Filtering.cells.genes = FALSE
if(Filtering.cells.genes){
  #st[["percent.mt"]] <- PercentageFeatureSet(st, pattern = "^Mt-")
  # Visualize QC metrics as a violin plot
  Idents(st) = st$condition
  
  VlnPlot(st, features = c("nCount_Spatial")) + 
    geom_hline( yintercept = c(2500, 5000)) + 
    scale_y_continuous(limits = c(0, 30000))
  ggsave(paste0(resDir, '/QCs_reseq_nUMI.counts.pdf'), width = 12, height = 8)
  
  VlnPlot(st, features = c("nFeature_Spatial")) + 
    geom_hline( yintercept = c(500, 1000)) + 
    scale_y_continuous(limits = c(0, 3000))
  ggsave(paste0(resDir, '/QCs_reseq_nFeatures_Spatial.pdf'), width = 12, height = 8)
  
  FeatureScatter(st, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial")
  
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
  
  SpatialFeaturePlot(st, features = features[2])
  #SpatialFeaturePlot(st, features = features[3])
  
}

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

features = rownames(st)[grep('MYH6|NPPA|CLU-AMEX60DD032706', rownames(st))]
FeaturePlot(st, features = features)

SpatialFeaturePlot(st, features = features[2])

save(design, varibleGenes, st, 
     file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered', species, version.analysis,'.Rdata'))

##########################################
# QCs with marker genes
##########################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered', species, version.analysis,'.Rdata'))

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
  
  markers = c("Ntn1", "Mfap5", "Ubb", "Spon2", "Sparc", "Comp", "Ncanb", "Cthrc1", "Gpr42", "Ffr2")
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
  p0 = SpatialFeaturePlot(st, features = xx[2], image.alpha = 0.5)
  p1 = SpatialDimPlot(st, image.alpha = 0.5)
  plot(wrap_plots(p0, p1, nrow = 2))
  
  pdfname = paste0(resDir, "/check_detected_celltypes_using_AdditionalMarkerGenes_Bern_", species, ".pdf")
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

#features = rownames(st)[grep('COL2A1|COL3A1|AMEX60DD013528|CIRBP', rownames(st))]
features = rownames(st)[grep('AMEX60DD028911|AMEX60DD002594|TGFBI|AMEX60DD029613', rownames(st))]

SpatialFeaturePlot(st, features = features[1], image.alpha = 0.5)

########################################################
########################################################
# Section : start the analysis 
# analyze the axolotl visium per condition
# To identify the boarder zone
########################################################
########################################################
##########################################
# loop over all conditions of st
##########################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered', species, '.Rdata'))
#st = Seurat::SplitObject(st, split.by = 'condition')

st$condition = factor(st$condition, levels = design$condition)

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
  
  
  ggsave(filename =  paste0(resDir, "/Spatial_patterning_", species, '_', cc[n], ".pdf"), width = 12, height = 8)
  
  features = rownames(st)[grep('MYH6|NPPA|CLU-AMEX60DD032706', rownames(st))]
  FeaturePlot(st, features = features)
  
  SpatialFeaturePlot(st, features = features[2])
  
  ##########################################
  # to check if there is border zone, manually select border zone and remote zones from images 
  ##########################################
  source('functions_Visium.R')
  aa = manual_selection_spots_image_Spata(aa, slice = cc[n])
  
  SpatialDimPlot(aa, group.by = 'segmentation', label = TRUE, label.size = 5)
  
  ggsave(filename =  paste0(resDir, "/manual_segmentation_", species, '_', cc[n], ".pdf"), width = 12, height = 8)
  
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
  
  ggsave(filename = paste0(resDir, '/Visium_Clustering_SptialDimPLot', cc[n], '.pdf'), width = 16, height = 8)
  
  # find marker of those clusters
  aa.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.3)
  
  aa.markers %>%
    group_by(cluster) %>%
    top_n(n = 15, wt = avg_log2FC) -> top10
  
  DoHeatmap(aa, features = top10$gene)
  
  ggsave(paste0(resDir, '/heatmap_markerGenes_', cc[n], '.pdf'), width = 12, height = 20)
  
  
}

########################################################
########################################################
# Section IV : cell type deconvolution for each time point with RCTD (not cell2location)
# using the snRNA-seq annotated by Elad
########################################################
########################################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered', species, '.Rdata'))
st$condition = factor(st$condition, levels = design$condition)

# refined subtypes by Elad
refs_file = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets_20221004_2.rds'
refs = readRDS(file = refs_file)
table(refs$subtypes)

# refs0 = readRDS(file ='/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets_20221004.rds')
## NK.cells.or.doublets cluster is weired
#refs = subset(refs, cells = colnames(refs)[grep('doubluets', refs$subtypes, invert = TRUE)])
#refs = subset(refs, cells = colnames(refs)[grep('Neuronal', refs$subtypes, invert = TRUE)])

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

## preapre the paramters for RCTD subtypes
refs$celltype_toUse = refs$subtypes
DefaultAssay(refs) = 'RNA'
DefaultAssay(st) = 'Spatial'
require_int_SpatialRNA = FALSE
RCTD_out = paste0(resDir, '/RCTD_subtype_out_v3.5')
max_cores = 32

# st = subset(st, condition == 'Amex_d4')

source('functions_Visium.R')
Run.celltype.deconvolution.RCTD(st, refs, 
                                require_int_SpatialRNA = require_int_SpatialRNA,
                                max_cores = max_cores,
                                RCTD_out = RCTD_out
)

source('analysis_RCTD_result.R')

########################################################
########################################################
# Section : cell-to-cell communication analysis
# 1) step: define spatial domain 
# 2) step: neighborhood enrichment analysis
# 3) step: ligand-receptor analysis
########################################################
########################################################
source('functions_Visium.R')
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered', species, '.Rdata'))
st$condition = factor(st$condition, levels = design$condition)

cat('visium conditions :\n')
print(table(st$condition))
cc = design$condition

VlnPlot(st, features = 'nFeature_Spatial', group.by = 'condition') +
  geom_hline(yintercept = c(200, 500, 1000, 2000))

ggsave(paste0(resDir, '/QCs_nFeatures_mergedReseq.pdf'), width = 12, height = 8)

VlnPlot(st, features = 'nFeature_SCT', group.by = 'condition') +
  geom_hline(yintercept = c(200, 500, 1000, 2000))
ggsave(paste0(resDir, '/QCs_nFeatures_SCT_mergedReseq.pdf'), width = 12, height = 8)


##########################################
# step 1) Spatial domain searching and potential define remote regions and border zone
# here using computational methods to define regions of interest or cell niches
##########################################
## import manually defined spatial domains
Import.manual.spatial.domains = FALSE
if(Import.manual.spatial.domains){
  require(SPATA2) # installation https://themilolab.github.io/SPATA2/articles/spata-v2-installation.html
  
  st$segmentation = NA
  
  inputDir = '/groups/tanaka/Collaborations/Jingkui-Elad/visium_axolotl_reseq/spata2_manual_regions/'
  #manual_selection_spots_image_Spata
  for(n in 1:length(cc))
  {
    # n = 1
    cat('slice -- ', cc[n], '\n')
    slice = cc[n]
    
    aa = readRDS(file = paste0(inputDir, 'visium_manual_segmentation_', slice, '.rds'))
    Idents(aa) = as.factor(aa$segmentation)
    SpatialDimPlot(aa)
    
    st$segmentation[match(colnames(aa), colnames(st))] = aa$segmentation
        
  }
  
}
st$segmentation = as.factor(st$segmentation)
Idents(st) = as.factor(st$segmentation)
SpatialDimPlot(st)

ggsave(paste0(resDir, '/Manual_segmentation_spata2_Elad.pdf'), width = 16, height = 6)

save(st, design, file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered_manualSegmentation', 
                               species, '.Rdata'))

## run bayesSpace to systematic spatial domain searching
source('functions_Visium.R')
run_bayesSpace(st, outDir = paste0(resDir, '/bayesSpace/'))

##########################################
# step 2) cell proximity analysis 
##########################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered_manualSegmentation', 
                   species, '.Rdata'))
source('functions_Visium.R')
outDir = paste0(resDir, '/neighborhood_test/')
RCTD_out = '../results/visium_axolotl_R12830_resequenced_20220308/RCTD_subtype_out_v3.5'

run_neighborhood_analysis(st, 
                          outDir = outDir,
                          RCTD_out = RCTD_out)


##########################################
# step 3) ligand-receptor-target prediction 
##########################################
source('functions_Visium.R')
refs = readRDS(file = paste0('../results/sc_multiome_R13591_intron.exon.20220729/Rdata/', 
                             'aa_annotated_no_doublets_Elad.rds'))

## NK.cells.or.doublets cluster is weired
refs = subset(refs, cells = colnames(refs)[grep('doubluets', refs$subtypes, invert = TRUE)])
refs = subset(refs, cells = colnames(refs)[grep('Neuronal', refs$subtypes, invert = TRUE)])

refs$celltypes = refs$subtypes

refs$celltypes[grep('CM_|CMs_|_CMs', refs$subtypes)] = 'CM'
refs$celltypes[grep('EC|EC_', refs$subtypes)] = 'EC'
refs$celltypes[grep('FB_', refs$subtypes)] = 'FB'
refs$celltypes[grep('B_cells', refs$subtypes)] = 'Bcell'
refs$celltypes[grep('Macrophages|_MF', refs$subtypes)] = 'Macrophages'
refs$celltypes[grep('Megakeryocytes', refs$subtypes)] = 'Megakeryocytes'
refs$celltypes[grep('RBC', refs$subtypes)] = 'RBC'


run_LIANA()

run_NicheNet()

########################################################
########################################################
# Section IV: spatial organization of cell types and genes  
# 
########################################################
########################################################
#source('functions_Visium.R')
#st = Find.SpatialDE(st)


