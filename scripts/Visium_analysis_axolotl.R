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
version.analysis = '_R12830_2020118'
dataDir = '../R12830_visium/nf_out'
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
check.QC.each.condtiion = TRUE

for(n in 1:nrow(design))
{
  # n = 1
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
  
  if(check.QC.each.condtiion){
    
    pdfname = paste0(resDir, '/QCs_gene_marker_check_', design$condition[n], '.pdf')
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
    
    dev.off()
    
    saveRDS(aa, file = paste0(RdataDir, 'seuratObject_design_st_', design$condition[n], '.rds'))
    
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

save(design, varibleGenes, st, file = paste0(RdataDir, 'seuratObject_design_variableGenes_', species, '.Rdata'))

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

st <- FindNeighbors(st, dims = 1:10)
st <- FindClusters(st, verbose = FALSE, resolution = 1.0)

st <- RunUMAP(st, dims = 1:20, n.neighbors = 20, min.dist = 0.05)
SpatialDimPlot(st, image.alpha = 0.5)

DimPlot(st, reduction = "umap", group.by = c("ident", "condition")) 

ggsave(filename = paste0(resDir, '/UMAP_all.timepoints_', species, '.pdf'), width = 16, height = 8)

features = rownames(st)[grep('MYH6|NPPA|CLU-AMEX60DD032706', rownames(st))]
FeaturePlot(st, features = features)

SpatialFeaturePlot(st, features = features[2])

save(design, varibleGenes, st, file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered', species, '.Rdata'))

load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered', species, '.Rdata'))

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
  markers = unique(c(markers, c('hmgb2', 'top2a')))
  
  markers = markers[which(markers != '')]
  markers = firstup(tolower(unique(markers)))
  markers = gsub(' ', '', markers)
  
  xx = c()
  for(n in 1:length(markers))
  {
    xx = c(xx, rownames(st)[grep(toupper(markers[n]), rownames(st))])
  }
  xx = unique(xx)
  markers = xx
  
  #markers[is.na(match(markers, rownames(st)))]
  #markers = markers[!is.na(match(markers, rownames(st)))]
  xx = rownames(st)[grep('CLU', rownames(st))]
  p0 = SpatialFeaturePlot(st, features = xx[2], image.alpha = 0.5)
  p1 = SpatialDimPlot(st, image.alpha = 0.5)
  plot(wrap_plots(p0, p1, nrow = 2))
  
  pdfname = paste0(resDir, "/check_detected_celltypes_using_AdditionalMarkerGenes_", species, ".pdf")
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
# Section : analyze the axolotl visium per condition
# To identify the boarder zone
########################################################
########################################################
##########################################
# loop over all conditions of st
##########################################
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered', species, '.Rdata'))
#st = Seurat::SplitObject(st, split.by = 'condition')

cat('visium conditions :\n')
print(table(design$condition))
cc = design$condition


for(n in 1:length(cc))
#for(n in 1:2)
{
  # n = 2
  aa = readRDS(file = paste0(RdataDir, 'seuratObject_design_st_', cc[n],  '.rds'))
  
  DefaultAssay(aa) <- "SCT"
  
  aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa)
  
  aa <- FindNeighbors(aa, dims = 1:10)
  
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 1.0)
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.05)
  
  p1 = DimPlot(aa, reduction = "umap", group.by = c("ident"), label = TRUE, label.size = 8)
  p2 <- SpatialDimPlot(aa, label = TRUE, label.size = 5)
  
  p1 + p2
  
  ggsave(filename = paste0(resDir, '/Visium_Clustering_SptialDimPLot', cc[n], '.pdf'), width = 16, height = 8)
  
  # find marker of those clusters
  aa.markers <- FindAllMarkers(aa, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.3)
  # saveRDS(aa.markers, file = paste0(RdataDir, 'Forte2020_logNormalize_allgenes_majorCellTypes_markerGenes.rds'))
  
  aa.markers %>%
    group_by(cluster) %>%
    top_n(n = 20, wt = avg_log2FC) -> top10
  
  DoHeatmap(aa, features = top10$gene)
  
  ggsave(paste0(resDir, '/heatmap_markerGenes_', mcells, '_subtypes.pdf'), width = 12, height = 26)
  
  
  
}


