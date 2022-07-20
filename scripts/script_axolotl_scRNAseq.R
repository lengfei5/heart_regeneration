##########################################################################
##########################################################################
# Project: heart regeneration
# Script purpose: analyze the scRNA-seq from scMultiome
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Jul 20 15:45:18 2022
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_R13591_20220720'

resDir = paste0("../results/sc_multiome", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R13591_axolotl_multiome'

source('functions_Visium.R')
library(pryr) # monitor the memory usage
require(ggplot2)
mem_used()

species = 'axloltl'

########################################################
########################################################
# Section I: import the scRNAseq data by kalisto
# check QCs 
########################################################
########################################################
design = data.frame(sampleID = seq(197249, 197253), 
                    condition = c(paste0('Amex_scRNA_d', c(0, 1, 4, 7, 14))), stringsAsFactors = FALSE)
varibleGenes = c()
check.QC.each.condition = TRUE

for(n in 1:nrow(design))
{
  # n = 1
  cat('-----------', design$condition[n], '-------------\n')
  
  # load nf output and process
  source('functions_Visium.R')
  
  aa = make_SeuratObj_visium(topdir = paste0(dataDir, '/', design$condition[n], '/'), 
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
