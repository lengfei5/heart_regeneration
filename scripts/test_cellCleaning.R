##########################################################################
##########################################################################
# Project: heart regeneration
# Script purpose: to clean the cells (but keep interesting subpopulation) and 
# to identify the time-specific population 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Nov 21 16:54:43 2022
##########################################################################
##########################################################################
## here is the starting point where there are 48867 cells
aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                           '_lognormamlized_pca_umap_v3.rds'))

# Elad's 37609 cells with 42 subtypes
xx = readRDS(file = paste0('/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_subtypes_final_20221117.rds'))

missed = which(is.na(match(colnames(aa), colnames(xx))))

aa$keep = 'keep'
aa$keep[missed] = 'missed'
aa$keep = as.factor(aa$keep)

# transfer the subtype labels
aa$subtypes = NA
mm = match(colnames(xx), colnames(aa))
aa$subtypes[mm] = as.character(xx$subtypes)

rm(xx)

# aggregate subtypes into major cell types 
aa$celltypes = as.character(aa$subtypes)

aa$celltypes[grep('CM_|CMs_|_CM|_CM_', aa$subtypes)] = 'CM'
aa$celltypes[grep('EC_|_EC', aa$subtypes)] = 'EC'
aa$celltypes[grep('FB_', aa$subtypes)] = 'FB'
aa$celltypes[grep('B_cells', aa$subtypes)] = 'B'
aa$celltypes[grep('T_cells', aa$subtypes)] = 'T'

aa$celltypes[grep('Macrophages|_MF|Macs', aa$subtypes)] = 'Macs'
aa$celltypes[grep('Megakeryocytes', aa$subtypes)] = 'Megakeryocytes'
aa$celltypes[grep('RBC', aa$subtypes)] = 'RBC'
aa$celltypes[grep('Neu_', aa$subtypes)] = 'Neutrophile'


p1 = VlnPlot(aa, features = 'nFeature_RNA', group.by = 'keep',
        y.max = 10000) +
  geom_hline(yintercept = c(500, 2500, 3000))

p2 = VlnPlot(aa, features = 'nCount_RNA', group.by = 'keep', y.max = 50000)
p3 = VlnPlot(aa, features = 'percent.mt', group.by = 'keep', y.max = 100)

p1|p2|p3

ggsave(filename = paste0(resDir, '/compare_keep_discarded_cells.pdf'), width = 20, height = 8)


p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'keep') 
p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes') 
p1 | p2

p3 = FeaturePlot(object = aa, features = 'nFeature_RNA')
p4 = FeaturePlot(aa, features = 'percent.mt') 
p3 | p4

(p1 | p2)/(p3 | p4)

ggsave(filename = paste0(resDir, '/compare_celltypes_nFeatures_pctMT_keep.vs.discarded_cells.pdf'), 
       width = 20, height = 16)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                          '_lognormamlized_pca_umap_keep.missed_subtypes.rds'))

########################################################
########################################################
# Section : we will keep the previous cell filtering (pct.MT < 60% & nFeatures.RNA > 200)  
# preditc doublets
########################################################
########################################################
library(DoubletFinder)
require(Seurat)

aa = readRDS(file = paste0(RdataDir, 'seuratObject_', species, version.analysis, 
                        '_lognormamlized_pca_umap_keep.missed_subtypes.rds'))

cc = unique(aa$condition)

for(n in 1:length(cc))
{
  # n = 1]
  cat(n, ' -- ', cc[n], '\n')
  subs <- subset(aa, condition == cc[n])
  
  subs <- FindVariableFeatures(subs, selection.method = "vst", nfeatures = 5000)
  subs <- ScaleData(subs)
  
  subs <- RunPCA(subs, features = VariableFeatures(object = subs), verbose = TRUE)
  
  subs <- FindNeighbors(subs, dims = 1:30)
  subs <- FindClusters(subs, resolution = 1)
  
  subs <- RunUMAP(subs, dims = 1:30)
  
  sweep.res.list_nsclc <- paramSweep_v3(subs)
  
  sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = FALSE)
  bcmvn_nsclc <- find.pK(sweep.stats_nsclc)
  
  pK <- bcmvn_nsclc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    select(pK) 
  
  pK <- as.numeric(as.character(pK[[1]]))
  annotations <- subs@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  
  nExp_poi <- round(0.076*nrow(subs@meta.data))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  subs <- doubletFinder_v3(subs, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj,  
                           reuse.pANN = FALSE, sct = FALSE)
  
  df_out = subs@meta.data
  subs$DF_out = df_out[, grep('DF.classification', colnames(df_out))]
  
  DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'DF_out',
          raster=FALSE)
  ggsave(filename = paste0(resDir, '/subs_doubletFinder_out_', cc[n], '.pdf'), 
         width = 12, height = 8)
  
  saveRDS(subs, file = paste0(RdataDir, 'subs_doubletFinder_out_', cc[n], '.rds'))
  
}


##########################################
# save the doubletFinder in the main table  
##########################################
cc = unique(aa$condition)
aa$DF_out = NA

for(n in 1:length(cc))
{
  # n = 1
  cat(n, '--', cc[n], '\n')
  subs = readRDS(file = paste0(RdataDir, 'subs_doubletFinder_out_', cc[n], '.rds'))
  aa$DF_out[match(colnames(subs), colnames(aa))] = subs$DF_out
  
}

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_doubletFinderOut_', 
                          species, version.analysis, '.rds'))

##########################################
# Visulize the doublet 
##########################################
aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)
aa <- ScaleData(aa)
aa <- RunPCA(aa, features = VariableFeatures(object = aa), verbose = FALSE)

ElbowPlot(aa, ndims = 30)

Idents(aa) = aa$condition
aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 30, min.dist = 0.1)

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'DF_out', raster=FALSE)
ggsave(filename = paste0(resDir, '/umap_doubletFinder_results.pdf'), width = 12, height = 8)

VlnPlot(aa, features = 'nCount_RNA', group.by = 'DF_out', pt.size = 0.1) +
  geom_hline(yintercept = c(25000, 12500), col = 'red')

ggsave(filename = paste0(resDir, '/nCounts_RNA_double.vs.singlet_doubletFinder_results.pdf'), 
       width = 12, height = 8)

as_tibble(data.frame(condition = aa$condition, group= aa$DF_out)) %>%
  group_by(condition, group) %>% tally() 

pcts = c()
for(n in 1:length(cc))
{
  # n =1
  pcts = c(pcts, length(which(aa$DF_out== 'Doublet' & aa$condition == cc[n]))/length(which(aa$condition == cc[n])))
  
}

data.frame(condition = cc, pct = pcts) %>%
  ggplot(aes(x = condition, y = pct, fill = condition)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none")  + 
  ggtitle('pct of doublets by DF ') + 
  theme(axis.text.x = element_text(angle = 90)) 

ggsave(filename = paste0(resDir, '/Percentages_doublet.vs.total_doubletFinder_results.pdf'), 
       width = 12, height = 8)

saveRDS(aa, file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_doubletFinderOut.v2_', 
                          species, version.analysis, '.rds'))