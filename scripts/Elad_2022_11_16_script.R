##########################################################################
##########################################################################
# Project: Heart regeneration 
# Script purpose: for Elad to check marker genes and annotate cell types and subtypes
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Aug  9 14:55:17 2022
##########################################################################
##########################################################################
#rm(list = ls())
require(Seurat)
#require(sctransform)
#library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)

##### AMEX60DD201054292.1 = MPO = Neu

# install.packages("Polychrome")
library("Polychrome")

# build-in color palette
Glasbey = glasbey.colors(32)
swatch(Glasbey)


#Color scheme for small_umap
#661CB0
#5D16A6
#9683EC
#AF91B2
#C79E77
#F7B801
#F4A001
#F18701
#F27103
#F35B04





install.packages("pals")
library(pals)
pal.bands(alphabet, alphabet2, cols25, glasbey, kelly, polychrome, 
          stepped, tol, watlington,
          show.names=FALSE)

install.packages(RColorBrewer)
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))


#####

SeuratObj = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/Rdata_spliced/'
resDir = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/res_Elad'
setwd("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome")

aa = readRDS(file = paste0(SeuratObj, 'seuratObject_axloltl_scRNAseq_R13591_20220720_lognormamlized_pca_umap.rds'))

DimPlot(aa_sub, label = TRUE, repel = TRUE, cols = "RColorBrewer") + ggtitle("scNuc (multiome)")
DimPlot(aa, label = TRUE, repel = TRUE, split.by = "condition") + ggtitle("scNuc (multiome)")
DimPlot(aa_sub, label = TRUE) + NoLegend()
DimPlot(aa_sub, group.by = "subtypes")

features = rownames(aa_sub)[grep('ARG1', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('DLGAP5', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

VlnPlot(aa, features = features)

#ggsave(filename = paste0(resDir, '/FeaturePlot_CD45.pdf'), width = 8, height = 6)
aa
features = rownames(aa)[grep('ITGA2B|MPL-AMEX60DD020059|ITGB3-AMEX60DD009754|PECAM|SELP', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('CD68', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('MYH6|ACTN2|NPPA|TNNT2|GATA4', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

DotPlot(Cluster_5_subest, features = features)$data[,c("features.plot", "id","pct.exp")]
features = rownames(aa)[grep('CD68|CD8A|CD74|CSF1R|ITGAM', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('MKI67|CCNB2|PCNA-|CDK1-', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


aa_var <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(aa), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(aa)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

write.csv(Prolif_CM_marker, "/groups/tanaka/People/current/Elad/Prolif_CM_marker.csv", quote = F)
##export table  

table(Marker$cluster)
##gives back number of values in the table

View(Marker_cluster_0)
#View the tabkle

aa_small <- RunUMAP(aa_small, dims = 1:10)
# runs UMAP

Marker_cluster_0 = FindAllMarkers(aa_small, min.pct = 0.6, min.diff.pct = 0.3)
# finds all DGE 


features = rownames(aa)[grep('COL8A1', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))


Cluster_5_subest <- subset(aa, idents = 5)
#Subseting cluster
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))



aa_small <- FindNeighbors(aa_small, dims = 1:10)

aa_small <- FindClusters(aa_small, resolution = 0.5)

DotPlot(Cluster_5_subest, features = features)$data[,c("features.plot", "id","pct.exp")]
#returns pct of cells expressing gene

pbmc <- RunUMAP(pbmc, dims = 1:10)

DotPlot(aa, features = features)$data[,c("features.plot", "id","pct.exp")]

Marker_cluster_0 = FindAllMarkers(aa_small, min.pct = 0.6, min.diff.pct = 0.3)


View(Marker_cluster_0)

CD45_subset <- subset(aa,  `PTPRC-AMEX60DD018338` > 1)
#backtics are for gene namesCD

DimPlot(CD45_subset, split.by = "condition")
#shows seperated by dataset

CD45_subset_3$condition <- factor(CD45_subset_3$condition,levels = unique(aa$condition))
#sort values according to factor, de facto rearaanged the appearcne of split graph to d0,1,4,7,14 instead of 1,14

FeaturePlot(aa, features = features[1], split.by = "condition")
#shows genes in a split way according to timepoint

DimPlot(Endo_0_subest, group.by = "condition")
#split way according to timepoint


###IF running Harmony - pipeline - NO PCA
Endo_0_subest_Harmony <- RunHarmony(Endo_0_subest, group.by.vars = "condition")

Endo_0_subest_Harmony <- FindNeighbors(Endo_0_subest_Harmony, dims = 1:10, reduction = "harmony")

aa_small <- FindClusters(aa_small, resolution = 0.5)

Endo_0_subest_Harmony <- RunUMAP(Endo_0_subest_Harmony, reduction = "harmony", dims = 1:20)


## doublet finder
library(DoubletFinder)

sweep.res.list_nsclc <- paramSweep_v3 (aa, Pcs = 1:20, sct = False)
sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = False)
#bcmvn_nsclc <- find pK(sweep.stats_nsclc)

ggplot(bcmvn_nsclc, aes(pK, BCmetric, group =1))+
  geom_point()+
  geom_line

#pK <- bcmvn_nsclc %>% filter(BCmetric == max()) %>% select select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- aa@meta.data$seurat_clusters
nExp_poi <- round(0.076*nrow(aa@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
aa_doublet <- doubletFinder_v3(aa, PCs = 1:25, pN = 0.25, pK = pK, nExp = nEXP_poi.adj, reuse.pANN = FALSE, sct = FALSE)
## below final valus calc
aa_doublet <- doubletFinder_v3(aa, PCs = 1:25, pN = 0.25, pK = 0.005, nExp = 3600, reuse.pANN = FALSE, sct = FALSE)


DimPlot(aa_doublet, reduction = 'umap', group.by = "DF.classifications_0.25_0.005_3600")

aa_singlet <- subset(aa_doublet,DF.classifications_0.25_0.005_3600 == "Singlet")



##############################################################
## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list_pbmc <- paramSweep_v3(pbmc.seurat.filtered, PCs = 1:20, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)

ggplot(bcmvn_pbmc, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()

pK <- bcmvn_pbmc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]]))


## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- pbmc.seurat.filtered@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.076*nrow(pbmc.seurat.filtered@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


# run doubletFinder 
pbmc.seurat.filtered <- doubletFinder_v3(pbmc.seurat.filtered, 
                                         PCs = 1:20, 
                                         pN = 0.25, 
                                         pK = pK, 
                                         nExp = nExp_poi.adj,
                                         reuse.pANN = FALSE, sct = FALSE)


# visualize doublets
DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = "DF.classifications_0.25_0.21_691")


# number of singlets and doublets
table(pbmc.seurat.filtered@meta.data$DF.classifications_0.25_0.21_691)
#####################################################################################
aa_day_0 <- subset(aa_gen, condition == "Amex_scRNA_d0")

aa_day_0 <- FindVariableFeatures(aa_day_0, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(aa_day_0)
aa_day_0 <- ScaleData(aa_day_0, features = all.genes)

aa_day_0 <- RunPCA(aa_day_0, features = VariableFeatures(object = aa_day_0))

aa_day_0 <- FindNeighbors(aa_day_0, dims = 1:30)


aa_day_0 <- FindClusters(aa_day_0, resolution = 0.5)

aa_day_0 <- RunUMAP(aa_day_0, dims = 1:30)

sweep.res.list_nsclc <- paramSweep_v3 (aa_day_0)
sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = FALSE)
bcmvn_nsclc <- find.pK(sweep.stats_nsclc)
pK <- bcmvn_nsclc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 


pK <- as.numeric(as.character(pK[[1]]))
annotations <- aa_day_0@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 


nExp_poi <- round(0.076*nrow(aa_day_0@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


aa_day_0 <- doubletFinder_v3(aa_day_0, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj,  reuse.pANN = FALSE, sct = FALSE)





###
aa_day_1 <- subset(aa_gen, condition == "Amex_scRNA_d1")

aa_day_1 <- FindVariableFeatures(aa_day_1, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(aa_day_1)
aa_day_1 <- ScaleData(aa_day_1, features = all.genes)

aa_day_1 <- RunPCA(aa_day_1, features = VariableFeatures(object = aa_day_1))

aa_day_1 <- FindNeighbors(aa_day_1, dims = 1:30)


aa_day_1 <- FindClusters(aa_day_1, resolution = 0.5)

aa_day_1 <- RunUMAP(aa_day_1, dims = 1:30)

sweep.res.list_nsclc <- paramSweep_v3 (aa_day_1)
sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = FALSE)
bcmvn_nsclc <- find.pK(sweep.stats_nsclc)
pK <- bcmvn_nsclc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 


pK <- as.numeric(as.character(pK[[1]]))
annotations <- aa_day_1@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 


nExp_poi <- round(0.076*nrow(aa_day_1@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


aa_day_1 <- doubletFinder_v3(aa_day_1, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj,  reuse.pANN = FALSE, sct = FALSE)

####

aa_day_4 <- subset(aa_gen, condition == "Amex_scRNA_d4")

aa_day_4 <- FindVariableFeatures(aa_day_4, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(aa_day_4)
aa_day_4 <- ScaleData(aa_day_4, features = all.genes)

aa_day_4 <- RunPCA(aa_day_4, features = VariableFeatures(object = aa_day_4))

aa_day_4 <- FindNeighbors(aa_day_4, dims = 1:30)


aa_day_4 <- FindClusters(aa_day_4, resolution = 0.5)

aa_day_4 <- RunUMAP(aa_day_4, dims = 1:30)

sweep.res.list_nsclc <- paramSweep_v3 (aa_day_4)
sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = FALSE)
bcmvn_nsclc <- find.pK(sweep.stats_nsclc)
pK <- bcmvn_nsclc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 


pK <- as.numeric(as.character(pK[[1]]))
annotations <- aa_day_4@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 


nExp_poi <- round(0.076*nrow(aa_day_4@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


aa_day_4 <- doubletFinder_v3(aa_day_4, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj,  reuse.pANN = FALSE, sct = FALSE)








aa_day_7 <- subset(aa_gen, condition == "Amex_scRNA_d7")

aa_day_7 <- FindVariableFeatures(aa_day_7, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(aa_day_7)
aa_day_7 <- ScaleData(aa_day_7, features = all.genes)

aa_day_7 <- RunPCA(aa_day_7, features = VariableFeatures(object = aa_day_7))

aa_day_7 <- FindNeighbors(aa_day_7, dims = 1:30)


aa_day_7 <- FindClusters(aa_day_7, resolution = 0.5)

aa_day_7 <- RunUMAP(aa_day_7, dims = 1:30)

sweep.res.list_nsclc <- paramSweep_v3 (aa_day_7)
sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = FALSE)
bcmvn_nsclc <- find.pK(sweep.stats_nsclc)
pK <- bcmvn_nsclc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 


pK <- as.numeric(as.character(pK[[1]]))
annotations <- aa_day_7@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 


nExp_poi <- round(0.076*nrow(aa_day_7@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


aa_day_7 <- doubletFinder_v3(aa_day_7, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj,  reuse.pANN = FALSE, sct = FALSE)



aa_day_14 <- subset(aa_gen, condition == "Amex_scRNA_d14")

aa_day_14 <- FindVariableFeatures(aa_day_14, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(aa_day_14)
aa_day_14 <- ScaleData(aa_day_14, features = all.genes)

aa_day_14 <- RunPCA(aa_day_14, features = VariableFeatures(object = aa_day_14))


aa_day_14 <- FindNeighbors(aa_day_14, dims = 1:30)

aa_day_14 <- FindClusters(aa_day_14, resolution = 0.5)

aa_day_14 <- RunUMAP(aa_day_14, dims = 1:30)

sweep.res.list_nsclc <- paramSweep_v3 (aa_day_14)
sweep.stats_nsclc <- summarizeSweep(sweep.res.list_nsclc, GT = FALSE)
bcmvn_nsclc <- find.pK(sweep.stats_nsclc)
pK <- bcmvn_nsclc %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 


pK <- as.numeric(as.character(pK[[1]]))
annotations <- aa_day_14@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 


nExp_poi <- round(0.076*nrow(aa_day_14@meta.data))  
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))


aa_day_14 <- doubletFinder_v3(aa_day_14, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj,  reuse.pANN = FALSE, sct = FALSE)

#####
aa_merge <- Merge(aa_day_0_singlet, Y = c(aa_day_1_singlet, aa_day_4_singlet, aa_day_7_singlet, aa_day_14_singlet))

##############
aa_gen_5000 <- subset(aa_gen, subset =  nFeature_RNA < 5000)
aa_gen_marker <- FindAllMarkers(aa_gen_5000, min.pct = 0.5, min.diff.pct =0.3 )

CD45_subset <- subset(aa_gen_5000,  `PTPRC-AMEX60DD018338` > 2)
CD45_subset_2 <- subset(aa_low_comp,  seurat_clusters %in% c(5,4,2))

CD45_subset <- FindVariableFeatures(CD45_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CD45_subset)
CD45_subset <- ScaleData(CD45_subset, features = all.genes)

CD45_subset <- RunPCA(CD45_subset, features = VariableFeatures(object = CD45_subset))

CD45_subset <- FindNeighbors(CD45_subset, dims = 1:10)


CD45_subset <- FindClusters(CD45_subset, resolution = 0.5)

CD45_subset <- RunUMAP(CD45_subset, dims = 1:10)

features = rownames(CD45_subset)[grep('PTPRC', rownames(CD45_subset))]
FeaturePlot(CD45_subset, features = features, cols = c('gray', 'red'))

CD45_subset_cluster0vs2.markers <- FindMarkers(CD45_subset, ident.1 = 0, ident.2 = c(2), min.pct = 0.3)
CD45_subset_cluster4.markers <- FindMarkers(CD45_subset, ident.1 = 4, ident.2 = c(2,0), min.pct = 0.3)
CD45_subset_cluster10.markers <- FindMarkers(CD45_subset, ident.1 = 10, min.pct = 0.25)
CD45_subset_cluster8.markers <- FindMarkers(CD45_subset, ident.1 = "B-Cells_2", min.pct = 0.25)
CD45_subset_cluster_doublets_1.markers <- FindMarkers(CD45_subset, ident.1 = "DOUBLETS", min.pct = 0.25)

new.cluster.ids <- c("Resident Macrophages", "B-Cells", "Macrophages/Monocytes", "Megakeryocytes", "Neutrophills", "T-Cells", "DOUBLETS_EC", "Proliferating_Mf/Mo", "Proliferating_B-Cells", "DOUBLETS_CM" ,"possible NK- PRF1 high")
###Magekeryocytes = MPL, CD42b (GP1BA), CD41 (ITGA2B)
names(new.cluster.ids) <- levels(CD45_subset)
CD45_subset <- RenameIdents(CD45_subset, new.cluster.ids)

########################################## Sub clestring by cluster number and not by PTPRC expression
CD45_subset_2 <- subset(aa_low_comp,  seurat_clusters %in% c(5,4,2))

CD45_subset_2 <- FindVariableFeatures(CD45_subset_2, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CD45_subset_2)
CD45_subset_2 <- ScaleData(CD45_subset_2, features = all.genes)

CD45_subset_2 <- RunPCA(CD45_subset_2, features = VariableFeatures(object = CD45_subset_2))

CD45_subset_2 <- FindNeighbors(CD45_subset_2, dims = 1:15)


CD45_subset_2 <- FindClusters(CD45_subset_2, resolution = 0.5)

CD45_subset_2 <- RunUMAP(CD45_subset_2, dims = 1:15)

CD45_subset_2_marker <- FindAllMarkers(CD45_subset_2, min.pct = 0.6, min.diff.pct = 0.3)

#cluster 0 = Megakeryocyte
#cluster 1 = Resident MF
#cluster 2 = 
#cluster 3 = B Cells
#cluster 4 = Recruited MF
#cluster 5 = 
#cluster 6 = Doublet_EC
#cluster 7 = Doublet_EC
#cluster 8 = Doublet_EC 
#cluster 9 = Doublet_EC
#cluster 10 = Doublet_EC
#cluster 11 = T_cells
#cluster 12 = Proliferating B Cells 
#cluster 13 = Proliferating MF
DimPlot(CD45_subset_2, label = TRUE, repel = TRUE)
features = rownames(CD45_subset_2)[grep('THBS1', rownames(CD45_subset_2))]
FeaturePlot(CD45_subset_2, features = features, cols = c('gray', 'red'))

new.cluster.ids <- c("Resident Macrophages", "B-Cells", "Macrophages/Monocytes", "Megakeryocytes", "Neutrophills", "T-Cells", "DOUBLETS_EC", "Proliferating_Mf/Mo", "Proliferating_B-Cells", "DOUBLETS_CM" ,"possible NK- PRF1 high")
###Magekeryocytes = MPL, CD42b (GP1BA), CD41 (ITGA2B)
names(new.cluster.ids) <- levels(CD45_subset)
CD45_subset <- RenameIdents(CD45_subset, new.cluster.ids)

########################################## Back to gen_cluster annotation

aa_marker <- FindAllMarkers(aa, min.pct = 0.6, min.diff.pct = 0.3)
DimPlot(aa, label = TRUE, repel = TRUE) + ggtitle("scNuc (multiome)")
#cluster 0 = EC
#cluster 1 = CM
#cluster 2 = EC_injury unique
#cluster 3 = CM_ NPPA high/ atrial or BZ
#cluster 4 = Proliferating EC
#cluster 5 = Macrophages
#cluster 6 = EC NOS3 high
#cluster 7 = Erythrocyte
#cluster 8 = EC
#cluster 9 = CM
#cluster 10 = Megakaryocytes
#cluster 11 = CM
#cluster 12 = FB_1
#cluster 13 = FB_2
#cluster 14 = B_cells
#cluster 15 = CM_ NPPA high/ atrial or BZ
#cluster 16 = Proliferating CM
#cluster 17 = CM
#cluster 18 = NOT SURE - Megakaryocytes
#cluster 19 = Proliferating CM + NPPA
#cluster 20 = EC
#cluster 21 = CM (LTBP2)
#cluster 22 = T Cells

#### aa_low_comp = lower number of clusters for sub-cluestring
aa_low_comp <- FindClusters(aa_low_comp, resolution = 0.015)
aa_low_comp <- FindClusters(aa_low_comp, resolution = 0.04) #for CD45 clustering

FB_cluster <- subset(aa_low_comp, idents = 3)

FB_cluster <- FindVariableFeatures(FB_cluster, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(FB_cluster)
FB_cluster <- ScaleData(FB_cluster, features = all.genes)

FB_cluster <- RunPCA(FB_cluster, features = VariableFeatures(object = FB_cluster))

FB_cluster <- FindNeighbors(FB_cluster, dims = 1:10)


FB_cluster <- FindClusters(FB_cluster, resolution = 0.2)

FB_cluster <- RunUMAP(FB_cluster, dims = 1:10)

DimPlot(FB_cluster, label = TRUE, repel = TRUE) 
FB_cluster.marker <-  FindAllMarkers(FB_cluster, min.pct = 0.6, min.diff.pct = 0.3)

##
FB_cluster <- JackStraw(FB_cluster, num.replicate = 100)
FB_cluster <- ScoreJackStraw(FB_cluster, dims = 1:20)
JackStrawPlot(FB_cluster, dims = 1:20)
##

#FB_SUB CLUSTER
#CLUSTER_0 = 
#CLUSTER_1 = 
#CLUSTER_2 = DOUBLET_EC
#CLUSTER_3 = DOUBLET_EC
#CLUSTER_4 = 
#CLUSTER_5 = 
#CLUSTER_6 = proliferating_FB
#CLUSTER_7 = DOUBLET_CM




####### ASCript for JINGKUI
CD45_subset$newId = Idents(CD45_subset)

aa$subtypes = NA
cell.sels = colnames(CD45_subset)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes[mm] = CD45_subset$newId


########## CD45_subset_3 - choosing clusters without 7 (eryth)
CD45_subset_3 <- subset(aa,  seurat_clusters %in% c(14,22,5,10,18))


CD45_subset_3 <- FindVariableFeatures(CD45_subset_3, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CD45_subset_3)
CD45_subset_3 <- ScaleData(CD45_subset_3, features = all.genes)

CD45_subset_3 <- RunPCA(CD45_subset_3, features = VariableFeatures(object = CD45_subset_3))

CD45_subset_3 <- FindNeighbors(CD45_subset_3, dims = 1:15)


CD45_subset_3 <- FindClusters(CD45_subset_3, resolution = 0.5)


CD45_subset_3 <- RunUMAP(CD45_subset_3, dims = 1:15)

CD45_subset_3_marker <- FindAllMarkers(CD45_subset_3, min.pct = 0.6, min.diff.pct = 0.3)

DimPlot(CD45_subset_3, label = TRUE, repel = TRUE)

features = rownames(CD45_subset_3)[grep('TCR', rownames(CD45_subset_3))]
FeaturePlot(CD45_subset_3, features = features,  cols = c('gray', 'red'))

VlnPlot(CD45_subset_3, features = c("PTPRC-AMEX60DD018338"))

#cluster 0 = Resident_MF
#cluster 1 = B_cells
#cluster 2 = Megakeryocytes
#cluster 3 = Macrophages
#cluster 4 = Doublets
#cluster 5 = Proliferating_Megakeryocytes
#cluster 6 = Doublets
#cluster 7 = Neutrophil
#cluster 8 = Doublets
#cluster 9 = T_cells
#cluster 10 = NK_cells? and doubluets
#cluster 11 = Proliferating_Macrophages
#cluster 12 = Neuronal
#cluster 13 = Proliferating_B_cells

Imm.cluster.ids <- c("Resident_MF", "B_cells", "Megakeryocytes", "Mono_Macrophages", "Doublets", "Proliferating_Megakeryocytes",
                     "Doublets", "Neutrophil", "Doublets", "T_cells", "NK_cells?_and_doubluets", "Proliferating_Mono_Macrophages", 
                     "Neuronal" , "Proliferating_B_cells")
names(Imm.cluster.ids) <- levels(CD45_subset_3)
CD45_subset_3 <- RenameIdents(CD45_subset_3, Imm.cluster.ids)

CD45_subset_3$newId = Idents(CD45_subset_3)


#aa$subtypes = NA
cell.sels = colnames(CD45_subset_3)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes[mm] = as.character(CD45_subset_3$newId)

DimPlot(CD45_subset_3, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()+ ggtitle("scNuc (multiome)")
DimPlot(CD45_subset_3, reduction = 'umap', group.by = "seurat_clusters", label = TRUE)


##### EC clusters 0,4,6,2,8,20

EC_subset <- subset(aa,  seurat_clusters %in% c(0,4,6,2,8,20))

EC_subset <- FindVariableFeatures(EC_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(EC_subset)
EC_subset <- ScaleData(EC_subset, features = all.genes)

EC_subset <- RunPCA(EC_subset, features = VariableFeatures(object = EC_subset))

EC_subset <- FindNeighbors(EC_subset, dims = 1:15)


EC_subset <- FindClusters(EC_subset, resolution = 0.15)

EC_subset <- RunUMAP(EC_subset, dims = 1:15)

EC_subset_marker <- FindAllMarkers(EC_subset, min.pct = 0.4, min.diff.pct = 0.2)

DimPlot(EC_subset, label = TRUE, repel = TRUE)
DimPlot(EC_subset, split.by = "condition")


features = rownames(EC_subset)[grep('VCAM1', rownames(EC_subset))]
FeaturePlot(EC_subset, features = features, cols = c('gray', 'red'))

#Cluster 0 = EC
#Cluster 1 = Injury_specific_EC
#Cluster 2 = EC_(VCAM1+)
#Cluster 3 = Proliferating ECs
#Cluster 4 = EC_(Wnt4+)
#Cluster 5 = EC_(LHX6+)


EC.cluster.ids <- c("EC", "Injury_specific_EC", "EC_(VCAM1+)", "Proliferating_ECs", "EC_(Wnt4+)", "EC_(LHX6+)")
names(EC.cluster.ids) <- levels(EC_subset)
EC_subset <- RenameIdents(EC_subset, EC.cluster.ids)

EC_subset$newId = Idents(EC_subset)

#aa$subtypes = NA
cell.sels = colnames(EC_subset)
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes[mm] = as.character(EC_subset$newId)

######## CM clusters 3,21,19,15,11,16,1,9,17

CM_subset <- subset(aa,  seurat_clusters %in% c(3,21,19,15,11,16,1,9,17))

CM_subset <- FindVariableFeatures(CM_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CM_subset)
CM_subset <- ScaleData(CM_subset, features = all.genes)

CM_subset <- RunPCA(CM_subset, features = VariableFeatures(object = CM_subset))

CM_subset <- FindNeighbors(CM_subset, dims = 1:15)


CM_subset <- FindClusters(CM_subset, resolution = 0.1)


CM_subset <- RunUMAP(CM_subset, dims = 1:15)

CM_subset_marker_2 <- FindAllMarkers(CM_subset, min.pct = 0.6, min.diff.pct = 0.3)

DimPlot(CM_subset, label = TRUE, repel = TRUE)
DimPlot(CM_subset, label = TRUE, repel = TRUE, split.by = "condition")


features = rownames(CM_subset)[grep('LHX6', rownames(CM_subset))]
FeaturePlot(CM_subset, features = features, cols = c('gray', 'red'))

#cluster 0 = CMs (ROBO2+) - ventricle
#cluster 1 = NPPA+_CMs - atrial CM
#cluster 2 = Doublets
#cluster 3 = Proliferating_CMs
#cluster 4 = CM (TBX5-, ERBB4+, FGF14+) - out flow track
#cluster 5 = Proliferating_CMs_NPPA+
#cluster 6 = Doublets
#cluster 7 = NPPA+_CMs
#cluster 8 = Doublets

CM.cluster.ids <- c("CM_ventricle", "CM_atrial", "Doublets", "Proliferating_CMs", "CM_OFT", "Proliferating_CMs_NPPA+", "Doublets", "NPPA+_CMs","Doublets")
names(CM.cluster.ids) <- levels(CM_subset)
CM_subset <- RenameIdents(CM_subset, CM.cluster.ids)




CM_subset$newId = Idents(CM_subset)

#aa$subtypes = NA
cell.sels = colnames(CM_subset)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes[mm] = as.character(CM_subset$newId)


######## FB clusters 12,13

FB_subset <- subset(aa,  seurat_clusters %in% c(12,13))

FB_subset <- FindVariableFeatures(FB_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(FB_subset)
FB_subset <- ScaleData(FB_subset, features = all.genes)

FB_subset <- RunPCA(FB_subset, features = VariableFeatures(object = FB_subset))

FB_subset <- FindNeighbors(FB_subset, dims = 1:10)


FB_subset <- FindClusters(FB_subset, resolution = 0.05)


FB_subset <- RunUMAP(FB_subset, dims = 1:10)

FB_subset_marker_1 <- FindAllMarkers(FB_subset, min.pct = 0.6, min.diff.pct = 0.3)

DimPlot(FB_subset, label = TRUE, repel = TRUE)
DimPlot(FB_subset, label = TRUE,  split.by = "condition")


features = rownames(FB_subset)[grep('EPHB1', rownames(FB_subset))]
FeaturePlot(FB_subset, features = features, cols = c('gray', 'red'))

#cluster 0 = FB_1
#cluster 1 = FB_2
#cluster 2 = FB_3
#cluster 3 = Doublets 

FB.cluster.ids <- c("FB_1", "FB_2", "FB_3", "Doublets")
names(FB.cluster.ids) <- levels(FB_subset)
FB_subset <- RenameIdents(FB_subset, FB.cluster.ids)
DimPlot(FB_subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



FB_subset$newId = Idents(FB_subset)

#aa$subtypes = NA
cell.sels = colnames(FB_subset)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes[mm] = as.character(FB_subset$newId)

###########################Eryothrocytes - 
Ery_subset <- subset(aa,  seurat_clusters %in% c(7))

Ery_subset <- FindVariableFeatures(Ery_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Ery_subset)
Ery_subset <- ScaleData(Ery_subset, features = all.genes)

Ery_subset <- RunPCA(Ery_subset, features = VariableFeatures(object = Ery_subset))

Ery_subset <- FindNeighbors(Ery_subset, dims = 1:10)


Ery_subset <- FindClusters(Ery_subset, resolution = 0.3)


Ery_subset <- RunUMAP(Ery_subset, dims = 1:10)

Ery_subset_marker_1 <- FindAllMarkers(Ery_subset, min.pct = 0.6, min.diff.pct = 0.3)

DimPlot(Ery_subset, reduction = "umap", label = TRUE, pt.size = 0.5) 

features = rownames(Ery_subset)[grep('HB-AM-AMEX60DD021531', rownames(Ery_subset))]
FeaturePlot(Ery_subset, features = features, cols = c('gray', 'red'))

#cluster 0 = RBC
#cluster 1 = Proliferating_RBC
#cluster 2 = Doublets
#cluster 3 =  RBC
#cluster 4 =  Doublets


Ery.cluster.ids <- c("RBC", "Proliferating_RBC", "Doublets","RBC", "Doublets")
names(Ery.cluster.ids) <- levels(Ery_subset)
Ery_subset <- RenameIdents(Ery_subset, Ery.cluster.ids)
DimPlot(Ery_subset, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



Ery_subset$newId = Idents(Ery_subset)

#aa$subtypes = NA
cell.sels = colnames(Ery_subset)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes[mm] = as.character(Ery_subset$newId)

features = rownames(aa)[grep('BUB1', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))


######

aa_sub <- subset(aa, cells = names(aa$subtypes)[aa$subtypes %in% c("CMs_(ROBO2+)", "NPPA+_CMs" ,"Proliferating_CMs","CM_(TBX5-,ERBB4+,FGF14+)" , "Proliferating_CMs_NPPA+","EC", "Injury_specific_EC"    ,    "EC_(VCAM1+)" , "Proliferating_ECs"    ,     "EC_(Wnt4+)" , "EC_(LHX6+)","RBC", "Proliferating_RBC" , "FB_1" ,     "FB_2"   ,   "FB_3","Resident_MF"     ,    "B_cells"   ,"Megakeryocytes" , "Mono_Macrophages" ,  "Proliferating_Megakeryocytes" , "Neutrophil", "T_cells"  ,      "NK_cells?_and_doubluets" ,"Proliferating_Mono_Macrophages" , "Neuronal"     ,     "Proliferating_B_cells" )])

aa_sub <- FindVariableFeatures(aa_sub, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(aa_sub)
aa_sub <- ScaleData(aa_sub, features = all.genes)

aa_sub <- RunPCA(aa_sub, features = VariableFeatures(object = aa_sub))

aa_sub <- FindNeighbors(aa_sub, dims = 1:25)


aa_sub <- FindClusters(aa_sub, resolution = 0.5)

aa_sub <- RunUMAP(aa_sub, dims = 1:25)


saveRDS(aa_sub, "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets.rds")

DimPlot(aa, group.by = "subtypes", label = TRUE, split.by = "condition") + NoLegend()
DimPlot(aa_sub, label = TRUE)

DimPlot(aa_sub, group.by = "subtypes", label = TRUE) + NoLegend()

features = rownames(aa_sub)[grep('NKX2-5', rownames(aa_sub))]
FeaturePlot(aa_sub, features = features ,cols = c('gray', 'red'))




#saveRDS(aa_sub, "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets.rds")

aa_cluster_5_mf <- FindMarkers(aa, ident.1 = 5, min.pct = 0.6, min.diff.pct = 0.3)

MF_5 <- subset(aa, seurat_clusters == c(5))

MF_5 <- FindVariableFeatures(MF_5, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(MF_5)
MF_5 <- ScaleData(MF_5, features = all.genes)

MF_5 <- RunPCA(MF_5, features = VariableFeatures(object = MF_5))

MF_5 <- FindNeighbors(MF_5, dims = 1:10)


MF_5 <- FindClusters(MF_5, resolution = 0.5)

MF_5 <- RunUMAP(MF_5, dims = 1:10)

DimPlot(MF_5, group.by = "subtypes", label = TRUE) + NoLegend()

DimPlot(MF_5, label = TRUE) + NoLegend()

MF_5_marker <- FindAllMarkers(MF_5, min.pct = 0.6, min.diff.pct = 0.3)

features = rownames(MF_5)[grep('TOP2A', rownames(MF_5))]
FeaturePlot(MF_5, features = features ,cols = c('gray', 'red'))

### EXPORTING A PDF WHICH IS EDITABLE IN AI
pdf("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_sub.pdf", width = 8)
DimPlot(aa_sub, group.by = "subtypes", label = TRUE) + NoLegend()
dev.off()
#####

VlnPlot(aa, group.by = "orig.ident",  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))


DimPlot(aa_sub, group.by = "seurat_clusters", label = TRUE) + NoLegend()

Ery_marker <- FindMarkers(aa_sub, ident.1 = 8, min.pct = 0.6, min.diff.pct = 0.3)



##### loading Rdataset of G2M and S genes
load("/groups/tanaka/People/current/Paco/Collaborations/Elad_Paco/CellCycleGenes.RData")


Ery_subset <- CellCycleScoring(Ery_subset, s.features = s.genes.new, g2m.features = g2m.genes.new, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(aa_cellcycle[[]])
RidgePlot(aa_cellcycle, features = c("AMEX60DD001288"), ncol = 2)

DimPlot(aa_cellcycle, group.by = "G2M.Score", label = TRUE) + NoLegend()

features = rownames(MF_5)[grep("G2M.Score", rownames(MF_5))]
FeaturePlot(MF_5, features = features ,cols = c('gray', 'red'))

FeaturePlot(MF_subset, c("G2M.Score", "S.Score"))



s.genes.new <- c()
for (i in s.genes) {
  s.genes.new <-c(s.genes.new, rownames(Ery_subset)[grep(i, rownames(Ery_subset))])
}

g2m.genes.new <- c()
for (i in g2m.genes) {
  g2m.genes.new <-c(g2m.genes.new, rownames(Ery_subset)[grep(i, rownames(Ery_subset))])
}

RidgePlot(aa_sub, group.by = "Phase", features = c("TOP2A-AMEX60DD010148"), ncol = 2)


###
BZ_CM <- ScaleData(BZ_CM, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(BZ_CM))


####

Prolif_CM <- subset(aa_sub, subtypes %in% c("Proliferating_CMs_NPPA+", "Proliferating_CMs"))

Prolif_CM <- FindVariableFeatures(Prolif_CM, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Prolif_CM)
Prolif_CM <- ScaleData(Prolif_CM, features = all.genes)

Prolif_CM <- RunPCA(Prolif_CM, features = VariableFeatures(object = Prolif_CM))

Prolif_CM <- FindNeighbors(Prolif_CM, dims = 1:10)


Prolif_CM <- FindClusters(Prolif_CM, resolution = 0.5)

Prolif_CM <- RunUMAP(Prolif_CM, dims = 1:10)

Prolif_CM_marker = FindAllMarkers(Prolif_CM, min.pct = 0.6, min.diff.pct = 0.3)


DimPlot(Prolif_CM, label = TRUE) 

DimPlot(Prolif_CM, split.by = "condition")
DimPlot(Prolif_CM)

Prolif_CM <- SetIdent(Prolif_CM, value = Prolif_CM$seurat_clusters)


features = rownames(Prolif_CM)[grep('RAI14', rownames(Prolif_CM))]

FeaturePlot(Prolif_CM, split.by = "condition", features = features, cols = c('gray', 'red'))
FeaturePlot(Prolif_CM, features = features, cols = c('gray', 'red'))

FeaturePlot(aa_sub, features = features, cols = c('gray', 'red'))

####

features = rownames(aa_sub2)[grep('PTPRC|PECAM|COL1A2|THBS1|ANK1-AMEX60DD002871|NCAM2|TNNI3K', rownames(aa_sub2))]
DotPlot(aa_sub2, group.by = "subtypes", features = features) + RotatedAxis()


FeaturePlot(aa_sub, features = features ,cols = c('gray', 'red'))


aa_sub2 <- subset(aa_sub, cells = names(aa_sub$subtypes)[aa_sub$subtypes != "NK_cells?_and_doubluets"]) 
General_markers = FindAllMarkers(aa_sub2, min.pct = 0.6, min.diff.pct = 0.3)

aa_sub2 <- SetIdent(aa_sub2, value = aa_sub2$subtypes)
DimPlot(aa_sub2)


FeaturePlot(aa_sub2, split.by = "condition", c("G2M.Score"))

####

BZ_CM <- subset(aa_sub2, subtypes %in% c("Proliferating_CMs"))

BZ_CM <- FindVariableFeatures(BZ_CM, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(BZ_CM)
BZ_CM <- ScaleData(BZ_CM, features = all.genes)

BZ_CM <- RunPCA(BZ_CM, features = VariableFeatures(object = BZ_CM))

BZ_CM <- FindNeighbors(BZ_CM, dims = 1:10)


BZ_CM <- FindClusters(BZ_CM, resolution = 0.3)

BZ_CM <- RunUMAP(BZ_CM, dims = 1:10)

BZ_CM_marker = FindAllMarkers(BZ_CM, min.pct = 0.6, min.diff.pct = 0.3)


DimPlot(BZ_CM, split.by = "condition",   label = TRUE) 
DimPlot(BZ_CM,   label = TRUE) 
DimPlot(aa_sub2,   label = TRUE) 


features = rownames(BZ_CM)[grep('ITGB1', rownames(BZ_CM))]
FeaturePlot(BZ_CM, split.by = "condition", features = features, cols = c('gray', 'red'))
FeaturePlot(BZ_CM, features = features, cols = c('gray', 'red'))

FeaturePlot(BZ_CM, features = features, cols = c('gray', 'red'))
FeaturePlot(aa_sub2, features = features, cols = c('gray', 'red'))

###
BZ_CM.cluster.ids <- c("Proliferating_CM1", "Proliferating_CM2", "Proliferating_CM_injury_transient" , "Proliferating_CM_injury", "Proliferating_CM3")
names(BZ_CM.cluster.ids) <- levels(BZ_CM)
BZ_CM <- RenameIdents(BZ_CM, BZ_CM.cluster.ids)
DimPlot(BZ_CM2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


BZ_CM$newId = Idents(BZ_CM)

#aa$subtypes = NA
cell.sels = colnames(BZ_CM)
mm = match(cell.sels, colnames(aa_sub2))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa_sub2$subtypes[mm] = as.character(BZ_CM$newId)

###


##### Ventricluar CM subsetting

Ven_CM <- subset(aa_sub2, subtypes %in% c("CMs_(ROBO2+)"))

Ven_CM <- FindVariableFeatures(Ven_CM, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Ven_CM)
Ven_CM <- ScaleData(Ven_CM, features = all.genes)

Ven_CM <- RunPCA(Ven_CM, features = VariableFeatures(object = BZ_CM))

Ven_CM <- FindNeighbors(Ven_CM, dims = 1:10)


Ven_CM <- FindClusters(Ven_CM, resolution = 0.2)

Ven_CM <- RunUMAP(Ven_CM, dims = 1:10)

Ven_CM_marker = FindAllMarkers(Ven_CM, min.pct = 0.6, min.diff.pct = 0.3)


DimPlot(Ven_CM, split.by = "condition",   label = TRUE) 
DimPlot(Ven_CM,   label = TRUE) 
DimPlot(aa_sub2,  group.by = "subtypes", label = TRUE) 
DimPlot(aa_sub2,  label = TRUE) 

features = rownames(Ven_CM)[grep('ACTC1', rownames(Ven_CM))]
FeaturePlot(Ven_CM, split.by = "condition", features = features, cols = c('gray', 'red'))
FeaturePlot(Ven_CM, features = features, cols = c('gray', 'red'))

FeaturePlot(Ven_CM, features = features, cols = c('gray', 'red'))
FeaturePlot(aa_sub2, features = features, cols = c('gray', 'red'))

###
Ven_CM.cluster.ids <- c("Ventricular_CM_ROBO2+", "Ventricular_CM_Cav3_1+", "Ventricular_CM_injury_specific")
names(Ven_CM.cluster.ids) <- levels(Ven_CM)
Ven_CM <- RenameIdents(Ven_CM, Ven_CM.cluster.ids)
##DimPlot(Ven_CM, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


Ven_CM$newId = Idents(Ven_CM)

#aa$subtypes = NA
cell.sels = colnames(Ven_CM)
mm = match(cell.sels, colnames(aa_sub2))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa_sub2$subtypes[mm] = as.character(Ven_CM$newId)

### renaming the clusters according to visium data


aa_sub2 <- RenameIdents(aa_sub2, 'CM_(TBX5-,ERBB4+,FGF14+)' = "CM_OFT")
aa_sub2 <- RenameIdents(aa_sub2, 'NPPA+_CMs' = "CM_Atria")
aa_sub2 <- RenameIdents(aa_sub2, 'Proliferating_CM1' = "Proliferating_CM",'Proliferating_CM2' = "Proliferating_CM",'Proliferating_CM_injury_transient' = "Proliferating_CM",  'Proliferating_CM_injury' = "Proliferating_CM2")
##version 2
aa_sub2 <- RenameIdents(aa_sub2, 'Proliferating_CM2' = "Proliferating_CM", 'Proliferating_CM3' = "Proliferating_CM" )

saveRDS(aa_sub2, "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets_20221004_2.rds")


####

aa_annotated_no_doublets_20221004$subtypes <- Idents(aa_annotated_no_doublets_20221004)


##### Macrophage  subsetting

MF_subset <- subset(aa, subtypes %in% c("Resident_MF", "Mono_Macrophages",  "Proliferating_Mono_Macrophages"))


MF_subset  <- FindVariableFeatures(MF_subset , selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(MF_subset )
MF_subset  <- ScaleData(MF_subset , features = all.genes)

#MF_subset <- ScaleData(MF_subset, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(MF_subset))


MF_subset  <- RunPCA(MF_subset , features = VariableFeatures(object = MF_subset))

MF_subset  <- FindNeighbors(MF_subset , dims = 1:10)


MF_subset  <- FindClusters(MF_subset , resolution = 0.5)

MF_subset  <- RunUMAP(MF_subset , dims = 1:10)

MF_subset_marker = FindAllMarkers(MF_subset , min.pct = 0.6, min.diff.pct = 0.3)

DimPlot(MF_subset)
DimPlot(MF_subset, split.by = "condition")

features = rownames(MF_subset)[grep('ACAN', rownames(MF_subset))]
FeaturePlot(MF_subset, features = features, cols = c('gray', 'red'))





DimPlot(MF_subset, split.by = "condition", group.by = "subtypes",  label = TRUE) 
DimPlot(MF_subset, split.by = "condition",   label = TRUE) 
DimPlot(MF_subset,   label = TRUE) 
DimPlot(aa_sub2,  group.by = "subtypes", label = TRUE) 
DimPlot(aa_sub2) 

features = rownames(aa)[grep('ACAN', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))
FeaturePlot(aa_sub2, features = features, cols = c('gray', 'red') , order = TRUE)

DefaultAssay(aa_sub2) <- "RNA"

FeaturePlot(Ven_CM, features = features, cols = c('gray', 'red'))
FeaturePlot(aa_sub2, features = features, cols = c('gray', 'red'))


####

###########################Eryothrocytes - 
Ery_subset <- subset(aa,  subtypes %in% c("RBC", "Proliferating_RBC"))

Ery_subset <- FindVariableFeatures(Ery_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Ery_subset)
Ery_subset <- ScaleData(Ery_subset, features = all.genes)

#Ery_subset <- RunPCA(Ery_subset, features = VariableFeatures(object = Ery_subset))
Ery_subset <- RunPCA(Ery_subset, approx=FALSE,  features = c(s.genes, g2m.genes))

Ery_subset <- FindNeighbors(Ery_subset, dims = 1:5)


Ery_subset <- FindClusters(Ery_subset, resolution = 0.1)


Ery_subset <- RunUMAP(Ery_subset, dims = 1:5)

Ery_subset_marker_1 <- FindAllMarkers(Ery_subset, min.pct = 0.4, min.diff.pct = 0.2)

DimPlot(Ery_subset, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "condition") 

features = rownames(Ery_subset)[grep('AMEX60DD008970', rownames(Ery_subset))]
FeaturePlot(Ery_subset, features = features, order = TRUE, cols = c('gray', 'red'))

#### G2M S
load("/groups/tanaka/People/current/Paco/Collaborations/Elad_Paco/CellCycleGenes.RData")

s.genes.new <- c()
for (i in s.genes) {
  s.genes.new <-c(s.genes.new, rownames(aa)[grep(i, rownames(aa))])
}

g2m.genes.new <- c()
for (i in g2m.genes) {
  g2m.genes.new <-c(g2m.genes.new, rownames(aa)[grep(i, rownames(aa))])
}

Ery_subset <- CellCycleScoring(Ery_subset, s.features = s.genes.new, g2m.features = g2m.genes.new, set.ident = TRUE)


Ery_subset <- ScaleData(Ery_subset, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Ery_subset))


aa <- CellCycleScoring(aa, s.features = s.genes.new, g2m.features = g2m.genes.new, set.ident = TRUE)

FeaturePlot(aa, c("G2M.Score"))

DimPlot(aa, group.by = "subtypes")


###

EC_subset <- subset(aa,  subtypes %in% c("EC", "Injury_specific_EC", "EC_(VCAM1+)", "Proliferating_ECs", "EC_(Wnt4+)", "EC_(LHX6+)"))

EC_subset <- FindVariableFeatures(EC_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(EC_subset)
EC_subset <- ScaleData(EC_subset, features = all.genes)

EC_subset <- RunPCA(EC_subset, features = VariableFeatures(object = EC_subset))


EC_subset <- FindNeighbors(EC_subset, dims = 1:10)


EC_subset <- FindClusters(EC_subset, resolution = 0.5)


EC_subset <- RunUMAP(EC_subset, dims = 1:10)

EC_subset_marker_1 <- FindAllMarkers(EC_subset, min.pct = 0.4, min.diff.pct = 0.2)

DimPlot(EC_subset, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(EC_subset, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "subtypes") 
DimPlot(EC_subset, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "condition") 

features = rownames(EC_subset)[grep('AKAP12', rownames(EC_subset))]
FeaturePlot(EC_subset, features = features, order = TRUE, cols = c('gray', 'red'))

VlnPlot(EC_subset, features = features)

VlnPlot(aa, features = features) + NoLegend()

EC_cluster4.markers <- FindMarkers(EC_subset, ident.1 = 4, ident.2 = c(9), min.pct = 0.4)



###
FB_subset <- subset(aa,  subtypes %in% c("FB_1", "FB_2", "FB_3"))

FB_subset <- FindVariableFeatures(FB_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(FB_subset)
FB_subset <- ScaleData(FB_subset, features = all.genes)

FB_subset <- RunPCA(FB_subset, features = VariableFeatures(object = FB_subset))


FB_subset <- FindNeighbors(FB_subset, dims = 1:10)


FB_subset <- FindClusters(FB_subset, resolution = 0.1)


FB_subset <- RunUMAP(FB_subset, dims = 1:10)

FB_subset_marker_1 <- FindAllMarkers(FB_subset, min.pct = 0.6, min.diff.pct = 0.3)

DimPlot(FB_subset, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(FB_subset, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "subtypes") 
DimPlot(FB_subset, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "condition")
DimPlot(FB_subset, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "seurat_clusters")
DimPlot(FB_subset, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "newId")

features = rownames(aa_small)[grep('ERBB', rownames(aa_small))]
FeaturePlot(aa_small,split.by = "condition", features = features, order = TRUE, cols = c('gray', 'red'))

VlnPlot(FB_subset, features = features)

VlnPlot(aa, features = features) + NoLegend()

DimPlot(aa, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(aa, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "subtypes")
#### 2 = doublets


FB.cluster.ids <- c("FB_1","FB_2","Doublet","FB_3","FB_4")
names(FB.cluster.ids) <- levels(FB_subset)
FB_subset <- RenameIdents(FB_subset, FB.cluster.ids)


FB_subset$subtypes = Idents(FB_subset)

#cell.sels = colnames(FB_subset)
##mm = match(cell.sels, colnames(aa))
#cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
#aa$subtypes[mm] <- as.character(FB_subset$subtypes)




cell.sels = colnames(FB_subset)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes = as.character(aa$subtypes)
aa$subtypes[mm] <- as.character(FB_subset$subtypes)
aa$subtypes = as.factor(aa$subtypes)

aa$subtypes <- Idents(aa)
#
FB_subset2 <- subset(aa, cells = names(aa$subtypes)[aa$subtypes %in% c("FB_1","FB_2","FB_3","FB_4")])

FB_subset2 <- FindVariableFeatures(FB_subset2, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(FB_subset2)
FB_subset2 <- ScaleData(FB_subset2, features = all.genes)

FB_subset2 <- RunPCA(FB_subset2, features = VariableFeatures(object = FB_subset2))


FB_subset2 <- FindNeighbors(FB_subset2, dims = 1:10)


FB_subset2 <- FindClusters(FB_subset2, resolution = 0.1)


FB_subset2 <- RunUMAP(FB_subset2, dims = 1:10)


DimPlot(FB_subset2, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(FB_subset2, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "subtypes") 
DimPlot(FB_subset2, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "condition", group.by = "subtypes")
DimPlot(FB_subset2, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "seurat_clusters")
DimPlot(FB_subset2, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "newId")


##
aa_sub <- subset(aa, cells = names(aa$subtypes)[aa$subtypes %in% c("Proliferating_CM" ,"Ventricular_CM_injury_specific",   "Ventricular_CM_ROBO2+", "Ventricular_CM_Cav3_1+","CM_Atria" ,"Proliferating_CMs","CM_OFT" , "Proliferating_CMs_NPPA+","EC", "Injury_specific_EC"    ,    "EC_(VCAM1+)" , "Proliferating_ECs"    ,     "EC_(Wnt4+)" , "EC_(LHX6+)","RBC", "Proliferating_RBC" , "FB_1" ,     "FB_2"   ,   "FB_3", "FB_4","Resident_MF"     ,    "B_cells"   ,"Megakeryocytes" , "Mono_Macrophages" ,  "Proliferating_Megakeryocytes" , "Neutrophil", "T_cells"   ,"Proliferating_Mono_Macrophages" , "Neuronal"     ,     "Proliferating_B_cells" )])




saveRDS(aa_sub, "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets_2022_10_17.rds")

########## aa low complexity for fig 2a




aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(aa)
aa_small <- ScaleData(aa, features = all.genes)

aa <- RunPCA(aa, features = VariableFeatures(object = aa_small))


aa <- FindNeighbors(aa, dims = 1:30)


aa <- FindClusters(aa, resolution = 0.05)


aa <- RunUMAP(aa, dims = 1:30)

DimPlot(aa)
DimPlot(aa, group.by = "seurat_clusters")
DimPlot(aa, group.by = "subtypes")
DimPlot(aa, split.by = "condition")

aa_small_for_subset <- subset(aa, seurat_clusters %in% c(3))

aa_small_for_subset <- FindVariableFeatures(aa_small_for_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(aa_small_for_subset)
aa_small <- ScaleData(aa_small_for_subset, features = all.genes)

aa_small_for_subset <- RunPCA(aa_small_for_subset, features = VariableFeatures(object = aa_small_for_subset))


aa_small_for_subset <- FindNeighbors(aa_small_for_subset, dims = 1:8)


aa_small_for_subset <- FindClusters(aa_small_for_subset, resolution = 0.2)


aa_small_for_subset <- RunUMAP(aa_small_for_subset, dims = 1:8)

DimPlot(aa_small_for_subset)
DimPlot(aa_small_for_subset, group.by = "seurat_clusters")
DimPlot(aa_small_for_subset, group.by = "subtypes")
DimPlot(aa_small_for_subset, split.by = "condition")





DimPlot(aa_small, cols = c())





features = rownames(aa_small)[grep('COL1A2', rownames(aa_small))]
FeaturePlot(aa_small, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))


FB.cluster.ids <- c("EC","CM","Blood","MYeloid","FB","FB", "Bcells")
names(FB.cluster.ids) <- levels(bb)
bb <- RenameIdents(bb, FB.cluster.ids)


Small_markers <- FindAllMarkers(bb, min.pct = 0.6, min.diff.pct = 0.3)

All_Markers = FindAllMarkers(aa, min.pct = 0.6, min.diff.pct = 0.3)



#####

CM_subset_ven <- subset(bb,  subtypes %in% c("Proliferating_CM", "Ventricular_CM_ROBO2+", "Ventricular_CM_Cav3_1+","Ventricular_CM_injury_specific"))
DimPlot(CM_subset)

CM_subset_ven <- FindVariableFeatures(CM_subset_ven, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CM_subset_ven)
CM_subset_ven <- ScaleData(CM_subset_ven, features = all.genes)

CM_subset_ven <- RunPCA(CM_subset_ven, features = VariableFeatures(object = CM_subset_ven))


CM_subset_ven <- FindNeighbors(CM_subset_ven, dims = 1:10)


CM_subset_ven <- FindClusters(CM_subset_ven, resolution = 0.25)


CM_subset_ven <- RunUMAP(CM_subset_ven, dims = 1:10)

DimPlot(CM_subset_ven)
DimPlot(CM_subset_ven, group.by = "seurat_clusters")
DimPlot(CM_subset_ven, group.by = "subtypes")
DimPlot(CM_subset_ven, split.by = "condition")

Idents(CM_subset_ven) <- CM_subset_ven$seurat_clusters

CM_subset_ven_marker = FindAllMarkers(CM_subset_ven, min.pct = 0.6, min.diff.pct = 0.3)
CM_subset_ven_pro_marker <- FindMarkers(CM_subset_ven, ident.1 = c(2), ident.2 = c(3), min.pct = 0.6)
CM_subset_ven_cav_marker <- FindMarkers(CM_subset_ven, ident.1 = c(1), ident.2 = c(4), min.pct = 0.6)


features = rownames(aa_small)[grep('CD3E|PECAM|PAX5|ITGAM|COL1A1|TNNI3K', rownames(aa_small))]
FeaturePlot(aa_small, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))

s.genes.new <- c()
for (i in s.genes) {
  s.genes.new <-c(s.genes.new, rownames(aa_small)[grep(i, rownames(aa_small))])
}

g2m.genes.new <- c()
for (i in g2m.genes) {
  g2m.genes.new <-c(g2m.genes.new, rownames(aa_small)[grep(i, rownames(aa_small))])
}

aa_small <- CellCycleScoring(aa_small, s.features = s.genes.new, g2m.features = g2m.genes.new, set.ident = TRUE)

FeaturePlot(aa_small, c("G2M.Score"), cols = c("#dee2e6", "#661CB0"))


CM_subset_ven.id <- c("Ven_CM","Ven_CM_Cav3_1","Ven_CM_Pro_1_IS","Ven_CM_Pro_2","Ven_CM_IS","Ven_CM_Pro_3")
names(CM_subset_ven.id) <- levels(CM_subset_ven)
CM_subset_ven <- RenameIdents(CM_subset_ven, CM_subset_ven.id )


CM_subset_ven$subtypes = Idents(CM_subset_ven)

saveRDS(aa_small, "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_small_for_major_clusters.rds")

#cell.sels = colnames(FB_subset)
##mm = match(cell.sels, colnames(aa))
#cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
#aa$subtypes[mm] <- as.character(FB_subset$subtypes)




cell.sels = colnames(CM_subset_ven)
mm = match(cell.sels, colnames(bb))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
bb$subtypes = as.character(bb$subtypes)
bb$subtypes[mm] <- as.character(CM_subset_ven$subtypes)
bb$subtypes = as.factor(bb$subtypes)

bb$subtypes <- Idents(bb)

##########
#########
aa_small <- aa

aa_small$subtypes -> Idents(aa_small)
aa_small <- RenameIdents(aa_small,  c('Proliferating_RBC' = "RBC"))


aa_small_markers <- FindAllMarkers(aa_small, min.pct = 0.6, min.diff.pct = 0.3)



features = rownames(aa_small)[grep('ANK1-AMEX60DD002868|TFR2-AMEX60DD01281|VIPR1-AMEX60DD021330|PAX5-AMEX60DD043936|ARHGAP24-AMEX60DD043883|CD79B-AMEX60DD009638|CMKLR1-AMEX60DD000718|MARCO-AMEX60DD056028|CSF1R-AMEX60DD030145|SLC6A4-AMEX60DD054806|KEL-AMEX60DD026372|TMTC1-AMEX60DD006845|CDH11-AMEX60DD015334|COL1A1-AMEX60DD009937|KRT19-AMEX60DD010100|SORBS2-AMEX60DD045413|CORIN-AMEX60DD045587|MYH7B-AMEX60DD028558|ITGB3-AMEX60DD009754|SERPINB10-AMEX60DD038401|SYTL2-AMEX60DD049645|NCAM2-AMEX60DD047429|LRRC75A-AMEX60DD054527|KCNQ5-AMEX60DD033595|GPR171-AMEX60DD001916|LCK-AMEX60DD006076|PTPN7-AMEX60DD008631', rownames(aa_small))]
Idents(aa_small) <- factor(aa_small@active.ident, sort(levels(aa_small@active.ident)))

features_RBC = rownames(aa_small)[grep('ANK1-AMEX60DD002868|TFR2-AMEX60DD01281|VIPR1-AMEX60DD021330', rownames(aa_small))]
features_EC = rownames(aa_small)[grep('SLC6A4-AMEX60DD054806|KEL-AMEX60DD026372|TMTC1-AMEX60DD006845', rownames(aa_small))]
features_FB = rownames(aa_small)[grep('CDH11-AMEX60DD015334|COL1A1-AMEX60DD009937|KRT19-AMEX60DD010100', rownames(aa_small))]
features_CM = rownames(aa_small)[grep('SORBS2-AMEX60DD045413|CORIN-AMEX60DD045587|MYH7B-AMEX60DD028558', rownames(aa_small))]
features_MK = rownames(aa_small)[grep('ITGB3-AMEX60DD009754|SERPINB10-AMEX60DD038401|PLEK-AMEX60DD035848', rownames(aa_small))]
features_Neur = rownames(aa_small)[grep('NCAM2-AMEX60DD047429|LRRC75A-AMEX60DD054527|KCNQ5-AMEX60DD033595', rownames(aa_small))]
features_TCELLS = rownames(aa_small)[grep('GPR171-AMEX60DD001916|LCK-AMEX60DD006076|PTPN7-AMEX60DD008631', rownames(aa_small))]
features_Myeloid = rownames(aa_small)[grep('MARCO-AMEX60DD056028|CSF1R-AMEX60DD030145|CMKLR1-AMEX60DD000718', rownames(aa_small))]
features_B_cells = rownames(aa_small)[grep('PAX5-AMEX60DD043936|ARHGAP24-AMEX60DD043883|CD79B-AMEX60DD009638', rownames(aa_small))]




DotPlot(aa_small, features = c( features_CM,features_FB,features_EC , features_TCELLS, features_B_cells, features_Myeloid,features_RBC , features_MK, features_Neur))+ RotatedAxis()

DotPlot(aa_small, features = features) + RotatedAxis()


write.csv(aa_small_markers, "/groups/tanaka/People/current/Elad/aa_small_markers", quote = F)

#### sorting the values
myLevels <- c("CM",
              "FB",
              "EC",
              "T_cells",
              "B_cells",
              "Myeloid",
              "RBC",
              "Megakaryocytes",
              "Neuronal"
)



factor(Idents(aa_small), levels= myLevels)
Idents(aa_small) <- factor(Idents(aa_small), levels= myLevels)


DotPlot(aa_small, features = c( features_CM,features_FB,features_EC , features_TCELLS, features_B_cells, features_Myeloid,features_RBC , features_MK, features_Neur)) + RotatedAxis()
############
######## CM clusters 3,21,19,15,11,16,1,9,17

CM_subset_gen <- subset(aa,  subtypes %in% c("Ventricular_CM_Cav3_1+", "CM_Atria",   "CM_OFT", "Ventricular_CM_ROBO2+" ))

CM_subset_gen <- FindVariableFeatures(CM_subset_gen, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CM_subset_gen)
CM_subset_gen <- ScaleData(CM_subset_gen, features = all.genes)

CM_subset_gen <- RunPCA(CM_subset_gen, features = VariableFeatures(object = CM_subset_gen))

CM_subset_gen <- FindNeighbors(CM_subset_gen, dims = 1:10)


CM_subset_gen <- FindClusters(CM_subset_gen, resolution = 0.1)


CM_subset_gen <- RunUMAP(CM_subset_gen, dims = 1:10)

CM_subset_marker_no_div <- FindAllMarkers(CM_subset_gen, min.pct = 0.6, min.diff.pct = 0.3)

CM.cluster.ids <- c("Ven_CM_Robo2", "CM_Atria", "Ven_CM_(Cav3_1)", "CM_Atria_Tagln", "CM_OFT", "CM_PM_(HCN4)")
names(CM.cluster.ids) <- levels(CM_subset_gen)
CM_subset_gen <- RenameIdents(CM_subset_gen, CM.cluster.ids)
CM_subset_gen$subtypes <- Idents(CM_subset_gen)

cell.sels = colnames(CM_subset_gen)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes = as.character(aa$subtypes)
aa$subtypes[mm] <- as.character(CM_subset_gen$subtypes)
aa$subtypes = as.factor(aa$subtypes)

DimPlot(CM_subset, label = TRUE, repel = TRUE) +NoLegend()
DimPlot(CM_subset,  repel = TRUE, split.by = "condition", group.by = "subtypes")
DimPlot(CM_subset,  repel = TRUE, group.by = "subtypes")

features = rownames(CM_subset)[grep('ADAMTS17', rownames(CM_subset))]
FeaturePlot(CM_subset, features = features, cols = c('gray', 'red'))

s.genes.new <- c()
for (i in s.genes) {
  s.genes.new <-c(s.genes.new, rownames(EC_subest_2)[grep(i, rownames(EC_subest_2))])
}

g2m.genes.new <- c()
for (i in g2m.genes) {
  g2m.genes.new <-c(g2m.genes.new, rownames(EC_subest_2)[grep(i, rownames(EC_subest_2))])
}

EC_subest_2 <- CellCycleScoring(EC_subest_2, s.features = s.genes.new, g2m.features = g2m.genes.new, set.ident = TRUE)

FeaturePlot(EC_subest_2, c("S.Score"), cols = c("#dee2e6", "#661CB0"))



CM.cluster.ids <- c("Ven_CM_Robo2", "CM_Atria", "Ven_CM_(Cav3_1)", "CM_Atria_Tagln", "CM_OFT", "CM_PM_(HCN4)")
names(CM.cluster.ids) <- levels(CM_subset)
CM_subset <- RenameIdents(CM_subset, CM.cluster.ids)




CM_subset$newId = Idents(CM_subset)


saveRDS(CM_subset, "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/CM_subset_no_div.rds")


features_Robo2 = rownames(CM_subset)[grep('ROBO2-AMEX60DD047486|JPH3-AMEX60DD017176|WSCD2-AMEX60DD000716', rownames(CM_subset))]
features_Cav = rownames(CM_subset)[grep('RBPMS-AMEX60DD043590|ERV3-1-AMEX60DD016594|SEPTIN9-AMEX60DD030744', rownames(CM_subset))]
features_OFT = rownames(CM_subset)[grep('LHCGR-AMEX60DD035890|CRISPLD2-AMEX60DD017214|THSD4-AMEX60DD014724', rownames(CM_subset))]
features_A1 = rownames(CM_subset)[grep('GFRA1-AMEX60DD053205|ADAMTS17-AMEX60DD003517|AGBL1-AMEX60DD004208', rownames(CM_subset))]
features_Tagln = rownames(CM_subset)[grep('NPPA-AMEX60DD051099|TAGLN-AMEX60DD053922|ADAM8-AMEX60DD052283', rownames(CM_subset))]
features_Ltbp2 = rownames(CM_subset)[grep('LTBP2-AMEX60DD011463|KCTD12-AMEX60DD022831|HCN4-AMEX60DD003883', rownames(CM_subset))]


CMmyLevels <- c("CM_Atria_Ltbp2",
                "CM_Atria_Tagln",
                "CM_Atria_1",
                "CM_OFT",
                "Ventricular_CM_Cav3_1",
                "Ven_CM_Robo2"
                
)




factor(Idents(CM_subset), levels= CMmyLevels)
Idents(CM_subset) <- factor(Idents(CM_subset), levels= CMmyLevels)


DotPlot(CM_subset, col.min = 0, features = c( features_Robo2,features_Cav,features_OFT,features_A1,features_Tagln,features_Ltbp2)) + RotatedAxis()
############

DimPlot(aa_small, cols = c(
  "#E60B5F",
  "#C8216F",
  "#AD377E",
  "#934D90",
  "#725F9B",
  "#5472A6",
  "#398AB9",
  "#1DA3CB",
  "#00BBE0") )

#


EC_subset <- subset(aa,  subtypes %in% c("EC", "EC_(Wnt4+)",   "EC_(VCAM1+)", "EC_(LHX6+)" , "Injury_specific_EC"))

EC_subset <- FindVariableFeatures(EC_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(EC_subset)
EC_subset <- ScaleData(EC_subset, features = all.genes)

EC_subset <- RunPCA(EC_subset, features = VariableFeatures(object = EC_subset))

EC_subset <- FindNeighbors(EC_subset, dims = 1:10)


EC_subset <- FindClusters(EC_subset, resolution = 0.19)


EC_subset <- RunUMAP(EC_subset, dims = 1:10)

EC_subset_marker_no_div <- FindAllMarkers(EC_subset, min.pct = 0.4, min.diff.pct = 0.2)


EC.cluster.ids <- c("EC", "EC_(VCAM1)", "EC", "EC_(WNT4)", "EC_INJ", "EC_(LHX6)")
names(EC.cluster.ids) <- levels(EC_subset)
EC_subset <- RenameIdents(EC_subset, EC.cluster.ids)

EC_subset$subtypes = Idents(EC_subset)


DimPlot(EC_subset, label = TRUE, repel = TRUE) +NoLegend()
DimPlot(EC_subset,  repel = TRUE, split.by = "condition", group.by = "subtypes")
DimPlot(EC_subset,  repel = TRUE, split.by = "condition")
DimPlot(EC_subset,  repel = TRUE, group.by = "subtypes")

features = rownames(EC_subset)[grep('PLAT-AMEX60DD002881', rownames(EC_subset))]
FeaturePlot(EC_subset, features = features, cols = c('gray', 'red'))

EC_subest_2 <- subset(EC_subset,  subtypes %in% c("EC", "EC_(VCAM1)", "EC_(WNT4)", "EC_(LHX6)"))


EC_subest_2 <- FindVariableFeatures(EC_subest_2, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(EC_subest_2)
EC_subest_2 <- ScaleData(EC_subest_2, features = all.genes)

EC_subest_2 <- RunPCA(EC_subest_2, features = VariableFeatures(object = EC_subest_2))

EC_subest_2 <- FindNeighbors(EC_subest_2, dims = 1:10)


EC_subest_2 <- FindClusters(EC_subest_2, resolution = 0.19)


EC_subest_2 <- RunUMAP(EC_subest_2, dims = 1:10)

DimPlot(EC_subest_2, label = TRUE, repel = TRUE) +NoLegend()
DimPlot(EC_subest_2,  repel = TRUE, split.by = "condition", group.by = "subtypes")
DimPlot(EC_subest_2,  repel = TRUE, split.by = "condition")
DimPlot(EC_subest_2,  repel = TRUE, group.by = "subtypes")

EC_subset_2_marker_no_div <- FindAllMarkers(EC_subest_2, min.pct = 0.4, min.diff.pct = 0.2)


features = rownames(aa)[grep('CRISPLD2', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

EC.cluster.ids <- c("EC", "EC_(NOS3)", "EC_(WNT4)", "EC_(VCAM1)", "EC","EC_(LHX6)")
names(EC.cluster.ids) <- levels(EC_subest_2)
EC_subest_2 <- RenameIdents(EC_subest_2, EC.cluster.ids)

DimPlot(EC_subest_2)

saveRDS(EC_subest_2, "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/EC_subset_2_2022_10_31.rds")

features_EC = rownames(EC_subest_2)[grep('LAMA2-AMEX60DD034687|KCNMB2-AMEX60DD030396|RADIL-AMEX60DD025680', rownames(CM_subset))]
features_NOS3 = rownames(EC_subest_2)[grep('NOS3-AMEX60DD026163|ADAM15-AMEX60DD014884|HSPA12B-AMEX60DD046249', rownames(CM_subset))]
features_WNT4 = rownames(EC_subest_2)[grep('KCNQ5-AMEX60DD033595|HEG1-AMEX60DD009667|WNT4-AMEX60DD052091', rownames(CM_subset))]
features_VCAM1 = rownames(EC_subest_2)[grep('NTRK3-AMEX60DD004199|VCAM1-AMEX60DD018837|CA8-AMEX60DD039785', rownames(CM_subset))]
features_LHX6 = rownames(EC_subest_2)[grep('PDE2A-AMEX60DD049847|ADGRF5-AMEX60DD032958|LHX6-AMEX60DD050778', rownames(CM_subset))]


ECmyLevels <- c( "EC",
                 "EC_(NOS3)",
                 "EC_(WNT4)",
                 "EC_(VCAM1)",
                 "EC_(LHX6)"
                 
)




factor(Idents(EC_subest_2), levels= ECmyLevels)
Idents(EC_subest_2) <- factor(Idents(EC_subest_2), levels= ECmyLevels)


#DotPlot(EC_subest_2, col.min = 0, features = c(features_EC,features_NOS3,features_WNT4,features_VCAM1,features_LHX6)) + RotatedAxis()
DotPlot(EC_subest_2, col.min = 0, features = c(features_LHX6,features_VCAM1,features_WNT4,features_NOS3,features_EC)) + RotatedAxis() & coord_flip()
############



FB_subset <- subset(aa,  subtypes %in% c("FB_1", "FB_2",   "FB_3", "FB_4" ))

FB_subset <- FindVariableFeatures(FB_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(FB_subset)
FB_subset <- ScaleData(FB_subset, features = all.genes)

FB_subset <- RunPCA(FB_subset, features = VariableFeatures(object = FB_subset))

FB_subset <- FindNeighbors(FB_subset, dims = 1:10)


FB_subset <- FindClusters(FB_subset, resolution = 0.1)


FB_subset <- RunUMAP(FB_subset, dims = 1:10)

FB_subset_marker_no_div <- FindAllMarkers(FB_subset, min.pct = 0.4, min.diff.pct = 0.2)


FB.cluster.ids <- c("FB_1", "FB_2", "FB_3", "FB_4", "FB_5", "FB_6")
names(FB.cluster.ids) <- levels(FB_subset)
FB_subset <- RenameIdents(FB_subset, FB.cluster.ids)

FB_subset$subtypes = Idents(FB_subset)


DimPlot(FB_subset, label = TRUE, repel = TRUE) +NoLegend()
DimPlot(FB_subset,  repel = TRUE, split.by = "condition", group.by = "subtypes")
DimPlot(FB_subset,  repel = TRUE, split.by = "condition")
DimPlot(FB_subset,  repel = TRUE, group.by = "subtypes")

features = rownames(FB_subset)[grep('PLAT-AMEX60DD002881', rownames(FB_subset))]
FeaturePlot(FB_subset, features = features, cols = c('gray', 'red'))





saveRDS(FB_subset, "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/FB_subset_2022_10_31.rds")

features_1 = rownames(FB_subset)[grep('EFEMP1-AMEX60DD035438|PKD1-AMEX60DD030949|PODXL-AMEX60DD008079', rownames(CM_subset))]
features_2 = rownames(FB_subset)[grep('TNXB-AMEX60DD010491|LAMA2-AMEX60DD034687|NGF-AMEX60DD008715', rownames(CM_subset))]
features_3 = rownames(FB_subset)[grep('ACAN-AMEX60DD004178|VWA2-AMEX60DD053181|CSPG5-AMEX60DD026060', rownames(CM_subset))]
features_4 = rownames(FB_subset)[grep('TFPI2-AMEX60DD022268|CXCL14-AMEX60DD028973|HAS1-AMEX60DD017885', rownames(CM_subset))]
features_5 = rownames(FB_subset)[grep('CKAP2L-AMEX60DD001288|IQGAP3-AMEX60DD016824|MASTL-AMEX60DD021858', rownames(CM_subset))]
features_6 = rownames(FB_subset)[grep('MPL-AMEX60DD020059|ADGRG1-AMEX60DD015471|PLEK-AMEX60DD035848', rownames(CM_subset))]


FBmyLevels <- c( "FB_1",
                 "FB_2",
                 "FB_3",
                 "FB_4",
                 "FB_5",
                 "FB_6"
                 
)




factor(Idents(FB_subset), levels= FBmyLevels)
Idents(FB_subset) <- factor(Idents(FB_subset), levels= FBmyLevels)


#DotPlot(EC_subest_2, col.min = 0, features = c(features_EC,features_NOS3,features_WNT4,features_VCAM1,features_LHX6)) + RotatedAxis()
DotPlot(FB_subset, col.min = 0, features = c(features_6,features_5,features_4,features_3,features_2,features_1)) + RotatedAxis() & coord_flip()
############


FB_subset_2 <- subset(FB_subset,  subtypes %in% c("FB_1", "FB_2",  "FB_3", "FB_5" ))

FB_subset_2 <- FindVariableFeatures(FB_subset_2, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(FB_subset_2)
FB_subset_2 <- ScaleData(FB_subset_2, features = all.genes)

FB_subset_2 <- RunPCA(FB_subset_2, features = VariableFeatures(object = FB_subset_2))

FB_subset_2 <- FindNeighbors(FB_subset_2, dims = 1:5)


FB_subset_2 <- FindClusters(FB_subset_2, resolution = 0.08)


FB_subset_2 <- RunUMAP(FB_subset_2, dims = 1:5)

FB_subset_2_marker_no_div <- FindAllMarkers(FB_subset_2, min.pct = 0.4, min.diff.pct = 0.2)


FB2.cluster.ids <- c("FB_1", "FB_2", "FB_3", "FB_4")
names(FB2.cluster.ids) <- levels(FB_subset_2)
FB_subset_2 <- RenameIdents(FB_subset_2, FB2.cluster.ids)

features_1 = rownames(FB_subset_2)[grep('EFEMP1-AMEX60DD035438|PKD1-AMEX60DD030949|CFB-AMEX60DD010523', rownames(CM_subset))]
features_2 = rownames(FB_subset_2)[grep('TNXB-AMEX60DD010491|LAMA2-AMEX60DD034687|NGF-AMEX60DD008715', rownames(CM_subset))]
features_3 = rownames(FB_subset_2)[grep('ACAN-AMEX60DD004178|VWA2-AMEX60DD053181|CSPG5-AMEX60DD026060', rownames(CM_subset))]
#features_4 = rownames(FB_subset_2)[grep('TFPI2-AMEX60DD022268|CXCL14-AMEX60DD028973|HAS1-AMEX60DD017885', rownames(CM_subset))]
features_4 = rownames(FB_subset_2)[grep('CKAP2L-AMEX60DD001288|IQGAP3-AMEX60DD016824|MASTL-AMEX60DD021858', rownames(CM_subset))]
#features_6 = rownames(FB_subset_2)[grep('MPL-AMEX60DD020059|ADGRG1-AMEX60DD015471|PLEK-AMEX60DD035848', rownames(CM_subset))]


FB2myLevels <- c( "FB_1",
                  "FB_2",
                  "FB_3",
                  "FB_4"
                  
)


factor(Idents(FB_subset_2), levels= FB2myLevels)
Idents(FB_subset_2) <- factor(Idents(FB_subset_2), levels= FB2myLevels)


#DotPlot(EC_subest_2, col.min = 0, features = c(features_EC,features_NOS3,features_WNT4,features_VCAM1,features_LHX6)) + RotatedAxis()
DotPlot(FB_subset_2, col.min = 0, features = c(features_4,features_3,features_2,features_1)) + RotatedAxis() & coord_flip()
############

saveRDS(FB_subset, "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/FB_SUB_ALL.rds")


######################### FIG 3 #######################
######################### FIG 3 #######################
######################### FIG 3 #######################

CM_subset <- subset(aa1,  subtypes %in% c("Ventricular_CM_Cav3_1+", "Ventricular_CM_ROBO2+ ", "Proliferating_CMs_NPPA+", "Proliferating_CM","Proliferating_CMs_NPPA+", "Ventricular_CM_injury_specific", "Prol_CM_3" ,"Ven_CM_(Cav3_1)"))

CM_subset <- FindVariableFeatures(CM_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CM_subset)
CM_subset <- ScaleData(CM_subset, features = all.genes)

CM_subset <- RunPCA(CM_subset, features = VariableFeatures(object = CM_subset))

CM_subset <- FindNeighbors(CM_subset, dims = 1:10)


CM_subset <- FindClusters(CM_subset, resolution = 0.3)

CM_subset <- RunUMAP(CM_subset, dims = 1:10)

CM_subset_marker <- FindAllMarkers(CM_subset, min.pct = 0.6, min.diff.pct = 0.3)


CM_subset$subtypes = Idents(CM_subset)


DimPlot(CM_subset, repel = TRUE) 
DimPlot(CM_subset,  repel = TRUE, split.by = "condition")
DimPlot(CM_subset,  repel = TRUE, split.by = "condition", group.by = "subtypes")
DimPlot(CM_subset,  repel = TRUE, group.by = "subtypes")

features = rownames(CM_subset)[grep('ADAMTS17', rownames(CM_subset))]
FeaturePlot(CM_subset, features = features, cols = c('gray', 'red'))


CM.cluster.ids <- c("CM_ven_(Robo2)", "CM_ven_(Cav3_1)", "CM_IS", "CM_Prol_1", "CM_Prol_2", "CM_Prol_3", "CM_Prol_IS")
names(CM.cluster.ids) <- levels(CM_subset)
CM_subset <- RenameIdents(CM_subset, CM.cluster.ids)

load("/groups/tanaka/People/current/Paco/Collaborations/Elad_Paco/CellCycleGenes.RData")

s.genes.new <- c()
for (i in s.genes) {
  s.genes.new <-c(s.genes.new, rownames(CM_subset)[grep(i, rownames(CM_subset))])
}

g2m.genes.new <- c()
for (i in g2m.genes) {
  g2m.genes.new <-c(g2m.genes.new, rownames(CM_subset)[grep(i, rownames(CM_subset))])
}

CM_subset <- CellCycleScoring(CM_subset, s.features = s.genes.new, g2m.features = g2m.genes.new, set.ident = TRUE)

CM_subset$seurat_clusters -> Idents(CM_subset)
CM_subset$subtypes <- Idents(CM_subset)

cell.sels = colnames(CM_subset)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes = as.character(aa$subtypes)
aa$subtypes[mm] <- as.character(CM_subset$subtypes)
aa$subtypes = as.factor(aa$subtypes)

saveRDS(CM_subset, "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/VEN_CM_subset.rds")


features_Robo2 = rownames(CM_subset)[grep('ROBO2-AMEX60DD047486|SCN5A-AMEX60DD021356|SCN2A-AMEX60DD021357', rownames(CM_subset))]
features_Cav = rownames(CM_subset)[grep('CACNA1G-AMEX60DD030985|AGBL1-AMEX60DD004208|SEPTIN9-AMEX60DD030744', rownames(CM_subset))]
features_Pro1 = rownames(CM_subset)[grep('ADAM8-AMEX60DD052283|AXL-AMEX60DD024484|RAI14-AMEX60DD041872', rownames(CM_subset))]
features_Pro2 = rownames(CM_subset)[grep('CENPF-AMEX60DD036088|ASPM-AMEX60DD018743|KIF2C-AMEX60DD019926', rownames(CM_subset))]
features_inj = rownames(CM_subset)[grep('ACAN-AMEX60DD007108|ADAMTS6-AMEX60DD042234|ARHGAP31-AMEX60DD047024', rownames(CM_subset))]
features_Pro3 = rownames(CM_subset)[grep('CLSPN-AMEX60DD005937|DIAPH3-AMEX60DD048907|CENPE-AMEX60DD044675', rownames(CM_subset))]


CMmyLevels <- c("Prol_CM_3",
                "Prol_CM_2",
                "Prol_CM_1_injury_specific",
                "CM_injury_specific",
                "Ventricular_CM_Cav3_1",
                "Ven_CM_Robo2"
                
)




factor(Idents(CM_subset), levels= CMmyLevels)
Idents(CM_subset) <- factor(Idents(CM_subset), levels= CMmyLevels)


DotPlot(CM_subset, col.min = 0, features = c( features_Robo2,features_Cav,features_inj,features_Pro1,features_Pro2,features_Pro3)) + RotatedAxis()
############

CM_subset_2 <- CM_subset

s.genes.new <- c()
for (i in s.genes) {
  s.genes.new <-c(s.genes.new, rownames(CM_subset_2)[grep(i, rownames(CM_subset_2))])
}

g2m.genes.new <- c()
for (i in g2m.genes) {
  g2m.genes.new <-c(g2m.genes.new, rownames(CM_subset_2)[grep(i, rownames(CM_subset_2))])
}
CM_subset_2 <- CellCycleScoring(CM_subset_2, s.features = s.genes.new, g2m.features = g2m.genes.new, set.ident = TRUE)

CM_subset_2 <- ScaleData(CM_subset_2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(CM_subset_2))

CM_subset_2 <- RunPCA(CM_subset_2, features = c(s.genes, g2m.genes))

CM.cluster.ids <- c("Ven_CM_Robo2", "Ventricular_CM_Cav3_1", "Prol_CM_1_injury_specific", "Prol_CM_2", "CM_injury_specific", "Prol_CM_3")
names(CM.cluster.ids) <- levels(CM_subset_2)
CM_subset_2 <- RenameIdents(CM_subset_2, CM.cluster.ids)


#CM_subset_2 <- SetIdent(CM_subset_2, value = CM_subset_2$seurat_clusters)

features_Robo2 = rownames(CM_subset_2)[grep('ROBO2-AMEX60DD047486|SCN5A-AMEX60DD021356|SCN2A-AMEX60DD021357', rownames(CM_subset_2))]
features_Cav = rownames(CM_subset_2)[grep('CACNA1G-AMEX60DD030985|AGBL1-AMEX60DD004208|SEPTIN9-AMEX60DD030744', rownames(CM_subset_2))]
features_Pro1 = rownames(CM_subset_2)[grep('ADAM8-AMEX60DD052283|AXL-AMEX60DD024484|RAI14-AMEX60DD041872', rownames(CM_subset_2))]
features_Pro2 = rownames(CM_subset_2)[grep('CENPF-AMEX60DD036088|ASPM-AMEX60DD018743|KIF2C-AMEX60DD019926', rownames(CM_subset_2))]
features_inj = rownames(CM_subset_2)[grep('ACAN-AMEX60DD007108|ADAMTS6-AMEX60DD042234|ARHGAP31-AMEX60DD047024', rownames(CM_subset_2))]
features_Pro3 = rownames(CM_subset_2)[grep('CLSPN-AMEX60DD005937|DIAPH3-AMEX60DD048907|CENPE-AMEX60DD044675', rownames(CM_subset_2))]


CMmyLevels <- c("Prol_CM_3",
                "Prol_CM_2",
                "Prol_CM_1_injury_specific",
                "CM_injury_specific",
                "Ventricular_CM_Cav3_1",
                "Ven_CM_Robo2"
                
)




factor(Idents(CM_subset_2), levels= CMmyLevels)
Idents(CM_subset_2) <- factor(Idents(CM_subset_2), levels= CMmyLevels)


DotPlot(CM_subset_2, col.min = 0, features = c( features_Robo2,features_Cav,features_inj,features_Pro1,features_Pro2,features_Pro3)) + RotatedAxis()
############
CM_subset_marker_after_reg <- FindAllMarkers(CM_subset_2, min.pct = 0.6, min.diff.pct = 0.3)

FeaturePlot(CM_subset_2, c("G2M.Score"), cols = c("#dee2e6", "#661CB0"))

CM_subset_2$CC.Difference <- CM_subset_2$S.Score - CM_subset_2$G2M.Score

CM_subset_2 <- ScaleData(CM_subset_2, vars.to.regress = "CC.Difference", features = rownames(CM_subset_2))

saveRDS(CM_subset_2, "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/VEN_CM_subset_after_regress_g2m.rds")

CM_subset_2 <- RunPCA(CM_subset_2, features = VariableFeatures(CM_subset_2), nfeatures.print = 10)
CM_subset_2 <- RunPCA(CM_subset_2, features = c(s.genes, g2m.genes))
DimPlot(CM_subset_2, reduction = "pca")

sumCM_subset_2 <- FindNeighbors(CM_subset_2, dims = 1:5)


CM_subset_2 <- FindClusters(CM_subset_2, resolution = 0.1)

CM_subset_2 <- RunUMAP(CM_subset_2, dims = 1:5)

DimPlot(CM_subset_2)


#########

CM_subset_3 <- subset(aa,  subtypes %in% c("Ventricular_CM_Cav3_1+", "Ventricular_CM_ROBO2+", "Proliferating_CM", "Ventricular_CM_injury_specific" ))

CM_subset_3 <- FindVariableFeatures(CM_subset_3, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(CM_subset_3)
CM_subset_3 <- ScaleData(CM_subset_3, features = all.genes)

CM_subset_3 <- RunPCA(CM_subset_3, features = VariableFeatures(object = CM_subset_3))

CM_subset_3 <- FindNeighbors(CM_subset_3, dims = 1:10)


CM_subset_3 <- FindClusters(CM_subset_3, resolution = 0.25)

CM_subset_3 <- RunUMAP(CM_subset_3, dims = 1:10)

CM.cluster.ids <- c("Ven_CM_Robo2", "Ventricular_CM_Cav3_1", "Prol_CM_1_injury_specific", "Prol_CM_2", "CM_injury_specific", "Prol_CM_3")
names(CM.cluster.ids) <- levels(CM_subset_3)
CM_subset_3 <- RenameIdents(CM_subset_3, CM.cluster.ids)

DimPlot(CM_subset_3, group.by = "newId")

CM_subset_3$newId -> CM_subset_2$newId
CM_subset_3$newId -> Idents(CM_subset_2)

DimPlot(CM_subset_2, group.by = "seurat_clusters")

load("/groups/tanaka/People/current/Paco/Collaborations/Elad_Paco/CellCycleGenes.RData")

s.genes.new <- c()
for (i in s.genes) {
  s.genes.new <-c(s.genes.new, rownames(CM_subset_2)[grep(i, rownames(CM_subset_2))])
}

g2m.genes.new <- c()
for (i in g2m.genes) {
  g2m.genes.new <-c(g2m.genes.new, rownames(CM_subset_2)[grep(i, rownames(CM_subset_2))])
}

CM_subset_2 <- CellCycleScoring(CM_subset_2, s.features = s.genes.new, g2m.features = g2m.genes.new, set.ident = TRUE)

CM_subset_2 <- ScaleData(CM_subset_2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(CM_subset_2))

CM_subset_2 <- RunPCA(CM_subset_2, features = VariableFeatures(CM_subset_2))
#####NEED TO USE .NEW!!!!!!!!!!

CM_subset_2 <- FindNeighbors(CM_subset_2, dims = 1:10)

CM_subset_2 <- FindClusters(CM_subset_2, resolution = 0.3)
CM_subset_2 <- RunUMAP(CM_subset_2, dims = 1:10)

DimPlot(CM_subset_2, split.by ="condition" )
DimPlot(CM_subset_2, reduction = "pca", group.by = "Phase")

DimPlot(CM_subset_2)
DimPlot(CM_subset_2, group.by = "newId")
DimPlot(CM_subset_2, group.by = "Phase")


CM_sub_marker <- FindAllMarkers(CM_subset_2, min.pct = 0.4, min.diff.pct = 0.2)


CM_subset_2 <- SetIdent(CM_subset_2, value = CM_subset_2$newId)

CM_subset_2$subtypes <- Idents(CM_subset_2)

#cell.sels = colnames(FB_subset)
##mm = match(cell.sels, colnames(aa))
#cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
#aa$subtypes[mm] <- as.character(FB_subset$subtypes)




cell.sels = colnames(CM_subset_2)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes_RNA = as.character(aa$subtypes_RNA)
aa$subtypes_RNA[mm] <- as.character(CM_subset_2$subtypes)
aa$subtypes_RNA = as.factor(aa$subtypes_RNA)

features = rownames(CM_subset_2)[grep('SLIT3', rownames(CM_subset_2))]
FeaturePlot(CM_subset_2, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))


#############
#MY_subset <- subset(aa,  subtypes %in% c("B_cells", "Neutrophil", "Proliferating_Mono_Macrophages", "Resident_MF","Mono_Macrophages", "Proliferating_B_cells", "T_cells" ))

MY_subset_2 <- subset(aa1,  subtypes %in% c("Mono_Macrophages", "Neutrophil","Proliferating_Mono_Macrophages", "Resident_MF" ))

#MY_subset <- subset(MY_subset,  subtypes %in% c("Mono/MF_(CD163L1+)", "B_cells_(MUSK+)", "B_Cells", "T_cells", "B_cells_(Prol)", "Mono/MF_(Prol)"))

MY_subset_2 <- FindVariableFeatures(MY_subset_2, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(MY_subset_2)
MY_subset_2 <- ScaleData(MY_subset_2, features = all.genes)

MY_subset_2 <- RunPCA(MY_subset_2, features = VariableFeatures(object = MY_subset_2))

MY_subset_2 <- FindNeighbors(MY_subset_2, dims = 1:10)


MY_subset_2 <- FindClusters(MY_subset_2, resolution = 0.3)

MY_subset_2 <- RunUMAP(MY_subset_2, dims = 1:10)

MY_subset_2_marker <- FindAllMarkers(MY_subset_2, min.pct = 0.4, min.diff.pct = 0.2)


features = rownames(MY_subset_2)[grep('ITGB3-AMEX60DD00975', rownames(MY_subset_2))]
FeaturePlot(MY_subset_2, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))

DimPlot(MY_subset_2)
DimPlot(MY_subset_2, group.by = "subtypes")
DimPlot(MY_subset_2, split.by = "condition")
DimPlot(MY_subset_2, split.by = "condition", group.by = "subtypes")

##MY.cluster.ids <- c("Mono/MF_(CD163L1+)", "B_cells_(MUSK+)", "Mono/MF_(PLBD1+)", "B_Cells", "T_cells", "B_cells_(Prol)", "T_cells", "Mono/MF_(Prol)")
#names(MY.cluster.ids) <- levels(MY_subset)
#MY_subset <- RenameIdents(MY_subset, MY.cluster.ids)

#MY_subset$subtypes -> MY_subset_2$newId_1


#####


features = rownames(MY_subset_2)[grep('IL1R1', rownames(MY_subset_2))]
FeaturePlot(MY_subset_2, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))

DimPlot(MY_subset_2, label = TRUE, group.by = "newId_1")+NoLegend()

DimPlot(MY_subset_2)
DimPlot(MY_subset_2, split.by = "condition")

MY.cluster.ids <- c("Mo/Macs_(FAXDC2)","Neu_(DYSF)","Mo/Macs_(SNX22)" ,"Mo/Macs_resident","Mo/Macs_Prol","doublet","Neu_(IL1R1)")
names(MY.cluster.ids) <- levels(MY_subset_2)
MY_subset_2 <- RenameIdents(MY_subset_2, MY.cluster.ids)

MY_subset_2$subtypes <- Idents(MY_subset_2)

cell.sels = colnames(MY_subset_2)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes = as.character(aa$subtypes)
aa$subtypes[mm] <- as.character(MY_subset_2$subtypes)
aa$subtypes = as.factor(aa$subtypes)

#####

FB_subset <- subset(aa,  subtypes %in% c("FB_1", "FB_2", "FB_3", "FB_4"))

FB_subset <- FindVariableFeatures(FB_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(FB_subset)
FB_subset <- ScaleData(FB_subset, features = all.genes)

FB_subset <- RunPCA(FB_subset, features = VariableFeatures(object = FB_subset))


FB_subset <- FindNeighbors(FB_subset, dims = 1:10)


FB_subset <- FindClusters(FB_subset, resolution = 0.2)


FB_subset <- RunUMAP(FB_subset, dims = 1:10)

FB_subset_marker_1 <- FindAllMarkers(FB_subset, min.pct = 0.6, min.diff.pct = 0.3)

DimPlot(FB_subset, label = TRUE, group.by = "subtypes")+NoLegend()

DimPlot(FB_subset)
DimPlot(FB_subset, split.by = "condition")

features = rownames(EC_subset)[grep('AMEX60DD008752', rownames(EC_subset))]
FeaturePlot(EC_subset, split.by = "condition", features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))


FB.cluster.ids <- c("FB_(PKD1)", "FB_(TNXB)", "FB_(VWA2)", "FB_IS_(TFPI2)", "FB_IS_(TNC)","FB_Prol", "doublet")
names(FB.cluster.ids) <- levels(FB_subset)
FB_subset <- RenameIdents(FB_subset, FB.cluster.ids)

Idents(FB_subset) -> FB_subset$subtypes

cell.sels = colnames(FB_subset)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes = as.character(aa$subtypes)
aa$subtypes[mm] <- as.character(FB_subset$subtypes)
aa$subtypes = as.factor(aa$subtypes)

features = rownames(FB_subset)[grep('COL1A2', rownames(FB_subset))]
FeaturePlot(FB_subset, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))


features = rownames(aa)[grep('BMP5|BMP10', rownames(aa))]
FeaturePlot(aa, features = features, split.by = "condition",  order = TRUE, cols = c("#dee2e6", "#661CB0"))

#AMEX60DD008752 = HMGA1?

EC_subset <- subset(aa1,  subtypes %in% c("EC", "EC_IS_(LOX)", "EC_(CEMIP)", "EC_(NOS3)", "EC_(WNT4)", "EC_Prol","EC_IS_(IARS1)","EC_IS_Prol","EC_(LHX6)"))

EC_subset <- FindVariableFeatures(EC_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(EC_subset)
EC_subset <- ScaleData(EC_subset, features = all.genes)

EC_subset <- RunPCA(EC_subset, features = VariableFeatures(object = EC_subset))

EC_subset <- FindNeighbors(EC_subset, dims = 1:10)


EC_subset <- FindClusters(EC_subset, resolution = 0.25)


EC_subset <- RunUMAP(EC_subset, dims = 1:10)

EC_subset$subtypes -> Idents(EC_subset)

DimPlot(EC_subset)
DimPlot(EC_subset, group.by = "subtypes")
DimPlot(EC_subset,  split.by = "condition")
DimPlot(EC_subset, group.by = "seurat_clusters")

EC_subset_marker <- FindAllMarkers(EC_subset, min.pct = 0.4, min.diff.pct = 0.2) 


features = rownames(aa)[grep('CSPG4', rownames(aa))]
FeaturePlot(aa,  features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))

load("/groups/tanaka/People/current/Paco/Collaborations/Elad_Paco/CellCycleGenes.RData")

s.genes.new <- c()
for (i in s.genes) {
  s.genes.new <-c(s.genes.new, rownames(EC_subset)[grep(i, rownames(EC_subset))])
}

g2m.genes.new <- c()
for (i in g2m.genes) {
  g2m.genes.new <-c(g2m.genes.new, rownames(EC_subset)[grep(i, rownames(EC_subset))])
}

EC_subset <- CellCycleScoring(EC_subset, s.features = s.genes.new, g2m.features = g2m.genes.new, set.ident = TRUE)


FeaturePlot(EC_subset, c("G2M.Score", "S.Score"), cols = c("#dee2e6", "#661CB0"))

Idents(EC_subset) <- EC_subset$seurat_clusters

EC.cluster.ids <- c("EC", "EC_Prol","EC_IS_(LOX)","EC_(NOS3)",  "EC_(CEMIP)", "EC_(WNT4)","EC_IS_(IARS1)","EC_IS_Prol","EC_(LHX6)")
names(EC.cluster.ids) <- levels(EC_subset)
EC_subset <- RenameIdents(EC_subset, EC.cluster.ids)


EC_subset$subtypes <- Idents(EC_subset)

#cell.sels = colnames(FB_subset)
##mm = match(cell.sels, colnames(aa))
#cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')


cell.sels = colnames(EC_subset)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes = as.character(aa$subtypes)
aa$subtypes[mm] <- as.character(EC_subset$subtypes)
aa$subtypes = as.factor(aa$subtypes)

saveRDS(EC_subset, "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/EC_subset_final_20221117.rds")

features_EC = rownames(EC_subset)[grep('SLCO2A1-AMEX60DD002307|KCNMB2-AMEX60DD030396|FMO2-AMEX60DD018379', rownames(CM_subset))]
features_EC_Prol = rownames(EC_subset)[grep('CENPF-AMEX60DD036088|ASPM-AMEX60DD018743|BRIP1-AMEX60DD054637', rownames(CM_subset))]
features_EC_IS_LOX = rownames(EC_subset)[grep('LOX-AMEX60DD042852|CILP-AMEX60DD003432|P2RY6-AMEX60DD049766', rownames(CM_subset))]
features_NOS3 = rownames(EC_subset)[grep('NOS3-AMEX60DD026163|ADAM15-AMEX60DD014884|HSPA12B-AMEX60DD046249', rownames(CM_subset))]
features_WNT4 = rownames(EC_subset)[grep('KCNQ5-AMEX60DD033595|HEG1-AMEX60DD009667|WNT4-AMEX60DD052091', rownames(CM_subset))]
features_CEMIP = rownames(EC_subset)[grep('CEMIP-AMEX60DD00429|VCAM1-AMEX60DD018837|PTH1R-AMEX60DD020537', rownames(CM_subset))]
features_LHX6 = rownames(EC_subset)[grep('PDE2A-AMEX60DD049847|UNC5B-AMEX60DD051903|LHX6-AMEX60DD050778', rownames(CM_subset))]
features_EC_IS_IARS1 = rownames(EC_subset)[grep('ACAN-AMEX60DD007108|SULF1-AMEX60DD03987|IARS1-AMEX60DD023979', rownames(CM_subset))]
features_EC_IS_Prol = rownames(EC_subset)[grep('CCN1-AMEX60DD051003|MICAL2-AMEX60DD004827|EXT1-AMEX60DD040379', rownames(CM_subset))]

EC.cluster.ids <- c("EC", "EC_Prol","EC_IS_(LOX)","EC_(NOS3)",  "EC_(CEMIP)", "EC_(WNT4)","EC_IS_(IARS1)","EC_IS_Prol","EC_(LHX6)")


ECmyLevels <- c("EC", "EC_(NOS3)","EC_(WNT4)","EC_(LHX6)",  "EC_(CEMIP)","EC_Prol", "EC_IS_(LOX)","EC_IS_(IARS1)","EC_IS_Prol")



factor(Idents(EC_subset), levels= ECmyLevels)
Idents(EC_subset) <- factor(Idents(EC_subset), levels= ECmyLevels)

dimp
#DotPlot(EC_subest_2, col.min = 0, features = c(features_EC,features_NOS3,features_WNT4,features_VCAM1,features_LHX6)) + RotatedAxis()
DotPlot(EC_subset, col.min = 0, features = c(features_EC_Prol,features_EC_IS_Prol,features_EC_IS_IARS1,features_EC_IS_LOX,features_CEMIP,features_LHX6,features_WNT4,features_NOS3,features_EC)) + RotatedAxis() & coord_flip()
############

########

LY_subset <- subset(aa,  subtypes %in% c("B_cells", "Proliferating_B_cells","T_cells"))


LY_subset <- FindVariableFeatures(LY_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(LY_subset)
LY_subset <- ScaleData(LY_subset, features = all.genes)

LY_subset <- RunPCA(LY_subset, features = VariableFeatures(object = LY_subset))

LY_subset <- FindNeighbors(LY_subset, dims = 1:10)


LY_subset <- FindClusters(LY_subset, resolution = 0.15)

LY_subset <- RunUMAP(LY_subset, dims = 1:10)

LY_subset_2_marker <- FindAllMarkers(LY_subset, min.pct = 0.4, min.diff.pct = 0.2)


features = rownames(aa)[grep('LPAR1', rownames(aa))]
FeaturePlot(aa, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))

DimPlot(LY_subset)
DimPlot(LY_subset, group.by = "subtypes")
DimPlot(LY_subset, split.by = "condition")
DimPlot(LY_subset, split.by = "condition", group.by = "subtypes")

LY.cluster.ids <- c("B_cells_(FOXO1)", "B_cells_(SIGLEC11)","T_cells","B_cells_Prol",  "T_cells", "doublet")
names(LY.cluster.ids) <- levels(LY_subset)
LY_subset <- RenameIdents(LY_subset, LY.cluster.ids)

LY_subset$subtypes <- Idents(LY_subset)

cell.sels = colnames(LY_subset)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes = as.character(aa$subtypes)
aa$subtypes[mm] <- as.character(LY_subset$subtypes)
aa$subtypes = as.factor(aa$subtypes)
######
Blo_subset <- subset(aa,  subtypes %in% c("Megakeryocytes", "Proliferating_RBC","RBC", "Proliferating_Megakeryocytes"))


Blo_subset <- FindVariableFeatures(Blo_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(Blo_subset)
Blo_subset <- ScaleData(Blo_subset, features = all.genes)

Blo_subset <- RunPCA(Blo_subset, features = VariableFeatures(object = Blo_subset))

Blo_subset <- FindNeighbors(Blo_subset, dims = 1:10)


Blo_subset <- FindClusters(Blo_subset, resolution = 0.3)

Blo_subset <- RunUMAP(Blo_subset, dims = 1:10)

Blo_subset_2_marker <- FindAllMarkers(Blo_subset, min.pct = 0.4, min.diff.pct = 0.2)


features = rownames(aa)[grep('MYH11', rownames(aa))]
FeaturePlot(aa, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))

DimPlot(Blo_subset)
DimPlot(Blo_subset, group.by = "subtypes")
DimPlot(Blo_subset, split.by = "condition")
DimPlot(Blo_subset, split.by = "condition", group.by = "subtypes")

Blo.cluster.ids <- c("B_cells_(FOXO1)", "B_cells_(SIGLEC11)","T_cells","B_cells_Prol",  "T_cells", "doublet")
names(Blo.cluster.ids) <- levels(Blo_subset)
Blo_subset <- RenameIdents(Blo_subset, Blo.cluster.ids)

Blo_subset$subtypes <- Idents(Blo_subset)

cell.sels = colnames(Blo_subset)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes = as.character(aa$subtypes)
aa$subtypes[mm] <- as.character(Blo_subset$subtypes)
aa$subtypes = as.factor(aa$subtypes)


aa <- RenameIdents(object = aa, `RBC_Prol` = "RBC_Prol")

aa$subtypes -> Idents(aa)

#######
CM_subset <- subset(aa,  subtypes %in%  c(  "CM_PM_(HCN4)", "CM_OFT","CM_Atria","CM_Atria_Tagln","CM_Atria_Ltbp2","CM_ven_(Robo2)", "CM_ven_(Cav3_1)", "CM_IS", "CM_Prol_1", "CM_Prol_2", "CM_Prol_3", "CM_Prol_IS"))
DimPlot(CM_subset)
CM_subset$subtypes -> Idents(CM_subset)

CM_subset <- FindVariableFeatures(CM_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CM_subset)
CM_subset <- ScaleData(CM_subset, features = all.genes)

CM_subset <- RunPCA(CM_subset, features = VariableFeatures(object = CM_subset))

CM_subset <- FindNeighbors(CM_subset, dims = 1:10)


CM_subset <- FindClusters(CM_subset, resolution = 0.3)

CM_subset <- RunUMAP(CM_subset, dims = 1:10)

CM_subset_marker <- FindAllMarkers(CM_subset, min.pct = 0.4, min.diff.pct = 0.2)

load("/groups/tanaka/People/current/Paco/Collaborations/Elad_Paco/CellCycleGenes.RData")

s.genes.new <- c()
for (i in s.genes) {
  s.genes.new <-c(s.genes.new, rownames(CM_subset)[grep(i, rownames(CM_subset))])
}

g2m.genes.new <- c()
for (i in g2m.genes) {
  g2m.genes.new <-c(g2m.genes.new, rownames(CM_subset)[grep(i, rownames(CM_subset))])
}

CM_subset <- CellCycleScoring(CM_subset, s.features = s.genes.new, g2m.features = g2m.genes.new, set.ident = TRUE)

CM_subset <- ScaleData(CM_subset, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(CM_subset))

CM_subset <- RunPCA(CM_subset, features = VariableFeatures(object = CM_subset))

CM_subset$subtypes -> Idents(CM_subset)

features = rownames(CM_subset)[grep('THSD4', rownames(CM_subset))]
FeaturePlot(CM_subset, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))

CM_Prol_IS.markers <- FindMarkers(CM_subset, ident.1 = "CM_Prol_IS", ident.2 = c("CM_Prol_1", "CM_Prol_2", "CM_Prol_3"), min.pct = 0.4, min.diff.pct = 0.2)
CM_IS.markers <- FindMarkers(CM_subset, ident.1 = "CM_IS", ident.2 = c("CM_ven_(Robo2)", "CM_ven_(Cav3_1)"), min.pct = 0.4, min.diff.pct = 0.2)

write.csv(CM_Prol_IS.markers, "/groups/tanaka/People/current/Elad/CM_Prol_IS.markers.csv", quote = F)
write.csv(CM_IS.markers, "/groups/tanaka/People/current/Elad/CM_IS.markers.csv", quote = F)

######

aa1 <- subset(aa,  subtypes %in%  c( 

"B_cells_(FOXO1)"   ,       "B_cells_(SIGLEC11)"   ,             "B_cells_Prol"     ,                "CM_Atria" ,

"CM_Atria_Tagln"   ,                     "CM_IS"          ,             "CM_OFT"         ,        "CM_PM_(HCN4)" ,

"CM_Prol_1"       ,             "CM_Prol_2"      ,              "CM_Prol_3"           ,        "CM_Prol_IS" ,

"CM_ven_(Cav3_1)"        ,       "CM_ven_(Robo2)"       ,                                      "EC" ,

"EC_(CEMIP)"      ,              "EC_(LHX6)"    ,                "EC_(NOS3)" ,                   "EC_(WNT4)" ,

"EC_IS_(IARS1)"    ,              "EC_IS_(LOX)"   ,                "EC_IS_Prol"  ,                    "EC_Prol" ,

"FB_(PKD1)"  ,                  "FB_(TNXB)"  ,                  "FB_(VWA2)" ,               "FB_IS_(TFPI2)" ,

"FB_IS_(TNC)"  ,                    "FB_Prol",               "Megakeryocytes"  ,           "Mo/Macs_(FAXDC2)" ,

"Mo/Macs_(SNX22)" ,                "Mo/Macs_Prol" ,            "Mo/Macs_resident"  ,                 "Neu_(DYSF)" ,

"Neu_(IL1R1)" ,                    "Neuronal", "Proliferating_Megakeryocytes",            "Proliferating_RBC",
 
"RBC"  ,                    "T_cells" ))

aa1 <- FindVariableFeatures(aa1, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(aa1)
aa1 <- ScaleData(aa1, features = all.genes)

aa1 <- RunPCA(aa1, features = VariableFeatures(object = aa1))

aa1 <- FindNeighbors(aa1, dims = 1:25)


#aa1 <- FindClusters(aa1, resolution = 0.3)

aa1 <- RunUMAP(aa1, dims = 1:25)

########

CMmyLevels <- c("CM_Atria",
                "CM_PM_(HCN4)",
                "CM_OFT",
                "CM_ven_(Robo2)",
                "CM_ven_(Cav3_1)",
                "CM_Prol_1",
                "CM_Prol_2",
                "CM_Prol_3",
                "CM_Atria_Tagln",
                "CM_IS",
                "CM_Prol_IS"
                
)

####
features_aCM = rownames(CM_subset)[grep('AGBL1-AMEX60DD004208|ADAMTS17-AMEX60DD003517|GFRA1-AMEX60DD053205', rownames(CM_subset))]
features_PM = rownames(CM_subset)[grep('LTBP2-AMEX60DD011463|HCN4-AMEX60DD003883|TENM2-AMEX60DD030344', rownames(CM_subset))]
features_OFT = rownames(CM_subset)[grep('LHCGR-AMEX60DD035890|FUT2-AMEX60DD017238|ADRA1B-AMEX60DD03002', rownames(CM_subset))]
features_ROBO = rownames(CM_subset)[grep('ROBO2-AMEX60DD047486|JPH3-AMEX60DD01717|STAC-AMEX60DD03796', rownames(CM_subset))]
features_CAV = rownames(CM_subset)[grep('CNTFR-AMEX60DD041547|ERV3-1-AMEX60DD016594|RBPMS-AMEX60DD043590', rownames(CM_subset))]
features_PR1 = rownames(CM_subset)[grep('CENPF-AMEX60DD036088|ASPM-AMEX60DD01874|ANLN-AMEX60DD03806', rownames(CM_subset))]
features_PR2 = rownames(CM_subset)[grep('SMC4-AMEX60DD001810|BUB1-AMEX60DD034434|IQGAP3-AMEX60DD01682', rownames(CM_subset))]
features_PR3 = rownames(CM_subset)[grep('CENPE-AMEX60DD044675|RTTN-AMEX60DD039311|DIAPH3-AMEX60DD048907', rownames(CM_subset))]
features_TAGLN = rownames(CM_subset)[grep('NPPA-AMEX60DD051099|PDCD6IP-AMEX60DD037944|TAGLN-AMEX60DD053922', rownames(CM_subset))]
features_IS = rownames(CM_subset)[grep('ACAN-AMEX60DD007108|ARHGAP31-AMEX60DD047024|NACC2-AMEX60DD050654', rownames(CM_subset))]
features_IS_PR = rownames(CM_subset)[grep('ANKRD1-AMEX60DD052557|CCN1-AMEX60DD051003|TLL1-AMEX60DD045228', rownames(CM_subset))]


CMmyLevels <- c("CM_Atria",
                "CM_PM_(HCN4)",
                "CM_OFT",
                "CM_ven_(Robo2)",
                "CM_ven_(Cav3_1)",
                "CM_Prol_1",
                "CM_Prol_2",
                "CM_Prol_3",
                "CM_Atria_Tagln",
                "CM_IS",
                "CM_Prol_IS"
                
)




factor(Idents(CM_subset), levels= CMmyLevels)
Idents(CM_subset) <- factor(Idents(CM_subset), levels= CMmyLevels)


DotPlot(CM_subset, col.min = 0, features = c( features_aCM,
                                                features_PM,
                                                features_OFT,
                                                features_ROBO,
                                                features_CAV,
                                                features_PR1,
                                                features_PR2,
                                                features_PR3,
                                                features_TAGLN,
                                                features_IS,
                                                features_IS_PR
                                                )) + RotatedAxis()
############
####
####
features_aCM = rownames(CM_subset)[grep('AGBL1-AMEX60DD004208|ADAMTS17-AMEX60DD003517|GFRA1-AMEX60DD053205', rownames(CM_subset))]
features_PM = rownames(CM_subset)[grep('LTBP2-AMEX60DD011463|HCN4-AMEX60DD003883|TENM2-AMEX60DD030344', rownames(CM_subset))]
features_OFT = rownames(CM_subset)[grep('LHCGR-AMEX60DD035890|FUT2-AMEX60DD017238|ADRA1B-AMEX60DD03002', rownames(CM_subset))]
features_ROBO = rownames(CM_subset)[grep('ROBO2-AMEX60DD047486|JPH3-AMEX60DD01717|STAC-AMEX60DD03796', rownames(CM_subset))]
features_CAV = rownames(CM_subset)[grep('CNTFR-AMEX60DD041547|ERV3-1-AMEX60DD016594|RBPMS-AMEX60DD043590', rownames(CM_subset))]
features_TAGLN = rownames(CM_subset)[grep('NPPA-AMEX60DD051099|PDCD6IP-AMEX60DD037944|TAGLN-AMEX60DD053922', rownames(CM_subset))]
features_IS = rownames(CM_subset)[grep('ACAN-AMEX60DD007108|ARHGAP31-AMEX60DD047024|NACC2-AMEX60DD050654', rownames(CM_subset))]
features_IS_PR = rownames(CM_subset)[grep('ANKRD1-AMEX60DD052557|CCN1-AMEX60DD051003|TLL1-AMEX60DD045228', rownames(CM_subset))]


CMmyLevels <- c("CM_Atria",
                "CM_Atria_Tagln",
                "CM_PM_(HCN4)",
                "CM_OFT",
                "CM_ven_(Robo2)",
                "CM_ven_(Cav3_1)",
                "CM_IS",
                "CM_Prol_IS"
)


CMmyLevels <- c("CM_Atria",
                "CM_PM_(HCN4)",
                "CM_OFT",
                "CM_ven_(Robo2)",
                "CM_ven_(Cav3_1)",
                "CM_Prol_1",
                "CM_Prol_2",
                "CM_Prol_3",
                "CM_Atria_Tagln",
                "CM_IS",
                "CM_Prol_IS"
                
)

factor(Idents(CM_subset), levels= CMmyLevels)
Idents(CM_subset) <- factor(Idents(CM_subset), levels= CMmyLevels)


DotPlot(CM_subset, col.min = 0, features = c(features_aCM,features_PM,
                                              features_OFT,
                                              features_ROBO,
                                              features_CAV,

                                              features_TAGLN,
                                              features_IS,
                                              features_IS_PR
)) + RotatedAxis()
############
####

factor(Idents(CM_subset_reg), levels= CMmyLevels)
Idents(CM_subset_reg) <- factor(Idents(CM_subset_reg), levels= CMmyLevels)

Idents(CM_subset) <- CM_subset$subtypes


DimPlot(CM_subset, split.by = "condition",  cols = c(
  
  
  "#4CC9F0",
  "#4AAFF0",
  "#49A2F0",
  "#4895EF",
  "#4361EE",
  "#414CDC",
  "#4042D3",
  "#3F37C9",
  "#941F56",
  "#AD1833",
  "#C61010"   
  ) )


#######
CM_subset_hs <- subset(aa,  subtypes %in%  c(  "CM_PM_(HCN4)", "CM_OFT","CM_Atria","CM_ven_(Robo2)", "CM_ven_(Cav3_1)"))
DimPlot(CM_subset_hs)
CM_subset_hs$subtypes -> Idents(CM_subset_hs)

CM_subset_hs <- FindVariableFeatures(CM_subset_hs, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CM_subset_hs)
CM_subset_hs <- ScaleData(CM_subset_hs, features = all.genes)

CM_subset_hs <- RunPCA(CM_subset_hs, features = VariableFeatures(object = CM_subset_hs))

CM_subset_hs <- FindNeighbors(CM_subset_hs, dims = 1:10)


CM_subset_hs <- FindClusters(CM_subset_hs, resolution = 0.3)

CM_subset_hs <- RunUMAP(CM_subset_hs, dims = 1:10)

CM_subset_hs_marker <- FindAllMarkers(CM_subset_hs, min.pct = 0.4, min.diff.pct = 0.2)

DimPlot(CM_subset_hs,  cols = c(
  
  
  "#4CC9F0",
  
  "#49A2F0",
  
  "#4361EE",
  
  
 
  "#941F56",

  "#C61010"   
) )

######################################### FB_HS

FB_subset_hs <- subset(aa,  subtypes %in%  c(  "FB_(PKD1)","FB_(TNXB)", "FB_(VWA2)"))
DimPlot(FB_subset_hs)
FB_subset_hs$subtypes -> Idents(FB_subset_hs)

FB_subset_hs <- FindVariableFeatures(FB_subset_hs, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(FB_subset_hs)
FB_subset_hs <- ScaleData(FB_subset_hs, features = all.genes)

FB_subset_hs <- RunPCA(FB_subset_hs, features = VariableFeatures(object = FB_subset_hs))

FB_subset_hs <- FindNeighbors(FB_subset_hs, dims = 1:10)


FB_subset_hs <- FindClusters(FB_subset_hs, resolution = 0.3)

FB_subset_hs <- RunUMAP(FB_subset_hs, dims = 1:10)
FB_subset_hs$subtypes -> Idents(FB_subset_hs)

FB_subset_hs_marker <- FindAllMarkers(FB_subset_hs, min.pct = 0.4, min.diff.pct = 0.2)

DimPlot(FB_subset_hs,  cols = c(
  
  
  "#4CC9F0",
  
  
  
  "#4361EE",
  
  
  
  "#941F56"
  
   
) )

features_1 = rownames(FB_subset_hs)[grep('EFEMP1-AMEX60DD035438|PKD1-AMEX60DD030949|CFB-AMEX60DD010523', rownames(FB_subset_hs))]
features_2 = rownames(FB_subset_hs)[grep('TNXB-AMEX60DD010491|LAMA2-AMEX60DD034687|NGF-AMEX60DD008715', rownames(FB_subset_hs))]
features_3 = rownames(FB_subset_hs)[grep('ACAN-AMEX60DD004178|VWA2-AMEX60DD053181|CSPG5-AMEX60DD026060', rownames(FB_subset_hs))]

FBmyLevels <- c( "FB_(PKD1)",
                 "FB_(TNXB)",
                 "FB_(VWA2)"
               
                 
)




factor(Idents(FB_subset_hs), levels= FBmyLevels)
Idents(FB_subset_hs) <- factor(Idents(FB_subset_hs), levels= FBmyLevels)


#DotPlot(EC_subest_2, col.min = 0, features = c(features_EC,features_NOS3,features_WNT4,features_VCAM1,features_LHX6)) + RotatedAxis()
DotPlot(FB_subset_hs, col.min = 0, features = c(features_3,features_2,features_1)) + RotatedAxis() & coord_flip()
############
######################################### EC_HS

EC_subset_hs <- subset(aa,  subtypes %in%  c(  "EC_(CEMIP)","EC", "EC_(LHX6)", "EC_(NOS3)", "EC_(WNT4)"))
DimPlot(EC_subset_hs)
EC_subset_hs$subtypes -> Idents(EC_subset_hs)

EC_subset_hs <- FindVariableFeatures(EC_subset_hs, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(EC_subset_hs)
EC_subset_hs <- ScaleData(EC_subset_hs, features = all.genes)

EC_subset_hs <- RunPCA(EC_subset_hs, features = VariableFeatures(object = EC_subset_hs))

EC_subset_hs <- FindNeighbors(EC_subset_hs, dims = 1:10)


EC_subset_hs <- FindClusters(EC_subset_hs, resolution = 0.3)

EC_subset_hs <- RunUMAP(EC_subset_hs, dims = 1:10)
EC_subset_hs$subtypes -> Idents(EC_subset_hs)

EC_subset_hs_marker <- FindAllMarkers(EC_subset_hs, min.pct = 0.4, min.diff.pct = 0.2)

DimPlot(EC_subset_hs,  cols = c(
  
  
  "#4CC9F0",
 
  "#49A2F0",
  
  
 
  "#4042D3",
 
  "#941F56",
  
  "#C61010" 
  
  
) )

features_EC = rownames(EC_subest_2)[grep('LAMA2-AMEX60DD034687|KCNMB2-AMEX60DD030396|RADIL-AMEX60DD025680', rownames(CM_subset))]
features_NOS3 = rownames(EC_subest_2)[grep('NOS3-AMEX60DD026163|ADAM15-AMEX60DD014884|HSPA12B-AMEX60DD046249', rownames(CM_subset))]
features_WNT4 = rownames(EC_subest_2)[grep('KCNQ5-AMEX60DD033595|HEG1-AMEX60DD009667|WNT4-AMEX60DD052091', rownames(CM_subset))]
features_VCAM1 = rownames(EC_subest_2)[grep('NTRK3-AMEX60DD004199|VCAM1-AMEX60DD018837|CA8-AMEX60DD039785', rownames(CM_subset))]
features_LHX6 = rownames(EC_subest_2)[grep('PDE2A-AMEX60DD049847|ADGRF5-AMEX60DD032958|LHX6-AMEX60DD050778', rownames(CM_subset))]

FBmyLevels <- c( "FB_(PKD1)",
                 "FB_(TNXB)",
                 "FB_(VWA2)"
                 
                 
)




factor(Idents(EC_subset_hs), levels= FBmyLevels)
Idents(EC_subset_hs) <- factor(Idents(EC_subset_hs), levels= FBmyLevels)


#DotPlot(EC_subest_2, col.min = 0, features = c(features_EC,features_NOS3,features_WNT4,features_VCAM1,features_LHX6)) + RotatedAxis()
DotPlot(EC_subset_hs, col.min = 0, features = c(features_3,features_2,features_1)) + RotatedAxis() & coord_flip()
############

##### color scheme

###Basic

"#4CC9F0",
"#4AAFF0",
"#49A2F0",
"#4895EF",
"#4361EE",
"#414CDC",
"#4042D3",
"#3F37C9",
"#941F56",
"#AD1833",
"#C61010"  

#Blues
"#E3F2FD",
"#BBDEFB",
"#90CAF9",
"#64B5F6",
"#42A5F5",
"#2196F3",
"#1E88E5",
"#1976D2",
"1565C0",
"0D47A1"

#Purples
"#DC97FF",
"#D283FF",
"#BD68EE",
"#AB51E3",
"#8B2FC9",
"#6818A5",
"#5A108F",
"#4A0A77",
"3C0663",
"310055"

#Reds
"#E01E37",
"#C71F37",
"#B21E35",
"A11D33",
"85182A",
"641220"





######



saveRDS(aa1, "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_subtypes_final_20221117.rds")


