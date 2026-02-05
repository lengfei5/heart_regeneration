##########################################################################
##########################################################################
# Project: Heart regeneration 
# Script purpose: manually annotate snRNA-seq data of axolotl heart with multiple iterations
# Usage example: 
# Author: Elad Bassat (elad.bassat.imba.oeaw.ac.at) and Jingkui Wang (jingkui.wang@imba.oeaw.ac.at)
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

#### install necessary R packages
# install.packages("Polychrome")
#library("Polychrome")
#install.packages("pals")
#install.packages(RColorBrewer)
# build-in color palette
#Glasbey = glasbey.colors(32)
#swatch(Glasbey)
library(pals)
pal.bands(alphabet, alphabet2, cols25, glasbey, kelly, polychrome, 
          stepped, tol, watlington,
          show.names=FALSE)


library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

SeuratObj = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/Rdata_spliced/'
resDir = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/res_Elad'

setwd("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome")


########################################################
########################################################
# Section I: improt the snRNA-seq Seurat object processed by Jingkui 
# general check and further processing
########################################################
########################################################
aa = readRDS(file = paste0(SeuratObj, 
                           'seuratObject_axloltl_scRNAseq_R13591_20220720_lognormamlized_pca_umap.rds'))

DimPlot(aa, label = TRUE, repel = TRUE, cols = "RColorBrewer") + ggtitle("scNuc (multiome)")
DimPlot(aa, label = TRUE, repel = TRUE, split.by = "condition") + ggtitle("scNuc (multiome)")
DimPlot(aa, label = TRUE) + NoLegend()
DimPlot(aa, group.by = "subtypes")

features = rownames(aa)[grep('KAZALD', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))
VlnPlot(aa, features = features) + NoLegend()

#ggsave(filename = paste0(resDir, '/FeaturePlot_CD45.pdf'), width = 8, height = 6)
aa
features = rownames(aa)[grep('ITGA2B|MPL-AMEX60DD020059|ITGB3-AMEX60DD009754|PECAM|SELP', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('HBEGF', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))
VlnPlot(aa, features = features) + NoLegend()

features = rownames(aa)[grep('MYH6|ACTN2|NPPA|TNNT2|GATA4', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('CD68|CD8A|CD74|CSF1R|ITGAM', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

features = rownames(aa)[grep('TBX18', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

aa_var <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(aa), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(aa)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

cluster5.markers <- FindMarkers(aa, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

##gives back number of values in the table
features = rownames(aa)[grep('COL8A1', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

Cluster_5_subest <- subset(aa, idents = 5)
#Subseting cluster
#aa_small <- FindNeighbors(aa_small, dims = 1:10)
#aa_small <- FindClusters(aa_small, resolution = 0.5)

DotPlot(Cluster_5_subest, features = features)$data[,c("features.plot", "id","pct.exp")]
#returns pct of cells expressing gene
#pbmc <- RunUMAP(pbmc, dims = 1:10)


##########################################
# further filtering doublets using DbouletFinder
# for each sample
##########################################
library(DoubletFinder)
aa_gen = aa;

### day0
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

aa_day_0 <- doubletFinder_v3(aa_day_0, PCs = 1:30, pN = 0.25, 
                             pK = pK, nExp = nExp_poi.adj,  
                             reuse.pANN = FALSE, sct = FALSE)


### day1
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

aa_day_1 <- doubletFinder_v3(aa_day_1, PCs = 1:30, 
                             pN = 0.25, pK = pK, nExp = nExp_poi.adj,  reuse.pANN = FALSE, sct = FALSE)

#### day 4
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

aa_day_4 <- doubletFinder_v3(aa_day_4, PCs = 1:30, 
                             pN = 0.25, pK = pK, nExp = nExp_poi.adj,  reuse.pANN = FALSE, sct = FALSE)



#### day 7
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


aa_day_7 <- doubletFinder_v3(aa_day_7, PCs = 1:30, pN = 0.25, 
                             pK = pK, nExp = nExp_poi.adj,  reuse.pANN = FALSE, sct = FALSE)

##### day 14
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


aa_day_14 <- doubletFinder_v3(aa_day_14, PCs = 1:30, 
                              pN = 0.25, pK = pK, nExp = nExp_poi.adj,  reuse.pANN = FALSE, sct = FALSE)


### merge singlet cells 
aa_day_0_singlet = subset(aa_day_0, 
                          cells = colnames(aa_day_0)[which(aa_day_0$DF.classifications_0.25_0.02_773 
                                                           == "Singlet")] )

aa_day_1_singlet = subset(aa_day_1, 
                          cells = colnames(aa_day_1)[which(aa_day_1$DF.classifications_0.25_0.01_606 
                                                          == "Singlet")] )

aa_day_4_singlet = subset(aa_day_4, 
                          cells = colnames(aa_day_4)[which(aa_day_4$DF.classifications_0.25_0.005_511 
                                                          == "Singlet")] )

aa_day_7_singlet = subset(aa_day_7, 
                          cells = colnames(aa_day_7)[which(aa_day_7$DF.classifications_0.25_0.005_610 
                                                           == "Singlet")] )

aa_day_14_singlet = subset(aa_day_14, 
                          cells = colnames(aa_day_14)[which(aa_day_14$DF.classifications_0.25_0.005_543 
                                                           == "Singlet")] )


aa_merge <- Merge(aa_day_0_singlet, 
                  Y = c(aa_day_1_singlet, aa_day_4_singlet, aa_day_7_singlet, aa_day_14_singlet))

saveRDS(aa_merge, file =  'SeuratObject_axolotl_snRNAseq_R13591_doublet_removal.rds')


########################################################
########################################################
# Section II: manual annotation of major clusters and subclusters: iteration I
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
########################################################
########################################################
aa = readRDS(file = 'SeuratObject_axolotl_snRNAseq_R13591_doublet_removal.rds')

DimPlot(aa, label = TRUE, repel = TRUE) + ggtitle("scNuc (multiome)")

aa_marker <- FindAllMarkers(aa, min.pct = 0.6, min.diff.pct = 0.3)


##########################################
# subclustering EC clusters 0,4,6,2,8,20
##########################################
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
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes[mm] = as.character(EC_subset$newId)

##########################################
# annotate subclusters of CM (clusters 3,21,19,15,11,16,1,9,17)
##########################################
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

CM.cluster.ids <- c("CM_ventricle", "CM_atrial", "Doublets", "Proliferating_CMs", 
                    "CM_OFT", "Proliferating_CMs_NPPA+", "Doublets", "NPPA+_CMs","Doublets")
names(CM.cluster.ids) <- levels(CM_subset)
CM_subset <- RenameIdents(CM_subset, CM.cluster.ids)

CM_subset$newId = Idents(CM_subset)

#aa$subtypes = NA
cell.sels = colnames(CM_subset)
mm = match(cell.sels, colnames(aa))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa$subtypes[mm] = as.character(CM_subset$newId)

##########################################
# subclustering FB clusters 12,13
##########################################
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


features = rownames(aa)[grep('HAPLN1', rownames(aa))]
FeaturePlot(aa, features = features, cols = c('gray', 'red'))

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

##########################################
# Eryothrocytes - 
##########################################
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


##########################################
# subset cells with annotations and save the object
##########################################
aa_sub <- subset(aa, cells = names(aa$subtypes)[aa$subtypes %in% c("CMs_(ROBO2+)", "NPPA+_CMs", 
                                                                   "Proliferating_CMs",
                                                                   "CM_(TBX5-,ERBB4+,FGF14+)" , 
                                                                   "Proliferating_CMs_NPPA+", 
                                                                   "EC", "Injury_specific_EC",
                                                                   "EC_(VCAM1+)", 
                                                                   "Proliferating_ECs",
                                                                   "EC_(Wnt4+)", 
                                                                   "EC_(LHX6+)", "RBC", 
                                                                   "Proliferating_RBC", 
                                                                   "FB_1" , "FB_2","FB_3",
                                                                   "Resident_MF",
                                                                   "B_cells",
                                                                   "Megakeryocytes", 
                                                                   "Mono_Macrophages",  
                                                                   "Proliferating_Megakeryocytes", 
                                                                   "Neutrophil", "T_cells",
                                                                   "NK_cells?_and_doubluets",
                                                                   "Proliferating_Mono_Macrophages", 
                                                                   "Neuronal",
                                                                   "Proliferating_B_cells")])

aa_sub <- FindVariableFeatures(aa_sub, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(aa_sub)
aa_sub <- ScaleData(aa_sub, features = all.genes)

aa_sub <- RunPCA(aa_sub, features = VariableFeatures(object = aa_sub))

aa_sub <- FindNeighbors(aa_sub, dims = 1:25)
aa_sub <- FindClusters(aa_sub, resolution = 0.5)

aa_sub <- RunUMAP(aa_sub, dims = 1:25)

DimPlot(aa, group.by = "subtypes", label = TRUE, split.by = "condition") + NoLegend()
DimPlot(aa_sub, label = TRUE)

DimPlot(aa_sub, group.by = "subtypes", label = TRUE) + NoLegend()

features = rownames(aa_sub)[grep('NKX2-5', rownames(aa_sub))]
FeaturePlot(aa_sub, features = features ,cols = c('gray', 'red'))

#saveRDS(aa_sub, "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets.rds")

### EXPORTING A PDF WHICH IS EDITABLE IN AI
pdf("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_sub.pdf", width = 8)
DimPlot(aa_sub, group.by = "subtypes", label = TRUE) + NoLegend()
dev.off()

########################################################
########################################################
# Section III: manual annotation of subclusters: iteration II 
# with cell cycle scores and visiuam spatial information
########################################################
########################################################
aa_sub = readRDS(file = paste0("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/",
                               "aa_annotated_no_doublets.rds"))

##### loading Rdataset of G2M and S genes
load("/groups/tanaka/People/current/Paco/Collaborations/Elad_Paco/CellCycleGenes.RData")


VlnPlot(aa_sub, group.by = "orig.ident",  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
DimPlot(aa_sub, group.by = "seurat_clusters", label = TRUE) + NoLegend()

s.genes.new <- c()
for (i in s.genes) {
  s.genes.new <-c(s.genes.new, rownames(CM_subset)[grep(i, rownames(CM_subset))])
}

g2m.genes.new <- c()
for (i in g2m.genes) {
  g2m.genes.new <-c(g2m.genes.new, rownames(CM_subset)[grep(i, rownames(CM_subset))])
}


aa <- CellCycleScoring(aa, s.features = s.genes.new, g2m.features = g2m.genes.new, set.ident = TRUE)

FeaturePlot(aa, c("G2M.Score"))
# 
# DimPlot(aa, group.by = "subtypes")

RidgePlot(aa_sub, group.by = "Phase", features = c("TOP2A-AMEX60DD010148"), ncol = 2)


##########################################
# CM 
##########################################
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

##export table
write.csv(Prolif_CM_marker, "/groups/tanaka/People/current/Elad/Prolif_CM_marker.csv", quote = F)


features = rownames(Prolif_CM)[grep('RAI14', rownames(Prolif_CM))]

FeaturePlot(Prolif_CM, split.by = "condition", features = features, cols = c('gray', 'red'))
FeaturePlot(Prolif_CM, features = features, cols = c('gray', 'red'))

FeaturePlot(aa_sub, features = features, cols = c('gray', 'red'))

####

features = rownames(aa_sub2)[grep('PTPRC|PECAM|COL1A2|THBS1|ANK1-AMEX60DD002871|NCAM2|TNNI3K', 
                                  rownames(aa_sub2))]
DotPlot(aa_sub2, group.by = "subtypes", features = features) + RotatedAxis()


FeaturePlot(aa_sub, features = features ,cols = c('gray', 'red'))


aa_sub2 <- subset(aa_sub, cells = names(aa_sub$subtypes)[aa_sub$subtypes != "NK_cells?_and_doubluets"]) 
General_markers = FindAllMarkers(aa_sub2, min.pct = 0.6, min.diff.pct = 0.3)

aa_sub2 <- SetIdent(aa_sub2, value = aa_sub2$subtypes)
DimPlot(aa_sub2)


FeaturePlot(aa_sub2, split.by = "condition", c("G2M.Score"))

##########################################
# border zone CM 
##########################################
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
BZ_CM.cluster.ids <- c("Proliferating_CM1", "Proliferating_CM2", 
                       "Proliferating_CM_injury_transient" , 
                       "Proliferating_CM_injury", "Proliferating_CM3")
names(BZ_CM.cluster.ids) <- levels(BZ_CM)
BZ_CM <- RenameIdents(BZ_CM, BZ_CM.cluster.ids)
DimPlot(BZ_CM2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


BZ_CM$newId = Idents(BZ_CM)

#aa$subtypes = NA
cell.sels = colnames(BZ_CM)
mm = match(cell.sels, colnames(aa_sub2))
cat(length(which(is.na(cell.sels))), '--', length(mm), '\n')
aa_sub2$subtypes[mm] = as.character(BZ_CM$newId)


##########################################
# Ventricluar CM subsetting
##########################################
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
aa_sub2 <- RenameIdents(aa_sub2, 
                        'Proliferating_CM1' = "Proliferating_CM",
                        'Proliferating_CM2' = "Proliferating_CM",
                        'Proliferating_CM_injury_transient' = "Proliferating_CM",  
                        'Proliferating_CM_injury' = "Proliferating_CM2")
##version 2
aa_sub2 <- RenameIdents(aa_sub2, 
                        'Proliferating_CM2' = "Proliferating_CM", 
                        'Proliferating_CM3' = "Proliferating_CM" )

saveRDS(aa_sub2, file= paste0("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/",
                        "aa_annotated_no_doublets_20221004_2.rds"))

aa_annotated_no_doublets_20221004$subtypes <- Idents(aa_annotated_no_doublets_20221004)

##########################################
# Macrophage  subsetting
##########################################
MF_subset <- subset(aa, subtypes %in% c("Resident_MF", 
                                        "Mono_Macrophages",  "Proliferating_Mono_Macrophages"))


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

##########################################
# Eryothrocytes - 
##########################################
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


Ery_subset <- ScaleData(Ery_subset, vars.to.regress = c("S.Score", "G2M.Score"), 
                        features = rownames(Ery_subset))


##########################################
# EC
##########################################
EC_subset <- subset(aa,  subtypes %in% c("EC", "Injury_specific_EC", 
                                         "EC_(VCAM1+)", "Proliferating_ECs", "EC_(Wnt4+)", "EC_(LHX6+)"))

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


##########################################
# FB 
##########################################
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
DimPlot(FB_subset2, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "condition", 
        group.by = "subtypes")
DimPlot(FB_subset2, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "seurat_clusters")
DimPlot(FB_subset2, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "newId")


##
aa_sub <- subset(aa, cells = names(aa$subtypes)[aa$subtypes %in% c("Proliferating_CM" ,
                                                                   "Ventricular_CM_injury_specific",  
                                                                   "Ventricular_CM_ROBO2+", 
                                                                   "Ventricular_CM_Cav3_1+",
                                                                   "CM_Atria",
                                                                   "Proliferating_CMs",
                                                                   "CM_OFT", 
                                                                   "Proliferating_CMs_NPPA+", "EC",
                                                                   "Injury_specific_EC",
                                                                   "EC_(VCAM1+)", 
                                                                   "Proliferating_ECs",
                                                                   "EC_(Wnt4+)", 
                                                                   "EC_(LHX6+)", "RBC", 
                                                                   "Proliferating_RBC", 
                                                                   "FB_1" , "FB_2", "FB_3", "FB_4",
                                                                   "Resident_MF",
                                                                   "B_cells",
                                                                   "Megakeryocytes",
                                                                   "Mono_Macrophages",
                                                                   "Proliferating_Megakeryocytes" , 
                                                                   "Neutrophil", 
                                                                   "T_cells",
                                                                   "Proliferating_Mono_Macrophages" , 
                                                                   "Neuronal",
                                                                   "Proliferating_B_cells")])



saveRDS(aa_sub, file = paste0("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/",
                              "aa_annotated_no_doublets_2022_10_17.rds"))

