##########################################################################
##########################################################################
# Project: Heart regeneration in axolotl 
# Script purpose: plots for figure with snRNA-seq data
# Usage example: 
# Author: Elad Bassat (elad.bassat.imba.oeaw.ac.at) and Jingkui Wang (jingkui.wang@imba.oeaw.ac.at)
# Date of creation: Thu Feb  5 14:14:02 2026
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
# Section I: 
# 
########################################################
########################################################
aa = readRDS(paste0("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/",
                     "aa_annotated_no_doublets_2022_10_17.rds"))

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

aa_small_for_subset <- FindVariableFeatures(aa_small_for_subset, 
                                            selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(aa_small_for_subset)
aa_small <- ScaleData(aa_small_for_subset, features = all.genes)

aa_small_for_subset <- RunPCA(aa_small_for_subset, 
                              features = VariableFeatures(object = aa_small_for_subset))

aa_small_for_subset <- FindNeighbors(aa_small_for_subset, dims = 1:8)

aa_small_for_subset <- FindClusters(aa_small_for_subset, resolution = 0.2)

aa_small_for_subset <- RunUMAP(aa_small_for_subset, dims = 1:8)

DimPlot(aa_small_for_subset)
DimPlot(aa_small_for_subset, group.by = "seurat_clusters")
DimPlot(aa_small_for_subset, group.by = "subtypes")
DimPlot(aa_small_for_subset, split.by = "condition")

features = rownames(aa_small)[grep('COL1A2', rownames(aa_small))]
FeaturePlot(aa_small, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))


FB.cluster.ids <- c("EC","CM","Blood","MYeloid","FB","FB", "Bcells")
names(FB.cluster.ids) <- levels(bb)
bb <- RenameIdents(bb, FB.cluster.ids)


Small_markers <- FindAllMarkers(bb, min.pct = 0.6, min.diff.pct = 0.3)

All_Markers = FindAllMarkers(aa, min.pct = 0.6, min.diff.pct = 0.3)


#####
CM_subset_ven <- subset(bb,  subtypes %in% c("Proliferating_CM", 
                                             "Ventricular_CM_ROBO2+", 
                                             "Ventricular_CM_Cav3_1+",
                                             "Ventricular_CM_injury_specific"))
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


features = rownames(aa)[grep('FCN3', rownames(aa))]
VlnPlot(aa, features = features)
FeaturePlot(CM_sub, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))


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


DotPlot(aa_small, features = c( features_CM,features_FB,features_EC , 
                                features_TCELLS, features_B_cells, features_Myeloid,features_RBC , features_MK, features_Neur)) + RotatedAxis()


############
######## CM clusters 3,21,19,15,11,16,1,9,17

CM_subset_gen <- subset(aa,  subtypes %in% c("Ventricular_CM_Cav3_1+", 
                                             "CM_Atria",   "CM_OFT", "Ventricular_CM_ROBO2+" ))

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

features_1 = rownames(FB_subset)[grep('EFEMP1-AMEX60DD035438|PKD1-AMEX60DD030949|PODXL-AMEX60DD008079', rownames(FB_subset))]
features_2 = rownames(FB_subset)[grep('TNXB-AMEX60DD010491|LAMA2-AMEX60DD034687|NGF-AMEX60DD008715', rownames(FB_subset))]
features_3 = rownames(FB_subset)[grep('ACAN-AMEX60DD004178|VWA2-AMEX60DD053181|CSPG5-AMEX60DD026060', rownames(FB_subset))]
features_4 = rownames(FB_subset)[grep('TFPI2-AMEX60DD022268|CXCL14-AMEX60DD028973|HAS1-AMEX60DD017885', rownames(FB_subset))]
features_5 = rownames(FB_subset)[grep('CKAP2L-AMEX60DD001288|IQGAP3-AMEX60DD016824|MASTL-AMEX60DD021858', rownames(FB_subset))]
features_6 = rownames(FB_subset)[grep('MPL-AMEX60DD020059|ADGRG1-AMEX60DD015471|PLEK-AMEX60DD035848', rownames(FB_subset))]


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

all.genes <- rownames(CM_subset_2)
CM_subset_2 <- ScaleData(CM_subset_2, features = all.genes)

CM_subset_2 <- RunPCA(CM_subset_2, features = VariableFeatures(object = CM_subset))

CM_subset <- FindNeighbors(CM_subset, dims = 1:10)


CM_subset <- FindClusters(CM_subset, resolution = 0.3)

CM_subset_2 <- RunUMAP(CM_subset_2, dims = 1:10)

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

EC_subset <- subset(aa,  subtypes %in% c("EC", "EC_IS_(LOX)", "EC_(CEMIP)", "EC_(NOS3)", "EC_(WNT4)", "EC_Prol","EC_IS_(IARS1)","EC_IS_Prol","EC_(LHX6)"))

EC_subset <- FindVariableFeatures(EC_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(EC_subset)
EC_subset <- ScaleData(EC_subset, features = all.genes)

EC_subset <- RunPCA(EC_subset, features = VariableFeatures(object = EC_subset))

EC_subset <- FindNeighbors(EC_subset, dims = 1:10)


#EC_subset <- FindClusters(EC_subset, resolution = 0.25)


EC_subset <- RunUMAP(EC_subset, dims = 1:10)

EC_subset$subtypes -> Idents(EC_subset)

DimPlot(EC_subset)
DimPlot(EC_subset, group.by = "subtypes")
DimPlot(EC_subset,  split.by = "condition")
DimPlot(EC_subset, group.by = "seurat_clusters")

EC_subset_marker <- FindAllMarkers(EC_subset, min.pct = 0.4, min.diff.pct = 0.2) 


features = rownames(EC_subset)[grep('EPHB4', rownames(EC_subset))]
FeaturePlot(EC_subset,  features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))

VlnPlot(EC_subset, features = features)

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
features_CAV = rownames(CM_subset)[grep('CNTFR-AMEX60DD041547|ERV3-1-AMEX60DD016594|HEY2-AMEX60DD034730', rownames(CM_subset))]
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



#####

CM_subset$celltype <- factor(Idents(CM_subset), levels = rev(CMmyLevels))
Idents(CM_subset) <- CM_subset$celltype

DotPlot(CM_subset,
        features = c(features_aCM,
                     features_PM,
                     features_OFT,
                     features_ROBO,
                     features_CAV,
                     features_PR1,
                     features_PR2,
                     features_PR3,
                     features_TAGLN,
                     features_IS,
                     features_IS_PR),
        col.min = 0,
        group.by = "celltype") + 
  RotatedAxis()

############
####
####
features_aCM = rownames(CM_subset)[grep('AGBL1-AMEX60DD004208|ADAMTS17-AMEX60DD003517|GFRA1-AMEX60DD053205', rownames(CM_subset))]
features_PM = rownames(CM_subset)[grep('LTBP2-AMEX60DD011463|HCN4-AMEX60DD003883|TENM2-AMEX60DD030344', rownames(CM_subset))]
features_OFT = rownames(CM_subset)[grep('LHCGR-AMEX60DD035890|FUT2-AMEX60DD017238|ADRA1B-AMEX60DD03002', rownames(CM_subset))]
features_ROBO = rownames(CM_subset)[grep('ROBO2-AMEX60DD047486|JPH3-AMEX60DD01717|STAC-AMEX60DD03796', rownames(CM_subset))]
features_CAV = rownames(CM_subset)[grep('CNTFR-AMEX60DD041547|ERV3-1-AMEX60DD016594|HEY2-AMEX60DD034730', rownames(CM_subset))]
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
Idents(CM_subset) <- CM_subset$subtypes

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

features_Robo2 = rownames(CM_subset_hs)[grep('ROBO2-AMEX60DD047486|JPH3-AMEX60DD017176|STAC-AMEX60DD03796', rownames(CM_subset_hs))]
features_Cav = rownames(CM_subset_hs)[grep('RBPMS-AMEX60DD043590|CNTFR-AMEX60DD041547|SEPTIN9-AMEX60DD030744', rownames(CM_subset_hs))]
features_OFT = rownames(CM_subset_hs)[grep('THSD4-AMEX60DD014724|CRISPLD2-AMEX60DD017214|GRIK2-AMEX60DD034010', rownames(CM_subset_hs))]
features_A1 = rownames(CM_subset_hs)[grep('GFRA1-AMEX60DD053205|ADAMTS17-AMEX60DD003517|AGBL1-AMEX60DD004208', rownames(CM_subset_hs))]
features_Ltbp2 = rownames(CM_subset_hs)[grep('LTBP2-AMEX60DD011463|TENM2-AMEX60DD030344|HCN4-AMEX60DD003883', rownames(CM_subset_hs))]


CMmyLevels <- c("CM_Atria",
                "CM_PM_(HCN4)",
                "CM_OFT",
                "CM_ven_(Cav3_1)",
                "CM_ven_(Robo2)"
                
)




factor(Idents(CM_subset_hs), levels= CMmyLevels)
Idents(CM_subset_hs) <- factor(Idents(CM_subset_hs), levels= CMmyLevels)


DotPlot(CM_subset_hs, col.min = 0, features = c( features_Robo2,features_Cav,features_OFT,features_Ltbp2,features_A1)) + RotatedAxis()
############


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


"#4CC9F0",



"#4361EE",



"#941F56"

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

features_1 = rownames(EC_subset_hs)[grep('LAMA2-AMEX60DD034687|KCNMB2-AMEX60DD030396|RADIL-AMEX60DD025680', rownames(EC_subset_hs))]
features_2 = rownames(EC_subset_hs)[grep('NOS3-AMEX60DD026163|ADAM15-AMEX60DD014884|HSPA12B-AMEX60DD046249', rownames(EC_subset_hs))]
features_3 = rownames(EC_subset_hs)[grep('KCNQ5-AMEX60DD033595|HEG1-AMEX60DD009667|WNT4-AMEX60DD052091', rownames(EC_subset_hs))]
features_4 = rownames(EC_subset_hs)[grep('PTH1R-AMEX60DD020537|SVEP1-AMEX60DD044566|NTRK3-AMEX60DD004199', rownames(EC_subset_hs))]
features_5 = rownames(EC_subset_hs)[grep('PDE2A-AMEX60DD049847|ADGRF5-AMEX60DD032958|LHX6-AMEX60DD050778', rownames(EC_subset_hs))]

ECmyLevels <- c( "EC",
                 "EC_(NOS3)",
                 "EC_(WNT4)",
                 "EC_(CEMIP)",
                 
                 "EC_(LHX6)"
                 
                 
                 
)




factor(Idents(EC_subset_hs), levels= ECmyLevels)
Idents(EC_subset_hs) <- factor(Idents(EC_subset_hs), levels= ECmyLevels)


#DotPlot(EC_subest_2, col.min = 0, features = c(features_EC,features_NOS3,features_WNT4,features_VCAM1,features_LHX6)) + RotatedAxis()
DotPlot(EC_subset_hs, col.min = 0, features = c(features_5,features_4,features_3,features_2,features_1)) + RotatedAxis() & coord_flip()
############

CM_subset_2 <- subset(aa,  subtypes %in%  c(  "CM_ven_(Robo2)", "CM_ven_(Cav3_1)", "CM_IS", "CM_Prol_1", "CM_Prol_3", "CM_Prol_IS"))
DimPlot(CM_subset_2)
CM_subset_2$subtypes -> Idents(CM_subset_2)

CM_subset_2 <- FindVariableFeatures(CM_subset_2, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CM_subset_2)
CM_subset_2 <- ScaleData(CM_subset_2, features = all.genes)

CM_subset_2 <- RunPCA(CM_subset_2, features = VariableFeatures(object = CM_subset_2))

CM_subset_2 <- FindNeighbors(CM_subset_2, dims = 1:10)


CM_subset_2 <- FindClusters(CM_subset_2, resolution = 0.3)
CM_subset_2$subtypes -> Idents(CM_subset_2)

CM_subset_2 <- RunUMAP(CM_subset_2, dims = 1:10)

CM_subset_2_marker <- FindAllMarkers(CM_subset_2, min.pct = 0.4, min.diff.pct = 0.2)

CM_subset_2$subtypes -> Idents(CM_subset_2)

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


CMmyLevels <- c(  "CM_ven_(Robo2)", "CM_ven_(Cav3_1)", "CM_IS","CM_Prol_IS", "CM_Prol_1", "CM_Prol_3")




factor(Idents(CM_subset_2), levels= CMmyLevels)
Idents(CM_subset_2) <- factor(Idents(CM_subset_2), levels= CMmyLevels)

DimPlot(CM_subset_2, group.by = "Phase",  cols = c(
  
  "#4CC9F0",
  
  "#C61010"  ,
  "#414CDC"
  
  
) )

features_Robo2 = rownames(CM_subset_2)[grep('ROBO2-AMEX60DD047486|SCN5A-AMEX60DD021356|HRH2-AMEX60DD027783', rownames(CM_subset_2))]
features_Cav = rownames(CM_subset_2)[grep('CACNA1G-AMEX60DD030985|AGBL1-AMEX60DD004208|CPA6-AMEX60DD039860', rownames(CM_subset_2))]
features_IS = rownames(CM_subset_2)[grep('ACAN-AMEX60DD007108|ADAMTS6-AMEX60DD042234|ARHGAP31-AMEX60DD047024', rownames(CM_subset_2))]
features_PIS = rownames(CM_subset_2)[grep('CACNA1D-AMEX60DD023703|NPPA-AMEX60DD051099|CTF1-AMEX60DD028289', rownames(CM_subset_2))]
features_Pro1 = rownames(CM_subset_2)[grep('CENPF-AMEX60DD036088|ASPM-AMEX60DD018743|ANLN-AMEX60DD03806', rownames(CM_subset_2))]
features_Pro3 = rownames(CM_subset_2)[grep('E2F1-AMEX60DD028458|DIAPH3-AMEX60DD048907|CENPE-AMEX60DD044675', rownames(CM_subset_2))]


CMmyLevels <- c(  "CM_Prol_3", "CM_Prol_1", "CM_Prol_IS","CM_IS", "CM_ven_(Cav3_1)", "CM_ven_(Robo2)")



factor(Idents(CM_subset_2), levels= CMmyLevels)
Idents(CM_subset_2) <- factor(Idents(CM_subset_2), levels= CMmyLevels)


DotPlot(CM_subset_2, col.min = 0, features = c( features_Robo2,features_Cav,features_IS,features_PIS,features_Pro1,features_Pro3)) + RotatedAxis()
############
############

CM_subset <- subset(aa,  subtypes %in%  c(  "CM_ven_(Robo2)", "CM_ven_(Cav3_1)", "CM_IS", "CM_Prol_1", "CM_Prol_3", "CM_Prol_IS"))
CM_subset$subtypes -> Idents(CM_subset)
DimPlot(CM_subset)


CM_subset <- FindVariableFeatures(CM_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CM_subset)
CM_subset <- ScaleData(CM_subset, features = all.genes)

CM_subset <- RunPCA(CM_subset, features = VariableFeatures(object = CM_subset))

CM_subset <- FindNeighbors(CM_subset, dims = 1:10)


CM_subset <- FindClusters(CM_subset, resolution = 0.3)
CM_subset$subtypes -> Idents(CM_subset)

CM_subset <- RunUMAP(CM_subset, dims = 1:10)

CM_subset_marker <- FindAllMarkers(CM_subset, min.pct = 0.4, min.diff.pct = 0.2)

CM_subset$subtypes -> Idents(CM_subset)

CM_IS_VS_ROBO2 <- FindMarkers(aa, ident.1 = "CM_IS", ident.2 = c("CM_ven_(Robo2)"), min.pct = 0.4, min.diff.pct = 0.2)

CMmyLevels <- c(  "CM_ven_(Robo2)", "CM_ven_(Cav3_1)", "CM_IS","CM_Prol_IS", "CM_Prol_1", "CM_Prol_3")




factor(Idents(CM_subset), levels= CMmyLevels)
Idents(CM_subset) <- factor(Idents(CM_subset), levels= CMmyLevels)

DimPlot(CM_subset,  cols = c(
  
  
  "#4CC9F0",
  
  "#49A2F0",
  "#941F56",
  
  "#C61010",  
  
  "#4361EE",
  
  "#3F37C9"
  
  
  
  
) )


FB_subset <- subset(aa,  subtypes %in%  c(  "FB_(PKD1)","FB_(TNXB)", "FB_(VWA2)", "FB_IS_(TFPI2)", "FB_IS_(TNC)", "FB_Prol"))
DimPlot(FB_subset)
FB_subset$subtypes -> Idents(FB_subset)

FB_subset <- FindVariableFeatures(FB_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(FB_subset)
FB_subset <- ScaleData(FB_subset, features = all.genes)

FB_subset <- RunPCA(FB_subset, features = VariableFeatures(object = FB_subset))

FB_subset <- FindNeighbors(FB_subset, dims = 1:10)


FB_subset <- FindClusters(FB_subset, resolution = 0.3)

FB_subset <- RunUMAP(FB_subset, dims = 1:10)
FB_subset$subtypes -> Idents(FB_subset)

FB_subset_marker <- FindAllMarkers(FB_subset, min.pct = 0.4, min.diff.pct = 0.2)

DimPlot(FB_subset, label = TRUE,  cols = c(
  
  "#4CC9F0",
  
  "#49A2F0",
  
  "#4361EE",
  
  "#4042D3",
  
  "#941F56",
  
  "#C61010"  
  
) )

DimPlot(FB_subset, split.by = "condition", pt.size = 0.2, cols = c(
  
  
  "#4CC9F0",
  
  "#49A2F0",
  "#4361EE",
  
  "#3F37C9",  
  
  "#941F56",
  #941F56
  "#C61010"
  
  #C61010
  
  
  
) )


features_1 = rownames(FB_subset)[grep('C3-AMEX60DD032076|SLC29A1-AMEX60DD036508|KRT19-AMEX60DD010100', rownames(FB_subset))]
features_2 = rownames(FB_subset)[grep('TNXB-AMEX60DD010491|ADAMTS15-AMEX60DD053555|POSTN-AMEX60DD049175', rownames(FB_subset))]
features_3 = rownames(FB_subset)[grep('ACAN-AMEX60DD004178|VWA2-AMEX60DD053181|CSPG5-AMEX60DD026060', rownames(FB_subset))]
features_4 = rownames(FB_subset)[grep('TFPI2-AMEX60DD022268|CXCL14-AMEX60DD028973|HAS1-AMEX60DD017885', rownames(FB_subset))]
features_5 = rownames(FB_subset)[grep('TNC-AMEX60DD050822|ERG-AMEX60DD047184|COL11A1-AMEX60DD018809', rownames(FB_subset))]
features_6 = rownames(FB_subset)[grep('ASPM-AMEX60DD018743|CENPF-AMEX60DD036088|IQGAP3-AMEX60DD016824', rownames(FB_subset))]


FBmyLevels <- c(  "FB_Prol","FB_IS_(TNC)", "FB_IS_(TFPI2)", "FB_(VWA2)", "FB_(TNXB)", "FB_(PKD1)")


factor(Idents(FB_subset), levels= FBmyLevels)
Idents(FB_subset) <- factor(Idents(FB_subset), levels= FBmyLevels)


#DotPlot(EC_subest_2, col.min = 0, features = c(features_EC,features_NOS3,features_WNT4,features_VCAM1,features_LHX6)) + RotatedAxis()
DotPlot(FB_subset, col.min = 0, features = c(features_1,features_2,features_3,features_4,features_5,features_6)) + RotatedAxis() 
############

#############
MY_subset <- subset(aa,  subtypes %in% c(
  
  "B_cells_(FOXO1)",
  "B_cells_(SIGLEC11)", 
  "B_cells_Prol",
  "Mo/Macs_(FAXDC2)",
  "Mo/Macs_(SNX22)",
  "Mo/Macs_Prol",
  "Mo/Macs_resident",
  "Neu_(DYSF)",
  "Neu_(IL1R1)",
  "T_cells"
))

#MY_subset_2 <- subset(aa1,  subtypes %in% c("Mono_Macrophages", "Neutrophil","Proliferating_Mono_Macrophages", "Resident_MF" ))

#MY_subset <- subset(MY_subset,  subtypes %in% c("Mono/MF_(CD163L1+)", "B_cells_(MUSK+)", "B_Cells", "T_cells", "B_cells_(Prol)", "Mono/MF_(Prol)"))

atria <- subset (aa, ident = "CM_Atria")

atria <- FindVariableFeatures(atria, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(atria)
atria <- ScaleData(atria, features = all.genes)

atria <- RunPCA(atria, features = VariableFeatures(object = atria))

atria <- FindNeighbors(atria, dims = 1:10)
atria$subtypes -> Idents(atria)

atria <- FindClusters(atria, resolution = 0.08)

atria <- RunUMAP(atria, dims = 1:10)

atria_markers <- FindAllMarkers(atria, min.pct = 0.4, min.diff.pct = 0.2)


features = rownames(aa)[grep('UCHL1', rownames(aa))]
FeaturePlot(aa, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))
VlnPlot(aa, features = features)

features = rownames(atria)[grep('CACNA1H', rownames(atria))]
FeaturePlot(atria, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))

DimPlot(MY_subset)
DimPlot(MY_subset, group.by = "subtypes")
DimPlot(MY_subset, split.by = "condition")
DimPlot(MY_subset, split.by = "condition", group.by = "subtypes")

DimPlot(MY_subset, split.by = "condition",  cols = c(
  "#ADE8F4",
  "#90E0EF",
  "#48CAE4",
  "#0096C7",
  
  "#0077B6",
  "#023E8A",
  "#03045E",
  "#941F56",
  "#AD1833",
  "#FB3640"  
  
  
  
  
  
) )
##MY.cluster.ids <- c("Mono/MF_(CD163L1+)", "B_cells_(MUSK+)", "Mono/MF_(PLBD1+)", "B_Cells", "T_cells", "B_cells_(Prol)", "T_cells", "Mono/MF_(Prol)")
#names(MY.cluster.ids) <- levels(MY_subset)
#MY_subset <- RenameIdents(MY_subset, MY.cluster.ids)

#MY_subset$subtypes -> MY_subset_2$newId_1

#############
MY_subset <- subset(aa,  subtypes %in% c(
  
  
  "Mo/Macs_(FAXDC2)",
  "Mo/Macs_(SNX22)",
  "Mo/Macs_Prol",
  "Mo/Macs_resident",
  "Neu_(DYSF)",
  "Neu_(IL1R1)"
  
))

#MY_subset_2 <- subset(aa1,  subtypes %in% c("Mono_Macrophages", "Neutrophil","Proliferating_Mono_Macrophages", "Resident_MF" ))

#MY_subset <- subset(MY_subset,  subtypes %in% c("Mono/MF_(CD163L1+)", "B_cells_(MUSK+)", "B_Cells", "T_cells", "B_cells_(Prol)", "Mono/MF_(Prol)"))

MY_subset <- FindVariableFeatures(MY_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(MY_subset)
MY_subset <- ScaleData(MY_subset, features = all.genes)

MY_subset <- RunPCA(MY_subset, features = VariableFeatures(object = MY_subset))

MY_subset <- FindNeighbors(MY_subset, dims = 1:10)
MY_subset$subtypes -> Idents(MY_subset)

#MY_subset <- FindClusters(MY_subset, resolution = 0.3)

MY_subset <- RunUMAP(MY_subset, dims = 1:10)

MY_subset_marker <- FindAllMarkers(MY_subset, min.pct = 0.4, min.diff.pct = 0.2)
#AMEX60DD054293 = EPX = MPL

features = rownames(MY_subset)[grep('AMEX60DD054293', rownames(MY_subset))]
FeaturePlot(MY_subset, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))

DimPlot(MY_subset)
DimPlot(MY_subset, group.by = "subtypes")
DimPlot(MY_subset, split.by = "condition")
DimPlot(MY_subset, split.by = "condition", group.by = "subtypes")

DimPlot(MY_subset,split.by = "condition",  cols = c(
  "#ADE8F4",
  
  "#48CAE4",
  "#0096C7",
  
  "#0077B6",
  
  
  "#941F56",
  "#AD1833",
  "#FB3640"  
  
) )

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

EC_subset_2 <- ScaleData(EC_subset, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(EC_subset))

EC_subset_2 <- RunPCA(EC_subset_2, features = VariableFeatures(object = EC_subset_2))

EC_subset_2 <- RunUMAP(EC_subset_2, dims = 1:10)

EC_subset$subtypes -> Idents(EC_subset)



features_1 = rownames(MY_subset)[grep('AXL-AMEX60DD024484|CD163L1-AMEX60DD01641|ADGRL3-AMEX60DD043519', rownames(MY_subset))]
features_2 = rownames(MY_subset)[grep('FAXDC2-AMEX60DD030063|PLBD1-AMEX60DD029125|PARVG-AMEX60DD006445', rownames(MY_subset))]
features_3 = rownames(MY_subset)[grep('LGALS9-AMEX60DD054072|SNX22-AMEX60DD003839|ITGAD-AMEX60DD028068', rownames(MY_subset))]
features_4 = rownames(MY_subset)[grep('NSD2-AMEX60DD046033|ASPM-AMEX60DD018743|SMC4-AMEX60DD00181', rownames(MY_subset))]
features_5 = rownames(MY_subset)[grep('DYSF-AMEX60DD046270|HVCN1-AMEX60DD000672|CFP-AMEX60DD020986', rownames(MY_subset))]
features_6 = rownames(MY_subset)[grep('ARG1-AMEX60DD034655|IL1R1-AMEX60DD048111|CSF3R-AMEX60DD00590', rownames(MY_subset))]


MYmyLevels <- c(
  
  "Neu_(IL1R1)",
  "Neu_(DYSF)",
  "Mo/Macs_Prol",
  "Mo/Macs_(SNX22)",
  "Mo/Macs_(FAXDC2)",
  "Mo/Macs_resident"
  
  
)


factor(Idents(MY_subset), levels= MYmyLevels)
Idents(MY_subset) <- factor(Idents(MY_subset), levels= MYmyLevels)


#DotPlot(EC_subest_2, col.min = 0, features = c(features_EC,features_NOS3,features_WNT4,features_VCAM1,features_LHX6)) + RotatedAxis()
DotPlot(MY_subset, col.min = 0, features = c(features_1,features_2,features_3,features_4,features_5,features_6)) + RotatedAxis() 
############

DimPlot (aa,  cols = c("B_cells_(FOXO1)" =  "#E3F2FD"   ,       "B_cells_(SIGLEC11)" = "#BBDEFB"  ,             "B_cells_Prol"  =  "#90CAF9"  ,           
                       
                       "CM_Atria"  = "#E01E37", "CM_Atria_Tagln" = "#C71F37" ,   "CM_IS" =    "#DA1E37"      , "CM_OFT" =    "#B41E35"     , 
                       "CM_PM_(HCN4)"=  "#A11D33",
                       
                       "CM_Prol_1" = "#931B2F",             "CM_Prol_2" = "#85182A" ,              "CM_Prol_3" = "#801729",        "CM_Prol_IS" = "#7A1627" ,
                       
                       "CM_ven_(Cav3_1)" = "#771626",       "CM_ven_(Robo2)" ="#641220" ,  
                       "EC" = "#007F5F", "EC_(CEMIP)" = "#2B9348"     , "EC_(LHX6)" = "#55A630"  ,  "EC_(NOS3)" = "#80B918",  "EC_(WNT4)" = "#AACC00",
                       
                       "EC_IS_(IARS1)"  = "#BFD200" ,              "EC_IS_(LOX)"= "#D4D700"  ,                "EC_IS_Prol" = "#E6E710",                    "EC_Prol"= "#E2E211",
                       
                       "FB_(PKD1)"=  "#CC5803",                  "FB_(TNXB)" ="#E2711D" ,                  "FB_(VWA2)"= "#FF9505",               "FB_IS_(TFPI2)" ="#FFB627",
                       
                       "FB_IS_(TNC)" = "#FFC04C",                    "FB_Prol"="#FFC971",               "Megakeryocytes" ="#EDC4B3" ,           "Mo/Macs_(FAXDC2)" ="#42A5F5",
                       
                       "Mo/Macs_(SNX22)"= "#64B5F6",                "Mo/Macs_Prol" ="#2196F3",            "Mo/Macs_resident"=  "#1E88E5",                 "Neu_(DYSF)"= "#1565C0",
                       
                       "Neu_(IL1R1)"="#0D47A1" ,                    "Neuronal"="#979dac",    "Proliferating_Megakeryocytes"="#DEAB90"  ,            "Proliferating_RBC"="#9D6B53",
                       
                       "RBC"= "#F3D5B5" ,                    "T_cells"= "#5c677d"))

####

DimPlot (MY_subset, label = TRUE,  cols = c("B_cells_(FOXO1)" =  "#CAF0F8"   ,       "B_cells_(SIGLEC11)" = "#ADE8F4"  ,             "B_cells_Prol"  =  "#90E0EF"  ,           
                                            
                                            "FB_IS_(TFPI2)" ="#48CAE4",
                                            
                                            "Mo/Macs_(FAXDC2)" ="#00B4D8",
                                            
                                            "Mo/Macs_(SNX22)"= "#0096C7",                "Mo/Macs_Prol" ="#0096C7",            "Mo/Macs_resident"=  "#0077B6",                 "Neu_(DYSF)"= "#023E8A",
                                            
                                            "Neu_(IL1R1)"="#03045E" ,                          "T_cells"= "#5c677d"))
#####

head(AverageExpression(CM_subset, pb.method = "aggregate" ,slot = "counts", assays = "RNA"))
head(PseudobulkExpression(CM_subset, slot = "counts", assays = "RNA"))




############

##### color scheme

###Basic

"#BBDEFB",
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


#light-blue - dark blue- purple - red - 10 colors
"#ADE8F4",
"#90E0EF",
"#48CAE4",
"#0096C7",

"#0077B6",
"#023E8A",
"#03045E",
"#941F56",
"#AD1833",
"#FB3640"  



######
features = rownames(EC_subset)[grep('FZD4', rownames(EC_subset))]
VlnPlot(EC_subset, features = features)
FeaturePlot(EC_subset, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))



Tagln_v_IS <- FindMarkers(CM_subset, ident.1 = "CM_Atria_Tagln", ident.2 = "CM_IS", min.pct = 0.4, min.diff.pct = 0.2)

features = rownames(aa)[grep('AXL', rownames(aa))]
FeaturePlot(aa,split.by = "condition" ,features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))

features = rownames(CM_subset_2)[grep('HMOX1', rownames(CM_subset))]
FeaturePlot(CM_subset_2, features = features, pt.size = 0.05, order = TRUE, cols = c("#dee2e6", "#C01E1C"))

VlnPlot(CM_subset, features = features)



features = rownames(CM_subset)[grep('NPPA|NPPB|UCHL1|MYH6|AMEX60DD056342|JUN|AMEX60DD001238|FOXP1', rownames(CM_subset))]
features = rownames(CM_subset)[grep('KAZAL', rownames(CM_subset))]
FeaturePlot(aa, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))
DoHeatmap(subset(CM_subset, downsample = 100), features = features, size = 3)

## view specific cells dimplot
DimPlot(aa, cells = names(aa$condition)[aa$condition == "Amex_scRNA_d0"])

number_cluster <- c()

for (i in levels(aa$condition)) {
  number_cluster <- rbind(number_cluster, summary(aa$subtypes[names(aa$condition)[aa$condition == i]]) )
}

rownames(number_cluster) <- levels(CM_subset$condition)

number_cluster <- c()

for (i in levels(CM_subset$condition)) {
  number_cluster <- cbind(number_cluster, summary(CM_subset$subtypes[names(CM_subset$condition)[CM_subset$condition == i]]) )
}


###########

CM_subset <- subset(aa,  subtypes %in%  c(  "CM_ven_(Robo2)", "CM_ven_(Cav3_1)"))
CM_subset$subtypes -> Idents(CM_subset)
DimPlot(CM_subset)


CM_subset <- FindVariableFeatures(CM_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CM_subset)
CM_subset <- ScaleData(CM_subset, features = all.genes)

CM_subset <- RunPCA(CM_subset, features = VariableFeatures(object = CM_subset))

CM_subset <- FindNeighbors(CM_subset, dims = 1:10)


#CM_subset <- FindClusters(CM_subset, resolution = 0.3)
CM_subset$subtypes -> Idents(CM_subset)

CM_subset <- RunUMAP(CM_subset, dims = 1:10)

CM_subset_marker <- FindAllMarkers(CM_subset, min.pct = 0.4, min.diff.pct = 0.2)

CM_subset$subtypes -> Idents(CM_subset)

CM_IS_VS_ROBO2 <- FindMarkers(aa, ident.1 = "CM_IS", ident.2 = c("CM_ven_(Robo2)"), min.pct = 0.4, min.diff.pct = 0.2)

CMmyLevels <- c(  "CM_ven_(Robo2)", "CM_ven_(Cav3_1)", "CM_IS","CM_Prol_IS", "CM_Prol_1", "CM_Prol_3")

#########
CM_subset1 <- subset(aa,  subtypes %in%  c(  "CM_PM_(HCN4)", "CM_OFT","CM_Atria","CM_ven_(Robo2)", "CM_ven_(Cav3_1)"))
DimPlot(CM_subset)
CM_subset$subtypes -> Idents(CM_subset)

CM_subset1 <- FindVariableFeatures(CM_subset1, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(CM_subset1)
CM_subset1 <- ScaleData(CM_subset1, features = all.genes)

CM_subset1 <- RunPCA(CM_subset1, features = VariableFeatures(object = CM_subset1))

CM_subset1 <- FindNeighbors(CM_subset1, dims = 1:10)



CM_subset1 <- RunUMAP(CM_subset, dims = 1:10)



DimPlot(CM_subset, split.by = "condition", group.by = "Phase", cols =c(
  
  "#DDE2E6",
  
  
  
  "#0096C9",
  
  
  
  "#B10A2E"))

FeaturePlot(CM_subset, features = "S.Score", order = TRUE, cols = c("#dee2e6", "#661CB0"))


######### CM_subset thresholding for cell cycle scoring

CM_subset$True_G2M <- 1*(CM_subset$G2M.Score>0.2)
CM_subset$True_G2M[CM_subset$True_G2M == 0] <- "Other"
CM_subset$True_G2M[CM_subset$True_G2M == 1] <- "G2M"
DimPlot(CM_subset, split.by = "condition", group.by = "True_G2M", cols = c("#661CB0", "#dee2e6"))
table(data.frame(CM_subset$True_G2M, CM_subset$condition))

CM_subset$True_S <- 1*(CM_subset$S.Score>0.13)
CM_subset$True_S[CM_subset$True_S == 0] <- "Other"
CM_subset$True_S[CM_subset$True_S == 1] <- "S_Phase"
DimPlot(CM_subset, split.by = "condition", group.by = "True_S", cols = c("#dee2e6", "#661CB0"))
table(data.frame(CM_subset$True_S, CM_subset$condition))

hist(CM_subset$G2M.Score, breaks = 60)
hist(CM_subset$S.Score, breaks = 60)

hist(aa$G2M.Score, breaks = 60)
hist(aa$S.Score, breaks = 60)


G2MPhaser <- function (G2M_Scores,S_Scores,G2M_Thresh=0, S_Thresh=0) {
  PhaseName <- rep('Other',length(G2M_Scores))
  PhaseName[(S_Scores > S_Thresh ) & (S_Scores > G2M_Scores)] <- 'S'
  PhaseName[(G2M_Scores > G2M_Thresh ) & (S_Scores < G2M_Scores)] <- 'G2M'
  return(PhaseName)
}


CM_subset$Phase_threshold <- G2MPhaser(CM_subset$G2M.Score,CM_subset$S.Score,0.22,0.13)
names(CM_subset$Phase_threshold) <- Cells(CM_subset)
DimPlot(CM_subset, group.by = "Phase_threshold", cols = c(
  
  "#B10A2E",
  
  
  
  "#DDE2E6",
  
  
  
  "#0096C9"), split.by = "condition")

aa$Phase_threshold <- G2MPhaser(aa$G2M.Score,aa$S.Score,0.22,0.13)
names(aa$Phase_threshold) <- Cells(aa)
DimPlot(aa, group.by = "Phase_threshold", cols = c("#ed7306", "#e5e5e5", "#671bb2"))

table(CM_subset$Phase_threshold)
#######

DimPlot(aa, cols = c("CM" = "#B4D034",  "EC" ="#C1D857",  "FB" =  "#9E7C29",  "Megakaryocytes"=  "#CE9B1E", "RBC" = "#FDBA13",  "Myeloid" =  "#EBC01C", "B_cells"="#D9C524", "T_cells"= "#3E3E3E",  "Neuronal" = "#6E5D34"))

##### TO look for markers across time within the same cluster
NEU_V_Neu <- FindMarkers(MY_subset, ident.1 = "Neu_(IL1R1)" , ident.2 = "Neu_(DYSF)"  ,min.pct = 0.6, min.diff.pct = 0.4)

PKD1_marker <- FindMarkers(FB_subset, ident.1 = (names(FB_subset$subtypes)[(FB_subset$subtypes == "FB_(PKD1)") & (FB_subset$condition %in% c("Amex_scRNA_d0"))]), ident.2 = (names(FB_subset$subtypes)[(FB_subset$subtypes == "FB_(PKD1)") & (FB_subset$condition %in% c("Amex_scRNA_d4","Amex_scRNA_d7","Amex_scRNA_d14"))]), min.pct = 0.4, min.diff.pct = 0.2)

EC_marker_lox <- FindMarkers(EC_subset, ident.1 = "EC_IS_(LOX)",  min.pct = 0.4, min.diff.pct = 0.2)
CM_IS_marker <- FindMarkers(CM_subset, ident.1 = "CM_IS",  min.pct = 0.4, min.diff.pct = 0.2)
CM_IS_marker_time <- FindMarkers(CM_subset, ident.1 = (names(CM_subset$subtypes)[(CM_subset$subtypes == "CM_IS") & (CM_subset$condition %in% c("Amex_scRNA_d1"))]), ident.2 = (names(CM_subset$subtypes)[(CM_subset$subtypes == "CM_IS") & (CM_subset$condition %in% c("Amex_scRNA_d14"))]), min.pct = 0.4, min.diff.pct = 0.2)
##################
#d1 = CM IS
CM_IS_day1_vs_cav31 <- FindMarkers(CM_subset, ident.1 = (names(CM_subset$subtypes)[(CM_subset$subtypes == "CM_IS") & (CM_subset$condition %in% c("Amex_scRNA_d1"))]), ident.2 = ("CM_ven_(Cav3_1)"), min.pct = 0.4, min.diff.pct = 0.2)
#d4 = Prol IS
CM_Prol_IS_day4_vs_IS_cav31 <- FindMarkers(CM_subset, ident.1 = (names(CM_subset$subtypes)[(CM_subset$subtypes == "CM_Prol_IS") & (CM_subset$condition %in% c("Amex_scRNA_d4"))]), ident.2 = c("CM_ven_(Cav3_1)", "CM_IS"), min.pct = 0.4, min.diff.pct = 0.2)
#d7 = Prol IS
CM_Prol_IS_day7_vs_IS_cav31 <- FindMarkers(CM_subset, ident.1 = (names(CM_subset$subtypes)[(CM_subset$subtypes == "CM_Prol_IS") & (CM_subset$condition %in% c("Amex_scRNA_d7"))]), ident.2 = c("CM_ven_(Cav3_1)", "CM_IS"), min.pct = 0.4, min.diff.pct = 0.2)
#d14 = CM IS
CM_IS_day14_vs_Prol_IS_cav31<- FindMarkers(CM_subset, ident.1 = (names(CM_subset$subtypes)[(CM_subset$subtypes == "CM_IS") & (CM_subset$condition %in% c("Amex_scRNA_d14"))]), ident.2 = c("CM_ven_(Cav3_1)", "CM_Prol_IS"), min.pct = 0.4, min.diff.pct = 0.2)

################## vs cav31
#d1 = CM IS
CM_IS_day1_vs_cav31 <- FindMarkers(CM_subset, ident.1 = (names(CM_subset$subtypes)[(CM_subset$subtypes == "CM_IS") & (CM_subset$condition %in% c("Amex_scRNA_d1"))]), ident.2 = ("CM_ven_(Cav3_1)"), min.pct = 0.4, min.diff.pct = 0.2)
#d4 = Prol IS
CM_Prol_IS_day4_vs_IS_cav31 <- FindMarkers(CM_subset, ident.1 = (names(CM_subset$subtypes)[(CM_subset$subtypes == "CM_Prol_IS") & (CM_subset$condition %in% c("Amex_scRNA_d4"))]), ident.2 = c("CM_ven_(Cav3_1)"), min.pct = 0.4, min.diff.pct = 0.2)
#d7 = Prol IS
CM_Prol_IS_day7_vs_IS_cav31 <- FindMarkers(CM_subset, ident.1 = (names(CM_subset$subtypes)[(CM_subset$subtypes == "CM_Prol_IS") & (CM_subset$condition %in% c("Amex_scRNA_d7"))]), ident.2 = c("CM_ven_(Cav3_1)"), min.pct = 0.4, min.diff.pct = 0.2)
#d14 = CM IS
CM_IS_day14_vs_Prol_IS_cav31<- FindMarkers(CM_subset, ident.1 = (names(CM_subset$subtypes)[(CM_subset$subtypes == "CM_IS") & (CM_subset$condition %in% c("Amex_scRNA_d14"))]), ident.2 = c("CM_ven_(Cav3_1)"), min.pct = 0.4, min.diff.pct = 0.2)


############
features = rownames(CM_subset)[grep('HEY2|TBX20|SMAD3|ETS1|E2F8|MYBL1', rownames(CM_subset))]
features = rownames(CM_subset)[grep('CD24', rownames(CM_subset))]
VlnPlot(CM_subset, features = features, pt.size = 0 )
FeaturePlot(CM_subset, features = features,  order = TRUE, cols = c("#dee2e6", "#661CB0"))


Prox <- subset(aa,  subtypes %in% c("CM_ven_(Robo2)", "FB_(PKD1)", "EC_IS_(LOX)", "CM_Prol_IS", "Mo/Macs_(FAXDC2)" ))
FeaturePlot(CM_subset, features = features,  order = TRUE, cols = c("#dee2e6", "#661CB0"))
DotPlot(CM_subset, col.min = 0 ,features = c(features)) + RotatedAxis()  
VlnPlot(CM_subset, features = features, pt.size = 0 )
RidgePlot(CM_subset, features = features, ncol = 2)

CMmyLevels <- c(  "CM_ven_(Robo2)", "CM_Prol_3","CM_Prol_1","CM_ven_(Cav3_1)", "CM_IS","CM_Prol_IS" )

Idents(CM_subset) <- factor(Idents(CM_subset), levels= CMmyLevels)

####### grab featuers from table 
ZFP_Ids <- c()
for (i in table[,1]){
  ZFP_Ids <- c(ZFP_Ids,rownames(aa)[grep(i, rownames(aa))])
}

DoHeatmap(subset(CM_subset, downsample = 100), features = ZFP_Ids, size = 3)

######### remove the gene identifier for "axo to human conversion" 


############### FB_subtype for paper writing

FB_subset <- subset(aa,  subtypes %in% c(  "FB_(PKD1)", "FB_(TNXB)", "FB_(VWA2)", "FB_IS_(TFPI2)"  , 
                                           "FB_IS_(TNC)"  ,  "FB_Prol"      ))

FB_subset <- FindVariableFeatures(FB_subset, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(FB_subset)
FB_subset <- ScaleData(FB_subset, features = all.genes)

FB_subset <- RunPCA(FB_subset, features = VariableFeatures(object = FB_subset))

FB_subset <- FindNeighbors(FB_subset, dims = 1:10)


FB_subset <- FindClusters(FB_subset, resolution = 0.1)


FB_subset <- RunUMAP(FB_subset, dims = 1:10)


FB_subset$subtypes = Idents(FB_subset)

FB_subset_marker <- FindAllMarkers(FB_subset, min.pct = 0.4, min.diff.pct = 0.2)

FB_subset_PKD1 <- subset(aa,  subtypes %in% c(  "FB_(PKD1)", "FB_IS_(TFPI2)" ))   


FB_subset_marker_PKD <- FindAllMarkers(FB_subset_PKD1, min.pct = 0.4, min.diff.pct = 0.2)

FB_subset_marker_PKD_2 <- FindMarkers(aa, ident.1 = "FB_(PKD1)", ident.2 = c("FB_IS_(TFPI2)"), min.pct = 0.4, min.diff.pct = 0.2)

features = rownames(CM_sub)[grep('GJA5|EPHB3|CACNA2D2', rownames(CM_sub))]
VlnPlot(CM_sub,features = features )
FeaturePlot(CM_sub, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))


########################



MY_subset <- subset(aa,  subtypes %in% c("Mo/Macs_(FAXDC2)","Neu_(DYSF)","Mo/Macs_(SNX22)" ,"Mo/Macs_resident","Mo/Macs_Prol","doublet","Neu_(IL1R1)"))

MY_subset <- FindVariableFeatures(MY_subset, selection.method = "vst", nfeatures = 2000)


all.genes <- rownames(MY_subset)
MY_subset <- ScaleData(MY_subset, features = all.genes)

MY_subset <- RunPCA(MY_subset, features = VariableFeatures(object = MY_subset))

MY_subset <- FindNeighbors(MY_subset, dims = 1:10)


MY_subset <- FindClusters(MY_subset, resolution = 0.1)


MY_subset <- RunUMAP(MY_subset, dims = 1:10)


MY_subset$subtypes = Idents(MY_subset)

features = rownames(aa)[grep('AMEX60DD025027', rownames(aa))]
FeaturePlot(aa, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))
VlnPlot(aa,features = features, pt.size = 0 ,y.max = 4)+ NoLegend()



#### LYZ AMEX60DD007762
#####MMP9 AMEX60DD027684

MY_subset_marker <- FindAllMarkers(MY_subset, min.pct = 0.4, min.diff.pct = 0.2)
###################
DimPlot(EC_subset, split.by = "condition", pt.size = 0.2, cols = c(
  "#4CC9F0",
  "#4AAFF0",
  "#4895EF",
  "#4361EE",
  
  
  "#941F56",
  "#AD1833",
  "#C61010",
  "#03045E"))

DimPlot(EC_subset, split.by = "condition", pt.size = 0.2, cols = c(
  "#4CC9F0",
  "#4CC9F0",
  "#4CC9F0",
  "#4CC9F0",
  
  "#4CC9F0",
  "#AD1833",
  "#3F37C9",
  
  "#4CC9F0",
  "#4CC9F0"
))

DimPlot(EC_subset,  pt.size = 0.2, cols = c(
  "#4CC9F0",
  "#4AAFF0",
  "#4895EF",
  "#4361EE",
  
  "#4CC9F0",
  "#AD1833",
  "#3F37C9",
  "#C61010",
  "#941F56",
  "#03045E"
))
#############

DimPlot (aa,split.by = "condition",  cols = c("B_cells_(FOXO1)" =  "#dee2e6"   ,       "B_cells_(SIGLEC11)" = "#dee2e6"  ,             "B_cells_Prol"  =  "#dee2e6"  ,           
                                              
                                              "CM_Atria"  = "#dee2e6", "CM_Atria_Tagln" = "#C71F37" ,   "CM_IS" =    "#DA1E37"      , "CM_OFT" =    "#dee2e6"     , 
                                              "CM_PM_(HCN4)"=  "#dee2e6",
                                              
                                              "CM_Prol_1" = "#dee2e6",             "CM_Prol_2" = "#dee2e6" ,              "CM_Prol_3" = "#dee2e6",        "CM_Prol_IS" = "#7A1627" ,
                                              
                                              "CM_ven_(Cav3_1)" = "#dee2e6",       "CM_ven_(Robo2)" ="#dee2e6" ,  
                                              "EC" = "#dee2e6", "EC_(CEMIP)" = "#dee2e6"     , "EC_(LHX6)" = "#dee2e6"  ,  "EC_(NOS3)" = "#dee2e6",  "EC_(WNT4)" = "#dee2e6",
                                              
                                              "EC_IS_(IARS1)"  = "#BFD200" ,              "EC_IS_(LOX)"= "#D4D700"  ,                "EC_IS_Prol" = "#E6E710",                    "EC_Prol"= "#dee2e6",
                                              
                                              "FB_(PKD1)"=  "#dee2e6",                  "FB_(TNXB)" ="#dee2e6" ,                  "FB_(VWA2)"= "#dee2e6",               "FB_IS_(TFPI2)" ="#FFB627",
                                              
                                              "FB_IS_(TNC)" = "#FFC04C",                    "FB_Prol"="#dee2e6",               "Megakeryocytes" ="#dee2e6" ,           "Mo/Macs_(FAXDC2)" ="#42A5F5",
                                              
                                              "Mo/Macs_(SNX22)"= "#64B5F6",                "Mo/Macs_Prol" ="#dee2e6",            "Mo/Macs_resident"=  "#dee2e6",                 "Neu_(DYSF)"= "#dee2e6",
                                              
                                              "Neu_(IL1R1)"="#0D47A1" ,                    "Neuronal"="#dee2e6",    "Proliferating_Megakeryocytes"="#dee2e6"  ,            "Proliferating_RBC"="#dee2e6",
                                              
                                              "RBC"= "#dee2e6" ,                    "T_cells"= "#dee2e6"))

###################



DimPlot (aa,  cols = c("B_cells_(FOXO1)" =  "#E3F2FD"   ,       "B_cells_(SIGLEC11)" = "#BBDEFB"  ,             "B_cells_Prol"  =  "#90CAF9"  ,           
                       
                       "CM_Atria"  = "#E01E37", "CM_Atria_Tagln" = "#dee2e6" ,   "CM_IS" =    "#dee2e6"      , "CM_OFT" =    "#B41E35"     , 
                       "CM_PM_(HCN4)"=  "#A11D33",
                       
                       "CM_Prol_1" = "#931B2F",             "CM_Prol_2" = "#85182A" ,              "CM_Prol_3" = "#801729",        "CM_Prol_IS" = "#dee2e6" ,
                       
                       "CM_ven_(Cav3_1)" = "#771626",       "CM_ven_(Robo2)" ="#641220" ,  
                       "EC" = "#007F5F", "EC_(CEMIP)" = "#2B9348"     , "EC_(LHX6)" = "#55A630"  ,  "EC_(NOS3)" = "#80B918",  "EC_(WNT4)" = "#AACC00",
                       
                       "EC_IS_(IARS1)"  = "#dee2e6" ,              "EC_IS_(LOX)"= "#dee2e6"  ,                "EC_IS_Prol" = "#dee2e6",                    "EC_Prol"= "#E2E211",
                       
                       "FB_(PKD1)"=  "#CC5803",                  "FB_(TNXB)" ="#E2711D" ,                  "FB_(VWA2)"= "#FF9505",               "FB_IS_(TFPI2)" ="#dee2e6",
                       
                       "FB_IS_(TNC)" = "#dee2e6",                    "FB_Prol"="#FFC971",               "Megakeryocytes" ="#EDC4B3" ,           "Mo/Macs_(FAXDC2)" ="#dee2e6",
                       
                       "Mo/Macs_(SNX22)"= "#dee2e6",                "Mo/Macs_Prol" ="#2196F3",            "Mo/Macs_resident"=  "#1E88E5",                 "Neu_(DYSF)"= "#1565C0",
                       
                       "Neu_(IL1R1)"="#dee2e6" ,                    "Neuronal"="#979dac",    "Proliferating_Megakeryocytes"="#DEAB90"  ,            "Proliferating_RBC"="#9D6B53",
                       
                       "RBC"= "#F3D5B5" ,                    "T_cells"= "#5c677d"))

####
DimPlot(aa, group.by = "condition", cols = c("Amex_scRNA_d0" =  "#7A1627" , "Amex_scRNA_d1" =  "#03045e","Amex_scRNA_d4" =  "#0077b6","Amex_scRNA_d7" =  "#48cae4","Amex_scRNA_d14" =  "#caf0f8")
)
####
RBC <- subset(aa,  subtypes %in% c("RBC", "Proliferating_RBC" ))
Idents(RBC) <- RBC$condition
RBC_marker_time <- FindAllMarkers(RBC, min.pct = 0.3, min.diff.pct = 0.2)

####
DimPlot (MY_subset, label = TRUE,  cols = c("B_cells_(FOXO1)" =  "#CAF0F8"   ,       "B_cells_(SIGLEC11)" = "#ADE8F4"  ,             "B_cells_Prol"  =  "#90E0EF"  ,           
                                            
                                            "FB_IS_(TFPI2)" ="#48CAE4",
                                            
                                            "Mo/Macs_(FAXDC2)" ="#00B4D8",
                                            
                                            "Mo/Macs_(SNX22)"= "#0096C7",                "Mo/Macs_Prol" ="#0096C7",            "Mo/Macs_resident"=  "#0077B6",                 "Neu_(DYSF)"= "#023E8A",
                                            
                                            "Neu_(IL1R1)"="#03045E" ,                          "T_cells"= "#5c677d"))
#####
DimPlot(CM_subset, split.by = "condition", cols = c(
  
  "#DDE2E6",
  
  
  
  "#0096C9",
  
  
  
  "#B10A2E"))


########## To subset from the metadata
day0_subset <- subset(aa, subset = condition=='Amex_scRNA_d0' )

### to get values in a table such as counting the number of cells per class
table(data.frame(CM_subset$Phase, CM_subset$condition))

## saverds
saveRDS(EC_subset_2, "/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/EC_subset_regress.rds")

write.csv(CM_subset_marker, "/groups/tanaka/People/current/Elad/CM_subset_marker.csv", quote = F)


#AMEX60DD030400 = KCNIP1
#AMEX60DD034639 = TMEM200A
#AMEX60DD047596 = NDP
#AMEX60DD025027 = non-annotated... new? gtra... tmem272 TM like


features = rownames(aa)[grep('FZD6', rownames(aa))]
VlnPlot(aa, features = features) + NoLegend()
FeaturePlot(EC_subset, features = features, order = TRUE, cols = c("#dee2e6", "#661CB0"))
