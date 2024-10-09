##########################################################################
##########################################################################
# Project: heart regeneration project
# Script purpose: process and analyze the scRNA-seq data and visium data from Kuppe_et_al_2022
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Feb 22 11:34:57 2024
##########################################################################
##########################################################################
rm(list = ls())

version.analysis = '_20231207'

resDir = paste0("../results/scRNAseq_human", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

# required libraries
library(data.table)
require(Seurat)
library(SeuratObject)
require(sctransform)
require(ggplot2)
library(dplyr)
library(patchwork)
require(tictoc)
library(DropletUtils)
library(edgeR)
library(future)
options(future.globals.maxSize = 160000 * 1024^2)
library(pryr) # monitor the memory usage

mem_used()
Normalization = 'lognormal' # ('lognormal or SCT')

dataDir = "../published_dataset/human/Kuppe_et_al_2022/processed_data_Robj/snRNA/"

# annot = import(paste0('/groups/tanaka/People/current/jiwang/Genomes/human/hg38/annotation/', 
#                       'Homo_sapiens.GRCh38.111.gtf'))
# 
# annot = data.frame(id = annot$gene_id, symbol = annot$gene_name, type = annot$gene_biotype)
# annot = annot[which(annot$type == 'protein_coding'), ]
# annot = annot[match(unique(annot$id), annot$id), ]
# 
# saveRDS(annot, file = paste0(RdataDir, 'human_gtf_ens_geneSymbol.rds'))

########################################################
########################################################
# Section I : Import the processed seurat object of scRNA-seq from
# Kuppe et al., 2022
# object were downloaded in https://cellxgene.cziscience.com/collections/8191c283-0816-424b-9b61-c3e1d6258a77
########################################################
########################################################
aa = readRDS(file = paste0(dataDir, 'all-snRNA.rds'))

p1 = DimPlot(aa, reduction = 'umap', group.by = 'cell_type_original', raster=FALSE, label = TRUE, repel = TRUE)

p2 = DimPlot(aa, reduction = 'umap', group.by = 'final_cluster', raster=FALSE, label = TRUE, repel = TRUE)

p1 + p2 

ggsave(paste0(resDir, '/Kupper2022_Umap_clusters_cellType.original.pdf'), 
       width = 16, height = 8)

p1 = DimPlot(aa, reduction = 'umap', group.by = 'patient_region_id', raster=FALSE, label = TRUE, repel = TRUE)
p2 = DimPlot(aa, reduction = 'umap', group.by = 'major_labl', raster=FALSE, label = TRUE, repel = TRUE)
p1 + p2

ggsave(paste0(resDir, '/Kupper2022_Umap_patienceID_major_labels.pdf'), 
       width = 16, height = 8)

saveRDS(aa, file = paste0(RdataDir, '/Kuppe2022_heart_all.rds'))

##########################################
# subset the snRNA-seq with the donor based on the metadata (41586_2022_5060_MOESM6_ESM)
# donors with days after infarction >30 were filtered
##########################################
Idents(aa) = aa$donor_id

aa = subset(aa, idents = c("P14", "P18", 'P20', 'P19', 'P4', 'P13', 'P11', 'P12', 'P5'), invert = TRUE)

DimPlot(aa, reduction = 'umap', group.by = 'donor_id', raster=FALSE, label = TRUE, repel = TRUE)

ggsave(paste0(resDir, '/Kupper2022_Umap_selectedDonor.pdf'), 
       width = 16, height = 12)


p1 = DimPlot(aa, reduction = 'umap', group.by = 'cell_type_original', raster=FALSE, label = TRUE, repel = TRUE)

p2 = DimPlot(aa, reduction = 'umap', group.by = 'donor_id', raster=FALSE, label = TRUE, repel = TRUE)

p1 + p2 

ggsave(paste0(resDir, '/Kupper2022_Umap_clusters_cellType.original_selectedDonor.pdf'), 
       width = 16, height = 8)


saveRDS(aa, file = paste0(RdataDir, '/Kuppe2022_heart_donorSelected.rds'))


########################################################
########################################################
# Section II: import and process the human visium data
# 
########################################################
########################################################
visiumDir = "../published_dataset/human/Kuppe_et_al_2022/processed_data_Robj/visium"

outDir = paste0(resDir, '/visium_data/')
if(!dir.exists(outDir)) dir.create(outDir)

st_file = list.files(path = visiumDir,
                        pattern = '*.rds', full.names = TRUE)

annot = readRDS(file = paste0(RdataDir, 'human_gtf_ens_geneSymbol.rds'))

##########################################
# manually select the control, BZ, RZ  
##########################################
st_file = st_file[grep('control_P|RZ_BZ_P|RZ_P|IZ_BZ_P2', st_file)]
st_file = st_file[grep('RZ_P9|RZ_BZ_P12', st_file, invert = TRUE)]

for(n in 1:length(st_file))
{
  # n = 2
  cc = gsub('.rds', '', basename(st_file[n]))
  cc = gsub('Visium-', '', cc)
  cat(n, ' -- ', cc, '\n')
  aa = readRDS(st_file[n])
  aa$sample = cc
  
  cat(nrow(aa), ' genes -- ', ncol(aa), ' spots\n')
  
  pdfname = paste0(outDir, '/QCs_gene_marker_check_', cc, '.pdf')
  pdf(pdfname, width=16, height = 8)
  
  p1 = VlnPlot(aa, features = c("nCounts_RNA", "nFeaturess_RNA", "percent.mt"), 
               ncol = 3, pt.size = 1.0)
  plot(p1)
  
  plot1 <- FeaturePlot(aa, features = "nCounts_RNA", reduction = 'spatial') 
  plot2 <- FeaturePlot(aa, features = "nFeaturess_RNA", reduction = 'spatial') 
  plot3 <- FeaturePlot(aa, features = "percent.mt", reduction = 'spatial') 
  plot(wrap_plots(plot1, plot2, plot3))
  
  markers = c('Nppa', 'Myh6', 'Lum', 'Itgb5', 'Vim', 'Col1a1')
  markers = toupper(markers)
  mm = match(annot$id[match(markers, annot$symbol)], rownames(aa))
  
  p4 = FeaturePlot(aa, features = rownames(aa)[mm], reduction = 'spatial', ncol = 3) 
  plot(p4)
  ##########################################
  # normalization 
  ##########################################
  # min_cells = 5 only use genes that have been detected in at least this many cells 
    
  # merge slices from different time points and 
  if(n == 1) {
    st = aa
  }else{
    st = merge(st, aa)
  }
  
  remove(aa)
  
  dev.off()
  
}

save(st, file = paste0(RdataDir,  'seuratObject_visium_humanHeart_selected_control_BZ_RZ.rds'))


########################################################
########################################################
# Section III: # add subtypes with the original subclusterd based on snRNA and scATAC
# original R objects from https://zenodo.org/records/7098004#.Y0P_LC0RoeY
# 
########################################################
########################################################
#library(org.Mm.eg.db)
library('rtracklayer')
library(GenomicRanges)
library('GenomicFeatures')

annot = readRDS(file = paste0(RdataDir, 'human_gtf_ens_geneSymbol.rds'))

aa = readRDS(file = paste0(RdataDir, '/Kuppe2022_heart_donorSelected.rds'))
aa$patient_region_id = droplevels(aa$patient_region_id)

aa$cell.id = sapply(colnames(aa), function(x) {unlist(strsplit(x, '_'))[1]})
aa$cell.id = paste0(aa$patient_region_id, '_', aa$cell.id)

# ax = readRDS(file = paste0())

DimPlot(aa, reduction = 'umap', group.by = 'cell_type_original', 
        raster=FALSE, label = TRUE, repel = TRUE)

annot_file = list.files(path = paste0('../published_dataset/human/Kuppe_et_al_2022/',
                        'processed_data_Robj/subtype_annot'),
                        pattern = '*.Rds', full.names = TRUE)

aa$annotation = NA

##########################################
# systematic annotate the subtypes 
##########################################
for(n in 1:length(annot_file))
{
  # n = 2
  celltype = gsub('_snRNA_snATAC.Rds', '', basename(annot_file[n]))
  
  outDir = paste0(resDir, '/', celltype)
  if(!dir.exists(outDir)) dir.create(outDir)
  
  cat(n, '--', basename(annot_file[n]), '-- cell type : ', celltype, '\n')
  xx = readRDS(annot_file[n])
  cat(nrow(xx), 'cells \n')
  print(table(xx$annotation))
  
  #xx$cell.id = sapply(colnames(xx), function(x) {unlist(strsplit(x, '#'))[2]})
  #xx$cell.id = paste0(xx$patient_region_id, '_', xx$cell.id)
  
  xx = DietSeurat(xx, data = TRUE, assays = 'RNA')
  
  subs = subset(aa, cells = colnames(aa)[which(aa$cell_type_original == celltype)])
  
  DimPlot(subs, group.by = 'cell_type_original', raster=FALSE, label = TRUE, repel = TRUE)
  
  subs <- FindVariableFeatures(subs, selection.method = 'vst', nfeatures = 2000) %>% 
    ScaleData() %>% RunPCA()
  
  subs <- RunUMAP(subs, dims = 1:20, reduction = 'pca', n.neighbors = 30, min.dist = 0.1)
  
  subs <- RunUMAP(subs, dims = 1:20, reduction = 'harmony', n.neighbors = 30, min.dist = 0.1,
                  reduction.name = "harmony_umap")
  
  p1 = DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'annotation', raster=FALSE)
  p2 =  DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'annotation', raster=FALSE, 
                reduction = 'harmony_umap')
  
  p1 + p2 
  
  #ElbowPlot(subs, ndims = 50)
  subs <- FindNeighbors(subs, dims = 1:20, reduction = 'harmony')
  subs <- FindClusters(subs, verbose = FALSE, algorithm = 3, resolution = 1.0)
  
 
  
  p1 = DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'cell_type_original', raster=FALSE)
  p2 =  DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  
  p1 + p2 
  ggsave(filename = paste0(outDir, '/umap_clusters_celltype_', celltype, '.pdf'), width = 16, height = 8)
  
  # pseudo-bulk by per donor per cell type
  pb <- AverageExpression(subs, return.seurat = TRUE, slot = 'counts', features = VariableFeatures(subs),
                          group.by = c("seurat_clusters"))
  clusters = pb@assays$RNA@data
  rm(pb)
  
  ref = AverageExpression(xx, return.seurat = TRUE, slot = 'counts',
                          group.by = c("annotation"))
  ref = ref@assays$RNA@data
  
  ggs = annot$symbol[match(rownames(clusters), annot$id)]
  
  gg_sels = unique(ggs[which(!is.na(ggs))])
  clusters = clusters[match(gg_sels, ggs), ]
  rownames(clusters) = gg_sels
  rm(ggs)
  
  ggs = intersect(rownames(clusters), rownames(ref))
  clusters = clusters[match(ggs, rownames(clusters)), ]
  ref = ref[match(ggs, rownames(ref)), ]
  
  keep = matrix(NA, nrow = ncol(clusters), ncol = 4)
  rownames(keep) = colnames(clusters)
  colnames(keep) = c('cluster.pearson', 'cluster.spearman',  
                     'cor.pearson', 'cor.spearman')
  
  for(m in 1:ncol(clusters))
  {
    # m = 1
    cat(m, ' -- ', colnames(clusters)[m], '\n')
    cors = cor(clusters[,m], ref, use = 'everything', method = 'pearson')
    keep[m, 1] = colnames(ref)[which.max(cors)]
    keep[m, 3] = cors[which.max(cors)]
    
    cors = cor(clusters[,m], ref, use = 'everything', method = 'spearman')
    keep[m, 2] = colnames(ref)[which.max(cors)]
    keep[m, 4] = cors[which.max(cors)]
    
  }
  
  keep = data.frame(keep)
  
  for(m in 1:nrow(keep))
  {
    # m = 1
    if(as.numeric(keep$cor.spearman[m]) > 0.7 & keep$cluster.pearson[m] == keep$cluster.spearman[m]){
      cat(rownames(keep)[m], '\n')
      cells = colnames(subs)[which(subs$seurat_clusters == rownames(keep)[m])]
      aa$annotation[match(cells, colnames(aa))] = keep$cluster.spearman[m]
    }
    
  }
  
  #DimPlot(aa, reduction = 'umap', group.by = 'cell_type_original', raster=FALSE, label = TRUE, repel = TRUE)
  #ggsave(filename = paste0(outDir, '/umap_celltype.pdf'), width = 12, height = 8)
  
  DimPlot(aa, reduction = 'umap', group.by = 'annotation', raster=FALSE, label = TRUE, repel = TRUE)
  ggsave(filename = paste0(outDir, '/umap_subtypes_for_celltype_', celltype, '.pdf'), width = 12, height = 8)
  
  saveRDS(aa, file = paste0(RdataDir, '/Kuppe2022_heart_donorSelected_subtypes_', celltype, '.rds'))
  
}

saveRDS(aa, file = paste0(RdataDir, '/Kuppe2022_heart_donorSelected_subtypes.added.rds'))


########################################################
########################################################
# Section IV: double check the subtypes  
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, '/Kuppe2022_heart_donorSelected_subtypes.added.rds'))

annot_file = list.files(path = paste0('../published_dataset/human/Kuppe_et_al_2022/',
                        'processed_data_Robj/subtype_annot'),
                        pattern = '*.Rds', full.names = TRUE)

##########################################
# manually add the batch information, those batch information were found in the subtype_annot files 
##########################################
aa$sample = droplevels(aa$sample)
aa$batch = 'B'
aa$batch[grep('CK1', aa$sample)] = 'A'

jj = which(aa$cell_type_original == 'Mast'| aa$cell_type_original == 'Cycling cells' |
             aa$cell_type_original == 'Adipocyte')

aa$annotation[jj] = as.character(aa$cell_type_original[jj])


p1 = DimPlot(aa, reduction = 'umap', group.by = 'cell_type_original', raster=FALSE, label = TRUE, repel = TRUE)
p2 = DimPlot(aa, reduction = 'umap', group.by = 'annotation', raster=FALSE, label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(resDir, '/umap_celltype_sutypes.pdf'), width = 16, height = 6)


##########################################
# ## double check CM
##########################################
source('functions_dataIntegration.R')

celltype = 'Cardiomyocyte'
#celltype = 'Fibroblast'

outDir = paste0(resDir, '/', celltype)
if(!dir.exists(outDir)) dir.create(outDir)
n = grep(celltype, annot_file)
cat(n, '--', basename(annot_file[n]), '-- cell type : ', celltype, '\n')

xx = readRDS(annot_file[n])
cat(nrow(xx), 'cells \n')
print(table(xx$annotation))

p1 = DimPlot(xx, group.by = 'annotation', reduction = 'umap')
p2 = DimPlot(xx, group.by = 'orig.ident', reduction = 'umap')
p3 = DimPlot(xx, group.by = 'patient', reduction = 'umap')
p4 = DimPlot(xx, group.by = 'batch', reduction = 'umap')

(p1 + p2)/(p3 + p4)

ggsave(filename = paste0(outDir, '/umap_subtypes_sample_patient_batch_', celltype, '.pdf'),
       width = 16, height = 12)


p1 = DimPlot(xx, group.by = 'annotation', reduction = 'umap_harmony')
p2 = DimPlot(xx, group.by = 'annotation', reduction = 'umap_harmony_v2')

p1 + p2

ggsave(filename = paste0(outDir, '/umap_clusters_celltype_originalIntegration_', celltype, '.pdf'),
       width = 16, height = 8)

rm(xx)

## code modified based on 
## https://github.com/saezlab/visium_heart/blob/master/sub-clustering/Cardiomyocyte/03_integrate_snRNA.ipynb
subs = subset(aa, cells = colnames(aa)[which(aa$cell_type_original == celltype)])

### we remove samples with less than 400 cells
df_cell_count <- as.data.frame(subs@meta.data) %>%
  group_by(sample) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

df_cell_count
df_cell_count <- subset(df_cell_count, count > 400)

subs <- subset(subs, sample %in% df_cell_count$sample)

DimPlot(subs, group.by = 'cell_type_original', raster=FALSE, label = TRUE, repel = TRUE)

subs  <-  subs %>% 
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, reduction = 'pca', n.neighbors = 30, min.dist = 0.1)

subs <- RunUMAP(subs, dims = 1:30, reduction = 'harmony', n.neighbors = 30, min.dist = 0.1,
                reduction.name = 'harmony_umap')

p1 = DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'annotation', raster=FALSE, reduction = 'umap')
p2 = DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'annotation', raster=FALSE, 
             reduction = 'harmony_umap')

p1 + p2 

ggsave(filename = paste0(outDir, '/umap_subtypes_umap.no.integration_original.harmony.umap_',
                         celltype, '.pdf'), width = 16, height = 8)


options(repr.plot.height = 5, repr.plot.width = 20)
p1 <- DimPlot(subs, reduction = "umap", group.by = "sample")
p2 <- DimPlot(subs, reduction = "umap", group.by = "patient_region_id")
p3 <- DimPlot(subs, reduction = "umap", group.by = "batch")
p4 <- DimPlot(subs, reduction = "umap", group.by = "patient_group", label = TRUE)

patchwork::wrap_plots(list(p1, p2, p3, p4), nrow = 2)

ggsave(filename = paste0(outDir, '/umap_samples_patient_region_batch_', celltype, '.pdf'), 
       width = 12, height = 8)

subs$sample = droplevels(subs$sample)
subs$patient_region_id = droplevels(subs$patient_region_id)
#subs$batch = droplevels(subs$batch)

subs <- harmony::RunHarmony(subs, 
                   group.by.vars = c("sample", "patient_region_id", "batch"),
                    reduction = "pca", 
                    max.iter.harmony = 30, 
                          dims.use = 1:30,
                          project.dim = FALSE,
                          plot_convergence = TRUE)

subs <- RunUMAP(subs, dims = 1:30, 
                reduction = 'harmony',
                reduction.name = "umap_harmony2",
                reduction.ke = 'umapharmony2_',
                verbose = FALSE,
                min.dist = 0.1)

p1 = DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'annotation', raster=FALSE, reduction = 'umap')
p2 = DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'annotation', raster=FALSE, 
             reduction = 'harmony_umap')
p3 = DimPlot(subs, label = TRUE, repel = TRUE, group.by = 'annotation', raster=FALSE, 
             reduction = 'umap_harmony2')

p1 + p2 + p3

ggsave(filename = paste0(outDir, '/umap_subtypes_umap.no.integration_original.harmony.umap_',
                         celltype, '.pdf'), width = 20, height = 8)

subs$condition = paste0(subs$sample, '_', subs$patient_region_id, '_', subs$batch)

ref.combined = IntegrateData_Seurat_RPCA(subs, group.by = 'condition', 
                                         redo.normalization.scaling = TRUE,
                                         correct.all = TRUE)

saveRDS(ref.combined, file = paste0(outDir, '/seuratObj_dataIntegration_seuratRPCA_', 
                                    celltype, '.rds'))

DimPlot(ref.combined, group.by = 'annotation', label = TRUE, repel = TRUE)

ggsave(filename = paste0(outDir, '/umap_subtypes_seuratRPCA_',
                         celltype, '.pdf'), width = 12, height = 8)


#ElbowPlot(subs, ndims = 50)
ref.combined <- FindNeighbors(ref.combined, dims = 1:20, reduction = 'pca')
ref.combined <- FindClusters(ref.combined, verbose = FALSE, algorithm = 3, resolution = 0.7)
ref.combined <- RunUMAP(ref.combined, dims = 1:30, reduction = 'pca', n.neighbors = 30, min.dist = 0.1,
                reduction.name = 'umap')

p1 = DimPlot(ref.combined, group.by = 'annotation', label = TRUE, repel = TRUE)
p2 = DimPlot(ref.combined, group.by = 'seurat_cluster', label = TRUE, repel = TRUE)

p1 + p2

ggsave(filename = paste0(outDir, '/umap_subtypes_seuratRPCA_clusters_',
                         celltype, '.pdf'), width = 16, height = 8)



