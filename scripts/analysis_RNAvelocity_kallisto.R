##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: prepare files for RNA velocity using nf-kallisto-velocity 
# https://github.com/csoneson/rna_velocity_quant
# https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Oct  7 14:09:17 2022
##########################################################################
##########################################################################
outDir = paste0(resDir, '/RNA_velocity_kallisto/')
system(paste0('mkdir -p ', outDir))

##########################################
# test kallisto annotation using BUSpaRse
##########################################
# test kallisto annotation file making 
if(make_annotationFile_kallisto){
  library(BUSpaRse)
  library(Seurat)
  #library(SeuratWrappers)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(AnnotationHub)
  library(zeallot) # For %<-% that unpacks lists in the Python manner
  library(DropletUtils)
  library(tidyverse)
  library(GGally) # For ggpairs
  library(velocyto.R)
  library(SingleR)
  library(scales)
  library(plotly)
  theme_set(theme_bw())
  
  ah <- AnnotationHub()
  query(ah, pattern = c("Ensembl", "97", "Mus musculus", "EnsDb"))
  
  # Get mouse Ensembl 97 annotation
  edb <- ah[["AH73905"]]
  
  require("ensembldb")
  
  get_velocity_files(edb, L = 91, Genome = BSgenome.Mmusculus.UCSC.mm10, 
                     out_path = paste0("/groups/tanaka/People/current/jiwang/Genomes/",
                                       "axolotl/Transcriptomics/kallisto_RNAvelocity_index/test"), 
                     isoform_action = "separate")
  
}

########################################################
########################################################
# Section : process the kallisto output and prepare the data for scvelo and cellrank
# 
########################################################
########################################################
aa =  readRDS(file = paste0("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/", 
                            "CM_subset_for_velocity.rds"))
aa$time = gsub('Amex_scRNA_', '', aa$condition)
aa$cell.ids = sapply(colnames(aa), function(x) unlist(strsplit(as.character(x), '-'))[1]) 
aa$cell.ids = paste0(aa$cell.ids, '_', aa$time)


levels = c("Amex_scRNA_d0", "Amex_scRNA_d1",
           "Amex_scRNA_d4", "Amex_scRNA_d7", 
           "Amex_scRNA_d14")
Idents(aa) = factor(aa$condition, levels = levels)

#DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

CMmyLevels <- c(  "CM_ven_(Robo2)", "CM_ven_(Cav3_1)", "CM_IS","CM_Prol_IS", "CM_Prol_1", "CM_Prol_3")
aa$subtypes = factor(aa$subtypes, levels = CMmyLevels)

cols = c("#4CC9F0", "#49A2F0",
         "#941F56", "#C61010",  
         "#4361EE",  "#3F37C9")

DimPlot(aa, dims = c(1, 2), label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE,
        cols = cols
)

ggsave(filename = paste0(outDir, 'UMAP_CMsubsets_forRNAvelocity.pdf'), width = 10, height = 8)

search_optimal_umap_parameters = FALSE
if(search_optimal_umap_parameters){
  source("functions_scRNAseq.R")
  explore.umap.params.combination(sub.obj = aa, resDir = outDir, 
                                  pdfname = 'axolotl_snRNA_RNAvelocity_umap_test.pdf',
                                  use.parallelization = FALSE,
                                  group.by = 'condition',
                                  nfeatures.sampling = c(1000, 2000, 3000, 5000),
                                  nb.pcs.sampling = c(10, 30, 50, 100), 
                                  n.neighbors.sampling = c(20, 30, 50, 100, 200),
                                  min.dist.sampling = c(0.1, 0.3)
                                  
  )
  
}



xx = aa
xx@reductions$umap@cell.embeddings = -xx@reductions$umap@cell.embeddings
DimPlot(xx, dims = c(1, 2), label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE,
        cols = cols
)

aa = xx

DimPlot(aa, dims = c(1, 2), label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE,
        cols = cols
)

ggsave(filename = paste0(outDir, 'UMAP_CMsubsets_forRNAvelocity.pdf'), width = 10, height = 8)

saveRDS(aa, file = paste0(RdataDir, 'CM_subset_for_velocity.rds'))

## this is how the umap was calculated
# CM_subset_2 = aa
# CM_subset_2 <- FindVariableFeatures(CM_subset_2, selection.method = "vst", nfeatures = 2000)
# 
# all.genes <- rownames(CM_subset_2)
# CM_subset_2 <- ScaleData(CM_subset_2, features = all.genes)
# 
# CM_subset_2 <- RunPCA(CM_subset_2, features = VariableFeatures(object = CM_subset_2))
# 
# CM_subset_2 <- FindNeighbors(CM_subset_2, dims = 1:10)
# 
# 
# CM_subset_2 <- FindClusters(CM_subset_2, resolution = 0.3)
# CM_subset_2$subtypes -> Idents(CM_subset_2)
# 
# CM_subset_2 <- RunUMAP(CM_subset_2, dims = 1:10)
# 
# DimPlot(CM_subset_2, dims = c(1, 2), label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE,
#         cols = cols
# ) 

##########################################
# Import abundances into R with tximeta
# https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/
##########################################
library(Seurat)
library(Seurat)
library(DropletUtils)
library(edgeR)
#library(BiocParallel)

source('functions_scRNAseq.R')

dataDir = "../R13591_axolotl_multiome/kallisto_velocity"
design = data.frame(sampleID = seq(197249, 197253), 
                    condition = c(paste0('Amex_scRNA_d', c(0, 1, 4, 7, 14))), 
                    stringsAsFactors = FALSE)

design$time = gsub('Amex_scRNA_', '', design$condition)

for(n in 1:nrow(design))
{
  # n = 1
  cat('-----------',n, ' : ', design$condition[n], '-------------\n')
  
  # load nf output and process
  topdir = paste0(dataDir, '/', design$condition[n], '/')
  
  # aa = make_SeuratObj_scRNAseq(topdir = topdir,
  #                              saveDir = paste0(resDir, '/', design$condition[n], '_', design$sampleID[n], '/'), 
  #                              changeGeneName.axolotl = TRUE, 
  #                              defaultDrops.only = TRUE)
  
  # unspliced and spliced 
  for(obj in c('unspliced', 'spliced')){
    # obj = 'spliced'
    exp = Matrix::readMM(paste0(topdir, obj, ".mtx")) #read matrix
    bc = read.csv(paste0(topdir, obj,  ".barcodes.txt"), header = F, stringsAsFactors = F)
    g = read.csv(paste0(topdir, obj, ".genes.txt"), header = F, stringsAsFactors = F)
    
    cat('change gene names \n')
    # change the gene names before making Seurat object
    annot = readRDS(paste0('/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                           'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_',
                           'curated.geneSymbol.toUse.rds'))
    
    mm = match(g$V1, annot$geneID)
    ggs = paste0(annot$gene.symbol.toUse[mm], '-',  annot$geneID[mm])
    g$V1[!is.na(mm)] = ggs[!is.na(mm)]
    
    dimnames(exp) = list(bc$V1, g$V1) # number added because of seurat format for barcodes
    count.data = Matrix::t(exp)
    rm(exp)
    
    # select only cells found in aa
    meta = data.frame(cell.id = paste0(colnames(count.data), '_', design$time[n]),
                      condition = design$condition[n])
    cell2keep = !is.na(match(meta$cell.id, aa$cell.ids)) 
    meta$cell2keep = cell2keep
    rownames(meta) = colnames(count.data)
    
    mm = match(rownames(count.data), rownames(aa))
    gene2keep = which(!is.na(mm))
    
    srat = CreateSeuratObject(counts = count.data[gene2keep, cell2keep],
                              meta.data = meta[cell2keep, ], 
                              min.cells = 5, 
                              min.features = 10)
    
    if(n == 1) {
      if(obj == 'spliced'){
        spliced = srat
      }else{
        unspliced = srat
      }
      
    }else{
      if(obj == 'spliced'){
        spliced = merge(spliced, srat)
      }else{
        unspliced = merge(unspliced, srat)
      }
    }
  }
}


save(design, spliced, unspliced, 
     file = paste0(RdataDir, 'seuratObject_spliced_unspliced_', species, version.analysis, '.Rdata'))

rm(srat)
rm(count.data)
rm(meta)
rm(g)


########################################################
########################################################
# Section : prepare the files for scVelo and cellRank
# 
########################################################
########################################################
aa = readRDS(file = paste0(RdataDir, 'CM_subset_for_velocity.rds'))
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE)
aa$celltypes = droplevels(aa$subtypes)
aa <- FindNeighbors(aa, dims = 1:20)
aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
aa$time = gsub('d', '', aa$time)

subsetting_further = FALSE
if(subsetting_further){
  
  #aa = subset(aa, cells = colnames(aa)[which(aa$subtypes == 'CM_IS'|aa$subtypes == "CM_Prol_IS")])
  
  ### select 4 major cell types
  #aa = subset(aa, cells = colnames(aa)[which(aa$subtypes != 'CM_Prol_1' & aa$subtypes != "CM_Prol_3" &
  #                                             aa$condition != "Amex_scRNA_d14")])
  aa = subset(aa, cells = colnames(aa)[which(aa$subtypes != 'CM_Prol_1' & aa$subtypes != "CM_Prol_3")])
  
  aa$celltypes = droplevels(aa$subtypes)
  
  # aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000)
  # aa <- ScaleData(aa)
  # aa <- RunPCA(aa, features = VariableFeatures(object = aa), weight.by.var = TRUE, verbose = FALSE)
  # ElbowPlot(aa, ndims = 50)
  # 
  # aa <- RunUMAP(aa, dims = 1:10, n.neighbors = 20, min.dist = 0.3)
  # DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  # 
  # p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE) 
  # p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  # 
  # p1 + p2
  
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
  aa <- ScaleData(aa)
  
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), weight.by.var = TRUE, verbose = FALSE)
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 20, min.dist = 0.3)
  #DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE) 
  p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  p1 + p2
  
  aa <- FindNeighbors(aa, dims = 1:20)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
  p3 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  
  p1 + p2+ p3
  
  aa$cluster = aa$seurat_clusters
  
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  
  
  regressOut_cellcycle = FALSE
  if(regressOut_cellcycle){
    # load s.genes and g2m.genes
    load("/groups/tanaka/People/current/Paco/Collaborations/Elad_Paco/CellCycleGenes.RData")
    
    geneIDs = get_geneID(rownames(aa))
    mm = match(geneIDs, s.genes)
    s.genes = rownames(aa)[which(!is.na(mm))]
    mm = match(geneIDs, g2m.genes)
    g2m.genes = rownames(aa)[which(!is.na(mm))]
    
    ggs = c(s.genes, g2m.genes)
    
    reRun.cellcycle_regression = FALSE
    if(reRun.cellcycle_regression){
      write.csv(ggs, file = paste0(outDir, "cellcycle_gene.csv"), quote = FALSE, row.names = FALSE,
                col.names = FALSE)
      
      aa <- CellCycleScoring(aa, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
      
      # view cell cycle scores and phase assignments
      head(aa[[]])
      
      DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
      
      
      aa <- ScaleData(aa, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(aa))
      
      saveRDS(aa,  file = paste0(RdataDir, 'CM_subset_4majoySubtypes_regressedOut.cellcycles.rds'))
      
    }
   
    aa = readRDS(file = paste0(RdataDir, 'CM_subset_4majoySubtypes_regressedOut.cellcycles.rds'))
    
    ## remove the cell cycle genes here
    aa = subset(aa, features = setdiff(rownames(aa), ggs))
    ## rerun the PCA and umap
    aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 2000)
    #aa <- ScaleData(aa)
    
    aa <- RunPCA(aa, features = VariableFeatures(object = aa), weight.by.var = TRUE, verbose = FALSE)
    aa <- RunUMAP(aa, dims = 1:30, n.neighbors = 20, min.dist = 0.3)
    #DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
    
    p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE) 
    p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'Phase', raster=FALSE)
    
    p1 + p2
    
    aa <- FindNeighbors(aa, dims = 1:20)
    aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
    p3 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
    
    p1 + p2+ p3
    
    aa$cluster = aa$seurat_clusters
    
    DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
    
    
  }
  
  Serach_for_good_umap_parameters = FALSE
  if(search_optimal_umap_parameters){
    source("functions_scRNAseq.R")
    explore.umap.params.combination(sub.obj = aa, resDir = outDir, 
                                    pdfname = 'axolotl_CMsubsets_RNAvelocity_umap_test.pdf',
                                    use.parallelization = FALSE,
                                    group.by = 'condition',
                                    nfeatures.sampling = c(1000, 2000, 3000, 5000),
                                    nb.pcs.sampling = c(20, 30, 50, 100), 
                                    n.neighbors.sampling = c(20, 30, 50, 100, 200),
                                    min.dist.sampling = c(0.1, 0.3)
                                    
    )
    
    aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 5000)
    aa <- ScaleData(aa)
    
    aa <- RunPCA(aa, features = VariableFeatures(object = aa), weight.by.var = TRUE, verbose = FALSE)
    aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 30, min.dist = 0.3)
    DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
    
    p1 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE) 
    p2 = DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
    
    p1 + p2
    
    
  }
  
  Test_Batch_correction_timepoint = FALSE
  if(Test_Batch_correction_timepoint){
    aa.list <- SplitObject(aa, split.by = "condition")
    
    aa.list <- lapply(X = aa.list, FUN = function(x) {
      x<- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
    })
    
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = aa.list)
    
    Use.softer.rpca = FALSE
    if(Use.softer.rpca){
      aa.list <- lapply(X = aa.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
      })
      
      aa.anchors <- FindIntegrationAnchors(object.list = aa.list, anchor.features = features, 
                                           reduction = "rpca")
    }else{
      aa.anchors <- FindIntegrationAnchors(object.list = aa.list, anchor.features = features)
    }
    
    aa.combined <- IntegrateData(anchorset = aa.anchors)
    
    DefaultAssay(aa.combined) <- "integrated"
    
    # Run the standard workflow for visualization and clustering
    aa.combined <- ScaleData(aa.combined, verbose = FALSE)
    aa.combined <- RunPCA(aa.combined, npcs = 50, verbose = FALSE)
    
    aa.combined$condition = factor(aa.combined$condition, levels = levels(aa$condition))
    
    p1 = DimPlot(aa.combined, label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE) 
    p2 = DimPlot(aa.combined, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
    p1 + p2
    
    
    aa.combined <- FindNeighbors(aa.combined, reduction = "pca", dims = 1:30)
    aa.combined <- FindClusters(aa.combined, resolution = 0.5)
    
    p3 = DimPlot(aa.combined, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
    aa.combined$cluster = aa.combined$seurat_clusters
    
    aa.combined <- RunUMAP(aa.combined, reduction = "pca", dims = 1:30, n.neighbors = 50, min.dist = 0.1)
    p2 = DimPlot(aa.combined, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
    p1 + p2+ p3
    
   
    
  }
  
}

##########################################
# load saved spliced and unspliced 
# and select the cells and genes to use
##########################################
load(file = paste0(RdataDir, 'seuratObject_spliced_unspliced_', species, version.analysis, '.Rdata'))
aa$cell.id = aa$cell.ids
aa$celltypes = aa$subtypes
DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE) 

CCA_batch = FALSE
if(CCA_batch){
  
  aa$cluster = aa.combined$cluster[match(colnames(aa), colnames(aa.combined))]
  embedding = aa.combined@reductions$umap@cell.embeddings
  aa[['umap']] = Seurat::CreateDimReducObject(embeddings=embedding, key='UMAP_', assay='RNA')
  aa[['pca']] = Seurat::CreateDimReducObject(embeddings=aa.combined@reductions$pca@cell.embeddings, 
                                           key='PC_', assay='RNA')
  
}

# swapping the unspliced and spliced matrix
#xx = spliced
#spliced = unspliced
#unspliced = xx 

## select features shares by spliced and unspliced
features = intersect(rownames(spliced), rownames(unspliced))
spliced = subset(spliced, features = features)
unspliced = subset(unspliced, features = features)

# subsetting aa with the same cell as spliced and unspliced
cells.shared = intersect(spliced$cell.id, aa$cell.id)
mnt = subset(aa, cells = colnames(aa)[match(cells.shared, aa$cell.id)])

counts = GetAssayData(spliced, slot = 'counts')
counts = counts[, match(mnt$cell.id, spliced$cell.id)]
colnames(counts) = colnames(mnt)

# subsetting mnt with the same genes as spliced and unspliced
DefaultAssay(mnt) = 'RNA'
mnt = subset(mnt, features = intersect(rownames(counts), rownames(mnt)))

counts = counts[match(rownames(mnt), rownames(counts)), ]
mnt[["spliced"]] <- CreateAssayObject(counts = counts)
rm(spliced)

counts = GetAssayData(unspliced, slot = 'counts')
counts = counts[, match(mnt$cell.id, unspliced$cell.id)]
colnames(counts) = colnames(mnt)
counts = counts[match(rownames(mnt), rownames(counts)), ]
mnt[["unspliced"]]<- CreateAssayObject(counts = counts)
rm(unspliced)

# try to save mulitple assays 
# https://github.com/mojaveazure/seurat-disk/issues/21

#library(SeuratData)
library(SeuratDisk)

VariableFeatures(mnt) = NULL
#mnt@assays$RNA@scale.data = NULL
#mnt@assays$RNA@data = NULL

DefaultAssay(mnt) = 'RNA'
mnt = DietSeurat(mnt, counts = TRUE, data = TRUE,
                 scale.data = FALSE,
                 features = rownames(mnt), 
                 assays = c('RNA', 'spliced', 'unspliced'), 
                 dimreducs = c('umap'), graphs = NULL, 
                 misc = TRUE
)


DefaultAssay(mnt) = 'RNA'
VariableFeatures(mnt)

Idents(mnt) = mnt$condition
mnt$condition = as.character(mnt$condition)
mnt$celltypes = as.character(mnt$subtypes)

mnt$celltypes[which(mnt$celltypes == "CM_ven_(Robo2)")] = "CM_ven_Robo2"
mnt$celltypes[which(mnt$celltypes == "CM_ven_(Cav3_1)")] = "CM_ven_Cav3_1"

#mnt = subset(mnt, downsample = 2000)

#saveDir = paste0("/Volumes/groups/tanaka/People/current/jiwang/projects/RA_competence/",
#                "results/scRNAseq_R13547_10x_mNT_20220813/RA_symetryBreaking/")
#saveFile = "RNAmatrix_umap_kalisto.velocity_spliced_unspliced_CMsutypes_v1.5.h5Seurat"
saveFile = 'RNAmatrix_umap_kalisto.velocity_spliced_unspliced_CMsutypesAll_2.1.h5Seurat'

SaveH5Seurat(mnt, filename = paste0(outDir, saveFile), 
             overwrite = TRUE)
Convert(paste0(outDir, saveFile), 
        dest = "h5ad", overwrite = TRUE)


##########################################
## test umap with spliced 
##########################################
Test_umap_use.only.splicedMatrix = FALSE
if(Test_umap_use.only.splicedMatrix){
  DefaultAssay(mnt) = 'spliced'
  DimPlot(mnt, label = TRUE, repel = TRUE, reduction = 'UMAP',
          group.by = 'condition', cols = cols_sel, raster=FALSE)
  
  
  mnt <- FindVariableFeatures(mnt, selection.method = "vst", nfeatures = 3000) # find subset-specific HVGs
  
  ## because the data was regressed and scaled already, only the HVGs were used to calculate PCA
  mnt = ScaleData(mnt)
  
  mnt <- RunPCA(mnt, features = VariableFeatures(object = mnt), verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(mnt, ndims = 50)
  
  Idents(mnt) = mnt$condition
  
  mnt <- RunUMAP(mnt, dims = 1:20, n.neighbors = 100, min.dist = 0.2)
  DimPlot(mnt, label = TRUE, repel = TRUE, group.by = 'condition',
          reduction = 'umap', cols = cols_sel, raster=FALSE)
  
  ggsave(filename = paste0(outDir, 'UMAP_splicedMatrix_',
                           'subsetting.RAsymmetryBreaking.onlyday3rep1.pdf'), width = 10, height = 8)
  
  DimPlot(mnt, label = TRUE, repel = TRUE, reduction = 'UMAP',
          group.by = 'condition', cols = cols_sel, raster=FALSE)
  
}

### double check the intron and exon ratios
xx = readRDS(file = "../results/Rdata/seuratObject_axloltl_scRNAseq_R13591_20220720_lognormamlized_pca_umap_v2.rds")


