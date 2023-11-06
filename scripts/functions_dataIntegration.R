##########################################################################
##########################################################################
# Project: RA competence and also generic functions
# Script purpose: single cell data integration
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Sep 18 11:17:50 2023
##########################################################################
##########################################################################
IntegrateData_Seurat_CCA = function()
{
  ref.list <- SplitObject(refs.merged, split.by = "dataset")
  
  rm(list = c('refs.merged', 'aa', 'srat')) # remove big seurat objects to clear memory
  
  # normalize and identify variable features for each dataset independently
  ref.list <- lapply(X = ref.list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize")
    x <- FindVariableFeatures(x, selection.method = "vst")
    
  })
  
  # select features that are repeatedly variable across datasets for integration run PCA on each
  # dataset using these features
  features <- SelectIntegrationFeatures(object.list = ref.list)
  
  ref.list <- lapply(X = ref.list, FUN = function(x) {
    x <- ScaleData(x, features = features.common, verbose = TRUE)
    x <- RunPCA(x, features = features, verbose = FALSE)
    
  })
  
  ref.anchors <- FindIntegrationAnchors(object.list = ref.list, 
                                        anchor.features = features, 
                                        #reference = c(2),
                                        #reduction = "cca", 
                                        reduction = 'rpca',
                                        #k.anchor = 5,
                                        dims = 1:50)
  
  rm(ref.list)
  
  # this command creates an 'integrated' data assay
  ref.combined <- IntegrateData(anchorset = ref.anchors, features.to.integrate = features.common, 
                                dims = 1:50) ## take ~100G memory
  
  rm(ref.anchors)
  
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(ref.combined) <- "integrated"
  
  #xx = DietSeurat(ref.combined, counts = TRUE, data = TRUE, scale.data = TRUE, assays = 'integrated')
  #xx@assays$integrated@counts = ref.combined@assays$RNA@counts
  #saveRDS(xx, file = paste0(outDir, 
  #                          '/Seurat.obj_mouseGastrulation_mNT_integrated.rds'))
  
  # Run the standard workflow for visualization and clustering
  #ref.combined = readRDS(file =paste0(outDir, 
  #                                    '/Seurat.obj_mouseGastrulation_mNT_integrated.rds'))
  
  
  ref.combined <- ScaleData(ref.combined, verbose = FALSE)
  ref.combined <- RunPCA(ref.combined, npcs = 50, verbose = FALSE)
  
  ElbowPlot(ref.combined, ndims = 50)
  
  kk = which(ref.combined$dataset == 'mNT') 
  ref.combined$celltype[kk] = paste0('mNT_', ref.combined$condition[kk])
  
  
  #ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:20)
  #ref.combined <- FindClusters(ref.combined, resolution = 0.2)
  ref.combined <- RunUMAP(ref.combined, reduction = "pca", dims = 1:50, n.neighbors = 50, 
                          min.dist = 0.2) 
  
  #DimPlot(ref.combined, reduction = "umap")
  
  # Visualization
  names(cols_sel) = paste0('mNT_', names(cols_sel))
  
  DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE,
          repel = TRUE, raster=FALSE, cols = c(cols_mouse, cols_sel))
  
  ggsave(paste0(outDir, '/Integration_dataset_celltypes.pdf'), 
         width = 14, height = 8)
  
  DimPlot(ref.combined, reduction = "umap", group.by = "dataset", raster=FALSE)
  ggsave(paste0(outDir, '/Integration_dataset.pdf'), 
         width = 10, height = 8)
  
  pdf(paste0(outDir, '/FeaturePlot_Markers.pdf'),
      width =10, height = 8, useDingbats = FALSE)
  
  ggs = c('Pax6', 'Foxa2', 'Pou5f1', 'Sox17', 'Sox1', 'Sox2')
  for(n in 1:length(ggs))
  {
    p1 = FeaturePlot(ref.combined, features = ggs[n], min.cutoff = 'q5')
    #FeaturePlot(ref.combined, features = 'Foxa2', min.cutoff = 'q5')
    #FeaturePlot(ref.combined, features = 'Sox17', min.cutoff = 'q5')
    plot(p1)
  }
  
  dev.off()
  
  DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE, split.by = 'dataset',
          repel = TRUE, raster=FALSE) + NoLegend()
  
  ggsave(paste0(outDir, '/Integration_celltypes_split.dataset.pdf'), 
         width = 24, height = 8)
  
  DimPlot(ref.combined, reduction = "umap", group.by = "stage", label = TRUE,
          repel = TRUE, raster=FALSE)
  
  ggsave(paste0(outDir, '/Integration_stage.pdf'), 
         width = 16, height = 8)
  
}

IntegrateData_Seurat_RPCA = function()
{
  ref.list <- SplitObject(refs.merged, split.by = "dataset")
  
  rm(list = c('refs.merged', 'aa', 'srat')) # remove big seurat objects to clear memory
  
  # normalize and identify variable features for each dataset independently
  ref.list <- lapply(X = ref.list, FUN = function(x) {
    x <- NormalizeData(x, normalization.method = "LogNormalize")
    x <- FindVariableFeatures(x, selection.method = "vst")
    
  })
  
  # select features that are repeatedly variable across datasets for integration run PCA on each
  # dataset using these features
  features <- SelectIntegrationFeatures(object.list = ref.list)
  
  ref.list <- lapply(X = ref.list, FUN = function(x) {
    x <- ScaleData(x, features = features.common, verbose = TRUE)
    x <- RunPCA(x, features = features, verbose = FALSE)
    
  })
  
  ref.anchors <- FindIntegrationAnchors(object.list = ref.list, 
                                        anchor.features = features, 
                                        #reference = c(2),
                                        #reduction = "cca", 
                                        reduction = 'rpca',
                                        k.anchor = 5,
                                        dims = 1:50)
  
  rm(ref.list)
  
  # this command creates an 'integrated' data assay
  ref.combined <- IntegrateData(anchorset = ref.anchors, features.to.integrate = features.common, 
                                dims = 1:50) ## take ~100G memory
  
  rm(ref.anchors)
  
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(ref.combined) <- "integrated"
  
  #xx = DietSeurat(ref.combined, counts = TRUE, data = TRUE, scale.data = TRUE, assays = 'integrated')
  #xx@assays$integrated@counts = ref.combined@assays$RNA@counts
  #saveRDS(xx, file = paste0(outDir, 
  #                          '/Seurat.obj_mouseGastrulation_mNT_integrated.rds'))
  
  # Run the standard workflow for visualization and clustering
  #ref.combined = readRDS(file =paste0(outDir, 
  #                                    '/Seurat.obj_mouseGastrulation_mNT_integrated.rds'))
  
  
  ref.combined <- ScaleData(ref.combined, verbose = FALSE)
  ref.combined <- RunPCA(ref.combined, npcs = 50, verbose = FALSE)
  
  ElbowPlot(ref.combined, ndims = 50)
  
  kk = which(ref.combined$dataset == 'mNT') 
  ref.combined$celltype[kk] = paste0('mNT_', ref.combined$condition[kk])
  names(cols_sel) = paste0('mNT_', names(cols_sel))
  
  
  #ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:20)
  #ref.combined <- FindClusters(ref.combined, resolution = 0.2)
  ref.combined <- RunUMAP(ref.combined, reduction = "pca", dims = 1:50, n.neighbors = 50, 
                          min.dist = 0.2) 
  
  #DimPlot(ref.combined, reduction = "umap")
  
  saveRDS(ref.combined, file = paste0(outDir, '/integrated_mNT_mouseGastrulation_SeuratRPCA.rds'))
  
  cols_used = c(cols_mouse, cols_sel)
  saveRDS(cols_used, file = paste0(outDir, '/integrated_mNT_mouseGastrulation_colorsUsed.rds'))
  
  # Visualization
  DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE,
          repel = TRUE, raster=FALSE, cols = c(cols_mouse, cols_sel))
  
  ggsave(paste0(outDir, '/Integration_dataset_celltypes.pdf'), 
         width = 14, height = 8)
  
  DimPlot(ref.combined, reduction = "umap", group.by = "dataset", raster=FALSE)
  ggsave(paste0(outDir, '/Integration_dataset.pdf'), 
         width = 10, height = 8)
  
  pdf(paste0(outDir, '/FeaturePlot_Markers.pdf'),
      width =10, height = 8, useDingbats = FALSE)
  
  ggs = c('Pax6', 'Foxa2', 'Pou5f1', 'Sox17', 'Sox1', 'Sox2')
  for(n in 1:length(ggs))
  {
    p1 = FeaturePlot(ref.combined, features = ggs[n], min.cutoff = 'q5')
    #FeaturePlot(ref.combined, features = 'Foxa2', min.cutoff = 'q5')
    #FeaturePlot(ref.combined, features = 'Sox17', min.cutoff = 'q5')
    plot(p1)
  }
  
  dev.off()
  
  DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE, split.by = 'dataset',
          repel = TRUE, raster=FALSE) + NoLegend()
  
  ggsave(paste0(outDir, '/Integration_celltypes_split.dataset.pdf'), 
         width = 24, height = 8)
  
  DimPlot(ref.combined, reduction = "umap", group.by = "stage", label = TRUE,
          repel = TRUE, raster=FALSE)
  
  ggsave(paste0(outDir, '/Integration_stage.pdf'), 
         width = 16, height = 8)
  
}

IntegrateData_runHarmony = function(aa, ref)
{
  library(harmony)
  library(Seurat)
  library(SeuratData)
  
  # ref = srat;
  # refs.merged = merge(aa, y = ref, add.cell.ids = c("mNT", "mouseGastrulation"), project = "RA_competence")
  
  refs.merged <- NormalizeData(refs.merged) %>% FindVariableFeatures() %>% ScaleData() %>% 
    RunPCA(verbose = FALSE)
  
  kk = which(refs.merged$dataset == 'mNT') 
  refs.merged$celltype[kk] = paste0('mNT_', refs.merged$condition[kk])
  names(cols_sel) = paste0('mNT_', names(cols_sel))
  
  refs.merged <- RunHarmony(refs.merged, 
                            group.by.vars = "dataset",
                            reduction = 'pca',
                            dims.use = c(1:50), 
                            nclust = 20,
                            reference_values = 'ref',
                            epsilon.harmony = -Inf,
                            max.iter.harmony = 20, 
                            verbose = TRUE,
                            plot_convergence = TRUE
                            )
  
  refs.merged <- RunUMAP(refs.merged, reduction = "harmony", dims = 1:30)
  #refs.merged <- FindNeighbors(refs.merged, reduction = "harmony", dims = 1:30) %>% FindClusters()
  
  DimPlot(refs.merged, group.by = c("celltype"), label = TRUE, repel = TRUE,
          raster=FALSE, cols = c(cols_mouse, cols_sel))
  
  ggsave(paste0(outDir, '/Integration_dataset_celltypes.pdf'), 
         width = 14, height = 8)
  
  DimPlot(refs.merged, group.by = c("dataset"), label = TRUE)
  ggsave(paste0(outDir, '/Integration_dataset.pdf'), 
         width = 10, height = 8)
  
  
  pdf(paste0(outDir, '/FeaturePlot_Markers.pdf'),
      width =10, height = 8, useDingbats = FALSE)
  
  ggs = c('Pax6', 'Foxa2', 'Pou5f1', 'Sox17', 'Sox1', 'Sox2')
  for(n in 1:length(ggs))
  {
    p1 = FeaturePlot(refs.merged, features = ggs[n], min.cutoff = 'q5')
    #FeaturePlot(ref.combined, features = 'Foxa2', min.cutoff = 'q5')
    #FeaturePlot(ref.combined, features = 'Sox17', min.cutoff = 'q5')
    plot(p1)
  }
  
  dev.off()
  
  ## a quick test
  Test_if_RunHarmony_works = FALSLE
  if(Test_if_RunHarmony_works){
    InstallData("pbmcsca")
    data("pbmcsca")
    pbmcsca <- NormalizeData(pbmcsca) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
    pbmcsca <- RunHarmony(pbmcsca, group.by.vars = "Method")
    
    pbmcsca <- RunUMAP(pbmcsca, reduction = "harmony", dims = 1:30)
    pbmcsca <- FindNeighbors(pbmcsca, reduction = "harmony", dims = 1:30) %>% FindClusters()
    DimPlot(pbmcsca, group.by = c("Method", "ident", "CellType"), ncol = 3)
    
  }
  
  ##########################################
  #  ## directly calling Harmony
  ## original code 
  # from http://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/docs/advanced.html
  ##########################################
  Run_Harmony = FALSE
  if(Run_Harmony){
    
    refs.merged = NormalizeData(refs.merged, normalization.method = "LogNormalize")
    refs.merged <- FindVariableFeatures(refs.merged, selection.method = "vst")
    refs.merged <- ScaleData(refs.merged, verbose = TRUE)
    refs.merged <- RunPCA(refs.merged, verbose = FALSE)
    
   
    #V <- harmony::cell_lines$scaled_pcs
    #V_cos <- cosine_normalize(V, 1)
    #meta_data <- harmony::cell_lines$meta_data  
    V = refs.merged@reductions$pca@cell.embeddings
    meta_data = as_tibble(refs.merged@meta.data)
    
    set.seed(1)
    harmony_embeddings <- harmony::HarmonyMatrix(
      data_mat = V, ## PCA embedding matrix of cells
      meta_data = meta_data, ## dataframe with cell labels
      #vars_use = 'dataset',
      #theta = 0.5, ## cluster diversity enforcement
      vars_use = 'dataset', ## variable to integrate out
      #npcs = 30,
      #nclust = 5, ## number of clusters in Harmony model
      max.iter.harmony = 10, ## stop after initialization
      lambda=0.2,
      return_object = FALSE, ## return the full Harmony model object
      do_pca = FALSE ## don't recompute PCs
      
    )
    
    source('functions_utility.R')
    p1 <- do_scatter(harmony_embeddings, meta_data, label_name = 'dataset') + 
      labs(title = 'Colored by dataset')
    # p2 <- do_scatter(harmony_embeddings, meta_data, 'cell_type') + 
    #   labs(title = 'Colored by cell type')
    # cowplot::plot_grid(p1, p2, nrow = 1)
    
    refs.merged[['harmony']] = Seurat::CreateDimReducObject(embeddings=harmony_embeddings,
                                                            key='HARMONY_',
                                                            assay='RNA')
    
    ref.combined <- RunUMAP(refs.merged, reduction = "pca", dims = 1:50, n.neighbors = 50, 
                            min.dist = 0.2) 
    
    kk = which(ref.combined$dataset == 'mNT') 
    ref.combined$celltype[kk] = paste0('mNT_', ref.combined$condition[kk])
    
    # Visualization
    names(cols_sel) = paste0('mNT_', names(cols_sel))
    
    DimPlot(ref.combined, reduction = "umap", group.by = "celltype", label = TRUE,
            repel = TRUE, raster=FALSE, cols = c(cols_mouse, cols_sel))
    
    
  }
  
}

IntegrateData_runFastMNN = function()
{
  library(Seurat)
  library(SeuratData)
  library(SeuratWrappers)
 
  # ref = srat;
  refs.merged = merge(aa, y = ref, add.cell.ids = c("mNT", "mouseGastrulation"), project = "RA_competence")
  
  refs.merged <- NormalizeData(refs.merged)
  refs.merged <- FindVariableFeatures(refs.merged)
  refs.merged <- RunFastMNN(object.list = SplitObject(refs.merged, split.by = "dataset"))
  
  refs.merged <- RunUMAP(refs.merged, reduction = "mnn", dims = 1:30)
  
  kk = which(refs.merged$dataset == 'mNT')
  refs.merged$celltype[kk] = paste0('mNT_', refs.merged$condition[kk])
  #cat(cols_sel)
  #names(cols_sel) = paste0('mNT_', names(cols_sel))
  
  
  DimPlot(refs.merged, group.by = c("celltype"), label = TRUE, repel = TRUE,
          raster=FALSE, cols = c(cols_mouse, cols_sel))
  ggsave(paste0(outDir, '/Integration_dataset_celltypes.pdf'), 
         width = 14, height = 8)
  
  DimPlot(refs.merged, group.by = c("dataset"), label = TRUE)
  ggsave(paste0(outDir, '/Integration_dataset.pdf'), 
         width = 10, height = 8)
  
  
  pdf(paste0(outDir, '/FeaturePlot_Markers.pdf'),
      width =10, height = 8, useDingbats = FALSE)
  
  ggs = c('Pax6', 'Foxa2', 'Pou5f1', 'Sox17', 'Sox1', 'Sox2')
  for(n in 1:length(ggs))
  {
    p1 = FeaturePlot(refs.merged, features = ggs[n], min.cutoff = 'q5')
    #FeaturePlot(ref.combined, features = 'Foxa2', min.cutoff = 'q5')
    #FeaturePlot(ref.combined, features = 'Sox17', min.cutoff = 'q5')
    plot(p1)
  }
  
  dev.off()
  
}

IntegrateData_runSCVI = function()
{
  library(SeuratDisk)
  
  mnt = refs.merged
  mnt <- NormalizeData(mnt)
  
  mnt <- FindVariableFeatures(mnt, nfeatures = 5000)
  mnt = subset(mnt, features = VariableFeatures(mnt))
  
  VariableFeatures(mnt) = NULL
  #mnt@assays$RNA@scale.data = NULL
  #mnt@assays$RNA@data = NULL
  
  DefaultAssay(mnt) = 'RNA'
  mnt = DietSeurat(mnt, counts = TRUE, data = TRUE,
                   scale.data = FALSE,
                   features = rownames(mnt), 
                   assays = c('RNA'), 
                   dimreducs = NULL, graphs = NULL, 
                   misc = TRUE
  )
  
  DefaultAssay(mnt) = 'RNA'
  VariableFeatures(mnt)
  
  #Idents(mnt) = mnt$condition
  #mnt = subset(mnt, downsample = 1500)
  
  saveFile = '/RNAmatrix_mouseGastrulation_mNT_HVG5k.h5Seurat'
  
  SaveH5Seurat(mnt, filename = paste0(outDir, saveFile), 
               overwrite = TRUE)
  Convert(paste0(outDir, saveFile), 
          dest = "h5ad", overwrite = TRUE)
  
}

##########################################
# two test functions
# - test to combine the adult heart datasets without integration
# - test to convert the batch corrected expression to counts (not good idea)
#  i.e. try to reversely calculated batch-corrected UMI counts using corrected gene expression matrix from Seurat 
##########################################
Combine.adult.mice.heart.without.integration = function(refs)
{
  # test refs without integration
  Test_refs_withoutIntegration = FALSE
  if(Test_refs_){
    DefaultAssay(refs) <- "RNA"
    refs = FindVariableFeatures(refs, selection.method = "vst", nfeatures = 3000)
    
    # Run the standard workflow for visualization and clustering
    refs <- ScaleData(refs, verbose = FALSE)
    refs <- RunPCA(refs, npcs = 30, verbose = FALSE)
    
    ElbowPlot(refs, ndims = 30)
    
    refs <- FindNeighbors(refs, reduction = "pca", dims = 1:20)
    refs <- FindClusters(refs, resolution = 0.5)
    
    refs <- RunUMAP(refs, reduction = "pca", dims = 1:30, n.neighbors = 50, min.dist = 0.05) 
    
    # Visualization
    p1 <- DimPlot(refs, reduction = "umap", group.by = "dataset")
    p2 <- DimPlot(refs, reduction = "umap", group.by = "annot.ref", label = TRUE,
                  repel = TRUE)
    p1 + p2 + ggsave(paste0(resDir, '/Forte2020_Ren2020_noCorrection_', Normalization, '.pdf'), 
                     width = 24, height = 10)
    
  }
  
}

Convert.batch.corrected.expression.matrix.to.UMIcount = function(refs){
  refs = readRDS(file = paste0(RdataDir, 
                               'SeuratObj_adultMiceHeart_refCombine_Forte2020.nonCM_Ren2020CM_cleanAnnot_logNormalize_v1.rds'))
  
  #jj = which(metadata$dataset == 'Ren2020')
  #aa = readRDS(file = paste0(RdataDir, 'Forte2020_logNormalize_allgenes.rds'))
  #cms = readRDS(file =  paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization_umap.rds'))
  p1 <- DimPlot(refs, reduction = "umap", group.by = "dataset")
  p2 <- DimPlot(refs, reduction = "umap", group.by = "celltype", label = TRUE,
                repel = TRUE)
  p1 + p2
  
  ggsave(filename = paste0(resDir, '/AdultMice_scRNAref_overView.pdf'),  width = 20, height = 8)
  
  VlnPlot(refs, features = c("nCount_RNA"), group.by = 'dataset')
  
  metadata = refs@meta.data   
  Ec = refs@assays$integrated@data
  
  counts = refs@assays$RNA@counts
  counts = counts[match(rownames(Ec), rownames(counts)), ]
  #E = refs@assays$RNA@data
  
  ss = colSums(counts)
  ccts = expm1(Ec)/10000
  ccts = t(t(ccts)*ss)
  ccts = round(ccts)
  
  metadata$nCount_RNA = ss
  
  aa <- CreateSeuratObject(counts = ccts, project = "adult", min.cells = 50, min.features = 500)
  aa = AddMetaData(aa, metadata, col.name = NULL) 
  
  aa <- NormalizeData(aa, normalization.method = "LogNormalize", scale.factor = 10000)
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
  
  plot1 <- VariableFeaturePlot(aa)
  
  top10 <- head(VariableFeatures(aa), 10)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot1 + plot2
  
  all.genes <- rownames(aa)
  aa <- ScaleData(aa, features = all.genes)
  
  aa <- RunPCA(aa, verbose = FALSE, weight.by.var = TRUE)
  ElbowPlot(aa, ndims = 30)
  
  aa <- FindNeighbors(aa, dims = 1:10)
  aa <- FindClusters(aa, verbose = FALSE, algorithm = 3, resolution = 0.5)
  
  aa = RunUMAP(aa, reduction = "pca", dims = 1:30, n.neighbors = 30, min.dist = 0.05) 
  
  # Visualization
  p1 <- DimPlot(aa, reduction = "umap", group.by = "dataset")
  p2 <- DimPlot(aa, reduction = "umap", group.by = "celltype", label = TRUE,
                repel = TRUE)
  p1 + p2 + ggsave(paste0(resDir, '/refCombined_correctUMIcounts_Forte2020_Ren2020__', Normalization, '.pdf'), 
                   width = 24, height = 10)
  
  refs = aa
  #aa = readRDS(file = paste0(RdataDir, 'Seurat.obj_adultMiceHeart_week0.week2_Ren2020_seuratNormalization.rds'))
  
  saveRDS(refs, file = paste0(RdataDir, 
                              'SeuratObj_adultMiceHeart_refCombine_Forte2020.nonCM_Ren2020CM_cleanAnnot_correctedUMIcounts_v1.rds'))
  
}


