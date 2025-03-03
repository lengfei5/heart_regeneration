##########################################################################
##########################################################################
# Project: RA competence and also generic functions
# Script purpose: single cell data integration
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Sep 18 11:17:50 2023
##########################################################################
##########################################################################
IntegrateData_Seurat_CCA = function(seuratObj, group.by = 'dataset', nfeatures = 2000,
                                    ndims = c(1:30), k.weight = 100,
                                    merge.order = NULL,
                                    redo.normalization.scaling = TRUE,
                                    correct.all = FALSE,
                                    reference = NULL)
{
  # seuratObj = aa; group.by = "condition"; ndims = c(1:30);reference = NULL;nfeatures = 2000
  ref.list <- SplitObject(seuratObj, split.by = group.by)
  
  # normalize and identify variable features for each dataset independently
  if(redo.normalization.scaling){
    ref.list <- lapply(X = ref.list, FUN = function(x) {
      x <- NormalizeData(x, normalization.method = "LogNormalize")
      x <- FindVariableFeatures(x, selection.method = "vst")
    })
    
    # select features that are repeatedly variable across datasets for integration run PCA on each
    # dataset using these features
    features <- SelectIntegrationFeatures(object.list = ref.list, nfeatures = nfeatures)
    
    ref.list <- lapply(X = ref.list, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = TRUE)
      x <- RunPCA(x, features = features, verbose = FALSE)
      
    })
    
  }else{
    ref.list <- lapply(X = ref.list, FUN = function(x) {
      x <- FindVariableFeatures(x, selection.method = "vst")
    })
    
    # select features that are repeatedly variable across datasets for integration run PCA on each
    # dataset using these features
    features <- SelectIntegrationFeatures(object.list = ref.list, nfeatures = nfeatures)
    
    ref.list <- lapply(X = ref.list, FUN = function(x) {
      x <- RunPCA(x, features = features, verbose = FALSE)
      
    })
    
  }
  
  
  ref.anchors <- FindIntegrationAnchors(object.list = ref.list, 
                                        anchor.features = features, 
                                        reference = reference,
                                        reduction = "cca",
                                        #k.anchor = 5,
                                        dims = ndims)
  
  # this command creates an 'integrated' data assay
  if(correct.all){
    feature2integrate = rownames(seuratObj)
  }else{
    feature2integrate = features
  }
  
  if(is.null(merge.order)){
    ref.combined <- IntegrateData(anchorset = ref.anchors, 
                                  features.to.integrate = feature2integrate, 
                                  dims = ndims, 
                                  k.weight = k.weight) 
  }else{
    ref.combined <- IntegrateData(anchorset = ref.anchors, 
                                  features.to.integrate = feature2integrate,
                                  sample.tree = merge.order,
                                  dims = ndims, k.weight = k.weight)
  }
  
  
  rm(ref.list)
  rm(ref.anchors)
  
  # specify that we will perform downstream analysis on the corrected data note that the
  # original unmodified data still resides in the 'RNA' assay
  DefaultAssay(ref.combined) <- "integrated"
  
  ref.combined <- ScaleData(ref.combined, verbose = FALSE)
  ref.combined <- RunPCA(ref.combined, npcs = 50, verbose = FALSE)
  
  ElbowPlot(ref.combined, ndims = 50)
  
  ref.combined <- RunUMAP(ref.combined, reduction = "pca", dims = 1:30, n.neighbors = 30, 
                          min.dist = 0.3) 
  DimPlot(ref.combined, reduction = "umap", group.by = group.by, label = TRUE,
          repel = TRUE, raster=FALSE)
  
  #ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:20)
  #ref.combined <- FindClusters(ref.combined, resolution = 0.2)
  #DimPlot(ref.combined, reduction = "umap")
  
  return(ref.combined)
  
}

IntegrateData_Seurat_RPCA = function(seuratObj, group.by = 'dataset', nfeatures = 2000,
                                     ndims = c(1:30), k.weight = 100,
                                     merge.order = NULL,
                                     redo.normalization.scaling = TRUE,
                                     correct.all = FALSE,
                                     reference = NULL,
                                     use.parallelization = FALSE)
{
  
  if(use.parallelization){
    library(future)
    cat('test if multicore is supported or not : \n')
    print(supportsMulticore())
    
    #cat('use future package to paralle \n')
    
    plan()
    # change the current plan to access parallelization
    plan("multicore", workers = 16)
    plan()
    
  }
  
  # seuratObj = aa; group.by = "species"; ndims = c(1:30); reference = NULL;nfeatures = 2000
  ref.list <- SplitObject(seuratObj, split.by = group.by)
  
  # normalize and identify variable features for each dataset independently
  if(redo.normalization.scaling){
    ref.list <- lapply(X = ref.list, FUN = function(x) {
      x <- NormalizeData(x, normalization.method = "LogNormalize")
      x <- FindVariableFeatures(x, selection.method = "vst")
    })
    
    # select features that are repeatedly variable across datasets for integration run PCA on each
    # dataset using these features
    features <- SelectIntegrationFeatures(object.list = ref.list, nfeatures = nfeatures)
    
    ref.list <- lapply(X = ref.list, FUN = function(x) {
      x <- ScaleData(x, features = features, verbose = TRUE)
      x <- RunPCA(x, features = features, verbose = FALSE)
      
    })
  }else{
    ref.list <- lapply(X = ref.list, FUN = function(x) {
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures)
    })
    
    # select features that are repeatedly variable across datasets for integration run PCA on each
    # dataset using these features
    features <- SelectIntegrationFeatures(object.list = ref.list, nfeatures = nfeatures)
    
    ref.list <- lapply(X = ref.list, FUN = function(x) {
      #x <- ScaleData(x, features = features, verbose = TRUE)
      x <- RunPCA(x, features = features, verbose = FALSE)
      
    })
  }
 
  
  ref.anchors <- FindIntegrationAnchors(object.list = ref.list, 
                                        anchor.features = features, 
                                        reference = reference,
                                        reduction = "rpca",
                                        #k.anchor = 5,
                                        dims = ndims)
  
  # this command creates an 'integrated' data assay
  # this command creates an 'integrated' data assay
  if(correct.all){
    feature2integrate = rownames(seuratObj)
  }else{
    feature2integrate = features
  }
  
  if(is.null(merge.order)){
    ref.combined <- IntegrateData(anchorset = ref.anchors, 
                                  features.to.integrate = feature2integrate, 
                                  dims = ndims,
                                  k.weight = k.weight) 
  }else{
    ref.combined <- IntegrateData(anchorset = ref.anchors, 
                                  features.to.integrate = feature2integrate,
                                  sample.tree = merge.order,
                                  dims = ndims,
                                  k.weight = k.weight)
  }
  
  rm(ref.list)
  rm(ref.anchors)
  
  DefaultAssay(ref.combined) <- "integrated"
  
  ref.combined <- ScaleData(ref.combined, verbose = FALSE)
  ref.combined <- RunPCA(ref.combined, npcs = 50, verbose = FALSE)
  
  ElbowPlot(ref.combined, ndims = 50)
  
  ref.combined <- RunUMAP(ref.combined, reduction = "pca", dims = 1:30, n.neighbors = 30, 
                          min.dist = 0.3) 
  DimPlot(ref.combined, reduction = "umap", group.by = group.by, label = TRUE,
          repel = TRUE, raster=FALSE)
  
  #ref.combined <- FindNeighbors(ref.combined, reduction = "pca", dims = 1:20)
  #ref.combined <- FindClusters(ref.combined, resolution = 0.2)
  #DimPlot(ref.combined, reduction = "umap")
  
  return(ref.combined)
  
  
}

##########################################
# customized RunFastMNN 
# original code from https://github.com/satijalab/seurat-wrappers/issues/83
# fix the issue of not returning corrected gene expression for all genes
##########################################
RunFastMNN_fixed <- function (object.list, 
                              assay = NULL, 
                              features = 2000, 
                              reduction.name = "mnn", 
                              reduction.key = "mnn_", 
                              reconstructed.assay = "mnn.reconstructed",
                              verbose = TRUE, ...) 
{
  
  if (!all(sapply(X = object.list, FUN = inherits, what = "Seurat"))) {
    stop("'object.list' must be a list of Seurat objects", 
         call. = FALSE)
  }
  
  if (length(x = object.list) < 2) {
    stop("'object.list' must contain multiple Seurat objects for integration", 
         call. = FALSE)
  }
  
  assay <- assay %||% DefaultAssay(object = object.list[[1]])
  features_common = rownames(object.list[[1]])
  
  for (i in 1:length(x = object.list)) {
    DefaultAssay(object = object.list[[i]]) <- assay
    if(i >1) {features_common = intersect(features_common, rownames(object.list[[i]]))}
  }
  message(paste(length(features_common),  "common features for batch correction "))
  
  
  if (is.numeric(x = features)) {
    if (verbose) {
      message(paste("Computing", features, "integration features"))
    }
    
    ## select features for each batch
    features <- SelectIntegrationFeatures(object.list = object.list, 
                                          nfeatures = features, 
                                          assay = rep(assay, length(object.list)))
  }
  
  objects.sce <- lapply(X = object.list, 
                        FUN = function(x, f) {return(as.SingleCellExperiment(x = subset(x = x, features = f)))},
                        f = features_common)
  
  integrated <- merge(x = object.list[[1]], y = object.list[2:length(x = object.list)])
  out <- do.call(what = batchelor::fastMNN, args = c(objects.sce, 
                                                     list(subset.row = features, 
                                                          ...)
  ))
  
  rownames(x = SingleCellExperiment::reducedDim(x = out)) <- colnames(x = integrated)
  colnames(x = SingleCellExperiment::reducedDim(x = out)) <- paste0(reduction.key, 
                                                    1:ncol(x = SingleCellExperiment::reducedDim(x = out)))
  integrated[[reduction.name]] <- CreateDimReducObject(embeddings = SingleCellExperiment::reducedDim(x = out), 
                                                       assay = DefaultAssay(object = integrated), key = reduction.key)
  
  
  # Add reconstructed matrix (gene x cell)
  integrated[[reconstructed.assay]] <- CreateAssayObject(
    data = as(object = SummarizedExperiment::assay(x = out), Class = "sparseMatrix"),
  )
  # Add variable features
  VariableFeatures(object = integrated[[reconstructed.assay]]) <- features
  Tool(object = integrated) <- S4Vectors::metadata(x = out)
  integrated <- LogSeuratCommand(object = integrated)
  return(integrated)
  
}

IntegrateData_runFastMNN = function(seuratObj, group.by = 'dataset', 
                                    assays = 'RNA', 
                                    nfeatures = 2000,
                                    ndims = c(1:30),
                                    merge.order = NULL,
                                    cos.norm = TRUE,
                                    correct.all = TRUE,
                                    reference = NULL)
{
  # seuratObj = aa; group.by = "dataset"; nfeatures = 3000;  merge.order = NULL; correct.all = TRUE
  # seuratObj[['integrated']] = NULL
  #seuratObj = ScaleData(seuratObj, )
  
  library(Seurat)
  library(SeuratData)
  library(SeuratWrappers)
  
  #DefaultAssay(seuratObj) = 'RNA'
  ## one error solved (see https://github.com/satijalab/seurat-wrappers/issues/126)
  seuratObj = DietSeurat(object = seuratObj, counts = TRUE, data = TRUE, assays = assays, 
                         dimreducs = c('pca'))
  
  ## By default, batches are merged in the user-supplied order in ..., i.e., 
  ## the first batch is merged with the second batch, the third batch is merged with the combined 
  ## first-second batch, the fourth batch is merged with the combined first-second-third batch and so on. 
  ## We refer to this approach as a progressive merge. When batch is supplied for a single object in ..., 
  ## the ordering of batches in a progressive merge is determined by the ordering of factor levels in batch.
  cat(' used the SeuratWrappers RunFastMNN \n')
  refs.merged <- SeuratWrappers::RunFastMNN(SplitObject(object = seuratObj, split.by = group.by), 
                                            features = nfeatures,
                                            merge.order = merge.order,
                                            correct.all = correct.all, 
                                            cos.norm = TRUE)
  
  # if(correct.all){
  #   cat(' used the modified version of RunFastMNN \n')
  #   refs.merged <- RunFastMNN_fixed(SplitObject(object = seuratObj, split.by = group.by), 
  #                                   features = nfeatures,
  #                                   merge.order = merge.order,
  #                                   correct.all = correct.all, 
  #                                   cos.norm = FALSE)
  # }
  
  refs.merged <- RunUMAP(refs.merged, reduction = "mnn", dims = 1:30)
  
  p1 = DimPlot(refs.merged, group.by = group.by, label = TRUE, repel = TRUE, raster=FALSE)
  plot(p1)
  
  return(refs.merged)
  
  
}

IntegrateData_runHarmony = function(seuratObj, group.by = 'dataset', nfeatures = 2000, 
                                    dims.use = c(1:50), 
                                    nclust = NULL,
                                    redo.normalization.hvg.scale.pca = TRUE,
                                    reference_values = NULL, 
                                    max.iter.harmony = 20)
{
  library(harmony)
  library(Seurat)
  library(SeuratData)
  library(SeuratWrappers)
  # ref = srat;
  # refs.merged = merge(aa, y = ref, add.cell.ids = c("mNT", "mouseGastrulation"), project = "RA_competence")
  if(redo.normalization.hvg.scale.pca){
    refs.merged <- NormalizeData(seuratObj) %>% FindVariableFeatures(nfeatures = nfeatures) %>% RunPCA(verbose = FALSE)
  }else{
    refs.merged = FindVariableFeatures(seuratObj, nfeatures = nfeatures) %>% RunPCA(verbose = FALSE)
  }
  rm(seuratObj)
  
  refs.merged <- RunHarmony(refs.merged,
                            group.by.vars = group.by,
                            reduction = 'pca',
                            dims.use = dims.use, 
                            nclust = nclust,
                            reference_values = reference_values,
                            max.iter.harmony = max.iter.harmony, 
                            verbose = TRUE,
                            plot_convergence = TRUE
                            )
  
  refs.merged <- RunUMAP(refs.merged, reduction = "harmony", dims = 1:30)
  #refs.merged <- FindNeighbors(refs.merged, reduction = "harmony", dims = 1:30) %>% FindClusters()
  
  DimPlot(refs.merged, group.by = group.by, label = TRUE, repel = TRUE, raster=FALSE)
  
  ## a quick test
  Test_if_RunHarmony_works = FALSE
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
  
  return(refs.merged)
  
}


runScanorama <- function(object, batch,  do.umap = F,  
                         nfeatures = 2000, vars2regress = "percent.mt", dim = 50L, verbose = T)
{
  cat("Checking python modules...\n")
  
  library(reticulate, quietly = T)
  
  if (!(py_module_available("scanorama"))){
    cat("Installing 'Scanorama' python module...\n")
    py_install("scanorama", pip = T)
    
  }
  
  scanorama <- import('scanorama')
  
  if (!("SCT" %in% names(object@assays))){
    object <- getMitoContent(object)
    object <- SCTransform(object, vars.to.regress = vars2regress, assay = "RNA", method = "glmGamPoi", verbose = F, vst.flavor = "v2")
  }
  
  merge.exists <- F
  if ("Seurat" %in% class(object)){
    
    object <- getMitoContent(object)
    so <- object;
    
    merge.exists <- T
    
    if (!("SCT" %in% names(so@assays))){
      miko_message("Running SCTransform...", verbose = verbose)
      so <- SCTransform(so, vars.to.regress = vars2regress, assay = "RNA", method = "glmGamPoi", verbose = F, vst.flavor = "v2")
    }
    
    
    object <- Seurat::SplitObject(so, split.by = batch)
    
  } else if ("list" %in% class(object)){
    
    object <- pbapply::pblapply(X = object, FUN = function(x){
      x <- getMitoContent(x)
      return(x)
    })
    
    for (i in 1:length(object)){
      if (!("SCT" %in% names(object[[i]]@assays))){
        miko_message(paste0("Running SCTransform for '", names(object)[i], "'..."), verbose = verbose)
        object[[i]] <- SCTransform(object[[i]], vars.to.regress = vars2regress, assay = "RNA", method = "glmGamPoi", verbose = F, vst.flavor = "v2")
      }
    }
    
  }
  
  miko_message("Selecting integration features...", verbose = verbose)
  so.features <- SelectIntegrationFeatures(object.list = object,
                                           assay = rep("SCT", length(object)),
                                           nfeatures = nfeatures,
                                           fvf.nfeatures = nfeatures)
  
  
  miko_message("Computing residuals for integration features...", verbose = verbose)
  for (i in 1:length(object)){
    object[[i]] <- GetResidual(object = object[[i]],features =  so.features, assay = "SCT", replace.value = F, verbose = F)
  }
  
  
  if (!merge.exists){
    miko_message("Merging seurat objects...", verbose = verbose)
    so <- tryCatch({
      merge(object[[1]], y = object[-1]) # 2-3x faster than alternative method
    }, error = function(e){
      return(mergeSeuratList(object))     # slower but robust
    }, silent = T)
  }
  
  datasets <- list()
  genes_list <- list()
  for (i in 1:length(object)) {
    datasets[[i]] <- as.matrix( Matrix::t((
      object[[i]]@assays$SCT@scale.data[so.features , ]
    )))
    genes_list[[i]] <- colnames(datasets[[i]])
  }
  
  # Integration.
  miko_message("Running Scanorama Integration...", verbose = verbose)
  suppressMessages({
    suppressWarnings({
      integrated.data <- scanorama$integrate(datasets, genes_list, dimred= dim)
    })
  })
  
  # df.embed <- rbind(integrated.data[[1]])
  df.embed <- NULL
  for (i in 1:length(integrated.data[[1]])){
    
    current_embed <- (integrated.data[[1]][[i]])
    
    if ("list" %in% class(object)){
      rownames(current_embed) <- colnames(object[[i]])
    }
    # rownames(current_embed) <- colnames(so.list[[i]])
    colnames(current_embed) <- paste0("PC_", seq(1, ncol(current_embed)))
    
    df.embed <- rbind(df.embed, current_embed)
    
  }
  
  if ("Seurat" %in% class(object)){
    rownames(current_embed) <- colnames(object)
  }
  
  so@assays[["integrated"]] <- CreateAssayObject(counts = so@assays[["SCT"]]@counts, min.cells = 0, min.features = 0)
  so@assays[["integrated"]]@key <- "integrated_"
  
  
  # set default assay
  DefaultAssay(so) <- "integrated"
  
  so[["pca"]]  <- CreateDimReducObject(embeddings = df.embed,
                                       key = "PC_",
                                       global  = F,
                                       loading = new(Class = "matrix"),
                                       assay = "integrated")
  # cluster.UMAP(so)
  
  if (do.umap){
    so <- RunUMAP(object = so, dims = 1:dim, verbose = verbose)
  }
  
  return(so)
  
  
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
# calculate the similarity of query and reference (without integration) 
# inspired by https://github.com/quadbio/simspec
# and https://pubmed.ncbi.nlm.nih.gov/33711282/
##########################################
calculate_similarity_query_ref = function(query,
                                          ref,
                                          assay_use = 'RNA',
                                          find_hvg = TRUE,
                                          nHVGs = 500,
                                          features_use = NULL,
                                          method = c("spearman", "pearson"),  
                                          group.by = 'celltype')
{
  # query = subset(aa, cells = colnames(aa)[which(aa$condition == 'day3_noRA')]); method = 'pearson'
  library(simspec)
  library(ggplot2)
  library(data.table)
  
  if(find_hvg){
    #if(length(VariableFeatures(ref)) != nHVGs){
    #  ref = FindVariableFeatures(ref, selection.method = 'vst', nfeatures = nHVGs)
    #}
    if(length(VariableFeatures(query)) != nHVGs) 
    {
      query = FindVariableFeatures(query, selection.method = 'vst', nfeatures = nHVGs)
    }
    
  }
  
  if(is.null(features_use)){
    candidats = unique(c(VariableFeatures(query)))
  }else{
    candidats = features_use
  }
  
  candidats = candidats[which(!is.na(match(candidats, intersect(rownames(query), rownames(ref)))))]
  cat(length(candidats), ' HVGs selected for similarity calculation \n ')
  
  if(assay_use  == 'RNA'){
    ref_avg = AverageExpression(object = ref, features = candidats, group.by = group.by, 
                                slot = "data")$RNA
    query_mat = query@assays$RNA@data
  }else{
    ref_avg = eval(parse(text = 
    paste0("AverageExpression(object = ref, features = candidats, group.by = group.by, slot = 'data')$", 
           assay_use)))
    query_mat = eval(parse(text = paste0("query@assays$", assay_use, "@data")))
  }
  
  
  query_mat = query_mat[match(rownames(ref_avg), rownames(query_mat)), ]
  
  if(method == 'pearson'){
    corr <- qlcMatrix::corSparse(query_mat, ref_avg)
  }else if(method == 'spearman'){
    ranked_query <- presto::rank_matrix(query_mat)$X_ranked
    ranked_ref <- presto::rank_matrix(ref_avg)$X_ranked
    corr <- qlcMatrix::corSparse(ranked_query, t(ranked_ref))
  }
  
  colnames(corr) = colnames(ref_avg)
  
  xx = reshape2::melt(corr, value.name = "value")
  
  p <- ggplot(xx, aes(x=Var2, y=value)) + 
    geom_violin(trim=FALSE)
  p1 = p + stat_summary(fun.data=mean_sdl, 
                        #mult= 0.5, 
                   geom="pointrange", color="red") + 
    #theme_classic() +
    labs(x = 'cell types (ref)', y = 'correlation') + 
    #scale_y_continuous(breaks = c(1:nb_chromStates)) +
    theme(axis.text.x = element_text(angle = 90, size = 10), 
          axis.text.y = element_text(angle = 0, size = 10), 
          axis.title =  element_text(size = 12),
          legend.text = element_text(size=12),
          legend.title = element_text(size = 14)
    ) + 
    #guides(fill=guide_legend(title="emission")) +
    ggtitle(method)
  
  #plot(p1)
  return(p1)
  
}

##########################################
# test the neighborhood-based correlation from 
# https://github.com/MarioniLab/RabbitGastrulation2022/blob/master/8-compare_species/compare_nhoods.ipynb
##########################################
calculate_neighbourhood_based_correlation = function()
{
  # Load packages 
  suppressPackageStartupMessages(library(scrabbitr)) 
  suppressPackageStartupMessages(library(SingleCellExperiment)) 
  suppressPackageStartupMessages(library(miloR)) 
  suppressPackageStartupMessages(library(DelayedArray)) 
  suppressPackageStartupMessages(library(Matrix)) 
  suppressPackageStartupMessages(library(ggraph)) 
  suppressPackageStartupMessages(library(igraph)) 
  suppressPackageStartupMessages(library(viridis)) 
  suppressPackageStartupMessages(library(gridExtra)) 
  suppressPackageStartupMessages(library(RColorBrewer)) 
  suppressPackageStartupMessages(library(jsonlite)) 
  suppressPackageStartupMessages(library(ggrastr)) 
  suppressPackageStartupMessages(library(ggridges)) 
  suppressPackageStartupMessages(library(ggalluvial)) 
  suppressPackageStartupMessages(library(ggrepel)) 
  # TODO: Remove source("../data-in/rabbit/load_rabbit.R") source("../scrabbitr/R/plot_utils.R") #temp   
  
  # Load rabbit data 
  r_data <- readRDS(file = '../data/RabbitGastrulation2022/r_sce.rds')
  r_data 
  
  # Load mouse data
  m_data <- readRDS("../data/RabbitGastrulation2022/embryo_sce.rds")
  m_data
  
  # TODO: Remove
  m_pcs <- read.csv("../data/RabbitGastrulation2022/pca_batch_corrected.csv")
  rownames(m_pcs) = m_pcs$cell
  m_pcs = m_pcs[, -1]
  reducedDim(m_data,"PCA") <- as.matrix(m_pcs)
  
  #m_meta <- read.csv("../data/RabbitGastrulation2022/metadata_cells.csv")
  #m_datacelltype.clustering
  
  #reducedDim(m_data,"UMAP") <- m_meta2[,c("BBKNN_UMAP1","BBKNN_UMAP2")]
  m_genes <- read.csv("../data/RabbitGastrulation2022/metadata_genes.csv")
  
  # Load rabbit-mouse one-to-one orthologs
  rm_orthologs <- read.table("../data-in/orthologs/mmusculus.tsv",sep="\t")
  rm_orthologs[1:5,]
  
  
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


