rm(list = ls())

version.analysis = '_R13591_atac_reseq_20221115'

resDir = paste0("../results/sc_multiome", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')


if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R14353_ax_snATAC_reseq'

#source('functions_scATAC.R')
#source('functions_scRNAseq.R')
#source('functions_Visium.R')

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(patchwork)
require(SeuratObject)

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
require(tictoc)

suppressPackageStartupMessages(library("argparse"))
library(topicmodels)
library(dplyr)
library(ggplot2)
library(Matrix)
library(here)
library(parallel)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

options(future.globals.maxSize = 80 * 1024^3)
set.seed(1234)
mem_used()

library(data.table)
require(cisTopic)

##########################################
# start the main function 
##########################################
srat_cr = readRDS(file = paste0(RdataDir, 
                                'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                'scATAC.merged.peaks.cr.',
                                '584K.annot_38280cells.rds'))
# normalize ATAC and UMAP
DefaultAssay(srat_cr) <- "ATAC"

Filter.cells.with.scATAC = FALSE
if(Filter.cells.with.scATAC){
  # quick filtering 
  srat_cr <- subset(
    x = srat_cr,
    subset = nCount_ATAC < 100000 &
      nCount_RNA < 25000 &
      nCount_ATAC > 200 &
      nCount_RNA > 1000 &
      nucleosome_signal < 6 &
      TSS.enrichment > 1
  )
}

Idents(srat_cr) = as.factor(srat_cr$condition)
#srat_cr = subset(srat_cr,  downsample = 1000)

srat_cr <- RunTFIDF(srat_cr)
srat_cr = FindTopFeatures(srat_cr, min.cutoff = 'q25')
srat_cr = subset(srat_cr, features = VariableFeatures(srat_cr, assay = 'ATAC'))

# Run LDA on count matrix -------------------------------------------------
cat("Running LDA \n")
count.mat = srat_cr@assays$ATAC@counts
#count.mat = count.mat[c(1:10000), ] # subset peaks 

# remove empty cols
print("Removing empty cells...")
print(dim(count.mat))

term_freq = rowSums(count.mat>0)/ncol(count.mat)
#peak.mean = apply()
jj = which(term_freq > 5/ncol(count.mat) & term_freq<10^-1) 
count.mat <- count.mat[jj, ] # filtered the terms with low (>5 cells) and high freq (<0.1)
print(dim(count.mat))

cols.empty <- colSums(count.mat) == 0
#rows.empty = rowSums(count.mat) == 0
count.mat <- count.mat[, !cols.empty]
print(dim(count.mat))

# binarize matrix
count.mat.orig <- count.mat
count.mat <- BinarizeMatrix(count.mat)
print(paste('Max count after binarizing', max(count.mat)))

topic.vec = unique(c(seq(50, 300, by=50), seq(20,  300, by = 20)))

tic("LDA running time")
if (length(topic.vec) > 1){
  print("Running multicore LDA for topics:")
  print(topic.vec)
  
  out.lda <- parallel::mclapply(topic.vec, function(nc)
  {topicmodels::LDA(x = t(count.mat), k = nc, method = "Gibbs", control=list(seed=0))}, 
  mc.cores = length(topic.vec)
  )
  
  # out.lda <- parallel::mclapply(topic.vec, function(nc)
  # {topicmodels::LDA(x = t(count.mat), k = nc, method = "Gibbs", control=list(seed=0))}, 
  # mc.cores = length(topic.vec)
  # )
  # 
  # out.lda <- parallel::mclapply(topic.vec, function(nc)
  # {topicmodels::LDA(x = t(count.mat), k = nc, method = "Gibbs", control=list(seed=0))}, 
  # mc.cores = length(topic.vec)
  # )
  
  #out.lda <- parallel::mclapply(topic.vec, function(nc){
  #  LDA(x = t(count.mat), k = nc, method = "Gibbs", control=list(seed=0))            
  #}, mc.cores = length(topic.vec))
  
  
}else{
  print("Running single LDA for topics:")
  print(topic.vec)
  # topic.vec = 20
  out.lda <- topicmodels::LDA(x = t(count.mat), 
                              k = topic.vec, 
                              method = "Gibbs", 
                              control=list(seed=0))
  
}

toc()

saveRDS(out.lda, file = paste0(RdataDir, 'test_LDA_saved_topFeatures.q25_filtered.termFreq_v7.rds'))

Process_LDA_results = FALSE
if(Process_LDA_results){ # perplexity(out.lda)
  
  #out.lda = readRDS(file = paste0(RdataDir, 'test_LDA_saved_v1.rds'))
  #out.lda = readRDS(file = paste0(RdataDir, 'test_LDA_saved_peaks.all_v4.rds'))
  #out.lda = readRDS(file = paste0(RdataDir, 'test_LDA_saved_topFeatures.q25_v5.rds'))
  #out.lda = readRDS(file = paste0(RdataDir, 'test_LDA_saved_topFeatures.q10_v5.rds'))
  out.lda = readRDS(file = paste0(RdataDir, 'test_LDA_saved_topFeatures.q25_filtered.termFreq_v7.rds'))
  
  #print("Saving LDA")
  #save(out.lda, count.mat, count.mat.orig, file = outpath)
  #print("Time elapsed after LDA")
  #print(Sys.time() - jstart)
  
  plot_umap_seurat = TRUE
  
  if(plot_umap_seurat){
    
    reduction='pca'
    
    #method ='Z-score'
    method = 'Probability'
    
    for(n in 1:length(out.lda))
    {
      # n = 4
      xx = out.lda[[n]]
      nb_topics = xx@k
      cat(n,  '-- nb of topics : ', nb_topics, '\n')
      
      tm.result <- posterior(out.lda[[n]])
      tm.result <- AddTopicToTmResult(tm.result)
      topics.mat = tm.result$topics
      rm(xx)
      
      # cistopicObject.reduced_space = t(cisTopic::modelMatSelection(cistopicObject,
      #                                                              target='cell',
      #                                                              method=method))
      colnames(topics.mat) = paste0('PC_', 1:ncol(topics.mat))
      dimensions = ncol(topics.mat)
      
      seurat_obj = subset(srat_cr, cells = rownames(topics.mat))
      
      topics.mat = topics.mat[match(colnames(seurat_obj), rownames(topics.mat)), ]
      
      #topics.mat <- scale(topics.mat, center = TRUE, scale = TRUE)
      
      if(method == "Z-score"){
        #topics.mat_zscore = log10(topics.mat)
        topics.mat_zscore <- apply(topics.mat, 2, scale, center=TRUE, scale=TRUE)
        rownames(topics.mat_zscore) = rownames(topics.mat)
        seurat_obj[['pca']] = Seurat::CreateDimReducObject(embeddings=topics.mat_zscore,
                                                           key='PC_',
                                                           assay='ATAC')
        
      }else{
        seurat_obj[['pca']] = Seurat::CreateDimReducObject(embeddings=topics.mat,
                                                           key='PC_',
                                                           assay='ATAC')
      }
      
      reduction = 'pca.l2'
      seurat_obj = seurat_obj %>%
        Seurat::L2Dim(reduction='pca') %>%
        Seurat::RunUMAP(#metric = "euclidean",
                        #metric = "cosine",
                        reduction = reduction, 
                        dims = 1:dimensions, 
                        n.neighbors = 20,
                        min.dist = 0.05)
      
      DimPlot(seurat_obj, reduction = 'umap', group.by = 'celltypes', label = TRUE, repel = TRUE) +
        #NoLegend() + 
        ggtitle(paste0('celltypes - nb.topics : ', nb_topics))
      
      p1 = DimPlot(seurat_obj, reduction = 'umap', group.by = 'celltypes', label = TRUE, repel = TRUE) +
        NoLegend() + ggtitle(paste0('celltypes - nb.topics : ', nb_topics))
      p2 = DimPlot(seurat_obj, reduction = 'umap', group.by = 'subtypes', label = TRUE, repel = TRUE) + 
        NoLegend() + ggtitle(paste0('subtypes - nb.topics : ', nb_topics))
      p1 + p2
      
      ggsave(paste0(resDir, '/UMAP_LDA_peaks.topFeatures.q25_filtered.termFreq_v7_methods.', method, '_', 
                    reduction, 
                    '_nb.topics.', nb_topics, '.pdf'), 
             width = 16, height = 6)
      
    }
    
    saveRDS(seurat_obj, file = paste0(RdataDir,
                                      'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                      'scATAC.merged.peaks.cr_',
                                      '584K.features_38247cells_200topics.rds'))
    
    #cistopicObject.seurat = cistopicObject.seurat %>%
    #  Seurat::FindNeighbors(reduction=reduction, nn.eps=0.25, dims=1:dimensions) %>%
    #  Seurat::FindClusters(reduction=reduction, n.start=20, resolution=resolution)

    #srat_cr[['Topic']] = Seurat::CreateDimReducObject(embeddings=cistopicObject.reduced_space,
    #                                                key='Topic_',
    #                                                   assay='ATAC')
    
    ## compare the umap from Seurat
    srat_cr <- RunTFIDF(srat_cr)
    srat_cr = FindTopFeatures(srat_cr, min.cutoff = 'q5')
    srat_cr <- RunSVD(srat_cr)
    
    DepthCor(srat_cr, n = 30)
    
    #cordat = DepthCor(srat_cr, reduction = "lsi", n = 30)$data
    #dims_use = cordat$Component[abs(cordat$counts)<0.3]
    
    dims_use = c(2:30)
    print(dims_use)
    
    srat_cr <- RunUMAP(object = srat_cr, reduction = 'lsi', dims = 2:30, n.neighbors = 30, min.dist = 0.1, 
                       reduction.name = "umap_lsi")
    
    DimPlot(object = srat_cr, label = TRUE, reduction = 'umap_lsi', group.by = 'celltypes') + NoLegend()
    ggsave(paste0(resDir, '/UMAP_seurat_tfidf.pdf'), width = 8, height = 6)
    
    
  }else{
    #philentropy::getDistMethods()
    #dists = philentropy::distance(topics.mat, method = "jensen-shannon")
    
    dists = sccore::jsDist(t(topics.mat))
    
    
    jsettings <- umap.defaults
    jsettings$n_neighbors <- 30
    jsettings$min_dist <- 0.1
    jsettings$random_state <- 123
    #umap.out <- umap(topics.mat, config = jsettings)
    #dat.umap.long <- data.frame(cell = rownames(umap.out$layout), 
    #                           umap1 = umap.out$layout[, 1], 
    #                            umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)
    
    dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)
    
    cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", 
                   "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
    
    ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) +
      geom_point() +
      theme_bw() +
      ggtitle(paste("LDA test")) +
      scale_color_manual(values = cbPalette) +
      theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.position = "right")
    
  }
  
}
