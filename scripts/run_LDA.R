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

options(future.globals.maxSize = 80000 * 1024^2)
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

# srat_cr = subset(srat_cr,  downsample = 500)

# Run LDA on count matrix -------------------------------------------------
cat("Running LDA \n")
count.mat = srat_cr@assays$ATAC@counts
#count.mat = count.mat[c(1:5000), ] # subset peaks 

# remove empty cols
print("Removing empty cells...")
print(dim(count.mat))
cols.empty <- colSums(count.mat) == 0
count.mat <- count.mat[, !cols.empty]
print(dim(count.mat))

# binarize matrix
count.mat.orig <- count.mat
count.mat <- BinarizeMatrix(count.mat)
print(paste('Max count after binarizing', max(count.mat)))

topic.vec = c(seq(6, 18, by = 2), seq(20,  100, by = 10))

tic("LDA running time")
if (length(topic.vec) > 1){
  print("Running multicore LDA for topics:")
  print(topic.vec)
  
  out.lda <- parallel::mclapply(topic.vec, function(nc)
    {topicmodels::LDA(x = t(count.mat), k = nc, method = "Gibbs", control=list(seed=0))}, 
    mc.cores = length(topic.vec)
  )

}else{
  print("Running single LDA for topics:")
  print(topic.vec)
  
  out.lda <- topicmodels::LDA(x = t(count.mat), k = topic.vec, method = "Gibbs", control=list(seed=0))
}

toc()

saveRDS(out.lda, file = paste0(RdataDir, 'test_LDA_saved_v1.rds'))

Process_LDA_results = FALSE
if(Process_LDA_results){# perplexity(out.lda)
  
  out.lda = readRDS(file = paste0(RdataDir, 'test_LDA_saved_v1.rds'))
  
  tm.result <- posterior(out.lda[[1]])
  tm.result <- AddTopicToTmResult(tm.result)
  topics.mat = tm.result$topics
  # save output
  #print("Saving LDA")
  #save(out.lda, count.mat, count.mat.orig, file = outpath)
  #print("Time elapsed after LDA")
  #print(Sys.time() - jstart)
  
  plot_umap_seurat = TRUE
  
  if(plot_umap_seurat){
    # method='Z-score'
    # #method = 'Probability'
    # cistopicObject.reduced_space = t(cisTopic::modelMatSelection(cistopicObject,
    #                                                              target='cell',
    #                                                              method=method))
    # 
    # colnames(cistopicObject.reduced_space) = paste0('PC_', 1:ncol(cistopicObject.reduced_space))
    # 
    # dimensions = ncol(cistopicObject.reduced_space)
    # 
    # #cistopicObject.seurat = run_dim_reduction(,
    # #                                          cistopicObject.reduced_space,
    # #                                          dims=1:120,
    # #                                          reduction='pca')
    # 
    # seurat_obj = Seurat::CreateSeuratObject(cistopicObject@binary.count.matrix, assay = 'peaks')
    # seurat_obj[['pca']] = Seurat::CreateDimReducObject(embeddings=cistopicObject.reduced_space, 
    #                                                    key='PC_', 
    #                                                    assay='peaks')
    # 
    # seurat_obj$subtypes = srat_cr$subtypes[match(colnames(seurat_obj), colnames(srat_cr))]
    # seurat_obj$celltypes = srat_cr$celltypes[match(colnames(seurat_obj), colnames(srat_cr))]
    # 
    # # seurat_obj = Seurat::L2Dim(seurat_obj, reduction='pca') 
    # seurat_obj = Seurat::RunUMAP(seurat_obj, reduction = "pca", dims = 1:120, n.neighbors = 30, 
    #                              min.dist = 0.05)
    # 
    # DimPlot(seurat_obj, reduction = 'umap', group.by = 'celltypes', label = TRUE, repel = TRUE) + 
    #   NoLegend()
    # 
    # ggsave(paste0(resDir, '/UMAP_DM_cisTopic_test.pdf'), width = 8, height = 6)
    # 
    # DimPlot(seurat_obj, reduction = 'umap', group.by = 'subtypes', label = TRUE, repel = TRUE) + NoLegend()
    
    #cistopicObject.seurat = cistopicObject.seurat %>%
    #  Seurat::FindNeighbors(reduction=reduction, nn.eps=0.25, dims=1:dimensions) %>%
    #  Seurat::FindClusters(reduction=reduction, n.start=20, resolution=resolution)
    
    #srat_cr[['Topic']] = Seurat::CreateDimReducObject(embeddings=cistopicObject.reduced_space, 
    #                                                key='Topic_', 
    #                                                   assay='ATAC')
    
  }else{
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
