##########################################################################
##########################################################################
# Project: Elad's heart regeneration
# Script purpose: functions of scRNA-seq data analysis
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Jul 20 16:03:21 2022
##########################################################################
##########################################################################
require(Seurat)
require(SeuratObject)
require(ggplot2)
require(tibble)
require(dplyr)
require(tictoc)
library(patchwork)

firstup <- function(x) {
  x = tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

get_geneName = function(xx)
{
  return(sapply(xx, function(x) {test = unlist(strsplit(as.character(x), '-')); 
  if(length(test) == 1) {
    test
  }else{
    test = test[-length(test)]; 
    paste0(test, collapse = '-')
  }
  }
  ))
}

get_geneID = function(xx)
{
  return(sapply(xx, function(x) {test = unlist(strsplit(as.character(x), '-')); return(test[length(test)])}))
}

########################################################
########################################################
# Section : import kalisto output to make 10x seurat object
# orignal code from Tomas' code
########################################################
########################################################
make_SeuratObj_scRNAseq = function(topdir = './', 
                                   saveDir = './results', 
                                   changeGeneName.axolotl = TRUE, 
                                   defaultDrops.only = FALSE,
                                   keyname = 'Amex_scRNA_d0',  
                                   QC.umi = FALSE)
{
  library(Seurat)
  library(DropletUtils)
  library(edgeR)
  library(BiocParallel)
  
  if(!dir.exists(saveDir)) dir.create(saveDir)
  
  # read in data
  # topdir = "${outbus}" # source dir
  exp = Matrix::readMM(paste0(topdir, "genecounts.mtx")) #read matrix
  bc = read.csv(paste0(topdir, "/genecounts.barcodes.txt"), header = F, stringsAsFactors = F)
  g = read.csv(paste0(topdir, "/genecounts.genes.txt"), header = F, stringsAsFactors = F)
  
  if(changeGeneName.axolotl){
    cat('change gene names \n')
    # change the gene names before making Seurat object
    annot = readRDS(paste0('/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                           'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))
    
    mm = match(g$V1, annot$geneID)
    ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
    g$V1[!is.na(mm)] = ggs[!is.na(mm)]
    
  }
  
  dimnames(exp) = list(paste0(bc$V1,"-1"), g$V1) # number added because of seurat format for barcodes
  count.data = Matrix::t(exp)
  
  #dimnames(exp) = list(paste0(bc\$V1,"-1"),g\$V1) # number added because of seurat format for barcodes
  #count.data = Matrix::t(exp)
  cat('-- filtering empty drops with UMI rank --  \n')
  
  # get emptyDrops and default cutoff cell estimates
  iscell_dd = defaultDrops(count.data, expected = 8000) # default cell estimate, similar to 10x cellranger
  cat(sum(iscell_dd, na.rm=TRUE), ' cell identified with default cellRanger method \n')
  
  if(!defaultDrops.only){
    eout = emptyDrops(count.data, lower = 100, BPPARAM = SerialParam())
    eout$FDR[is.na(eout$FDR)] = 1
    iscell_ed = eout$FDR<=0.01
    cat(sum(iscell_ed, na.rm=TRUE), ' cell identified with emptyDrops \n')
   
  }else{
    iscell_ed = rep(NA, length(iscell_dd))
  }
  
  meta = data.frame(row.names = paste0(bc$V1,"-1"),
                    iscell_dd = iscell_dd, iscell_ed = iscell_ed)
  
  # plot rankings for number of UMI
  br.out <- barcodeRanks(count.data)
  
  pdf(paste0(saveDir, "UMIrank.pdf"), height = 6, width =10, useDingbats = FALSE)
  
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  if(!defaultDrops.only) {abline(v = sum(iscell_ed), col = 'darkgreen', lwd = 2.0)}
  abline(v = sum(iscell_dd), col = 'darkblue', lwd = 2.0)
  abline(v = c(3000, 5000, 8000, 10000, 12000), col = 'gray')
  text(x = c(3000, 5000, 8000, 10000, 12000), y =10000, labels = c(3000, 5000, 8000, 10000, 12000), col = 'red')
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
         legend=c("knee", "inflection"))
  
  dev.off()
  
  if(QC.umi){
    # plot rankings for number of UMI
   
    # UMI duplication
    # umi = read.table("${umic}", sep = "\t", header = F, stringsAsFactors = F)
    # sumUMI = c()
    # sumi = sum(umi\$V4)
    
    cat('umi duplication check \n')
    umi = read.table(paste0(topdir, "umicount.txt"), sep = "\t", header = F, stringsAsFactors = F)
    colnames(umi) = c('cell.bc', 'umi', 'kallisto.seqIndex', 'counts')
    
    # for(i in 0:250){ sumUMI = c(sumUMI, sum(umi\$V4[umi\$V4>i])/sumi) }
    # pdf("${params.samplename}_UMIduplication.pdf", height = 3.5, width = 7, useDingbats = F)
    # par(mfrow = c(1,2))
    # plot(sumUMI, ylim = c(0,1), pch = 20, col = "grey30", ylab = "% of total reads",
    #      xlab = "More than xx UMI", main = "${params.samplename}")
    # diffUMI = sumUMI[-length(sumUMI)] - sumUMI[-1]
    # plot(diffUMI, ylim = c(0,0.2), pch = 20, col = "grey30", ylab = "Change in % of total reads",
    #      xlab = "More than xx UMI", main = "${params.samplename}")
    # dev.off()
    
    sumUMI = c()
    sumi = sum(umi$counts)
    
    for(i in 0:250){
      sumUMI = c(sumUMI, sum(umi$counts[umi$counts>i])/sumi)
    }
    
    counts.umi = c()
    for(i in 1:100){
      counts.umi = c(counts.umi, length(which(umi$counts == i))/nrow(umi))
    }
    
    pdf(paste0(saveDir, "UMIduplication.pdf"), height = 6, width = 12, useDingbats = F)
    
    par(mfrow = c(1,1))
    plot(1:100, counts.umi, ylim = c(0, 1), col = 'gray30', 
         ylab = '% of umi', xlab = 'duplication nb of umi', pch = 20)
    
    par(mfrow = c(1,2))
    
    plot(sumUMI, ylim = c(0,1), pch = 20, col = "grey30", ylab = "% of total reads",
         xlab = "More than xx UMI")
    
    diffUMI = sumUMI[-length(sumUMI)] - sumUMI[-1]
    plot(diffUMI, ylim = c(0,0.2), pch = 20, col = "grey30", ylab = "Change in % of total reads",
         xlab = "More than xx UMI")
    
    dev.off()
    
  }
  
  # create Seurat object
  ## we're only keeping what might potentially be a cell (by DD or ED)
  if(defaultDrops.only){
    srat = CreateSeuratObject(counts = count.data[, iscell_dd],
                              meta.data = meta[iscell_dd,], 
                              min.cells = 10, 
                              min.features = 50)
  }else{
    srat = CreateSeuratObject(counts = count.data[,iscell_dd | iscell_ed],
                              meta.data = meta[iscell_dd | iscell_ed,], 
                              min.cells = 10, 
                              min.features = 50)
    
  }
  
  amb_prop = estimateAmbience(count.data)[rownames(srat@assays$RNA@meta.features)]
  srat@assays$RNA@meta.features = data.frame(row.names = rownames(srat@assays$RNA@meta.features),
                                              "ambient_prop" = amb_prop)
  
  # get MT% (genes curated from NCBI chrMT genes)
  # mtgenes = c("COX1", "COX2", "COX3", "ATP6", "ND1", "ND5", "CYTB", "ND2", "ND4",
  #            "ATP8", "MT-CO1", "COI", "LOC9829747")
  # mtgenes = c(mtgenes, paste0("MT", mtgenes), paste0("MT-", mtgenes))
  # 
  # mt_features = rownames(srat)[get_geneName(g[,1]) %in% mtgenes]
  # srat = PercentageFeatureSet(srat, col.name = "percent.mt", assay = "RNA",
  #                            features = mt_features)
  
  saveRDS(srat, file = paste0(saveDir, "srat.RDS"))
  
  return(srat)
  
}


##########################################
# gene filtering using scran 
##########################################
geneFiltering.scran = function(srat = aa, gg.mt, gg.rb)
{
  # srat = aa;
  require(SingleCellExperiment)
  library(scran)
  library(scater)
  library(scuttle)
  library(Seurat)
  library(SeuratObject)
  
  sce <- as.SingleCellExperiment(srat)
  rm(srat)
  
  fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
  plotHighestExprs(sce, n=50, exprs_values = "counts")
  
  #fontsize <- theme(axis.text=element_text(size=16), axis.title=element_text(size=16))
  #plotHighestExprs(sce, n=30) + fontsize
  
  ave.counts <- calculateAverage(sce, assay.type = "counts")
  
  hist(log10(ave.counts), breaks=100, main="", col="grey80",
       xlab=expression(Log[10]~"average count"))
  
  num.cells <- nexprs(sce, byrow=TRUE)
  
  smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
                xlab=expression(Log[10]~"average count"))
  
  # detected in >= 5 cells, ave.counts >=5 but not too high
  genes.to.keep <- num.cells > 50 & ave.counts >= 10^-2  & ave.counts <10^2  
  summary(genes.to.keep)
  
  # remove mt and ribo genes
  genes.to.keep = genes.to.keep & ! rownames(sce) %in% gg.mt & ! rownames(sce) %in% gg.rb
  summary(genes.to.keep)
  
  sce <- sce[genes.to.keep, ]
  
  aa = Seurat::as.Seurat(sce)
  rm(sce);
  
  metadata = aa@meta.data
  library(Seurat)
  
  aa = CreateSeuratObject(counts = GetAssayData(object = aa, assay = "RNA", slot = "counts"),
                            meta.data = metadata[, c(1,3:10)])
  
  Idents(aa) = factor(aa$condition, levels = levels)
  VlnPlot(aa, features = 'nFeature_RNA', y.max = 10000)
  VlnPlot(aa, features = 'nCount_RNA', y.max = 100000)
  
  saveRDS(aa, file = paste0(RdataDir, 'seuratObject_merged_cellFiltered_doubletFinderOut.v2_geneFiltered.15kGene_', 
                             species, version.analysis, '.rds'))
  
  return(aa)
  
}

########################################################
########################################################
# Section : normalization methods
# 
########################################################
########################################################
Normalize_with_scran = function(aa, method = c('scran'))
{
  require(SingleCellExperiment)
  library(scran)
  library(scuttle)
  library(scRNA.seq.funcs)
  # sce =  as.SingleCellExperiment(subs)
  sce <- as.SingleCellExperiment(aa)
  
  if(method == 'scran'){
    clusters <- quickCluster(sce)
    sce <- computeSumFactors(sce, clusters=clusters)
    #summary(sizeFactors(scc))
    sce <- logNormCounts(sce)
    
  }
  
  if(method == 'sf'){
    expr_mat = counts(sce)
    sizefactors = calc_sf(expr_mat = expr_mat)
    sce <- logNormCounts(sce)
  }
  
  if(method == 'UQ'){
    
  }
  if(method == 'RLE'){
    
  }
  
  aa = Seurat::as.Seurat(sce)
  
  return(aa)
  
}



### test SC3
findClusters_SC3 = function(aa)
{
  #load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata')) 
  library(SC3)
  require(M3Drop)
  require(scran)
  require(SingleCellExperiment)
  
  DefaultAssay(aa) = 'Spatial'
  sce <- as.SingleCellExperiment(aa)
  
  expr_matrix =  exp(logcounts(sce))
  
  Brennecke_HVG <- BrenneckeGetVariableGenes(
    expr_mat = expr_matrix,
    spikes = NA,
    fdr = 0.2,
    minBiolDisp = 0.2
  )
  
  #HVG_genes <- Brennecke_HVG
  
  ## another method to identify by Kolodziejczyk AA, Kim JK, Tsang JCH et al. (2015)
  assay(sce, "normcounts") <- exp(logcounts(sce))
  means <- rowMeans(normcounts(sce))
  cv2 <- apply(normcounts(sce), 1, var)/means^2
  dm.stat <- DM(means, cv2)
  #head(dm.stat)
  
  DM_HVG = names(dm.stat)[which(dm.stat>0.3)]
  #DM_HVG = DM_HVG[which()]
  #dev.off()
  
  sce.HVG.Brenneck = sce[rownames(sce)%in%Brennecke_HVG, ]
  sce.HVG.DM = sce[rownames(sce)%in%DM_HVG, ]
  
  save(sce, sce.HVG.Brenneck, sce.HVG.DM, file=paste0(RdataDir, version.DATA, 
                                                      '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))
  #HVG_genes <- Brennecke_HVG$Gene
  #load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))
  
  sce.sels = sce.HVG.Brenneck
  #sce.sels = sce.HVG.DM
  sce = sce.sels
  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)
  # remove features with duplicated names
  sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
  
  plotPCA(
    sce,
    run_args = list(exprs_values = "logcounts"),
    colour_by = "total_counts",
    size_by = "total_features_by_counts"
  )
  
  ### test SC3 clustering method
  sce = sc3_estimate_k(sce)
  metadata(sce)$sc3$k_estimation
  
  sce <- sc3(sce.sels, gene_filter = FALSE, ks = 2:30, biology = TRUE, n_cores = 6)
  #rowData(sce)$feature_symbol
  #sce <- sc3(sce, ks = 2, gene_filter = TRUE, biology = TRUE)
  #deng <- plotTSNE(deng, rand_seed = 1, return_SCE = TRUE)
  
  col_data <- colData(sce)
  head(col_data[ , grep("sc3_", colnames(col_data))])
  plotPCA(
    sce, 
    colour_by = "sc3_10_clusters", 
    size_by = "sc3_10_log2_outlier_score"
  )
  
  plotPCA(
    sce, 
    colour_by = "sc3_3_clusters",
    shape_by = "celltypes",
    size_by = "total_features_by_counts"
  )
  
  row_data <- rowData(sce)
  head(row_data[ , grep("sc3_", colnames(row_data))])
  
  plotFeatureData(
    sce, 
    aes(
      x = sc3_3_markers_clusts, 
      y = sc3_3_markers_auroc, 
      colour = sc3_3_markers_padj
    )
  )
  
  set.seed(1)
  plotTSNE(sce, colour_by="sc3_6_clusters", size_by = "total_features_by_counts",
           run_args = list(perplexity=20)) 
  
  #sc3_plot_consensus(sce, k = 3)
  sc3_plot_consensus(
    sce, k = 3, 
    show_pdata = c(
      "celltypes", 
      "sc3_3_clusters", 
      "log10_total_features_by_counts",
      "log10_total_counts",
      
      "sc3_3_log2_outlier_score"
    )
  )
  
  sc3_plot_cluster_stability(sce, k = 6)
  
  sc3_plot_de_genes(sce, k = 3, p.val = 0.2)
  
  sc3_plot_markers(sce, k = 3)
  
}


##########################################
# doublet detection with DoubletFinder for scATAC-seq
# the old version of function is from 
# https://github.com/lengfei5/scRNAseq_MS_lineage_dev/blob/development/scripts/scATAC_functions.R
##########################################
detect_doubletCell_scRNAseq = function(seuratObj)
{
  # seuratObj = aa
  cc = unique(seuratObj$condition)
  seuratObj$DF_out = NA
  print(cc)
  
  for(n in 1:length(cc))
  {
    # n = 1
    cat(' ------------------------\n')
    cat(n, ' --- ', cc[n], '---\n')
    cat(' ------------------------\n')
    
    subs <- subset(seuratObj, condition == cc[n])
    
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
    
    seuratObj$DF_out[match(colnames(subs), colnames(seuratObj))] = subs$DF_out
    
  }
  
  return(seuratObj)
  
}


detect_doubletCell_scATAC= function(seurat.cistopic)
{
  ##########################################
  # check the distribution of fragment for cluster 7
  ##########################################
  #source.my.script('scATAC_functions.R')
  filtered_mtx_dir = paste0("../output_cellranger.ce11_scATACpro/filtered_matrix_peaks_barcodes")
  tenx.bmat = load_tenx_atac(paste0(filtered_mtx_dir, '/matrix.mtx'), 
                             paste0(filtered_mtx_dir, '/peaks.bed'), 
                             paste0(filtered_mtx_dir, '/barcodes.tsv'))
  
  pdfname = paste0(resDir, "/detect_doublet/distribution_counts_each_clusters.pdf")
  pdf(pdfname, width=12, height = 8)
  par(cex =1.0, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  par(mfrow=c(4,4))
  for(n in c(1:16)){
    # n = 1
    cat('cluster -- ', n, '\n')
    mm = match(colnames(seurat.cistopic)[which(seurat.cistopic$peaks_snn_res.0.8 == n)], 
               colnames(tenx.bmat))
    bmat.cluster = as.numeric(tenx.bmat[, mm])
    xx = bmat.cluster[which(bmat.cluster<=5)]
    hist(xx, main = paste0('histogrma of cluster ', n))
    #abline(v = log2(1+2)+1, col = 'red')
  }
  
  dev.off()
  
  ##########################################
  # detect doublets using published methods 
  ##########################################
  library(DoubletFinder)
  DefaultAssay(seurat.cistopic) = 'RNA'
  
  seu_kidney <- ScaleData(seurat.cistopic)
  seu_kidney <- FindVariableFeatures(seu_kidney, selection.method = "vst", nfeatures = 2000)
  seu_kidney <- RunPCA(seu_kidney, verbose = TRUE)
  seu_kidney <- RunUMAP(seu_kidney, dims = 1:10)
  
  sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- seu_kidney@meta.data$peaks_snn_res.0.8
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*length(seu_kidney$peaks_snn_res.0.8))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, 
                                 sct = FALSE)
  
  seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, 
                                 nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_816", sct = FALSE)
  
  xx = data.frame(seu_kidney$peaks_snn_res.0.8, seu_kidney$DF.classifications_0.25_0.09_770)
  colnames(xx) = c('cluster', 'doublet')
  yy = table(xx$cluster)
  yy = rbind(yy, table(xx$cluster[which(xx$doublet=='Singlet')]))
  yy = rbind(yy, table(xx$cluster[which(xx$doublet=='Doublet')]))
  rownames(yy) = c('total', 'singlet', 'doublet')
  
  barplot(yy[c(1,3),], main="nb of total and doublet",
          xlab="cluster index", col=c("darkblue","red"),
          legend = rownames(yy)[c(1,3)], beside=TRUE)
  
  xx = seurat.cistopic
  xx$Singlet = seu_kidney$DF.classifications_0.25_0.09_770
  #xx = subset(seurat.cistopic, cells = which(seu_kidney$DF.classifications_0.25_0.09_770 == 'Singlet') )
  
  nb.pcs = 80;
  n.neighbors = 50; min.dist = 0.1;
  
  xx <- RunUMAP(object = xx, reduction = 'pca', 
                dims = 1:nb.pcs, 
                n.neighbors = n.neighbors, min.dist = min.dist)
  
  DimPlot(xx, label = TRUE, pt.size = 1.5, label.size = 5,  repel = TRUE, split.by = 'Singlet') + 
    NoLegend()
  
  
}

##########################################
# explore the umap parameter combination 
##########################################
explore.umap.params.combination = function(sub.obj,
                                           resDir = './',
                                          pdfname = 'UMAP_explore_parameterCombinations.pdf',
                                          group.by = 'seurat_clusters', 
                                          cols = NULL,
                                          with_legend = TRUE,
                                          weight.by.var = TRUE,
                                          nfeatures.sampling = c(5000, 8000, 10000),
                                          nb.pcs.sampling = c(20, 30, 50), 
                                          n.neighbors.sampling = c(30, 50),
                                          min.dist.sampling = c(0.1, 0.3), 
                                          spread.sampling = 1,
                                          use.parallelization = FALSE
)
{
  options(future.globals.maxSize = 80000 * 1024^2)
  # 
  if(use.parallelization){
    cat('use future package to paralle \n')
    library(future)
    plan()
    # change the current plan to access parallelization
    plan("multicore", workers = 16)
    plan()
    
  }
  
  if (length(dev.list()!=0)) {dev.off()}
  
  pdfname = paste0(resDir,'/', pdfname); while (!is.null(dev.list()))  dev.off();
  
  pdf(pdfname, width=12, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  for(nfeatures in nfeatures.sampling)
  {
    # nfeatures = nfeatures.sampling[1]
    cat('------------- nfeatures - ', nfeatures, '\n')
    sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures, 
                                    verbose = FALSE)
    # not scale every time
    # sub.obj = ScaleData(sub.obj, verbose = FALSE)
    sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, 
                      npcs = max(50, nb.pcs.sampling), 
                      weight.by.var = weight.by.var)
    
    for(nb.pcs in nb.pcs.sampling)
    {
      for(n.neighbors in n.neighbors.sampling)
      {
        for(min.dist in min.dist.sampling)
        {
          for(spread in spread.sampling){
            cat('--- nb.pcs - ', nb.pcs, ', n.neighbors - ', n.neighbors, 
                ', min.dist - ', min.dist, ', spread - ', spread,  '\n')
            # nfeatures = 5000;
            # nb.pcs = 50 # nb of pcs depends on the considered clusters or ids 
            # n.neighbors = 50;
            # min.dist = 0.05; spread = 1;
            sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", 
                               dims = 1:nb.pcs, 
                               spread = spread, n.neighbors = n.neighbors, 
                               min.dist = min.dist, verbose = FALSE)
            
            if(with_legend){
              if(is.null(cols)){
                pp = DimPlot(sub.obj, group.by = group.by, reduction = 'umap', label = TRUE, label.size = 6, 
                             repel = TRUE) + 
                  ggtitle(paste0('nfeatures - ', nfeatures,  ', nb.pcs - ', nb.pcs, ', n.neighbors - ', n.neighbors, 
                                 ', min.dist - ', min.dist, ', spread - ', spread))
              }else{
                pp = DimPlot(sub.obj, group.by = group.by, reduction = 'umap', label = TRUE, label.size = 6, 
                             repel = TRUE, cols = cols) + 
                  ggtitle(paste0('nfeatures - ', nfeatures,  ', nb.pcs - ', nb.pcs, ', n.neighbors - ', n.neighbors, 
                                 ', min.dist - ', min.dist, ', spread - ', spread))
              }
             
            }else{
              if(is.null(cols)){
                pp = DimPlot(sub.obj, group.by = group.by, reduction = 'umap', label = TRUE, label.size = 6, 
                             repel = TRUE) + 
                  NoLegend() + 
                  ggtitle(paste0('nfeatures - ', nfeatures,  ', nb.pcs - ', nb.pcs, ', n.neighbors - ', n.neighbors, 
                                 ', min.dist - ', min.dist, ', spread - ', spread))
              }else{
                pp = DimPlot(sub.obj, group.by = group.by, reduction = 'umap', label = TRUE, label.size = 6, 
                            repel = TRUE, cols = cols) + 
                  NoLegend() + 
                  ggtitle(paste0('nfeatures - ', nfeatures,  ', nb.pcs - ', nb.pcs, ', n.neighbors - ', n.neighbors, 
                                 ', min.dist - ', min.dist, ', spread - ', spread))
              }
             
            }
            
            plot(pp)
          }
        }
      }
      
    }
    
  }
  
  dev.off()
  
  ## change back to single core
  if(use.parallelization){
    library(future)
    plan()
    # change the current plan to access parallelization
    plan("multicore", workers = 1)
    plan()
    
  }
  
}
