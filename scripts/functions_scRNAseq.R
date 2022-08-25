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
  return(sapply(xx, function(x) {test = unlist(strsplit(as.character(x), '-')); test = test[-length(test)]; 
  paste0(test, collapse = '-')}))
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
# doublet detection with DoubletFinder
# the old version of function is from 
# https://github.com/lengfei5/scRNAseq_MS_lineage_dev/blob/development/scripts/scATAC_functions.R
##########################################
detect_doubletCell= function(seurat.cistopic)
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
  options(future.globals.maxSize = 20000 * 1024^2)
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
    sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    sub.obj = ScaleData(sub.obj, verbose = FALSE)
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
            sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, 
                               spread = spread, n.neighbors = n.neighbors, 
                               min.dist = min.dist, verbose = FALSE)
            
            if(with_legend){
              pp = DimPlot(sub.obj, group.by = group.by, reduction = 'umap', label = TRUE, label.size = 6, 
                           pt.size = 2, repel = TRUE) + 
                ggtitle(paste0('nfeatures - ', nfeatures,  ', nb.pcs - ', nb.pcs, ', n.neighbors - ', n.neighbors, 
                               ', min.dist - ', min.dist, ', spread - ', spread))
            }else{
              pp = DimPlot(sub.obj, group.by = group.by, reduction = 'umap', label = TRUE, label.size = 6, 
                           pt.size = 2, repel = TRUE) + 
                NoLegend() + 
                ggtitle(paste0('nfeatures - ', nfeatures,  ', nb.pcs - ', nb.pcs, ', n.neighbors - ', n.neighbors, 
                               ', min.dist - ', min.dist, ', spread - ', spread))
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


########################################################
########################################################
# Section : subclustering and annotation
# 
########################################################
########################################################
subclustering_manual.annotation = function(aa)
{
  ##########################################
  # transfer the cell labels from spliced data
  ##########################################
  Transferring_cellLabels_splicedData = FALSE
  if(Transferring_cellLabels_splicedData){
    RdataDir_spliced = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/results/Rdata/'
    bb = readRDS(file = paste0(RdataDir_spliced, 'seuratObject_axloltl_scRNAseq_splicedOnly_manualAnnotation.rds'))
    mm = match(colnames(aa), colnames(bb))
    aa$celltypes = bb$celltypes[mm]
    aa$subtypes = bb$subtypes[mm]
    
    rm(bb)
    
    aa$clusters = aa$seurat_clusters
    
    saveRDS(aa, file = paste0(RdataDir, 'seuratObject_', species, version.analysis, '_lognormamlized_pca_umap_manualAnnot.rds'))
    
    ##########################################
    # some cleaning for each cluster ??
    ##########################################
    
  }
  
  DimPlot(aa, label = TRUE, repel = TRUE) + ggtitle("snRNAseq (multiome)")
  p1 = FeaturePlot(aa, features = 'nFeature_RNA')
  p2 = FeaturePlot(aa, features = 'percent.mt')
  
  p1 | p2
  
  ggsave(filename = paste0(resDir, '/FeaturesPlot_nFeatures_pertMT.pdf'), width = 14, height = 6)
  
  # annotation for all cells
  Idents(aa) = aa$subtypes
  
  DimPlot(aa, label = TRUE, group.by = 'subtypes', repel = TRUE) + ggtitle("first manual annotation")
  ggsave(filename = paste0(resDir, '/first_test_umap_v2_manualAnnot.pdf'), width = 10, height = 8)
  
  aa$condition = factor(aa$condition, levels = paste0('Amex_scRNA_d', c(0, 1, 4, 7, 14)))
  aa$subtypes[which(is.na(aa$subtypes))] = 'unknow'
  
  DimPlot(aa, label = TRUE, repel = TRUE, split.by = 'condition') + NoLegend() +  
    ggtitle("first manual annotation")
  
  ggsave(filename = paste0(resDir, '/first_test_umap_v2_manualAnnot_byCondition.pdf'), width = 25, height = 8)
  
  # aa$subtypes = aa$celltypes
  # aa$clusters = aa$seurat_clusters
  # markers.coarse = readRDS(file = paste0(RdataDir_spliced, 'top10_markerGenes_coarseCluster.rds'))
  
  ##########################################
  # ## subclustering CMs
  ##########################################
  celltype.sels = 'CM'
  sub.obj = subset(aa, cells = colnames(aa)[!is.na(match(aa$celltypes, celltype.sels))])
  
  #sub.obj = SCTransform(aa, ncells = 3000, assay = "RNA", verbose = FALSE, 
  #            variable.features.n = 3000, return.only.var.genes = TRUE, vst.flavor = "v2")
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
  sub.obj = ScaleData(sub.obj)
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 30)
  
  nb.pcs = 20 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 20; min.dist = 0.05;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, 
                     n.neighbors = n.neighbors,
                     min.dist = min.dist)
  
  # sub.obj$clusters = sub.obj$seurat_clusters
  p1 = DimPlot(sub.obj, group.by = 'clusters', reduction = 'umap', label = TRUE, label.size = 5) +
    ggtitle(celltype.sels)
  
  features = rownames(aa)[grep('MYH6|ACTN2|NPPA|TNN|GATA4', rownames(aa))]
  FeaturePlot(aa, features = features, cols = c('gray', 'red'))
  
  sub.obj <- FindNeighbors(sub.obj, dims = 1:20)
  sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 0.7)
  
  p2 = DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 6) +
    ggtitle(paste0(celltype.sels, ' -- subclusters'))
  
  markers = FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3)
  
  markers %>%
    filter(!str_detect(gene, '^(AMEX|LOC|N/A)')) %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC) -> top10
  
  xx = subset(sub.obj, downsample = 500)
  #tops.common = c(markers.coarse$gene[!is.na(match(markers.coarse$cluster, c(1, 2, 3, 6, 13, 14, 15)))])
  
  p3 = DoHeatmap(xx, features = top10$gene) + NoLegend()
  p3
  #ggsave(filename = paste0(resDir, '/first_test_clusterMarkers_v2.pdf'), width = 10, height = 30)
  
  pdfname = paste0(resDir, '/subclustering_associatedMarkerGenes_', celltype.sels, '_v2.pdf')
  pdf(pdfname, width=16, height = 16)
  
  p1
  p2
  p3 
  dev.off()
  
  ## save the subclustering labels 
  cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 2)]
  mm = which(!is.na(match(colnames(aa), cell.sels)))
  cat(length(cell.sels), '--', length(mm), '\n')
  aa$subtypes[mm] = "CM.Endo.doublet"
  
  cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters != 2)]
  mm = which(!is.na(match(colnames(aa), cell.sels)))
  cat(length(cell.sels), '--', length(mm), '\n')
  
  
  ##########################################
  # subclustering FBs
  ##########################################
  celltype.sels = 'FB'
  sub.obj = subset(aa, cells = colnames(aa)[!is.na(match(aa$celltypes, celltype.sels))])
  
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
  
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 30)
  
  nb.pcs = 20 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 20; min.dist = 0.05;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, 
                     n.neighbors = n.neighbors,
                     min.dist = min.dist)
  
  
  # sub.obj$clusters = sub.obj$seurat_clusters
  p1 = DimPlot(sub.obj, group.by = 'clusters', reduction = 'umap', label = TRUE, label.size = 5) +
    ggtitle(celltype.sels)
  
  features = rownames(aa)[grep('VIM|COL1A2|FSTL1|POSTN', rownames(aa))]
  FeaturePlot(aa, features = features, cols = c('gray', 'red'))
  
  
  sub.obj <- FindNeighbors(sub.obj, dims = 1:20)
  sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 0.7)
  
  p2 = DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 5) +
    ggtitle(paste0(celltype.sels, ' -- subclusters'))
  
  p1 + p2
  
  markers = FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3)
  
  markers %>%
    filter(!str_detect(gene, '^(AMEX|LOC|N/A)')) %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC) -> top10
  
  xx = subset(sub.obj, downsample = 500)
  #tops.common = c(markers.coarse$gene[!is.na(match(markers.coarse$cluster, c(1, 2, 3, 6, 13, 14, 15)))])
  
  p3 = DoHeatmap(xx, features = top10$gene) + NoLegend()
  p3
  
  pdfname = paste0(resDir, '/subclustering_associatedMarkerGenes_', celltype.sels, '_v2.pdf')
  pdf(pdfname, width=18, height = 16)
  
  p1
  p2
  p3 
  dev.off()
  
  ## save the subclustering labels 
  cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 1)]
  mm = which(!is.na(match(colnames(aa), cell.sels)))
  cat(length(cell.sels), '--', length(mm), '\n')
  aa$subtypes[mm] = "FB/Endo"
  
  cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 8)]
  mm = which(!is.na(match(colnames(aa), cell.sels)))
  cat(length(cell.sels), '--', length(mm), '\n')
  aa$subtypes[mm] = "FB.CM.doublet"
  
    
  ##########################################
  # subclustering immune cells 
  ##########################################
  #celltype.sels = c('immune', 'blood', 'Macrophage', 'B')
  cluster.sels = c(24, 22, 18, 17, 9, 11, 20, 21, 13)
  sub.obj = subset(aa, cells = colnames(aa)[!is.na(match(aa$clusters, cluster.sels))])
  
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 5000)
  sub.obj = ScaleData(sub.obj)
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 30)
  
  nb.pcs = 20 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 20; min.dist = 0.1;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, 
                     n.neighbors = n.neighbors,
                     min.dist = min.dist)
  
  # sub.obj$clusters = sub.obj$seurat_clusters
  p1 = DimPlot(sub.obj, group.by = 'clusters', reduction = 'umap', label = TRUE, label.size = 5)
  p11 = DimPlot(sub.obj, group.by = 'subtypes', reduction = 'umap', label = TRUE, repel = TRUE, label.size = 5) + 
    NoLegend()
  
  p1
  p11
  
  features = rownames(sub.obj)[grep('PTPRC|CD68|CD8A|CD74|CSF1R|ITGAM', rownames(sub.obj))]
  FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))
  
  sub.obj <- FindNeighbors(sub.obj, dims = 1:20)
  sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 0.3)
  
  p2 = DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 5)
  
  p1 + p2
  
  markers = FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.3)
  
  markers %>%
    filter(!str_detect(gene, '^(AMEX|LOC|N/A)')) %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC) -> top10
  
  xx = subset(sub.obj, downsample = 500)
  #tops.common = c(markers.coarse$gene[!is.na(match(markers.coarse$cluster, c(1, 2, 3, 6, 13, 14, 15)))])
  
  p3 = DoHeatmap(xx, features = top10$gene) + NoLegend()
  p3
  
  pdfname = paste0(resDir, '/subclustering_associatedMarkerGenes_', celltype.sels, '_v2.pdf')
  pdf(pdfname, width=18, height = 16)
  
  p1
  p11
  p2
  p3 
  dev.off()
  
  ## save the subclustering labels
  Save.subtype.Annot = FALSE
  if(Save.subtype.Annot){
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 0)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = "Macrophage"
    
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 1)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = 'blood'
    
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 2)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = 'Thrombocytes/Megakaryocytes'
    
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 3)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = 'Bcell'
    
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 4)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = 'notSure.immune.CM.doublet'
    
    
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 5)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = 'immune/MF'
    
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 6)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = 'blood/neutrophil'
    
    cell.sels = colnames(sub.obj)[which(sub.obj$seurat_clusters == 7)]
    mm = which(!is.na(match(colnames(aa), cell.sels)))
    cat(length(cell.sels), '--', length(mm), '\n')
    aa$subtypes[mm] = 'Tcell'
    
  }
  
  ##########################################
  # subclustering Endo
  ##########################################
  celltype.sels = c('Endo')
  sub.obj = subset(aa, cells = colnames(aa)[!is.na(match(aa$celltypes, celltype.sels))])
  
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
  
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 30)
  
  nb.pcs = 20 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 20; min.dist = 0.05;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, 
                     n.neighbors = n.neighbors,
                     min.dist = min.dist)
  
  # sub.obj$clusters = sub.obj$seurat_clusters
  p1 = DimPlot(sub.obj, group.by = 'clusters', reduction = 'umap', label = TRUE, label.size = 5) +
    ggtitle(celltype.sels)
  
  #features = rownames(sub.obj)[grep('PTPRC|CD68|CD8A|CD74|CSF1R|ITGAM', rownames(sub.obj))]
  #FeaturePlot(sub.obj, features = features, cols = c('gray', 'red'))
  
  sub.obj <- FindNeighbors(sub.obj, dims = 1:10)
  sub.obj <- FindClusters(sub.obj, verbose = FALSE, algorithm = 3, resolution = 0.5)
  
  p2 = DimPlot(sub.obj, reduction = 'umap', label = TRUE, label.size = 5) +
    ggtitle(paste0(celltype.sels, ' -- subclusters'))
  
  p1 + p2
  
  markers = FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
  
  markers %>%
    filter(!str_detect(gene, '^(AMEX|LOC|N/A)')) %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC) -> top10
  
  xx = subset(sub.obj, downsample = 500)
  #tops.common = c(markers.coarse$gene[!is.na(match(markers.coarse$cluster, c(1, 2, 3, 6, 13, 14, 15)))])
  
  p3 = DoHeatmap(xx, features = top10$gene) + NoLegend()
  p3
  
  pdfname = paste0(resDir, '/subclustering_associatedMarkerGenes_', celltype.sels, '_v2.pdf')
  pdf(pdfname, width=18, height = 16)
  
  p1
  p2
  p3 
  dev.off()
  
  
}
