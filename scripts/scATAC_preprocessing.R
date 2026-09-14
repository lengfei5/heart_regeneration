##########################################################################
##########################################################################
# Project: Aleks' single-cell project
# Script purpose: mainly cell filtering for scATAC-seq data for peak- and window-basded matrix
# and also does crude clustering for windown-based matrix
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Mar 25 11:41:30 2020
##########################################################################
##########################################################################
source.my.script <- function(name.of.function){
  tryCatch(path <- rstudioapi::getSourceEditorContext()$path,
           error = function(e){
             install.packages("rstudioapi")
             path <-  rstudioapi::getSourceEditorContext()$path})
  source.path <- sub(basename(path), "", path)
  source(paste0(source.path,name.of.function))
}

## set up the paths for the data and results
tryCatch(path <- rstudioapi::getSourceEditorContext()$path, 
         error = function(e){
           install.packages("rstudioapi")
           path <-  rstudioapi::getSourceEditorContext()$path})
source.path <- sub(basename(path), "", path)

source.my.script('scATAC_functions.R')
if(!require("DropletUtils")) BiocManager::install('DropletUtils')
library(DropletUtils)
library(data.table)
library(Matrix)
library(tictoc)

## in this script the results are saved in the 
scATACpro_output_path = '~/workspace/imp/scRNAseq_MS_lineage_dev/output_cellranger.ce11_scATACpro'

setwd(scATACpro_output_path)

raw_mtx_file = paste0(scATACpro_output_path, '/raw_matrix/MACS2/matrix.mtx')
filtered_mtx_dir = paste0(scATACpro_output_path, '/filtered_matrix_peaks_barcodes')
raw_mtx_dir = dirname(raw_mtx_file)

if(!dir.exists(filtered_mtx_dir)){dir.create(filtered_mtx_dir)}

Processing.window.based.matrix = FALSE
########################################################
########################################################
# Section : prepare annotation files
# 
########################################################
########################################################
# # args = commandArgs(T)
Generate.enhancer.annotation = FALSE
if(Generate.enhancer.annotation){
  enhancers = read.table('/Volumes/groups/cochella/jiwang/annotations/promoter_enhancer_annotation-elife-37344-fig2-data1-v2.txt', header = TRUE)
  enhancers = enhancers[ which(enhancers$annot == 'putative_enhancer'), ]
  enhancers = data.frame(enhancers$chrom_ce11, enhancers$start_ce11, enhancers$end_ce11, enhancers$associated_gene_id, stringsAsFactors = FALSE)
  enhancers$scores = '.'
  enhancers$strand = '.'
  write.table(enhancers, file = '/Volumes/groups/cochella/jiwang/annotations/ce11_enhancers_Jaenes_et_al_elife_2018.bed',
                        row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
}

########################################################
########################################################
# Section : Cell filtering using DropletUtils
#  prepare filtered matrix, barcodes and peaks
########################################################
########################################################
mat = readMM(raw_mtx_file)
features = fread(paste0(raw_mtx_dir, '/features.txt'), header = F)
barcodes = fread(paste0(raw_mtx_dir, '/barcodes.txt'), header = F)
rownames(mat) = features$V1
colnames(mat) = barcodes$V1
# metrics = read.csv('res_counts/outs/singlecell.csv', header = TRUE)

set.seed(2019)
tic('emptyDrops running :')
cell.out <- emptyDrops(mat, BPPARAM = SerialParam())
toc()

PLOT.barcode.rank = TRUE
if(PLOT.barcode.rank){
  
  pdfname = paste0(filtered_mtx_dir, "/EmptyDrop_barcode_filtering.pdf")
  pdf(pdfname, width=14, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  br.out <- barcodeRanks(mat)
  
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))

  cell.out
  
  is.cell <- cell.out$PValue < 0.0001
  sum(is.cell, na.rm=TRUE)
  
  table(Limited=cell.out$Limited, Significant=is.cell)
  plot(cell.out$Total, -cell.out$LogProb, col=ifelse(is.cell, "red", "black"),
       xlab="Total UMI count", ylab="-Log Probability", log='x')
  
  dev.off()
  
}

filter.out <- cell.out[complete.cases(cell.out), ]

saveRDS(filter.out, file = paste0(filtered_mtx_dir, '/EmptyDrop_obj.rds'))


fdr = 0.01
#filter.out = filter.out[filter.out$FDR <= fdr, ]
kk = which(filter.out$Total>10^3 &  
              filter.out$PValue < 0.001)
length(kk)

filter.out = filter.out[kk, ]
select.cells = rownames(filter.out)
#length(select.cells)

out_mat = mat[, colnames(mat) %in% select.cells]
barcodes = colnames(out_mat)

peaks = rownames(out_mat)

## for peak.bed file
peaks = t(sapply(peaks, FUN = function(x){
  ccord = unlist(strsplit(as.character(x), "-"))
  if(!startsWith(ccord[1], 'chr')){
    return(c(paste0('chr', ccord[1]), ccord[-1]))
  }else{
    return(ccord)
  }  
}))

#metrics$barcode = sapply(metrics$barcode, function(x) gsub('-1', '', x))
#mm = match(barcodes, metrics$barcode)
#metrics = metrics[]

## save filtered matrix, peaks and barcodes
system(paste('mkdir -p', filtered_mtx_dir))

writeMM(out_mat, file = paste0(filtered_mtx_dir, '/matrix.mtx'))  
write.table(barcodes, file = paste0(filtered_mtx_dir, '/barcodes.tsv'), sep = '\t', 
            row.names = F, quote = F, col.names = F)
write.table(peaks, file = paste0(filtered_mtx_dir, '/peaks.bed'), sep = '\t',
            row.names = F, quote = F, col.names = F)

########################################################
########################################################
# Section :  lsi and lsi_log transform for window-based filtered matrix
# 
########################################################
########################################################
if(Processing.window.based.matrix){
  source.my.script('scATAC_functions.R')
  tenx.bmat = load_tenx_atac(paste0(filtered_mtx_dir, '/matrix.mtx'), 
                             paste0(filtered_mtx_dir, '/peaks.bed'), 
                             paste0(filtered_mtx_dir, '/barcodes.tsv'))
  
  if(Processing.window.based.matrix){
    # select the top 20k features for window-based matrix
    ss = Matrix::rowSums(tenx.bmat)
    ss.rank = order(ss)
    sels = ss.rank <=20000
    tenx.bmat = tenx.bmat[sels, ]
    
  }else{
    # select features showing in > 50 cells
    tenx.bmat = filter_features(tenx.bmat, cells=50)
  }
  
  
  # Binarize the matrix for consistency
  tenx.bmat@x[tenx.bmat@x > 1] = 1
  
  tenx.seurat.lsi = lsi_workflow(tenx.bmat, dims=2:50, log_scale_tf=FALSE, 
                                 reduction='pca', resolution=0.8)
  tenx.seurat.lsi_log = lsi_workflow(tenx.bmat, dims=2:50, log_scale_tf=TRUE, 
                                     reduction='pca.l2', resolution=0.8)
  
  plot_clustering_comparison(tenx.seurat.lsi,
                             tenx.seurat.lsi_log,
                             reduction='umap',
                             description1='LSI',
                             description2='LSI logTF',
                             cluster_column1='peaks_snn_res.0.8',
                             cluster_column2='peaks_snn_res.0.8')
  
  
  ##########################################
  # select the transformation methods (lsi or lsi_log)
  # tune the resolution for clustering 
  # and save the crude clusters for making bam files and the peak calling in each cluster 
  ##########################################
  p1 = DimPlot(object = tenx.seurat.lsi, label = TRUE, reduction = 'umap') + NoLegend()
  p2 = DimPlot(object = tenx.seurat.lsi_log, label = TRUE, reduction = 'umap') + NoLegend()
  
  p1 + p2
  
  ## refine parameters
  tenx.seurat = tenx.seurat.lsi # select lsi transformation
  
  tenx.seurat = FindNeighbors(object = tenx.seurat, 
                                  reduction='pca', nn.eps=0.25, dims=2:50)
  tenx.seurat = FindClusters(object = tenx.seurat, 
                                 n.start=20, resolution=1.5,
                                 algorithm = 3)
  
  nb.pcs = 40; n.neighbors = 20; min.dist = 0.2;
  tenx.seurat <- RunUMAP(object = tenx.seurat, 
                             reduction = 'pca', 
                             dims = 2:nb.pcs, 
                             n.neighbors = n.neighbors, 
                             min.dist = min.dist)
  DimPlot(object = tenx.seurat, label = TRUE, reduction = 'umap') + 
    NoLegend()
  
  table(tenx.seurat$peaks_snn_res.1.5)
  tenx.seurat = SetIdent(tenx.seurat, value='peaks_snn_res.1.5')
  
  ##########################################
  # generate bigwigs for each cluster for visualization
  ##########################################
  barcode.cluster = data.frame(Barcode = colnames(tenx.seurat),
                               Cluster = tenx.seurat$peaks_snn_res.1.5, 
                               stringsAsFactors = FALSE)
  # manually merge small clusters with nb.cells <200
  barcode.cluster$Cluster[which(barcode.cluster$Cluster == '14')] = '0'
  
  write.table(barcode.cluster, file = paste0(filtered_mtx_dir, '/barcodes_and_clusters.txt'), 
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
  
  
  ########################################################
  ########################################################
  # Section : Merge peaks for peaks files in crude clusters 
  # 
  ########################################################
  ########################################################
  Merge.peaks.for.each.clusters = FALSE
  if(Merge.peaks.for.each.clusters){
    
    library(bedr)
    library(data.table)
    library(Signac) # blacklist_ce11 is included in this package
    
    peak_cluster_dir = '/Volumes/sshftps/scratch/jiwang/scATAC_pro_test/output/preprocess_bin2kb/peaks_by_cluster'
    
    files = dir(peak_cluster_dir)
    files = files[grepl(files, pattern = "narrowPeak")]
    
    peaks = NULL
    for(file0 in files){
      peaks = rbind(peaks, fread(paste0(peak_cluster_dir, '/', file0)))
    }
    
    # xx = peaks[!grepl(V1, pattern = 'random', ignore.case = T)]
    peaks = peaks[!grepl(V1, pattern = 'random', ignore.case = T)]
    peaks = peaks[!grepl(V1, pattern = 'Un', ignore.case = T)]
    peaks = peaks[!grepl(V1, pattern = 'M', ignore.case = T)]
    
    regions = paste0(peaks$V1, ':', peaks$V2, '-', peaks$V3)
    
    blacklist = as.data.frame(blacklist_ce11)
    blackregions = paste0(blacklist$seqnames, ':', 
                          blacklist$start, '-', 
                          blacklist$end)
    
    #xx = is.valid.region(regions, throw.error = TRUE, check.chr = FALSE)
    a.sort   <- bedr.sort.region(regions, check.chr = FALSE,
                                 chr.to.num = c('chrI' = 1, 
                                                'chrII' = 2, 
                                                'chrIII' = 3, 
                                                'chrIV' = 4,
                                                'chrV' = 5, 
                                                'chrX' = 6))
    b.sort   <- bedr.sort.region(blackregions, check.chr = FALSE,
                                 chr.to.num = c('chrI' = 1, 
                                                'chrII' = 2, 
                                                'chrIII' = 3, 
                                                'chrIV' = 4,
                                                'chrV' = 5, 
                                                'chrX' = 6))
    
    a.merged <- bedr.merge.region(a.sort, distance = 0, check.chr = FALSE)
    
    
    if (check.binary("bedtools")) {
      #a <- bedr(engine = "bedtools", input = list(i = a), method = "sort", params = "");
      #b <- bedr(engine = "bedtools", input = list(i = b), method = "sort", params = "");
      #d <- bedr.subtract.region(a,b);
      dd = bedr.subtract.region(a.merged, b.sort, 
                                remove.whole.feature = TRUE,
                                check.chr = FALSE)
    }
    
    
    
    dd = data.table('peak' = dd)
    dd[, 'chr' := unlist(strsplit(peak, ':'))[1], by = peak]
    dd[, 'range' := unlist(strsplit(peak, ':'))[2], by = peak]
    dd[, 'start' := unlist(strsplit(range, '-'))[1], by = peak]
    dd[, 'end' := unlist(strsplit(range, '-'))[2], by = peak]
    
    dd[, c('peak', 'range') := NULL]
    
    write.table(dd, file = paste0(peak_cluster_dir, '/merged_peaks.bed'), quote = F, 
                row.names = F, col.names = F, sep = '\t')
    
  }
  
}
