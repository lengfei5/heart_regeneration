##########################################################################
##########################################################################
# Project: heart regeneration   
# Script purpose: to have peak consensus with macs2/3
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Sep  2 10:07:15 2022
##########################################################################
##########################################################################
source('functions_scRNAseq.R')


for(n in 2:nrow(design))
{
  # n = 1
  cat('----------- : ', n, ':',  design$condition[n], '-------------\n')
  
  # load nf output and process
  topdir = paste0(dataDir, '/multiome_', design$timepoint[n], '/outs')
  
  counts <- Read10X_h5(filename = paste0(topdir, "/filtered_feature_bc_matrix.h5"))
  fragpath <- paste0(topdir, "/atac_fragments.tsv.gz")
  metadata <- read.csv(
    file = paste0(topdir, '/per_barcode_metrics.csv'),
    header = TRUE,
    row.names = 1
  )
  
  chrom_assay <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = granges_axolotl,
    min.cells = 10,
    min.features = 20
  )
  
  bb <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    
  )
  rm(chrom_assay, counts)
  
  bb$condition = design$condition[n]
  
  mm = match(colnames(bb), metadata$gex_barcode)
  metadata = metadata[, c(1,2, 3, 21:30)]
  
  bb = AddMetaData(bb, metadata = metadata)
  
  bb$pct_reads_in_peaks <- bb$atac_peak_region_fragments / bb$atac_fragments * 100
  bb$pct_usable_fragments = bb$atac_fragments/bb$atac_raw_reads
  
  DefaultAssay(bb) <- "ATAC"
  Idents(bb) = bb$condition
  
  # bb <- NucleosomeSignal(bb)
  # bb <- TSSEnrichment(bb, fast = FALSE)
 
  # bb$high.tss <- ifelse(bb$TSS.enrichment > 5, 'High', 'Low')
  
  #pdfname = paste0(resDir, '/scATAC_QCs_cellRangerOutput_', design$condition[n],  version.analysis,  '.pdf')
  #pdf(pdfname, width=12, height = 8)
  #TSSPlot(bb, group.by = 'high.tss') + NoLegend()
  
  VlnPlot(bb, features = "nCount_ATAC", ncol = 1, y.max = 5000) +
    geom_hline(yintercept = c(500, 1000))
  
  # VlnPlot(bb, 
  #         features = c("nCount_ATAC", "atac_fragments", "pct_reads_in_peaks"), 
  #         ncol = 3)
  # 
  # cat(median(bb$nCount_ATAC), '--', median(bb$atac_fragments), ' -- ', median(bb$pct_reads_in_peaks), ' \n')
  
  # VlnPlot(
  #   object = bb,
  #   features = c("TSS.enrichment", "nucleosome_signal"),
  #   ncol = 2,
  #   pt.size = 0
  # )
  
  # quick filtering 
  bb <- subset(
    x = bb,
    subset = nCount_ATAC > 200 &
      nCount_ATAC < 10000 
      #pct_reads_in_peaks > 15 &
      #blacklist_ratio < 0.05 &
      #nucleosome_signal < 4 &
      #TSS.enrichment > 1
  )
  
  bb = RunTFIDF(bb)
  bb = FindTopFeatures(bb, min.cutoff = "q5")
  bb = RunSVD(bb)
  
  DepthCor(bb, n = 30)
  #cordat = DepthCor(bb, n = 30)$data
  dims_use = c(2:30)
  print(dims_use)
  
  bb = FindNeighbors(object = bb, reduction = 'lsi', dims = dims_use, 
                     force.recalc = T, graph.name = "thegraph")
  bb = FindClusters(object = bb, verbose = FALSE, algorithm = 3, 
                    graph.name = "thegraph", resolution = 1)
  
  bb <- RunUMAP(object = bb, reduction = 'lsi', dims = dims_use)
  DimPlot(object = bb, label = TRUE) + NoLegend()
  
  bc_clusters = data.frame(bc = colnames(bb), 
                           cluster = paste0("cluster_", bb$seurat_clusters), 
  stringsAsFactors = FALSE)
  
  cell.nbs = table(bc_clusters$cluster)
  cell.nbs = cell.nbs[which(cell.nbs>200)]
  bc_clusters = bc_clusters[!is.na(match(bc_clusters$cluster, names(cell.nbs))), ]
  
  write.table(bc_clusters, file = paste0(topdir, '/cellBarcode_clusters_round_1.txt'), 
              sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  rm(bb)
  
}

##########################################
# round 1 peak combining from clusters of each samples 
##########################################
library(ChIPseeker)
library(rtracklayer)
library("ChIPpeakAnno")
library("ggplot2")
library("GenomicFeatures")


pval.cutoff = 5 # the 10^- seems to be better when checking peak overlapping in replicates

for(m in 1:nrow(design))
{
  # m = 1
  cat(m, '--', design$condition[m], '\n')
  
  peakDir = paste0(dataDir, '/multiome_', design$timepoint[m], '/outs/calledPeaks/macs2')
  peak.files = list.files(path = peakDir,
                          pattern = '*_peaks.xls', full.names = TRUE)
  cat(length(peak.files), ' peak files\n')
  
  for(n in 1:length(peak.files)) 
  {
    p = readPeakFile(peak.files[n], as = "GRanges");
    #eval(parse(text = paste0("p = pp.", k)));
    with.p.values = "X.log10.pvalue." %in% colnames(mcols(p))
    if(with.p.values) {
      p <- p[mcols(p)[,"X.log10.pvalue."] > pval.cutoff];
      p = reduce(p);
      #peaks10= c(peaks10, p10);
    }else{ 
      cat("no p values conlumn found for -- ", design.matrix$file.name[k], "\n");
      PLOT.p10 = FALSE;
    }
    cat('--- n = ', n, ' - peak numbers:  ', length(p), '\n')
    #p = reduce(p)
    if(m == 1 & n == 1) {
      peaks = p
    }else{
      peaks = union(peaks, p)
    }
  }
}

saveRDS(peaks, file = paste0(RdataDir, 'merged_peaks_macs2_firstRound_800K.rds'))

##########################################
# used the first round of merged peak to redo clustering  
##########################################
combined.peaks = readRDS(file = paste0(RdataDir, 'merged_peaks_macs2_firstRound_273K.rds'))

# Create a unified set of peaks to quantify in each dataset
# Filter out bad peaks based on length
peakwidths = width(combined.peaks)
combined.peaks = combined.peaks[peakwidths > 100]
# there is a problem with coordinates starting at 0 for some reason...
combined.peaks = restrict(combined.peaks, start = 1)
combined.peaks

library(future)
plan()
plan("multicore", workers = 1)
plan()


## Make all individual Seurats with common peaks
for(n in 2:nrow(design))
{
  # n = 1
  cat('----------- : ', n, ':',  design$condition[n], '-------------\n')
  
  # load nf output and process
  topdir = paste0(dataDir, '/multiome_', design$timepoint[n], '/outs')
  
  counts <- Read10X_h5(filename = paste0(topdir, "/filtered_feature_bc_matrix.h5"))
  
  fragpath <- paste0(topdir, "/atac_fragments.tsv.gz")
  frags_l = CreateFragmentObject(path = fragpath, cells = colnames(counts$Peaks))
 
  # slow step takes 16 mins without parall computation
  feat = FeatureMatrix(fragments = frags_l, features = combined.peaks, 
                       cells = colnames(counts$Peaks))
  
  chrom_assay <- CreateChromatinAssay(
    counts = feat,
    sep = c(":", "-"),
    fragments = frags_l,
    annotation = granges_axolotl,
    min.cells = 0,
    min.features = 0 
  )
  
  bb <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "ATAC",
    
  )
  rm(chrom_assay)
  
  bb$condition = design$condition[n]
  
  DefaultAssay(bb) <- "ATAC"
  Idents(bb) = bb$condition
  
  VlnPlot(bb, features = "nCount_ATAC", ncol = 1, y.max = 5000) +
    geom_hline(yintercept = c(500, 1000))
  
  # quick filtering 
  bb <- subset(
    x = bb,
    subset = nCount_ATAC > 100 &
      nCount_ATAC < 10000 
    #pct_reads_in_peaks > 15 &
    #blacklist_ratio < 0.05 &
    #nucleosome_signal < 4 &
    #TSS.enrichment > 1
  )
  
  bb = RunTFIDF(bb, method = 3)
  bb = FindTopFeatures(bb, min.cutoff = "q5")
  bb = RunSVD(bb)
  
  DepthCor(bb, n = 30)
  #cordat = DepthCor(bb, n = 30)$data
  dims_use = c(2:30)
  print(dims_use)
  
  bb = FindNeighbors(object = bb, reduction = 'lsi', dims = dims_use, 
                     force.recalc = T, graph.name = "thegraph")
  bb = FindClusters(object = bb, verbose = FALSE, algorithm = 3, 
                    graph.name = "thegraph", resolution = 1)
  
  bb <- RunUMAP(object = bb, reduction = 'lsi', dims = dims_use)
  DimPlot(object = bb, label = TRUE) + NoLegend()
  
  bc_clusters = data.frame(bc = colnames(bb), 
                           cluster = paste0("cluster_", bb$seurat_clusters), 
                           stringsAsFactors = FALSE)
  
  cell.nbs = table(bc_clusters$cluster)
  cell.nbs = cell.nbs[which(cell.nbs>200)]
  bc_clusters = bc_clusters[!is.na(match(bc_clusters$cluster, names(cell.nbs))), ]
  
  write.table(bc_clusters, file = paste0(topdir, '/cellBarcode_clusters_round_2.txt'), 
              sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  rm(bb)
  
}

##########################################
# round II peak combining from clusters of each samples 
##########################################
library(ChIPseeker)
library(rtracklayer)
library("ChIPpeakAnno")
library("ggplot2")
library("GenomicFeatures")

pval.cutoff = 5 # the 10^-6 seems to be better when checking peak overlapping in replicates

for(m in 1:nrow(design))
{
  # m = 1
  cat(m, '--', design$condition[m], '\n')
  
  peakDir = paste0(dataDir, '/multiome_', design$timepoint[m], '/outs/calledPeaks/macs2')
  peak.files = list.files(path = peakDir,
                          pattern = '*_peaks.xls', full.names = TRUE)
  cat(length(peak.files), ' peak files\n')
  
  for(n in 1:length(peak.files)) 
  {
    p = readPeakFile(peak.files[n], as = "GRanges");
    #eval(parse(text = paste0("p = pp.", k)));
    with.p.values = "X.log10.pvalue." %in% colnames(mcols(p))
    if(with.p.values) {
      p <- p[mcols(p)[,"X.log10.pvalue."] > pval.cutoff];
      p = reduce(p);
      #peaks10= c(peaks10, p10);
    }else{ 
      cat("no p values conlumn found for -- ", design.matrix$file.name[k], "\n");
      PLOT.p10 = FALSE;
    }
    cat('--- n = ', n, ' - peak numbers:  ', length(p), '\n')
    #p = reduce(p)
    if(m == 1 & n == 1) {
      peaks = p
    }else{
      peaks = union(peaks, p)
    }
  }
}

length(peaks)

saveRDS(peaks, file = paste0(RdataDir, 'merged_peaks_macs2_secondRound_800K_pval5.rds'))




