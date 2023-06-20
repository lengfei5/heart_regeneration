##########################################################################
##########################################################################
# Project: heart regeration
# Script purpose: test trajectory
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Apr 18 13:07:41 2023
##########################################################################
##########################################################################
source('functions_scATAC.R')
source('functions_scRNAseq.R')
source('functions_Visium.R')
#source('functions_scRNAseq.R')

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(patchwork)
require(SeuratObject)
library(data.table)
source('functions_scATAC.R')
library(ArchR)

library(JASPAR2020)
library(TFBSTools)
library(chromVAR)

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
require(tictoc)
library(future)
require(pheatmap)
require(RColorBrewer)
options(future.globals.maxSize = 80 * 1024^3)

mem_used()

###  import annotation and metadata
set.seed(1234)
library(BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M)
library(ballgown)
gtf_axolotl = paste0("/groups/tanaka/People/current/jiwang/scripts/axolotl_multiome/r_package/", 
                     "AmexT_v47.FULL_corr_chr_cut.gtf")

granges_axolotl = ballgown::gffReadGR(gtf_axolotl)

# adding a gene biotype, as that's necessary for TSS metaprofile
granges_axolotl$gene_biotype = "protein_coding"

species = 'axloltl_scATAC'

outDir = paste0(resDir, '/CM_trajectory_test/')
system(paste0('mkdir -p ', outDir))


##########################################
# import the multiome data
##########################################
aa = readRDS(file = paste0('../results/sc_multiome_R13591_atac_reseq_20221115/Rdata/',
                           'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                           'scATAC.merged.peaks.cr_filtered_umap.lsi',
                           '584K.features_37680cells_umap.topics_updated.umap.subtypes_celltypes.rds'))

aa$time = gsub('Amex_', '', aa$condition)
aa$cell.ids = sapply(colnames(aa), function(x) unlist(strsplit(as.character(x), '-'))[1]) 
aa$cell.ids = paste0(aa$cell.ids, '_', aa$time)


# identify the DARs using the celltypes 
DefaultAssay(aa) <- 'ATAC'

aa = subset(aa, cells = colnames(aa)[which(aa$celltypes != 'Neuronal')])
Idents(aa) = aa$celltypes

motif_tf = readRDS(file = paste0('../results/sc_multiome_R13591_atac_reseq_20221115/Rdata/', 
                                 'motif_to_tfs_pfm_JASPAR2020_CORE_vertebrate_v1.rds'))
chromvar = readRDS(file = paste0('../results/sc_multiome_R13591_atac_reseq_20221115/Rdata/', 
                                 'atac_seuratObject_motifClass_chromVAR_v3.rds'))
DefaultAssay(chromvar) <- 'chromvar'
ss = colSums(chromvar@assays$chromvar@data)
length(which(is.na(ss)))
data = chromvar@assays$chromvar@data
data[which(is.na(data))] = 0
chromvar@assays$chromvar@data = data

##########################################
# subset CM sels 
##########################################
celltype_sel = 'CM'

#sub_obj = subset(aa, cells = colnames(aa)[which(aa$celltypes == celltype_sel)])
#sub_obj$subtypes = droplevels(sub_obj$subtypes)

## make sure the CM subtypes are the latest version
ref =  readRDS(file = paste0("/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/", 
                             "CM_subset_for_velocity.rds"))

mm = match(colnames(ref), colnames(aa))
cat(length(which(is.na(mm))), ' cells lost in the updated CM annotations \n')

cells_keep = intersect(colnames(ref), colnames(aa))
cat(length(cells_keep), ' cells to consider \n')

sub_obj = subset(aa, cells = cells_keep)

sub_obj$subtypes = ref$subtypes[match(colnames(sub_obj), colnames(ref))]

sub_obj$subtypes = droplevels(sub_obj$subtypes)

Idents(sub_obj) <- factor(sub_obj$subtypes)

DimPlot(ref, group.by = 'subtypes')
DimPlot(sub_obj, group.by = 'subtypes')

umap.embedding = ref@reductions$umap@cell.embeddings
umap.embedding = umap.embedding[match(colnames(sub_obj), rownames(umap.embedding)), ]

sub_obj[['umap']] = Seurat::CreateDimReducObject(embeddings=umap.embedding,
                                                 key='UMAP_',
                                                 assay='RNA')
rm(umap.embedding)
DefaultAssay(sub_obj) <- 'RNA'

DimPlot(sub_obj, group.by = 'subtypes')

saveRDS(sub_obj, file = paste0(outDir, 'CM_subset_RNA_ATAC_4trajectory.pseudotime.rds'))

########################################################
########################################################
# Section I : estimate pseudotime using scanpy or dpt
# 
########################################################
########################################################
Estimate_pseudotime_Scanpy = FALSE
if(Estimate_pseudotime_Scanpy){
  subtype_version = '_injurySubtypes'
  
  sel_subtypes = c('CM_IS', 'CM_ven_(Cav3_1)', 'CM_Prol_IS')
  # sub_obj = subset(sub_obj, cells = colnames(sub_obj)[which(!is.na(match(sub_obj$subtypes, sel_subtypes)))])
  
  sub_obj = subset(sub_obj, cells = colnames(sub_obj)[which(!is.na(match(sub_obj$subtypes, sel_subtypes)))])
  sub_obj$subtypes = droplevels(sub_obj$subtypes)
  
  table(sub_obj$subtypes)
  
  refine_subpopulation_for_trajectory = FALSE
  if(refine_subpopulation_for_trajectory){
    
    #subtypes_sel = c("EC", "EC_IS_(IARS1)", "EC_IS_(LOX)", "EC_IS_Prol")
    subtypes_sel = setdiff(unique(sub_obj$subtypes), c("EC_IS_(IARS1)", "EC_IS_(LOX)", "EC_IS_Prol"))
    
    sub_obj = subset(sub_obj,  cells = colnames(sub_obj)[which(!is.na(match(sub_obj$subtypes, subtypes_sel)))])
    
    sub_obj$subtypes = droplevels(sub_obj$subtypes)
    
    ## redo the clustering in case needed in the downstream analysis
    sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000)
    sub_obj <- ScaleData(sub_obj)
    sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), weight.by.var = TRUE, 
                      verbose = FALSE)
    ElbowPlot(sub_obj, ndims = 50)
    
    sub_obj <- RunUMAP(sub_obj, dims = 1:10, n.neighbors = 30, min.dist = 0.1)
    
    DimPlot(sub_obj, group.by = 'subtypes', label = TRUE, repel = TRUE)
    
  }else{
    ## redo the clustering in case needed in the downstream analysis
    sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 1000)
    sub_obj <- ScaleData(sub_obj)
    sub_obj <- RunPCA(sub_obj, features = VariableFeatures(object = sub_obj), weight.by.var = TRUE, 
                      verbose = FALSE)
    ElbowPlot(sub_obj, ndims = 50)
    
    sub_obj <- RunUMAP(sub_obj, dims = 1:10, n.neighbors = 30, min.dist = 0.3)
    DimPlot(sub_obj, group.by = 'subtypes', label = TRUE, repel = TRUE)
    
  }
  
  sub_obj <- FindNeighbors(sub_obj, dims = 1:10)
  sub_obj <- FindClusters(sub_obj, verbose = FALSE, algorithm = 3, resolution = 0.3)
  
  p1 = DimPlot(sub_obj, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()
  p2 = DimPlot(sub_obj, label = TRUE, repel = TRUE, group.by = 'seurat_clusters', raster=FALSE)
  
  p1 + p2
  
  sub_obj$time = gsub('d', '', sub_obj$time)
  sub_obj$clusters = sub_obj$seurat_clusters
  
  ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_subset_', celltype_sel, subtype_version,
                           '_reclustered.pdf'), 
         height = 5, width = 14)
  
  # DefaultAssay(sub_obj) <- 'RNA'
  # p1 = DimPlot(sub_obj, label = TRUE, group.by = 'subtypes',  repel = TRUE) + NoLegend()
  # 
  # DefaultAssay(sub_obj) <- 'ATAC'
  # p2 = DimPlot(sub_obj, label = TRUE, reduction = 'umap_topics',
  #              group.by = 'subtypes',  repel = TRUE) + NoLegend()
  # p1 + p2
  # 
  # ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_subset_', celltype_sel,  subtype_version, '.pdf'), 
  #        height = 6, width = 14)
  # 
  # DefaultAssay(sub_obj) <- 'RNA'
  # DimPlot(sub_obj, label = TRUE, group.by = 'subtypes', split.by = 'condition', repel = TRUE) + NoLegend()
  # 
  # ggsave(filename = paste0(outDir, '/multiome_snRNA_scATAC_subset_', celltype_sel, subtype_version, 
  #                          '_bytimePoint.pdf'), 
  #        height = 5, width = 20)
  
  
  ##########################################
  # preapre the spliced and unspliced matrix   
  ##########################################
  DefaultAssay(sub_obj) = 'RNA'
  
  #saveRDS(sub_obj, file = paste0(outDir, 'seuratObj_multiome_snRNA_scATAC_',  celltype_sel, 
  #                               subtype_version, '.rds'))
  
  source('utility_velocity.R')
  mnt = preapre_dataFile_for_RNAvelocity_PAGA(seuratObj = sub_obj)
  
  Idents(mnt) = mnt$condition
  table(mnt$condition)
  
  #mnt = subset(mnt, downsample = 2000)
  
  saveFile = paste0('RNAmatrix_umap_kalisto.velocity_spliced_unspliced_',
                    'CM_subtypes', subtype_version, '_timepoints.all_downsample.h5Seurat')
  
  SaveH5Seurat(mnt, filename = paste0(outDir, saveFile), 
               overwrite = TRUE)
  Convert(paste0(outDir, saveFile), 
          dest = "h5ad", overwrite = TRUE)
  
  
}

##########################################
# 3D DM to visualize the trajectory
##########################################
Test_Pseudotime_DPT_Slingshot = FALSE
if(Test_Pseudotime_DPT_Slingshot){
  library(slingshot, quietly = FALSE)
  library(destiny, quietly = TRUE)
  library(mclust, quietly = TRUE)
  library(scater)
  library(SingleCellExperiment)
  library(scran)
  library(RColorBrewer)
  
  DefaultAssay(sub_obj) = 'RNA'
  sub_obj_diet = DietSeurat(sub_obj, 
                            counts = TRUE, 
                            data = TRUE,
                            scale.data = TRUE,
                            features = rownames(sub_obj), 
                            assays = c('RNA'), 
                            dimreducs = c('umap'), 
                            graphs = NULL, 
                            misc = TRUE)
  sce = as.SingleCellExperiment(sub_obj_diet)
  
  rm(sub_obj_diet)
  
  dec <- modelGeneVar(sce)
  
  #plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
  #curve(metadata(dec)$trend(x), col="blue", add=TRUE)
  
  nb_features = 3000; 
  top.hvgs <- getTopHVGs(dec, n=nb_features)
  
  sce <- runPCA(sce, subset_row=top.hvgs, ncomponents = 100)
  # reducedDimNames(sce)
  ll.pca = reducedDim(sce, 'PCA')[, c(1:50)]
  
  n_neighbors = 100; n_eigs = 50; sigma = 'global';
  
  tic()
  dm <- DiffusionMap(ll.pca, sigma = sigma, k = n_neighbors, n_eigs = n_eigs, distance = 'euclidean')
  toc()
  
  cells = names(dm$DC1)
  metadata = sub_obj@meta.data
  dcs = data.frame(dm@eigenvectors, stringsAsFactors = FALSE)
  
  dcs = dcs[match(rownames(metadata), cells), ]
  
  dcs = as.matrix(dcs)
  
  sub_obj[["DC"]] <- CreateDimReducObject(embeddings = as.matrix(dcs), key = "DC_", 
                                          assay = DefaultAssay(sub_obj))
  
  rm(metadata)
  
  #sub_obj = RunUMAP(sub_obj, reduction = "DC", dims = 1:30, n.neighbors = 100, min.dist = 0.3, 
  # metric = 'euclidean',
  #             reduction.name = "dc_umap")
  #DimPlot(sub_obj, reduction = 'dc_umap', label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE)
  
  ## try to make 3d for DC1, DC2 and DC3
  ## example from https://plotly-r.com/d-charts.html
  library(plotly)
  dcs = as.data.frame(dcs)
  dcs$subtypes = sub_obj$subtypes[match(rownames(dcs), colnames(sub_obj))]
  dcs$condition = sub_obj$condition[match(rownames(dcs), colnames(sub_obj))]
  
  
  plot_ly(data.frame(dcs), x = ~DC1, y = ~DC2, z = ~DC3, size = 3) %>%
    add_markers(color = ~ subtypes)
  
  plot_ly(data.frame(dcs), x = ~DC1, y = ~DC2, z = ~DC4, size = 3) %>%
    add_markers(color = ~ subtypes)
  
  plot_ly(data.frame(dcs), x = ~DC1, y = ~DC3, z = ~DC4, size = 3) %>%
    add_markers(color = ~ subtypes)
  
  
  p1 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
               group.by = 'subtypes')
  p2 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 3), label = TRUE, repel = TRUE,
               group.by = 'subtypes')
  
  p3 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 4), label = TRUE, repel = TRUE,
               group.by = 'subtypes')
  
  p4 = DimPlot(sub_obj, reduction = 'DC', dims = c(3, 4), label = TRUE, repel = TRUE,
               group.by = 'subtypes')
  
  (p1 / p2)|(p3 / p4)
  
  ggsave(filename = paste0(outDir, 'trajectory_test_', celltype_sel, subtype_version, '.pdf'), 
         width = 12, height = 16)
  
  p1 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
               group.by = 'condition')
  p2 = DimPlot(sub_obj, reduction = 'DC', dims = c(1, 3), label = TRUE, repel = TRUE,
               group.by = 'condition')
  p1/p2
  
  DimPlot(sub_obj, reduction = 'DC', dims = c(1, 2), label = TRUE, repel = TRUE,
          group.by = 'subtypes', split.by = 'condition')
  ggsave(filename = paste0(outDir, 'trajectory_test_', celltype_sel, subtype_version, '_DC_split.by.time.pdf'), 
         width = 6, height = 18)
  
  save(sub_obj, dcs, file = paste0(outDir, 'trajectory_test_DM_', 
                                   celltype_sel, subtype_version, 
                                   '.Rdata'))
}


########################################################
########################################################
# Section II: use the scanpy pseudotime 
# - gene expression (in parituclar TFs or SPs), 
# - enhancers, 
# - motif activity along the pseudotime
########################################################
########################################################
sub_obj = readRDS(file = paste0(outDir, 'CM_subset_RNA_ATAC_4trajectory.pseudotime.rds'))

##########################################
# add pseudotime from scanpy
##########################################
sub_obj$pst = NA 
sub_obj$branch = NA
xx = read.csv(file = paste0(outDir, 
                            'CM_subtypes_injurySubtypes_timepoints.all_downsample/scanpy_dpt_pseudotime.csv'))

mm = match(xx$X, colnames(sub_obj))
sub_obj$branch[mm] = 'injury'
sub_obj$pst[mm] = xx$dpt_pseudotime

xx = read.csv(file = paste0(outDir, 
                            'CM_subtypes_uninjurySubtypes_timepoints.all_downsample/scanpy_dpt_pseudotime.csv'))

mm = match(xx$X, colnames(sub_obj))
sub_obj$branch[mm] = 'uninjury'
sub_obj$pst[mm] = xx$dpt_pseudotime

sub_obj$pst = as.numeric(sub_obj$pst)

p1 = DimPlot(sub_obj, group.by = 'subtypes')
p2 = DimPlot(sub_obj, group.by = 'branch')
p3 = FeaturePlot(sub_obj, features = 'pst')

p1 + p2 +p3

ggsave(filename = paste0(outDir, 'dpt_pseudotime_CM.pdf'), 
       width = 18, height = 6)

saveRDS(sub_obj, file = paste0(outDir, 'CM_subset_RNA_ATAC_pseudotime.scanpy.rds'))


##########################################
# test trajectory genes along pseudotime  
# 1) tradeseq (not sure it will work)
##########################################
aa = readRDS(file = paste0(outDir, 'CM_subset_RNA_ATAC_pseudotime.scanpy.rds'))

Test_TradeSeq = FALSE
if(Test_TradeSeq){
  library(tradeSeq)
  library(RColorBrewer)
  library(SingleCellExperiment)
  library(Seurat)
  library(SeuratDisk)
  library(SeuratObject)
  library(slingshot)
  palette(brewer.pal(8, "Dark2"))
  
  # data(countMatrix, package = "tradeSeq")
  # counts <- as.matrix(countMatrix)
  # rm(countMatrix)
  # data(crv, package = "tradeSeq")
  # data(celltype, package = "tradeSeq")
  # 
  # set.seed(7)
  # pseudotime <- slingPseudotime(crv, na = FALSE)
  # cellWeights <- slingCurveWeights(crv)
  # 
  
  # assign cell weights to two trajectories
  index_cells = which(!is.na(aa$branch))
  
  cellWeights = matrix(NA, ncol = 2, nrow = length(index_cells))
  rownames(cellWeights) = colnames(aa)[index_cells]
  colnames(cellWeights) = c('injury', 'uninjury')
  
  pseudotime = matrix(NA, ncol = 2, nrow = nrow(cellWeights))
  rownames(pseudotime) = rownames(cellWeights)
  colnames(pseudotime) = colnames(cellWeights)
  
  cell_injury = colnames(aa)[which(aa$branch == 'injury')]
  cell_uninjury = colnames(aa)[which(aa$branch == 'uninjury')]
  
  kk1 = match(cell_injury, rownames(cellWeights))
  cellWeights[kk1, 1] = 1; cellWeights[kk1, 2] = 0
  pseudotime[kk1, 1] = aa$pst[match(cell_injury, colnames(aa))] # keep the pseudo computed before
  #pseudotime[kk1, 2] = pseudotime[kk1, 1]
  pseudotime[kk1, 2] = 0
  
  kk2 = match(cell_uninjury, rownames(cellWeights))
  cellWeights[kk2, 1] = 0; cellWeights[kk2, 2] = 1
  pseudotime[kk2, 2] = aa$pst[match(cell_uninjury, colnames(aa))] # keep the pseudo computed before
  #pseudotime[kk2, 1] = pseudotime[kk2, 2]
  pseudotime[kk2, 1] = 0
  
  rm(list =c('kk1', 'kk2'))
  
  counts = GetAssayData(aa, slot = 'counts', assay = 'RNA')
  counts = counts[, match(rownames(pseudotime), colnames(counts))]
  
  #save(counts, pseudotime, cellWeights, 
  #     file = paste0(outDir, '/counts_pseudotime_cellWeights_for_tradeSeq.Rdata'))
  #load(file = paste0(outDir, '/counts_pseudotime_cellWeights_for_tradeSeq.Rdata'))
  
  #RNGversion("3.5.0")
  counts <- as.matrix(counts)
    
  ### Downstream of any trajectory inference method using pseudotime and cell weights
  # slow, but still ok. with 60K cells, it takes 10-15 mins for each knot.
  require(tictoc)
  library(BiocParallel)
  parallel::detectCores()
  
  BPPARAM <- BiocParallel::bpparam()
  BPPARAM # lists current options
  BPPARAM$workers <- 1 # use 2 cores
  
  tic()
  set.seed(7)
  #index_sub = sample(c(1:ncol(counts)), 10000)
  #set.seed(7)
  icMat <- evaluateK(counts = counts,
                      pseudotime = pseudotime,
                      cellWeights = cellWeights,
                      k=seq(5, 100, by = 5),
                      nGenes = 300,
                      verbose = TRUE,
                      plot = TRUE,
                      parallel=FALSE,
                     BPPARAM = BPPARAM)
  
  saveRDS(icMat, file = paste0(outDir, 'tradeSeqDE_evaluateK_res_v2.rds'))

  toc()
  
  
  pdfname = paste0(outDir, "Knot_number_estimation_genes300.pdf")
  pdf(pdfname, width=12, height = 10)
  par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  #icMat = readRDS(file = paste0(outDir, 'tradeSeqDE_evaluateK_res.rds'))
  #icMat = readRDS(file = paste0(outDir, 'tradeSeqDE_evaluateK_res_k3.12_100genes_v7.rds'))
  #icMat = readRDS(file = paste0(outDir, 'tradeSeqDE_evaluateK_res_k3.12_200genes_v6.rds'))
  icMat = readRDS(file = paste0(outDir, 'tradeSeqDE_evaluateK_res_k3.12_300genes_v5.rds'))
  icMat = readRDS(file = paste0(outDir, 'tradeSeqDE_evaluateK_res_v1.rds'))
  icMat = data.frame(icMat)
  colnames(icMat) = gsub('k..', '', colnames(icMat))
  
  icDev = icMat - apply(icMat, 1, mean)
  boxplot(icDev, ylim = c(-50, 50))
  
  plot(as.numeric(as.character(colnames(icMat))), apply(icMat, 2, mean, na.rm = TRUE), type = 'b')
  
  icRel =  icMat/apply(icMat, 1, max, na.rm = TRUE)
  plot(as.numeric(as.character(colnames(icMat))), apply(icRel, 2, mean, na.rm = TRUE), type = 'b')
  
  dev.off()
  
  # downsample the cells and also the genes of interest to test fitGAM
  # library(BiocParallel)
  # BPPARAM <- BiocParallel::bpparam()
  # BPPARAM # lists current options
  # BPPARAM$workers <- 16 # use 2 cores
  # 
  # 
  # set.seed(7)
  # nb_cells = 10000
  # set.seed(7)
  # index_sub = sample(c(1:ncol(counts)), nb_cells)
  # 
  # candidates = readRDS(file = paste0(outDir, 'DElist_3512genes_pairwiseComaprison.rds'))
  # genes.sel = match(candidates, rownames(counts))
  # genes.sel = genes.sel[which(!is.na(genes.sel))]
  # length(genes.sel)
  # 
  # nb.knots = 6;
  # 
  # 
  # tic()
  # set.seed(7)
  # sce <- fitGAM(counts = counts[, index_sub], 
  #               pseudotime = pseudotime[index_sub, ], 
  #               cellWeights = cellWeights[index_sub, ],
  #               genes = genes.sel,
  #               nknots = nb.knots, 
  #               verbose = TRUE, 
  #               parallel=TRUE, 
  #               BPPARAM = BPPARAM
  #               )
  # toc()
  #save(sce, file = paste0(RdataDir, 'fitGAM_output_tradeSeq_v3.Rdata'))
  
  ## reload the fitGam results
  load(paste0(outDir, 'tradeSeqDE_fitGAM_output_cellnb.10000.Rdata'))
  table(rowData(sce)$tradeSeq$converged)
  candidates = readRDS(file = paste0(outDir, 'DElist_3512genes_pairwiseComaprison.rds'))
  load(file = paste0(outDir, '/counts_pseudotime_cellWeights_for_tradeSeq.Rdata'))
  
  counts = counts[match(candidates, rownames(counts)), match(colnames(sce), colnames(counts))] 
  pseudotime = pseudotime[match(colnames(sce), colnames(counts)), ] 
  cellWeights = cellWeights[match(colnames(sce), colnames(counts)), ] 
  
  plotGeneCount(curve = pseudotime, 
                counts = counts,
                #clusters = apply(slingClusterLabels(crv), 1, which.max),
                models = sce)
  earlyDERes <- earlyDETest(sce, knots = c(1, 2))
  
  oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
  earlyDERes = earlyDERes[oEarly, ]
  head(rownames(earlyDERes))
  
  plotSmoothers(sce, counts, gene = rownames(earlyDERes)[oEarly][1])
  
  plotSmoothers(sce, counts, gene = 'Zfp703', nPoints = 1000)
  plotSmoothers(sce, counts, gene = 'Cyp26a1', nPoints = 1000)
  
  assoRes <- associationTest(sce)
  head(assoRes)
  
  startRes <- startVsEndTest(sce)
  
  oStart <- order(startRes$waldStat, decreasing = TRUE)
  sigGeneStart <- names(sce)[oStart[3]]
  plotSmoothers(sce, counts[, index_sub], gene = sigGeneStart)
  
  endRes <- diffEndTest(sce)
  o <- order(endRes$waldStat, decreasing = TRUE)
  sigGene <- names(sce)[o[1]]
  plotSmoothers(sce, counts[, index_sub], sigGene)
  
  plotGeneCount(crv, counts[, index_sub], gene = sigGene)
  
  patternRes <- patternTest(sce)
  oPat <- order(patternRes$waldStat, decreasing = TRUE)
  head(rownames(patternRes)[oPat])
  
  plotSmoothers(sce, counts[, index_sub], gene = rownames(patternRes)[oPat][4])
  
  ##########################################
  # test how to define the pseudotime for tradeseq 
  ##########################################
  # For reproducibility
  # palette(brewer.pal(8, "Dark2"))
  # data(countMatrix, package = "tradeSeq")
  # counts <- as.matrix(countMatrix)
  # rm(countMatrix)
  # data(celltype, package = "tradeSeq")
  # 
  # set.seed(22)
  # library(monocle3) # unable to install monocle3
  # 
  # # Create a cell_data_set object
  # cds <- new_cell_data_set(counts, cell_metadata = pd,
  #                          gene_metadata = data.frame(gene_short_name = rownames(counts),
  #                                                     row.names = rownames(counts)))
  # # Run PCA then UMAP on the data
  # cds <- preprocess_cds(cds, method = "PCA")
  # cds <- reduce_dimension(cds, preprocess_method = "PCA",
  #                         reduction_method = "UMAP")
  # 
  # # First display, coloring by the cell types from Paul et al
  # plot_cells(cds, label_groups_by_cluster = FALSE, cell_size = 1,
  #            color_cells_by = "cellType")
  # 
  # # Running the clustering method. This is necessary to the construct the graph
  # cds <- cluster_cells(cds, reduction_method = "UMAP")
  # # Visualize the clusters
  # plot_cells(cds, color_cells_by = "cluster", cell_size = 1)
  # 
  # # Construct the graph
  # # Note that, for the rest of the code to run, the graph should be fully connected
  # cds <- learn_graph(cds, use_partition = FALSE)
  # 
  # # We find all the cells that are close to the starting point
  # cell_ids <- colnames(cds)[pd$cellType ==  "Multipotent progenitors"]
  # closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  # closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # closest_vertex <- closest_vertex[cell_ids, ]
  # closest_vertex <- as.numeric(names(which.max(table(closest_vertex))))
  # mst <- principal_graph(cds)$UMAP
  # root_pr_nodes <- igraph::V(mst)$name[closest_vertex]
  # 
  # # We compute the trajectory
  # cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)
  # plot_cells(cds, color_cells_by = "pseudotime")
  
}


##########################################
# 
##########################################




