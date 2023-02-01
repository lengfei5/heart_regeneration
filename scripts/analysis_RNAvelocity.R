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
outDir = paste0(resDir, '/RNA_velocity/')
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

##########################################
# test eisaR and alevin 
##########################################
use_eisaR_alevein = FALSE
if(use_eisaR_alevein){
  suppressPackageStartupMessages({
    library(Biostrings)
    library(BSgenome)
    library(eisaR)
    library(GenomicFeatures)
    library(SummarizedExperiment)
    library(tximeta)
    library(rjson)
    library(reticulate)
    library(SingleCellExperiment)
    library(scater)
  })
  
  # gtf and fa files of mouse using eisaR
  # we load the eisaR package and extract a GRanges object 
  # containing the genomic coordinates of each annotated transcript and intron. 
  # In this example, we use the ‘separate’ approach to define introns separately for each transcript, 
  # and add a flank length of 90nt to each intron.
  
  genomeDir = '/groups/tanaka/People/current/jiwang/Genomes/axolotl/Transcriptomics/alevin_velocity/'
  
  annotDir = '/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'
  gtf = paste0(annotDir, 'AmexT_v47.chr.FINAL.FULL.gtf')
  #annot = import(paste0(annotDir, 'AmexT_v47.release.gtf'))
  #txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf")
  
  #gtf <- "/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/Mus_musculus.GRCm38.87.gtf"
  genome.file = "/groups/tanaka/People/current/jiwang/Genomes/axolotl/AmexG_v6.DD.corrected.round2.chr.fa"
  
  ##########################################
  # clean a bit the gtf for transcript name 
  ##########################################
  Clean_transcript_id = FALSE
  if(Clean_transcript_id){
    library('rtracklayer')
    library(GenomicRanges)
    library('GenomicFeatures')
    
    gtf = paste0(annotDir, 'AmexT_v47.release.gtf')
    annot = rtracklayer::import(paste0(annotDir, 'AmexT_v47.release.gtf'))
    
    xx = annot[grep('chr', seqnames(annot))]
    seqlevels(xx) <- seqlevelsInUse(xx)
    
    write.table(seqlevels(xx), file = paste0(genomeDir, 'AmexT_v47_genome.chrnames.txt'), 
                sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    #sels = which(!is.na(xx$transcript_id))
    
    transcripts = xx$transcript_id
    transcripts = sapply(transcripts, function(x) {test = unlist(strsplit(as.character(x), '[|]')); 
    test[length(test)]})
    xx$transcript_id = transcripts
    
    rtracklayer::export(xx, paste0(genomeDir, 'AmexT_v47.release_chr_transcript.id.cleaned.gtf'))
  }
  
  
  gtf = paste0(genomeDir, 'AmexT_v47.release_chr_transcript.id.cleaned.gtf')
  
  grl <- eisaR::getFeatureRanges(
    gtf = gtf,
    featureType = c("spliced", "intron"), 
    intronType = "separate", 
    flankLength = 90L, 
    joinOverlappingIntrons = FALSE, 
    verbose = TRUE
  )
  
  grl[4:6]
  
  eisaR::exportToGtf(
    grl, 
    filepath = paste0(genomeDir, "AmexT_v47.annotation.expanded.gtf")
  )
  
  head(metadata(grl)$corrgene)
  
  write.table(
    metadata(grl)$corrgene, 
    file = paste0(genomeDir, "AmexT_v47.annotation.expanded.features.tsv"),
    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
  )
  
  df <- eisaR::getTx2Gene(
    grl, filepath = paste0(genomeDir, "AmexT_v47.annotation.expanded.tx2gene.tsv")
  )
  
  genome <- Biostrings::readDNAStringSet(genome.file)
  
  names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)
  seqs <- GenomicFeatures::extractTranscriptSeqs(
    x = genome, 
    transcripts = grl
  )
  
  Biostrings::writeXStringSet(
    seqs, filepath = paste0(genomeDir,  "AmexT_v47.annotation.expanded.fa")
  )
  
  
}

########################################################
########################################################
# Section : process the alevin output and prepare the data for scvelo and cellrank
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

DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)

CMmyLevels <- c(  "CM_ven_(Robo2)", "CM_ven_(Cav3_1)", "CM_IS","CM_Prol_IS", "CM_Prol_1", "CM_Prol_3")
aa$subtypes = factor(aa$subtypes, levels = CMmyLevels)

cols = c("#4CC9F0", "#49A2F0",
         "#941F56", "#C61010",  
         "#4361EE",  "#3F37C9")

DimPlot(aa, dims = c(1, 2), label = TRUE, repel = TRUE, group.by = 'subtypes', raster=FALSE,
        cols = cols
)

ggsave(filename = paste0(outDir, 'UMAP_CMsubsets_forRNAvelocity.pdf'), width = 10, height = 8)

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
# aa = readRDS(file = paste0(RdataDir, 
#                            'seuratObject_RA.symmetry.breaking_doublet.rm_mt.ribo.filtered_regressout.nCounts_',
#                            'cellCycleScoring_annot.v2_newUMAP_clusters_time_',
#                            species, version.analysis, '.rds'))

subsetting_further = FALSE
if(subsetting_further){
  
  aa = subset(aa, cells = colnames(aa)[which(aa$subtypes == 'CM_IS'|aa$subtypes == "CM_Prol_IS")])
  
  aa <- FindVariableFeatures(aa, selection.method = "vst", nfeatures = 3000)
  # 
  aa <- ScaleData(aa)
  # 
  aa <- RunPCA(aa, features = VariableFeatures(object = aa), weight.by.var = TRUE, verbose = FALSE)
  # 
  ElbowPlot(aa, ndims = 50)
  
  aa <- RunUMAP(aa, dims = 1:20, n.neighbors = 20, min.dist = 0.3)
  DimPlot(aa, label = TRUE, repel = TRUE, group.by = 'condition', raster=FALSE)
  
  source("functions_scRNAseq.R")
  explore.umap.params.combination(sub.obj = aa, resDir = outDir, 
                                  pdfname = 'axolotl_snRNA_RNAvelocity_umap_test.pdf',
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
  
  
}


## load saved spliced and unspliced 
load(file = paste0(RdataDir, 'seuratObject_spliced_unspliced_', species, version.analysis, '.Rdata'))
aa$cell.id = aa$cell.ids

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
mnt = subset(mnt, features = rownames(counts))

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
#saveFile = "RNAmatrix_umap_kalisto.velocity_spliced_unspliced_CMsutypes_v1.2.h5Seurat"
saveFile = 'RNAmatrix_umap_kalisto.velocity_spliced_unspliced_CMsutypes_injurySpec_v2.2.h5Seurat'

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


