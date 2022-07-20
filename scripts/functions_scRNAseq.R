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
  return(sapply(xx, function(x) unlist(strsplit(as.character(x), '_'))[1]))
}

get_geneID = function(xx)
{
  return(sapply(xx, function(x) {test = unlist(strsplit(as.character(x), '_')); return(test[length(test)])}))
  
}

########################################################
########################################################
# Section : import kalisto output to make 10x seurat object
# orignal code from Tomas' code
########################################################
########################################################
make_SeuratObj_scRNAseq(topdir = './', saveDir = './results', changeGeneName = TRUE, keyname = 'Amex_scRNA_d0',
                        QC.umi = TRUE)
{
  library(Seurat)
  library(DropletUtils)
  library(edgeR)
  
  if(!dir.exists(saveDir)) dir.create(saveDir)
  
  # read in data
  # topdir = "${outbus}" # source dir
  exp = Matrix::readMM(paste0(topdir, "genecounts.mtx")) #read matrix
  bc = read.csv(paste0(topdir, "/genecounts.barcodes.txt"), header = F, stringsAsFactors = F)
  g = read.csv(paste0(topdir, "/genecounts.genes.txt"), header = F, stringsAsFactors = F)
  
  if(changeGeneName){
    cat('change gene names \n')
    # change the gene names before making Seurat object
    annot = readRDS(paste0('/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                           'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse_manual_v1.rds'))
    
    
    mm = match(g$V1, annot$geneID)
    ggs = paste0(annot$gene.symbol.toUse[mm], '_',  annot$geneID[mm])
    g$V1[!is.na(mm)] = ggs[!is.na(mm)]
    
  }
  
  dimnames(exp) = list(paste0(bc$V1,"-1"),g$V1) # number added because of seurat format for barcodes
  count.data = Matrix::t(exp)
  
  #dimnames(exp) = list(paste0(bc\$V1,"-1"),g\$V1) # number added because of seurat format for barcodes
  #count.data = Matrix::t(exp)
  if(QC.umi){
    # plot rankings for number of UMI
    cat('plot UMI rank \n')
    
    # get emptyDrops and default cutoff cell estimates
    iscell_dd = defaultDrops(count.data, expected = 5000)
    eout = emptyDrops(count.data, lower = 100)
    eout$FDR[is.na(eout$FDR)] = 1
    
    iscell_ed = eout$FDR<=0.01
    meta = data.frame(row.names = paste0(bc$V1,"-1"),
                      iscell_dd = iscell_dd, iscell_ed = iscell_ed)
    
    # plot rankings for number of UMI
    br.out <- barcodeRanks(count.data)
    
    pdf(paste0(saveDir, "UMIrank.pdf"), height = 8, width = 8, useDingbats = FALSE)
    
    plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
    
    o <- order(br.out$rank)
    lines(br.out$rank[o], br.out$fitted[o], col="red")
    abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
    abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
    legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"),
           legend=c("knee", "inflection"))
    
    dev.off()
    
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
  srat = CreateSeuratObject(counts = count.data[,iscell_dd | iscell_ed],
                            meta.data = meta[iscell_dd | iscell_ed,])
  
  amb_prop = estimateAmbience(count.data)[rownames(srat@assays$RNA@meta.features)]
  srat@assays$RNA@meta.features = data.frame(row.names = rownames(srat@assays$RNA@meta.features),
                                              "ambient_prop" = amb_prop)
  
  # get MT% (genes curated from NCBI chrMT genes)
  mtgenes = c("COX1", "COX2", "COX3", "ATP6", "ND1", "ND5", "CYTB", "ND2", "ND4",
             "ATP8", "MT-CO1", "COI", "LOC9829747")
  mtgenes = c(mtgenes, paste0("MT", mtgenes), paste0("MT-", mtgenes))
  
  mt_features = rownames(srat)[get_geneName(g[,1]) %in% mtgenes]
  srat = PercentageFeatureSet(srat, col.name = "percent.mt", assay = "RNA",
                             features = mt_features)
  
  saveRDS(srat, file = paste0(saveDir, "srat.RDS"))
  
  return(srat)
  
}
