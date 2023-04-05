##########################################################################
##########################################################################
# Project: 
# Script purpose: double check axolotl annotation for CellOracle
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Mar  3 14:31:08 2023
##########################################################################
##########################################################################
rm(list = ls())

annotDir = '/groups/tanaka/People/current/jiwang/Genomes/axolotl/CellOracle_annot/'
resDir = annotDir

##########################################
# collect gene id and gene symobls 
##########################################
aa = rtracklayer::import(paste0(annotDir, 'AmexT_v47.release.gtf'))
xx = aa
xx = as.data.frame(xx)

aa = aa[which(aa$type == 'gene')]
aa = as.data.frame(aa)
aa = aa[, c(1:5, 10:11)]

# import manually collected gene symbols
annot = readRDS(paste0(annotDir,
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr',
                       '_curated.geneSymbol.toUse.rds'))
annot = annot[which(!is.na(annot$gene.symbol.toUse)), ]

mm = match(aa$gene_id, annot$geneID)

aa$gene_name[which(is.na(mm))] = NA
aa$gene_name[which(!is.na(mm))] = annot$gene.symbol.toUse[mm[which(!is.na(mm))]]

kk = c()
for(n in 1:nrow(aa))
{
  #n = 1
  strands = as.character(xx$strand[which(xx$gene_id == aa$gene_id[n])])
  if(length(unique(strands)) > 1){
    cat(n, '--', aa$gene_id[n], '\n')
    
  }
  
}

#aa$gene_name = gsub("\\s*\\([^\\)]+\\)", "", aa$gene_name[2])
write.table(aa, file = paste0(resDir, '/amexT_v47_4CellOracle.txt'), 
            sep = '\t', col.names = TRUE, row.names = FALSE, 
            quote = FALSE)

##########################################
# double check the annotation for CellOracle
##########################################
annotDir = '/groups/tanaka/People/current/jiwang/Genomes/axolotl/CellOracle_annot/'

library('rtracklayer')
library(GenomicRanges)
library('GenomicFeatures')

gtf = paste0(annotDir, 'AmexT_v47.release.gtf')
annot = import(paste0(annotDir, 'AmexT_v47.release.gtf'))

genes = read.table(file = paste0(annotDir, 'amexT_v47_hoxPatch.txt'), sep = '\t', header = TRUE)

mm = match(genes$gene_id, annot$gene_id)

kk = match(annot$gene_id, genes$gene_id)
jj = which(is.na(kk))

jj0 = jj[grep('^chr', seqnames(annot)[jj])]


kk = which(!is.na(kk))

geneids = annot$gene_id
head(grep('AMEX60DD20', geneids))

