##########################################################################
##########################################################################
# Project:
# Script purpose: mapping transcripts to genes
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Jan 11 16:55:06 2022
##########################################################################
##########################################################################
inputDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/Transcriptomics/'
annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'

aa = read.delim(file = paste0(inputDir, 'AmexT_v47_FULL_transcripts.all.txt'), header = FALSE)
aa$V1 = gsub('>', '', aa$V1)

annot = readRDS(file = paste0(annotDir, 'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID.rds'))
aa$gene = annot$geneID[match(aa$V1, annot$transcritID)]

write.table(aa, file = paste0(inputDir, 'AmexT_v47_transcripts_genes_t2g.txt'), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')


##########################################
# prepare gene_transcript and intron sequences for kallisto 
# for scNuc sequencing and also RNA velocity
##########################################
inputDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/Transcriptomics/'
annotDir = '/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/'

## 181985 transcripts in total in current annotation
aa = read.table(file = paste0(inputDir, 'AmexT_v47_transcripts_genes_t2g.txt'), sep = '\t')

annot = readRDS(file = paste0(annotDir, 
        'AmexT_v47_transcriptID_transcriptCotig_geneSymbol.nr_geneSymbol.hs_geneID_gtf.geneInfo_gtf.transcriptInfo.rds'))

ids = unique(annot$geneID_gtf.gene)
kk = match(ids, annot$geneID_gtf.gene)

ggs = data.frame(chrom = annot$chr_gene[kk],start = annot$start_gene[kk], end = annot$end_gene[kk], 
                 names =annot$geneID_gtf.gene[kk], scores = 0, 
                 strand = annot$strand_transcript[kk], stringsAsFactors = FALSE)

ggs$names = paste0(ggs$names, '.0')
write.table(ggs, file = paste0(inputDir, 'AmexT_v47_genes.bed'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

mapping = data.frame(ggs$names, gsub('.0$', '', ggs$names))
colnames(mapping) = colnames(aa)

mapping = rbind(aa, mapping)
mapping = mapping[order(mapping$V2), ]

write.table(mapping, file = paste0(inputDir, 'AmexT_v47_unspliced.spliced.transcripts_genes_t2g.txt'), 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)


##########################################
# make intron annotation  
# test intron sequence extract code from 
# https://github.com/csoneson/rna_velocity_quant/blob/master/scripts/extractIntronSeqs.R
# it works quite well actually
##########################################
library('rtracklayer')
library(GenomicRanges)
library('GenomicFeatures')

gtf = paste0(annotDir, 'AmexT_v47.release.gtf')
annot = import(paste0(annotDir, 'AmexT_v47.release.gtf'))

txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf")

# genes = which()

grl <- GenomicFeatures::exonsBy(txdb, by = "gene")

## Collapse the exons of each gene
grl <- GenomicRanges::reduce(grl)

grl <- BiocGenerics::setdiff(range(grl), grl)

flanklength = 50
grl <- grl + flanklength

gr <- unlist(grl)

## Add -I{X} to names
xx <- gsub("\\.I\\.", ".I", make.unique(paste0(names(gr), ".I")))

names(gr) <- xx
rm(xx)
start(gr)[which(start(gr)<=1)] = 1

export(gr, con = paste0(inputDir, 'AmexT_v47_introns.bed'), format = 'bed')

mapping = data.frame(names(gr), sapply(names(gr), function(x) {unlist(strsplit(as.character(x), '[.]'))[1]}))

write.table(mapping, file = paste0(inputDir, 'AmexT_v47_introns_genes_t2g.txt'), 
            sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)



########################################################
########################################################
# Section : 
# 
########################################################
########################################################

