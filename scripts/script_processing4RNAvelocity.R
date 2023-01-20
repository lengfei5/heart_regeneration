##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: prepare files and matrix for scVelo
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Jan  9 14:15:39 2023
##########################################################################
##########################################################################

########################################################
########################################################
# Section : processing for spliced and unspliced matrix quantification using kallisto
# the same code as RA competence
########################################################
########################################################
rm(list = ls())

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
  library(BUSpaRse)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(AnnotationHub)
})

#devtools::install_github("BUStools/BUSpaRse")
#require(Biostrings)

gtf_file = '/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/AmexT_v47.FULL.gtf'
transcriptome_amex6 <- paste0("/groups/tanaka/People/current/jiwang/Genomes/axolotl/Transcriptomics/Kallisto_Index", 
                              "/AmexT_v47.FULL.fa")
genome_file = paste0("/groups/tanaka/People/current/jiwang/Genomes/axolotl", 
                     "/AmexG_v6.DD.corrected.round2.chr.fa")

outDir = paste0("/groups/tanaka/People/current/jiwang/Genomes/axolotl/Transcriptomics/",
                  "kallisto_RNAvelocity_index/test")
system(paste0('mkdir -p ', outDir))

## test BUSpaRse but did not work; original code from https://bustools.github.io/BUS_notebooks_R/velocity.html
Use.BUSpaRse.to.prepare.annotation = FALSE
if(Use.BUSpaRse.to.prepare.annotation){
  ah <- AnnotationHub()
  query(ah, pattern = c("Ensembl", "97", "Mus musculus", "EnsDb"))
  edb <- ah[["AH73905"]]
  require("ensembldb")
  
  get_velocity_files(edb, L = 91, Genome = BSgenome.Mmusculus.UCSC.mm10, 
                     out_path = outDir, 
                     isoform_action = "separate")
  
  
  get_velocity_files(gtf_file,
                     L = 51,
                     Genome = genome_amex6,
                     Transcriptome = transcriptome_amex6,
                     gene_version = 'amex6',
                     transcript_version = 'v47',
                     out_path = paste0("/groups/tanaka/People/current/jiwang/Genomes/axolotl/Transcriptomics/",
                                       "kallisto_RNAvelocity_index"),
                     isoform_action = "separate")
  
  
}

##########################################
# preapre the intron and transcript reference fasta files
# this part actually has been done previously by 
# https://github.com/lengfei5/axolotl_annotation/blob/main/scripts/AmexT_transcript2gene.R
##########################################
# gtf and fa files of mouse using eisaR
# we load the eisaR package and extract a GRanges object 
# containing the genomic coordinates of each annotated transcript and intron. 
# In this example, we use the ‘separate’ approach to define introns separately for each transcript, 
# and add a flank length of 90nt to each intron.
#genomeDir = '/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/'
#gtf <- "/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/Mus_musculus.GRCm38.87.gtf"
#genome.file = "/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/Mus_musculus.GRCm38.dna.toplevel.fa"
Test_eisaR = FALSE
if(Test_eisaR){
  grl <- eisaR::getFeatureRanges(
    gtf = gtf_file,
    featureType = c("spliced", "intron"), 
    intronType = "separate", 
    flankLength = 90L, 
    joinOverlappingIntrons = FALSE, 
    verbose = TRUE
  )
  
  grl[4:6]
  
  #genome <- Biostrings::readDNAStringSet(genome.file)
  genome_amex6 <- Biostrings::readDNAStringSet(genome_file)
  
  names(genome_amex6) <- sapply(strsplit(names(genome_amex6), " "), .subset, 1)
  seqs <- GenomicFeatures::extractTranscriptSeqs(
    x = genome_amex6, 
    transcripts = grl
  )
  
  Biostrings::writeXStringSet(
    seqs, filepath = paste0(genomeDir,  "Mus.GRCm38.87.annotation.expanded.fa")
  )
  
  eisaR::exportToGtf(
    grl, 
    filepath = paste0(outDir, "/Amex47.annotation.expanded.gtf")
  )
  
  head(metadata(grl)$corrgene)
  
  write.table(
    metadata(grl)$corrgene, 
    file = paste0(outDir, "/Amex47.annotation.expanded.features.tsv"),
    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
  )
  
  df <- eisaR::getTx2Gene(
    grl, filepath = paste0(outDir, "/Amex47.annotation.expanded.tx2gene.tsv")
  )
  
}
