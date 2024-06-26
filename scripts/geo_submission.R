##########################################################################
##########################################################################
# Project: heart regeneration
# Script purpose: process metadata for GEO submission
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Nov 14 13:16:08 2023
##########################################################################
##########################################################################

metadata_Dir = '/groups/tanaka/People/current/jiwang/projects/heart_regeneration/GEO_submission_visium_mice/'

design = openxlsx::read.xlsx(paste0('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/',
                                    'R11934_visium_mice/Visium_Sample_Info.xlsx'))

design = design[, c(7,15)]
colnames(design) = c('sample', 'condition')

design$condition = gsub('Visium ', '', design$condition)
design$condition = gsub('Neonate hearts', 'neonate', design$condition)
design$condition = gsub(' D', '_d', design$condition)

raw = read.table(file = paste0(metadata_Dir, 'md5sum_raw.file'), header = FALSE)
raw = raw[grep('fastq.gz', raw$V2), ]

raw = raw[, c(2, 1)]
raw[, 1] = basename(raw[,1])

write.csv2(raw, file = paste0(metadata_Dir, 'md5sum_raw.csv'), row.names = FALSE)

processed = read.table(file = paste0(metadata_Dir, 'md5sum_processed.file'), header = FALSE)
processed = processed[, c(2, 1)]
processed[, 1] = basename(processed[,1])
write.csv2(processed, file = paste0(metadata_Dir, 'md5sum_processed.csv'), row.names = FALSE)

design = data.frame(design)

design$raw_R1 = NA
design$raw_R2 = NA
design$raw_I1 = NA
design$matrix = NA
design$barcode = NA
design$feature = NA
design$image = NA
for(n in 1:nrow(design))
{
  # n = 1
  jj1 = grep(design$sample[n], raw[,1])
  design$raw_R1[n] = raw$V2[jj1[grep('R1', raw$V2[jj1])]]
  design$raw_R2[n] = raw$V2[jj1[grep('R2', raw$V2[jj1])]]
  design$raw_I1[n] = raw$V2[jj1[grep('I1', raw$V2[jj1])]]
  jj2 = grep(design$sample[n], processed[,1])
  design$matrix[n] = processed$V2[jj2[grep('matrix', processed$V2[jj2])]]
  design$barcode[n] = processed$V2[jj2[grep('barcode', processed$V2[jj2])]]
  design$feature[n] = processed$V2[jj2[grep('feature', processed$V2[jj2])]]
  design$image[n] = processed$V2[jj2[grep('tif', processed$V2[jj2])]]
}

write.csv2(design, file = paste0(metadata_Dir, 'metadata_raw_processed.csv'), row.names = FALSE)

