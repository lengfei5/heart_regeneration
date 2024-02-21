##########################################################################
##########################################################################
# Project: main script for cross species analysis 
# Script purpose: 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Feb 20 15:26:14 2024
##########################################################################
##########################################################################
rm(list = ls())

library(pryr) # monitor the memory usage
require(ggplot2)
#library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
library(circlize)
library(RColorBrewer)
require(scran)
require(scater)
library(liana)
source('functions_scRNAseq.R')
source('functions_Visium.R')
source('functions_cccInference.R')

version.analysis = '_20240220'

resDir = paste0("../results/cross_species", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

mem_used()

########################################################
########################################################
# Section I: cross-species co-embedding of scRNA-seq
# 
########################################################
########################################################

########################################################
########################################################
# Section II: cross-species cell type correlation and differential analysis
# 
########################################################
########################################################

########################################################
########################################################
# Section III : Differential LR analysis in ST 
# 
########################################################
########################################################

outDir = paste0(resDir, '/LR_differentialAnalysis/')
if(!dir.exists(outDir)) dir.create(outDir)

##########################################
# test the LIANA + tensor_cell2cell
# original code : https://saezlab.github.io/liana/articles/liana_cc2tensor.html
##########################################
library(tidyverse, quietly = TRUE)
library(SingleCellExperiment, quietly = TRUE)
library(reticulate, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(liana, quietly = TRUE)
library(ExperimentHub, quietly = TRUE)

eh <- ExperimentHub()
# Get Data
(sce <- eh[["EH2259"]])


# basic feature filtering
sce <- sce[rowSums(counts(sce) >= 1) >= 5, ]

# basic outlier filtering
qc <- scater::perCellQCMetrics(sce)

# remove cells with few or many detected genes
ol <- scater::isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol]

# Remove doublets
sce <- sce[, sce$multiplets=="singlet"]

# Set rownames to symbols
rownames(sce) <- rowData(sce)$SYMBOL

# log-transform
sce <- scuttle::logNormCounts(sce)

# Create a label unique for every sample
sce$context <- paste(sce$stim, sce$ind, sep="|")

# Plot
sce %>%
  get_abundance_summary(sample_col = "context",
                        idents_col = "cell", 
                        min_cells = 10, # min cells per sample
                        min_samples = 3, # min samples
                        min_prop = 0.2 # min prop of samples
  ) %>%
  plot_abundance_summary()

# filter non abundant celltypes
sce <- liana::filter_nonabundant_celltypes(sce,
                                           sample_col = "context",
                                           idents_col = "cell")

# Run LIANA by sample
sce <- liana_bysample(sce = sce,
                      sample_col = "context",
                      idents_col = "cell",
                      method = "sca", # we use SingleCellSignalR's score alone
                      expr_prop = 0, # expression proportion threshold
                      inplace=TRUE, # saves inplace to sce
                      return_all = FALSE # whether to return non-expressed interactions 
)
