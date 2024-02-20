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
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
library(circlize)
library(RColorBrewer)
require(scran)
require(scater)
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


