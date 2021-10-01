##########################################################################
##########################################################################
# Project: Elad's heart regeneration 
# Script purpose: the common functions used for adult, neonadal mice and axolotl
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Oct  1 10:37:45 2021
##########################################################################
##########################################################################
require(Seurat)
require(SeuratObject)
require(ggplot2)
require(tibble)
require(dplyr)
library(patchwork)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

