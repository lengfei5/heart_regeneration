##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: utility functions for RNA velocity prepartion
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Apr 27 11:56:00 2023
##########################################################################
##########################################################################
process_spliced_unspliced_kallisto = function(aa)
{
  #a$time = gsub('Amex_', '', aa$condition)
  #aa$cell.ids = sapply(colnames(aa), function(x) unlist(strsplit(as.character(x), '-'))[1]) 
  #aa$cell.ids = paste0(aa$cell.ids, '_', aa$time)
  
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
      cat(nrow(srat), ' cells in ', design$condition[n], '\n')
      
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
       file = paste0('../results/Rdata/', 
                     'seuratObject_spliced_unspliced_celltypes.all.Rdata'))
  
  #rm(srat)
  #rm(count.data)
  #rm(meta)
  #rm(g)
  
}


preapre_dataFile_for_RNAvelocity_PAGA = function(seuratObj)
{
  # seuratObj = sub_obj
  seuratObj$cell.id = seuratObj$cell.ids
  seuratObj$celltypes = seuratObj$subtypes
  DefaultAssay(seuratObj) = 'RNA'
  
  cat('-- Drop other assays besides RNA -- \n')
  seuratObj = DietSeurat(seuratObj, counts = TRUE, data = TRUE,
                   scale.data = TRUE,
                   features = rownames(seuratObj), 
                   assays = c('RNA'), 
                   dimreducs = c('umap'), 
                   graphs = NULL, 
                   misc = TRUE
  )
  
  #p1 = DimPlot(seuratObj, label = TRUE, repel = TRUE, group.by = 'celltypes', raster=FALSE) 
  #plot(p1)
  
  # import the processed spliced and unspliced seurat object
  load(file = paste0('../results/Rdata/', 
                     'seuratObject_spliced_unspliced_celltypes.all.Rdata'))
  
  cat('-- select features shares by spliced and unspliced -- \n') 
  features = intersect(rownames(spliced), rownames(unspliced))
  spliced = subset(spliced, features = features)
  unspliced = subset(unspliced, features = features)
  
  cat('-- subsetting seuratObj with the same cell as spliced and unspliced -- \n') 
  cells.shared = intersect(spliced$cell.id, seuratObj$cell.id)
  mnt = subset(seuratObj, cells = colnames(seuratObj)[match(cells.shared, seuratObj$cell.id)])
  
  counts = GetAssayData(spliced, slot = 'counts')
  counts = counts[, match(mnt$cell.id, spliced$cell.id)]
  colnames(counts) = colnames(mnt)
  
  cat('-- subsetting mnt with the same genes as spliced and unspliced -- \n') 
  DefaultAssay(mnt) = 'RNA'
  mnt = subset(mnt, features = intersect(rownames(counts), rownames(mnt)))
  
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
  
  cat('-- reduce the integrated seurat object -- \n') 
  DefaultAssay(mnt) = 'RNA'
  mnt = DietSeurat(mnt, counts = TRUE, data = TRUE,
                   scale.data = FALSE,
                   features = rownames(mnt), 
                   assays = c('RNA', 'spliced', 'unspliced'), 
                   dimreducs = c('umap'), graphs = NULL, 
                   misc = TRUE
  )
  
  
  DefaultAssay(mnt) = 'RNA'
  #VariableFeatures(mnt)
  
  Idents(mnt) = mnt$condition
  mnt$condition = as.character(mnt$condition)
  mnt$celltypes = as.character(mnt$subtypes)
  
  return(mnt)
  
  
}