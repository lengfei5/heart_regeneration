##########################################################################
##########################################################################
# Project:
# Script purpose: functions of cell-cell communication inference
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Oct 17 09:41:55 2022
##########################################################################
##########################################################################
########################################################
########################################################
# Section : ligand-receptor-TFs analysis
# 
########################################################
########################################################
run_LIANA_defined_celltype = function(subref,
                                      celltypes,
                                      receivers = NULL,
                                      additionalLabel = '_fixedCelltypes')
{
  source('functions_scRNAseq.R')
  
  # subref = refs
  sce <- as.SingleCellExperiment(subref)
  colLabels(sce) = as.factor(sce$celltypes)
  rownames(sce) = toupper(get_geneName(rownames(sce)))
  
  ave.counts <- calculateAverage(sce, assay.type = "counts")
  
  #hist(log10(ave.counts), breaks=100, main="", col="grey80",
  #     xlab=expression(Log[10]~"average count"))
  
  num.cells <- nexprs(sce, byrow=TRUE)
  #smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
  #              xlab=expression(Log[10]~"average count"))
  
  # detected in >= 5 cells, ave.counts >=5 but not too high
  genes.to.keep <- num.cells > 20 & ave.counts >= 10^-4  & ave.counts <10^2  
  summary(genes.to.keep)
  
  sce <- sce[genes.to.keep, ]
  
  
  ## run the liana wrap function by specifying resource and methods
  # Resource currently included in OmniPathR (and hence `liana`) include:
  show_resources()
  # Resource currently included in OmniPathR (and hence `liana`) include:
  show_methods()
  
  liana_test <- liana_wrap(sce,  
                           method = c("natmi", "connectome", "logfc", "sca", "cytotalk"  
                                      #'cellphonedb'
                                      ),
                           resource = c("Consensus", 'CellPhoneDB', "OmniPath", "LRdb",
                                        "CellChatDB",  "CellTalkDB"), 
                           assay.type = "logcounts", 
                           idents_col = 'celltypes')
  
  # Liana returns a list of results, each element of which corresponds to a method
  # liana_test %>% glimpse
  
  # We can aggregate these results into a tibble with consensus ranks
  liana_test <- liana_test %>%
    liana_aggregate(resource = 'Consensus')
  
  saveRDS(liana_test, file = paste0(outDir, '/res_lianaTest_Consensus', additionalLabel, '.rds'))
  
  # liana_test = readRDS(file = paste0(outDir, '/res_lianaTest_Consensus_day1.rds'))
  if(is.na(receivers)){ # loop over all cell type candidates
    receivers = celltypes
  }
  
  for(m in 1:length(receivers)){
    # m = 1
    #liana_test %>%
    #  liana_dotplot(source_groups = celltypes[n],
    #                target_groups = celltypes,
    #                ntop = ntop)
    for(ntop in c(100, 200, 300, 500)){
      liana_test %>%
        liana_dotplot(source_groups = celltypes,
                      target_groups = receivers[m],
                      ntop = ntop)
      #liana_test_save =  liana_test %>% filter()
      #  liana_dotplot(source_groups = celltypes,
      #                target_groups = receivers[m],
      #                ntop = ntop)
      ggsave(filename = paste0(outDir, '/liana_LR_prediction_recieveCell', 
                               additionalLabel, 
                               '_receiverCells.', receivers[m], 
                               '_ntop.', ntop, '.pdf'), 
             width = 30, height = min(c(10*ntop/20, 50)), limitsize = FALSE)
    }
   
  }
  
  pdfname = paste0(outDir, '/liana_celltype_communication_freqHeatmap', additionalLabel, '.pdf')
  pdf(pdfname, width=20, height = 8)
  
  liana_trunc <- liana_test %>%
    # only keep interactions concordant between methods
    filter(aggregate_rank <= 0.01) # this can be FDR-corr if n is too high
  
  heat_freq(liana_trunc)
  
  dev.off()
  
}

##########################################
# main function of LIANA
##########################################
# original code from https://saezlab.github.io/liana/articles/liana_tutorial.html
run_LIANA = function(refs,
                     celltypes = NULL,
                     celltypes_timeSpecific = NULL,
                     timepoint_specific = FALSE,
                     receiver_cells = list(c('CM_IS')),
                     outDir = '../results/Ligand_Receptor_analysis',
                     ntop = 50,
                     RUN.CPDB.alone = FALSE)
{
  require(tidyverse)
  require(magrittr)
  require(liana)
  require(scran)
  require(scater)
  require(DelayedMatrixStats)
  require(genefilter)
  require(matrixStats)
  require(MatrixGenerics)
  require(sparseMatrixStats)
  
  ## double check timepoint_specific and celltype 
  if(!timepoint_specific){
    cat(' -- no timepoint to consider -- \n')
    if(is.null(celltypes)) {
      stop('no cell types selected found \n')
    }else{
      if(is.vector(celltypes) != TRUE | length(celltypes) ==0){
        stop('a vector of cell types expected \n')
      }else{
        mm = match(celltypes, refs$celltypes)
        if(length(which(is.na(mm)))>0){
          stop('some selected cell types not found in refs -- ', celltypes[which(is.na(mm))])
        }else{
          cat(celltypes, '\n')
          cat('all selected cell types were found in the refs \n')
        }
      }
    }
  }else{
    cat('-- time specific to consider -- \n')
    if(is.null(celltypes_timeSpecific)){
      stop('no list of time-specific cell types found \n')
    }else{
      if(!is.list(celltypes_timeSpecific) | length(celltypes_timeSpecific) == 0){
        stop('a list of time-specific cell types expected \n')
      }else{
        cat(length(celltypes_timeSpecific), 'time points specified \n')
        cat(names(celltypes_timeSpecific), '\n')
        
        for(n in 1:length(celltypes_timeSpecific))
        {
          celltypes = celltypes_timeSpecific[[n]]
          mm = match(celltypes, refs$celltypes)
          if(length(which(is.na(mm)))>0){
            stop(names(celltypes_timeSpecific)[n], 
                 ': some selected cell types not found in refs -- ', celltypes[which(is.na(mm))])
          }else{
            cat(names(celltypes_timeSpecific)[n], celltypes,'\n')
            cat('--all selected cell types were found in the refs \n')
          }
        }
      }
    }
    
  }
  
  # prepre for the LIANA loop
  system(paste0('mkdir -p ', outDir))
  Idents(refs) = as.factor(refs$celltypes)
  
  if(!timepoint_specific){
    subref = subset(refs, cells = colnames(refs)[!is.na(match(refs$celltypes, celltypes))])
    subref$celltypes = droplevels(as.factor(subref$celltypes))
    table(subref$celltypes)
    #subref = subset(x = subref, downsample = 1000)
    cat('celltype to consider -- ', names(table(subref$celltypes)), '\n')
    
    Idents(subref) = subref$celltypes
    
    # Run liana
    # liana_test <- liana_wrap(testdata, method = 'cellphonedb', resource = 'CellPhoneDB')
    run_LIANA_defined_celltype(sburef,
                               celltypes = celltypes,
                               receiver_cells = receiver_cells,
                               additionalLabel = '_fixedCelltypes')
    
  }else{
    for(n in 1:length(celltypes_timeSpecific))
    {
      # n = 1
      celltypes = celltypes_timeSpecific[[n]]
      timepoint = names(celltypes_timeSpecific)[n]
      receivers = receiver_cells[[n]]
      
      cat(n, '-- time point : ', timepoint, '\n')
      cat('---- celltype : ', celltypes, '\n')
      cat('---- receiver : ', receivers, "\n")
      
      celltypes_sel = unique(c(celltypes, receivers))
      
      subref = subset(refs, cells = colnames(refs)[!is.na(match(refs$celltypes, celltypes_sel))])
      subref$celltypes = droplevels(as.factor(subref$celltypes))
      table(subref$celltypes)
      #subref = subset(x = subref, downsample = 1000)
      
      cat('celltype to consider -- ', names(table(subref$celltypes)), '\n')
      
      Idents(subref) = subref$celltypes
      
      # Run liana
      run_LIANA_defined_celltype(subref,
                                 celltypes = celltypes,
                                 receivers = receivers,
                                 additionalLabel = paste0('_', timepoint))

    }
    
      
  }
  
  # RUN CPDB alone
  if(RUN.CPDB.alone){
    cpdb_test <- liana_wrap(sce,
                            method = 'cellphonedb',
                            resource = c('CellPhoneDB'),
                            permutation.params = list(nperms=1000,
                                                      parallelize=FALSE,
                                                      workers=4), 
                            expr_prop=0.05)
    
    # Plot toy results
    # identify interactions of interest
    cpdb_int <- cpdb_test %>%
      # only keep interactions with p-val <= 0.05
      filter(pvalue <= 0.05) %>% # this reflects interactions `specificity`
      # then rank according to `magnitude` (lr_mean in this case)
      rank_method(method_name = "cellphonedb",
                  mode = "magnitude") %>%
      # keep top 20 interactions (regardless of cell type)
      distinct_at(c("ligand.complex", "receptor.complex")) %>%
      head(ntop)
    
    # Plot toy results
    cpdp_res = cpdb_test %>%
      # keep only the interactions of interest
      inner_join(cpdb_int, 
                 by = c("ligand.complex", "receptor.complex")) %>%
      # invert size (low p-value/high specificity = larger dot size)
      # + add a small value to avoid Infinity for 0s
      mutate(pvalue = -log10(pvalue + 1e-10)) 
    
    for(n in 1:length(celltypes)){
      # n = 1
      cpdp_res %>% 
        liana_dotplot(source_groups = celltypes[n],
                      target_groups = celltypes,
                      specificity = "pvalue",
                      magnitude = "lr.mean",
                      show_complex = TRUE,
                      size.label = "-log10(p-value)", 
                      ntop = ntop)
      
      ggsave(filename = paste0(outDir, '/CellPhoneDB_LR_prediction_senderCell_', celltypes[n], '.pdf'), 
             width = 14, height = min(c(10*ntop/20, 50)), limitsize = FALSE)
      
    }
  }
  
}

assembly.liana.plot = function()
{
  liana_test = readRDS(file = paste0(outDir, '/res_lianaTest_Consensus.rds'))
  ntop = 80
  
  for(n in 1:length(celltypes)){
    # n = 1
    #liana_test %>%
    #  liana_dotplot(source_groups = celltypes[n],
    #                target_groups = celltypes,
    #                ntop = ntop)
    liana_test %>%
      liana_dotplot(source_groups = celltypes,
                    target_groups = celltypes[n],
                    ntop = ntop)
    
    ggsave(filename = paste0(outDir, '/liana_LR_prediction_recieveCell_', celltypes[n], 
                             '_ntop', ntop, '.pdf'), 
           width = 16, height = 6*ntop/20, limitsize = FALSE)
    
  }
  
}


########################################################
########################################################
# Section II : Differential NicheNet analysis
# 
########################################################
########################################################
run_Diff_NicheNet = function(refs = refs, 
                             timepoint_specific = TRUE,
                             celltypes = NULL,
                             celltypes_BZ_timeSpecific = NULL, 
                             celltypes_RZ_timeSpecific = NULL,
                             receivers_BZ_timeSpecific = list(c('CM_IS')),
                             receivers_RZ_timeSpecific = list(c('CM_ROBO2')), 
                             outDir = outDir,
                             ntop = c(50, 100, 150, 200))
{
  library(pryr) # monitor the memory usage
  require(ggplot2)
  library(nichenetr)
  library(Seurat) # please update to Seurat V4
  library(tidyverse)
  library(circlize)
  library(RColorBrewer)
  require(scran)
  require(scater)
  library(ggrepel)
  library(cowplot)
  source('functions_scRNAseq.R')
  source('functions_Visium.R')
  
  ## double check timepoint_specific and celltype 
  if(!timepoint_specific){
    cat(' -- no timepoint to consider -- \n')
    if(is.null(celltypes)) {
      stop('no cell types selected found \n')
    }else{
      if(is.vector(celltypes) != TRUE | length(celltypes) ==0){
        stop('a vector of cell types expected \n')
      }else{
        mm = match(celltypes, refs$celltypes)
        if(length(which(is.na(mm)))>0){
          stop('some selected cell types not found in refs -- ', celltypes[which(is.na(mm))])
        }else{
          cat(celltypes, '\n')
          cat('all selected cell types were found in the refs \n')
        }
      }
    }
  }else{
    cat('-- time specific to consider -- \n\n')
    
    if(is.null(celltypes_BZ_timeSpecific)| is.null(celltypes_RZ_timeSpecific)){
      stop('no list of time-specific cell types BZ or RZ found \n')
    }else{
      if(!is.list(celltypes_BZ_timeSpecific) | length(celltypes_BZ_timeSpecific) == 0 |
         !is.list(celltypes_RZ_timeSpecific) | length(celltypes_RZ_timeSpecific) == 0){
        stop('lists of time-specific BZ/RZ cell types expected \n')
        
      }else{
        cat(length(celltypes_BZ_timeSpecific), 'time points specified for BZ  \n')
        cat(names(celltypes_BZ_timeSpecific), '\n\n')
        
        cat(length(celltypes_RZ_timeSpecific), 'time points specified for RZ \n')
        cat(names(celltypes_RZ_timeSpecific), '\n\n')
        
        if(length(celltypes_BZ_timeSpecific) != length(celltypes_RZ_timeSpecific)){
          stop('not the same time points specified for BZ and RZ')
        } 
        
        for(n in 1:length(celltypes_BZ_timeSpecific))
        {
          celltypes = celltypes_BZ_timeSpecific[[n]]
          mm = match(celltypes, refs$celltypes)
          if(length(which(is.na(mm)))>0){
            stop(names(celltypes_BZ_timeSpecific)[n], 
                 ': some selected cell types in BZ not found in refs  -- ', celltypes[which(is.na(mm))])
          }else{
            cat(names(celltypes_BZ_timeSpecific)[n], celltypes,'\n')
            cat('--all selected cell types of BZ were found in the refs \n')
          }
        }
        
        for(n in 1:length(celltypes_RZ_timeSpecific))
        {
          celltypes = celltypes_RZ_timeSpecific[[n]]
          mm = match(celltypes, refs$celltypes)
          if(length(which(is.na(mm)))>0){
            stop(names(celltypes_RZ_timeSpecific)[n], 
                 ': some selected cell types in RZ not found in refs  -- ', celltypes[which(is.na(mm))])
          }else{
            cat(names(celltypes_RZ_timeSpecific)[n], celltypes,'\n')
            cat('--all selected cell types of RZ were found in the refs \n')
          }
        }
        
        for(n in length(receivers_BZ_timeSpecific))
        {
          # n = 1
          receivers = receivers_BZ_timeSpecific[[n]]
          receivers.ctl =  receivers_RZ_timeSpecific[[n]]
          
          mm = match(receivers, celltypes_BZ_timeSpecific[[n]])
          if(any(is.na(mm))){stop(receivers, ' not found in the celltypes_BZ') }
          
          if(length(receivers.ctl)!= 1) stop('only one control receiver is allowed at each time point')
          mm.ctl = match(receivers.ctl, celltypes_RZ_timeSpecific[[n]])
          if(is.na(mm.ctl)){stop(receivers.ctl, ' not found in the celltype_RZ')}
          
        }
      }
    }
  }
  
  
  # prepre for the LIANA loop
  system(paste0('mkdir -p ', outDir))
  Idents(refs) = as.factor(refs$celltypes)
  
  ##########################################
  # step 0):  load the Nichenet data 
  ##########################################
  # NicheNet’s ligand-target prior model
  ligand_target_matrix = readRDS(paste0(dataPath_nichenet,  "ligand_target_matrix.rds"))
  ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
  
  # ligand-receptor network, Putative ligand-receptor links from NicheNet
  lr_network = readRDS(paste0(dataPath_nichenet, "lr_network.rds"))
  lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% 
    distinct(ligand, receptor, bonafide)
  
  # If wanted, users can remove ligand-receptor interactions that were predicted based on 
  # protein-protein interactions and 
  # only keep ligand-receptor interactions that are described in curated databases.
  # lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  head(lr_network)
  
  # ## weighted integrated networks 
  weighted_networks = readRDS(paste0(dataPath_nichenet,  "weighted_networks.rds"))
  head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
  head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
  
  weighted_networks_lr = weighted_networks$lr_sig %>% 
    dplyr::rename(ligand = from, receptor = to) %>%
    inner_join(lr_network %>% 
                 distinct(ligand,receptor), by = c("ligand","receptor"))
  
  
  ##########################################
  # loop over the time points
  ##########################################
  for(n in 1:length(celltypes_BZ_timeSpecific))
  {
    # n = 2
    timepoint = names(celltypes_BZ_timeSpecific)[n]
    celltypes_BZ = celltypes_BZ_timeSpecific[[n]]
    celltypes_RZ = celltypes_RZ_timeSpecific[[n]]
    
    receivers_BZ = receivers_BZ_timeSpecific[[n]]
    receivers_RZ = receivers_RZ_timeSpecific[[n]]
    
    cat(n, '-- ', timepoint, '\n')
    cat('BZ cell types : \n ', celltypes_BZ, '\n')
    cat('RZ cell types : \n ', celltypes_RZ, '\n')
    
    cat(' -- receiver cell of BZ niche :  ', receivers_BZ, '\n')
    cat(' -- receiver cell of RZ niche :  ', receivers_RZ, '\n')
    
    ##########################################
    # from here the arguments are 
    # - BZ cell types
    # - RZ cell types
    # - receiver cell of BZ niche
    # - receiver cell of RZ niche
    ##########################################
    celltypes_sel = unique(c(celltypes_BZ, receivers_BZ, 
                             celltypes_RZ, receivers_RZ)
    )
    Idents(refs) = as.factor(refs$celltypes)
    subref = subset(refs, cells = colnames(refs)[!is.na(match(refs$celltypes, celltypes_sel))])
    subref$celltypes = droplevels(as.factor(subref$celltypes))
    table(subref$celltypes)
    
    #subref$celltypes = as.character(subref$celltypes)
    cat('celltype to consider -- ', names(table(subref$celltypes)), '\n')
    #Idents(subref) = as.factor(subref$celltypes)
    table(subref$celltypes)
    
    subref = subset(x = subref, downsample = 3000) # downsample the CM and EC for the sake of speed
    table(subref$celltypes)
    Idents(subref) = subref$celltypes
    
    ##########################################
    # Gene filtering, in particular filtering gene without UMI counts or lowly expressed genes
    ##########################################
    Gene.filtering.preprocessing = TRUE
    if(Gene.filtering.preprocessing){
      sce <- as.SingleCellExperiment(subref)
      ggs = rownames(sce)
      
      colLabels(sce) = as.factor(sce$celltypes)
      rownames(sce) = get_geneName(rownames(sce))
      
      ave.counts <- calculateAverage(sce, assay.type = "counts")
      
      #hist(log10(ave.counts), breaks=100, main="", col="grey80",
      #     xlab=expression(Log[10]~"average count"))
      
      num.cells <- nexprs(sce, byrow=TRUE)
      #smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
      #              xlab=expression(Log[10]~"average count"))
      
      # detected in >= 10 cells, ave.counts cutoff are both very lose
      genes.to.keep <- num.cells > 10 & ave.counts >= 10^-5  & ave.counts <10^3  
      summary(genes.to.keep)
      sce <- sce[genes.to.keep, ]
      ggs = ggs[genes.to.keep]
      
      rownames(sce) = make.names(rownames(sce), unique = TRUE)
      
      geneDup = data.frame(gene = ggs, gg.uniq = rownames(sce), stringsAsFactors = FALSE)
      geneDup$geneSymbol = get_geneName(geneDup$gene)
      
      kk = which(geneDup$geneSymbol != geneDup$gg.uniq)
      
      subref = as.Seurat(sce, counts = "counts", data = "logcounts")
      rm(sce)
      Idents(subref) = as.factor(subref$celltypes)
      
      # saveRDS(geneDup, paste0(outDir, '/geneSymbol_duplication_inLRanalysis.rds'))
      
    }
    
    ##########################################
    # overview of scRNA-seq data 
    ##########################################
    check.QC.subref = FALSE
    if(check.QC.subref){
      subref$celltype = subref$celltypes
      
      # rerun the processing and umap
      subref <- FindVariableFeatures(subref, selection.method = "vst", nfeatures = 3000)
      subref <- ScaleData(subref)
      
      subref <- RunPCA(subref, verbose = FALSE, features = VariableFeatures(object = subref), weight.by.var = TRUE,
                       reduction.key = 'pca_')
      ElbowPlot(subref, ndims = 30)
      
      subref <- RunUMAP(subref, dims = 1:30, n.neighbors = 30, min.dist = 0.1, reduction = 'pca', 
                        reduction.key = 'umap_')
      DimPlot(subref, group.by = "celltype", label = TRUE, reduction = 'umap') 
      
      
      geneDup = readRDS(paste0(outDir, '/geneSymbol_duplication_inLRanalysis.rds'))
      kk = which(geneDup$geneSymbol != geneDup$gg.uniq)
      ggs = unique(geneDup$geneSymbol[kk])
      cat(length(ggs), ' genes have symbols issue \n')
      ligand_dup = ggs[ggs %in% lr_network$ligand]
      receptor_dup = ggs[ggs %in% lr_network$receptor]
      
      ligand_dup = geneDup[!is.na(match(geneDup$geneSymbol, ligand_dup)), ]
      ligand_dup = ligand_dup[order(ligand_dup$geneSymbol), ]
      
      receptor_dup = geneDup[!is.na(match(geneDup$geneSymbol, receptor_dup)), ]
      receptor_dup = receptor_dup[order(receptor_dup$geneSymbol), ]
      
      write.csv2(ligand_dup, file = paste0(outDir, '/geneDuplication_issue_ligands.csv'))
      write.csv2(receptor_dup, file = paste0(outDir, '/geneDuplication_issue_receptors.csv'))
      saveRDS(subref, file = paste0(outDir, '/seuratObject_snRNAseq_subset_for_NicheNet.rds')) 
      
    }
    
    ##########################################
    # double check the subref before running nichenet
    ##########################################
    #subref = readRDS(file = paste0(outDir, '/seuratObject_snRNAseq_subset_for_NicheNet.rds'))
    #subref = subset(x = subref, downsample = 1000) 
    table(subref$celltypes)
    subref$celltype = subref$celltypes
    seurat_obj = SetIdent(subref, value = "celltype")
    
    table(seurat_obj$celltypes)
    #seurat_obj$celltype = as.factor(seurat_obj$celltypes)
    Idents(seurat_obj) = seurat_obj$celltype
    rm(subref)
    
    ##########################################
    # Step 1:  Define the niches/microenvironments of interest
    ##########################################
    seurat_obj$celltype = as.factor(seurat_obj$celltype)
    table(seurat_obj$celltype)
    #celltypes_all = levels(seurat_obj$celltype)
    #celltypes_BZ = celltypes_all[grep('_BZ', celltypes_all)]
    #celltypes_noBZ = setdiff(celltypes_all, celltypes_BZ)
    #rm(celltypes_all)
    cat('cell types in BZ : ', celltypes_BZ, '\n')
    cat('cell types in RZ : ', celltypes_RZ, '\n')
    
    # receiver_cells = 'CM_BZ'
    cat('--- receiver cells in BZ: ', receivers_BZ, ' ---\n')
    cat('--- receiver cells in RZ: ', receivers_RZ, ' ---\n')
    
    lfc_cutoff = 0.15
    expression_pct = 0.10
    top_n_target = 200
    
    for(receiver in receivers_BZ)
    {
      # receiver = receivers_BZ[1]
      receiver_ctl = receivers_RZ;
      
      cat('-- define niches wiht receiver ', receiver, ' and control ', receiver_ctl, '--\n')
      
      # define the niches
      niches = list(
        "BZ_niche" = list(
          #"sender" = setdiff(celltypes_BZ, receiver_cells),
          "sender" = celltypes_BZ, # including autocrine
          "receiver" = receiver),
        "RZ_niche" = list(
          #"sender" = setdiff(celltypes_noBZ, gsub('_BZ',  '', receiver_cells)),
          "sender" = celltypes_RZ, 
          "receiver" = receiver_ctl)
      )
      
      res = prioritize_ligands_between_niches(seurat_obj, niches,
                                              expression_pct = expression_pct,
                                              lfc_cutoff = lfc_cutoff,
                                              top_n_target = top_n_target)
      
      saveRDS(res, file = paste0(outDir, '/nichenet_prioritization_tables_output',
                            '_receiverBZ.', receiver,
                            '_receiverRZ.', receiver_ctl, 
                            '_timepoint.', timepoint, '.rds'))
      
      # res = readRDS(file = paste0(outDir, '/nichenet_prioritization_tables_output',
      #                             '_receiverBZ.', receivers_BZ_timeSpecific,
      #                             '_receiverRZ.', receivers_RZ_timeSpecific, 
      #                             '_timepoint.', timepoint, '.rds'))
      
      prioritization_tables = res[[1]]
      output = res[[2]]
      
      ##########################################
      # 8). Visualization of the Differential NicheNet output 
      ##########################################
      top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
        dplyr::select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
        group_by(ligand) %>% 
        top_n(1, prioritization_score) %>% 
        ungroup() %>% 
        dplyr::select(ligand, receptor, niche) %>% 
        dplyr::rename(top_niche = niche)
      
      top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
        dplyr::select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
        group_by(ligand, receptor) %>% 
        top_n(1, prioritization_score) %>% 
        ungroup() %>% 
        dplyr::select(ligand, receptor, niche) %>% 
        dplyr::rename(top_niche = niche)
      
      for(ntop in c(50, 100, 150, 200))
      {
        # ntop = 50
        cat('ntop -- ', ntop, '\n')
        ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
          dplyr::select(niche, sender, receiver, ligand, prioritization_score) %>% 
          group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% 
          ungroup() %>% distinct() %>% 
          inner_join(top_ligand_niche_df) %>% 
          filter(niche == top_niche) %>% 
          group_by(niche) %>% 
          top_n(n = ntop, prioritization_score) %>% 
          ungroup() # get the top50 ligands per niche
        
        receiver_oi = receiver
        
        filtered_ligands = ligand_prioritized_tbl_oi %>% 
          filter(receiver == receiver_oi) %>% 
          pull(ligand) %>% unique()
        
        prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
          filter(ligand %in% filtered_ligands) %>% 
          dplyr::select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
          distinct() %>% 
          inner_join(top_ligand_receptor_niche_df) %>% 
          group_by(ligand) %>% 
          filter(receiver == receiver_oi) %>% 
          top_n(2, prioritization_score) %>% ungroup() 
        
        lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, 
                                                 prioritized_tbl_oi, 
                                                 prioritization_tables$prioritization_tbl_ligand_receptor, 
                                                 plot_legend = FALSE, 
                                                 heights = NULL, widths = NULL)
        lfc_plot
        
        ggsave(paste0(outDir, '/Ligand_receptors_LFC', 
                      '_receiverBZ.', receiver,
                      '_receiverRZ.', receiver_ctl, 
                      '_timepoint.', timepoint,
                      '_ntop.', ntop, '.pdf'), 
               width = 12, height = 12*ntop/50, limitsize = FALSE)
        
        
        exprs_activity_target_plot = 
          make_ligand_activity_target_exprs_plot(receiver_oi, 
                                                 prioritized_tbl_oi,  
                                                 prioritization_tables$prioritization_tbl_ligand_receptor,  
                                                 prioritization_tables$prioritization_tbl_ligand_target, 
                                                 output$exprs_tbl_ligand,  
                                                 output$exprs_tbl_target, 
                                                 lfc_cutoff, 
                                                 ligand_target_matrix, 
                                                 plot_legend = FALSE, 
                                                 heights = NULL, widths = NULL)
        
        exprs_activity_target_plot$combined_plot
        ggsave(paste0(outDir, '/Combined_plots_ligand_noFiltering',
                      '_receiverBZ.', receiver,
                      '_receiverRZ.', receiver_ctl, 
                      '_timepoint.', timepoint,
                      '_ntop.', ntop, '.pdf'), 
               width = 60, height = 12*ntop/50, limitsize = FALSE)
        
        ## Circos plot of prioritized ligand-receptor pairs
        Make.Circos.plot = FALSE
        if(Make.Circos.plot){
          filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% 
            top_n(15, prioritization_score) %>% pull(ligand) %>% unique()
          
          prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
            filter(ligand %in% filtered_ligands) %>% 
            dplyr::select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
            distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% 
            group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(2, prioritization_score) %>% ungroup() 
          
          colors_sender = brewer.pal(n = prioritized_tbl_oi$sender %>% unique() %>% sort() %>% 
                                       length(), name = 'Spectral') %>% 
            magrittr::set_names(prioritized_tbl_oi$sender %>% unique() %>% sort())
          colors_receiver = c("lavender")  %>% magrittr::set_names(prioritized_tbl_oi$receiver %>% 
                                                                     unique() %>% sort())
          
          circos_output = make_circos_lr(prioritized_tbl_oi, colors_sender, colors_receiver)
          
        }
        
        
      }
      
    }
    
  }
  
  
  
}


prioritize_ligands_between_niches = function(seurat_obj, 
                                             niches,
                                             assay_oi = "RNA",
                                             expression_pct = 0.10,
                                             lfc_cutoff = 0.15,
                                             include_spatial_info_sender = FALSE,
                                             include_spatial_info_receiver = FALSE,
                                             top_n_target = 250
                                             )
{
  ##########################################
  # step 0): Nichenet data loaded already 
  ##########################################
  # NicheNet’s ligand-target prior model
  ligand_target_matrix = readRDS(paste0(dataPath_nichenet,  "ligand_target_matrix.rds"))
  ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
  
  # ligand-receptor network, Putative ligand-receptor links from NicheNet
  lr_network = readRDS(paste0(dataPath_nichenet, "lr_network.rds"))
  lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
  lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% 
    distinct(ligand, receptor, bonafide)
  
  # If wanted, users can remove ligand-receptor interactions that were predicted based on 
  # protein-protein interactions and 
  # only keep ligand-receptor interactions that are described in curated databases.
  # lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  head(lr_network)
  
  # ## weighted integrated networks 
  weighted_networks = readRDS(paste0(dataPath_nichenet,  "weighted_networks.rds"))
  head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
  head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
  
  weighted_networks_lr = weighted_networks$lr_sig %>% 
    dplyr::rename(ligand = from, receptor = to) %>%
    inner_join(lr_network %>% 
                 distinct(ligand,receptor), by = c("ligand","receptor"))
  
  
  ##########################################
  # step 2. Calculate differential expression between the niches 
  # issue solved because find markers for all genes, should run only for receptors and ligands
  ##########################################
  # only ligands important for sender cell types
  DE_sender = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), 
                                 niches = niches, 
                                 type = "sender", 
                                 assay_oi = assay_oi)
  
  # only receptors now, later on: DE analysis to find targets
  DE_receiver = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), 
                                   niches = niches, 
                                   type = "receiver", 
                                   assay_oi = assay_oi)
  
  DE_sender = DE_sender %>% 
    mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), 
                               ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
  DE_receiver = DE_receiver %>% 
    mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), 
                               ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
  
  # expected percentage of cells expressing the genes 
  DE_sender_processed = process_niche_de(DE_table = DE_sender, 
                                         niches = niches, 
                                         expression_pct = expression_pct, 
                                         type = "sender")
  DE_receiver_processed = process_niche_de(DE_table = DE_receiver, 
                                           niches = niches, 
                                           expression_pct = expression_pct, 
                                           type = "receiver")
  
  # Combine sender-receiver DE based on L-R pairs:
  # As mentioned above, DE of ligands from one sender cell type is determined be calculating DE 
  # between that cell type, and all the sender cell types of the other niche. 
  # To summarize the DE of ligands of that cell type we have several options: 
  # we could take the average LFC, but also the minimum LFC compared to the other niche. 
  # We recommend using the minimum LFC, because this is the strongest specificity measure of ligand expression, 
  # because a high min LFC means that a ligand is more strongly expressed in the cell type of niche 1 
  # compared to all cell types of niche 2 (in contrast to a high average LFC, 
  # which does not exclude that one or more cell types in niche 2 also strongly express that ligand).
  specificity_score_LR_pairs = "min_lfc"
  DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, 
                                                  DE_receiver_processed, 
                                                  lr_network, 
                                                  specificity_score = specificity_score_LR_pairs)
  
  
  ##########################################
  # Step 3) Optional: Calculate differential expression between the different spatial regions 
  ##########################################
  # if not spatial info to include: put this to false 
  
  # user adaptation required on own dataset
  spatial_info = tibble(celltype_region_oi = "CAF_High", 
                        celltype_other_region = "myofibroblast_High", 
                        niche =  "pEMT_High_niche", celltype_type = "sender") 
  specificity_score_spatial = "lfc"
  
  # this is how this should be defined if you don't have spatial info (mock spatial info)
  if(include_spatial_info_sender == FALSE & include_spatial_info_receiver == FALSE){
    spatial_info = tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% 
      mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
  } 
  
  
  if(include_spatial_info_sender == TRUE){
    sender_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% 
                                               subset(features = lr_network$ligand %>% unique()), 
                                             spatial_info = spatial_info %>% filter(celltype_type == "sender"))
    sender_spatial_DE_processed = process_spatial_de(DE_table = sender_spatial_DE, 
                                                     type = "sender", 
                                                     lr_network = lr_network, 
                                                     expression_pct = expression_pct, 
                                                     specificity_score = specificity_score_spatial)
    
    # add a neutral spatial score for sender celltypes in which the spatial is not known / not of importance
    sender_spatial_DE_others = get_non_spatial_de(niches = niches, 
                                                  spatial_info = spatial_info, 
                                                  type = "sender", 
                                                  lr_network = lr_network)
    sender_spatial_DE_processed = sender_spatial_DE_processed %>% bind_rows(sender_spatial_DE_others)
    
    sender_spatial_DE_processed = sender_spatial_DE_processed %>% 
      mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
    
  }else{
    #add a neutral spatial score for all sender celltypes (for none of them, spatial is relevant in this case)
    sender_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, 
                                                     type = "sender", lr_network = lr_network)
    sender_spatial_DE_processed = sender_spatial_DE_processed %>% 
      mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  
    
  }
  
  if(include_spatial_info_receiver == TRUE){
    receiver_spatial_DE = calculate_spatial_DE(seurat_obj = seurat_obj %>% 
                                                 subset(features = lr_network$receptor %>% unique()), 
                                               spatial_info = spatial_info %>% 
                                                 filter(celltype_type == "receiver"))
    receiver_spatial_DE_processed = process_spatial_de(DE_table = receiver_spatial_DE, 
                                                       type = "receiver", 
                                                       lr_network = lr_network, 
                                                       expression_pct = expression_pct, 
                                                       specificity_score = specificity_score_spatial)
    
    # add a neutral spatial score for receiver celltypes in which the spatial is not known / not of importance
    receiver_spatial_DE_others = get_non_spatial_de(niches = niches, 
                                                    spatial_info = spatial_info, 
                                                    type = "receiver", 
                                                    lr_network = lr_network)
    receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% bind_rows(receiver_spatial_DE_others)
    
    receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% 
      mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
    
  }else{
    # # add a neutral spatial score for all receiver celltypes (for none of them, spatial is relevant in this case)
    receiver_spatial_DE_processed = get_non_spatial_de(niches = niches, 
                                                       spatial_info = spatial_info, 
                                                       type = "receiver", 
                                                       lr_network = lr_network)
    receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% 
      mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  }
  
  
  ##########################################
  # Step 4). Calculate ligand activities and infer active ligand-target links
  ##########################################
  ## first define the target genes by DE
  # It is always useful to check the number of genes in the geneset before doing the ligand activity analysis. 
  # We recommend having between 20 and 1000 genes in the geneset of interest, 
  # and a background of at least 5000 genes for a proper ligand activity analysis. 
  # If you retrieve too many DE genes, it is recommended to use a higher lfc_cutoff threshold. 
  # We recommend using a cutoff of 0.15 if you have > 2 receiver cells/niches to compare and 
  # use the min_lfc as specificity score. 
  # If you have only 2 receivers/niche, we recommend using a higher threshold (such as using 0.25). 
  # If you have single-cell data like Smart-seq2 with high sequencing depth, 
  #we recommend to also use higher threshold.
  # Smart-seq2 data and only 2 niches to compare, 
  # so we will use a stronger LFC threshold to keep less DE genes, but more trustworthy ones.
  # lfc_cutoff = 0.15 # recommended for 10x as min_lfc cutoff. 
  specificity_score_targets = "min_lfc"
  
  # here DE all genes in the receiver cells, not only for ligand and receptors as before
  DE_receiver_targets = calculate_niche_de_targets(seurat_obj = seurat_obj, 
                                                   niches = niches, 
                                                   lfc_cutoff = lfc_cutoff, 
                                                   expression_pct = expression_pct, 
                                                   assay_oi = assay_oi)
  
  DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, 
                                                             niches = niches, 
                                                             expression_pct = expression_pct, 
                                                             specificity_score = specificity_score_targets)
  
  background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
  geneset_niche1 = DE_receiver_processed_targets %>% 
    filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & 
             target_present == 1) %>% pull(target) %>% unique()
  geneset_niche2 = DE_receiver_processed_targets %>% 
    filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & 
             target_present == 1) %>% pull(target) %>% unique()
  
  # Good idea to check which genes will be left out of the ligand activity analysis 
  # (=when not present in the rownames of the ligand-target matrix).
  # If many genes are left out, this might point to some issue in the gene naming 
  # (eg gene aliases and old gene symbols, bad human-mouse mapping)
  geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
  geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
  
  length(geneset_niche1)
  length(geneset_niche2)
  
  niche_geneset_list = list(
    "BZ_niche" = list(
      "receiver" = niches[[1]]$receiver,
      "geneset" = geneset_niche1,
      "background" = background),
    "noBZ_niche" = list(
      "receiver" = niches[[2]]$receiver,
      "geneset" = geneset_niche2 ,
      "background" = background)
  )
  
  # get_ligand_activites_targets calling the function predict_ligand_activities for each niche
  ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, 
                                                            ligand_target_matrix = ligand_target_matrix, 
                                                            top_n_target = top_n_target)
  
  head(ligand_activities_targets)
  ligand_activities_targets %>% arrange(-activity) %>% filter(receiver %in% niches$BZ_niche$receiver)
  
  ##########################################
  # step 5. Calculate (scaled) expression of ligands, receptors and targets 
  # across cell types of interest (log expression values and expression fractions)
  # we will calculate average (scaled) expression, and fraction of expression, of 
  # ligands, receptors, and target genes across all cell types of interest
  ##########################################
  features_oi = union(lr_network$ligand, lr_network$receptor) %>% 
    union(ligand_activities_targets$target) %>% setdiff(NA)
  
  require(dplyr)
  # save dotplot 
  dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()), 
                                             features = features_oi, assay = assay_oi))
  
  exprs_tbl = dotplot$data %>% as_tibble()
  exprs_tbl = exprs_tbl %>% 
    dplyr::rename(celltype = id, gene = features.plot, expression = avg.exp, 
                  expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
    mutate(fraction = fraction/100) %>% as_tibble() %>% 
    dplyr::select(celltype, gene, expression, expression_scaled, fraction) %>% 
    distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))
  
  exprs_tbl_ligand = exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% 
    dplyr::rename(sender = celltype, ligand = gene, ligand_expression = expression, 
                  ligand_expression_scaled = expression_scaled, ligand_fraction = fraction)
  
  exprs_tbl_receptor = exprs_tbl %>% 
    filter(gene %in% lr_network$receptor) %>% 
    dplyr::rename(receiver = celltype, receptor = gene, receptor_expression = expression, 
                  receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
  exprs_tbl_target = exprs_tbl %>% 
    filter(gene %in% ligand_activities_targets$target) %>% 
    dplyr::rename(receiver = celltype, target = gene, target_expression = expression, 
                  target_expression_scaled = expression_scaled, target_fraction = fraction)
  
  exprs_tbl_ligand = exprs_tbl_ligand %>% 
    mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% 
    mutate(ligand_fraction_adapted = ligand_fraction) %>% 
    mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% 
    mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))
  
  exprs_tbl_receptor = exprs_tbl_receptor %>% 
    mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled)) %>% 
    mutate(receptor_fraction_adapted = receptor_fraction) %>% 
    mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct) %>% 
    mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))
  
  ##########################################
  # step 6. Expression fraction and receptor
  # score ligand-receptor interactions based on expression strength of the receptor, 
  # in such a way that we give higher scores to the most strongly expressed receptor of a certain ligand, 
  # in a certain celltype. 
  # This will not effect the rank of individual ligands later on, 
  # but will help in prioritizing the most important receptors per ligand 
  # (next to other factors regarding the receptor - see later).
  ##########################################
  exprs_sender_receiver = lr_network %>% 
    inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
    inner_join(exprs_tbl_receptor, by = c("receptor")) %>% 
    inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))
  
  ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% 
    group_by(ligand, receiver) %>% 
    mutate(rank_receptor_expression = dense_rank(receptor_expression), 
           rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% 
    mutate(ligand_scaled_receptor_expression_fraction = 0.5* ((rank_receptor_fraction / max(rank_receptor_fraction)) + 
                                                                ((rank_receptor_expression / max(rank_receptor_expression))) ))  %>% 
    distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% 
    distinct() %>% ungroup() 
  
  ##########################################
  # step 7. Prioritization of ligand-receptor and ligand-target links 
  ##########################################
  prioritizing_weights = c(
    "scaled_ligand_score" = 5, # # niche-specific expression of the ligand: Recommended 5 (between 1-5)
    # scaled ligand expression in one sender compared to the other cell types in the dataset
    "scaled_ligand_expression_scaled" = 1,
    "ligand_fraction" = 1, # Ligands expressed in a smaller fraction of cells of cell type than cutoff (default: 0.10)
    "scaled_ligand_score_spatial" = 0, 
    # receptor DE score: niche-specific expression, Recommended 0.5 (>=0.5 and lower than "scaled_ligand_score")
    "scaled_receptor_score" = 1, 
    "scaled_receptor_expression_scaled" = 0.5, # Recommended weight: 0.5
    # Receptors that are expressed in a smaller fraction of cells of a cell type than exprs_cutoff(default: 0.10) 
    # will get a lower ranking, proportional to their fraction, Recommended weight: 1. 
    "receptor_fraction" = 1, 
    # this factor let us give higher weights to the most highly expressed receptor of a ligand in the receiver.
    # Recommended value: 1 (minimum: 0.5)
    "ligand_scaled_receptor_expression_fraction" = 1,
    "scaled_receptor_score_spatial" = 0,
    # Absolute ligand activity: to further prioritize ligand-receptor pairs based on their predicted effect of 
    # the ligand-receptor interaction on the gene expression in the receiver cell type - 
    # prioritizing_weights argument: "scaled_activity". Recommended weight: 0, 
    # unless absolute enrichment of target genes is of specific interest.
    "scaled_activity" = 0, 
    # Normalized ligand activity: to further prioritize ligand-receptor pairs based on their predicted effect of
    # the ligand-receptor interaction on the gene expression in the receiver cell type, Recommended weight: >=1.
    "scaled_activity_normalized" = 4,
    "bona_fide" = 1)
  
  output = list(DE_sender_receiver = DE_sender_receiver, 
                ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, 
                sender_spatial_DE_processed = sender_spatial_DE_processed, 
                receiver_spatial_DE_processed = receiver_spatial_DE_processed,
                ligand_activities_targets = ligand_activities_targets, 
                DE_receiver_processed_targets = DE_receiver_processed_targets, 
                exprs_tbl_ligand = exprs_tbl_ligand,  
                exprs_tbl_receptor = exprs_tbl_receptor, 
                exprs_tbl_target = exprs_tbl_target
  )
  
  prioritization_tables = get_prioritization_tables(output, prioritizing_weights)
  
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
    filter(receiver == niches[[1]]$receiver) %>% 
    head(10)
  
  prioritization_tables$prioritization_tbl_ligand_target %>% 
    filter(receiver == niches[[1]]$receiver) %>% 
    head(10)
  
  return(list(prioritization_tables, output))
  
  
}

########################################################
########################################################
# Section III : combine LIANA and NicheNet 
# 
########################################################
########################################################
##########################################
# Combine Liana and NicheNet  
# # original code from https://saezlab.github.io/liana/articles/liana_nichenet.html
##########################################
run_liana_nitchenet = function()
{
  library(tidyverse)
  library(liana)
  library(nichenetr)
  library(Seurat)
  library(ggrepel)
  library(cowplot)
  options(timeout=180) # required for downloading single-cell expression on slow connection
  
  # single-cell expression matrix described in Puram et al. 2017
  hnscc_expression <-  readRDS(url("https://zenodo.org/record/3260758/files/hnscc_expression.rds"))
  expression <- hnscc_expression$expression
  sample_info <- hnscc_expression$sample_info
  colnames(sample_info) <- make.names(colnames(sample_info))
  
  # filter samples based on vignette's information and add cell type
  tumors_remove <-  c("HN10", "HN", "HN12", "HN13", "HN24", "HN7", "HN8", "HN23")
  sample_info <- sample_info %>%
    subset( !(tumor %in% tumors_remove) & Lymph.node == 0) %>%
    # fix some cell type identity names
    mutate(cell_type = ifelse(classified..as.cancer.cell == 1, "Tumor", non.cancer.cell.type)) %>%
    subset(cell_type %in% c("Tumor", "CAF"))
  
  # cell ID as rownames
  rownames(sample_info) <- sample_info$cell
  
  # subset expression to selected cells
  expression <- expression[sample_info$cell, ]
  
  # model weights
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
  
  # gene set of interest
  geneset_oi <- read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), 
                         col_types = cols(), col_names = "gene") %>%
    pull(gene) %>%
    .[. %in% rownames(ligand_target_matrix)]
  
  # create seurat object
  seurat_object <- Seurat::CreateAssayObject(counts = expm1(t(expression))) %>%
    Seurat::CreateSeuratObject(., meta.data = sample_info) %>%
    Seurat::NormalizeData()
  
  # set cell identity to cell type
  Idents(seurat_object) <- seurat_object@meta.data$cell_type
  
  liana_results <- liana_wrap(seurat_object) %>%
    liana_aggregate()
  
  
  # filter results to cell types of interest
  caf_tumor_results <- liana_results %>%
    subset(source == "CAF" & target == "Tumor")
  
  # filter results to top N interactions
  n <- 50
  top_n_caf_tumor <- caf_tumor_results %>%
    arrange(aggregate_rank) %>%
    slice_head(n = n) %>%
    mutate(id = fct_inorder(paste0(ligand, " -> ", receptor)))
  
  # visualize median rank
  top_n_caf_tumor %>%
    ggplot(aes(y = aggregate_rank, x = id)) +
    geom_bar(stat = "identity") +
    xlab("Interaction") + ylab("LIANA's aggregate rank") +
    theme_cowplot() +
    theme(axis.text.x = element_text(size = 8, angle = 60, hjust = 1, vjust = 1))
  
  
  # get ligands and filter to those included in NicheNet's ligand-target matrix
  ligands <- unique(top_n_caf_tumor$ligand)
  ligands <- ligands[ligands %in% colnames(ligand_target_matrix)]
  ligands
  
  background_genes <- expression[sample_info$cell[sample_info$cell_type == "Tumor"], ] %>%
    apply(2,function(x){10*(2**x - 1)}) %>%
    apply(2,function(x){log2(mean(x) + 1)}) %>%
    .[. >= 4] %>%
    names()
  
  nichenet_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_genes,
    ligand_target_matrix = ligand_target_matrix, potential_ligands = ligands
  )
  
  # prepare data for visualization
  vis_liana_nichenet <- top_n_caf_tumor %>%
    inner_join(nichenet_activities, by = c("ligand" = "test_ligand")) %>%
    arrange(pearson) %>%
    mutate(ligand = fct_inorder(ligand))
  
  # prepare NicheNet figure
  nichenet_scores_plot <- vis_liana_nichenet %>%
    group_by(ligand) %>%
    summarize(pearson = mean(pearson)) %>%
    ggplot(aes(y = ligand, x = pearson)) +
    geom_bar(stat = "identity") +
    ggtitle("NicheNet") +
    xlab("Pearson's score") +
    theme_cowplot() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.line.y = element_line(color = "white"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
  
  # prepare LIANA figure
  liana_receptor_heatmap <- vis_liana_nichenet %>%
    ggplot(aes(y = ligand, x = receptor, fill = aggregate_rank)) +
    geom_tile() +
    theme_cowplot() +
    ggtitle("LIANA") +
    ylab("Ligand") + xlab("Receptor") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_line(colour = "gray", linetype = 2),
          legend.position = "left")
  
  # combine plots
  plot_grid(liana_receptor_heatmap, nichenet_scores_plot,
            align = "h", nrow = 1, rel_widths = c(0.8,0.3))
  
}
