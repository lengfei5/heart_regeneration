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

##########################################
# main function of LIANA
##########################################
# original code from https://saezlab.github.io/liana/articles/liana_tutorial.html
run_LIANA = function(refs,
                     timepoint_specific = TRUE,
                     include_autocrine = FALSE,
                     celltypes = NULL,
                     celltypes_timeSpecific = NULL,
                     #receiver_cells = list(c('CM_IS')),
                     outDir = '../results/Ligand_Receptor_analysis',
                     species = 'ax6',
                     ntop = 100,
                     plotGeneExpression.LR = FALSE,
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
  
  # celltypes_timeSpecific = celltypes_BZ_timeSpecific[[1]]
  cat('output directory : ', outDir,  '\n')
  system(paste0('mkdir -p ', outDir))
  
  
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
        #cat(length(celltypes_timeSpecific), 'subtype pairs specified \n')
        cat(names(celltypes_timeSpecific), ' as receivers \n')
        
        for(n in 1:length(celltypes_timeSpecific))
        {
          celltypes = unique(c(celltypes_timeSpecific[[n]], names(celltypes_timeSpecific)[n]))
          mm = match(celltypes, refs$celltypes)
          if(length(which(is.na(mm)))>0){
            stop(names(celltypes_timeSpecific)[n], 
                 ': some selected cell types not found in refs -- ', celltypes[which(is.na(mm))])
          }else{
            cat('-----------------\n')
            cat("receiver : ", names(celltypes_timeSpecific)[n],'\n')
            cat("senders : ", celltypes, '\n')
            cat('--all selected cell types were found in the refs --\n')
          }
        }
      }
    }
    
  }
  
  # prepre for the LIANA loop
  Idents(refs) = as.factor(refs$celltypes)
  
  if(!timepoint_specific){
    # celltypes = rownames(pairs)
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
      #timepoint = names(celltypes_timeSpecific)[n]
      receivers = names(celltypes_timeSpecific)[n]
      
      cat('---------- run LIANA for ', receivers, ' -----------\n')
      cat('---------- senders : ', celltypes, '\n')
      #cat('---- receivers : ', receivers, "\n")
      
      
      if(include_autocrine) {
        cat('-- autocrine is considerede -- \n')
        celltypes = unique(c(celltypes, receivers))
      }
      
      celltypes_sel = unique(c(celltypes, receivers))
      
      subref = subset(refs, cells = colnames(refs)[!is.na(match(refs$celltypes, celltypes_sel))])
      subref$celltypes = droplevels(as.factor(subref$celltypes))
      table(subref$celltypes)
      #subref = subset(x = subref, downsample = 1000)
      
      #cat('celltype to consider -- ', names(table(subref$celltypes)), '\n')
      
      Idents(subref) = subref$celltypes
      
      # Run liana
      run_LIANA_defined_celltype(subref,
                                 celltypes = celltypes,
                                 receivers = receivers,
                                 additionalLabel = paste0('_', receivers),
                                 outDir = outDir, 
                                 species = species,
                                 Test.geneExpresion = plotGeneExpression.LR)

    }
  
  }
  
  
}


run_LIANA_defined_celltype = function(subref,
                                      celltypes,
                                      receivers = NULL,
                                      additionalLabel = '_fixedCelltypes',
                                      outDir = outDir,
                                      species = species,
                                      Test.geneExpresion = TRUE)
{
  source('functions_scRNAseq.R')
  
  sce <- as.SingleCellExperiment(subref)
  colLabels(sce) = as.factor(sce$celltypes)
  
  if(species == 'ax6'){
    rownames(sce) = toupper(get_geneName(rownames(sce)))
  }else{
    rownames(sce) = toupper(rownames(sce))
  }
 
  
  ave.counts <- calculateAverage(sce, assay.type = "counts")
  
  #hist(log10(ave.counts), breaks=100, main="", col="grey80",
  #     xlab=expression(Log[10]~"average count"))
  
  num.cells <- nexprs(sce, byrow=TRUE)
  smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
                xlab=expression(Log[10]~"average count"))
  
  # detected in >= 5 cells, ave.counts >=5 but not too high
  genes.to.keep <- num.cells > 20 & ave.counts >= 10^-4  & ave.counts <10^2  
  summary(genes.to.keep)
  
  sce <- sce[genes.to.keep, ]
  
  ## run the liana wrap function by specifying resource and methods
  # Resource currently included in OmniPathR (and hence `liana`) include:
  show_resources()
  # Resource currently included in OmniPathR (and hence `liana`) include:
  show_methods()
  
  if(min(exec("logcounts", sce)) < 0){
    xx = logcounts(sce)
    xx[which(xx<0)] = 0
    xx = Matrix(xx, sparse = TRUE)
    logcounts(sce) = xx
    rm(xx)
  }
  
  liana_test <- liana_wrap(sce,  
                           # method = c("natmi", "connectome", "logfc", "sca", "cytotalk"),
                           method = c("natmi", "connectome", "logfc", "sca"),
                           #resource = c("Consensus", 'CellPhoneDB', "OmniPath", "LRdb", "CellChatDB",  
                           # "CellTalkDB"), 
                           resource = c("Consensus", 'CellPhoneDB', "CellChatDB"),
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
  
  df_test = liana_test %>% filter(target %in% receivers & source %in% celltypes) %>% as.data.frame() 
  write.table(df_test, file = paste0(outDir, '/res_lianaTest_Consensus', additionalLabel, '.txt'), 
              sep = '\t', quote = FALSE)
  
  
  if(Test.geneExpresion){
    
    system(paste0('mkdir -p ', outDir, '/examples_plotted'))
    
    DimPlot(subref, group.by = 'subtypes', label = TRUE, repel = TRUE) + NoLegend()
    ggsave(filename = paste0(outDir, '/examples_plotted/snRNA_subtypes_labeled', additionalLabel, '.pdf'),
           width = 10, height = 8)
    
    geneNames = get_geneName(rownames(subref))
    
    for(mm in 1:100)
    {
      cat('-- example number : ', mm, '\n')
      #FeaturePlot(refs, features = rownames(refs)[grep('VEGFC|VIPR2', rownames(refs))])
      features = c(df_test$ligand.complex[mm], df_test$receptor.complex[mm])
      features = sapply(features, function(x) unlist(strsplit(as.character(x), '_')))
      
      features = rownames(subref)[!is.na(match(geneNames, features))]
      FeaturePlot(subref, features = features)
      
      ggsave(filename = paste0(outDir, '/examples_plotted/featurePlots_ligand_receptor_', additionalLabel, 
                               '_top_', mm, 
                               "_" ,df_test$source[mm], 
                               '_to_', df_test$target[mm],  '.pdf'),
             width = 12, height = 10)
    }
    
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
# aggreate LIANA outputs and visiualize them
##########################################
aggregate_output_LIANA = function(liana_out)
{
  # liana_out = outDir
  cat('liana output folder --- ', liana_out, '\n')
  xlist = list.files(path = liana_out, pattern = '*.txt', full.names = TRUE)
  xlist = xlist[grep('allPairs_LIANA', xlist, invert = TRUE)]
  # for(n in 1:length(xlist))
  # {
  #   # n = 2
  #   cat(n, ' -- ', basename(xlist[n]), '\n')
  #   test = read.table(xlist[n], header = TRUE)
  #   
  #   ## aggregate_rank interpretation: https://saezlab.github.io/liana/articles/liana_tutorial.html
  #   ## RRA scores can be interpreted as p-values and 
  #   ## interactions which are ranked consistently higher than random are assigned low scores/p-values.
  #   test = data.frame(test[, c(1:6, 
  #                              which(colnames(test) == 'natmi.edge_specificity'), 
  #                              which(colnames(test) == 'sca.LRscore'))], stringsAsFactors = FALSE)
  #   
  #   write.table(test, file = paste0(liana_out, '/',  gsub('.txt', '_simplied.txt', basename(xlist[n]))),
  #               quote = FALSE, col.names = TRUE, row.names = FALSE)
  #   
  #   
  # }
  
  cat('concate the tables \n')
  res = c()
  for(n in 1:length(xlist))
  {
    # n = 1
    test = read.table(file = xlist[n])
    test = test[, c(1:5, which(colnames(test) == 'natmi.edge_specificity'), 
                    which(colnames(test) == 'sca.LRscore'))]
    colnames(test)[1:4] = c('sender', 'receiver', 'ligand', 'receptor')
    
    kk = which(test$ligand == 'GAS6' & test$receptor == 'AXL')[1]
    if(!is.na(kk)){
      cat('receiver: ', unique(test$receiver), '--',
          'sender: ', test$sender[kk], 
          'rank :', kk, '\n')
      #test = test[1:ntop, c(1:5)]
    }
    
    res = rbind(res, test)
  }
  
  saveRDS(res, file = paste0(liana_out, '/test_LianaOut_for_circoPlot.rds'))
  
  return(res)
  
  
}

assembly.liana.plot = function(output)
{
  ##########################################
  # Test some other packages
  ##########################################
  require(cellcall)
  
  res = readRDS(file = paste0(resDir, '/test_LianaOut_for_circoPlot.rds'))
  #colors_use <- RColorBrewer::brewer.pal(n=length(unique(c(res$sender, res$receiver))),"Set2")
  
  subtypes = unique(c(res$sender, res$receiver))
  cell_color = data.frame(color = distinctColorPalette(length(subtypes)), stringsAsFactors = FALSE)
  rownames(cell_color) <- subtypes
  
  ViewInterCircos(object = mt, font = 2, cellColor = cell_color, 
                  lrColor = c("#F16B6F", "#84B1ED"),
                  arr.type = "big.arrow",arr.length = 0.04,
                  trackhight1 = 0.05, 
                  slot="expr_l_r_log2_scale",
                  linkcolor.from.sender = TRUE,
                  linkcolor = NULL, 
                  gap.degree = 2,
                  order.vector=c('ST', "SSC", "SPGing", "SPGed"),
                  trackhight2 = 0.032, 
                  track.margin2 = c(0.01,0.12), 
                  DIY = FALSE)
  
  
  
  
  ## test celltalker
  library(celltalker)
  xx = readRDS(file = paste0(resDir, '/test_LianaOut_for_circoPlot.rds'))
  xx$sender = gsub('_', '.', xx$sender)
  xx$receiver = gsub('_', '.', xx$receiver)
  xx$interaction_pairs = paste0(xx$sender, '_', xx$receiver)
  xx$interaction = paste0(xx$ligand, '_', xx$receptor)
  xx = xx[,c(8, 9, 7, 1:5)]
  colnames(xx)[3:5] = c('value', 'cell_type1', 'cell_type2')
  
  top_stats_xx <- xx %>% as_tibble() %>%
    #mutate(fdr=p.adjust(p_val,method="fdr")) %>%
    #filter(fdr<0.05) %>%
    group_by(cell_type2) %>%
    top_n(-50, aggregate_rank) %>%
    ungroup()
  
  #top_stats_xx = as_tibble(xx)
  subtypes = unique(c(top_stats_xx$cell_type1, top_stats_xx$cell_type2))
  colors_use <- distinctColorPalette(length(subtypes))
  
  
  svg(paste0(resDir, "/ligand_target_circos_day4_test.svg"), width = 10, height = 10)
  
  source("functions_cccInference.R")
  circos_plot_customized(ligand_receptor_frame=top_stats_xx,
                         cell_group_colors=colors_use,
                         ligand_color="#84B1ED",
                         receptor_color="#F16B6F",
                         cex_outer=0.7,
                         cex_inner=0.1,
                         link.lwd=0.1, 
                         arr.length=0.1, 
                         arr.width=0.05
  )
  dev.off()
  
  
  
  
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

##########################################
# circo plot from connectome 
##########################################
#' CircosPlot
#'
#' Plotting function to make Circos plots using the circlize package, following the vignette by the Saeys Lab at: https://github.com/saeyslab/nichenetr/blob/master/vignettes/circos.md Note that this plotting type is incompatible with edges where the ligand and the receptor are the exact same gene.
#'
#' @param connectome A connectomic object, ideally filtered to only edges of interest.
#' @param weight.attribute Column to use to define edgeweights for network analysis. 'weight_sc' or 'weight_norm'. Defaults to 'weight_sc'. If 'weight_sc', function will automatically filter at min.z = 0 to remove negative source/sink values.
#' @param cols.use Optional. Colors for plotting nodes.
#' @param min.z Minimum z-score for ligand and receptor.
#' @param lab.cex Text size for gene names
#' @param balanced.edges Edges in this plot can change thickness along their length. This parameter decides whether to scale edges by a single edgeweight (chosen in weight.attribute) or by the separate cell-specific ligand and receptor values.  Default balanced (TRUE).  If FALSE, the edges will expand or contract to join ligand weight to receptor weight.
#' @param edge.color.by.source Default TRUE - edges will be colored by their source cell type. If false, edges will be colored by receiving cell instead.
#' @param small.gap Default 1. Amount of distance between sectors.  If the number of edges is very large, this will have to be reduced in size.
#' @param big.gap Default 10. Amount of distance between the source cells and the target cells (top and bottom arc of graph).  If the number of edges is very large, this can be reduced in size in addition to 'small.gap'
#' @param title Character string for title of plot.
#' @param ... Arguments passed to FilterConnectome
#' @export

my_CircosPlot <- function(connectome,
                       weight.attribute = 'weight_sc',
                       cols.use = NULL,
                       min.z = NULL,
                       lab.cex = 1,
                       balanced.edges = T,
                       edge.color.by.source = T,
                       small.gap = 1,
                       big.gap = 10,
                       title = NULL,...)
{
  # connectome = test; weight.attribute = 'weight_norm'; sources.include = cells.of.interest;
  # targets.include = cells.of.interest;min.z = NULL;cols.use = cell_color;edge.color.by.source = T;
  # balanced.edges = T;small.gap = 1;big.gap = 10;lab.cex = 1;
  
  library(tidyverse)
  library(circlize)
  library(dplyr)
  library(scales)
  library(ComplexHeatmap)
  
  # If (weight.attribute != 'weight_norm'){
  #if (weight.attribute == 'weight_sc' & is.null(min.z)){
  #  connectome <- FilterConnectome(connectome, remove.na = T,min.z = 0,...)
  #}else{
  #  connectome <- FilterConnectome(connectome,remove.na = T,min.z = min.z,...)
  #}
  
  #}
  # Pull the dataframe of interest for plotting and format with weight as third column
  connectome$lig.stash <- as.character(connectome$ligand)
  connectome$rec.stash <- as.character(connectome$receptor)
  df <- data.frame(connectome %>% select(ligand,receptor))
  df$ligand <- make.unique(as.character(df$ligand))
  df$receptor <- make.unique(as.character(df$receptor))
  #df$weight <- connectome[,weight.attribute]
  temp <- connectome[,!colnames(connectome) %in% colnames(df)]
  df <- cbind(df,temp)
  
  # Squash ligands back together to single name if they are duplicates (different edges on same cell type)
  for (i in 1:length(unique(df$lig.stash))){
    temp <- subset(df,lig.stash == unique(df$lig.stash)[i])
    for (j in 1:length(unique(temp$source))){
      temp2 <- subset(temp,source == unique(temp$source)[j])
      dummy <- paste(rep(' ',j-1), collapse = '') 
      # Add number of spaces corresponding to number of unique sources
      df[rownames(temp2),]$ligand <- paste(as.character(temp2$lig.stash),dummy,sep='')
    }
    #if(length(unique(temp$source)) == 1){
    #  df[rownames(temp),]$ligand <- as.character(temp$lig.stash)
    #}
  }
  
  # Squash receptors back together to single name if they are duplicates 
  # (different edges on same cell type)
  for (i in 1:length(unique(df$rec.stash))){
    temp <- subset(df,rec.stash == unique(df$rec.stash)[i])
    for (j in 1:length(unique(temp$target))){
      temp2 <- subset(temp,target == unique(temp$target)[j])
      dummy <- paste(rep(' ',j-1),collapse = '') 
      # Add number of spaces corresponding to number of unique targets
      df[rownames(temp2),]$receptor <- paste(as.character(temp2$rec.stash),dummy,sep='')
    }
    #if(length(unique(temp$target)) == 1){
    #  df[rownames(temp),]$receptor <- as.character(temp$rec.stash)
    #}
  }
  
  # Squash ligands back together, by cell type, if they are expressed on multiple cell types
  #temp <- subset(df,ligand != lig.stash) # this is a problem
  #if (nrow(temp)>0){
  #  for (i in 1:length(unique(temp$source))){
  #    temp2 <- subset(temp,source == unique(temp$source)[i])
  #    dummy <- paste(rep(' ',i),collapse = '') # Add number of spaces corresponding to number of unique sources
  #    df[rownames(temp2),]$ligand <- paste(as.character(temp2$lig.stash),dummy,sep='')
  #  }
  #}
  # Squash receptors back together, by cell type, if they are expressed on multiple cell types
  #temp <- subset(df,receptor != rec.stash) # this is a problem
  #if (nrow(temp)>0){
  #  for (i in 1:length(unique(temp$target))){
  #    temp2 <- subset(temp,target == unique(temp$target)[i])
  #    dummy <- paste(rep(' ',i),collapse = '') # Add number of spaces corresponding to number of unique targets
  #    df[rownames(temp2),]$receptor <- paste(as.character(temp2$rec.stash),dummy,sep='')
  #  }
  #}
  
  #Establish ordering (order) so that genes are grouped nicely by celltype
  source.order <- df[order(df$source), ]
  target.order <- df[order(df$target), ]
  source.order.un <- unique(source.order[,c('ligand','source')])
  target.order.un <- unique(target.order[,c('receptor','target')])
  
  source.order$id <- 1:nrow(source.order)
  target.order$id <- 1:nrow(target.order)
  source.order.un$id <- 1:nrow(source.order.un)
  target.order.un$id <- 1:nrow(target.order.un)
  
  sector.order.un <- c(as.character(source.order.un$ligand),
                       as.character(target.order.un$receptor))
  
  # Coloring setup
  if (is.null(cols.use)){
    nodes <- as.character(unique(union(df$source,df$target)))
    cols.use <- hue_pal()(length(nodes))
    names(cols.use) <- nodes
    cols.use <- data.frame(cols.use)
    cols.use$cell <- rownames(cols.use)
  }else{
    cols.use <- data.frame(cols.use)
    cols.use$cell <- rownames(cols.use)
  }
  
  
  # Map to get ligand colorings (edges)
  map <- base::merge(source.order, cols.use, by.x = "source", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  lig.cols.edges <- as.character(map$cols.use)
  names(lig.cols.edges) <- map$ligand
  
  # Map to get receptor colorings (edges) # this does not work
  map <- base::merge(target.order, cols.use, by.x = "target", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  rec.cols.edges <- as.character(map$cols.use)
  names(rec.cols.edges) <- map$receptor
  
  # Map to get ligand colorings (sectors)
  map <- base::merge(source.order.un, cols.use, by.x = "source", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  lig.cols.sect <- as.character(map$cols.use)
  names(lig.cols.sect) <- map$ligand
  
  # Map to get receptor colorings (sectors)
  map <- base::merge(target.order.un, cols.use, by.x = "target", by.y = "cell", all = FALSE)
  map <- map[order(map$id), ]
  rec.cols.sect <- as.character(map$cols.use)
  names(rec.cols.sect) <- map$receptor
  
  # Make sector colors (grid.cols)
  sectors <- c(source.order.un$ligand, target.order.un$receptor)
  sector.cols <- c(as.character(lig.cols.sect), as.character(rec.cols.sect))
  
  counts_sectors = table(sectors)
  lig_rep = names(counts_sectors[which(counts_sectors >1)])
  if(length(lig_rep) > 0){
    cat('genes annotated as both ligand and receptor : ', lig_rep, '\n')
    #stop()
  }
  
  # Plotting
  # Decide edge order and edge color order
  if (edge.color.by.source == T){
    edge.color <- lig.cols.edges
    df.plot <- source.order
  }else{
    edge.color <- rec.cols.edges
    df.plot <- target.order
  }
  
  # Decide weight attributes and balanced vs. not
  if (weight.attribute == 'weight_norm'){
    if (balanced.edges == T){
      df.plot <- df.plot[,c('ligand','receptor','weight_norm')]
    }else{
      df.plot <- df.plot[,c('ligand','receptor','ligand.expression','recept.expression')]
    }
  }
  
  if (weight.attribute == 'weight_sc'){
    if (balanced.edges == T){
      df.plot <- df.plot[,c('ligand','receptor','weight_sc')]
    }else{
      df.plot <- df.plot[,c('ligand','receptor','ligand.scale','recept.scale')]
    }
  }
  
  #if (weight.attribute == 'score'){
  #  if (balanced.edges == T){
  #    df.plot <- df.plot[,c('ligand','receptor','score')]
  #  }else{
  #    df.plot <- df.plot[,c('ligand','receptor','ligand.norm.lfc','recept.norm.lfc')]
  #  }
  #}
  
  circos.clear()
  #circos.par(gap.degree = gap.degree)
  chordDiagram(df.plot,
               order = sector.order.un,
               col = edge.color,
               grid.col = sector.cols,
               directional = 1,
               direction.type = "arrows",
               link.arr.type = "big.arrow",
               annotationTrack = "grid",
               preAllocateTracks = 1,
               small.gap = small.gap,
               big.gap = big.gap)
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .01, sector.name, facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5),cex = lab.cex)
    #circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)
  # Make and add legend
  legend <- Legend(at = as.character(unique(union(df$source,df$target))),
                   type = "grid",
                   legend_gp = gpar(fill = as.character(cols.use[as.character(unique(union(df$source,df$target))),]$cols.use)),
                   title_position = "topleft",
                   title = "Cell Type")
  draw(legend, x = unit(20, "mm"), y = unit(20, "mm"), just = c("left", "bottom"))
  if(!is.null(title)){title(title)}
  
  p1.base <- recordPlot()
  return(p1.base)
  
  
}


########################################################
########################################################
# Section II : Differential NicheNet analysis
# 
########################################################
########################################################
run_Diff_NicheNet = function(refs = refs, 
                             timepoint_specific = TRUE,
                             include_autocrine = TRUE,
                             celltypes = NULL,
                             celltypes_BZ_specificDay = NULL, 
                             celltypes_RZ_specificDay = NULL,
                             #receivers_BZ_specificDay = list(c('CM_IS')),
                             #receivers_RZ_specificDay = list(c('CM_ROBO2')), 
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
  
  # celltypes_BZ_specificDay = celltypes_BZ_timeSpecific[[n]]; celltypes_RZ_specificDay = celltypes_RZ_timeSpecific[[n]]
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
    
    if(is.null(celltypes_BZ_specificDay)| is.null(celltypes_RZ_specificDay)){
      stop('no list of time-specific cell types BZ or RZ found \n')
    }else{
      if(!is.list(celltypes_BZ_specificDay) | length(celltypes_BZ_specificDay) == 0 |
         !is.list(celltypes_RZ_specificDay) | length(celltypes_RZ_specificDay) == 0){
        stop('lists of time-specific BZ/RZ cell types expected \n')
        
      }else{
        #cat(length(celltypes_BZ_specificDay), 'time points specified for BZ  \n')
        cat( 'receivers in BZ :', names(celltypes_BZ_specificDay), '\n\n')
        
        #cat(length(celltypes_RZ_specificDay), 'time points specified for RZ \n')
        cat('receivers in RZ :', names(celltypes_RZ_specificDay), ' \n')
        
        if(length(celltypes_BZ_specificDay) != length(celltypes_RZ_specificDay)){
          stop('not the same time points specified for BZ and RZ')
        } 
        
        # double check the celltypes defined in BZ
        for(n in 1:length(celltypes_BZ_specificDay))
        {
          celltypes = unique(c(celltypes_BZ_specificDay[[n]], names(celltypes_BZ_specificDay[n])))
          mm = match(celltypes, refs$celltypes)
          if(length(which(is.na(mm)))>0){
            stop(names(celltypes_BZ_specificDay)[n], 
                 ': some selected cell types in BZ not found in refs  -- ', celltypes[which(is.na(mm))])
          }else{
            cat(celltypes,'\n')
            cat('--all selected cell types of BZ were found in the refs \n')
          }
        }
        
        for(n in 1:length(celltypes_RZ_specificDay))
        {
          celltypes = unique(c(celltypes_RZ_specificDay[[n]], names(celltypes_RZ_specificDay)[n]))
          mm = match(celltypes, refs$celltypes)
          if(length(which(is.na(mm)))>0){
            stop(names(celltypes_RZ_specificDay)[n], 
                 ': some selected cell types in RZ not found in refs  -- ', celltypes[which(is.na(mm))])
          }else{
            cat(celltypes,'\n')
            cat('--all selected cell types of RZ were found in the refs \n')
          }
        }
        # for(n in length(receivers_BZ_specificDay))
        # {
        #   # n = 1
        #   receivers = receivers_BZ_specificDay[[n]]
        #   receivers.ctl =  receivers_RZ_specificDay[[n]]
        #   
        #   mm = match(receivers, celltypes_BZ_specificDay[[n]])
        #   if(any(is.na(mm))){stop(receivers, ' not found in the celltypes_BZ') }
        #   
        #   if(length(receivers.ctl)!= 1) stop('only one control receiver is allowed at each time point')
        #   mm.ctl = match(receivers.ctl, celltypes_RZ_specificDay[[n]])
        #   if(is.na(mm.ctl)){stop(receivers.ctl, ' not found in the celltype_RZ')}
        #   
        # }
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
  # loop over the receiver cells
  ##########################################
  for(n in 1:length(celltypes_BZ_specificDay))
  {
    # n = 1
    celltypes_BZ = celltypes_BZ_specificDay[[n]]
    celltypes_RZ = celltypes_RZ_specificDay[[n]]
    
    receivers_BZ = names(celltypes_BZ_specificDay)[n]
    receivers_RZ = names(celltypes_RZ_specificDay)[n]
    
    if(include_autocrine){
      cat('-- autocrine is considerede -- \n')
      celltypes_BZ = unique(c(celltypes_BZ, receivers_BZ))
      celltypes_RZ = unique(c(celltypes_RZ, receivers_RZ))
    }
    
    #timepoint = receivers_BZ # in case this variable still remains
    
    #cat(n, '-- ', timepoint, '\n')
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
    
    set.seed(0)
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
                            #'_timepoint.', timepoint, 
                            '.rds'))
      
      # res = readRDS(file = paste0(outDir, '/nichenet_prioritization_tables_output',
      #                             '_receiverBZ.', receivers_BZ_specificDay,
      #                             '_receiverRZ.', receivers_RZ_specificDay, 
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
                      #'_timepoint.', timepoint,
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
                      #'_timepoint.', timepoint,
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
  cat('find DE targets in receiver cell, may take some time \n')
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
  ligand_activities_targets %>% arrange(-activity) %>% filter(receiver %in% niches$BZ_niche$receiver) %>%
    filter(ligand == 'ITGA4')
  ligand_activities_targets %>% arrange(-activity) %>% filter(receiver %in% niches$RZ_niche$receiver) %>%
    filter(ligand == 'ITGA4')
  
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
    filter(receiver == niches[[2]]$receiver) %>% 
    head(10)
  
  Test_by_plotting = FALSE
  if(Test_by_plotting){
    
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
    
    ntop = 50
    cat('ntop -- ', ntop, '\n')
    
   xx =  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
      filter(sender == 'RBC') %>% 
      dplyr::group_by(ligand_receptor, niche) %>%
      filter(ligand_receptor %in% c('FGF17--FGFR1', 'ITGA4--VCAM1', 'PTDSS1--ERBB2', 'WNT7B--LRP5')) %>% 
      data.frame()
    
    ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
      dplyr::select(niche, sender, receiver, ligand, prioritization_score) %>% 
      group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% 
      ungroup() %>% distinct() %>% 
      inner_join(top_ligand_niche_df) %>% 
      filter(niche == top_niche) %>% 
      group_by(niche) %>% 
      top_n(n = ntop, prioritization_score) %>% 
      ungroup() # get the top50 ligands per niche
    
    #  select top ligand-receptor pairs for cell population of interest
    # (here, we will take the top 2 scoring receptors per prioritized ligand)
    receiver_oi = receiver
    filtered_ligands = ligand_prioritized_tbl_oi %>% 
      filter(receiver == receiver_oi) %>% 
      pull(ligand) %>% unique()
    
    prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% 
      filter(ligand %in% filtered_ligands) %>% 
      dplyr::select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
      distinct() %>% 
      inner_join(top_ligand_receptor_niche_df) %>% 
      group_by(ligand_receptor) %>%
      #group_by(ligand) %>% 
      #filter(receiver == receiver_oi) %>% 
      top_n(2, prioritization_score) %>% ungroup() 
    
    
    lfc_plot = make_ligand_receptor_lfc_plot(receiver_oi, 
                                             prioritized_tbl_oi, 
                                             prioritization_tables$prioritization_tbl_ligand_receptor, 
                                             plot_legend = FALSE, 
                                             heights = NULL, widths = NULL)
    lfc_plot
    
    
  }
  
  return(list(prioritization_tables, output))
  
  
}


make_ligand_receptor_lfc_plot_customized = function(receiver_oi, 
                                                    prioritized_tbl_oi, 
                                                    prioritization_tbl_ligand_receptor, 
                                                    plot_legend = TRUE, 
                                                    heights = NULL, 
                                                    widths = NULL)
{
  # prioritization_tbl_ligand_receptor = prioritization_tables$prioritization_tbl_ligand_receptor
  
  filtered_ligand_receptors = prioritized_tbl_oi %>% pull(ligand_receptor) %>% unique()
  
  ordered_ligand_receptors = prioritization_tbl_ligand_receptor %>% 
    filter(ligand_receptor %in% filtered_ligand_receptors) %>% 
    select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% 
    distinct() %>% 
    group_by(ligand_receptor) %>% 
    summarise(prioritization_score = max(prioritization_score)) %>% 
    inner_join(prioritization_tbl_ligand_receptor %>% 
                 select(niche, sender, ligand, receptor, ligand_receptor, prioritization_score) %>% 
                 distinct()) %>% 
    arrange(sender, prioritization_score)
  ordered_ligand_receptors_max_ligand_score = prioritization_tbl_ligand_receptor %>% 
    filter(ligand_receptor %in% filtered_ligand_receptors) %>% 
    select(niche, sender, ligand, prioritization_score) %>% 
    distinct() %>% 
    group_by(ligand) %>% 
    summarise(prioritization_score_ligand = max(prioritization_score)) %>% 
    inner_join(prioritization_tbl_ligand_receptor %>% 
                 select(niche, sender, ligand, prioritization_score) %>% 
                 distinct()) %>% 
    arrange(sender, prioritization_score_ligand) %>% distinct()
  
  ordered_ligand_receptors = ordered_ligand_receptors %>% 
    inner_join(ordered_ligand_receptors_max_ligand_score) %>% 
    arrange(sender, prioritization_score_ligand, prioritization_score)
  ordered_ligand_receptors = ordered_ligand_receptors %>% 
    mutate(ligand_receptor_ordered = factor(ligand_receptor, ordered = T, 
                                            levels = unique(ordered_ligand_receptors$ligand_receptor))) %>% 
    distinct(ligand_receptor, ligand, receptor, ligand_receptor_ordered, niche) %>% 
    dplyr::rename(niche_prior = niche)
  
  plot_data = prioritization_tbl_ligand_receptor %>% inner_join(ordered_ligand_receptors)
  
  p_lig_lfc = plot_data %>% 
    mutate(scores_test = abs(scaled_activity_normalized) * abs(ligand_score)) %>%
    #filter(receiver %in% niches[[1]]$receiver) %>%
    #filter(niche %in% names(niches)[1]) %>% 
    #ggplot(aes(sender, ligand_receptor_ordered, fill = ligand_score)) +
    ggplot(aes(sender, ligand_receptor_ordered, fill = scores_test)) +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + 
    labs(fill = "Ligand:\nmin LFC vs\nother niches") + 
    ylab("Prioritized Ligand-Receptor pairs") + xlab("Ligand LFC\n in Sender")
  
  max_lfc = max(abs(plot_data$ligand_score) %>% max(), abs(plot_data$receptor_score) %>% max())
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "OrRd"), 
                                           #%>% 
                                          #   rev(),
                                          # values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  
                                           #limits = c(-1*max_lfc, max_lfc),
                                           limits = c(0, max_lfc),
                                           )
  
  p_lig_lfc = p_lig_lfc + custom_scale_fill
  p_lig_lfc
  
  p_lig_lfc_ctl = plot_data %>% 
    filter(receiver %in% niches[[2]]$receiver) %>%
    filter(niche %in% names(niches)[2]) %>% 
    ggplot(aes(sender, ligand_receptor_ordered, fill = ligand_score)) +
    geom_tile(color = "black") +
    facet_grid(~niche, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + 
    labs(fill = "Ligand:\nmin LFC vs\nother niches")  + 
    ylab("Prioritized Ligand-Receptor pairs") + xlab("Ligand LFC\n in Sender") +
    custom_scale_fill
  
  design = "A#B"
  p_L_niches = patchwork::wrap_plots(A = p_lig_lfc, 
                                     B = p_lig_lfc_ctl + ylab(""), 
                                     nrow = 1, guides = "collect", 
                                     design = design, 
                                     widths = c(plot_data$sender %>% unique() %>% length(), 1, 
                                                plot_data$receiver %>% unique() %>% length() +0.5))
  p_L_niches
  
  p_rec_lfc = plot_data %>%
    ggplot(aes(receiver, ligand_receptor_ordered, fill = receptor_score)) +
    geom_tile(color = "black") +
    facet_grid(~receiver, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 0),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.75, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 7, color = "black", face = "bold"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid"),
      strip.background.y = element_blank(),
      strip.text.y = element_blank()
    ) + labs(fill = "Receptor:\nmin LFC vs\nother niches")  + xlab("Receptor LFC\n in Receiver")
  max_lfc = max(abs(plot_data$ligand_score) %>% max(), abs(plot_data$receptor_score) %>% max())
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p_rec_lfc = p_rec_lfc + custom_scale_fill
  p_rec_lfc
  
  design = "A#B"
  p_LR_pair = patchwork::wrap_plots(A = p_lig_lfc, B = p_rec_lfc + ylab(""), nrow = 1, guides = "collect", design = design, widths = c(plot_data$sender %>% unique() %>% length(), 1 ,plot_data$receiver %>% unique() %>% length() +0.5))
  p_LR_pair
  
  
  
}



##########################################
# post-analysis of NicheNet 
##########################################
extract_tables_from_res_Diff_NicheNet = function(outDir)
{
  saved_list = list.files(path = outDir, 
                          pattern = '*.rds', full.names = TRUE)
  
  for(n in 1:length(saved_list))
  {
    # n = 1
    cat(n, ' -- ', basename(saved_list[n]), '\n')
    res = readRDS(file = saved_list[n])
    
    prioritization_tables = res[[1]]
    output = res[[2]]
    
    file_name = gsub('.rds','', basename(saved_list[n]))
    out_Res = paste0(outDir, '/', file_name, '/')
    system(paste0('mkdir -p ', out_Res))
    
    #table_targets = output$ligand_activities_targets
    write.table(output$ligand_activities_targets, 
                file = paste0(out_Res, 'ligand_activities_targets.txt'), 
                sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    
    write.table(prioritization_tables$prioritization_tbl_ligand_receptor, 
                file = paste0(out_Res, 'prioritization_tbl_ligand_receptor.txt'), 
                sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
    
    write.table(prioritization_tables$prioritization_tbl_ligand_target, 
                file = paste0(out_Res, 'prioritization_tbl_ligand_target.txt'), 
                sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
  }
  
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


########################################################
########################################################
# Section : LR interaction visualization
# 
########################################################
########################################################
##########################################
# original code from https://github.com/ShellyCoder/cellcall/blob/master/R/ViewInterCircos.R
##########################################
#' plot circle graph with communication profile
#' @param object a Cellwave objects
#' @param font the size of font
#' @param cellColor a color dataframe, rownames is cell type, value is color
#' @param lrColor a color vector denotes the color of ligand and receptor, containing two elements, default is c('#D92E27', "#35C6F4")
#' @param order.vector default is null, a celltype vector with the order you want in the circle graph
#' @param trackhight1 Height of the outer track
#' @param trackhight2 Height of the inner track
#' @param linkcolor.from.sender logical value, whether the color of line correspond with color of sender cell
#' @param linkcolor one color you want link to be, only if parameter linkcolor.from.sender=FALSE
#' @param arr.type Type of the arrows, default value is big.arrow There is an additional option triangle
#' @param arr.length Length of the arrows, measured in 'cm'. If arr.type is set to big.arrow, the value is percent to the radius of the unit circle.
#' @param DIY logical value, if TRUE, the parameter object should be a dataframe, and set slot="expr_l_r_log2_scale". otherwise object should be a Cellwave objects.
#' @param slot plot the graph with the data of specific slot
#' @param gap.degree between two neighbour sectors. It can be a single value or a vector. If it is a vector, the first value corresponds to the gap after the first sector.
#' @param track.margin2 affect current track
#' @importFrom grid pushViewport unit upViewport viewport gpar
#' @importFrom graphics plot.new par
#' @importFrom gridBase gridOMI
#' @importFrom circlize circos.clear circos.par circos.initialize circos.trackPlotRegion circos.link get.cell.meta.data highlight.sector colorRamp2 rand_color
#' @importFrom stringr str_split
#' @importFrom magrittr %>% set_colnames
#' @importFrom dplyr filter
#' @importFrom ComplexHeatmap draw Legend packLegend
#' @export

ViewInterCircos_customized <- function(object, font = 2, cellColor ,lrColor = NULL, order.vector = NULL,
                            trackhight1 = 0.05, linkcolor.from.sender = TRUE, linkcolor = NULL,
                            arr.type = "big.arrow",arr.length = 0.04, DIY = FALSE, 
                            gap.degree = NULL,
                            trackhight2 = 0.032, track.margin2 = c(0.01,0.12),
                            slot="expr_l_r_log2_scale")
{
  
  require(ggplot2)
  library(tidyverse)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(gridBase)
  library(gridExtra)
  require(ComplexHeatmap)
  
  # object = mt; slot = "expr_l_r_log2_scale"; order.vector=c('ST', "SSC", "SPGing", "SPGed")
  plot.new()
  circle_size = unit(1, "snpc") # snpc unit gives you a square region
  pushViewport(viewport(x = 0, y = 0.5, width = circle_size, 
                        height = circle_size,
                        just = c("left", "center")))
  par(omi = gridOMI(), new = TRUE)
  
  circos.clear()  
  
  # cell.padding: the padding between sector and other sector
  if(is.null(gap.degree)){
    circos.par(canvas.xlim =c(-1.1,1.1),canvas.ylim = c(-1.1,1.1),cell.padding = c(0.01,0,0.01,0))
  }else{
    circos.par(canvas.xlim =c(-1.1,1.1),canvas.ylim = c(-1.1,1.1),cell.padding = c(0.01,0,0.01,0), 
               gap.degree=gap.degree)
  }
  
  # library(stringr)
  if(DIY){
    a <- colSums(object)
  }else{
    a <- colSums(object@data[[slot]])
  }
  
  if(sum(a>0)==0){
    stop("There is no predicted L-R in @data$expr_l_r_log2_scale, looser paramter might be of help.")  
  }
  
  b <- stringr::str_split(names(a), "-", simplify = T)
  c <- data.frame(b, stringsAsFactors = FALSE)
  c$x3 <- as.numeric(a)
  test <- c
  colnames(test) <- c( "cell_from", "cell_to", "n")
  
  test$id1 <- paste("sender", test$cell_from, sep = "_", 1:nrow(test))
  test$id2 <- paste("recevier", test$cell_to, sep = "_", 1:nrow(test))
  
  a_tmp <- test[,c(1,4)]
  colnames(a_tmp) <- c('celltype', 'id')
  a_tmp$arrowType <-  "sender"
  
  b_tmp <- test[,c(2,5)]
  colnames(b_tmp) <- c('celltype', 'id')
  b_tmp$arrowType <-  "recevier"
  
  ab_tmp <- rbind(a_tmp,b_tmp)
  
  ab_tmp <- data.frame(ab_tmp, stringsAsFactors = FALSE)
  ab_tmp <- ab_tmp[order(ab_tmp$celltype, ab_tmp$id, decreasing = TRUE),]
  
  sector_id <- ab_tmp$id
  fa = ab_tmp$id
  
  if(is.null(order.vector)){
    fa = factor(fa,levels = fa)
  }else{
    fa.df <- str_split(fa, "_", simplify = T) %>% as.data.frame() %>% 
      magrittr::set_colnames(c('sender_or_receiver', 'clltype', "index_number"))
    my.levels <- do.call(rbind,lapply(order.vector, function(x){
      fa.df[which(fa.df$clltype==x),]
    })) %>% apply(1, function(x){
      paste(x[1],x[2],x[3],sep='_')
    }) %>% unlist %>% as.character()
    
    fa = factor(fa,levels = my.levels)
  }
  
  circos.initialize(factors = fa, xlim = c(0,1)) # 
  
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = trackhight1,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector.index = get.cell.meta.data('sector.index')
      xlim = get.cell.meta.data('xlim')
      ylim = get.cell.meta.data('ylim')
      # print(sector.index)
      # print(xlim)
      # print(ylim)
    }
  )
  
  if(is.null(cellColor)){
    cell_color <- data.frame(color = rand_color(length(unique(ab_tmp$celltype)), luminosity = "light"), stringsAsFactors = FALSE)
    rownames(cell_color) <- unique(ab_tmp$celltype)
  }else{
    cell_color <- cellColor
  }
  
  cell_type <- unique(ab_tmp$celltype)
  for(i in 1:length(cell_type)){
    myCell_Type <-  cell_type[i]
    mySector_id <- as.character(ab_tmp[ab_tmp$celltype==myCell_Type,'id'])
    myColor <- as.character(cell_color[myCell_Type,1])
    highlight.sector(mySector_id, track.index = 1,
                     text = myCell_Type, text.vjust = -1,niceFacing = T, font = font, col = myColor)
    
  }
  
  circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = trackhight2,
    bg.border = NA,
    track.margin = track.margin2,
    panel.fun = function(x, y) {
      sector.index = get.cell.meta.data('sector.index')
      xlim = get.cell.meta.data('xlim')
      ylim = get.cell.meta.data('ylim')
    }
  )
  
  if(is.null(lrColor)){
    ligand_receptor <- c('#D92E27', "#35C6F4")  #color of sender and receiver
  }else{
    ligand_receptor <- lrColor
  }
  
  for(i in 1:length(cell_type)){
    myCell_Type <-  cell_type[i]
    mySector_id_tmp <- ab_tmp[ab_tmp$celltype==myCell_Type,c('arrowType', 'id')]
    my_ligand <- dplyr::filter(mySector_id_tmp, arrowType == 'sender')
    my_receptor <- dplyr::filter(mySector_id_tmp, arrowType == 'recevier')
    
    if(length(as.character(my_ligand[,'id']))>0){
      highlight.sector(as.character(my_ligand[,'id']), track.index = 2,
                       text = '', niceFacing = F,col = ligand_receptor[1],text.col = 'white')
      
    }
    
    if(length(as.character(my_receptor[,'id']))>0){
      highlight.sector(as.character(my_receptor[,'id']), track.index = 2,
                       text = '', niceFacing = F,col = ligand_receptor[2],text.col = 'white')
      
    }
    
  }
  
  
  test$weighted_n <- (test$n-min(test$n))/(max(test$n)-min(test$n))
  test <- test[order(test$weighted_n, decreasing = F),]
  # print(test$weighted_n)
  # test$weighted_n <- as.numeric(test$n/sum(test$n))
  
  min_n <- min(test$weighted_n)
  max_n <- max(test$weighted_n)
  mean_n <- mean(min_n, max_n)
  
  if(linkcolor.from.sender){
    
    color.function <- apply(cell_color, 1, function(x){
      col_fun = colorRamp2(c(0, 1), c("#FFFFFF", x))
    })
    
    for(i in 1:nrow(test)){
      my_Line_color = as.character(cell_color[test[i,1],1])
      print(test$weighted_n[i])
      circos.link(sector.index1 = test[i,'id1'],
                  point1 = c(0,1),
                  sector.index2 = test[i,'id2'],
                  point2 = c(0,1),
                  directional = 1,
                  arr.type = arr.type,
                  # arr.width = 0,
                  arr.length = arr.length,
                  col=color.function[[test[i,'cell_from']]](test$weighted_n[i])
      )
      
    }
    
    upViewport()
    list.obj <- list()
    # discrete
    lgd_points = Legend(at = c('ligand', "receptor"), type = "points", gap = unit(2, "mm"),
                        legend_gp = gpar(col = lrColor), title_position = "topleft",
                        title = "")
    # discrete
    lgd_lines = Legend(at = rownames(cell_color), type = "lines", gap = unit(2, "mm"),
                       legend_gp = gpar(col = cell_color$color, lwd = 2), title_position = "topleft",
                       title = "Cell type")
    
    list.obj <- c(list.obj, list(lgd_lines))
    list.obj <- c(list.obj, list(lgd_points))
    # # continuous
    for (f in color.function) {
      lgd_links = Legend(at = c(0, 0.5, 1), col_fun = f, gap = unit(1, "mm"),direction="horizontal",
                         title_position = "topleft", title = "")
      list.obj <- c(list.obj, lgd_links)
    }
    
    lgd_list_vertical = packLegend(list = list.obj,
                                   gap = unit(2, "mm") # ligend distance
    )
    
    draw(lgd_list_vertical, x = circle_size, just = "left")
  }else{
    
    if(is.null(linkcolor)){
      col_fun = colorRamp2(c(0, 1), linkcolor)
    }else{
      col_fun = colorRamp2(c(0, 1), c("#FFFFFF", "#f349eb"))
    }
    
    for(i in 1:nrow(test)){
      my_Line_color = as.character(cell_color[test[i,1],1])
      print(test$weighted_n[i])
      circos.link(sector.index1 = test[i,'id1'],
                  point1 = c(0,1),
                  sector.index2 = test[i,'id2'],
                  point2 = c(0,1),
                  directional = 1,
                  arr.type = arr.type,
                  # arr.width = 0,
                  arr.length = arr.length,
                  col=col_fun(test$weighted_n[i])
      )
      
    }
    upViewport()
    
    # discrete
    lgd_points = Legend(at = c('ligand', "receptor"), type = "points", gap = unit(2, "mm"),
                        legend_gp = gpar(col = lrColor), title_position = "topleft",
                        title = "")
    # discrete
    lgd_lines = Legend(at = rownames(cell_color), type = "lines", gap = unit(2, "mm"),
                       legend_gp = gpar(col = cell_color$color, lwd = 2), title_position = "topleft",
                       background="#FFFFFF", title = "Cell type")
    # # continuous
    lgd_links = Legend(at = c(0, 0.5, 1), col_fun = col_fun, gap = unit(2, "mm"), direction="vertical",
                       background="#FFFFFF", title_position = "topleft", title = "Adjusted Score")
    
    lgd_list_vertical = packLegend(lgd_lines, lgd_points, lgd_links,
                                   gap = unit(4, "mm") # ligend distance
    )
    
    draw(lgd_list_vertical, x = circle_size, just = "left")
  }
  
  
}

##########################################
# plot code from celltalker
# https://github.com/arc85/celltalker/blob/master/R/circos_plot.R
##########################################
circos_plot_customized = function(ligand_receptor_frame,
                        cell_group_colors,
                        ligand_color="blue",
                        receptor_color="red",
                        cex_outer=0.5,
                        cex_inner=0.4,
                        link.lwd=0.5, 
                        arr.length=0.2, 
                        arr.width=(3*0.1)/2) {
  
  # ligand_receptor_frame=top_stats_xx; cell_group_colors=colors_use;
  # Bind variables
  cell_type1 <- lig <- cell_type2 <- rec <- classes <- ranges <-
    max_range <- to_class <- to_rec <- lig_rec <- ordered_lig_rec <- type <-
    lig.rec <- to.class <- to.rec <- ordered.lig.rec <- NULL
  
  # Reformat data
  part1 <- ligand_receptor_frame %>%
    mutate(lig=sapply(strsplit(interaction,split="_"),function(x) x[[1]])) %>%
    mutate(rec=sapply(strsplit(interaction,split="_"),function(x) x[[2]])) %>%
    select(cell_type1,lig) %>%
    distinct() %>%  
    mutate(type="lig")
  part2 <- ligand_receptor_frame %>%
    mutate(lig=sapply(strsplit(interaction,split="_"),function(x) x[[1]])) %>%
    mutate(rec=sapply(strsplit(interaction,split="_"),function(x) x[[2]])) %>%
    select(cell_type2,rec) %>%
    distinct() %>%  
    mutate(type="rec")
  colnames(part1) <- colnames(part2) <- c("classes","lig.rec","type")
  
  part12 <- rbind(part1,part2) %>%
    group_by(classes) %>%
    group_split()
  
  part12 <- lapply(part12,function(x) {
    
    x <- x %>%
      mutate(ordered.lig.rec=paste(type,lig.rec,sep="_")) %>%
      mutate(ranges=as.numeric(as.factor(ordered.lig.rec))) %>%
      select(-ordered.lig.rec)
    
  })
  
  part12 <- do.call(rbind,part12)
  
  to.join <- ligand_receptor_frame %>%
    mutate(lig=sapply(strsplit(interaction,split="_"),function(x) x[[1]])) %>%
    mutate(rec=sapply(strsplit(interaction,split="_"),function(x) x[[2]])) %>%
    select(cell_type1,cell_type2,lig,rec)
  
  colnames(to.join)[1:3] <- c("classes","to.class","lig.rec")
  
  part3 <- part12
  
  joined <- left_join(part3,to.join,by=c("classes","lig.rec"))
  joined$to.rec <- NA
  
  for(i in 1:nrow(joined)) 
  {
    sub.group <- joined[i,]
    sub.joined <- joined %>% filter(classes==sub.group$to.class)
    joined$to.rec[i] <- sub.joined[match(sub.group$rec,sub.joined$lig.rec),"ranges"] %>% pull()
    
  }
  
  final.construct <- joined
  
  # Repair single class
  single.class <- final.construct %>%
    group_by(classes) %>%
    summarize(max_range=max(ranges)) %>%
    filter(max_range==1) %>%
    pull(classes)
  
  if (!length(single.class)==0) {
    
    for (i in 1:length(single.class)) 
    {
      row.add <- final.construct[final.construct$classes==single.class[i],][1,]
      row.add$ranges <- 2
      
      final.construct <- rbind(final.construct,row.add)
      final.construct <- final.construct %>%
        arrange(classes)
    }
    
  }
  
  final.construct <- final.construct %>%
    arrange(classes,ranges)
  
  circos.clear()
  circos.par(gap.degree= 2, track.margin=c(0,0.2))
  circos.initialize(factors=final.construct$classes, x=final.construct$ranges)
  
  suppressMessages({
    circos.track(ylim = c(0, 1), track.height=0.05, panel.fun = function(x, y) {
      circos.rect(CELL_META$cell.xlim[1], 
                  CELL_META$cell.ylim[1],
                  CELL_META$cell.xlim[2],
                  CELL_META$cell.ylim[2],
                  col=cell_group_colors[CELL_META$sector.numeric.index])
      circos.text(CELL_META$xcenter, 
                  y = 2.5, 
                  CELL_META$sector.index,
                  #facing = "downward",
                  facing = "clockwise", 
                  niceFacing = TRUE,
                  cex=cex_outer,
                  adj = c(0.0, 0))
    })
  })
  
  ## Build interior track with ligand/receptors colors and gene labels
  circos.track(ylim = c(0, 1), track.height= 0.1, bg.border="white")
  
  # Define multiplers for each sector
  final.construct2 <- final.construct %>%
    select(classes,lig.rec,ranges,type) %>%
    distinct() %>%
    arrange(classes,ranges)
  ref.tab <- unname(table(final.construct2$classes))
  sec.multi <- (ref.tab-1)/ref.tab
  names(sec.multi) <- names(table(final.construct2$classes))
  
  # Loop to construct all sectors
  # Ligands first
  # Split into list of sectors
  int.types.list <- final.construct2 %>%
    group_split(classes)
  
  names(int.types.list) <- sapply(int.types.list,function(x) x$classes[1])
  int.types.list.multi <- int.types.list.individ <- list("NA")
  int.types.list.multi <- int.types.list[!names(int.types.list) %in% single.class]
  int.types.list.individ <- int.types.list[names(int.types.list) %in% single.class]
  
  for (i in 1:length(int.types.list.multi)) 
  {
    for (a in 1:nrow(int.types.list.multi[[i]])) 
    {
      if (a==1) {
        sec.multi.use <- sec.multi[names(sec.multi)==int.types.list.multi[[i]]$classes[1]]
        
        suppressMessages({
          circos.rect(2, 0, 
                      2+sec.multi.use*a, 
                      1,
                      sector.index=int.types.list.multi[[i]]$classes[a],
                      col=ifelse(int.types.list.multi[[i]]$type[a]=="lig",
                                 ligand_color,
                                 receptor_color),
                      track.index = 2)
          circos.text(1+sec.multi.use*a/4, 2,
                      sector.index=int.types.list.multi[[i]]$classes[a],
                      labels=int.types.list.multi[[i]]$lig.rec[a],
                      track.index = 2,
                      facing = "clockwise", 
                      niceFacing = TRUE,
                      cex=cex_inner)
        })
        
      } else {
        sec.multi.use <- sec.multi[names(sec.multi)==int.types.list.multi[[i]]$classes[1]]
        suppressMessages({
          circos.rect(2+sec.multi.use*(a-1), 0, 
                      2+sec.multi.use*a,
                      1,
                      sector.index=int.types.list.multi[[i]]$classes[a],
                      col=ifelse(int.types.list.multi[[i]]$type[a]=="lig",ligand_color,receptor_color),
                      track.index = 2)
          circos.text(1+sec.multi.use*a-sec.multi.use/4, 2,
                      sector.index=int.types.list.multi[[i]]$classes[a],
                      labels=int.types.list.multi[[i]]$lig.rec[a],
                      track.index = 2,
                      facing = "clockwise", niceFacing = TRUE,
                      cex=cex_inner)
        })
        
      }
    }
  }
  
  if (length(int.types.list.individ)>0) {
    
    for (i in 1:length(int.types.list.individ)) {
      
      circos.rect(1,0,2,1,sector.index=int.types.list.individ[[i]]$classes[1],
                  col=ifelse(int.types.list.individ[[i]]$type[1]=="lig",ligand_color,receptor_color),
                  track.index = 2)
      circos.text(1.5,4,sector.index=int.types.list.individ[[i]]$classes[1],
                  labels=int.types.list.individ[[i]]$lig.rec[1],
                  track.index = 2,
                  facing="downward",cex=cex_inner)
      
    }
    
  }
  
  
  ## Draw links
  final.construct3 <- joined %>%
    select(classes,lig.rec,ranges,to.class,to.rec) %>%
    distinct()
  
  split.construct <- final.construct3 %>%
    split(.$classes)
  
  final.construct3 <- lapply(split.construct,function(x) {
    class.length <- length(unique(x$ranges))
    if (class.length==1) {
      x[,"ranges"] <- 1.5
      x
    } else {
      x
    }
  }) %>%
    do.call(rbind,.)
  
  int.types.list <- final.construct3 %>%
    group_split(classes)
  
  names(int.types.list) <- sapply(int.types.list,function(x) x$classes[1])
  
  
  for (i in 1:length(int.types.list)) 
  {
    for (a in 1:nrow(int.types.list[[i]])) 
    {
      target <- which(!is.na(match(names(int.types.list),int.types.list[[i]]$to.class[[a]])))
      
      if(length(target)==0) {
        
      }else{
        if (!int.types.list[[i]]$to.class[[a]] %in% single.class) {
          circos.link(int.types.list[[i]]$classes[a], 
                      1+sec.multi[i]*int.types.list[[i]]$ranges[a]-sec.multi[i]/2,
                      int.types.list[[i]]$to.class[[a]], 
                      1+sec.multi[target]*int.types.list[[i]]$to.rec[a]-sec.multi[target]/2,
                      0.43, 0.43, directional=1, 
                      lwd=link.lwd, 
                      arr.length=arr.length, 
                      arr.width=arr.width)
          
        } else {
          circos.link(int.types.list[[i]]$classes[a], 
                      1+sec.multi[i]*int.types.list[[i]]$ranges[a]-sec.multi[i]/2,
                      int.types.list[[i]]$to.class[[a]], 
                      1.5,
                      0.43, 0.43, 
                      directional=1, 
                      lwd=link.lwd, 
                      arr.length=arr.length, 
                      arr.width=arr.width)
        }
        
        
      }
      
    }
    
  }
  
}

