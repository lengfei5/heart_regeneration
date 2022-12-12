
rm(list = ls())

species = 'axolotl'
version.analysis = '_R12830_resequenced_20220308'
dataDir = '../R12830_visium_reseqenced/nf_out'
resDir = paste0("../results/visium_axolotl", version.analysis)
RdataDir = paste0('../results/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)


library(pryr) # monitor the memory usage
require(ggplot2)
library(nichenetr)
library(Seurat) # please update to Seurat V4
library(tidyverse)
library(circlize)
library(RColorBrewer)
require(scran)
require(scater)
library(nichenetr)
library(Seurat)
library(tidyverse)
library(knitr)
library(tidyr)
library(reshape2)
library(gridExtra)
library(ComplexHeatmap)
library(clusterProfiler)
library(plotly)
library(ggalluvial)
library(stringr)

options(stringsAsFactors = F)

source('functions_scRNAseq.R')
source('functions_Visium.R')
dataPath_nichenet = '../data/NicheNet/'

mem_used()

##########################################
# load processed scRNA-seq and visium data
##########################################
# load ST data with additional region annotations
load(file = paste0(RdataDir, 'seuratObject_design_variableGenes_umap.clustered_manualSegmentation', 
                   species, '.Rdata')) # visium data

table(st$segmentation, st$condition)

## snRNA-seq reference  
refs = readRDS(file = paste0(RdataDir, 'RCTD_refs_subtypes_final_20221117.rds'))
refs$subtypes = refs$celltype_toUse # clean the special symbols
refs$celltypes = refs$celltype_toUse
#refs_file = '/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_annotated_no_doublets_20221004_2.rds'
#refs = readRDS(file = refs_file)
table(refs$subtypes)

# subtype time-specificity 
condition.specific_celltypes = readRDS(paste0(RdataDir, 'RCTD_refs_condition_specificity.rds'))

celltypes_timeSpecific = list(day1 = c('EC', 'EC_NOS3', 'EC_IS_IARS1', 'FB_IS_TFPI2', 'Mo.Macs_SNX22', 'Neu_IL1R1', 
                                       'CM_IS', "RBC"),
                              day4 = c('EC_IS_LOX', 'EC_IS_Prol', 'Mo.Macs_SNX22', 'Neu_DYSF', 'CM_IS', 
                                       'CM_Prol_IS', 'RBC'),
                              day7 = c('EC_IS_LOX', 'EC_IS_Prol', 'Mo.Macs_FAXDC2', 'Neu_DYSF', 'Neu_IL1R1', 'CM_IS', 
                                       'CM_Prol_IS', 'RBC'),
                              day14 = c('EC_IS_LOX', 'EC_IS_Prol', 'FB_PKD1', 'Neu_DYSF', 'CM_IS', 'Megakeryocytes', 
                                        'RBC')
)

timepoint_specific = TRUE

# test first one receiver with its control 
receiver_cells = 'CM_IS'
control_cells = "CM_ven_Robo2"
#refs = readRDS(file = refs_file)
table(refs$celltypes)

########################################################
########################################################
# Section : test default NicheNet with Seurat object 
# 
########################################################
########################################################
outDir = paste0(resDir, '/Ligand_Receptor_analysis/NicheNet_v3_timeSpecific_receiverCell.CM_IS')
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

ligands = lr_network %>% pull(ligand) %>% unique()
receptors = lr_network %>% pull(receptor) %>% unique()

# ## weighted integrated networks 
weighted_networks = readRDS(paste0(dataPath_nichenet,  "weighted_networks.rds"))
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

weighted_networks_lr = weighted_networks$lr_sig %>% 
  dplyr::rename(ligand = from, receptor = to) %>%
  inner_join(lr_network %>% 
               distinct(ligand,receptor), by = c("ligand","receptor"))

##########################################
# we directly consider time-specific 
# loop over the cell types for each time points
# here we are using the basic function predict_ligand_activities in NicheNet
# Shoval is also using the same one, more flexible to adapt for data 
# original code from https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_steps.md
##########################################
for(n in 1:length(celltypes_timeSpecific)) # loop over each time point
{
  # n = 1
  celltypes = celltypes_timeSpecific[[n]]
  timepoint = names(celltypes_timeSpecific)[n]
  cat(n, '-- ', timepoint, '\n')
  cat(celltypes, '\n')
  
  
  celltypes_sel = unique(c(celltypes, receiver_cells, control_cells)) 
  
  subref = subset(refs, cells = colnames(refs)[!is.na(match(refs$celltypes, celltypes_sel))])
  subref$celltypes = droplevels(as.factor(subref$celltypes))
  table(subref$celltypes)
  
  #subref = subset(x = subref, downsample = 1000)
  cat('celltype to consider -- ', names(table(subref$celltypes)), '\n')
  Idents(subref) = subref$celltypes
  
  # process and clean the gene names
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
    
    #saveRDS(geneDup, paste0(outDir, '/geneSymbol_duplication_inLRanalysis.rds'))
    
  }
  
  # step 1:  Define a “sender/niche” cell population and a “receiver/target” cell population 
  # present in your expression data and determine which genes are expressed in both populations
  if(is.na(receiver_cells)){ # loop over all cell type candidates
    receiver_cells = celltypes
  }
  
  # to save the result for each time point
  expressed_genes = c()
  res = tibble(sample=character(), sender=character(), receiver=character(),  
               test_ligand=character(), auroc=double(), aupr=double(), 
               pearson = double(),  rank=integer())
  
  for(receiver in receiver_cells) # loop over receiver cells 
  {
    # specify receiver
    # receiver = receiver_cells[1]
    cat('-- receiver : ', receiver, '-- \n')
    expressed_genes_receiver = get_expressed_genes(receiver, subref, pct = 0.10)
    background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
    
    # step 2: Define a gene set of interest: these are the genes in the “receiver/target” cell population 
    # that are potentially affected by ligands expressed by interacting cells 
    # (e.g. genes differentially expressed upon cell-cell interaction)
    DE_table_receiver = FindMarkers(object = subref, 
                                    ident.1 = receiver, 
                                    ident.2 = control_cells, 
                                    min.pct = 0.10) %>% rownames_to_column("gene")
    
    geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
    geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
    
    ## sender
    sender_celltypes = celltypes
    #loop over the same receiver for sender cells, i.e. incl. self-activation (shoval's example)
    for (sender in sender_celltypes){
      # sender = sender_celltypes[1]
      cat('-- sender : ', sender, '-- \n')
      expressed_genes_sender = get_expressed_genes(sender, subref, pct=0.1)
      
      # Step 3: Define a set of potential ligands
      # these are ligands that are expressed by the “sender/niche” cell population and 
      # bind a (putative) receptor expressed by the “receiver/target” population
      # Get potential ligands in our datasets
      expressed_ligands = intersect(ligands, expressed_genes_sender)
      expressed_receptors = intersect(receptors,expressed_genes_receiver)
            
      potential_ligands = lr_network %>% filter(ligand %in% expressed_ligands & receptor %in% expressed_receptors) %>%
        pull(ligand) %>% unique()
      
      # Step 4: Perform NicheNet’s ligand activity analysis on the gene set of interest
      # rank the potential ligands based on the presence of their target genes in the gene set of interest 
      # (compared to the background set of genes)     
      # the main function used : predict_ligand_activities
      ligand_activities = predict_ligand_activities(geneset = geneset_oi,
                                                    background_expressed_genes = background_expressed_genes,
                                                    ligand_target_matrix = ligand_target_matrix, 
                                                    potential_ligands = potential_ligands) 
      ligand_activities = ligand_activities %>% 
        arrange(-pearson) %>% dplyr::mutate(sender=sender, receiver=receiver, 
                                            rank = rank(dplyr::desc(pearson)))
      
      # pearson correlation between the target expresson and prediciton as ligand activity scores 
      res = rbind(res, ligand_activities) 
    }
    
    # list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seuratObj, 0.10) # lapply to get the expressed genes of every sender cell type separately here
    # expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
  }
  
  # processing after the ligand activity prediction
  res = res %>% dplyr::mutate(interaction=paste0(sender,"-",receiver))
  #res = res %>% filter(test_ligand %in% growth.factors) %>% dplyr::mutate(interaction=paste0(sender,"-",receiver))
  
  res %>% DT::datatable()
  expressed_genes = unique(expressed_genes)
  
  ##########################################
  # ligand visualization
  ##########################################
  res = res %>% dplyr::mutate(interaction=paste0(sender,"-",receiver))
  
  # Plot ligand activities scores in general
  p1 = res %>% ggplot(aes(x=pearson)) + 
    geom_histogram(color="black", fill="darkorange") + 
    geom_vline(aes(xintercept=0.10), color="red", linetype="dashed", size=1) +  
    labs(x="ligand activity (Pearson Correlation)", y = "# ligands") +  
    theme_classic() + ggtitle("Ligands activity scores")
  plot(p1)
  
  ggsave(filename = paste0(outDir, '/ligand_pearson_correlation_', timepoint, '.pdf'), 
         width = 10, height = 8)
  
  # Plot ligand activities scores per interaction
  p2 = res %>% ggplot(aes(x=pearson)) + 
    geom_histogram(color="black", fill="darkorange") + 
    geom_vline(aes(xintercept=0.1), color="red", linetype="dashed", size=1) +  
    labs(x="ligand activity (Pearson Correlation)", y = "# ligands") +  
    theme_classic() + 
    ggtitle("Ligands activity scores") + 
    facet_wrap(.~interaction, )
  plot(p2)
  ggsave(filename = paste0(outDir, '/ligand_pearson_correlation_individualInteraction_', timepoint, '.pdf'), 
         width = 12, height = 10)
  
  # Visualize ligand expression along different days
  ligand.expression.df = AverageExpression(subref, features =unique(res$test_ligand))[['RNA']]
  pseudo_signal = ceiling(log2(range(ligand.expression.df[ligand.expression.df>0]))[1])-1
  ligand.expression.df = log2(2^pseudo_signal +ligand.expression.df)[sort(rownames(ligand.expression.df)), 
                                                                     sort(colnames(ligand.expression.df))]
  write.csv2(inner_join(res, ligand.expression.df %>% as.data.frame %>% rownames_to_column('test_ligand')), 
             file = paste0(outDir, '/predicted_ligand_activity_avExpr_for_receivers_senders_', timepoint, '.csv'),
             row.names = FALSE)
  
  #p2 = Heatmap(ligand.expression.df, cluster_columns = F, column_title = "Ligand Expression", name="Expression (log2)", row_order = row_order(p1))
  
  # Visualization top ligand  per interaction
  p1 = res %>% top_n(500, pearson) %>% transmute(test_ligand=test_ligand, x=interaction, pearson=pearson) %>% 
    spread(x, pearson)  %>% 
    replace(is.na(.), 0) %>% 
    ungroup() %>% 
    column_to_rownames('test_ligand')
  p2 = ligand.expression.df[match(rownames(p1), rownames(ligand.expression.df)), ]
  
  #ha.col = list(Interaction=setNames(c(brewer.pal(9, 'YlOrRd')), colnames(p1)))
  ha.col = c(brewer.pal(9, 'YlOrRd'))[1:length(colnames(p1))]
  names(ha.col) = colnames(p1)
  ha.col = as.list(ha.col)
  ha = HeatmapAnnotation(Interaction=colnames(p1), annotation_name_side = 'left', col=ha.col)
  
  p3 = p1 %>% Heatmap(
    column_title = "Average Ligand pearson correlation\n(NAs replaced with zeros)", 
    cluster_columns = F, show_row_names = F, show_column_names = F, 
    top_annotation = ha, col=brewer.pal(4, 'Oranges'), 
    name="Ligand Score(pearson)")
  
  ha = colnames(p2)
  ha = HeatmapAnnotation(Cluster=ha, col=ha.col)
  p4 = Heatmap(p2, 
               cluster_columns = F, show_column_names = F,
               column_title = "Ligand Expression", name="log10(Expression)", 
               row_order = row_order(p3), 
               top_annotation = ha)
  
  pdfname = paste0(outDir, '/heatmap_top_ligand_pearson_correlation_individualInteraction_', timepoint, '.pdf')
  pdf(pdfname, width=16, height = 25)
  
  pp = p3 + p4
  plot(pp)
  dev.off()
  
}
  
  
########################################################
########################################################
# Section : test code not used currently
# 
########################################################
########################################################


  best_upstream_ligands = ligand_activities %>% 
    top_n(30, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
  DotPlot(subref, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
  # Now, we want to rank the ligands based on their ligand activity. 
  # In our validation study, we showed that the pearson correlation coefficient (PCC) 
  # between a ligand’s target predictions and the observed transcriptional response was 
  # the most informative measure to define ligand activity. 
  # Therefore, we will rank the ligands based on their pearson correlation coefficient.
  
  # show histogram of ligand activity scores to check the selected cutoff 
  # by looking at the distribution of the ligand activity values. 
  # Here, we show the ligand activity histogram 
  #(the score for the 20th ligand is indicated via the dashed line).
  p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
    geom_histogram(color="black", fill="darkorange")  + 
    # geom_density(alpha=.1, fill="orange") +
    geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, pearson) %>% pull(pearson))), 
               color="red", linetype="dashed", size=1) + 
    labs(x="ligand activity (PCC)", y = "# ligands") +
    theme_classic()
  p_hist_lig_activity
  
  
  # Step 5: Infer target genes of top-ranked ligands and visualize in a heatmap
  active_ligand_target_links_df = best_upstream_ligands %>% 
    lapply(get_weighted_ligand_target_links, 
           geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% 
    bind_rows() %>% drop_na()
  
  active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df,
                                                                   ligand_target_matrix = ligand_target_matrix, 
                                                                   cutoff = 0.33)
  
  order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% 
    rev() %>% make.names()
  order_targets = active_ligand_target_links_df$target %>% unique() %>% 
    intersect(rownames(active_ligand_target_links)) %>% make.names()
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() 
  # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() 
  # make.names() for heatmap visualization of genes like H2-T23
  
  vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  
  p_ligand_target_network = vis_ligand_target %>% 
    make_heatmap_ggplot("Prioritized ligands","Predicted target genes", 
                        color = "purple",legend_position = "top", 
                        x_axis_position = "top",
                        legend_title = "Regulatory potential")  + 
    theme(axis.text.x = element_text(face = "italic")) + 
    scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
  p_ligand_target_network
  
  
  # Follow-up analysis 1: Ligand-receptor network inference for top-ranked ligands
  # get the ligand-receptor network of the top-ranked ligands
  lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% 
    distinct(from,to)
  best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
  
  lr_network_top_df_large = weighted_networks_lr %>% 
    filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  
  lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
  lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% 
    magrittr::set_rownames(lr_network_top_df$to)
  
  dist_receptors = dist(lr_network_top_matrix, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
  
  vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
  p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% 
    make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", 
                        x_axis_position = "top",legend_title = "Prior interaction potential")
  p_ligand_receptor_network
  
  # Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor interactions documented 
  # in literature and publicly available databases
  lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
  receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()
  
  lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% 
    inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
  lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% 
    inner_join(lr_network_top_df_large, by = c("from","to"))
  
  lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
  lr_network_top_matrix_strict = lr_network_top_df_strict %>% 
    select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)
  
  dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
  hclust_receptors = hclust(dist_receptors, method = "ward.D2")
  order_receptors = hclust_receptors$labels[hclust_receptors$order]
  
  dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
  hclust_ligands = hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
  
  order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
  order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))
  
  vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
  rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
  colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
  
  p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% 
    make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred",
                        x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
  p_ligand_receptor_network_strict
  
  
  
  # Add log fold change information of ligands from sender cells
  # DE analysis for each sender cell type
  # this uses a new nichenetr function - reinstall nichenetr if necessary!
  DE_table_all = Idents(subref) %>% levels() %>% 
    intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = subref, 
                                           condition_colname = "celltypes", 
                                           condition_oi = condition_oi, 
                                           condition_reference = condition_reference, 
                                           expression_pct = 0.10, celltype_col = NULL) %>% reduce(full_join) 
  # use this if cell type labels are the identities of your Seurat object -- 
  # if not: indicate the celltype_col properly
  DE_table_all[is.na(DE_table_all)] = 0
  
  # Combine ligand activities with DE information
  ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
  ligand_activities_de[is.na(ligand_activities_de)] = 0
  
  # make LFC heatmap
  lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
  rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()
  
  order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
  vis_ligand_lfc = lfc_matrix[order_ligands,]
  
  colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()
  
  p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
  p_ligand_lfc
  
  # change colors a bit to make them more stand out
  p_ligand_lfc = p_ligand_lfc + scale_fill_gradientn(colors = c("midnightblue","blue", "grey95", "grey99","firebrick1","red"),values = c(0,0.1,0.2,0.25, 0.40, 0.7,1), limits = c(vis_ligand_lfc %>% min() - 0.1, vis_ligand_lfc %>% max() + 0.1))
  p_ligand_lfc
  
  
  # Follow-up analysis 2: Visualize expression of top-predicted ligands and 
  # their target genes in a combined heatmap
  library(RColorBrewer)
  library(cowplot)
  library(ggpubr)
  
  ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% 
    magrittr::set_rownames(ligand_activities$test_ligand)
  
  vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% 
    magrittr::set_colnames("Pearson")
  
  p_ligand_pearson = vis_ligand_pearson %>% 
    make_heatmap_ggplot("Prioritized CAF-ligands","Ligand activity", 
                        color = "darkorange",
                        legend_position = "top", 
                        x_axis_position = "top", 
                        legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)")
  p_ligand_pearson
  
  # Prepare expression of ligands in fibroblast per tumor
  expression_df_CAF = expression[CAF_ids,order_ligands] %>% 
    data.frame() %>% rownames_to_column("cell") %>% as_tibble() %>% 
    inner_join(sample_info %>% select(cell,tumor), by =  "cell")
  
  aggregated_expression_CAF = expression_df_CAF %>% group_by(tumor) %>% select(-cell) %>% summarise_all(mean)
  aggregated_expression_df_CAF = aggregated_expression_CAF %>% select(-tumor) %>% t() %>% 
    magrittr::set_colnames(aggregated_expression_CAF$tumor) %>% data.frame() %>% rownames_to_column("ligand") %>% as_tibble() 
  
  aggregated_expression_matrix_CAF = aggregated_expression_df_CAF %>% select(-ligand) %>% 
    as.matrix() %>% magrittr::set_rownames(aggregated_expression_df_CAF$ligand)
  
  order_tumors = c("HN6","HN20","HN26","HN28","HN22","HN25","HN5","HN18","HN17","HN16") 
  # this order was determined based on the paper from Puram et al. Tumors are ordered according to p-EMT score.
  vis_ligand_tumor_expression = aggregated_expression_matrix_CAF[order_ligands,order_tumors]
  
  library(RColorBrewer)
  color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
  p_ligand_tumor_expression = vis_ligand_tumor_expression %>% 
    make_heatmap_ggplot("Prioritized CAF-ligands","Tumor", color = color[100],
                        legend_position = "top", x_axis_position = "top", 
                        legend_title = "Expression\n(averaged over\nsingle cells)") + 
    theme(axis.text.y = element_text(face = "italic"))
  p_ligand_tumor_expression
  
  # Prepare expression of target genes in malignant cells per tumor
  expression_df_target = expression[malignant_ids,geneset_oi] %>% data.frame() %>% 
    rownames_to_column("cell") %>% as_tibble() %>% inner_join(sample_info %>% select(cell,tumor), by =  "cell") 
  
  aggregated_expression_target = expression_df_target %>% group_by(tumor) %>% 
    select(-cell) %>% summarise_all(mean)
  
  aggregated_expression_df_target = aggregated_expression_target %>% select(-tumor) %>% t() %>% 
    magrittr::set_colnames(aggregated_expression_target$tumor) %>% 
    data.frame() %>% rownames_to_column("target") %>% as_tibble() 
  
  aggregated_expression_matrix_target = aggregated_expression_df_target %>% 
    select(-target) %>% as.matrix() %>% magrittr::set_rownames(aggregated_expression_df_target$target)
  
  vis_target_tumor_expression_scaled = aggregated_expression_matrix_target %>% t() %>% 
    scale_quantile() %>% .[order_tumors,order_targets]
  
  p_target_tumor_scaled_expression = vis_target_tumor_expression_scaled  %>% 
    make_threecolor_heatmap_ggplot("Tumor","Target", low_color = color[1], 
                                   mid_color = color[50], mid = 0.5, 
                                   high_color = color[100], 
                                   legend_position = "top", x_axis_position = "top" , 
                                   legend_title = "Scaled expression\n(averaged over\nsingle cells)") +
    theme(axis.text.x = element_text(face = "italic"))
  p_target_tumor_scaled_expression
  
  
  ## Inferring ligand-to-target signaling paths
  ligand_tf_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_tf_matrix.rds"))
  sig_network = readRDS(url("https://zenodo.org/record/3260758/files/signaling_network.rds"))
  gr_network = readRDS(url("https://zenodo.org/record/3260758/files/gr_network.rds"))
  
  ligands_all = "TGFB3" # this can be a list of multiple ligands if required
  targets_all = c("TGFBI","LAMC2","TNC")
  
  active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, 
                                                       ligands_all = ligands_all, 
                                                       targets_all = targets_all, 
                                                       weighted_networks = weighted_networks)
  
  # For better visualization of edge weigths: normalize edge weights to 
  # make them comparable between signaling and gene regulatory interactions
  active_signaling_network_min_max = active_signaling_network
  active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% 
    mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
  active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% 
    mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
  
  graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, 
                                                    ligands_all = ligands_all, 
                                                    targets_all = targets_all, 
                                                    sig_color = "indianred", 
                                                    gr_color = "steelblue")
  
  # To render the graph: uncomment following line of code
  DiagrammeR::render_graph(graph_min_max, layout = "tree")
  
  
  ## Read in the expression data of interacting cells
  seuratObj = readRDS(paste0(dataPath_nichenet,  "seuratObj.rds"))
  seuratObj@meta.data %>% head()
  
  sels =c(#which(refs$celltype == 'CM')[1:1000], 
    which(refs$celltype == 'FB')[1:1000], 
    which(refs$celltype == 'prolife.Mphage'))
  subref = subset(refs, cells = colnames(refs)[sels])
  Idents(subref) = subref$celltype
  subref@meta.data$celltype %>% table()
  
  # note that the number of cells of some cell types is very low and should preferably be higher for a real application
  seuratObj@meta.data$celltype %>% table() 
  
  DimPlot(seuratObj, reduction = "tsne")
  
  seuratObj@meta.data$aggregate %>% table()
  
  DimPlot(seuratObj, reduction = "tsne", group.by = "aggregate")
  
  ## Perform the NicheNet analysis
  # indicated cell types should be cell class identities
  # check via: 
  # seuratObj %>% Idents() %>% table()
  nichenet_output = nichenet_seuratobj_aggregate(
    seurat_obj = seuratObj, 
    receiver = "CD8 T", 
    condition_colname = "aggregate", 
    condition_oi = "LCMV", 
    condition_reference = "SS", 
    sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"),
    
    ligand_target_matrix = ligand_target_matrix, 
    lr_network = lr_network, 
    weighted_networks = weighted_networks, 
    organism = "mouse")
  
  ## Interpret the NicheNet analysis output
  nichenet_output$ligand_activities
  
  nichenet_output$top_ligands
  nichenet_output$ligand_expression_dotplot
  
  nichenet_output$ligand_differential_expression_heatmap
  
  nichenet_output$ligand_target_heatmap
  
  nichenet_output$ligand_target_heatmap + 
    scale_fill_gradient2(low = "whitesmoke",  high = "royalblue", breaks = c(0,0.0045,0.009)) + 
    xlab("anti-LCMV response genes in CD8 T cells") + 
    ylab("Prioritized immmune cell ligands")
  
  nichenet_output$ligand_activity_target_heatmap
  
  features = rownames(st)[grep('BMP7', rownames(st))]
  SpatialFeaturePlot(st,  features = features)
  
  FeaturePlot(refs, features = rownames(refs)[grep('BMP7', rownames(refs))])

}