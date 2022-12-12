
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
require()
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
outDir = paste0(resDir, '/Ligand_Receptor_analysis/NicheNet_v3.2_clusterDE_timeSpecific_receiverCell.CM_IS')
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
#lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% 
#  distinct(ligand, receptor, bonafide)

# If wanted, users can remove ligand-receptor interactions that were predicted based on 
# protein-protein interactions and 
# only keep ligand-receptor interactions that are described in curated databases.
# lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
head(lr_network)

#ligands = lr_network %>% pull(ligand) %>% unique()
#receptors = lr_network %>% pull(receptor) %>% unique()

# ## weighted integrated networks 
weighted_networks = readRDS(paste0(dataPath_nichenet,  "weighted_networks.rds"))
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

weighted_networks_lr = weighted_networks$lr_sig %>% 
  #dplyr::rename(ligand = from, receptor = to) %>%
  inner_join(lr_network %>% 
               distinct(from, to), by = c("from","to"))

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
    
    hist(log10(ave.counts), breaks=100, main="", col="grey80",
         xlab=expression(Log[10]~"average count"))
    
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
  
  # to save the result for each time point
  expressed_genes = c()
  res = tibble(sample=character(), sender=character(), receiver=character(),  
               test_ligand=character(), auroc=double(), aupr=double(), 
               pearson = double(),  rank=integer())
  
  # step 1:  Define a “sender/niche” cell population and a “receiver/target” cell population 
  # present in your expression data and determine which genes are expressed in both populations
  if(is.na(receiver_cells)){ # loop over all cell type candidates
    receiver_cells = celltypes
  }
  
  for(receiver in receiver_cells) # loop over receiver cells 
  {
    # specify receiver
    # receiver = receiver_cells[1]
    cat('-- receiver : ', receiver, '-- \n')
    #expressed_genes_receiver = get_expressed_genes(receiver, subref, pct = 0.10)
    #background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
    
    # step 2: Define a gene set of interest: these are the genes in the “receiver/target” cell population 
    # that are potentially affected by ligands expressed by interacting cells 
    # (e.g. genes differentially expressed upon cell-cell interaction)
    #DE_table_receiver = FindMarkers(object = subref, 
    #                                ident.1 = receiver, 
    #                                ident.2 = control_cells, 
    #                                min.pct = 0.10) %>% rownames_to_column("gene")
    
    #geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
    #geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
    
    ## sender
    #sender_celltypes = celltypes
    
    ## aggregate function from nichenet
    nichenet_output = nichenet_seuratobj_cluster_de(seurat_obj = subref, 
                                                    receiver_affected = receiver, 
                                                    receiver_reference = control_cells, 
                                                    sender = celltypes, 
                                                    ligand_target_matrix, 
                                                    lr_network,
                                                    weighted_networks, 
                                                    expression_pct = 0.10, 
                                                    lfc_cutoff = 0.25, 
                                                    geneset = "DE", 
                                                    filter_top_ligands = TRUE, 
                                                    top_n_ligands = 100, 
                                                    top_n_targets = 100, 
                                                    cutoff_visualization = 0.33, 
                                                    organism = "human",
                                                    verbose = TRUE, assay_oi = NULL)
    
    nichenet_output$ligand_activities
    # res = rbind(res, nichenet_output$ligand_activities) 
    
    nichenet_output$top_ligands
    p1 = nichenet_output$ligand_expression_dotplot
    plot(p1) 
    ggsave(filename = paste0(outDir, '/ligand_expression_dotplot_', timepoint, '.pdf'), 
           width = 12, height = 20)
    
    nichenet_output$ligand_differential_expression_heatmap
    
    p2 = nichenet_output$ligand_target_heatmap + 
      scale_fill_gradient2(low = "whitesmoke",  high = "royalblue")
    plot(p2) 
    ggsave(filename = paste0(outDir, '/ligand_target_heatmap_', timepoint, '.pdf'), 
           width = 25, height = 20)
    
    nichenet_output$top_targets
    
    DotPlot(subref %>% subset(idents = c(receiver, control_cells)), 
            features = nichenet_output$top_targets %>% rev())
    
    p3 = nichenet_output$ligand_activity_target_heatmap
    plot(p3)
    
    ggsave(filename = paste0(outDir, '/ligand_activity_target_heatmap_', timepoint, '.pdf'), 
           width = 27, height = 20)
    
    
  }
  
  ##########################################
  # combine receiver cells to make one plot for one time point
  ##########################################
  Combine.multiple.receiver_cells = FALSE 
  if(Combine.multiple.receiver_cells){
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
    res %>% ggplot(aes(x=pearson)) + 
      geom_histogram(color="black", fill="darkorange") + 
      geom_vline(aes(xintercept=0.10), color="red", linetype="dashed", size=1) +  
      labs(x="ligand activity (Pearson Correlation)", y = "# ligands") +  
      theme_classic() + ggtitle("Ligands activity scores")
    
    # Plot ligand activities scores per interaction
    res %>% ggplot(aes(x=pearson)) + 
      geom_histogram(color="black", fill="darkorange") + 
      geom_vline(aes(xintercept=0.1), color="red", linetype="dashed", size=1) +  
      labs(x="ligand activity (Pearson Correlation)", y = "# ligands") +  
      theme_classic() + 
      ggtitle("Ligands activity scores") + 
      facet_wrap(.~interaction, )
    
    # Visualize ligand expression along different days
    ligand.expression.df = AverageExpression(subref, features =unique(res$test_ligand))[['RNA']]
    pseudo_signal = ceiling(log2(range(ligand.expression.df[ligand.expression.df>0]))[1])-1
    ligand.expression.df = log2(2^pseudo_signal +ligand.expression.df)[sort(rownames(ligand.expression.df)), 
                                                                       sort(colnames(ligand.expression.df))]
    write.csv2(inner_join(res, ligand.expression.df %>% as.data.frame %>% rownames_to_column('test_ligand')), 
               file = paste0(outDir, '/predicted_ligand_activity_avExpr_for_receivers_senders_', timepoint, '.csv'),
               row.names = FALSE)
    
    #p2 = Heatmap(ligand.expression.df, cluster_columns = F, column_title = "Ligand Expression", name="Expression (log2)", row_order = row_order(p1))
    
    # Visualization per interaction
    p1 = res %>% transmute(test_ligand=test_ligand, x=interaction, pearson=pearson) %>% 
      #spread(x, pearson)  %>% 
      replace(is.na(.), 0) %>% 
      ungroup() %>% 
      column_to_rownames('test_ligand')
    
    ha.col = list(Interaction=setNames(c(brewer.pal(9, 'YlOrRd')), colnames(p1)))
    ha = HeatmapAnnotation(Interaction=colnames(p1), annotation_name_side = 'left', col=ha.col)
    
    p1 = p1 %>% Heatmap(
      column_title = "Average Ligand pearson correlation\n(NAs replaced with zeros)", 
      cluster_columns = F, show_row_names = F, show_column_names = F, 
      top_annotation = ha, col=brewer.pal(4, 'Oranges'), 
      name="Ligand Score(pearson)")
    ha = colnames(ligand.expression.df)
    ha = HeatmapAnnotation(Cluster=ha, col=ha.col)
    p2 = Heatmap(ligand.expression.df, 
                 cluster_columns = F, show_column_names = F,
                 column_title = "Ligand Expression", name="log10(Expression)", 
                 row_order = row_order(p1), 
                 top_annotation = ha)
    
    p1+p2
    
  }
  
}

