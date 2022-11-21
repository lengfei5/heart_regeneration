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
run_LIANA = function(refs,
                     celltypes = c('Mono_Macrophages', 'Proliferating_CM', 'Neutrophil', 'Injury_specific_EC'),
                     outDir = '../results/Ligand_Receptor_analysis',
                     ntop = 50,
                     RUN.CPDB.alone = FALSE
) # original code from https://saezlab.github.io/liana/articles/liana_tutorial.html
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
  
  system(paste0('mkdir -p ', outDir))
  # liana_path <- system.file(package = "liana")
  # testdata <- readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
  # 
  # testdata %>% glimpse()
  
  Idents(refs) = as.factor(refs$celltypes)
  subref = subset(refs, cells = colnames(refs)[!is.na(match(refs$celltypes, celltypes))])
  subref$celltypes = droplevels(as.factor(subref$celltypes))
  table(subref$celltypes)
  #subref = subset(x = subref, downsample = 1000)
  cat('celltype to consider -- ', names(table(subref$celltypes)), '\n')
  
  # sels =c(which(refs$celltypes == 'CM')[1:1000], 
  #         which(refs$celltypes == 'FB')[1:1000], 
  #         which(refs$celltypes == 'Macrophages')[1:1000], 
  #         which(refs$celltypes == 'Neutrophil'))
  
  #subref = subset(refs, cells = colnames(refs)[sels])
  
  #rownames(subref) = toupper(rownames(subref))
  
  Idents(subref) = subref$celltypes
  
  # Run liana
  # liana_test <- liana_wrap(testdata, method = 'cellphonedb', resource = 'CellPhoneDB')
  source('functions_scRNAseq.R')
  sce <- as.SingleCellExperiment(subref)
  colLabels(sce) = as.factor(sce$celltypes)
  rownames(sce) = toupper(get_geneName(rownames(sce)))
  
  ave.counts <- calculateAverage(sce, assay.type = "counts")
  
  hist(log10(ave.counts), breaks=100, main="", col="grey80",
       xlab=expression(Log[10]~"average count"))
  
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
  
  liana_test <- liana_wrap(sce,  
                           method = c("natmi", "connectome", "logfc", "sca", "cytotalk",  
                                      'cellphonedb'),
                           resource = c("Consensus", 'CellPhoneDB', "OmniPath", "LRdb",
                                        "CellChatDB",  "CellTalkDB"), 
                           assay.type = "logcounts", 
                           idents_col = 'celltypes')
  
  # Liana returns a list of results, each element of which corresponds to a method
  liana_test %>% glimpse
  
  # We can aggregate these results into a tibble with consensus ranks
  liana_test <- liana_test %>%
    liana_aggregate(resource = 'Consensus')
  
  saveRDS(liana_test, file = paste0(outDir, '/res_lianaTest_Consensus.rds'))
  
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
           width = 16, height = min(c(10*ntop/20, 50)), limitsize = FALSE)
    
  }
  
  pdfname = paste0(outDir, '/liana_celltype_communication_freqHeatmap.pdf')
  pdf(pdfname, width=12, height = 8)
  
  liana_trunc <- liana_test %>%
    # only keep interactions concordant between methods
    filter(aggregate_rank <= 0.01) # this can be FDR-corr if n is too high
  
  heat_freq(liana_trunc)
  
  dev.off()
  
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
