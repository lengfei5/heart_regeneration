##########################################################################
##########################################################################
# Project: heart regeneration
# Script purpose: make plots
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Jan 23 16:03:26 2023
##########################################################################
##########################################################################
make_plot_umap_scRNA = FALSE
if(make_plot_umap_scRNA){
  annots = readRDS('/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds')
  annots = annots[, c(1:7, 16)]
  
  write.table(annots, file = paste0("/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/", 
                                    "geneAnnotation_curated.geneSymbol.toUse.txt"), 
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  ##########################################
  # CM umap  
  ##########################################
  
  aa1 <- subset(aa,  subtypes %in%  c( 
    
    "B_cells_(FOXO1)", "B_cells_(SIGLEC11)",             "B_cells_Prol"     ,               
    "CM_Atria" ,
    
    "CM_Atria_Tagln", 
    "CM_IS"          ,             
    "CM_OFT"         ,        "CM_PM_(HCN4)" ,
    
    "CM_Prol_1"       ,             "CM_Prol_2"      ,              "CM_Prol_3"           ,       
    "CM_Prol_IS" ,
    
    "CM_ven_(Cav3_1)"        ,       "CM_ven_(Robo2)"       ,                                      "EC" ,
    
    "EC_(CEMIP)"      ,              "EC_(LHX6)"    ,                "EC_(NOS3)" ,                   "EC_(WNT4)" ,
    
    "EC_IS_(IARS1)"    ,              "EC_IS_(LOX)"   ,                "EC_IS_Prol"  ,                    "EC_Prol" ,
    
    "FB_(PKD1)"  ,                  "FB_(TNXB)"  ,                  "FB_(VWA2)" ,               "FB_IS_(TFPI2)" ,
    
    "FB_IS_(TNC)"  ,                    "FB_Prol",               "Megakeryocytes"  ,           "Mo/Macs_(FAXDC2)" ,
    
    "Mo/Macs_(SNX22)" ,                "Mo/Macs_Prol" ,            "Mo/Macs_resident"  ,                 "Neu_(DYSF)" ,
    
    "Neu_(IL1R1)" ,                    "Neuronal", "Proliferating_Megakeryocytes",            "Proliferating_RBC",
    
    "RBC"  ,                    "T_cells" ))
  
  
}


##########################################
# test linear discriminant analysis (LDA) 
##########################################
Test_LDA = FALSE
if(Test_LDA){
  library(MASS)
  library(ggplot2)
  
  #attach iris dataset to make it easy to work with
  attach(iris)
  
  #view structure of dataset
  str(iris)
  
  #scale each predictor variable (i.e. first 4 columns)
  iris[1:4] <- scale(iris[1:4])
  
  #find mean of each predictor variable
  apply(iris[1:4], 2, mean)
  
  #find standard deviation of each predictor variable
  apply(iris[1:4], 2, sd) 
  
  #make this example reproducible
  set.seed(1)
  
  #Use 70% of dataset as training set and remaining 30% as testing set
  sample <- sample(c(TRUE, FALSE), nrow(iris), replace=TRUE, prob=c(0.7,0.3))
  train <- iris[sample, ]
  test <- iris[!sample, ] 
  
  #fit LDA model
  model <- lda(Species~., data=train)
  
  #view model output
  model
  #define data to plot
  
  #use LDA model to make predictions on test data
  predicted <- predict(model, train)
  
  head(predicted$class)
  #view posterior probabilities for first six observations in test set
  head(predicted$posterior)
  
  #view linear discriminants for first six observations in test set
  head(predicted$x)
  
  lda_plot <- cbind(train, predicted$x)
  
  #create plot
  ggplot(lda_plot, aes(LD1, LD2)) +
    geom_point(aes(color = Species))
  
}


##########################################
# plot heatmap of features along trajectories
##########################################
get_smooth_curve_spline = function(x, t, newt, downsample = TRUE)
{
  # x = as.numeric(cds[1, jj]); t = Pseudotime; newt = pseudot_comomon;
  if(downsample){
    nb_t = min(5000, length(t))
    nn = sample(1:length(t), size = nb_t, replace = FALSE)
    t = t[nn]
    x = x[nn]
  }
  
  fit_sel = smooth.spline(t, x, df = 3)
  
  #plot(Pseudotime, cds_sel, cex = 0.5)
  #lines(fit_sel, col = 'red', lwd =2.0)
  newx = predict(fit_sel, newt)
  return(newx$y)
  #VGAM::vglm(~sm.ns(Pseudotime, df=3), family = 'gaussian', data = cds_sel)
  
}

## original functions from Monocle
## https://github.com/cole-trapnell-lab/monocle-release/blob/master/R/plotting.R
plot_features_branched_heatmap <- function(seuratObj = aa,
                                           feature = 'gene',
                                        feature_subset = candidates,
                                        nbCell_condition = 50,
                                        Get.Smooth.Curve = TRUE, 
                                        scale_max=3, 
                                        scale_min=-3, 
                                        hmcols = NULL, 
                                        hclust_method = "ward.D2", 
                                        num_clusters = 6,
                                        cluster_rows = TRUE,
                                        add_annotation_row = NULL,
                                        show_rownames = TRUE) 
{
  # seuratObj = chromvar_CM;nbCell_condition = 50;scale_max=3; scale_min=-3; Get.Smooth.Curve = TRUE;
  # feature_subset = feature_sels;hclust_method = "ward.D2";num_clusters = 6;hmcols = NULL;
  library(VGAM) # an example code from https://online.stat.psu.edu/stat504/lesson/8/8.2/8.2.2
  library(MASS)
  library(tidyr)
  require(colorRamps)
  require(pheatmap)
  
  ## subsampling cells to make the heatmap symmetric and also avoid too many cells in some cell types
  subsampling.cells_perCondition = FALSE
  if(subsampling.cells_perCondition){
    
    subs = seuratObj
    Idents(subs) = as.factor(subs$subtypes)
    subs = subset(x = subs, downsample = 1000)
    
    Idents(subs) = as.factor(subs$branch)
    subs = subset(x = subs, downsample = min(table(subs$branch)))
    
    cell_ds = colnames(subs)
    
    saveRDS(cell_ds, file = paste0(outDir, 'saved_cells_downsampled.rds'))
    rm(subs)
    
  }
  
  ## select the corresponding matrix of selected features
  if(feature == 'gene') cds <- seuratObj@assays$RNA@scale.data
  if(feature == 'motif_activity') cds = seuratObj@assays$chromvar@data
  rownames(cds) = rownames(seuratObj)
  cds = cds[which(!is.na(match(rownames(cds), feature_subset))), ]
  
  ## smooth the gene expression along the trajectory
  if(Get.Smooth.Curve){
    cat(' -- smoothing the single cell data for subsampled cells -- \n')
    
    jj = which(seuratObj$branch == 'injury')
    Pseudotime = as.numeric(seuratObj$pst[jj])
    names = colnames(seuratObj)[jj]
    #kk_BrachA = jj
    #kk_BrachA = grep('_RA$|_RA.rep1', subs$condition)
    #kk_BrachA = kk_BrachA[order(subs$pseudot[kk_BrachA])]
    #pseudot_BrachA = subs$pseudot[kk_BrachA]
    newt = Pseudotime[order(Pseudotime)]
    names = names[order(Pseudotime)]
    BranchA_exprs <- t(apply(cds[,jj], 1, get_smooth_curve_spline, t = Pseudotime, newt = newt))
    colnames(BranchA_exprs) = names
    
    jj = which(seuratObj$branch == 'uninjury')
    Pseudotime = as.numeric(seuratObj$pst[jj])
    names = colnames(seuratObj)[jj]
    newt = Pseudotime[order(Pseudotime)]
    names = names[order(Pseudotime)]
    BranchB_exprs <- t(apply(cds[ ,jj], 1, get_smooth_curve_spline, t = Pseudotime, 
                             newt = newt))
    colnames(BranchB_exprs) = names
    
    save(BranchA_exprs, BranchB_exprs, file = paste0(outDir, 'smoothed_geneExpr_forHeatmap.Rdata'))
    
  }
  
  # reload the downsampled cells
  cell_ds = readRDS(file = paste0(outDir, 'saved_cells_downsampled.rds'))
  cat(length(cell_ds), ' pre-selected \n')
  cell_ds = intersect(cell_ds, colnames(seuratObj))
  cat(length(cell_ds), ' selected to use \n')
  
  ## use only cells saved in the previously downsampled
  load(file = paste0(outDir, 'smoothed_geneExpr_forHeatmap.Rdata'))
  BranchA_exprs = BranchA_exprs[, !is.na(match(colnames(BranchA_exprs), cell_ds))]
  BranchB_exprs = BranchB_exprs[, !is.na(match(colnames(BranchB_exprs), cell_ds))]
  col_gap_ind <- c(ncol(BranchB_exprs))
  
  heatmap_matrix <- cbind(BranchB_exprs[, ncol(BranchB_exprs):1], 
                        BranchA_exprs)
  
  if(feature == 'motif_activity'){
    scale_max=2; scale_min=-2;
  }
  heatmap_matrix=heatmap_matrix[!apply(heatmap_matrix, 1, sd)==0,]
  #if(feature == 'gene') 
  heatmap_matrix=Matrix::t(scale(Matrix::t(heatmap_matrix),center=TRUE))
  heatmap_matrix=heatmap_matrix[is.na(row.names(heatmap_matrix)) == FALSE, ]
  heatmap_matrix[is.nan(heatmap_matrix)] = 0
  heatmap_matrix[heatmap_matrix>scale_max] = scale_max
  heatmap_matrix[heatmap_matrix<scale_min] = scale_min
  
  saveRDS(heatmap_matrix, file = paste0(outDir, "/heatmap_matrix_forPlot.rds"))
  #heatmap_matrix_ori <- heatmap_matrix
  #heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[, 1]) & 
  # is.finite(heatmap_matrix[, col_gap_ind]), ] #remove the NA fitting failure genes for each branch 
  
  heatmap_matrix = readRDS(file = paste0(outDir, "/heatmap_matrix_forPlot.rds"))
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  
  exp_rng <- range(heatmap_matrix) #bks is based on the expression range
  
  bks <- seq(exp_rng[1] - 0.1, exp_rng[2] + 0.1, by=0.1)
  if(is.null(hmcols)) {
    hmcols <- blue2green2red(length(bks) - 1)
  }
  
  num_clusters = 4
  
  # prin  t(hmcols)
  ph <- pheatmap(heatmap_matrix, 
                 useRaster = T,
                 cluster_cols=FALSE, 
                 cluster_rows=TRUE, 
                 show_rownames=F, 
                 show_colnames=F, 
                 #scale="row",
                 clustering_distance_rows=row_dist, #row_dist
                 clustering_method = hclust_method,
                 cutree_rows=num_clusters,
                 silent=TRUE,
                 #filename=NA,
                 breaks=bks,
                 color=hmcols
                 #color=hmcols#
  )
  
  cells = colnames(heatmap_matrix)
  annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)),
                               pseudot = seuratObj$pst[match(cells, colnames(seuratObj))],
                               subtype = seuratObj$subtypes[match(cells, colnames(seuratObj))],
                               branch = seuratObj$branch[match(cells, colnames(seuratObj))])
  
  #write.csv2(annotation_row, file = paste0(outDir, 'gene_clusters.csv'), row.names = TRUE, quote = FALSE)
  
  #colnames(annotation_col) <- "Cell Type"  
  #if(!is.null(add_annotation_col)) {
  #  annotation_col <- cbind(annotation_col, add_annotation_col[fData(cds[row.names(annotation_col), ])$gene_short_name, 1])  
  #}
  #branch_colors = cols_sel[grep("day3_RA.rep2", names(cols_sel), invert = TRUE)]
  #annotation_colors=list("condition"=branch_colors)
  #names(annotation_colors$`condition`) = c('Pre-branch', branch_labels)
  
  #names(branch_colors) <- c("Pre-branch", branch_labels[1], branch_labels[2])
  #annotation_colors=list("Cell Type"=branch_colors)
  #names(annotation_colors$`Cell Type`) = c('Pre-branch', branch_labels)
  
  if(feature == 'motif_activity'){
    motif_tf = readRDS(file = paste0('../results/sc_multiome_R13591_atac_reseq_20221115/Rdata/', 
                                     'motif_to_tfs_pfm_JASPAR2020_CORE_vertebrate_v1.rds'))
    feature_label = motif_tf$name[match(rownames(heatmap_matrix), motif_tf$motif)]
    row_ann_labels = motif_tf$name[match(rownames(annotation_row), motif_tf$motif)]
  }else{
    feature_label <- row.names(heatmap_matrix)
    row_ann_labels <- row.names(annotation_row)
  }
  
  row.names(heatmap_matrix) <- feature_label
  row.names(annotation_row) <- row_ann_labels
  
  
  ##########################################
  # heatmap plotting
  ##########################################
  Plot_use_heatmap = TRUE
  if(Plot_use_heatmap){
    if(feature == 'gene') myColor = hmcols
    if(feature == 'motif_activity') myColor = inlmisc::GetColors(length(bks) +1, scheme = "BuRd")
    
    # all DE genes
    pheatmap(heatmap_matrix[, ], #ph$tree_row$order
             useRaster = T,
             cluster_cols=FALSE, 
             cluster_rows=TRUE, 
             show_rownames=FALSE,
             show_colnames=FALSE, 
             scale='none',
             clustering_distance_rows=row_dist, #row_dist
             clustering_method = hclust_method, #ward.D2
             cutree_rows=num_clusters,
             # cutree_cols = 2,
             annotation_row=annotation_row,
             annotation_col=annotation_col,
             #annotation_colors=annotation_colors,
             gaps_col = col_gap_ind,
             treeheight_row = 30, 
             breaks=bks,
             fontsize = 6,
             color=myColor, 
             border_color = NA,
             silent=TRUE, 
             filename=paste0(outDir, "/expression_pseudotime_pheatmap_DEtrajectory_", feature, ".pdf"),
             width = 5, height = 8
    )
    
    pheatmap(heatmap_matrix[, ], #ph$tree_row$order
             useRaster = T,
             cluster_cols=FALSE, 
             cluster_rows=TRUE, 
             show_rownames=TRUE,
             show_colnames=FALSE, 
             scale='none',
             clustering_distance_rows=row_dist, #row_dist
             clustering_method = hclust_method, #ward.D2
             cutree_rows=num_clusters,
             # cutree_cols = 2,
             annotation_row=annotation_row,
             annotation_col=annotation_col,
             #annotation_colors=annotation_colors,
             gaps_col = col_gap_ind,
             treeheight_row = 50, 
             breaks=bks,
             fontsize = 6,
             color=myColor, 
             border_color = NA,
             silent=TRUE, 
             filename=paste0(outDir, "/expression_pseudotime_pheatmap_DEtrajectory_", feature, "withNames.pdf"),
             width = 8, height = 12
    )
    
    ##########################################
    # DE TFs and signaling pathways 
    ##########################################
    tfs = readRDS(file = paste0('../../RA_competence/data/annotations/curated_human_TFs_Lambert.rds'))
    tfs = unique(tfs$`HGNC symbol`)
    
    ggs = get_geneName(rownames(heatmap_matrix))
    sels = which(!is.na(match(ggs, tfs)))
    
    #sels = which(!is.na(match(rownames(heatmap_matrix), targets)))
    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix[sels,])))/2)
    row_dist[is.na(row_dist)] <- 1
    annotation_rowSel = data.frame(Cluster = annotation_row[sels, ])
    rownames(annotation_rowSel) = rownames(annotation_row)[sels]
    
    pheatmap(heatmap_matrix[sels, ], #ph$tree_row$order
             useRaster = T,
             cluster_cols=FALSE, 
             cluster_rows=TRUE, 
             show_rownames=FALSE,
             show_colnames=FALSE, 
             scale='none',
             clustering_distance_rows=row_dist, #row_dist
             clustering_method = hclust_method, #ward.D2
             #cutree_rows=num_clusters,
             # cutree_cols = 2,
             annotation_row=annotation_rowSel,
             annotation_col=annotation_col,
             #annotation_colors=annotation_colors,
             gaps_col = col_gap_ind,
             treeheight_row = 30, 
             breaks=bks,
             fontsize = 6,
             color=hmcols, 
             border_color = NA,
             silent=TRUE, 
             filename=paste0(outDir, "/expression_pseudotime_pheatmap_allDEgenes_TF.pdf"),
             width = 5, height = 12
    )
    
    pheatmap(heatmap_matrix[sels, ], #ph$tree_row$order
             useRaster = T,
             cluster_cols=FALSE, 
             cluster_rows=TRUE, 
             show_rownames=TRUE,
             show_colnames=FALSE, 
             scale='none',
             clustering_distance_rows=row_dist, #row_dist
             clustering_method = hclust_method, #ward.D2
             cutree_rows=num_clusters,
             # cutree_cols = 2,
             annotation_row=annotation_rowSel,
             annotation_col=annotation_col,
             annotation_colors=annotation_colors,
             gaps_col = col_gap_ind,
             treeheight_row = 30, 
             breaks=bks,
             fontsize = 4,
             color=hmcols, 
             border_color = NA,
             silent=TRUE, 
             filename=paste0(outDir, "/expression_pseudotime_pheatmap_allDEgenes_TF.SP",
                             "_withgeneNames.pdf"),
             width = 8, height = 20
    )
    
  }
  
  
}

