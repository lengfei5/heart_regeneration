##########################################################################
##########################################################################
# Project:
# Script purpose: LIANA hacking hepl functions
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Mar  7 15:31:03 2022
##########################################################################
##########################################################################
row_scale <- function(mat){
  col_means = rowMeans(mat,
                       na.rm = TRUE) # Get the column means
  col_sd = MatrixGenerics::rowSds(mat,
                                  center = col_means,
                                  na.rm = TRUE) # Get the column sd
  
  # return scaled mat
  return(as.matrix((mat - col_means) / col_sd))
}

fast_mean <- function(mat){
  if(class(mat)=="dgCMatrix"){
    sum(mat@x)/(as.numeric(nrow(mat)) * as.numeric(ncol(mat)))
  } else{
    Matrix::mean(mat)
  }
}

compute_pem_scores <- function(sce, assay.type = "logcounts") {
  
  # extract normalized data
  norm_data <- t(as.matrix(sce@assays@data[[assay.type]]))
  
  # expm1(x) computes exp(x) - 1 accurately also for |x| << 1
  exp_data <- exp(norm_data)
  
  # calculate mean per cell using expm1 data
  means_per_cell <- exp_data %>%
    aggregate(., list(colLabels(sce)), FUN = mean) %>%
    tibble::as_tibble() %>%
    dplyr::rename(celltype = Group.1) %>%
    pivot_longer(-celltype, names_to = "gene") %>%
    tidyr::pivot_wider(names_from = celltype,
                       id_cols = gene,
                       values_from = value) %>%
    column_to_rownames("gene")
  
  # calculate PEMs from the mean matrix
  pem_out <- pems_from_means(means_per_cell)
  
  return(pem_out)
  
}

#' PEM from mean matrix
#'
#' Computes the Preferential Expression Measure (PEM) scores from a matrix containing
#' genes as rows and cell types as samples. The values in the matrix are supposed
#' to contain the average expression values of the count matrix per cell, obtained
#' from the exponential of the log-transformed count matrix.
#'
#' @param means_per_cell A matrix containing the mean expression per cell type
#'
#' @return A matrix with the computed PEM scores
#'
pems_from_means <- function(means_per_cell) {
  
  # calculate sum of means for rows and columns
  sums_per_cell <- colSums(means_per_cell)
  sums_per_gene <- rowSums(means_per_cell)
  total_sums <- sum(sums_per_cell)
  
  # will lead to bugs later, could check downstream
  if (total_sums == Inf) {
    message("Inf introduced into PEM score, is the data log1p tranformed?")
    return(1)
  }
  
  # what proportion of this cell type's rowmean sum accounts for the whole?
  cell_type_prop <- sums_per_cell / total_sums
  
  # gene_proportions
  gene_prop <- sapply(cell_type_prop, function(x) x * sums_per_gene)
  
  # calculate PEM scores
  pem_out <- log10(means_per_cell / gene_prop)
  
  # set negative and NAs PEMs to 0
  pem_out[pem_out < 0 | is.na(pem_out)] <- 0
  
  return(pem_out)
  
}

get_log2FC <- function(sce,
                       assay.type){
  
  # normalize counts across libraries
  if(assay.type == 'count'){
    sce <- scater::logNormCounts(sce,
                                log=FALSE,
                                 assay.type=assay.type)
  }
  
  # iterate over each possible cluster leaving one out
  levels(colLabels(sce)) %>%
    map(function(subject){
      # Subject (i.e. target) Cluster avg
      subject_avg <-
        scater::calculateAverage(subset(sce,
                                        select = colLabels(sce)==subject),
                                 assay.type = assay.type
        ) %>%
        as_tibble(rownames = "gene") %>%
        dplyr::rename(subject_avg = value)
      
      # All other cells average
      loso_avg <-
        scater::calculateAverage(subset(sce,
                                        select = !(colLabels(sce) %in% subject)),
                                 assay.type = assay.type
        ) %>%
        as_tibble(rownames = "gene") %>%
        dplyr::rename(loso_avg = value)
      
      # Join avg and calculate FC
      left_join(subject_avg, loso_avg, by="gene") %>%
        mutate(avg_log2FC =
                 log2((subject_avg + 1)) - log2((loso_avg + 1))) %>%
        select(gene, avg_log2FC)
      
    }) %>% setNames(levels(colLabels(sce))) %>%
    enframe(name = "cell") %>%
    unnest(value)
  
}

#' Helper Function to join DEG stats to LR
#'
#' @param cluster_markers dataframe with DE stats for a cluster
#' @param entity Transmitter or Receiver vector passed as tibble
#' @param source_target whether this is the source or target cluster
#'
#' @return A tibble with stats for receivers or transmitters per cluster
#'
#' @noRd
ligrec_degformat <- function(cluster_markers,
                             entity,
                             source_target){
  cluster_markers %>%
    left_join(entity, ., by = "gene") %>%
    na.omit() %>%
    {
      if(source_target=="source"){
        dplyr::select(
          .,
          ligand = gene,
          ligand.pval = p.value,
          ligand.FDR = FDR,
          ligand.stat = stat
        )
      }else if(source_target=="target"){
        dplyr::select(
          .,
          receptor = gene,
          receptor.pval = p.value,
          receptor.FDR = FDR,
          receptor.stat = stat
        )
      } else{
        stop("Incorrect entity!")
      }
    } %>%
    distinct() %>%
    na.omit()
}
#' Join Expression per Cluster
#'
#' @param lr_res LR formatted DE results from \link{ligrec_degformat}
#' @param means Gene avg expression per cluster
#' @param source_target target or source cell
#' @param entity ligand or receptor
#' @param type type of mean to join (count or scaled)
#'
#' @importFrom magrittr %>% %<>%
#'
#' @return Returns the Average Expression Per Cluster
#'
#' @noRd
join_means <- function(lr_res,
                       means,
                       source_target,
                       entity,
                       type,
                       pb = NULL){
  
  if(!is.null(pb)){
    pb$tick()$print()
  }
  
  entity.avg <- sym(str_glue("{entity}.{type}"))
  
  means %<>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "avg") %>%
    dplyr::rename({{ source_target }} := cell,
                  {{ entity }} := gene,
                  {{ entity.avg }} := avg)
  
  lr_res %>%
    left_join(means, by=c(source_target, entity))
}


#' Join Expression per Cluster
#'
#' @param lr_res LR formatted DE results from \link{ligrec_degformat}
#' @param means Gene avg expression per cluster
#' @param entity ligand or receptor
#'
#' @return Returns the Summed Average Expression Per Cluster
#'
#' @noRd
join_sum_means <- function(lr_res, means, entity){
  
  entity.expr = sym(str_glue("{entity}.sum"))
  
  sums <- means %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    rowwise("gene") %>%
    mutate(sum.means = sum(c_across(where(is.numeric)))) %>%
    select(gene, sum.means) %>%
    dplyr::rename({{ entity }} := gene,
                  {{ entity.expr }} := sum.means)  %>%
    distinct() %>%
    ungroup()
  
  lr_res %>%
    left_join(sums, by=c(entity))
}


#' Helper Function to join log2FC dataframe to LR_res
#'
#' @param lr_res LR formatted DE results from \link{ligrec_degformat}
#' @param logfc_df obtained via \link{get_log2FC}
#' @param source_target target or source cell
#' @param entity ligand or receptor
#'
#' @noRd
join_log2FC <- function(lr_res,
                        logfc_df,
                        source_target,
                        entity){
  
  entity.fc = sym(str_glue("{entity}.log2FC"))
  
  logfc <- logfc_df %>%
    dplyr::rename(
      {{ source_target }} := cell,
      {{ entity }} := gene,
      {{ entity.fc }} := avg_log2FC)  %>%
    distinct()
  
  lr_res %>%
    left_join(logfc, by=c(entity, source_target))
  
}


#' Helper Function to 'decomplexify' ligands and receptors into individual subunits
#'
#' @param resource a ligrec resource
#'
#' @param columns columns to separate and pivot long (e.g. genesymbol or uniprot),
#' `source_genesymbol` and `target_genesymbol` by default
#'
#' @return returns a longer tibble with complex subunits on seperate rows
#'
#' @details takes any number of columns, and assumes `_` as sep.
#'
#' @export
decomplexify <- function(resource,
                         columns = c("source_genesymbol",
                                     "target_genesymbol")){
  columns %>%
    map(function(col){
      sep_cols <- c(str_glue("col{rep(1:5)}"))
      col.complex <- str_glue("{col}_complex")
      
      resource <<- resource %>%
        mutate({{ col.complex }} :=
                 resource[[str_glue("{col}")]]) %>%
        separate(col,
                 into = sep_cols,
                 sep = "_",
                 extra = "drop",
                 fill = "right") %>%
        pivot_longer(cols = all_of(sep_cols),
                     values_to = col,
                     names_to = NULL) %>%
        tidyr::drop_na(col) %>%
        distinct() %>%
        mutate_at(.vars = c(col),
                  ~str_replace(., "COMPLEX:", ""))
    })
  return(resource)
}


test.liana.wrap.step.by.step = function()
{
  ## run liana_wrap step by step:
  sce <- liana_prep(sce,
                    idents_col = 'celltype',
                    verbose = verbose)
  
  
  method = 'cellphonedb'; resource = 'CellPhoneDB'; verbose = TRUE; assay.type = 'logcounts'
  method %<>% stringr::str_to_lower()
  
  resource %<>% select_resource # if null OmniPath
  
  ####### !! did not work 
  lr_resuits = liana_pipe(sce, op_resource = resource$CellPhoneDB, assay = 'integrated', assay.type = 'logcounts')
  
  # hack the function liana_pipe
  op_resource = resource
  ### this whole chunk needs to move to liana_wrap
  # Resource Format
  transmitters <- op_resource$CellPhoneDB$source_genesymbol %>%
    as_tibble() %>%
    select(gene = value)
  receivers <- op_resource$CellPhoneDB$target_genesymbol %>%
    as_tibble() %>%
    select(gene = value)
  entity_genes = union(transmitters$gene,
                       receivers$gene)
  
  # calculate global_mean required for SCA
  global_mean <- fast_mean(exec(assay.type, sce))
  
  # Filter `sce` to only include ligand receptor genes
  # and any cells which don't contain any expressed LR genes
  sce <- sce[toupper(rownames(sce)) %in% entity_genes,
             Matrix::colSums(counts(sce)) > 0]
  
  source('LIANA_helperFunctions.R')
  # Scale genes across cells
  sce@assays@data[["scaledata"]] <- row_scale(exec(assay.type, sce))
  
  source('LIANA_helperFunctions.R')
  # calculate PEM scores
  pem_scores <- compute_pem_scores(sce = sce, assay.type = assay.type)
  
  # Get Log2FC
  logfc_df <- get_log2FC(sce, "logcounts")
  
  test.type = "wilcox"
  pval.type = "all"
  # Find Markers and Format
  cluster_markers <- scran::findMarkers(sce,
                                        groups = colLabels(sce),
                                        direction = "any",
                                        full.stats = TRUE,
                                        test.type = test.type,
                                        pval.type = pval.type,
                                        assay.type = assay.type) %>%
    pluck("listData") %>%
    map2(., names(.), function(cluster, cluster_name){
      cluster %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>%
        as_tibble() %>%
        select(gene, p.value, FDR, stat = summary.stats)
    })
  
  # Get all Possible Cluster pair combinations
  pairs <- expand_grid(source = unique(colLabels(sce)),
                       target = unique(colLabels(sce)))
  
  # Get DEGs to LR format
  lr_res <- pairs %>%
    pmap(function(source, target){
      source_stats <- ligrec_degformat(cluster_markers[[source]],
                                       entity = transmitters,
                                       source_target = "source")
      target_stats <- ligrec_degformat(cluster_markers[[target]],
                                       entity = receivers,
                                       source_target = "target")
      
      op_resource %>%
        select(ligand = source_genesymbol,
               receptor = target_genesymbol) %>%
        left_join(source_stats, by = "ligand") %>%
        left_join(target_stats, by = "receptor") %>%
        na.omit() %>%
        distinct() %>%
        mutate(source = source,
               target = target)
    }) %>%
    bind_rows()
  
  
  source('LIANA_helperFunctions.R')
  
  # Join Expression Means
  lr_res %<>%
    join_means(means = means,
               source_target = "source",
               entity = "ligand",
               type = "expr") %>%
    join_means(means = means,
               source_target = "target",
               entity = "receptor",
               type = "expr") %>%
    join_means(means = scaled,
               source_target = "source",
               entity = "ligand",
               type = "scaled") %>%
    join_means(means = props,
               source_target = "target",
               entity = "receptor",
               type = "prop") %>%
    join_means(means = props,
               source_target = "source",
               entity = "ligand",
               type = "prop") %>%
    join_means(means = scaled,
               source_target = "target",
               entity = "receptor",
               type = "scaled") %>%
    
    join_sum_means(means = means,
                   entity = "ligand") %>%
    join_sum_means(means = means,
                   entity = "receptor") %>%
    # Join PEM scores
    join_means(means = pem_scores,
               source_target = "target",
               entity = "receptor",
               type = "pem") %>%
    join_means(means = pem_scores,
               source_target = "source",
               entity = "ligand",
               type = "pem") %>%
    # logFC
    join_log2FC(logfc_df,
                source_target = "source",
                entity="ligand") %>%
    join_log2FC(logfc_df,
                source_target = "target",
                entity="receptor") %>%
    # Global Mean
    mutate(global_mean = global_mean)
  
  liana_message("LIANA: LR summary stats calculated!",
                verbose = verbose
  )
  
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
  
  # Resource currently included in OmniPathR (and hence `liana`) include:
  show_resources()
  
  # Resource currently included in OmniPathR (and hence `liana`) include:
  show_methods()
  
  liana_path <- system.file(package = "liana")
  testdata <- readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))
  
  testdata %>% glimpse()
  
  sels =c(which(refs$celltype == 'CM')[1:1000], 
          which(refs$celltype == 'FB')[1:1000], 
          which(refs$celltype == 'prolife.Mphage'))
  
  subref = subset(refs, cells = colnames(refs)[sels])
  
  #rownames(subref) = toupper(rownames(subref))
  
  Idents(subref) = subref$celltype
  
  # Run liana
  # liana_test <- liana_wrap(testdata, method = 'cellphonedb', resource = 'CellPhoneDB')
  
  sce <- as.SingleCellExperiment(subref)
  colLabels(sce) = as.factor(sce$celltype)
  rownames(sce) = toupper(rownames(sce))
  
  raw = counts(sce)
  
  expr = logcounts(sce)
  #logcounts(sce) = exp(expr)
  ss1 = apply(raw, 1, sum)
  ss2 = apply(expr, 1, sum)
  cc1 = apply(raw, 2, sum)
  cc2 = apply(expr, 2, sum)
  
  sce = sce[which(ss1>10 & ss2> 10), which(cc1>10 & cc2 > 0)]
  
  ## working now !!!!
  liana_test <- liana_wrap(sce,  
                           resource = c("Consensus", 'CellPhoneDB',  "CellChatDB",  "CellTalkDB"), 
                           assay.type = "logcounts")
  
  #> Warning in .filter_sce(sce): 3465 genes and/or 0 cells were removed as they had
  #> no counts!
  
  # Liana returns a list of results, each element of which corresponds to a method
  liana_test %>% glimpse
  
  # We can aggregate these results into a tibble with consensus ranks
  liana_test <- liana_test %>%
    liana_aggregate(resource = 'Consensus')
  
  liana_test %>% 
    filter(source =="FB") %>%
    top_n(n = 25, wt = aggregate_rank) %>%
    liana_dotplot(source_groups = c("FB"),
                  target_groups = c("CM", "prolife.Mphage"))
  
  
  liana_trunc <- liana_test %>%
    # only keep interactions concordant between methods
    filter(aggregate_rank <= 0.01) # this can be FDR-corr if n is too high
  
  heat_freq(liana_trunc)
  
  # run any method of choice
  # Load Sce testdata
  # sce <- readRDS(file.path(liana_path , "testdata", "input", "testsce.rds"))
  
  # RUN CPDB alone
  # cpdb_test <- liana_wrap(sce,
  #                         method = 'cellphonedb',
  #                         resource = c('CellPhoneDB'),
  #                         permutation.params = list(nperms=100,
  #                                                   parallelize=FALSE,
  #                                                   workers=4))
  # 
  # # Plot toy results
  # cpdb_test %>%
  #   # filter(pvalue <= 0.05) %>% # only keep interactions with p-val <= 0.05
  #   # invert size (low p-value/high specificity = larger dot size)
  #   # + add a small value to avoid Infinity for 0s
  #   mutate(pvalue = -log10(pvalue + 1e-10)) %>% 
  #   liana_dotplot(source_groups = c("c"),
  #                 target_groups = c("c", "a", "b"),
  #                 specificity = "pvalue",
  #                 magnitude = "lr.mean",
  #                 show_complex = TRUE)
  # 
  # # Run liana re-implementations with the CellPhoneDB resource
  # complex_test <- liana_wrap(testdata,
  #                            method = c('natmi', 'sca', 'logfc'),
  #                            resource = c('CellPhoneDB'))
  # #> Warning in .filter_sce(sce): 3465 genes and/or 0 cells were removed as they had
  # #> no counts!
  # 
  # complex_test %>% liana_aggregate()
  # 
  
  
}
