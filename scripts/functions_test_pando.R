# Get variables from object
library(sparseMatrixStats)
object = srat_cr_plus

peak_to_gene_method = c('Signac')
upstream = 10e+6; downstream = 10e+6
#extend = 1000000,
only_tss = FALSE
peak_to_gene_domains = NULL
parallel = FALSE
tf_cor = 0.1
peak_cor = 0.
aggregate_rna_col = NULL
aggregate_peaks_col = NULL
method = c('glm', 'glmnet', 'cv.glmnet', 'brms', 'xgb', 'bagging_ridge', 'bayesian_ridge')
interaction_term = ':'
adjust_method = 'fdr'

params <- Params(object)
motif2tf <- NetworkTFs(object)

if (is.null(motif2tf)){
  stop('Motif matches have not been found. Please run `find_motifs()` first.')
}

gene_annot <- Signac::Annotation(object[[params$peak_assay]])

if (is.null(gene_annot)){
  stop('Please provide a gene annotation for the ChromatinAssay.')
}

# Select target genes for GRN inference
if (is.null(genes)){
  genes <- VariableFeatures(object, assay=params$rna_assay)
  if (is.null(genes)){
    stop('Please provide a set of features or run `FindVariableFeatures()`')
  }
}

# Get assay data or summary
if (is.null(aggregate_rna_col)){
  gene_data <- Matrix::t(Seurat::GetAssayData(object, assay=params$rna_assay))
  gene_groups <- TRUE
} else {
  gene_data <- GetAssaySummary(
    object,
    assay = params$rna_assay,
    group_name = aggregate_rna_col,
    verbose = FALSE
  )
  gene_groups <- object@meta.data[[aggregate_rna_col]]
}

if (is.null(aggregate_peaks_col)){
  peak_data <- Matrix::t(Seurat::GetAssayData(object, assay=params$peak_assay))
  peak_groups <- TRUE
} else {
  peak_data <- GetAssaySummary(
    object,
    assay = params$peak_assay,
    group_name = aggregate_peaks_col,
    verbose = FALSE
  )
  peak_groups <- object@meta.data[[aggregate_peaks_col]]
}

# Select genes to use by intersecting annotated genes with all
# detected genes in the object
features <- intersect(gene_annot$gene_name, genes) %>%
  intersect(rownames(GetAssay(object, params$rna_assay)))
gene_annot <- gene_annot[gene_annot$gene_name%in%features, ]

# Get regions
regions <- NetworkRegions(object)
peak_data <- peak_data[, regions@peaks]
colnames(peak_data) <- rownames(regions@motifs@data)
peaks2motif <- regions@motifs@data

# Find candidate regions near gene bodies
if (is.null(peak_to_gene_domains)){
  #log_message('Selecting candidate regulatory regions near genes', verbose=verbose)
  peaks_near_gene <- find_peaks_near_genes(
    peaks = regions@ranges,
    method = peak_to_gene_method,
    genes = gene_annot,
    upstream = upstream,
    downstream = downstream,
    only_tss = only_tss
  )
  
} else {
  log_message('Selecting candidate regulatory regions in provided domains', verbose=verbose)
  peaks_near_gene <- find_peaks_near_genes(
    peaks = regions@ranges,
    method = 'Signac',
    genes = peak_to_gene_domains,
    upstream = 0,
    downstream = 0,
    only_tss = FALSE
  )
}

peaks2gene <- aggregate_matrix(t(peaks_near_gene), groups=colnames(peaks_near_gene), fun='sum')

# Select peaks passing criteria
peaks_at_gene <- as.logical(colMaxs(peaks2gene))
peaks_with_motif <- as.logical(rowMaxs(peaks2motif*1))

# Subset data to good peaks
peaks_use <- peaks_at_gene & peaks_with_motif
peaks2gene <- peaks2gene[, peaks_use, drop=FALSE]
peaks2motif <- peaks2motif[peaks_use, , drop=FALSE]
peak_data <- peak_data[, peaks_use, drop=FALSE]

log_message('Preparing model input', verbose=verbose)
tfs_use <- colnames(motif2tf)
motif2tf <- motif2tf[, tfs_use, drop=FALSE]

log_message('Fitting models for ', length(features), ' target genes' , verbose=verbose)
# Loop through features and fit models/run CV for each
names(features) <- features
model_fits <- map_par(features, function(g){
  
  # g = features[1]
  # Select peaks near gene
  if (!g%in%rownames(peaks2gene)){
    log_message('Warning: ', g, ' not found in EnsDb', verbose=verbose==2)
    return()
  }
  gene_peaks <- as.logical(peaks2gene[g, ])
  
  if (sum(gene_peaks)==0){
    log_message('Warning: No peaks found near ', g, verbose=verbose==2)
    return()
  }
  
  # Select peaks correlating with target gene expression
  g_x <- gene_data[gene_groups, g, drop=FALSE]
  peak_x <- peak_data[peak_groups, gene_peaks, drop=FALSE]
  peak_g_cor <- as(sparse_cor(peak_x, g_x), 'generalMatrix')
  peak_g_cor[is.na(peak_g_cor)] <- 0
  peaks_use <- rownames(peak_g_cor)[abs(peak_g_cor[, 1]) > peak_cor]
  
  if (length(peaks_use)==0){
    log_message('Warning: No correlating peaks found for ', g, verbose=verbose==2)
    return()
  }
  peak_x <- peak_x[, peaks_use, drop=FALSE]
  peak_motifs <- peaks2motif[gene_peaks, , drop=FALSE][peaks_use, , drop=FALSE]
  
  # Select TFs with motifs in peaks
  gene_peak_tfs <- map(rownames(peak_motifs), function(p){
    # p = rownames(peak_motifs)[1]
    x <- as.logical(peak_motifs[p, ])
    peak_tfs <- colMaxs(motif2tf[x, , drop=FALSE])
    peak_tfs <- colnames(motif2tf)[as.logical(peak_tfs)]
    peak_tfs <- setdiff(peak_tfs, g)
    return(peak_tfs)
  })
  
  names(gene_peak_tfs) <- rownames(peak_motifs)
  
  # Check correlation of peaks with target gene
  gene_tfs <- purrr::reduce(gene_peak_tfs, union)
  tf_x <- gene_data[gene_groups, gene_tfs, drop=FALSE]
  tf_g_cor <- as(sparse_cor(tf_x, g_x), 'generalMatrix')
  tf_g_cor[is.na(tf_g_cor)] <- 0
  tfs_use <- rownames(tf_g_cor)[abs(tf_g_cor[, 1]) > tf_cor]
  if (length(tfs_use)==0){
    log_message('Warning: No correlating TFs found for ', g, verbose=verbose==2)
    return()
  }
  tf_g_corr_df <- as_tibble(
    tf_g_cor[unique(tfs_use), , drop=F],
    rownames='tf',
    .name_repair='check_unique'
  ) %>%
    rename('tf'=1, 'corr'=2)
  
  # Filter TFs and make formula string
  frml_string <- map(names(gene_peak_tfs), function(p){
    peak_tfs <- gene_peak_tfs[[p]]
    peak_tfs <- peak_tfs[peak_tfs%in%tfs_use]
    if (length(peak_tfs)==0){
      return()
    }
    peak_name <- str_replace_all(p, '-', '_')
    tf_name <- str_replace_all(peak_tfs, '-', '_')
    formula_str <- paste(
      paste(peak_name, interaction_term, tf_name, sep=' '), collapse = ' + ')
    return(list(tfs=peak_tfs, frml=formula_str))
  })
  frml_string <- frml_string[!map_lgl(frml_string, is.null)]
  if (length(frml_string)==0){
    log_message('Warning: No valid peak:TF pairs found for ', g, verbose=verbose==2)
    return()
  }
  
  target <- str_replace_all(g, '-', '_')
  model_frml <- as.formula(
    paste0(target, ' ~ ', paste0(map(frml_string, function(x) x$frml),  collapse=' + '))
  )
  
  # Get expression data
  nfeats <- sum(map_dbl(frml_string, function(x) length(x$tfs)))
  gene_tfs <- purrr::reduce(map(frml_string, function(x) x$tfs), union)
  gene_x <- gene_data[gene_groups, union(g, gene_tfs), drop=FALSE]
  model_mat <- as.data.frame(cbind(gene_x, peak_x))
  if (scale) model_mat <- as.data.frame(scale(as.matrix(model_mat)))
  colnames(model_mat) <- str_replace_all(colnames(model_mat), '-', '_')
  
  log_message('Fitting model with ', nfeats, ' variables for ', g, verbose=verbose==2)
  result <- try(fit_model(
    model_frml,
    data = model_mat,
    method = method,
    ...
  ), silent=TRUE)
  if (any(class(result)=='try-error')){
    log_message('Warning: Fitting model failed for ', g, verbose=verbose)
    log_message(result, verbose=verbose==2)
    return()
  } else {
    result$gof$nvariables <- nfeats
    result$corr <- tf_g_corr_df
    return(result)
  }
}, verbose=verbose, parallel=parallel)


model_fits <- model_fits[!map_lgl(model_fits, is.null)]
if (length(model_fits)==0){
  log_message('Warning: Fitting model failed for all genes.', verbose=verbose)
}

coefs <- map_dfr(model_fits, function(x) x$coefs, .id='target')
coefs <- format_coefs(coefs, term=interaction_term, adjust_method=adjust_method)
corrs <- map_dfr(model_fits, function(x) x$corr, .id='target')
if (nrow(coefs)>0){
  coefs <- suppressMessages(left_join(coefs, corrs))
}
gof <- map_dfr(model_fits, function(x) x$gof, .id='target')

params <- list()
params[['method']] <- method
params[['family']] <- family
params[['dist']] <- c('upstream'=upstream, 'downstream'=downstream)
params[['only_tss']] <- only_tss
params[['interaction']] <- interaction_term
params[['tf_cor']] <- tf_cor
params[['peak_cor']] <- peak_cor

network_obj <- new(
  Class = 'Network',
  features = features,
  coefs = coefs,
  fit = gof,
  params = params
)
object@grn@networks[[network_name]] <- network_obj
object@grn@active_network <- network_name
return(object)