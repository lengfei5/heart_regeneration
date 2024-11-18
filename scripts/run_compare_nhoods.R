##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: test the neighbourhood analysis 
# original code from https://github.com/lengfei5/RabbitGastrulation2022/blob/master/
# 8-compare_species/compare_nhoods.ipynb
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Oct 30 11:11:52 2024
##########################################################################
##########################################################################
rm(list = ls())

# Load packages
suppressPackageStartupMessages(library(scrabbitr))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(miloR))
suppressPackageStartupMessages(library(DelayedArray))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggraph))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(ggalluvial))
suppressPackageStartupMessages(library(ggrepel))
library(scuttle)
library(Seurat)
source('functions_scRNAseq.R')
source('functions_Visium.R')
source('functions_cccInference.R')


version.analysis = '_20240220'

resDir = paste0("../results/cross_species", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

outDir = paste0(resDir, '/neighborhoods_correlations/')
system(paste0('mkdir -p ', outDir))

##########################################
# load ax and mm data and process them 
##########################################
ax = readRDS(file = paste0(RdataDir, 'ax_scRNAseq.rds'))
nm = readRDS(file = paste0(RdataDir, 'nm_scRNAseq.rds'))

ax = as.SingleCellExperiment(ax, assay = 'RNA')
nm = as.SingleCellExperiment(nm, assay = 'integrated')

reducedDim(nm,"PCA") <- reducedDim(nm, "MNN")
xx = logcounts(nm)
xx[which(xx<0)] = 0
logcounts(nm) = xx

ax$celltype = ax$celltypes
nm$celltype = nm$subtype
# # Load rabbit data
# dataDir = '../published_dataset/scripts/'
# 
# r_data <- readRDS(paste0(dataDir, 'r_sce.rds'))
# r_data
# 
# m_data = readRDS(paste0('../../RA_competence/results/dataset_scRNAseq_MouseGastrulationData/Rdata/', 
#                      'EmbryoAtlasData_all36sample.rds'))
# m_data
# 
# m_data <- logNormCounts(m_data)
# assayNames(m_data)
# sels = which(!is.na(reducedDim(m_data,"pca.corrected")[,1]))
# m_data = m_data[,sels]
# 
# # check mouse rowData
# head(rowData(m_data))
# head(rowData(r_data))

##########################################
# prepare the one-to-one orthologs 
##########################################
# Load rabbit-mouse one-to-one orthologs
# library(biomaRt)
# listEnsembl()
# ensembl <- useEnsembl(biomart = "genes")
# ensembl
# datasets <- listDatasets(ensembl)
# head(datasets)
# 
# rm_orthologs = getEnsemblHomologs(ref_species = 'Oryctolagus cuniculus',
#                                   ref_genes = rownames(r_data),
#                                   query_species = 'mmusculus')

# rm_orthologs <- read.csv(paste0(dataDir, "mammalian_orthologs_20231113.csv"))
# rm_orthologs[1:5, ]
# 
# jj = which(colnames(rm_orthologs) == "Rabbit_EnsemblID" | colnames(rm_orthologs) == "Mouse_EnsemblID")
# rm_orthologs = rm_orthologs[,jj]
# rm_orthologs = rm_orthologs[, c(2, 1)]
# 
# rm_orthologs = rm_orthologs[which(!is.na(rm_orthologs[, 1]) & !is.na(rm_orthologs[, 2])),]
# colnames(rm_orthologs) = c('ref', 'query')
# 
# rm_orthologs = rm_orthologs[match(unique(rm_orthologs$ref), rm_orthologs$ref), ]
# 
# rownames(rm_orthologs) = rm_orthologs$ref
# rm_orthologs[1:5, ]

an_orthologs = data.frame(ref = rownames(ax), query = rownames(ax))
rownames(an_orthologs) = an_orthologs$ref
an_orthologs$query = sapply(an_orthologs$query, 
                            function(x){firstup(unlist(strsplit(as.character(x), '-'))[1])})

jj = which(!is.na(match(an_orthologs$query, rownames(nm))))
an_orthologs = an_orthologs[jj, ]

counts = table(an_orthologs$query)
gg_uniq = names(counts)[which(counts == 1)]
jj2 = which(!is.na(match(an_orthologs$query, gg_uniq)))

an_orthologs = an_orthologs[jj2, ]


##########################################
# compute the milo graph
##########################################
reLoad_computated_r_milo_m_miol = TRUE
if(reLoad_computated_r_milo_m_miol)
{
  r_milo = readRDS(paste0('../published_dataset/scripts/Rdata/', "r_milo.rds"))
  m_milo = readRDS(paste0('../published_dataset/scripts/Rdata/', "m_milo.rds"))
  
}else{
  # Compute rabbit neighbourhoods
  a_milo <- Milo(ax)
  
  a_milo <- buildGraph(a_milo, k=30, d=50, reduced.dim="PCA")
  
  a_milo <- makeNhoods(a_milo, prop=0.05, k=30, d=50, refined=T, reduced_dims="PCA")
  
  a_milo <- buildNhoodGraph(a_milo)
  
  a_milo
  
  # Export miloR object
  writeMM(a_milo@nhoods, paste0(outDir, "r_nhoods.mtx"))
  saveRDS(a_milo, paste0(outDir, "r_milo.rds"))
  
  
  options(repr.plot.width = 12, repr.plot.height = 5, repr.plot.res = 300)
  p1 <- scrabbitr::plotNhoodSizeHist(a_milo, colour="blue")
  
  p1
  
  p2 <- plotNhoodGraph(a_milo, size_range=c(0.1,3), node_stroke=0.1) + 
    scale_fill_viridis(name = "Nhood size", option = "viridis", direction = 1) 
  #grid.arrange(p1, p2, nrow=1)
  
  p1  + p2
  ggsave(paste0(outDir, "ax_milo_nhood_size_hist_size_graph.pdf"), width=8, height=4, dpi=300)
  
  # Compute mouse neighbourhoods
  n_milo <- Milo(nm)
  
  n_milo <- buildGraph(n_milo, k=30, d=50, reduced.dim="PCA")
  
  n_milo <- makeNhoods(n_milo, prop=0.05, k=30, d=50, refined=T, reduced_dims="PCA")
  
  n_milo <- buildNhoodGraph(n_milo)
  
  
  # Plot nhoods and size distribution
  options(repr.plot.width = 12, repr.plot.height = 5, repr.plot.res = 300)
  
  p1 <- scrabbitr::plotNhoodSizeHist(n_milo, colour="red")
  #ggsave("../plots/compare_nhoods/m_milo_nhood_size_hist.pdf", p1, width=4, height=4, dpi=300)
  
  p2 <- plotNhoodGraph(n_milo, size_range=c(0.1,3),node_stroke=0.1) + 
    scale_fill_viridis(name = "Nhood size", option = "viridis", direction=1)
  #ggsave("../plots/compare_nhoods/m_nhood_size_graph.pdf", p2, width=6, height=5, dpi=300)
  
  #grid.arrange(p1, p2, nrow = 1)
  p1  + p2
  ggsave(paste0(outDir, "nm_milo_nhood_size_hist_size_graph.pdf"), width=8, height=4, dpi=300)
  
  saveRDS(n_milo, paste0(outDir, "n_milo.rds"))
  
}

#reducedDim(m_milo,"UMAP") <- reducedDim(m_milo,"umap")
#reducedDim(m_milo,"PCA") <- reducedDim(m_milo,"pca.corrected")

# check mouse rowData
head(rowData(a_milo))
head(rowData(n_milo))

# Add mouse colData
head(colData(a_milo))
#table(m_milo$stage)
#table(m_milo$dissection)
#table(m_milo$somite_count)

head(colData(n_milo))
#table(r_milo$somite_count)
#table(r_milo$dissection)
head(n_milo$celltype)
#a_milo$celltype = a_milo$celltypes
#n_milo$celltype = n_milo$subtype

## Run neighbourhood comparison pipeline
# Run pipeline
out <- scrabbitr::calcNhoodSim(a_milo, n_milo, an_orthologs, 
                               sim_preprocessing="gene_spec", 
                               sim_measure="pearson",
                               hvg_join_type="intersection", 
                               max_hvgs=2000, 
                               #r_exclude = r_exclude, 
                               #m_exclude = m_exclude,
                               export_dir = outDir, 
                               verbose = TRUE)


# Load exported results
#nhood_sim <- as.matrix(fread(paste0(dataDir, "test_out/", "nhood_sim.tsv"), sep="\t"), rownames=1)

#r_vals <- fread("../data-out/compare_nhoods/nhood_comparisons/intersect_2000_pearson_w_gspec_exclude_cc_mito/r_vals.tsv", sep="\t")
#m_vals <- fread("../data-out/compare_nhoods/nhood_comparisons/intersect_2000_pearson_w_gspec_exclude_cc_mito/m_vals.tsv", sep="\t")
#out <- list(r_vals = r_vals, m_vals = m_vals, nhood_sim = nhood_sim)
#nhood_sim[1:5,1:5]

out$nhood_sim[1:5, 1:5]

# Extract neighbourhood graph
a_graph <- nhoodGraph(a_milo)
n_graph <- nhoodGraph(n_milo)


# Add nhood attributes to igraph
a_nhoodIDs <- as.numeric(vertex_attr(a_graph)$name) 
a_indCells <- colnames(a_milo)[a_nhoodIDs]

V(a_graph)$cell_name <- a_indCells
V(a_graph)$celltype <- colData(a_milo)[a_indCells, "celltype"]

n_nhoodIDs <- as.numeric(vertex_attr(n_graph)$name) 
n_indCells <- colnames(n_milo)[n_nhoodIDs]

V(n_graph)$cell_name <- n_indCells
V(n_graph)$celltype <- colData(n_milo)[n_indCells, "celltype"]


# Calculate maximum correlations  
a_maxNhoods <- getMaxMappings(out$nhood_sim, 1, long_format=FALSE) # rabbit-mouse
n_maxNhoods <- getMaxMappings(out$nhood_sim, 2, long_format=FALSE) # mouse-rabbit
df_simFilt <- rbind(a_maxNhoods, n_maxNhoods)                                                  

options(repr.plot.width = 18, repr.plot.height = 8, repr.plot.res = 300)
p1 <- plotNhoodMaxSim(a_milo, a_maxNhoods)
p2 <- plotNhoodMaxSim(n_milo, n_maxNhoods)
p1 + p2

ggsave(paste0(outDir, "n_milo_a_milo_max_corr.pdf"), width=10, height=4, dpi=300)


##########################################
# highlight the distribution of correlation of celltypes
##########################################
celltypes <- unique(c(unique(colData(a_milo)$celltype), unique(colData(n_milo)$celltype)))
cols =  randomcoloR::distinctColorPalette(k = length(celltypes), altCol = FALSE, runTsne = FALSE)
cols = cols[length(cols):1]
names(cols) = celltypes

# Highlight 
options(repr.plot.width = 5.2, repr.plot.height = 5, repr.plot.res = 300)
#a_milo$isCM <- ifelse(a_milo$celltypes %in% exe_celltypes,"Extra-embryonic", "Embryonic")

p <- plotNhoodSimGroups(a_milo, a_maxNhoods$sim, 
                        group_by = "celltypes", 
                        xlabel="Correlation", ylabel="subtypes", 
                        colour_by="celltypes", 
                        group_colours = cols,
                        size=0.15, 
                        rel_min_height=0.001, show_rank = FALSE
                        )

p <- p + 
  theme(text = element_text(size=14),
        axis.text = element_text(size=12), 
        axis.text.y = element_text(size=6),
        axis.ticks = element_line(size = 1),
        panel.grid.minor = element_line(size = 0.1), 
        panel.grid.major = element_line(size = 0.2))  +
  NoLegend()
p

ggsave(paste0(outDir, "a_nhoods_celltype_ranked_.pdf"), p, 
       width=5, height=8, dpi=300)


# ## similiary, we can look at the distribution of max correlation across samples or stages
# options(repr.plot.width = 6, repr.plot.height = 2.5, repr.plot.res = 300)
# 
# sample_pal <- scrabbitr::getDiscreteColours(unique(a_milo$sample))
# names(sample_pal) <- unique(a_milo$sample)
# 
# p <- scrabbitr::plotNhoodSimGroups(a_milo, a_maxNhoods$sim, xlabel="Correlation", 
#                                    group_by="sample", group_colours = sample_pal,
#                                    ylabel="Sample", type="boxplot", 
#                                    orientation = "horizontal", size=0.1, rel_min_height=0.001) + 
#   theme(axis.text.y = element_text(size=12),
#         axis.text.x = element_text(size=12)) + theme_linedraw()
# p
# ggsave(paste0(outDir, "sample_boxplot.pdf"), p, 
#        width=6, height=2.5, dpi=300)
# 
# options(repr.plot.width = 4, repr.plot.height = 2, repr.plot.res = 300)
# 
# p <- plotNhoodSimGroups(a_milo, a_maxNhoods$sim, xlabel="Correlation", 
#                         group_by="stage",
#                         group_colours=scrabbitr::getStageColours(), ylabel="Stage", 
#                         size=0.1, rel_min_height=0.001) + 
#   theme(axis.text.y = element_text(size=12),
#         axis.text.x = element_text(size=12)) + theme_linedraw()
# p
# 
# ggsave(paste0(outDir, "stage_ridgeline.pdf"), p, width=4, height=2, dpi=300)

##########################################
# hightlight the projection between species
##########################################
options(repr.plot.width = 18, repr.plot.height = 8, repr.plot.res = 300)

celltypes <- unique(c(unique(colData(a_milo)$celltype), unique(colData(n_milo)$celltype)))
#celltypes = celltypes[grep('FB|CM|EC', celltypes)]

p_all <- plotTrajMappings(a_milo, n_milo, df_simFilt, 
                         group="celltype", 
                         groups=celltypes, 
                         dimred="UMAP", 
                         #colour_by = cols, 
                         #rotate=90,
                         offset=c(0, 50), reflect.X=FALSE, reflect.Y=FALSE,
                         line_alpha=0.02,
                         edge_alpha=0.01, 
                         legend_pos="right") +
  guides(fill = guide_legend(override.aes = list(size=6))) + ggtitle("all")

p_all

ggsave(paste0(outDir, "axolotl_nm_mapping_test_v10.pdf"), width=18, height=10, dpi=300)


#### CM only
celltypes <- unique(colData(a_milo)$celltypes)
cols =  randomcoloR::distinctColorPalette(k = length(celltypes), altCol = FALSE, runTsne = FALSE)
cols = cols[length(cols):1]
names(cols) = celltypes

CM_celltypes <- c(unique(colData(a_milo)$celltypes), unique(colData(n_milo)$celltype))
CM_celltypes = CM_celltypes[grep('CM', CM_celltypes)]

p_cm <- plotTrajMappings(a_milo, n_milo, df_simFilt, 
                          group="celltype", 
                          groups=CM_celltypes, 
                          dimred="UMAP", 
                          #colour_by = 'celltypes',
                          offset=c(10,10),line_alpha=0.2,edge_alpha=0.01, 
                          reflect.X=FALSE, legend_pos="right") +
  guides(fill = guide_legend(override.aes = list(size=6))) + ggtitle("CM")

p_cm

ggsave(paste0(outDir, "gut_mapping.pdf"), p_gut, width=18, height=8, dpi=300)

##########################################
# projection 
##########################################
meso_celltypes <- c("Sclerotome","Dermomyotome", "Somitic mesoderm", "Anterior somitic tissues",
                    "Posterior somitic tissues","Cranial mesoderm","Anterior Primitive Streak",
                    "Caudal epiblast","Paraxial mesoderm","Intermediate mesoderm","Epiblast",
                    "Primitive Streak","Presomitic mesoderm","NMPs","NMPs/Mesoderm-biased",
                    "Caudal mesoderm", "Nascent mesoderm")
r_celltypes = unique(r_milo$celltype)
m_celltypes = unique(m_milo$celltype)
match(meso_celltypes, r_celltypes)
match(meso_celltypes, m_celltypes)

# Mesoderm
EC_celltypes <- c(unique(colData(a_milo)$celltypes), unique(colData(n_milo)$celltype))
EC_celltypes = EC_celltypes[grep('EC', EC_celltypes)]

p_ec <- plotTrajMappings(a_milo, n_milo, df_simFilt, 
                           group="celltype", 
                           groups=EC_celltypes, 
                           dimred="UMAP", offset=c(20,20), reflect.X=FALSE,
                           line_alpha=0.1,edge_alpha=0.01, 
                           legend_pos="right") + guides(fill = guide_legend(override.aes = list(size=6))) +
  ggtitle("EC")

p_ec

ggsave(paste0(outDir, "EC_mapping.pdf"), width=18, height=8, dpi=300)

#ggsave("../plots/compare_nhoods/blood_mapping.pdf", p_blood, width=18, height=8, dpi=300)

# Blood
FB_celltypes <- c(unique(colData(a_milo)$celltypes), unique(colData(n_milo)$celltype))
FB_celltypes = FB_celltypes[grep('FB', FB_celltypes)]

p_fb <- plotTrajMappings(a_milo, n_milo, df_simFilt, 
                            group="celltype", groups=FB_celltypes, 
                            dimred="UMAP", offset=c(15,-15),
                         line_alpha=0.03,edge_alpha=0.01, 
                            reflect.X=FALSE, legend_pos="right") + 
  guides(fill = guide_legend(override.aes = list(size=6))) + ggtitle("FB")

p_fb

ggsave(paste0(outDir, "FB_mapping.pdf"), width=18, height=8, dpi=300)


#ggsave("../plots/compare_nhoods/blood_mapping.pdf", p_blood, width=18, height=8, dpi=300)

all_celltypes <- c(unique(colData(a_milo)$celltypes), unique(colData(n_milo)$celltype))

p_all <- plotTrajMappings(a_milo, n_milo, df_simFilt, 
                         group="celltype", groups=all_celltypes, 
                         dimred="UMAP", offset=c(0,0),
                         line_alpha=0.03,edge_alpha=0.01, 
                         reflect.X=TRUE, legend_pos="right") + 
  guides(fill = guide_legend(override.aes = list(size=6))) + ggtitle("all")

p_all

ggsave(paste0(outDir, "all_mapping.pdf"), width=18, height=12, dpi=300)


amnion_celltypes <- c("Epiblast","Primitive Streak", 'Non-neural ectoderm 2', "Non-neural ectoderm 4",
                      "Non-neural ectoderm 1", "Non-neural ectoderm 3","Non-neural ectoderm 5",
                      "Placodal ectoderm", "Amnion 3", "Surface ectoderm", "Pharyngeal endoderm")

p_amnion <- plotTrajMappings(a_milo, n_milo, df_simFilt, 
                             group="celltype", groups=amnion_celltypes, 
                             dimred="UMAP",
                             offset=c(0,0),line_alpha=0.1,edge_alpha=0.01, 
                             reflect.X=FALSE, legend_pos="right") +
  guides(fill = guide_legend(override.aes = list(size=6))) + ggtitle("Amnion")



#ggsave("../plots/compare_nhoods/amnion_mapping.pdf", p_amnion, width=18, height=8, dpi=300)

allantois_celltypes <- c('Epiblast', "Primitive Streak", "Nascent mesoderm","Lateral plate mesoderm", 
                         "Allantois")

p_allantois <- plotTrajMappings(a_milo, n_milo, df_simFilt, 
                                group="celltype", groups=allantois_celltypes, 
                                dimred="FA", offset=c(0.3,0.5),
                                reflect.X=F,reflect.Y=T,line_alpha=0.1,
                                edge_alpha=0.01, legend_pos="right") +
  guides(fill = guide_legend(override.aes = list(size=6))) + ggtitle("Allantois")

p_allantois
#ggsave("../plots/compare_nhoods/allantois_mapping_fa.pdf", p_allantois, width=18, height=8, dpi=300)

# options(repr.plot.width = 18, repr.plot.height = 14, repr.plot.res = 300)
# grid.arrange(p_gut, p_meso, p_amnion, p_allantois, nrow=2)
# 
# 
# a_meso <- scrabbitr::subsetMiloGroups(a_milo, "celltype", meso_celltypes)
# n_meso <- scrabbitr::subsetMiloGroups(n_milo, "celltype", meso_celltypes)
# 
# a_meso <- runUMAP(a_meso, n_neighbors=50, min_dist=0.2, dimred = "PCA")
# n_meso <- runUMAP(n_meso, n_neighbors=50, min_dist=0.2, dimred = "pca.corrected")
# 
# options(repr.plot.width = 12, repr.plot.height = 6.5, repr.plot.res = 300)
# 
# p <- plotNhoodMappings(a_meso, n_meso, df_simFilt, "UMAP", "celltype",
#                        offset=c(0,15), reflect.X=TRUE, line_alpha=0.1, edge_alpha=0, legend_pos="right") +
#   guides(fill = guide_legend(override.aes = list(size=1))) + ggtitle("Mesoderm")
# 
# p
# 
# ggsave(paste0(outDir, "mapping_embeddings_recomputed_knngrpah.pdf"), p, width=16, height=12, dpi=300)

options(repr.plot.width = 12, repr.plot.height = 6.5, repr.plot.res = 300)

p <- plotNhoodMappings(a_milo, n_milo, df_simFilt, dimred = "UMAP", 
                       offset=c(0,15), reflect.X=TRUE, line_alpha=0.1, edge_alpha=0, legend_pos="right") +
  guides(fill = guide_legend(override.aes = list(size=1))) + ggtitle("Mesoderm")

p

