##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: test Pando for GRN inference 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Mar 24 13:34:26 2023
##########################################################################
##########################################################################
rm(list = ls())
version.analysis = '_R13591_atac_reseq_20221115'

resDir = paste0("../results/sc_multiome", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../R14353_ax_snATAC_reseq'

#source('functions_scATAC.R')
#source('functions_scRNAseq.R')
#source('functions_Visium.R')
#source('functions_scRNAseq.R')

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(patchwork)
require(SeuratObject)
library(data.table)

get_geneName = function(xx)
{
  return(sapply(xx, function(x) {test = unlist(strsplit(as.character(x), '-')); 
  if(length(test) == 1) {
    test
  }else{
    test = test[-length(test)]; 
    paste0(test, collapse = '-')
  }
  }
  ))
}

library(pryr) # monitor the memory usage
require(ggplot2)
require(dplyr)
require(stringr)
require(tidyr)
require(tictoc)
library(future)
options(future.globals.maxSize = 80 * 1024^3)
set.seed(1234)
mem_used()

##########################################
# import annotation and metadata
##########################################
library(BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M)

library(ballgown)
gtf_axolotl = paste0("/groups/tanaka/People/current/jiwang/scripts/axolotl_multiome/r_package/", 
                     "AmexT_v47.FULL_corr_chr_cut.gtf")

granges_axolotl = ballgown::gffReadGR(gtf_axolotl)
# adding a gene biotype, as that's necessary for TSS metaprofile
granges_axolotl$gene_biotype = "protein_coding"

# with need to add the "proper" gene name
# basedir = "/links/groups/treutlein/USERS/tomasgomes/gene_refs/axolotl/Amex_T_v47/"
# gene_name_match = read.table(paste0(basedir, "AmexT_v47.FULL_t2g_note.txt"), sep = "\t")[,2:3]
# gene_name_match = gene_name_match[!duplicated(gene_name_match$V2), ]
# rownames(gene_name_match) = gene_name_match$V2
# newgenenames = gene_name_match[granges_axolotl$gene_id,2]
# granges_axolotl$gene_name = newgenenames

species = 'axloltl_scATAC'


########################################################
########################################################
# Section : start the data processing 
# 
########################################################
########################################################
srat_cr = readRDS(file = paste0(RdataDir, 
                                'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                'scATAC.merged.peaks.cr_filtered_lsi.umap_',
                                '584K.annot_38280cells.rds'))

Idents(srat_cr) = as.factor(srat_cr$celltypes)

srat_cr = subset(srat_cr, idents = 'CM')

## remove empty cells
count.mat = srat_cr@assays$ATAC@counts
cols.empty <- colSums(count.mat) == 0
srat_cr = subset(srat_cr, cells = colnames(srat_cr)[!cols.empty]) 

saveRDS(srat_cr, file = paste0(RdataDir, 
                                'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                'scATAC.merged.peaks.cr_filtered_lsi.umap_',
                                '584K.annot_38280cells_CM.test.pando.rds'))


##########################################
# original code from
# https://github.com/lengfei5/pallium_evo/blob/main/analysis/GRN_analysis/PandoBasic.r
##########################################
srat_cr = readRDS(file = paste0(RdataDir, 
                                'seuratObj_multiome_snRNA.annotated.normalized.umap_',
                                'scATAC.merged.peaks.cr_filtered_lsi.umap_',
                                '584K.annot_38280cells_CM.test.pando.rds'))

# libraries
library(Signac)
library(BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M)
library(Pando)
library(ggplot2)
library(igraph)

DefaultAssay(srat_cr) <- 'ATAC'

# plotting most important labels
DimPlot(srat_cr, reduction = "umap", group.by = "subtypes", label = T, repel = TRUE) + NoLegend()
DimPlot(srat_cr, reduction = "umap", group.by = "celltypes", label = T, repel = TRUE) + NoLegend()
#DimPlot(srat_cr, reduction = "umap", group.by = "pred_regions_top", label = T)


# Split ATAC into peaks from 3 types of regions CellRanger peaks
DefaultAssay(srat_cr) = "ATAC"
transc = Annotation(srat_cr)[Annotation(srat_cr)$type=="transcript",] # get transcripts
transc = granges(transc[width(transc)>600,])

# restrict fixes issues with negative genome coordinates
proms = restrict(promoters(transc, upstream=5000, downstream=500), start = 1) 
proms$type = "promoter"
genebodies = setdiff(transc, proms)
genebodies$type = "genebody"
feat = c(proms, genebodies)

closest_f = ClosestFeature(srat_cr, regions = rownames(srat_cr), annotation = feat)
closest_f$type = ifelse(closest_f$distance!=0, "distal", closest_f$type)

for(tt in unique(closest_f$type))
{
  p = closest_f$query_region[closest_f$type==tt]
  i = paste0("ATAC_", tt)
  
  srat_cr[[i]] = CreateChromatinAssay(srat_cr@assays$ATAC@counts[p,])
  
  srat_cr = RunTFIDF(srat_cr, assay = i)
  srat_cr = FindTopFeatures(srat_cr, assay=i, min.cutoff = 300)
  print(length(srat_cr[[i]]@var.features))
  
}

# Basic Pando run
DefaultAssay(srat_cr) = "RNA"

# Get motif data
data(motifs)

# define genes
genes = VariableFeatures(srat_cr, assay = "RNA")
#genes = unique(mk_ct_comp$feature)
length(genes)

# define peaks
srat_cr = FindTopFeatures(srat_cr, assay = "ATAC_promoter", min.cutoff = 450)
srat_cr = FindTopFeatures(srat_cr, assay = "ATAC_distal", min.cutoff = 450)
srat_cr = FindTopFeatures(srat_cr, assay = "ATAC_genebody", min.cutoff = 450)
regions = unique(c(srat_cr@assays$ATAC_promoter@var.features,
                   srat_cr@assays$ATAC_distal@var.features,
                   srat_cr@assays$ATAC_genebody@var.features))

isreg = rownames(srat_cr@assays$ATAC@meta.features) %in% regions
regions = srat_cr@assays$ATAC@ranges[isreg,]
length(regions)

# Initiate GRN object and select candidate regions
srat_cr_plus = initiate_grn(srat_cr, 
                            peak_assay = "ATAC", 
                            rna_assay = "RNA",
                            regions = regions)

# manual correct the motif-to-TFs mapping
utils::data(motif2tf, envir = environment())
motif2tf <- motif2tf %>% dplyr::select('motif'=1,'tf'=2)

DefaultAssay(srat_cr) = "RNA"
ggs = get_geneName(rownames(srat_cr))
ggs_tf = unique(intersect(ggs, motif2tf$tf))

motif_tfs = c()
for(n in 1:length(ggs_tf))
{
  # n = 1
  cat(n, '--', ggs_tf[n], '\n')
  jj = which(ggs == ggs_tf[n])
  kk = which(motif2tf$tf == ggs_tf[n])
  
  for(j in jj)
  {
    for(k in kk)
    {
      motif_tfs = rbind(motif_tfs, c(motif2tf$motif[k], rownames(srat_cr)[j]))
    }
  }
    
}

motif_tfs = data.frame(motif_tfs)
colnames(motif_tfs) = c('motif', 'tf')

# Scan candidate regions for TF binding motifs
srat_cr_plus = find_motifs(srat_cr_plus, 
                           pfm = motifs, 
                           motif_tfs = motif_tfs,
                           genome = BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M)

DefaultAssay(srat_cr_plus) = "ATAC"

# Infer gene regulatory network
library(doParallel)
library(dplyr)
library(tidyr) 
library(purrr)
cl = parallel::makeCluster(8)
doParallel::registerDoParallel(cl)

library(tictoc)
library(sparseMatrixStats)
#tic()
srat_cr_plus = infer_grn(srat_cr_plus, 
                         upstream = 10e+6, 
                         downstream = 10e+6,
                         only_tss = FALSE,
                         parallel = FALSE,
                         peak_to_gene_method = 'Signac')

toc()


# Print inferred coefficients
coef(srat_cr_plus)

# save
saveRDS(srat_cr_plus, 
        file = paste0(RdataDir, "multiome_integATAC_SCT_plus.RDS"))
write.csv(coef(srat_cr_plus), 
          file = paste0(dir, "results/multiome/srat_cr_plus_coef.csv"), 
          col.names = T, row.names = F, quote = F)


# plot global networks
pando_atac_coef = coef(srat_cr_plus)
pando_atac_coef = pando_atac_coef[pando_atac_coef$padj<=0.05,]

fff = pando_atac_coef$estimate>0.01 & pando_atac_coef$target %in% pando_atac_coef$tf
dat = as.data.frame(pando_atac_coef[fff,c(1,2,5,9)])

igdat = igraph::graph.data.frame(dat[,1:2])
xxx = igraph::betweenness(igdat, directed = T)
ppp = igraph::page_rank(igdat, directed = T)
ddd = igraph::degree(igdat, mode = "out")
eee = ego_size(igdat, mode = "out", order = 100)
iii = ego_size(igdat, mode = "in", order = 100)
plot(xxx, eee, pch = 21, cex = (ddd+1)/4)
V(igdat)$name[eee==12]

netdat = network::as.network(unique(dat[,1:2]))
set.seed(1)
plot_df = ggnetwork(netdat, layout = "fruchtermanreingold", cell.jitter = 1, component_wise = T)

mk_cluster = lapply(unique(mk_ct_comp[mk_ct_comp$logFC>0.5,"feature"]), function(x) as.character(mk_ct_comp$group[which(mk_ct_comp$logFC==max(mk_ct_comp$logFC[mk_ct_comp$feature==x]))]))
names(mk_cluster) = unique(mk_ct_comp[mk_ct_comp$logFC>0.5,"feature"])
mk_cluster = unique(reshape2::melt(mk_cluster)[,2:1])
colnames(mk_cluster) = c("gene", "cluster")
plot_df = merge(plot_df, mk_cluster, by.x = "vertex.names", by.y = "gene", all.x = T)
plot_df$cluster[is.na(plot_df$cluster)] = "other genes"
plot_df$highlevel = unlist(lapply(strsplit(plot_df$cluster, "_"), function(x) x[1]))
plot_df$highlevel[is.na(plot_df$highlevel)] = "other genes"

labdf = head(dat[order(abs(dat$estimate), decreasing = T),], 1000)
ggplot(plot_df, aes(x = x, y = y, fill = highlevel)) +
  geom_edges(arrow = arrow(length = unit(6, "pt")), colour = "grey35",
             mapping = aes(xend = xend, yend = yend)) +
  geom_nodes() +
  geom_label(data = unique(plot_df[plot_df$vertex.names %in% c(labdf$tf, labdf$target),]), 
             mapping = aes(label = vertex.names), size=2.7, label.padding = unit(0.1, "lines"))+
  theme_blank()+
  theme(legend.position = "bottom")

ggplot(plot_df, aes(x = x, y = y, fill = cluster)) +
  geom_edges(arrow = arrow(length = unit(6, "pt")), colour = "grey35",
             mapping = aes(xend = xend, yend = yend)) +
  guides(fill = guide_legend(nrow = 2))+
  geom_nodes() +
  geom_label(data = unique(plot_df[plot_df$vertex.names %in% c(labdf$tf, labdf$target),]), 
             mapping = aes(label = vertex.names), size=2.7, label.padding = unit(0.1, "lines"))+
  theme_blank()+
  theme(legend.position = "bottom")
