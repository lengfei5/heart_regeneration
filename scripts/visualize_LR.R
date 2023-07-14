##########################################################################
##########################################################################
# Project: heart regeneration 
# Script purpose: visualize ligand-receptor analysis
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Jul 11 17:09:43 2023
##########################################################################
##########################################################################

########################################################
########################################################
# Section : test Circos plot from NicheNet 
# https://github.com/saeyslab/nichenetr/blob/master/vignettes/circos.md
########################################################
########################################################

library(nichenetr)
library(tidyverse)
library(circlize)

hnscc_expression = readRDS(url("https://zenodo.org/record/3260758/files/hnscc_expression.rds"))
expression = hnscc_expression$expression
sample_info = hnscc_expression$sample_info # contains meta-information about the cells

tumors_remove = c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")

CAF_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `non-cancer cell type` == "CAF") %>% pull(cell)
endothelial_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `non-cancer cell type` == "Endothelial") %>% pull(cell)
malignant_ids = sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) & `classified  as cancer cell` == 1) %>% pull(cell)

expressed_genes_CAFs = expression[CAF_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_endothelial = expression[endothelial_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()
expressed_genes_malignant = expression[malignant_ids,] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% names()

ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns


pemt_geneset = readr::read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"), col_names = "gene") %>% pull(gene) %>% .[. %in% rownames(ligand_target_matrix)] # only consider genes also present in the NicheNet model - this excludes genes from the gene list for which the official HGNC symbol was not used by Puram et al.
head(pemt_geneset)


background_expressed_genes = expressed_genes_malignant %>% .[. %in% rownames(ligand_target_matrix)]
head(background_expressed_genes)

lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands_CAFs = intersect(ligands,expressed_genes_CAFs)
expressed_ligands_endothelial = intersect(ligands,expressed_genes_endothelial)
expressed_ligands = union(expressed_ligands_CAFs, expressed_genes_endothelial)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_malignant)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
head(potential_ligands)


ligand_activities = predict_ligand_activities(geneset = pemt_geneset, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)


ligand_activities %>% arrange(-aupr_corrected) 

best_upstream_ligands = ligand_activities %>% top_n(20, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
head(best_upstream_ligands)

best_upstream_ligands %>% intersect(expressed_ligands_CAFs) 

best_upstream_ligands %>% intersect(expressed_ligands_endothelial)


# lot of overlap between both cell types in terms of expressed ligands
# therefore, determine which ligands are more strongly expressed in which of the two
ligand_expression_tbl = tibble(
  ligand = best_upstream_ligands, 
  CAF = expression[CAF_ids,best_upstream_ligands] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}),
  endothelial = expression[endothelial_ids,best_upstream_ligands] %>% apply(2,function(x){10*(2**x - 1)}) %>% apply(2,function(x){log2(mean(x) + 1)}))

CAF_specific_ligands = ligand_expression_tbl %>% filter(CAF > endothelial + 2) %>% pull(ligand)
endothelial_specific_ligands = ligand_expression_tbl %>% filter(endothelial > CAF + 2) %>% pull(ligand)
general_ligands = setdiff(best_upstream_ligands,c(CAF_specific_ligands,endothelial_specific_ligands))


ligand_type_indication_df = tibble(
  ligand_type = c(rep("CAF-specific", times = CAF_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length()),
                  rep("Endothelial-specific", times = endothelial_specific_ligands %>% length())),
  ligand = c(CAF_specific_ligands, general_ligands, endothelial_specific_ligands))


active_ligand_target_links_df = best_upstream_ligands %>% 
  lapply(get_weighted_ligand_target_links,geneset = pemt_geneset, 
         ligand_target_matrix = ligand_target_matrix, n = 250) %>% 
  bind_rows()

active_ligand_target_links_df = active_ligand_target_links_df %>% 
  mutate(target_type = "p_emt") %>% 
  inner_join(ligand_type_indication_df) 
# if you want ot make circos plots for multiple gene sets, 
# combine the different data frames and differentiate which target belongs to which gene set via the target type

cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.66)

active_ligand_target_links_df_circos = active_ligand_target_links_df %>% 
  filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), 
                            active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), 
                            active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% 
  filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

## Prepare the circos visualization: give each segment of ligands and targets a specific color and order
grid_col_ligand =c("General" = "lawngreen",
                   "CAF-specific" = "royalblue",
                   "Endothelial-specific" = "gold")
grid_col_target =c(
  "p_emt" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

# extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% mutate(ligand = paste(ligand," "))

circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% 
  mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% 
  mutate(transparency = 1-weight) %>% .$transparency 

# Prepare the circos visualization: order ligands and targets
target_order = circos_links$target %>% unique()
ligand_order = c(CAF_specific_ligands,general_ligands,endothelial_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)

width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, 
      times = (circos_links %>% filter(ligand_type == "CAF-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, 
      times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, 
      times = (circos_links %>% filter(ligand_type == "Endothelial-specific") %>% 
                 distinct(ligand) %>% nrow() -1)), 
  width_ligand_target,
  rep(width_same_cell_same_target_type, 
      times = (circos_links %>% filter(target_type == "p_emt") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
)

circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) 

circos.clear()


circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #