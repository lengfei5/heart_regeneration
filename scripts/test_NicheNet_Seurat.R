########################################################
########################################################
# Section : not Used codes, test default NicheNet with Seurat object 
# 
########################################################
########################################################
# Step 3: Define a set of potential ligands
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% 
  filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()


# Step 4: Perform NicheNet’s ligand activity analysis on the gene set of interest
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities

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