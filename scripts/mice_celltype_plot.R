##########################################
# neonatal mice
##########################################
aa = readRDS(file = paste0('ref_scRNAseq_neonatalMice_clean.v1.2.rds')) # specify the path

DefaultAssay(aa) = 'RNA'

kk = which(aa$celltype == 'Immune')
aa$celltype[kk] = aa$subtype[kk]
aa$celltype[which(aa$celltype == 'B cells')] = 'B'
aa$celltype[which(aa$celltype == 'T cells')] = 'T'
aa$celltype[grep('Macrophage', aa$celltype)] = "Macrophage"

## specify colors
require(RColorBrewer)
cols = brewer.pal(9,"Set1")
cols = c(cols[c(1:3, 5, 6, 9)], brewer.pal(8,"Set2")[-c(2, 4)])
aa$celltype = factor(aa$celltype, levels = c('CM', "EC", 'FB',"B", 'Macrophage', 'T', 
                                             'EPI', 'DC-like', "Monocyte", "Pericyte", 'SMC', "Gra" 
))

cols = c("#BBDEFB",
         "#4CC9F0",
         "#4AAFF0",
         "#49A2F0",
         "#4895EF",
         "#4361EE",
         "#414CDC",
         "#4042D3",
         "#3F37C9",
         "#941F56",
         "#AD1833",
         "#C61010")
cols = cols[length(cols):1]


p1 = DimPlot(aa, group.by = 'celltype', label = TRUE, repel = TRUE, cols = cols)

p2 = FeaturePlot(aa, features = 'Axl') +  
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))

p1 + p2


##########################################
# adult mice
##########################################
refs = readRDS(file = paste0('/ref_scRNAseq_adultMice_clean.v1.rds'))

## color to specify
require(RColorBrewer)
cols = brewer.pal(9,"Set1")
refs$celltype = factor(refs$celltype, levels = c('CM', "EC", 'FB', 'GN', "B", 'MHCII.Mphage', 
                                                 'Mphage.MCT', 'prolife.Mphage', 'NK.T'))

p1 = DimPlot(refs, reduction = 'umap', group.by = 'celltype', raster = T,shuffle= T, pt.size = 2, 
             label = TRUE, repel = TRUE, 
             cols = cols)


#DimPlot(refs, reduction = 'umap', group.by = 'subtype',raster = T,shuffle= T, pt.size = 2, 
#       label = TRUE, repel = TRUE)
p2 = FeaturePlot(refs, features = 'Axl') +  
  scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))

p1 + p2

