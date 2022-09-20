##########################################
# TFs and ligand receptor annotations
##########################################
annot = readRDS(paste0('/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))

## import tf annotation
tfs = readRDS(file = paste0('/groups/tanaka/People/current/jiwang/Databases/motifs_TFs/TFs_annot/',
                            'curated_human_TFs_Lambert.rds'))
tfs = unique(tfs$`HGNC symbol`)

## ligand-receptor from nichenet 
## https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_activity_geneset.md
dataPath_nichenet = '../data/NicheNet/'

ligand_target_matrix = readRDS(paste0(dataPath_nichenet,  "ligand_target_matrix.rds"))
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

lr_network = readRDS(paste0(dataPath_nichenet, "lr_network.rds"))
head(lr_network)

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

annotated_tfs = intersect(tfs, annot$gene.symbol.toUse)
annotated_ligands = intersect(ligands, annot$gene.symbol.toUse)
annotated_receptors = intersect(receptors, annot$gene.symbol.toUse)

##########################################
# Elad's snRNA-seq
##########################################
Check.TFs.LR.Elad.multiome = FALSE
if(Check.TFs.LR.Elad.multiome){
  ggs = get_geneName(rownames(aa))
  
  ## filter lowly expressed and highly expressed genes
  Filtering.lowly.highly.expressed.genes = FALSE
  if(Filtering.lowly.highly.expressed.genes){
    require(SingleCellExperiment)
    library(scran)
    library(scater)
    library(scuttle)
    library(Seurat)
    library(SeuratObject)
    
    sce <- as.SingleCellExperiment(aa)
    # rm(aa)
    
    fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
    plotHighestExprs(sce, n=50, exprs_values = "counts")
    
    ave.counts <- calculateAverage(sce, assay.type = "counts")
    
    hist(log10(ave.counts), breaks=100, main="", col="grey80",
         xlab=expression(Log[10]~"average count"))
    
    num.cells <- nexprs(sce, byrow=TRUE)
    hist(log10(num.cells), breaks = 100)
    
    smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
                  xlab=expression(Log[10]~"average count"))
    
    # detected in >= 5 cells, ave.counts >=5 but not too high
    genes.to.keep <- num.cells > 100 & ave.counts >= 10^-2  & ave.counts <10^1  
    summary(genes.to.keep)
    
    # remove mt and ribo genes
    #genes.to.keep = genes.to.keep & ! rownames(sce) %in% gg.mt & ! rownames(sce) %in% gg.rb
    summary(genes.to.keep)
    
    sce <- sce[genes.to.keep, ]
    ave.counts <- calculateAverage(sce, assay.type = "counts")
    sce = sce[order(ave.counts), ]
    
    ggs = get_geneName(rownames(sce))
    
  }
  
  
  expressed_tfs = intersect(tfs, ggs)
  expressed_ligands = intersect(ligands, ggs)
  expressed_receptors = intersect(receptors, ggs)
  
  cat('TFs -- ', length(tfs), '(database)', '-- ', 
      length(annotated_tfs), ' (annotated in axolotl) --', 
      length(expressed_tfs), ' annotated in axolotl\n' )
  
  cat('ligands -- ', length(ligands), '(database)', '-- ', 
      length(annotated_ligands), ' (annotated in axolotl) --', 
      length(expressed_ligands), ' annotated in axolotl\n' )
  
  cat('receptor -- ', length(receptors), '(database)', '-- ', 
      length(annotated_receptors), ' (annotated in axolotl) --', 
      length(expressed_receptors), ' annotated in axolotl\n' )
  
  #lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
  #head(lr_network_expressed)
  
  res = data.frame(tfs = c(length(tfs), length(annotated_tfs), length(expressed_tfs)),
                   ligands = c(length(ligands), length(annotated_ligands), length(expressed_ligands)), 
                   receptors = c(length(receptors), length(annotated_receptors), length(expressed_receptors)))
  rownames(res) = c('database', 'annotated.axolotl', 'detected.snRNA')
  
  write.csv(res, file = paste0(resDir, '/tfs_ligands_receptors_coverage_snRNA.csv'), 
            row.names = TRUE)
  
  ## double check individual profile of TFs, ligands and receptors
  ggs = get_geneName(rownames(sce))
  mm = match(ggs, expressed_tfs)
  kk = which(!is.na(mm) & ggs != '')
  
  
  subs = sce[kk, ]
  ave.counts <- calculateAverage(subs, assay.type = "counts")
  num.cells <- nexprs(subs, byrow=TRUE)
  
  res = data.frame(gene = row.names(subs), ave.counts = ave.counts, number.cells = num.cells, 
                   stringsAsFactors = FALSE)
  res = res[order(-res$ave.counts), ]
  
  write.csv(res, file = paste0(resDir, '/tfs_detected_aveCounts_numCells.csv'), 
            row.names = TRUE)
  
  
  Idents(aa) = aa$subtypes
  
  for(n in 1:length(kk))
    #for(n in 1:10)
  {
    # n = 1
    cat(n, '\n')
    p1 = FeaturePlot(aa, features = rownames(sce)[kk[n]])
    p2 = VlnPlot(aa, features = rownames(sce)[kk[n]]) + NoLegend()
    
    p1 + p2
    
    ggsave(filename = paste0(resDir, '/tfs/CoverageTest_', rownames(sce)[kk[n]], '.pdf'),
           width = 10, height = 16)
    
  }
  
  
  ggs = get_geneName(rownames(sce))
  mm = match(ggs, expressed_ligands)
  kk = which(!is.na(mm) & ggs != '')
  
  subs = sce[kk, ]
  ave.counts <- calculateAverage(subs, assay.type = "counts")
  num.cells <- nexprs(subs, byrow=TRUE)
  
  res = data.frame(gene = row.names(subs), ave.counts = ave.counts, number.cells = num.cells, 
                   stringsAsFactors = FALSE)
  res = res[order(-res$ave.counts), ]
  
  write.csv(res, file = paste0(resDir, '/ligands_detected_aveCounts_numCells.csv'), 
            row.names = TRUE)
  
  ggs = get_geneName(rownames(sce))
  kk = which(!is.na(match(ggs, expressed_ligands)))
  Idents(aa) = aa$subtypes
  
  for(n in 1:length(kk))
    #for(n in 1:10)
  {
    # n = 1
    cat(n, '\n')
    p1 = FeaturePlot(aa, features = rownames(sce)[kk[n]])
    p2 = VlnPlot(aa, features = rownames(sce)[kk[n]]) + NoLegend()
    
    p1 + p2
    
    ggsave(filename = paste0(resDir, '/ligands/CoverageTest_', rownames(sce)[kk[n]], '.pdf'),
           width = 10, height = 16)
    
  }
  
  ggs = get_geneName(rownames(sce))
  mm = match(ggs, expressed_receptors)
  kk = which(!is.na(mm) & ggs != '')
  
  subs = sce[kk, ]
  ave.counts <- calculateAverage(subs, assay.type = "counts")
  num.cells <- nexprs(subs, byrow=TRUE)
  
  res = data.frame(gene = row.names(subs), ave.counts = ave.counts, number.cells = num.cells, 
                   stringsAsFactors = FALSE)
  res = res[order(-res$ave.counts), ]
  
  write.csv(res, file = paste0(resDir, '/receptors_detected_aveCounts_numCells.csv'), 
            row.names = TRUE)
  
  
  kk = which(!is.na(match(ggs, expressed_receptors)))
  Idents(aa) = aa$subtypes
  
  for(n in 1:length(kk))
    #for(n in 1:10)
  {
    # n = 1
    cat(n, '\n')
    p1 = FeaturePlot(aa, features = rownames(sce)[kk[n]])
    p2 = VlnPlot(aa, features = rownames(sce)[kk[n]]) + NoLegend()
    
    p1 + p2
    
    ggsave(filename = paste0(resDir, '/receptors/CoverageTest_', rownames(sce)[kk[n]], '.pdf'),
           width = 10, height = 16)
    
  }
  
}

##########################################
# Lust et al. brain data to compare 
# start with multiome and do the same analysis as before 
# and compare the snRNA-seq and multiome with the rank of TFs, ligands and receptors
##########################################
lust = readRDS(file = paste0("../data/Lust_brain/all_nuclei_clustered_highlevel_anno.rds"))

DimPlot(lust, dims = c(1,2), reduction = 'umap_harmony', group.by = 'high_level_anno')

srat_mtm = subset(lust, chem == 'multiome')
srat_sn = subset(lust, chem == 'v3.1')

rm(lust)


## filter lowly expressed and highly expressed genes
Filtering.lowly.highly.expressed.genes = FALSE
if(Filtering.lowly.highly.expressed.genes){
  require(SingleCellExperiment)
  library(scran)
  library(scater)
  library(scuttle)
  library(Seurat)
  library(SeuratObject)
  
  #sce <- as.SingleCellExperiment(srat_mtm)
  sce <- as.SingleCellExperiment(srat_sn)
  ggs = rownames(sce)
  # rm(aa)
  
  fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
  plotHighestExprs(sce, n=50, exprs_values = "counts")
  
  ave.counts <- calculateAverage(sce, assay.type = "counts")
  
  hist(log10(ave.counts), breaks=100, main="", col="grey80",
       xlab=expression(Log[10]~"average count"))
  
  num.cells <- nexprs(sce, byrow=TRUE)
  hist(log10(num.cells), breaks = 100)
  
  smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
                xlab=expression(Log[10]~"average count"))
  
  
  # detected in >= 5 cells, ave.counts >=5 but not too high
  genes.to.keep <- num.cells > 100 & ave.counts >= 10^-3  & ave.counts <10^2  
  summary(genes.to.keep)
  
  # remove mt and ribo genes
  #genes.to.keep = genes.to.keep & ! rownames(sce) %in% gg.mt & ! rownames(sce) %in% gg.rb
  summary(genes.to.keep)
  
  sce <- sce[genes.to.keep, ]
  ave.counts <- calculateAverage(sce, assay.type = "counts")
  
  sce = sce[order(ave.counts), ]
  
  ggs = rownames(sce)
  
}

expressed_tfs = intersect(tfs, ggs)
expressed_ligands = intersect(ligands, ggs)
expressed_receptors = intersect(receptors, ggs)

cat('TFs -- ', length(tfs), '(database)', '-- ', 
    length(annotated_tfs), ' (annotated in axolotl) --', 
    length(expressed_tfs), ' annotated in axolotl\n' )

cat('ligands -- ', length(ligands), '(database)', '-- ', 
    length(annotated_ligands), ' (annotated in axolotl) --', 
    length(expressed_ligands), ' annotated in axolotl\n' )

cat('receptor -- ', length(receptors), '(database)', '-- ', 
    length(annotated_receptors), ' (annotated in axolotl) --', 
    length(expressed_receptors), ' annotated in axolotl\n' )

#lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
#head(lr_network_expressed)

res = data.frame(tfs = c(length(tfs), length(annotated_tfs), length(expressed_tfs)),
                 ligands = c(length(ligands), length(annotated_ligands), length(expressed_ligands)), 
                 receptors = c(length(receptors), length(annotated_receptors), length(expressed_receptors)))
rownames(res) = c('database', 'annotated.axolotl', 'detected.multiome')

res = t(res)

res = cbind(res, c(length(expressed_tfs), length(expressed_ligands), length(expressed_receptors)))
colnames(res)[ncol(res)] = 'detected.sn'

write.csv(res, file = paste0(resDir, '/tfs_ligands_receptors_coverage_multiome.vs.snRNA_LustBrain.csv'), 
          row.names = TRUE)


##########################################
# compare explicitly TFs
##########################################
ggs = rownames(sce)
mm = match(ggs, unique(c(expressed_tfs, expressed_ligands, expressed_receptors)))
kk = which(!is.na(mm) & ggs != '')

subs = sce[kk, ]
ave.counts <- calculateAverage(subs, assay.type = "counts")
num.cells <- nexprs(subs, byrow=TRUE)

subs_sn = as.SingleCellExperiment(srat_sn)
subs_sn = subs_sn[match(rownames(subs), rownames(subs_sn)), ]
ave.counts_sn <- calculateAverage(subs_sn, assay.type = "counts")
num.cells_sn <- nexprs(subs_sn, byrow=TRUE)


smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
              xlab=expression(Log[10]~"average count"))

smoothScatter(log10(ave.counts_sn), num.cells_sn, ylab="Number of cells",
              xlab=expression(Log[10]~"average count"))

keep = data.frame(ave.counts = c(ave.counts, ave.counts_sn), 
                  num.cells = c(num.cells/ncol(subs), num.cells_sn/ncol(subs_sn)), 
                  condition = c(rep('multiome', length(ave.counts)), rep('sn', length(ave.counts_sn))))
keep$ave.counts = log10(keep$ave.counts)

xx = data.frame(multiome = log10(ave.counts), sn = log10(ave.counts_sn))
xx = data.frame(xx, gene = rownames(xx))
xx$group = NA
xx$group[!is.na(match(xx$gene, expressed_tfs))] = 'tfs'
xx$group[!is.na(match(xx$gene, expressed_ligands))] = 'ligands'
xx$group[!is.na(match(xx$gene, expressed_receptors))] = 'receptors'

p1 = ggplot(data = xx, aes(x = multiome, y = sn, color = group)) +
  geom_point(size = 0.7) + 
  geom_abline(slope = 1,  intercept = 0, col = 'red', lwd = 1.) +
  geom_abline(slope = 1, intercept = 0.38, lwd = 1.) +
  theme_bw() + 
  ggtitle('log10 (average umi )') + NoLegend()

xx = data.frame(multiome =num.cells/ncol(subs) , sn = num.cells_sn/ncol(subs_sn))
xx = data.frame(xx, gene = rownames(xx))
xx$group = NA
xx$group[!is.na(match(xx$gene, expressed_tfs))] = 'tfs'
xx$group[!is.na(match(xx$gene, expressed_ligands))] = 'ligands'
xx$group[!is.na(match(xx$gene, expressed_receptors))] = 'receptors'

p2 = ggplot(data = xx, aes(x = multiome, y = sn, color = group)) +
  geom_point(size = 0.7) + 
  geom_abline(slope = 1,  intercept = 0, col = 'red') +
  theme_bw() + 
  ggtitle('detected in % of cell ')

p1 | p2

ggsave(filename = paste0(resDir, '/multiome_vs._snRNA_LustBrainData.pdf'), width = 20, height = 8)


res = data.frame(gene = row.names(subs), ave.counts_multiome = ave.counts, number.cells_multiome = num.cells, 
                 pct.cells_multiome = num.cells/ncol(subs),
                 ave.counts_sn = ave.counts_sn, number.cells_sn = num.cells_sn, 
                 pct.cells_sn = num.cells_sn/ncol(subs_sn),
                 stringsAsFactors = FALSE)

res = res[order(-res$ave.counts_multiome), ]

write.csv(res, file = paste0(resDir, '/tfs_ligands_receptors_detected_aveCounts_numCells_LustBrain.csv'), 
          row.names = TRUE)


for(n in 1:nrow(res))
#for(n in 1:10)
{
  # n = 1
  cat(n, '--', res$gene[n], '\n')
  p1 = FeaturePlot(srat_mtm, reduction = 'umap_harmony', features = res$gene[n]) +
    ggtitle(paste0(res$gene[n],  '-- multiome'))
  p2 = FeaturePlot(srat_sn, reduction = 'umap_harmony', features = res$gene[n]) +
    ggtitle(paste0(res$gene[n],  '-- sn'))
  
  p1 | p2
  
  ggsave(filename = paste0(resDir, '/coverageTest_LustBrain/CoverageTest_rank_', n, '_', res$gene[n], '.pdf'),
         width = 20, height = 8)
  
}

