##########################################
# test the TFs and ligand receptor coverages 
##########################################
annot = readRDS(paste0('/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr_curated.geneSymbol.toUse.rds'))

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

##########################################
# Lust et al. brain data to compare  
##########################################
lust = readRDS(file = paste0("../data/Lust_brain/all_nuclei_clustered_highlevel_anno.rds"))

DimPlot(lust, dims = c(1,2), reduction = 'umap_harmony', group.by = 'high_level_anno')

srat_mtm = subset(lust, chem == 'multiome')
srat_sn = subset(lust, chem == 'v3.1')

rm(lust)







