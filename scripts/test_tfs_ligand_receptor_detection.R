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
  
  ggs = get_geneName(rownames(sce))
  
}


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
