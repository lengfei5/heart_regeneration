

# srat_cr = subset(srat_cr,  downsample = 500)
jstart <- Sys.time() 
suppressPackageStartupMessages(library("argparse"))
library(topicmodels)
library(dplyr)
library(ggplot2)
library(Matrix)
library(here)
library(parallel)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument('inpath', metavar='INFILE',
                    help='.RData where count mat is in count.dat$counts or .rds object to count mat')
parser$add_argument('outdir', metavar='OUTDIR',
                    help='Out directory')
parser$add_argument("-t", "--topics", metavar='Comma sep string', required=TRUE,
                    help='CSV of topics to iterate')
parser$add_argument("-b", "--binarizemat", action="store_true", default=FALSE,
                    help="Binarize matrix")
parser$add_argument("-n", "--projname", metavar='Name of project', default="MyProj",
                    help="Name of project for naming pdf and Robj output. 
                    Make this meaningful otherwise it will overwrite projects!")
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")
parser$add_argument("--SkipPlots", action="store_true", default=FALSE,
                    help="Do not make plots, default FALSE")
parser$add_argument("--RemoveDupRows", action="store_true", default=FALSE,
                    help="Remove duplicated rows, default FALSE")
parser$add_argument("--RemoveEmptyCells", action="store_true", default=FALSE,
                    help="Remove empty cols, default FALSE")
parser$add_argument("--SkipMeanVar", action="store_true", default=FALSE,
                    help="Do not make plots, default FALSE")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

setwd(here())
print(paste("Work directory: ", getwd()))

# Run LDA on count matrix -------------------------------------------------
cat("Running LDA \n")
count.mat = srat_cr@assays$ATAC@counts

# remove empty cols
print("Removing empty cells...")
print(dim(count.mat))
cols.empty <- colSums(count.mat) == 0
count.mat <- count.mat[, !cols.empty]
print(dim(count.mat))


# binarize matrix
count.mat.orig <- count.mat
count.mat <- BinarizeMatrix(count.mat)
print(paste('Max count after binarizing', max(count.mat)))

topic.vec = c(5, 10)
tic()
if (length(topic.vec) > 1){
  print("Running multicore LDA for topics:")
  print(topic.vec)
  
  out.lda <- parallel::mclapply(topic.vec, function(nc)
    {topicmodels::LDA(x = t(count.mat), k = nc, method = "Gibbs", control=list(seed=0))}, 
    mc.cores = length(topic.vec)
  )
  
} else {
  print("Running single LDA for topics:")
  print(topic.vec)
  
  out.lda <- topicmodels::LDA(x = t(count.mat), k = topic.vec, method = "Gibbs", control=list(seed=0))
  
}

toc()

tm.result <- posterior(out.lda)
tm.result <- AddTopicToTmResult(tm.result)
topics.mat = tm.result$topics
# save output
#print("Saving LDA")
#save(out.lda, count.mat, count.mat.orig, file = outpath)
#print("Time elapsed after LDA")
#print(Sys.time() - jstart)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)

dat.umap.long <- data.frame(cell = rownames(umap.out$layout), 
                            umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 
                                                                                  2], stringsAsFactors = FALSE)
dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings = jsettings)

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", 
               "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point() +
  theme_bw() +
  ggtitle(paste("LDA test")) +
  scale_color_manual(values = cbPalette) +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "right") 

ggplot(dat.umap2, aes(x = umap1, y = umap2, color = louvain)) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = cbPalette) +
  ggtitle("From TSS-TES") +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "bottom")

print("Time elapsed after tuning")
print(Sys.time() - jstart)



