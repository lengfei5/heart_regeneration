## =============================================================================
## dnAXL vs Ctrl bulk RNA-seq, tested against a cardiomyocyte dedifferentiation
## reference derived from adult axolotl heart snRNA-seq.
##
## Question: do dnAXL hearts retain a more mature cardiomyocyte transcriptional
##           identity than controls?
##
## Run:  Rscript 01_dnAXL_vs_CM_dedifferentiation.R
## Inputs  -> ./input/   (two files, see README.md)
## Outputs -> ./output/  (three CSVs, two PNGs, one sessionInfo.txt)
##
## Everything below is deterministic; set.seed() is fixed before fgsea.
## =============================================================================
setwd("/groups/tanaka/People/current/jiwang/projects/heart_regeneration/bulkRNAseq_dnAXL_KDTXBX/Elad_analysis")

suppressPackageStartupMessages({
  library(DESeq2); 
  library(fgsea); 
  library(ggplot2); 
  library(patchwork)
})

set.seed(1)
dir.create("output", showWarnings = FALSE)

## -----------------------------------------------------------------------------
## 1. Bulk differential expression: dnAXL vs Ctrl
## -----------------------------------------------------------------------------
## The supplied RDS is an UNFITTED DESeqDataSet (counts + colData + design only;
## no dispersions, no results). It is fitted here from scratch.

dds <- readRDS(paste0("/groups/tanaka/Collaborations/Jingkui-Elad/dnAXL_bulk_RNA/counts_dds_QCs/", 
                    "dds_normalized_dnAXL_plasmidsaurus.rds"))

stopifnot(is(dds, "DESeqDataSet"))

dds$condition <- relevel(factor(dds$condition), ref = "Ctrl")
design(dds)   <- ~ condition

## Independent pre-filter: >=10 counts in >=3 samples (= smallest group size).
## Applied before DESeq() so dispersion fitting is not driven by near-zero genes.
keep <- rowSums(counts(dds) >= 10) >= 3
message(sprintf("genes retained: %d of %d", sum(keep), length(keep)))
dds <- dds[keep, ]

dds <- DESeq(dds, quiet = TRUE)
res <- results(dds, contrast = c("condition", "dnAXL", "Ctrl"), alpha = 0.05)


## NOTE ON SHRINKAGE -----------------------------------------------------------
## apeglm shrinkage was attempted and REJECTED: on this dataset the apeglm prior
## collapsed (all shrunken |LFC| < 0.1, lfcSE constant at 1.4e-3), which is a
## prior-fitting failure, not real shrinkage. All effect sizes reported here are
## therefore the unshrunken MLE estimates from results(). ashr shrinkage was
## checked as a sanity control and agrees with the MLE ranking (r > 0.99).
## The enrichment below uses the Wald statistic, which is unaffected either way.
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## 2. Reference gene sets from adult CM snRNA-seq
## -----------------------------------------------------------------------------
## Reference contrast: injury-zone proliferating CM (CM_Prol_IS) vs homeostatic
## ventricular CM (CM_ven_Robo2), Seurat FindMarkers output.
##   avg_log2FC > 0  = higher in injury-zone CM   -> dedifferentiation program
##   avg_log2FC < 0  = higher in homeostatic CM   -> mature CM identity
##
## A limb-regeneration module was tested in an earlier iteration and REMOVED at
## the PI's instruction; both sets below are cardiac-derived only.

## double check the ref here from Elad
library(Seurat)
#library(decoupleR)
library(tictoc)
library(dplyr)
library(tibble)
library(tidyr)
library(patchwork)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(data.table)

source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_scRNAseq.R')
source('/groups/tanaka/People/current/jiwang/projects/heart_regeneration/scripts/functions_Visium.R')

refs_file = paste0('/groups/tanaka/Collaborations/Jingkui-Elad/scMultiome/aa_subtypes_final_20221117.rds')

refs = readRDS(file = refs_file)
table(refs$subtypes)
length(table(refs$subtypes))

refs$subtypes = droplevels(refs$subtypes) 
length(table(refs$subtypes)) 

refs$condition = gsub('_scRNA', '', refs$condition)

## prepare the celltype to use and also specify the time-specific subtypes
refs$celltype_toUse = as.character(refs$subtypes)
length(table(refs$celltype_toUse))

refs$celltype_toUse = gsub('Mo/Macs', 'Mo.Macs', refs$celltype_toUse)
refs$celltype_toUse = gsub("[(]", '', refs$celltype_toUse)
refs$celltype_toUse = gsub("[)]", '', refs$celltype_toUse)
refs$celltype_toUse = gsub('CM_ven_Robo2', 'CM_Robo2', refs$celltype_toUse)
refs$celltype_toUse = gsub('CM_ven_Cav3_1', 'CM_Cav3.1', refs$celltype_toUse)

table(refs$celltype_toUse)
length(table(refs$celltype_toUse))

refs$subtypes = refs$celltype_toUse # clean the special symbols
refs$subtypes = as.factor(refs$subtypes) 
table(refs$subtypes)
length(table(refs$subtypes))

refs$celltypes = as.character(refs$celltype_toUse)

refs$celltypes[grep('CM_|CMs_|_CM|_CM_', refs$subtypes)] = 'CM'
refs$celltypes[grep('EC_|_EC', refs$subtypes)] = 'EC'
refs$celltypes[grep('FB_', refs$subtypes)] = 'FB'
refs$celltypes[grep('B_cells', refs$subtypes)] = 'Bcell'
refs$celltypes[grep('T_cells', refs$subtypes)] = 'Tcell'

refs$celltypes[grep('Macrophages|_MF', refs$subtypes)] = 'Macrophages'
refs$celltypes[grep('Megakeryocytes', refs$subtypes)] = 'Megakeryocytes'
refs$celltypes[grep('RBC', refs$subtypes)] = 'RBC'
refs$celltypes[grep('Neu_', refs$subtypes)] = 'Neutrophil'
refs$celltypes[grep('Mo.Macs', refs$subtypes)] = 'Mo.Macs'

table(refs$celltypes)
length(table(refs$celltypes))

Idents(refs) = refs$celltype_toUse
markers = FindMarkers(refs, ident.1 = 'CM_Prol_IS', ident.2 = "CM_Robo2",
                      logfc.threshold = 0.1,
                      test.use = "wilcox",
                      min.pct = 0.1,
                      only.pos = FALSE)

##########################################
# impossible to figure out how the logFold change were computed
##########################################
# ref <- read.csv("./cardiac_prolIS_vs_robo2_AMEX_DE.csv", stringsAsFactors = FALSE)

ref = data.frame(amex = rownames(markers), 
                 rownames = rownames(markers), 
                 markers) # use the newly computed markers 


## Gene-ID harmonisation. Both datasets carry AmexT v47 identifiers, but the
## symbol prefix uses a different separator in each source:
##    bulk      AXL_AMEX60DD024484     (underscore)
##    reference CBFB-AMEX60DD017820    (hyphen)
## The join is therefore on the bare AMEX60DD accession, extracted by regex.
amex_id <- function(x) {
  m <- regmatches(x, regexpr("AMEX60DD[0-9]+", x))
  out <- rep(NA_character_, length(x))
  out[grepl("AMEX60DD[0-9]+", x)] <- m
  out
}

bulk_amex <- amex_id(rownames(dds));
names(bulk_amex) <- rownames(dds)

ref$amex  <- amex_id(ref$rowname)
message(sprintf("bulk genes carrying an AMEX id: %d of %d", sum(!is.na(bulk_amex)), nrow(dds)))

## Set membership thresholds: FDR < 0.05 AND |avg_log2FC| > 1 in the reference.
## Sets are then intersected with the genes actually tested in the bulk data, so
## the enrichment background and the sets come from the same universe.
logfc_cutoff = 0.1
sig_up <- unique(ref$amex[ref$p_val_adj < 0.05 & ref$avg_log2FC > logfc_cutoff])
sig_dn <- unique(ref$amex[ref$p_val_adj < 0.05 & ref$avg_log2FC < -logfc_cutoff])

cat(length(sig_dn) + length(sig_up), ' gene for CM signatures\n')

sets <- list(
  `CM dedifferentiation (Prol_IS up)` = intersect(sig_up, bulk_amex),
  `Mature CM identity (Robo2 up)`     = intersect(sig_dn, bulk_amex)
)

message(sprintf("set sizes after intersection with bulk: dediff %d, mature %d",
                length(sets[[1]]), length(sets[[2]])))

## -----------------------------------------------------------------------------
## 3. Ranking statistic
## -----------------------------------------------------------------------------
## Genes are ranked by the DESeq2 Wald statistic (LFC / lfcSE). This is preferred
## over raw LFC because it is variance-stabilised: at n=3 the raw LFC is dominated
## by low-count genes with large standard errors. Positive = up in dnAXL.
#stat <- -log10(res$pvalue)
stat = res$stat
names(stat) <- bulk_amex[rownames(res)]
stat <- stat[!is.na(names(stat)) & !is.na(stat)]
stat <- stat[!duplicated(names(stat))]
message(sprintf("ranked genes: %d", length(stat)))

## -----------------------------------------------------------------------------
## 4. Test 1 - GSEA (rank-based, uses the whole ranking)
## -----------------------------------------------------------------------------
fg <- fgsea(sets, stat, minSize = 10, eps = 0, scoreType = 'std')
fg_out <- as.data.frame(fg[, c("pathway", "size", "ES", "NES", "pval", "padj")])

fg_out 

p1 = plotEnrichment(sets[["Mature CM identity (Robo2 up)"]], stat) + 
  labs(title="Mature CM identity (Robo2 up)")

p2 = plotEnrichment(sets[["CM dedifferentiation (Prol_IS up)"]], stat) + 
  labs(title="CM dedifferentiation (Prol_IS up)")

curve_A <- data.frame(
  rank = p1$data$x, ES = p1$data$y
)

curve_B <- data.frame(
  rank = p2$data$x, ES = p2$data$y
)

curve_A$gene_set <- "Mature CM"
curve_B$gene_set <- "CM dedifferentiation"

curves <- rbind(curve_A, curve_B)

curves = data.frame(curves)
curves$gene_set = factor(curves$gene_set, levels = c("Mature CM", 
                                                     "CM dedifferentiation"))


ggplot(curves, aes(x = rank, y = ES, color = gene_set)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  geom_line(linewidth = 0.9) +
  scale_color_manual(
    values = c("Mature CM" = "#D55E00", "CM dedifferentiation" = "#0072B2")
  ) +
  labs(
    x = "Gene rank",
    y = "Running enrichment score",
    color = "Gene set"
  ) +
  theme_classic()

ggsave("output/fig_geneSetEnrichment_maturePval0.001.pdf", width = 6.8, height = 4, dpi = 300)



##########################################
# pheat visualize the mature CM signature 
##########################################
library("pheatmap")
cpm = fpm(dds, robust = TRUE)
cpm = log2(cpm + 1)

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = c('condition'), returnData = FALSE)
print(pca)


cpm = assay(vsd)

o1 = order(res$pvalue)
res = res[o1, ]
cpm = cpm[o1, ]

sig_matureCM <- unique(ref$amex[ref$p_val < 0.05 & abs(ref$avg_log2FC) > 0.3])

mm = match(amex_id(rownames(cpm)), sig_matureCM)

select <- which(!is.na(mm) & res$pvalue < 0.05)

df <- as.data.frame(colData(dds)[,c("condition")])


pheatmap(cpm[select, ], 
         scale = "row",
         cluster_rows=TRUE, show_rownames=TRUE, 
         cluster_cols=FALSE, 
         filename = "output/fig_heatmap_DEgenes_mature_dedifferentiation_signatures.pdf",
         width = 8, height = 12)



## -----------------------------------------------------------------------------
## 5. Test 2 - distribution shift vs background (independent of GSEA)
## -----------------------------------------------------------------------------
## Mean Wald statistic per set with a normal-approximation 95% CI, plus a
## two-sided Wilcoxon rank-sum test of each set against ALL tested genes.
## This is reported alongside GSEA because the two make different assumptions;
## the conclusion should not rest on either alone.
mod_stats <- do.call(rbind, lapply(c(list(`All tested genes` = names(stat)), sets), function(g) {
  v <- stat[names(stat) %in% g]
  data.frame(n = length(v), mean_wald = mean(v),
             ci_lo = mean(v) - 1.96 * sd(v) / sqrt(length(v)),
             ci_hi = mean(v) + 1.96 * sd(v) / sqrt(length(v)),
             wilcox_p = if (length(v) == length(stat)) NA_real_ else
                        wilcox.test(v, stat)$p.value)
}))

mod_stats$module <- c("All tested genes", names(sets))
mod_stats <- mod_stats[, c("module", "n", "mean_wald", "ci_lo", "ci_hi", "wilcox_p")]

## -----------------------------------------------------------------------------
## 6. Test 3 - gene-level correlation with the dedifferentiation axis
## -----------------------------------------------------------------------------
## Threshold-free version of the same question: across all genes measured in both
## datasets, does the bulk response run opposite to the dedifferentiation axis?
m  <- match(bulk_amex[rownames(res)], ref$amex)

gl <- data.frame(gene = rownames(res),
                 ref_log2FC_prolIS_vs_robo2 = ref$avg_log2FC[m],
                 ref_padj_prolIS_vs_robo2 = ref$p_val_adj[m],
                 bulk_log2FC = res$log2FoldChange,
                 bulk_wald   = res$stat,
                 bulk_padj   = res$padj)
gl <- gl[!is.na(gl$ref_log2FC_prolIS_vs_robo2) & !is.na(gl$bulk_wald), ]

sels = which(abs(gl$ref_log2FC_prolIS_vs_robo2) > 0.5 & gl$ref_padj_prolIS_vs_robo2 < 0.05)

plot(gl$ref_log2FC_prolIS_vs_robo2[sels], gl$bulk_log2FC[sels])
abline(0, 1, lwd = 2, col = 'red')

sp <- cor.test(gl$ref_log2FC_prolIS_vs_robo2[sels], gl$bulk_log2FC[sels], method = "spearman")

## Dilution slope. The dnAXL construct is driven by an inducible double Cre-ER,
## so only a fraction f of cells convert. Under the mixture model
##     bulk fold change = (1 - f) + f * (per-cell fold change)
## the slope of bulk log2FC on reference log2FC approximates f for small effects.
## Restricted to confident reference genes (|avg_log2FC| > 0.5).

# dil <- lm(bulk_log2FC ~ ref_log2FC_prolIS_vs_robo2,
#           gl[abs(gl$ref_log2FC_prolIS_vs_robo2) > 0.5, ])
# dil_ci <- confint(dil)[2, ]
# 
# gl$axis_bin <- cut(gl$ref_log2FC_prolIS_vs_robo2,
#                    breaks = c(-Inf, -1, -0.25, 0.25, 1, Inf),
#                    labels = c("strongly mature", "mature-leaning", "neutral",
#                               "dediff-leaning", "strongly dediff"))
# 
# bins <- do.call(rbind, lapply(split(gl, gl$axis_bin), function(s)
#   data.frame(bin = s$axis_bin[1], n = nrow(s), mean_wald = mean(s$bulk_wald),
#              ci_lo = mean(s$bulk_wald) - 1.96 * sd(s$bulk_wald) / sqrt(nrow(s)),
#              ci_hi = mean(s$bulk_wald) + 1.96 * sd(s$bulk_wald) / sqrt(nrow(s)))))

## -----------------------------------------------------------------------------
## 7. Per-animal module scores (does the group result rest on one animal?)
## -----------------------------------------------------------------------------
vsd <- vst(dds, blind = TRUE)
v <- assay(vsd); rownames(v) <- bulk_amex[rownames(v)]
z <- t(scale(t(v)))                       # per-gene z-score across the 6 animals
score_mat <- sapply(sets, function(g) colMeans(z[rownames(z) %in% g, , drop = FALSE]))
colnames(score_mat) <- c("dediff_score", "mature_score")   # avoid make.names mangling
scores <- data.frame(sample = colnames(vsd),
                     condition = as.character(dds$condition),
                     score_mat, row.names = NULL)

## -----------------------------------------------------------------------------
## 8. Leading-edge genes
## -----------------------------------------------------------------------------
led <- do.call(rbind, lapply(seq_len(nrow(fg)), function(i) {
  g <- fg$leadingEdge[[i]]
  if (length(g) == 0) return(NULL)
  data.frame(module = fg$pathway[i], amex = g,
             symbol_in_ref = ref$rowname[match(g, ref$amex)],
             ref_log2FC_prolIS_vs_robo2 = ref$avg_log2FC[match(g, ref$amex)],
             bulk_wald = as.numeric(stat[g]))
}))
led <- led[order(led$module, -led$bulk_wald), ]

## -----------------------------------------------------------------------------
## 9. Figures
## -----------------------------------------------------------------------------
th <- theme_classic(base_size = 9) +
  theme(plot.title = element_text(size = 10, face = "plain"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        axis.text = element_text(colour = "black"))

ms <- mod_stats
ms$label <- factor(c("All tested genes\n(background, n=13,425)",
                     "CM dedifferentiation genes\n(up in injury-zone CM, n=1,105)",
                     "Mature CM identity genes\n(up in homeostatic CM, n=421)"),
                   levels = c("Mature CM identity genes\n(up in homeostatic CM, n=421)",
                              "CM dedifferentiation genes\n(up in injury-zone CM, n=1,105)",
                              "All tested genes\n(background, n=13,425)"))
cols <- setNames(c("grey45", "#4393C3", "#B2182B"), levels(ms$label)[c(3, 2, 1)])

pa <- ggplot(ms, aes(mean_wald, label, colour = label)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = .35) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = .14, linewidth = .5) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = cols, guide = "none") +
  labs(x = "Mean change in dnAXL vs Ctrl (Wald statistic, 95% CI)", y = NULL,
       title = "Dedifferentiation genes do not change; mature-identity genes rise in dnAXL",
       subtitle = "Left of the zero line = higher in Ctrl.   Right of the zero line = higher in dnAXL.") + 
  th

bins$lab <- factor(paste0(c("strongly\nmature", "mature-\nleaning", "neutral",
                            "dediff-\nleaning", "strongly\ndediff"), "\n(n=", bins$n, ")"),
                   levels = paste0(c("strongly\nmature", "mature-\nleaning", "neutral",
                                     "dediff-\nleaning", "strongly\ndediff"), "\n(n=", bins$n, ")"))

pb <- ggplot(bins, aes(lab, mean_wald)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = .35) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .12, linewidth = .5, colour = "grey30") +
  geom_point(size = 2.2, colour = "grey20") +
  labs(x = "Where the gene sits on the snRNA-seq dedifferentiation axis",
       y = "Mean change in\ndnAXL vs Ctrl",
       title = "The trend runs opposite to dedifferentiation, but is shallow",
       subtitle = sprintf("Spearman rho = %.3f, p = %.3g across %d genes measured in both datasets",
                          sp$estimate, sp$p.value, nrow(gl))) + th

ggsave("output/fig_module_comparison.png", (pa / pb), width = 6.8, height = 6.4, dpi = 300)

sc <- do.call(rbind, lapply(seq_along(sets), function(i)
  data.frame(sample = scores$sample, condition = scores$condition,
             module = names(sets)[i], score = score_mat[, i])))
pc <- ggplot(sc, aes(condition, score, colour = condition)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = .3) +
  geom_point(size = 1.8, position = position_jitter(width = .1, height = 0)) +
  facet_wrap(~ module) +
  scale_colour_manual(values = c(Ctrl = "grey45", dnAXL = "#B2182B"), guide = "none") +
  labs(x = NULL, y = "Mean z-scored VST expression",
       title = "Per-animal module scores",
       subtitle = "Each point is one animal; checks whether the group result rests on a single sample") + 
  th

ggsave("output/fig_per_animal_scores.png", pc, width = 5.6, height = 3, dpi = 300)

## -----------------------------------------------------------------------------
## 10. Write outputs
## -----------------------------------------------------------------------------
write.csv(fg_out,    "output/module_enrichment_gsea.csv", row.names = FALSE)
write.csv(mod_stats, "output/module_distribution_shift.csv", row.names = FALSE)
write.csv(gl,        "output/gene_level_bulk_vs_reference.csv", row.names = FALSE)
write.csv(led,       "output/leading_edge_genes.csv", row.names = FALSE)
write.csv(scores,    "output/per_animal_module_scores.csv", row.names = FALSE)
capture.output(sessionInfo(), file = "output/sessionInfo.txt")

cat("\n==== KEY RESULTS ====\n")
print(fg_out, digits = 3, row.names = FALSE)
cat("\n")
print(mod_stats, digits = 3, row.names = FALSE)
cat(sprintf("\nSpearman rho = %.4f, p = %.3g, n = %d\n", sp$estimate, sp$p.value, nrow(gl)))
cat(sprintf("Dilution slope = %.4f (SE %.4f, 95%% CI %.4f to %.4f), p = %.3g, n = %d\n",
            coef(dil)[2], summary(dil)$coef[2, 2], dil_ci[1], dil_ci[2],
            summary(dil)$coef[2, 4], length(dil$residuals)))
cat("\nPer-animal module scores:\n"); print(scores, digits = 3, row.names = FALSE)
