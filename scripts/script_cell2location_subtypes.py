import sys
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.io import mmread
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

# solve the jax and jaxlib issue from https://github.com/google/jax/issues/5501 for mac pro M1
import cell2location
import scvi
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns

## output folders
results_folder = '../results/visium_axolotl_R12830_resequenced_20220308/cell2location_subtype_out'
# create paths and names to results folders for reference regression and cell2location models
ref_run_name = f'{results_folder}/reference_signatures'
run_name = f'{results_folder}/cell2location_map'

## process snRNA-seq 
original_counts = np.round(csr_matrix(mmread("../data/snRNAseq_countMatrix.mtx")))
original_meta = pd.read_csv("../data/snRNAseq_countMatrix_metadata.csv", index_col = 0)
original_genes = pd.read_csv("../data/snRNAseq_countMatrix_gene.csv", index_col = 0)

adata_ref = ad.AnnData(original_counts, obs=original_meta, var = original_genes)

# filter the object
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

adata_ref = adata_ref[:, selected].copy()

# prepare anndata for the regression model
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        batch_key=None,
                        # cell type, covariate used for constructing signatures
                        labels_key='subtypes',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        categorical_covariate_keys=['condition']
                       )

# create the regression model
from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)

# Use all data for training (validation not implemented yet, train_size=1)
# takes ~ 3.5h for 175 epochs, but should increase to 200/250 epochs
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=False)

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': False}
)

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc_subtypes_v1.h5ad"
adata_ref.write(adata_file)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']

### load visium data 
vis = sc.read_loom("../data/visium_Amex_all.loom") # all 4 slices merged
vis.X = np.round(vis.X)

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(vis.var_names, inf_aver.index)
vis = vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=vis, 
                                                 categorical_covariate_keys=['condition'],
                                                 continuous_covariate_keys = ["nCount_Spatial"])

# create and train the model
mod = cell2location.models.Cell2location(
    vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=5,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
#mod.view_anndata_setup()

mod.train(max_epochs=20000, # between 15000 to 30000
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=False)

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
vis = mod.export_posterior(
    vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': False}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
vis.write(adata_file)

# add 5% quantile, representing confident cell abundance, 'at least this amount is present',
# to adata.obs with nice names for plotting
vis.obs[vis.uns['mod']['factor_names']] = vis.obsm['q05_cell_abundance_w_sf']
vis.obs.to_csv(results_folder+"/predictions_cell2loc.csv")

