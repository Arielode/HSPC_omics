#~/bin/Python-3.7/bin/python3
# Using scanpy to analyze scRNA-seq data
# input data is restored in data matrix format
###

import scanpy as sc
import pandas as pd
import numpy as np

def filterCellsByGroup(data, feature, fd=3):
    cells_remove = []
    group = np.unique(data.obs[feature])
    for gp in group:
        tmpData = data[data.obs[feature] == gp, :]
        tmpData.var['mt'] = tmpData.var_names.str.startswith('MT-')
        tmpData.var['rbs'] = tmpData.var_names.str.startswith(('RPS', 'RPL'))
        sc.pp.calculate_qc_metrics(tmpData, qc_vars=['mt', 'rbs'], percent_top=None, log1p=False, inplace=True)
        gene_counts_mean = np.mean(tmpData.obs.n_genes_by_counts)  # average gene number detected
        gene_counts_sd = np.std(tmpData.obs.n_genes_by_counts)  # standard deviation of gene numbers
        mt_mean = np.mean(tmpData.obs.pct_counts_mt)  # average mt-genes
        mt_sd = np.std(tmpData.obs.pct_counts_mt)  # standard deviation of mt-genes
        rbs_mean = np.mean(tmpData.obs.pct_counts_rbs)  # average rbs-genes
        rbs_sd = np.std(tmpData.obs.pct_counts_rbs)  # standard deviation of rbs-genes
        cells_remove += tmpData.obs_names[(tmpData.obs.n_genes_by_counts < (gene_counts_mean - fd * gene_counts_sd)) | (
                tmpData.obs.n_genes_by_counts > (gene_counts_mean + fd * gene_counts_sd))].to_list()
        cells_remove += tmpData.obs_names[tmpData.obs.pct_counts_mt > mt_mean + fd * mt_sd].to_list()
        cells_remove += tmpData.obs_names[tmpData.obs.pct_counts_rbs > rbs_mean + fd * rbs_sd].to_list()
    return cells_remove

# setting
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
#sc.logging.print_header()
#sc.settings.set_figure_params(dpi=80, facecolor='white')
input_file = 'data/RNA.h5ad'
results_file = 'write/RNA.h5ad'  # the file that will store the analysis results

adata = sc.read_h5ad(input_file)  #1271 ¡Á 19836
print(adata)
cells_rm = filterCellsByGroup(data=adata, feature='celltype', fd=3)
adata = adata[[(i not in cells_rm) for i in adata.obs_names], :]  # 1216 ¡Á 19836
print(adata)

# highly expressed genes
# sc.pl.highest_expr_genes(adata, n_top=40)

# annotate the group of mitochondrial genes as 'mt'
# annotate the group of ribosome genes as 'rbs'
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['rbs'] = adata.var_names.str.startswith(('RPS', 'RPL'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'rbs'], percent_top=None, log1p=False, inplace=True)
# cell cycle score calculation
cell_cycle_genes = [x.strip() for x in open('./data/regev_lab_cell_cycle_genes.txt')]
#cell_cycle_genes = [x.strip() for x in open('./data/cyclegenes.txt')] # modified by zhucy
s_genes = cell_cycle_genes[:43]  # s phage genes
g2m_genes = cell_cycle_genes[43:]  # g2m phase genes
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# a violin plot of some of the computed quality measures:
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rbs'], jitter=0.4, multi_panel=True,show=False,save="quality.features.ng")

# # distribution of measures:
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',show=False,save="distribution-counts-mt-pct.png")
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',show=False,save="distribution-counts-ngene.png")

# remove ribosome and mitochondrial genes
adata = adata[:, (~ adata.var['rbs'])]
adata = adata[:, (~ adata.var['mt'])]
print(adata)  # 1216 ¡Á 19722

## filter genes expressed in less than 3 cells
sc.pp.filter_genes(adata, min_cells=10)
print(adata)   # 1216 ¡Á 18526 
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

# normalize data
sc.pp.normalize_total(adata, target_sum=1e4)
# logarithmize data
sc.pp.log1p(adata)

# identify highly-variable genes
sc.pp.highly_variable_genes(adata,n_top_genes=2000)
adata.raw = adata

# filtering out HVGs 
adata = adata[:, adata.var.highly_variable]  

# Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed.
sc.pp.regress_out(adata, ['total_counts', 'S_score', 'G2M_score'])

# Scale each gene to unit variance. Clip values exceeding standard deviation 10.
sc.pp.scale(adata, max_value=10)

# # Principal component analysis
sc.tl.pca(adata, svd_solver='arpack')

# Inspect contribution of PCs
sc.pl.pca_variance_ratio(adata, log=True)

# integration with BBKNN
sc.external.pp.bbknn(adata, batch_key='individual', n_pcs=17)
for i in np.arange(1.0,1.9,0.05):
   j = round(i,3)
   sc.tl.leiden(adata, resolution=j, key_added='leiden'+str(j))
   
color_cluster = ['leiden' + str(round(i, 3)) for i in np.arange(1.0, 1.9, 0.05)]
color_cluster.extend(['lastleiden', 'celltype','phase','individual'])

res=1.75
dist=0.3
spread=2.0
out=str(res)+'_'+str(dist)+'_'+str(spread)
sc.tl.leiden(adata, resolution=res)
sc.tl.umap(adata, spread=spread, min_dist=dist)
sc.pl.umap(adata, color=['leiden', 'celltype', 'individual', 'lastleiden'], legend_loc='on data', show=False, save=out+'.png')
adata.write(out+'.h5ad')
sc.pl.umap(adata, color=['leiden', 'individual', 'celltype', 'phase'], legend_loc='on data',save=out+'.pdf')

