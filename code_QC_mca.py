
# Quality control (QC) for MCA data #
######################################


# Packages
import scanpy




# Pmport AnnData object (test data)
dades_subset = scanpy.read("MCA1.1_adata_Peripheral.h5ad")


# Shape
print("Print shape:")
print(dades_subset.shape)


# Genes that yield the highest fraction of counts in each single cell, 
# across all cells.
scanpy.pl.highest_expr_genes(dades_subset, n_top=20, )



#  assemble some information about mitochondrial gene
dades_subset.var['mt'] = dades_subset.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
scanpy.pp.calculate_qc_metrics(dades_subset, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


# Violin plot of some of the computed quality measures
# - the number of genes expressed in the count matrix
# - the total counts per cell
# - the percentage of counts in mitochondrial genes

scanpy.pl.violin(dades_subset, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)


# Scatter plot:
scanpy.pl.scatter(dades_subset, x='total_counts', y='n_genes_by_counts')
