#!/usr/bin/env python
# coding: utf-8
from configparser import ConfigParser

import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as an
import numpy as np

from mftools.cellgene import create_scanpy_object
from mftools.fileio import MerfishAnalysis


config:ConfigParser = ConfigParser()
config.read('config.ini')
ioconf = config["IO Options"]
csconf = config["2 Label Transfer"] # TODO: rename to 1_Cell_Segmentation

MERSCOPE_DIR = ioconf['msdir']
EXPERIMENT_NAME = ioconf['experiment']
REFERENCE_DATA_PATH = "/THIS IS NOT A REAL PATH"
MER_RAWDATA_DIR = "data"
MER_OUTPUT_DIR = "output"

CELLTYPE_KEY = "SOMETHING ELSE HERE"
DATSET_KEY = 'org_datset'

# Output variables
OUTPUT_DIRECTORY = ""
HI_DEF_DPI = 1000
# TODO: Scale function to get spatial embedding dots close to actuall cell size


image_dataset_path = f"{MERSCOPE_DIR}/{MER_RAWDATA_DIR}/{EXPERIMENT_NAME}/" # Used to load imageset
expiriment_out_path = f"{MERSCOPE_DIR}/{MER_OUTPUT_DIR}/{EXPERIMENT_NAME}/" # Used to find barcodes 

output = MerfishAnalysis(expiriment_out_path)

# Read in data
mfdata = create_scanpy_object(output)
refdata = sc.read(REFERENCE_DATA_PATH)
print(mfdata)
print(refdata)


var_intersect = mfdata.var_names.intersection(refdata.var_names)
mfdata = mfdata[:, var_intersect].copy()
refdata = refdata[:, var_intersect].copy()


sc.pp.filter_cells(mfdata, min_genes=5)
sc.pp.filter_genes(mfdata, min_cells=3)
sc.pp.normalize_total(mfdata)
sc.pp.log1p(mfdata)


adata = an.concat([mfdata, refdata], label=DATSET_KEY, join='outer')


sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)



sc.pl.umap(adata, color=DATSET_KEY)





sc.external.pp.harmony_integrate(adata, key=DATSET_KEY, max_iter_harmony=20)





sc.pp.neighbors(adata, use_rep="X_pca_harmony")
sc.tl.umap(adata)





sc.pl.umap(adata, color=DATSET_KEY)





sc.pl.umap(adata, color=CELLTYPE_KEY)





from sklearn.neighbors import KNeighborsClassifier
traindata = adata[adata.obs[DATSET_KEY] == '0'] 
nn = KNeighborsClassifier(n_jobs=16)
nn.fit(traindata.obsm["X_pca_harmony"], traindata.obs[CELLTYPE_KEY])
pred = nn.predict(adata[adata.obs[DATSET_KEY] == '0'].obsm["X_pca_harmony"])





adata.obs.loc[adata.obs[DATSET_KEY] == '0', CELLTYPE_KEY] = pred
adata = adata[adata.obs[DATSET_KEY] == '0']




sc.pl.umap(adata, color=CELLTYPE_KEY)





# Attach celltype and age labels to each cell
# Celltype label
mfdata.obs["celltype"] = list(adata[adata.obs[DATSET_KEY] == "0"].obs[CELLTYPE_KEY])

sc.set_figure_params(dpi=HI_DEF_DPI)
sc.pl.embedding(mfdata, basis="X_spatial", color="celltype", frameon=False, title="")#, groups=["Kupffer cell"])
sc.pl.embedding(mfdata, basis="X_spatial", color="age", frameon=False, title="")#, groups=["Kupffer cell"])

for celltype in mfdata.obs["celltype"].unique():
    sc.pl.embedding(mfdata[(mfdata.obsm["X_spatial"][:,0] > -1500) & (mfdata.obsm["X_spatial"][:,0] < 1500) & (mfdata.obsm["X_spatial"][:,1] > -1500) & (mfdata.obsm["X_spatial"][:,1] < 1500)], basis="X_spatial", color="celltype", frameon=False, title=celltype, groups=[celltype], legend_loc="off")


adata.write_h5ad("data_raw/M126to127_test2.h5ad")

