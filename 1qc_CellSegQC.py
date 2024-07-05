#!/usr/bin/env python
# coding: utf-8
from configparser import ConfigParser
from pathlib import Path

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


from mftools.fileio import MerfishAnalysis
from mftools.cellgene import create_scanpy_object
from mftools.plotting import fov_show


### TODO: THIS has been copied from 1_CellSeg.py, make sure to replace it here as well 
config:ConfigParser = ConfigParser()
config.read('config.ini')
ioconf = config["IO Options"]
csconf = config["Cellpose"] # TODO: rename to 1_Cell_Segmentation

# Filepath naming (Ren9)
MERSCOPE_DIR = ioconf['msdir']
EXPERIMENT_NAME = ioconf['experiment']

# Script params
ZSLICE = csconf['zslice']
CHANNEL = csconf['channel']
MODEL= ""
TEST_FOV = 150

# PATH PARAMS 
# CHANGE ONLY IF NEEDED: these should be standardized in the future and creation of these dirs automated from
# globalized config file
MER_RAWDATA_DIR = "data"
MER_OUTPUT_DIR = "output"
POSITIONS = "settings/positions.csv" # TODO: for testing, change later
CELLPOSE_DIR = f"./cellpose_{EXPERIMENT_NAME}" # TODO: for testing, change later 

# Path assembly TODO: remove this when we can trust MERBOT to generate these for
image_dataset_path = f"{MERSCOPE_DIR}/{MER_RAWDATA_DIR}/{EXPERIMENT_NAME}/" # Used to load imageset
# expiriment_out_path = f"{MERSCOPE_DIR}/{MER_OUTPUT_DIR}/{EXPERIMENT_NAME}/"
cellpose_path = CELLPOSE_DIR # Used to load masks
positions_path = f"{image_dataset_path}/{POSITIONS}" # Used to load positions
fig_outpath = f"./fig_dump_{EXPERIMENT_NAME}"

print(f"Looking for images in {image_dataset_path}")
print(f"Looking for positions in {positions_path}")


output = MerfishAnalysis(cellpose_path)
imageset = ImageDataset(image_dataset_path)
seg = CellSegmentation(imagedata=imageset, output=output, channel=CHANNEL, zslice=3)
adata = create_scanpy_object(output)



# Fix FOV coordinates:
positions = pd.read_csv(POSITIONS, header=None)

gx = 220 * ((2048 - adata.obs["fov_x"]) / 2048) + positions.iloc[adata.obs["fov"]][1].values
gy = 220 * ((2048 - adata.obs["fov_y"]) / 2048) + positions.iloc[adata.obs["fov"]][0].values

adata.obs["global_x"] = gx
adata.obs["global_y"] = gy
adata.obsm["X_spatial"] = adata.obs[["global_x", "global_y"]].to_numpy()

# TODO: this is useful, put this in MFtools
def rotate_spatial_coordinates(adata):
    X_spatial = adata.obsm['X_spatial']
    rotation_matrix = np.array([[0, -1],
                                    [1, 0]])
    X_spatial_rotated = X_spatial.dot(rotation_matrix)
    adata.obsm['X_spatial'] = X_spatial_rotated

rotate_spatial_coordinates(adata)


sc.pl.embedding(adata, basis="spatial", frameon=False)


from mftools.fileio import ImageDataset
from mftools.segmentation import CellSegmentation
import matplotlib.pyplot as plt
import numpy as np


dapi = imageset.load_image(channel=CHANNEL, fov=TEST_FOV, max_projection=True)
mask = seg[TEST_FOV]
plt.figure(dpi=150)
plt.imshow(dapi, cmap="gray")
plt.contour(mask, [x+0.5 for x in np.unique(mask)], colors="tab:blue")
plt.axis("off")



# TODO: is there somewhere smarter to place this where it can be reused? or accessed more easily?
# Maybe some sort of merfish experiment summary variable
def segqc_summary(adata, show=True):
    nCells = len(adata.obs)
    nTSCP = adata.obs['total_counts'].sum()
    medTSCP_Cell = adata.obs['total_counts'].median()
    medGenes_Cell = adata.obs['n_genes_by_counts'].median()

    if show:
        print(f"nCells = {nCells}")
        print(f"nTSCP = {nTSCP}")
        print(f"medTSCP_Cell = {medTSCP_Cell}")
        print(f"medGenes_Cell = {medGenes_Cell}")


# calculate quality control metrics
sc.pp.calculate_qc_metrics(adata, percent_top=(50, 100, 200, 300), inplace=True)

def segqc_plots(adata):
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))

    axs[0].set_title(f"{EXPERIMENT_NAME} - Total transcripts per cell")
    sns.histplot(
        adata.obs["total_counts"],
        kde=False,
        ax=axs[0],
    )

    axs[1].set_title(f"{EXPERIMENT_NAME} - Unique transcripts per cell")
    sns.histplot(
        adata.obs["n_genes_by_counts"],
        kde=False,
        ax=axs[1],
    )

    axs[2].set_title(f"{EXPERIMENT_NAME} - Volume of segmented cells")
    sns.histplot(
        adata.obs["volume"],
        kde=False,
        ax=axs[2],
    )

    fig.savefig(f"{fig_outpath}/before")

sc.pp.filter_cells(adata, min_genes=5)
sc.pp.filter_genes(adata, min_cells=3)





import matplotlib 
samples = [adataR001, adataR002, adataR003, adataR004, adataR005, adataR006, adataR007, adataR008, adataR009, adataR010, adataR011, adataR012, adataR013]
names = ['Ren001', 'Ren002', 'Ren003', 'Ren004', 'Ren005', 'Ren006', 'Ren007', 'Ren008', 'Ren009', 'Ren010', 'Ren011', 'Ren012', 'Ren013']
for n, adata in enumerate(samples):
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))

    axs[0].set_title(f"{names[n]} - Total transcripts per cell")
    sns.histplot(
        adata.obs["total_counts"],
        kde=False,
        ax=axs[0],
    )

    axs[1].set_title(f"{names[n]} - Unique transcripts per cell")
    sns.histplot(
        adata.obs["n_genes_by_counts"],
        kde=False,
        ax=axs[1],
    )

    axs[2].set_title(f"{names[n]} - Volume of segmented cells")
    sns.histplot(
        adata.obs["volume"],
        kde=False,
        ax=axs[2],
    )


# # QC Metrics




samples = [adataR001, adataR002, adataR003, adataR004, adataR005, adataR006, adataR007, adataR008, adataR009, adataR010, adataR011, adataR012, adataR013]
names = ['REN001', 'REN002', 'REN003', 'REN004', 'REN005', 'REN006', 'REN007', 'REN008', 'REN009', 'REN010', 'REN011', 'REN012', 'REN013']

for n, sample in enumerate(samples):
    nCells = len(sample.obs)
    nTSCP = sample.obs['total_counts'].sum()
    medTSCP_Cell = sample.obs['total_counts'].median()
    medGenes_Cell = sample.obs['n_genes_by_counts'].median()
    # print(f"{names[n]} - nCells = {nCells}")
    # print(f"{names[n]} - nTSCP = {nTSCP}")
    # print(f"{names[n]} - medTSCP_Cell = {medTSCP_Cell}")
    print(f"{names[n]} - medGenes_Cell = {medGenes_Cell}")





import scanpy as sc
# adata = sc.read_h5ad('/mnt/merfish18/BICAN/objects/BICAN_R1toR17_pre.h5ad')
# sc.settings.verbosity = 3
# sc.pp.filter_cells(adata, min_genes=5)
# sc.pp.filter_genes(adata, min_cells=3)

adataR001 = adata[adata.obs['batch'] == 'REN001']
adataR002 = adata[adata.obs['batch'] == 'REN002']
adataR003 = adata[adata.obs['batch'] == 'REN003']
adataR004 = adata[adata.obs['batch'] == 'REN004']
adataR005 = adata[adata.obs['batch'] == 'REN005']
adataR006 = adata[adata.obs['batch'] == 'REN006']
adataR007 = adata[adata.obs['batch'] == 'REN007']
adataR008 = adata[adata.obs['batch'] == 'REN008']
adataR009 = adata[adata.obs['batch'] == 'REN009']
adataR010 = adata[adata.obs['batch'] == 'REN010']
adataR011 = adata[adata.obs['batch'] == 'REN011']
adataR012 = adata[adata.obs['batch'] == 'REN012']
adataR013 = adata[adata.obs['batch'] == 'REN013']
adataR014 = adata[adata.obs['batch'] == 'REN014']
adataR015 = adata[adata.obs['batch'] == 'REN015']
adataR016 = adata[adata.obs['batch'] == 'REN016']
adataR017 = adata[adata.obs['batch'] == 'REN017']


samples = [adataR001, adataR002, adataR003, adataR004, adataR005, adataR006, adataR007, adataR008, adataR009, adataR010, adataR011, adataR012, adataR013, adataR014, adataR015, adataR016, adataR017]
names = ['REN001', 'REN002', 'REN003', 'REN004', 'REN005', 'REN006', 'REN007', 'REN008', 'REN009', 'REN010', 'REN011', 'REN012', 'REN013', 'REN014', 'REN015', 'REN016', 'REN017']

for n, sample in enumerate(samples):
    nCells = len(sample.obs)
    nTSCP = sample.obs['total_counts'].sum()
    medTSCP_Cell = sample.obs['total_counts'].median()
    medGenes_Cell = sample.obs['n_genes_by_counts'].median()
    # print(f"{names[n]} - nCells = {nCells}")
    # print(f"{names[n]} - nTSCP = {nTSCP}")
    # print(f"{names[n]} - medTSCP_Cell = {medTSCP_Cell}")
    print(f"{names[n]} - medGenes_Cell = {medGenes_Cell}")


# # Control Genes




import pandas as pd # PERCENTAGE

# Create a dictionary to store the results
results = {'Sample': [], 'Gene': [], 'Percentage': []}
samples = [adataR001, adataR002, adataR003, adataR004, adataR005, adataR006, adataR007, adataR008, adataR009, adataR010, adataR011, adataR012, adataR013, adataR014, adataR015, adataR016, adataR017]
sample_names = ['R001', 'R002', 'R003', 'R004', 'R005', 'R006', 'R007', 'R008', 'R009', 'R010', 'R011', 'R012', 'R013', 'R014', 'R015', 'R016', 'R017']
genes = ['UBR2', 'PAFAH1B1','RBM5', 'PSMD1', 'RBM6', 'SRSF11', 'LUC7L2', 'SNX14', 'PRKDC', 'YME1L1', 'PTCD3', 'SENP5', 'GOSR1', 'PSPC1', 'N4BP2L2', 'YLPM1', 'CUL5', 'SMARCAD1', 'OPA1', 'DCAF10']
# Iterate over each sample
for n, sample in enumerate(samples):  # Assuming samples are named Sample1, Sample2, ..., Sample9
    # Load or access each AnnData object
    adata = sample  # Replace this with your method of loading/accessing AnnData
    
    # Iterate over 20 different genes
    for gene_name in genes:  # Replace with your gene names
        # Check if the gene is present in the AnnData object
        if gene_name in adata.var_names:
            # Get the index of the gene in the var_names array
            gene_index = list(adata.var_names).index(gene_name)
            
            # Count the number of cells expressing the gene
            cells_with_gene = sum(adata.X[:, gene_index] > 0)
            
            # Calculate the percentage
            percentage_cells_with_gene = (cells_with_gene / adata.n_obs) * 100
            
            # Append the results to the dictionary
            results['Sample'].append(f'Sample{sample_names[n]}')
            results['Gene'].append(gene_name)
            results['Percentage'].append(percentage_cells_with_gene)
        else:
            # If gene not found, append NaN
            results['Sample'].append(f'Sample{sample_names[n]}')
            results['Gene'].append(gene_name)
            results['Percentage'].append(float('nan'))

# Convert the dictionary to a DataFrame
df = pd.DataFrame(results)

# Display the DataFrame
print(df)

pivot_df = df.pivot(index='Gene', columns='Sample', values='Percentage')
pivot_df = pivot_df.reindex(genes)


plt.figure(figsize=(10, 6))
sns.heatmap(pivot_df, cmap='viridis', annot=True, fmt=".2f", vmin=0, vmax=70, xticklabels=['R1-UWA.7648.CX22.Pu', 'R2-UWA.7648.CX22.Pu', 'R3-UWA.7648.CX17.CaH', 'R4-UWA.7648.CX20.CaB', 'R5-UWA.7648.CX26.EC', 'R6-UWA.7648.CX21.EC', 'R7-UWA.7648.CX23.THM.01', 'R8-UWA.7648.CX24.M1C', 'R9-UWA.7648.CX10.A46', 'R10-UWA.7648.CX41.V1C', 'R11-UWA.7648.CX31.PPHC', 'R12-UWA.7648.CX04.A10(m)','R13-UWA.7648.CX07.A9', 'R14-UWA.7648.CX22.Pu.R1', 'R15-UWA.7648.CX18.A8', 'R16-UWA.7648.CX17.A44', 'R17-UWA.7648.CX23.THM.01.R1'])
plt.xlabel('Samples')
plt.xticks(rotation=45, ha='right')
plt.ylabel('Control Genes')
plt.title('Percentage of Cells Expressing Control Genes across Samples - Ren Lab')
plt.show()





df.to_csv('BICAN_internalgenes.csv')





adataR001.obs_vector('PAFAH1B1').mean()





import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

adata_objects = [adataR001, adataR002, adataR003, adataR004, adataR005, adataR006, adataR007, adataR008, adataR009, adataR010, adataR011, adataR012, adataR013, adataR014, adataR015, adataR016, adataR017]
sample_names = ['R001', 'R002', 'R003', 'R004', 'R005', 'R006', 'R007', 'R008', 'R009', 'R010', 'R011', 'R012', 'R013', 'R014', 'R015', 'R016', 'R017']
control_genes = ['UBR2', 'PAFAH1B1','RBM5', 'PSMD1', 'RBM6', 'SRSF11', 'LUC7L2', 'SNX14', 'PRKDC', 'YME1L1', 'PTCD3', 'SENP5', 'GOSR1', 'PSPC1', 'N4BP2L2', 'YLPM1', 'CUL5', 'SMARCAD1', 'OPA1', 'DCAF10']

def calculate_mean_expression(control_genes, adata_objects):
    mean_expression_values = []
    for gene in control_genes:
        expression_values = []
        for adata in adata_objects:
            expression_values.append(adata.obs_vector(gene).mean())
        mean_expression_values.append(expression_values)
    return mean_expression_values

# Function to display mean expression values as a table
def display_table(control_genes, mean_expression_values):
    df = pd.DataFrame(mean_expression_values, index=control_genes)
    df.columns = ['REN{}'.format(i+1) for i in range(len(mean_expression_values[0]))]
    print(df)

# Function to display mean expression values as a heatmap
def display_heatmap(control_genes, mean_expression_values):
    plt.figure(figsize=(10, 6))
    sns.heatmap(mean_expression_values, cmap='viridis', annot=True, vmin=0, vmax=3, xticklabels=['R1-UWA.7648.CX22.Pu', 'R2-UWA.7648.CX22.Pu', 'R3-UWA.7648.CX17.CaH', 'R4-UWA.7648.CX20.CaB', 'R5-UWA.7648.CX26.EC', 'R6-UWA.7648.CX21.EC', 'R7-UWA.7648.CX23.THM.01', 'R8-UWA.7648.CX24.M1C', 'R9-UWA.7648.CX10.A46', 'R10-UWA.7648.CX41.V1C', 'R11-UWA.7648.CX31.PPHC', 'R12-UWA.7648.CX04.A10(m)','R13-UWA.7648.CX07.A9', 'R14-UWA.7648.CX22.Pu.R1', 'R15-UWA.7648.CX18.A8', 'R16-UWA.7648.CX17.A44', 'R17-UWA.7648.CX23.THM.01.R1'])
    plt.xlabel('Samples')
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('Control Genes')
    plt.title('Mean Expression of Control Genes across Samples - Ren Lab')
    plt.show()

control_genes = control_genes


# Calculate mean expression values
mean_expression_values = calculate_mean_expression(control_genes, adata_objects)

# Display mean expression values as a table
display_table(control_genes, mean_expression_values)

# Display mean expression values as a heatmap
display_heatmap(control_genes, mean_expression_values)


# ## Wei's results




control_genes = ['UBR2', 'PAFAH1B1','RBM5', 'PSMD1', 'RBM6', 'SRSF11', 'LUC7L2', 'SNX14', 'PRKDC', 'YME1L1', 'PTCD3', 'SENP5', 'GOSR1', 'PSPC1', 'N4BP2L2', 'YLPM1', 'CUL5', 'SMARCAD1', 'OPA1', 'DCAF10']





Ecker_20mean = pd.read_csv('/home/jeolness/Documents/References/Collab_data/ic20.mean.csv')
Ecker_20pcnt = pd.read_csv('/home/jeolness/Documents/References/Collab_data/ic20.pcnt.csv')
Ecker_20mean.set_index(Ecker_20mean.columns[0], inplace=True)
Ecker_20pcnt.set_index(Ecker_20pcnt.columns[0], inplace=True)


Ecker_20mean_reordered = Ecker_20mean.reindex(control_genes)
Ecker_20pcnt_reordered = Ecker_20pcnt.reindex(control_genes)
Ecker_20pcnt_reordered = Ecker_20pcnt_reordered.multiply(100)


Ecker_20pcnt_reordered






import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(10, 6))
sns.heatmap(Ecker_20pcnt_reordered, cmap='viridis', annot=True, vmin=0, vmax=70, xticklabels=['Ecker: UWA.7648.CX21.EC-500GP', 'Ecker: UWA.7648.CX26.EC-500GP', 'Ecker: UWA.7648.CX17.CaH-', 'Ecker: UWA.7648.CX20.CaB', 'Ecker: UWA.7648.CX23.THM.01-500GP'])
plt.xlabel('Samples')
plt.xticks(rotation=45, ha='right')
plt.ylabel('Control Genes')
plt.title('Percentage of Cells Expressing Control Genes across Samples - Ecker Lab')
plt.show()





import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(10, 6))
sns.heatmap(Ecker_20mean_reordered, cmap='viridis', annot=True, vmin=0, vmax=3, xticklabels=['Ecker: UWA.7648.CX21.EC-500GP', 'Ecker: UWA.7648.CX26.EC-500GP', 'Ecker: UWA.7648.CX17.CaH-', 'Ecker: UWA.7648.CX20.CaB', 'Ecker: UWA.7648.CX23.THM.01-500GP'])
plt.xlabel('Samples')
plt.xticks(rotation=45, ha='right')
plt.ylabel('Control Genes')
plt.title('Mean Expression of Control Genes across Samples - Ecker Lab')
plt.show()


# # Clustering




import scanpy as sc
adata = sc.read_h5ad('/mnt/merfish18/BICAN/objects/BICAN_R1toR17_pre.h5ad')





sc.settings.verbosity = 3
sc.pp.filter_cells(adata, min_genes=5)
sc.pp.filter_genes(adata, min_cells=3)





adataR001 = adata[adata.obs['batch'] == 'REN001']
adataR002 = adata[adata.obs['batch'] == 'REN002']
adataR003 = adata[adata.obs['batch'] == 'REN003']
adataR004 = adata[adata.obs['batch'] == 'REN004']
adataR005 = adata[adata.obs['batch'] == 'REN005']
adataR006 = adata[adata.obs['batch'] == 'REN006']
adataR007 = adata[adata.obs['batch'] == 'REN007']
adataR008 = adata[adata.obs['batch'] == 'REN008']
adataR009 = adata[adata.obs['batch'] == 'REN009']
adataR010 = adata[adata.obs['batch'] == 'REN010']
adataR011 = adata[adata.obs['batch'] == 'REN011']
adataR012 = adata[adata.obs['batch'] == 'REN012']
adataR013 = adata[adata.obs['batch'] == 'REN013']
adataR014 = adata[adata.obs['batch'] == 'REN014']
adataR015 = adata[adata.obs['batch'] == 'REN015']
adataR016 = adata[adata.obs['batch'] == 'REN016']
adataR017 = adata[adata.obs['batch'] == 'REN017']





pu = adataR001.concatenate([adataR002, adataR014], batch_categories=['REN001', 'REN002', 'REN014'])
cah = adataR003
cab = adataR004
thm = adataR007.concatenate([adataR017], batch_categories=['REN007', 'REN017'])
ec = adataR005.concatenate([adataR006], batch_categories=['REN005', 'REN006'])
m1c = adataR008
a46 = adataR009
v1c = adataR010
pphc = adataR011
a10 = adataR012
a9 = adataR013
a8 = adataR015
a44 = adataR016





regions = [pu, cah, cab, ec, thm, m1c, a46, v1c, pphc, a10, a9, a8, a44]

for region in regions:
    # sc.pp.highly_variable_genes(region, flavor="seurat_v3", n_top_genes=4000)
    sc.pp.normalize_total(region, inplace=True)
    sc.pp.log1p(region)
    sc.pp.pca(region)
    sc.pp.neighbors(region, n_neighbors=30, metric='correlation')
    sc.tl.umap(region, min_dist=0.001, spread=3.0)
    sc.tl.leiden(region)





# Leiden Clustering Visualization

regions = [pu, cah, cab, ec, thm, m1c, a46, v1c, pphc, a10, a9, a8, a44]
names = ['Pu', 'CaH', 'CaB', 'EC', 'Thm', 'M1C', 'A46', 'V1C', 'PPHC', 'A10', 'A9', 'A8', 'A44']

for n, region in enumerate(regions):
    sc.pl.umap(
        region,
        color=[
            "total_counts",
            "n_genes_by_counts",
            "leiden",
        ],
        wspace=0.4,
        title=f"{names[n]}, Total Counts"
    )





adataR001 = pu[pu.obs['batch'] == 'REN001']
adataR002 = pu[pu.obs['batch'] == 'REN002']
adataR003 = cah[cah.obs['batch'] == 'REN003']
adataR004 = cab[cab.obs['batch'] == 'REN004']
adataR005 = ec[ec.obs['batch'] == 'REN005']
adataR006 = ec[ec.obs['batch'] == 'REN006']
adataR007 = thm[thm.obs['batch'] == 'REN007']
adataR008 = m1c[m1c.obs['batch'] == 'REN008']
adataR009 = a46[a46.obs['batch'] == 'REN009']
adataR010 = v1c[v1c.obs['batch'] == 'REN010']
adataR011 = pphc[pphc.obs['batch'] == 'REN011']
adataR012 = a10[a10.obs['batch'] == 'REN012']
adataR013 = a9[a9.obs['batch'] == 'REN013']
adataR014 = pu[pu.obs['batch'] == 'REN014']
adataR015 = a8[a8.obs['batch'] == 'REN015']
adataR016 = a44[a44.obs['batch'] == 'REN016']
adataR017 = thm[thm.obs['batch'] == 'REN017']





adataR001





import numpy as np
def rotate_spatial_coordinates(adata):
    X_spatial = adata.obsm['X_spatial']
    rotation_matrix = np.array([[0, -1],
                                    [1, 0]])
    X_spatial_rotated = X_spatial.dot(rotation_matrix)
    adata.obsm['X_spatial'] = X_spatial_rotated

rotate_spatial_coordinates(adataR004)





import squidpy as sq
samples = [adataR001, adataR002, adataR003, adataR004, adataR005, adataR006, adataR007, adataR008, adataR009, adataR010, adataR011, adataR012, adataR013, adataR014, adataR015, adataR016, adataR017]
names = ['REN001', 'REN002', 'REN003', 'REN004', 'REN005', 'REN006', 'REN007', 'REN008', 'REN009', 'REN010', 'REN011', 'REN012', 'REN013', 'REN014', 'REN015', 'REN016', 'REN017']

for n, sample in enumerate(samples):
    sq.pl.spatial_scatter(
        sample,
        shape=None,
        color=[
            "leiden",
        ],
        wspace=0.4,
        spatial_key='X_spatial',
        frameon=False,
        figsize=[5, 10],
        title=f"{names[n]}"
    )





adata = adataR001.concatenate([adataR002, adataR003, adataR004, adataR005, adataR006, adataR007, adataR008, adataR009, adataR010, adataR011, adataR012, adataR013, adataR014, adataR015, adataR016, adataR017], batch_categories=['REN001', 'REN002', 'REN003', 'REN004', 'REN005', 'REN006', 'REN007', 'REN008', 'REN009', 'REN010', 'REN011', 'REN012', 'REN013', 'REN014', 'REN015', 'REN016', 'REN017'])





adata.write('/mnt/merfish18/BICAN/objects/BICAN_R1toR17_post.h5ad')


# # Correlation of genes in Ren2 and Ren8/Ren9




adataR002 = adata[adata.obs['batch'] == 'REN002']
adataR008 = adata[adata.obs['batch'] == 'REN008']
adataR009 = adata[adata.obs['batch'] == 'REN009']





genes002 = adataR002.X.toarray()
genes008 = adataR008.X.toarray()
genes009 = adataR009.X.toarray()





adata002_mean = genes002.mean(axis=0)
adata008_mean = genes008.mean(axis=0)
adata009_mean = genes009.mean(axis=0)





import scanpy as sc
import pandas as pd

import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import scipy
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression

# Pearson correlation
correlation_coefficient, _ = pearsonr(adata002_mean, adata008_mean)# , adata2_mean.flatten())
# Plot scatterplot
plt.scatter(adata002_mean, adata008_mean, alpha=0.7)
plt.xlabel('Human Brain MERFISH Datasets R002')
plt.ylabel('Human Brain MERFISH Datasets R008')

# Regression line
slope, intercept = np.polyfit(adata002_mean, adata008_mean, 1)
plt.plot(adata002_mean, slope * adata002_mean + intercept, color='red', linestyle=':', label='Regression Line')
# # Equality line (y = x)
# plt.plot(adata1_mean, adata1_mean, color='gray', linestyle='--', label='Equality Line (y=x)')
plt.title(f'Pearson Correlation: {correlation_coefficient:.2f}')
plt.legend()
plt.show()





import scanpy as sc
import pandas as pd

import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import scipy
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression

# Pearson correlation
correlation_coefficient, _ = pearsonr(adata002_mean, adata009_mean)# , adata2_mean.flatten())
# Plot scatterplot
plt.scatter(adata002_mean, adata009_mean, alpha=0.7)
# plt.title(f'Pearson Correlation: {correlation_coefficient:.2f}')
plt.xlabel('Human Brain MERFISH Datasets R002')
plt.ylabel('Human Brain MERFISH Datasets R009')

# Regression line
slope, intercept = np.polyfit(adata002_mean, adata009_mean, 1)
plt.plot(adata002_mean, slope * adata002_mean + intercept, color='red', linestyle=':', label='Regression Line')
# # Equality line (y = x)
# plt.plot(adata1_mean, adata1_mean, color='gray', linestyle='--', label='Equality Line (y=x)')
plt.title(f'Pearson Correlation: {correlation_coefficient:.2f}')
plt.legend()
plt.show()





import scanpy as sc
import pandas as pd

import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import scipy
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression

# Pearson correlation
correlation_coefficient, _ = pearsonr(adata008_mean, adata009_mean)# , adata2_mean.flatten())
# Plot scatterplot
plt.scatter(adata008_mean, adata009_mean, alpha=0.7)
# plt.title(f'Pearson Correlation: {correlation_coefficient:.2f}')
plt.xlabel('Human Brain MERFISH Datasets R008')
plt.ylabel('Human Brain MERFISH Datasets R009')

# Regression line
slope, intercept = np.polyfit(adata008_mean, adata009_mean, 1)
plt.plot(adata008_mean, slope * adata008_mean + intercept, color='red', linestyle=':', label='Regression Line')
# # Equality line (y = x)
# plt.plot(adata1_mean, adata1_mean, color='gray', linestyle='--', label='Equality Line (y=x)')
plt.title(f'Pearson Correlation: {correlation_coefficient:.2f}')
plt.legend()
plt.show()

