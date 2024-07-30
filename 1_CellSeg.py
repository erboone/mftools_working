#!/usr/bin/env python
# coding: utf-8

import os, sys
from configparser import ConfigParser
from pathlib import Path
from itertools import chain

import numpy as np
from skimage.measure import regionprops # This has been moved to plotting
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.axes import Axes 
from cellpose import models

from mftools.fileio import ImageDataset, MerfishAnalysis
from mftools.plotting import fov_show
from mftools.segmentation import CellSegmentation
from mftools.barcodes import assign_to_cells, link_cell_ids, create_cell_by_gene_table
from mftools.cellgene import create_scanpy_object
from mftools.mfexperiment import MerscopeExperiment

# ----- Load Config ----- #
# TODO: change this to work with config files generated from merbot
config:ConfigParser = ConfigParser()
config.read('config.ini')
ioconf = config["IO Options"]
csconf = config["1_CellSeg"]

MERSCOPE_DIR = ioconf['msdir']
EXPERIMENT_NAME = ioconf['experiment']

# PATH PARAMS 
# CHANGE ONLY IF NEEDED: these should be standardized in the future and creation of these dirs automated from
# globalized config file
MER_RAWDATA_DIR = "data"
MER_OUTPUT_DIR = "output"
CELLPOSE_DIR = f"./cellpose_{EXPERIMENT_NAME}" # TODO: for testing, change later
MASKS_DIR = "masks" 

# Script params
ZSLICE = csconf['zslice']
CHANNEL = csconf['channel']
MODEL= ""
TEST_FOVS = [150, 300, 500]

# Path assembly TODO: remove this when we can trust MERBOT to generate these for
# image_dataset_path = f"{MERSCOPE_DIR}/{MER_RAWDATA_DIR}/{EXPERIMENT_NAME}/" # Used to load imageset
# expiriment_out_path = f"{MERSCOPE_DIR}/{MER_OUTPUT_DIR}/{EXPERIMENT_NAME}/" # Used to find barcodes 
# cellpose_out_path = f"./{CELLPOSE_DIR}/" # used to save final cellpose output
# masks_out_path = f"{cellpose_out_path}{MASKS_DIR}/" # used to save masks

# print(f"Looking for images in {image_dataset_path}")
# print(f"Looking for barcodes in {expiriment_out_path}")
# print(f"Writing segmenatation output {cellpose_out_path}; {masks_out_path}")


# # Check directory structure
# needed_dirs = [image_dataset_path, expiriment_out_path, cellpose_out_path, masks_out_path]
# print(needed_dirs)
# for path in needed_dirs:
#     os.makedirs(path, exist_ok=True)

#     if os.path.exists(path):
#         if not os.access(path, os.W_OK) or not os.access(path, os.R_OK):
#             raise RuntimeError(f"You do not have read/write permissions in directory \"{path}\"")
#     else:
#         raise RuntimeError(f"Attempted to create \"{path}\" but failed. Check that you have permission.")
    

# Use file path like format ("/mnt/merfish12/MERSCOPE/merfish_raw_data/202401261424_20240126M134UWA7648CX22PuS6_VMSC10102
experiment = MerscopeExperiment(MERSCOPE_DIR, EXPERIMENT_NAME)
e = experiment

fig:Figure; axs:list[list[Axes]]
# fig, axs = plt.subplots(2, len(TEST_FOVS))

# TODO: come back to this; implement plotting and saving a couple of test segmentations
for i, fov in enumerate(TEST_FOVS):
    fig, ax = plt.subplots(1, 1)
    fov_show(seg=e.seg, imgs=e.imgs, fov=fov, show=False, ax=ax)
    fig.savefig(f"quality_control/test_{fov}")


exit()
# either load or create meatadata, then save if not already
metadata = e.seg.metadata
e.save_cell_metadata(metadata)


# load barcodes, assign to cells, save to disk

barcodes = e.load_barcode_table()
assign_to_cells(barcodes, e.seg)
link_cell_ids(barcodes, e.seg.linked_cells)
e.save_barcode_table(barcodes)

barcodes = e.load_barcode_table()
cbgtab = create_cell_by_gene_table(barcodes)
cbgtab.index = cbgtab.index.astype(int)
e.save_cell_by_gene_table(cbgtab)
# TODO: Create and write the scanpy object somewhere.