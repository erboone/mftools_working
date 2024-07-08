#!/usr/bin/env python
# coding: utf-8

import os, sys
from configparser import ConfigParser
from pathlib import Path

import numpy as np
from skimage.measure import regionprops # This has been moved to plotting
import matplotlib.pyplot as plt
from matplotlib.axes import Axes 
from cellpose import models

from mftools.fileio import ImageDataset, MerfishAnalysis
from mftools.plotting import fov_show
from mftools.segmentation import CellSegmentation
from mftools.barcodes import assign_to_cells, link_cell_ids, create_cell_by_gene_table
from mftools.cellgene import create_scanpy_object

# ----- Load Config ----- #
# TODO: change this to work with config files generated from merbot
config:ConfigParser = ConfigParser()
config.read('config.ini')
ioconf = config["IO Options"]
csconf = config["2 Label Transfer"] # TODO:

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
TEST_FOV = 150

# Path assembly TODO: remove this when we can trust MERBOT to generate these for
image_dataset_path = f"{MERSCOPE_DIR}/{MER_RAWDATA_DIR}/{EXPERIMENT_NAME}/" # Used to load imageset
expiriment_out_path = f"{MERSCOPE_DIR}/{MER_OUTPUT_DIR}/{EXPERIMENT_NAME}/" # Used to find barcodes 
cellpose_out_path = f"./{CELLPOSE_DIR}/" # used to save final cellpose output
masks_out_path = f"{cellpose_out_path}{MASKS_DIR}/" # used to save masks

print(f"Looking for images in {image_dataset_path}")
print(f"Looking for barcodes in {expiriment_out_path}")
print(f"Writing segmenatation output {cellpose_out_path}; {masks_out_path}")


# Check directory structure
needed_dirs = [image_dataset_path, expiriment_out_path, cellpose_out_path, masks_out_path]
print(needed_dirs)
for path in needed_dirs:
    os.makedirs(path, exist_ok=True)

    if os.path.exists(path):
        if not os.access(path, os.W_OK) or not os.access(path, os.R_OK):
            raise RuntimeError(f"You do not have read/write permissions in directory \"{path}\"")
    else:
        raise RuntimeError(f"Attempted to create \"{path}\" but failed. Check that you have permission.")
    

# Use file path like format ("/mnt/merfish12/MERSCOPE/merfish_raw_data/202401261424_20240126M134UWA7648CX22PuS6_VMSC10102
output = MerfishAnalysis(cellpose_out_path)
original = MerfishAnalysis(expiriment_out_path)
imageset = ImageDataset(image_dataset_path)
seg = CellSegmentation(imagedata=imageset, output=output, channel=CHANNEL, zslice=ZSLICE)

fov_show(seg, imageset, TEST_FOV)

# either load or create meatadata, then save if not already
metadata = seg.metadata
output.save_cell_metadata(metadata)


# load barcodes, assign to cells, save to disk

barcodes = original.load_barcode_table()
assign_to_cells(barcodes, seg)
link_cell_ids(barcodes, seg.linked_cells)
output.save_barcode_table(barcodes)

barcodes = output.load_barcode_table()
cbgtab = create_cell_by_gene_table(barcodes)
cbgtab.index = cbgtab.index.astype(int)
output.save_cell_by_gene_table(cbgtab)
# TODO: Create and write the scanpy object somewhere.