import unittest
from mftools import mfexperiment
from mftools.segmentation import CellSegmentation
import scanpy as sc
import numpy as np
from pathlib import Path # used from mfexperiment (helps me with the paths)

class MSEperiment_Test(unittest.TestCase): 
    

    def __init__(self, cls, res, *args, **kwargs):
        super(MSEperiment_Test, self).__init__(*args, **kwargs)
        self.e = cls(res.rootdir, res.name, res.region)


    def test__create_scanpy_object(self):
        adata:sc.AnnData = self.e.create_scanpy_object()
        self.assertTrue(isinstance(adata, mfexperiment.MerData))
        self.assertTrue(adata.shape[0] > 0 and adata.shape[1] > 0)
        self.assertFalse(adata.obs.isnull().values.any())
        self.assertFalse(adata.var.isnull().values.any())
        self.assertTrue('X_spatial' in adata.obsm)
        self.assertTrue(len(adata.obsm['X_spatial']) == adata.shape[0])
        self.assertTrue(all(~adata.obs['fov'].isna()))
        self.assertTrue(all(~adata.obs['global_x'].isna()))
        self.assertTrue(all(~adata.obs['global_y'].isna()))
        self.assertTrue(all(~adata.obs['global_volume'].isna()))
        self.assertTrue(np.all(np.mod(adata.X, 1) == 0))


    def test__images(self):
        pass

        # Things to check:
        # check that X is counts
        # check that X has no NaN
        # check status of obs
        
    def test__fileschema(self):
        pass

    def test__segmentation(self):
        pass

    def test__images(self):
        pass

# Creating a class for unit tests for a regular merscope experiment
# Inherits the other tests from MSEperiment_Test
class Test_Regular_MSExperiment(MSEperiment_Test): 
    
    # manually creating a directory to use for the test
    def setUp(self):
        # note -> change the test_regular_exp to an ACTUAL experiment (maybe?)
        self.test_dir = Path("test_regular_exp") 
        self.e.savepath = self.test_dir # making sure I use this path in the test

        # TODO: create a dataframe (or use our own to load in to check)

    # testing necessary functions -> cell_by_gene (should work regardless)
    def test_cell_by_gene(self):
        self.e.save_cell_by_gene_table() # add the dataframe we're using in here
        loaded = self.e.load_cell_by_gene_table() # loading the dataframe
        self.assertTrue((self.savepath / "cell_by_gene.csv").exists()) # making sure this exists
        # TODO: I need to make sure the contents of the load are correct too


    # test an interesting function (adjust_spatial_coordinates)
    # note: only in regular merscope experiment (important!)
    def test_adjust_spatial_coordinates(self):
        adata:sc.AnnData = self.e.adjust_spatial_coordinates() # wrong -> fix

    # checking for the metadata - could be helpful
    def test_cell_metadata(self):
        pass


if __name__ == "__main__":
    var = 'string'

    t = newtest()
    
    print(t.foo)
    print(t.bar)
    msexp_classe


