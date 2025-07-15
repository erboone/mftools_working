import unittest
from mftools import mfexperiment
from mftools.segmentation import CellSegmentation
import scanpy as sc
import numpy as np

class MSEperiment_Test(unittest.TestCase): 
    

    def __init__(self, cls, res, *args, **kwargs):
        super(MSEperiment_Test, self).__init__(*args, **kwargs)
        self.e = cls(res.rootdir, res.name, res.region)


    def __test__create_scanpy_object(self):
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



if __name__ == "__main__":
    test = MSEperiment_Test()
    unittest.main()