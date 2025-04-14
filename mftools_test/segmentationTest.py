import unittest
from mftools import mfexperiment
from mftools.segmentation import CellSegmentation
import scanpy as sc
import numpy as np

TEST_FOV = 150

class Segmentation_Test(unittest.TestCase): 
    

    def __init__(self, cls, res, *args, **kwargs):
        super(Segmentation_Test, self).__init__(*args, **kwargs)
        self.e = cls(res.rootdir, res.name, res.region,
                     alt_paths={
                         'cellpose':'/home/erboone/mftools/mftools_test/_test',
                         'masks':'/home/erboone/mftools/mftools_test/_test/masks'
                        }
                    )

    def test__segmentation(self):
        seg:CellSegmentation = self.e.seg
        mask = seg.segment_fov(TEST_FOV)
        self.assertEqual(mask.ndim, 3)
        self.assertEqual(mask.shape[1], mask.shape[2])
        self.assertTrue(mask.shape[0], mask.shape[1])



if __name__ == "__main__":
    var = 'string'

    t = newtest()
    
    print(t.foo)
    print(t.bar)
    msexp_classes

