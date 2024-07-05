from abc import ABC
from segmentation import CellSegmentation
from fileio import ImageDataset, MerfishAnalysis


class _AbsMFExperiment(ABC):
    """
    The goal of this class is to serve as an abstract wrapper for a merfish experiment class. 
    It should contain all the information needed to automatically resolve and check paths to the data. 
    I would like some sort of general solution that allows us to quickly update the pathing schema when we
    want to apply the abstract class to a new experiment format (i.e. different schema for MERSCOPE and Xenium etc.)
    """

    # TODO: improve this system
    MER_RAWDATA_DIR = "data"
    MER_OUTPUT_DIR = "output"
    CELLPOSE_DIR = "cellpose"
    MASKS_DIR = "masks" 


    def __init__(self, location:str, name:str):

        image_dataset_path = f"{location}/{self.MER_RAWDATA_DIR}/{name}/" # Used to load imageset
        expiriment_out_path = f"{location}/{self.MER_OUTPUT_DIR}/{name}/" # Used to find barcodes 
        cellpose_out_path = f"{expiriment_out_path}/{self.CELLPOSE_DIR}/" # used to save final cellpose output
        masks_out_path = f"{cellpose_out_path}{self.MASKS_DIR}/" # used to save masks

        self.original = MerfishAnalysis(expiriment_out_path)
        self.output = MerfishAnalysis(masks_out_path)
        self.imgset = ImageDataset(image_dataset_path)
        self.seg = CellSegmentation(imagedata=self.imgset, channel='dapi', zslice=3)

    @property
    def seg(self):
        """Getter for the CellSegmentation object for this MerfishExperiment.
        If object is not initialized, resolve relevant paths, initialize, and
        return, otherwise, just return.
        """

    @property.setter
    def cellsegmentation(self):
        # TODO: decide if this is needed
        """Setter for the cell segmentation object for this MerfishExperiment.
        """

    @property
    def imgs(self):
        """Getter for the ImageDataset object for this MerfishExperiment.
        If object is not initialized, resolve relevant paths, initialize, and
        return, otherwise, just return.
        """
    @property
    def imgs(self):
        """Getter for the ImageDataset object for this MerfishExperiment.
        If object is not initialized, resolve relevant paths, initialize, and
        return, otherwise, just return.
        """
        




if __name__ == "__main__":
    print("this works")