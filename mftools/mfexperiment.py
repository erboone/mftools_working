from abc import ABC, abstractmethod
from pathlib import Path
import pandas as pd

from .segmentation import CellSegmentation
from .fileio import ImageDataset, MerfishAnalysis, MerscopeSchema, _AbsExperimentSchema

class _AbsMFExperiment(ABC):
    """
    The goal of this class is to serve as an abstract wrapper for a merfish experiment class. 
    It should contain all the information needed to automatically resolve and check paths to the data. 
    I would like some sort of general solution that allows us to quickly update the pathing schema when we
    want to apply the abstract class to a new experiment format (i.e. different schema for MERSCOPE and Xenium etc.)
    """

    _schema_class = _AbsExperimentSchema
    _segmentator_class = CellSegmentation
    _imageset_class = ImageDataset

    def __init__(self, 
            root:str,
            name:str,
            alt_paths:dict={},
            seg_kwargs:dict={},
            img_kwargs:dict={}
        ):

        # setting properties
        self.root = root
        self.name = name
        self.files = self._schema_class(root, name, schema_mod=alt_paths)
        self.savepath = Path(self.files['cellpose'])
        # Properties set using @property.setter; see below
        # def namespace
        self._segmentator_instance:CellSegmentation = self._segmentator_class
        self._imageset_instance:ImageDataset = None
        # These set the objects based of info in self.files; supplimented
        # by kwargs.
        # Note: For now, imgs must be defined before seg.
        self.imgs_set(**img_kwargs)
        self.seg_set(**seg_kwargs)        

    # TODO: Decide if this is need.__s
    # @property
    # def files(self):
    #     """_summary_
    #     Object built to retrieve files using a known mfexperiment schema
    #     """
    #     pass
    
    # @files.setter
    # def files(self):
    #     pass

    def seg_get(self, **kwargs):
        """Getter for the CellSegmentation object for this MerfishExperiment.
        If object is not initialized, resolve relevant paths, initialize, and
        return, otherwise, just return.
        """
        # This check may not be needed anymore
        if self._segmentator_instance is None:
            self.seg.__setattr__()
        return self._segmentator_instance

    def seg_set(self, **kwargs):
        print("Setting segmentor")
        # TODO: Change this once implementation of CellSegmentation has been updated
        """Setter for the cell segmentation object for this MerfishExperiment.
        """
        self._segmentator_instance = self._segmentator_class(
            mask_folder = self.files['masks'],
            output = MerfishAnalysis(self.files['output']),
            imagedata=self.imgs,
            **kwargs
        )

    seg = property(seg_get, seg_set)

    def imgs_get(self):
        """Getter for the ImageDataset object for this MerfishExperiment.
        If object is not initialized, resolve relevant paths, initialize, and
        return, otherwise, just return.
        """
        return self._imageset_instance


    # @imgs.setter
    def imgs_set(self, **kwargs):
        print('Setting imageset')
        """Getter for the ImageDataset object for this MerfishExperiment.
        If object is not initialized, resolve relevant paths, initialize, and
        return, otherwise, just return.
        """
        self._imageset_instance = self._imageset_class(
            self.files['data'],
            **kwargs
        )

    imgs = property(imgs_get, imgs_set)

    def __load_dataframe(self, name: str, add_region: bool) -> pd.DataFrame:
        filename = self.savepath / name
        if filename.exists():
            return pd.read_csv(filename, index_col=0)
        # Check if this is a multi-region MERSCOPE experiment
        if list(self.root.glob("region_*")):
            region_dfs = []
            for region in list(self.root.glob("region_*")):
                num = str(region).rsplit("_", maxsplit=1)[-1]
                dataframe = pd.read_csv(region / name, index_col=0)
                if add_region:
                    dataframe["region"] = num
                region_dfs.append(dataframe)
            dataframe = pd.concat(region_dfs)
            dataframe.to_csv(filename)
            return dataframe
        raise FileNotFoundError(filename)

    def save_cell_metadata(self, celldata: pd.DataFrame) -> None:
        celldata.to_csv(self.savepath / "cell_metadata.csv")

    def load_cell_metadata(self) -> pd.DataFrame:
        return self.__load_dataframe("cell_metadata.csv", add_region=True)

    def has_cell_metadata(self) -> bool:
        return Path(self.savepath, "cell_metadata.csv").exists()

    def save_linked_cells(self, links) -> None:
        with open(self.savepath / "linked_cells.txt", "w", encoding="utf8") as f:
            for link in links:
                print(repr(link), file=f)

    def load_linked_cells(self):
        links = []
        with open(self.savepath / "linked_cells.txt", encoding="utf8") as f:
            for line in f:
                links.append(eval(line))
        return links

    def save_barcode_table(self, barcodes, dask=False) -> None:
        if dask:
            barcodes.to_csv(self.savepath / "detected_transcripts")
        else:
            barcodes.to_csv(self.savepath / "detected_transcripts.csv")

    def load_barcode_table(self) -> pd.DataFrame:
        return self.__load_dataframe("detected_transcripts.csv", add_region=True)

    def save_cell_by_gene_table(self, cellbygene) -> None:
        cellbygene.to_csv(self.savepath / "cell_by_gene.csv")

    def load_cell_by_gene_table(self) -> pd.DataFrame:
        return self.__load_dataframe("cell_by_gene.csv", add_region=False)
    
    @abstractmethod
    def create_scanpy_object():
        pass
        

class MerscopeExperiment(_AbsMFExperiment):

    _schema_class = MerscopeSchema

    def __init__(self, 
            root:str,
            name:str,
            alt_paths:dict={},
            seg_kwargs:dict={},
            img_kwargs:dict={}
        ):
        super().__init__(root, name, alt_paths=alt_paths, seg_kwargs=seg_kwargs,
                          img_kwargs=img_kwargs)

    def create_scanpy_object():
        pass


if __name__ == "__main__":
    print("this works")