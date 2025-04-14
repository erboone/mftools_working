from abc import ABC, abstractmethod
from pathlib import Path
from glob import glob
import warnings
import json

import pandas as pd
import scanpy as sc
import numpy as np
from scanpy import AnnData

from .segmentation import CellSegmentation
from .fileio import ImageDataset, MerfishAnalysis
from .fileio import _AbsExperimentSchema, MerscopeSchema, SmallMerscopeSchema, XeniumSchema 
from .scanpy_wrap import MerData
from .stats import stats

warnings.filterwarnings("ignore")

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

    # def __new__(cls,
    #         root:str,
    #         name:str,
    #         alt_paths:dict={},
    #         seg_kwargs:dict={},
    #         img_kwargs:dict={}
    #     ):
    #     pass

    def __init__(self, 
            root:str,
            name:str,
            reg:str=None,
            alt_paths:dict={},
            seg_kwargs:dict={},
            img_kwargs:dict={}
        ):

        # setting properties
        self.root = root
        self.name = name
        self.reg = reg
        self.files = self._schema_class(root, name, reg, schema_mod=alt_paths)
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
        # print("Setting segmentor")
        # TODO: Change this once implementation of CellSegmentation has been updated
        """Setter for the cell segmentation object for this MerfishExperiment.
        """
        if self._segmentator_class is None:
            self._segmentator_instance = object
        else:    
            self._segmentator_instance = self._segmentator_class(
                mask_folder = self.files['masks'],
                output = MerfishAnalysis(self.files['cellpose']),
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
        # print('Setting imageset')
        """Getter for the ImageDataset object for this MerfishExperiment.
        If object is not initialized, resolve relevant paths, initialize, and
        return, otherwise, just return.
        """
        if self._imageset_class is None:
            _imageset_instance = object
        else:
            self._imageset_instance = self._imageset_class(
                self.files['data'],
                **kwargs
            )

    imgs = property(imgs_get, imgs_set)

    def __load_dataframe(self, name: str, add_region: bool) -> pd.DataFrame:
        location = self.files['output']
        filename = Path(f'{location}/{name}')
        if filename.exists():
            return pd.read_csv(filename, index_col=0)
        # Check if this is a multi-region MERSCOPE experiment
        print()
        if list(glob(f'{location}/region_*')):
            region_dfs = []
            for region in list(glob(f'{location}/region_*')):
                region = Path(region)
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
    
    def assign_to_cells(self, barcodes, masks, drifts=None, transpose=False, flip_x=True, flip_y=True):
        for fov in tqdm(np.unique(barcodes["fov"]), desc="Assigning barcodes to cells"):
            group = barcodes.loc[barcodes["fov"] == fov]
            if drifts is not None:
                xdrift = drifts.loc[fov]["X drift"]
                ydrift = drifts.loc[fov]["Y drift"]
            else:
                xdrift, ydrift = 0, 0
            x = (group["x"] + xdrift).round() // config.get("scale")
            y = (group["y"] + ydrift).round() // config.get("scale")
            if flip_x:
                x = 2048 - x
            if flip_y:
                y = 2048 - y
            if transpose:
                x, y = y, x
            x = x.clip(upper=2047).astype(int)
            y = y.clip(upper=2047).astype(int)
            if len(masks[fov].shape) == 3:
                # TODO: Remove hard-coding of scale
                z = (group["z"].round() / 6.333333).astype(int)
                barcodes.loc[barcodes["fov"] == fov, "cell_id"] = masks[fov][z, x, y] + 10000 * fov
            else:
                barcodes.loc[barcodes["fov"] == fov, "cell_id"] = masks[fov][x, y] + 10000 * fov
        barcodes.loc[barcodes["cell_id"] % 10000 == 0, "cell_id"] = 0
        barcodes["cell_id"] = barcodes["cell_id"].astype(int)
        stats.set("Barcodes assigned to cells", len(barcodes[barcodes["cell_id"] != 0]))
        stats.set(
            "% barcodes assigned to cells",
            stats.get("Barcodes assigned to cells") / len(barcodes),
        )
        return barcodes
    
    @abstractmethod
    def create_scanpy_object(self):
        pass
        
class MerscopeExperiment(_AbsMFExperiment):

    _schema_class = MerscopeSchema
    _segmentator_class = CellSegmentation
    _imageset_class = ImageDataset

    def __init__(self, 
            root:str,
            name:str,
            reg:str=None,
            alt_paths:dict={},
            seg_kwargs:dict={},
            img_kwargs:dict={}
        ):
        super().__init__(root, name, reg, alt_paths=alt_paths, seg_kwargs=seg_kwargs,
                          img_kwargs=img_kwargs)
        try:
            with open(Path(self.files['settings']) / 'microscope_parameters.json', 'r') as file:
                self.settings = json.load(file)
        except FileNotFoundError:
                self.settings = {'image_dimensions': (204, 204)}

    def create_scanpy_object(self, name=None, positions=None, codebook=None, keep_empty_cells=True, return_adata=False):
        from .cellgene import adjust_spatial_coordinates
        merscope_ana = MerfishAnalysis(self.files['cellpose'])
        
        cellgene = merscope_ana.load_cell_by_gene_table()
        celldata = merscope_ana.load_cell_metadata()
        cellgene.index = cellgene.index.astype(str)
        celldata.index = celldata.index.astype(str)
        if keep_empty_cells:
            empty_cells = pd.DataFrame(
                [
                    pd.Series(data=0, index=cellgene.columns, name=cellid)
                    for cellid in celldata.index.difference(cellgene.index)
                ]
            )
            cellgene = pd.concat([cellgene, empty_cells])
        blank_cols = np.array(["notarget" in col or "blank" in col.lower() for col in cellgene])
        adata = sc.AnnData(cellgene.loc[:, ~blank_cols], dtype=np.int32)
        adata.obsm["X_blanks"] = cellgene.loc[:, blank_cols].to_numpy()
        adata.uns["blank_names"] = list(blank_cols)
        if "global_x" in celldata:
            adata.obsm["X_spatial"] = np.array(
                celldata[["global_x", "global_y"]].reindex(index=adata.obs.index)
            )
        elif "center_x" in celldata:
            adata.obsm["X_spatial"] = np.array(
                celldata[["center_x", "center_y"]].reindex(index=adata.obs.index)
            )
        if "fov_x" in celldata:
            adata.obsm["X_local"] = np.array(celldata[["fov_x", "fov_y"]].reindex(index=adata.obs.index))
        for column in celldata.columns:
            adata.obs[column] = celldata[column]
        adata.obs["fov"] = adata.obs["fov"].astype(str)
        adata.layers["counts"] = adata.X
        sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
        adata.obs["blank_counts"] = adata.obsm["X_blanks"].sum(axis=1)
        adata.obs["misid_rate"] = (adata.obs["blank_counts"] / len(adata.uns["blank_names"])) / (
            adata.obs["total_counts"] / len(adata.var_names)
        )
        adata.obs["counts_per_volume"] = adata.obs["total_counts"] / adata.obs["volume"]
        if codebook:
            adata.varm["codebook"] = codebook.set_index("name").loc[adata.var_names].filter(like="bit").to_numpy()
            for bit in range(adata.varm["codebook"].shape[1]):
                adata.obs[f"bit{bit+1}"] = adata[:, adata.varm["codebook"][:, bit] == 1].X.sum(axis=1)
        if positions:
            adata.uns["fov_positions"] = positions.to_numpy()
        if name:
            adata.uns["dataset_name"] = name
        else:
            adata.uns["dataset_name"] = merscope_ana.root.name

        adjust_spatial_coordinates(adata, flip_horizontal=True)

        if not return_adata:
            merdata = MerData(
                X=adata.X,
                obs=adata.obs,
                var=adata.var,
                uns=adata.uns,
                obsm=adata.obsm,
                varm=adata.varm,
                layers=adata.layers
                )
            return merdata
        else:
            return adata
    
class SmallMerscopeExperiment(_AbsMFExperiment):

    _schema_class = SmallMerscopeSchema
    _segmentator_class = None
    _imageset_class = None

    

    def __init__(self,
            root:str,
            name:str,
            alt_paths:dict={},
            seg_kwargs:dict={},
            img_kwargs:dict={}
        ):
        
        super().__init__(root, name, alt_paths=alt_paths, 
                        seg_kwargs=seg_kwargs,
                        img_kwargs=img_kwargs)
        
    def create_scanpy_object(self):
        # Load data
        data_p = self.files['data']
        cbgt = sc.read_csv(f'{data_p}/cell_by_gene.csv', first_column_names=True)
        cells = pd.read_csv(f'{data_p}/cell_metadata.csv')
        
        # Maniplating into raw counts data
        cbgt.obs = cells
        drop_cols = ['EntityID', 'anisotropy', 'transcript_count', 'perimeter_area_ratio', 'solidity', 'SPP1_raw', 'SPP1_high_pass', 'TAGLN_raw', 'TAGLN_high_pass', 'SFTPB_raw', 'SFTPB_high_pass', 'Cellbound2_raw', 'Cellbound2_high_pass', 'Cellbound3_raw', 'Cellbound3_high_pass', 'DAPI_raw', 'DAPI_high_pass', 'UMOD_raw', 'UMOD_high_pass', 'PolyT_raw', 'PolyT_high_pass', 'Cellbound1_raw', 'Cellbound1_high_pass', 'ACTA2_raw', 'ACTA2_high_pass', 'OLFM4_raw', 'OLFM4_high_pass']
        cbgt.obs_names = cbgt.obs['EntityID'].rename('cell').astype('string')
        cbgt.obs.drop(drop_cols, axis=1, inplace=True, errors='ignore')
        cbgt = cbgt[:, ~cbgt.var_names.str.contains('Blank-')]
        cbgt.obsm['X_spatial'] = np.array(cbgt.obs[['center_x','center_y']])
        adata = cbgt
        return adata

    def seg_set(self, **kwargs):
        # print("Setting segmentor")
        # TODO: Change this once implementation of CellSegmentation has been updated
        """Setter for the cell segmentation object for this MerfishExperiment.
        """
        self._segmentator_instance = object
    
    def img_set(self, **kwargs):
        # print("Setting segmentor")
        # TODO: Change this once implementation of CellSegmentation has been updated
        """Setter for the cell segmentation object for this MerfishExperiment.
        """
        self._segmentator_instance = object


class XeniumExperiment(_AbsMFExperiment):
    # This class mostly exists to access files
    # The existing image access class is heavily biased to MERSCOPE
    # Adding Xenium will be a lot of work.
    _schema_class = XeniumSchema
    _segmentator_class = object

    def __init__(self, 
            root:str,
            name:str,
            alt_paths:dict={},
            seg_kwargs:dict={},
            img_kwargs:dict={'data_organization':'/mnt/merfish15/MERSCOPE/data/202407221121_20240722M176BICANRen22_VMSC10002/dataorganization.csv'}
        ):
        super().__init__(root, name, alt_paths=alt_paths, seg_kwargs=seg_kwargs,
                          img_kwargs=img_kwargs)

    def seg_set(self, **kwargs):
        # print("Setting segmentor")
        # TODO: Change this once implementation of CellSegmentation has been updated
        """Setter for the cell segmentation object for this MerfishExperiment.
        """
        self._segmentator_instance = object

    def create_scanpy_object(self):
        xenium_ad = sc.read_10x_h5(f"{self.files['reseg']}/cell_feature_matrix.h5")
        df = pd.read_csv(f"{self.files['reseg']}/cells.csv.gz", index_col=0)
        xenium_ad.obsm["X_spatial"] = df.loc[xenium_ad.obs_names][["x_centroid", "y_centroid"]].to_numpy()
        return xenium_ad

    # def seg_get(self, **kwargs):
    #     """Getter for the CellSegmentation object for this MerfishExperiment.
    #     If object is not initialized, resolve relevant paths, initialize, and
    #     return, otherwise, just return.
    #     """
    #     # This check may not be needed anymore
    #     return None 

    # def seg_set(self, **kwargs):
    #     print("Setting segmentor")
    #     # TODO: Change this once implementation of CellSegmentation has been updated
    #     """Setter for the cell segmentation object for this MerfishExperiment.
    #     """
    #     self._segmentator_instance = None
    # seg = property(seg_get, seg_set)

if __name__ == "__main__":
    print("this works")