from itertools import product
from re import match

import scanpy as sc
from scanpy import AnnData

 #=-=#=-=#=-=#=-=#=-=#=-=#=- Constants -=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#

QUALITY_CONTROL_COLUMNS_RE =  [
    'total_*_by_*', 
    'total_*', 
    'pct_*_in_top_*_*', 
    'total_*_*', 
    'pct_*_*', 
    'total_*', 
    'n_genes_by_*', 
    'mean_*', 
    'n_cells_by_*', 
    'pct_dropout_by_*'
]

REQ_LOG = [
    'step_name',
    'step_hash',
    'function',
    'commit_id',
    'normalization'
]

REQ_COL = [
    '_BATCH', '_CELLTYPE',
    'fov_y', 'fov_x', 'fov', # fov info
    'global_x', 'global_y', 
    'volume' # for QC filtering later
]

#=-=#=-=#=-=#=-=#=-=#=-=#=- Error msg. -=#=-=#=-=#=-=#=-=#=-=#=-=#=-=#
UNPACK_SETTER_ERR = "Pass an iterable with two items; only the first two will be used"
MISSING_LOG_ERR = """You are missing data required to keep track of changes file"""
MISSING_COL_ERR = """You are missing columns required to identify rows in the future"""
FOUND_QC_ERR = "Multiple columns in MerData object match the output of method 'calculate_qc_metrics'. Please save this elsewhere"

class MerData(AnnData):
    """
    This is a wrapper for the standard scanpy Anndata object to implement
    some features that will make standardization in our pipeline easier. 

    Args:
        AnnData (_type_): _description_
    """

    # Cannonized keys, type format to avoid collision
    BATCH_KEY = '_BATCH'
    CTYPE_KEY = '_CELLTYPE'
    LOG_KEY = 'current_log'
    HISTORY_KEY = 'past_logs'

    def __init__(self, adata:AnnData=None, *args, **kwargs):
        if adata is not None:
            kwargs = {
                'X': adata.X ,
                'obs': adata.obs,
                'var': adata.var,
                'uns': adata.uns,
                'obsm': adata.obsm,
                'varm': adata.varm,
                'layers': adata.layers
            }
    
        super().__init__(*args, **kwargs)
        self.uns[self.HISTORY_KEY] = [] # Never change the order of this list
        self.uns[self.LOG_KEY]  = {}
    
    @property
    # Minor alias history  
    def history(self):
        return self.uns[self.HISTORY_KEY]

    @property
    # Minor alias for log to make it less unwieldy 
    def log(self):
        return self.uns[self.LOG_KEY]
    
    @log.setter
    def log(self, x):
        try:
            a, b = x
        except:
            raise ValueError(UNPACK_SETTER_ERR)
        self.uns[self.LOG_KEY][a] = b
    
    def _save_log(self) -> None:
        try:
            self.history.append(
                '/n'.join([f'{l}: {self.log[l]}' for l in REQ_LOG])
            )

        except KeyError as e:
            raise KeyError(MISSING_LOG_ERR) from e


    def check_requirements(self):
        
        # check metadata
        try:
            for l in REQ_LOG:
                self.log[l]

        except KeyError as e:
            raise KeyError(MISSING_LOG_ERR) from e
        
        # check columns
        try:
            for c in REQ_COL:
                temp = self.obs[c]
                # TODO: include some checks for the data in the columns?
                # if len(set(temp)) < 1:

        except KeyError as e:
            raise KeyError(MISSING_COL_ERR) from e
        
        # TODO: This doesn't work
        i = 0
        for re, s in product(QUALITY_CONTROL_COLUMNS_RE, self.obs.columns):
            if match(re, s):
                i += 1
            if i > 4:
                raise RuntimeError(FOUND_QC_ERR)
        
        return True


    # Check validity, then write. DO NOT alter data.
    def safe_write(self, *args, **kwargs):
        """Required information
        Log keys:
            'step_name',
            'step_hash',
            'function',
            'commit_id',
            'normalization'

        Obs columns:
            '_BATCH', 
            'fov_y', 'fov_x', 'fov', # fov info
            'global_x', 'global_y', 
            'volume'                 # for QC filtering later
            **NO calculate_qc_metrics in obs**
        Raises:
            KeyError: _description_
            KeyError: _description_
            RuntimeError: _description_
        """
        if self.check_requirements():
            self._save_log()
            sc.write(adata=self, *args, **kwargs)
    

if __name__ == "__main__":
    pass    
    


