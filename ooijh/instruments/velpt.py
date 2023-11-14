from datetime import datetime
import multiprocessing
import numpy as np
import xarray as xr

from ooijh.core import KDATA
from ooijh.drops import DROP_VELPT_VARS, drop_qartod_test_vars


class VELPT(KDATA):
    def __init__(self, site, node, instrument = 'VELPT', stream = 'velpt', 
                 begin_datetime = datetime(2014,1,1), end_datetime = datetime(2040,12,31,23,59,59),
                process_qartod_vars = True, drop_qartod_test_vars = True):
        super().__init__(site.upper(), node.upper(), instrument.upper(), stream.lower(), begin_datetime, end_datetime)
        
    
    def process(self):
        # COMBINE DATA GOES HERE
        ds_list = self.open_datasets()
        ds_list = [ds.drop_vars(DROP_VELPT_VARS, errors = 'ignore') for ds in ds_list]
        if self.process_qartod is True:
            ds_list = [self.nan_by_qartod(ds, self.nan_flags) for ds in ds_list]
        if self.drop_qartod is True:
            with multiprocessing.Pool(len(ds_list)) as pool:
                ds_list = pool.map(drop_qartod_test_vars, ds_list)
        ds = self.combine_data(ds_list)
        
        
        
        # ADDITIONAL CLEANUP GOES HERE
  
        
        # ADDITIONAL PROCESSING GOES HERE

        
        ds = ds[sorted(ds.data_vars)] #Sort variables alphabetically because it is less obtrusive.
        return ds
        
        
            