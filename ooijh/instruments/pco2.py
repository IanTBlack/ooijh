from datetime import datetime
import gsw
import multiprocessing
import numpy as np
import xarray as xr

from ooijh.core import KDATA
from ooijh.drops import DROP_PCO2_VARS, drop_qartod_test_vars


class PCO2(KDATA):
    def __init__(self, site, node, instrument = 'PCO2', stream = 'pco2', 
                 begin_datetime = datetime(2014,1,1), end_datetime = datetime(2040,12,31,23,59,59),
                process_qartod_vars = True, drop_qartod_test_vars = True):
        super().__init__(site.upper(), node.upper(), instrument.upper(), stream.lower(), begin_datetime, end_datetime)
        
    
    def process(self):
        # COMBINE DATA GOES HERE
        ds_list = self.open_datasets()
        ds_list = [ds.drop_vars(DROP_PCO2_VARS, errors = 'ignore') for ds in ds_list]
        if self.process_qartod is True:
            ds_list = [self.nan_by_qartod(ds, self.nan_flags) for ds in ds_list]
        if self.drop_qartod is True:
            with multiprocessing.Pool(len(ds_list)) as pool:
                ds_list = pool.map(drop_qartod_test_vars, ds_list)
        ds = self.combine_data(ds_list)
        
        
        
        # ADDITIONAL CLEANUP GOES HERE

        
        
        # ADDITIONAL PROCESSING GOES HERE

        
        
        return ds
        
        
            