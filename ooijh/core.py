from datetime import datetime
import fsspec
import multiprocessing
import numpy as np
import os
import re
import xarray as xr

from ooijh.drops import DROP_METHODS, DROP_STREAMS, DROP_SYS_VARS, drop_qc_test_vars

_USER_DIR = os.path.expanduser('~')


class KDATA():
    
    BASE_DIR = f"{_USER_DIR}/ooi/kdata/"

    def __init__(self, site: str, node: str, instrument: str, stream: str,
                 begin_datetime: datetime, end_datetime: datetime, 
                 process_qartod: bool = True , nan_flags: list = [4,9], drop_qartod: bool = True):
        
        self.__local = fsspec.filesystem('file')
        self.bdt = begin_datetime
        self.edt = end_datetime
        self._usnims = (site, node, instrument, stream)
        self.files = self.find_files()
        self.process_qartod = process_qartod
        self.drop_qartod = drop_qartod        
        self.nan_flags = nan_flags
        
        
    def list_all_datasets(self, data_only: bool = True):
        folders = self.__local.glob(self.BASE_DIR + '*')
        rds = sorted([i.replace(self.BASE_DIR,'') for i in folders])
        if data_only is True:
            for _method in DROP_METHODS:
                rds = [rd for rd in rds if _method not in rd]
            for _stream in DROP_STREAMS:
                rds = [rd for rd in rds if _stream not in rd]
        return rds
    
    
    def find_datasets(self, site, node, instrument, stream, data_only: bool = True):
        rds = self.list_all_datasets(data_only = data_only)
        rds = [rd for rd in rds if site in rd and node in rd and instrument in rd and stream in rd]
        return rds
    
    
    def list_dataset_files(self, rd, filetype = '.nc'):
        files = sorted(self.__local.glob(self.BASE_DIR + rd + f'/*{filetype}'))
        dtrfiles =[]
        for file in files:
            if len(re.findall(rd,file)) != 2: #Drop ancillary datasets.
                continue

            file_begin, file_end = [datetime.strptime(dt,'%Y%m%dT%H%M%S') for dt in re.findall("(\d{8}T\d{6})",file)]
            if (file_begin <= self.bdt <= file_end or file_begin <= self.edt <= file_end or 
                self.bdt <= file_begin <= self.edt or self.bdt <= file_end <= self.edt):            
                dtrfiles.append(file)        
        
        
        return dtrfiles
    
    def find_files(self):
        rds = self.find_datasets(*self._usnims)
        files = [] 
        for rd in rds:
            rdfiles = self.list_dataset_files(rd)
            files = files + rdfiles
        return files
            
    
    def import_file(self, filepath):
        ds = xr.open_dataset(filepath)
        ds = ds.swap_dims({'obs':'time'})
        ds = ds.drop_vars(['obs'], errors = 'ignore')
        ds = ds.drop_vars(DROP_SYS_VARS, errors = 'ignore') # Remove confusing OOI backend variables.
        ds = drop_qc_test_vars(ds) # Remove custom qc test variables. Deprecated???
        
        if 'lat' not in ds.data_vars and 'latitude' not in ds.data_vars: #If lat is not a variable, make lat and lon variables from attribute data.
            ds['latitude'] = ds.lat
            ds['longitude'] = ds.lon

        ds = ds.sortby('time')
        return ds
    
    def open_datasets(self):
        with multiprocessing.Pool(len(self.files)) as pool:
            ds_list = pool.map(self.import_file, self.files)
        return ds_list
    
    def combine_data(self, ds_list):
        try:
            ds = xr.combine_by_coords(ds_list, combine_attrs = 'drop')
            ds = ds.drop_duplicates(list(ds.dims))
        except:
            ds = xr.concat(ds_list, dim = 'time', combine_attrs = 'drop')
        ds = ds.sortby('time')
        ds = ds.drop_duplicates(list(ds.dims))
        ds = ds.sel(time = slice(self.bdt.strftime('%Y-%m-%d %H:%M:%S'),self.edt.strftime('%Y-%m-%d %H:%M:%S')))
        return ds
    
    
    def nan_by_qartod(self, ds, nan_flags = [4,9], variables = 'all'):
        if variables == 'all':
            qecs = [v for v in list(ds.data_vars) if 'qartod_executed' in v]
        else:
            qecs = [v + '_qartod_executed' for v in variables]
        for qec in qecs:
            qdc = qec.split('_qartod_executed')[0]
            qrc = qdc + '_qartod_results'   
            ds[qdc] = ds[qdc].where(ds[qec] != 1 & ~ds[qrc].isin(nan_flags), np.nan)
            # if handle_confusion is True:
            #     ds[qdc] = ds[qdc].where(~ds[qec].isin([3, 4]), np.nan)
        return ds