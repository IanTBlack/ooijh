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
                 process_qartod: bool = True , nan_flags: list = [4,9], drop_qartod: bool = True) -> None:

        """
        A CAF class for finding and import OOI data from the kdata directory .

        :param site: An 8-character OOI designator for a site.
            It is recommended to use the full site designator.
        :param node: A 5-character designator for a node. Partial strings acceptable.
        :param instrument: A 12-character designator for an instrument. Partial strings acceptable.
        :param stream: A variable character designator for an instrument stream. Partial strings acceptable.
        :param begin_datetime: Beginning of date range as a datetime object.
        :param end_datetime: Ending of date range as a datetime object.
        :param process_qartod: If True, variables with associated qartod tests will be NaNed if they fall within the
            suppled nan_flags.
        :param nan_flags: Determines which data to NaN.
        :param drop_qartod: Drop the qartod test variables to clean up space.
        """

        self.__local = fsspec.filesystem('file')
        self.bdt = begin_datetime
        self.edt = end_datetime
        self._usnims = (site, node, instrument, stream)
        self.files = self.find_files()
        self.process_qartod = process_qartod
        self.nan_flags = nan_flags
        self.drop_qartod = drop_qartod        


    def list_all_datasets(self, data_only: bool = True) -> list:
        """
        List all datasets in the kdata directory.

        :param data_only: If True, excludes engineering/calibration datasets.
        :return: A list of reference designators.
        """

        folders = self.__local.glob(self.BASE_DIR + '*')
        rds = sorted([i.replace(self.BASE_DIR,'') for i in folders])
        if data_only is True:
            for _method in DROP_METHODS:
                rds = [rd for rd in rds if _method not in rd]
            for _stream in DROP_STREAMS:
                rds = [rd for rd in rds if _stream not in rd]
        return rds
    
    
    def find_datasets(self, site: str, node: str, instrument: str, stream: str, data_only: bool = True) -> list:
        """
        Find datasets if they contain supplied strings.

        :param site: An 8-character OOI designator for a site.
            It is recommended to use the full site designator.
        :param node: A 5-character designator for a site. Partial strings acceptable.
        :param instrument: A 12-character designator for an instrument. Partial strings acceptable.
        :param stream: A variable character designator for an instrument stream. Partial strings acceptable.
        :param data_only: If True, excludes engineering/calibration datasets.
        :return:
        """

        rds = self.list_all_datasets(data_only = data_only)
        rds = [rd for rd in rds if site in rd and node in rd and instrument in rd and stream in rd]
        return rds
    
    
    def list_dataset_files(self, rd: str, file_ext: str = '.nc') -> list:
        """
        List dataset files if they contain the corresponding reference designator and file extension.
        :param rd: The reference designator of interest.
        :param file_ext: The file extension. OOI kdata science data comes in netcdfs (.nc).
        :return:
        """
        
        files = sorted(self.__local.glob(self.BASE_DIR + rd + f'/*{file_ext}'))
        dtrfiles =[]
        for file in files:
            if len(re.findall(rd,file)) != 2: # Drop ancillary datasets.
                continue
            file_begin, file_end = [datetime.strptime(dt,'%Y%m%dT%H%M%S') for dt in re.findall("(\d{8}T\d{6})", file)]
            if (file_begin <= self.bdt <= file_end or file_begin <= self.edt <= file_end or 
                self.bdt <= file_begin <= self.edt or self.bdt <= file_end <= self.edt):            
                dtrfiles.append(file)        
        return dtrfiles
    
    
    def find_files(self) -> list:
        """
        A wrapper for finding files based on instance inputs.

        :return: A list of filepaths that match user inputs.
        """
        rds = self.find_datasets(*self._usnims)
        files = [] 
        for rd in rds:
            rdfiles = self.list_dataset_files(rd)
            files = files + rdfiles
        return files
            
    
    def import_file(self, filepath: os.path.abspath) -> xr.Dataset:
        """
        Import a file and perform minor cleaning.

        :param filepath: An absolute filepath to the file of interest.
        :return: An xarray dataset containing slightly cleaned data.
        """

        ds = xr.open_dataset(filepath)
        ds = ds.swap_dims({'obs':'time'})
        ds = ds.drop_vars(['obs'], errors = 'ignore')
        ds = ds.drop_vars(DROP_SYS_VARS, errors = 'ignore') # Remove confusing OOI backend variables.
        ds = drop_qc_test_vars(ds) # Remove custom qc test variables. Deprecated???

        # If lat is not a variable, make lat and lon variables from attribute data.
        # Primarily intended for the computation of absolute salinity by the ctd module.
        if 'lat' not in ds.data_vars and 'latitude' not in ds.data_vars and 'lat' not in ds.coords:
            ds['latitude'] = ds.lat
            ds['longitude'] = ds.lon
        elif 'lat' in ds.coords:
            ds['latitude'] = ds.lat
            ds['longitude'] = ds.lon
            ds = ds.drop(['lat','lon'])
            
            
        ds = ds.sortby('time')
        return ds
    
    
    def open_datasets(self) -> list:
        """
        Import and open datasets found in the files attribute.

        :return: A list of xarray datasets.
        """
        with multiprocessing.Pool(len(self.files)) as pool:
            ds_list = pool.map(self.import_file, self.files)
        return ds_list
    
    
    def combine_data(self, ds_list: list) -> xr.Dataset:
        """
        Combine a list of xarray datasets into a single dataset.

        :param ds_list: A list of datasets.
        :return: A combined or concatenated dataset, sliced by the instance begin_datetime and end_datetime.
        """
        try: # Try to merge by coordinates.
            ds = xr.combine_by_coords(ds_list, combine_attrs = 'drop')
            ds = ds.drop_duplicates(list(ds.dims))
        except: # In some instances, merging my coordinates doesn't work, so concatenation may work, but is slower.
            ds = xr.concat(ds_list, dim = 'time', combine_attrs = 'drop')
        ds = ds.sortby('time')
        ds = ds.drop_duplicates(list(ds.dims))
        ds = ds.sel(time = slice(self.bdt.strftime('%Y-%m-%d %H:%M:%S'),self.edt.strftime('%Y-%m-%d %H:%M:%S')))
        return ds
    
    
    def nan_by_qartod(self, ds: xr.Dataset, nan_flags: list = [4,9], variables: str or list = 'all') -> xr.Dataset:
        """
        NaN QARTOD data if it did not meet user supplied primary flag criteria.

        :param ds: The input xarray dataset.
        :param nan_flags: A list of flags the user wants to NaN.
        :param variables: A string indicating 'all' QARTOD variables or a list of specific qartod variable.
        :return: The output xarray dataset, with data NaNed if it did not meet criteria.
        """
        if variables == 'all':
            qecs = [v for v in list(ds.data_vars) if 'qartod_executed' in v]
        else:
            qecs = [v + '_qartod_executed' for v in variables]
        for qec in qecs:
            qdc = qec.split('_qartod_executed')[0]
            qrc = qdc + '_qartod_results'   
            #ds[qdc] = ds[qdc].where(ds[qec] != 1 & ~ds[qrc].isin(nan_flags), np.nan) # Might be deprecated if OOI does not change executed flags to single digits.
            ds[qdc] = ds[qdc].where(~ds[qrc].isin(nan_flags), np.nan)
        return ds