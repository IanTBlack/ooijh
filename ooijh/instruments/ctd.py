from datetime import datetime
import gsw
import multiprocessing
import os
import xarray as xr

from ooijh.core import KDATA
from ooijh.drops import DROP_CTD_VARS, drop_qartod_test_vars


class CTD(KDATA):
    def __init__(self, site: str, node: str, instrument: str = 'CTD', stream: str = 'ctd', 
                 begin_datetime: datetime = datetime(2014,1,1), end_datetime: datetime = datetime(2040,12,31,23,59,59),
                 process_qartod: bool = True, nan_flags: list = [4,9], drop_qartod: bool = True):

        """
        A class for obtaining preprocessed and processed CTD data.

        :param site: An 8-character OOI designator for a site.
            It is recommended to use the full site designator.
        :param node: A 5-character designator for a site. Partial strings acceptable.
        :param instrument: A 12-character designator for an instrument. Partial strings acceptable.
        :param stream: A variable character designator for an instrument stream. Partial strings acceptable.
        :param begin_datetime: Beginning of date range as a datetime object.
        :param end_datetime: Ending of date range as a datetime object.
        :param process_qartod: If True, variables with associated qartod tests will be NaNed if they fall within the
            suppled nan_flags.
        :param nan_flags: Determines which data to NaN.
        :param drop_qartod: Drop the qartod test variables to clean up space.
        """

        super().__init__(site.upper(), node.upper(), instrument.upper(), stream.lower(), 
                         begin_datetime, end_datetime, 
                         process_qartod, nan_flags, drop_qartod)
    
    
    def preprocess(self, ds_list: list) -> list:
        """
        The preprocess function is for operations or actions that may be dataset dependent.

        :param ds_list: A list of xarray datasets.
        :return: A list of preprocessed xarray datasets.
        """
        ds_list = [ds.drop_vars(DROP_CTD_VARS, errors = 'ignore') for ds in ds_list]
        if self.process_qartod is True:
            ds_list = [self.nan_by_qartod(ds, self.nan_flags) for ds in ds_list]
        if self.drop_qartod is True:
            with multiprocessing.Pool(os.cpu_count()-1) as pool:
                ds_list = pool.map(drop_qartod_test_vars, ds_list)
        return ds_list
    
    
    def process(self, ds: xr.Dataset) -> xr.Dataset:
        """
        The process function is for operations or actions that can be performed after multiple datasets of the same
            type are combined.

        :param ds: The input xarray dataset.
        :return: The output processed xarray dataset.
        """

        try:
            ds['sea_water_absolute_salinity'] = (['time'],gsw.SA_from_SP(ds.sea_water_practical_salinity.values, 
                                                               ds.sea_water_pressure.values, 
                                                               ds.longitude.values, ds.latitude.values))
        except:
            ds['sea_water_absolute_salinity'] = (['time'],gsw.SA_from_SP(ds.sea_water_practical_salinity.values, 
                                                               ds.sea_water_pressure.values, 
                                                               ds.lon.values, ds.lat.values))
            
        ds['sea_water_conservative_temperature'] = (['time'], gsw.CT_from_t(ds.sea_water_absolute_salinity.values, 
                                                                 ds.sea_water_temperature.values, 
                                                                 ds.sea_water_pressure.values))
        ds['sea_water_density'] = (['time'], gsw.density.rho(ds.sea_water_absolute_salinity.values, 
                                                  ds.sea_water_conservative_temperature.values, 
                                                  ds.sea_water_pressure.values)) #Overwrites OOI sea_water_density derivation.
        

        return ds
        
        
    def get_data(self) -> xr.Dataset:
        """
        The get_data function performs the preprocess and process functions and serves up data to the end user.

        :return: An xarray dataset containing the data of interest.
        """

        ds_list = self.open_datasets() # Open all datasets.
        ds_list = self.preprocess(ds_list)
        ds = self.combine_data(ds_list)
        ds = self.process(ds)
        ds = ds[sorted(ds.data_vars)] # Sort variables alphabetically because it is less obtrusive.
        return ds
        

