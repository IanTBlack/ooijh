from datetime import datetime
import gsw
import multiprocessing
import numpy as np
import xarray as xr

from ooijh.core import KDATA
from ooijh.drops import DROP_CTD_VARS, drop_qartod_test_vars


class CTD(KDATA):
    def __init__(self, site, node, instrument = 'CTD', stream = 'ctd', 
                 begin_datetime = datetime(2014,1,1), end_datetime = datetime(2040,12,31,23,59,59),
                process_qartod_vars = True, drop_qartod_test_vars = True):
        super().__init__(site.upper(), node.upper(), instrument.upper(), stream.lower(), begin_datetime, end_datetime)
    
    
    def preprocess(self, ds_list):
        ds_list = self.open_datasets()
        ds_list = [ds.drop_vars(DROP_CTD_VARS, errors = 'ignore') for ds in ds_list]
        if self.process_qartod is True:
            ds_list = [self.nan_by_qartod(ds, self.nan_flags) for ds in ds_list]
        if self.drop_qartod is True:
            with multiprocessing.Pool(len(ds_list)) as pool:
                ds_list = pool.map(drop_qartod_test_vars, ds_list)
        return ds_list
    
  
    
    def process(self, ds):
        # Remove bad data.
        ds['sea_water_practical_salinity'] = xr.where((ds.sea_water_practical_salinity < 2) & (ds.sea_water_practical_salinity > 2) , np.nan, ds.sea_water_practical_salinity)
        ds['sea_water_temperature'] = xr.where((ds.sea_water_temperature < -5) & (ds.sea_water_temperature > 50) , np.nan, ds.sea_water_temperature)
       
    
        # Compute additional data products.
        ds['sea_water_absolute_salinity'] = gsw.SA_from_SP(ds.sea_water_practical_salinity, ds.sea_water_pressure, ds.longitude, ds.latitude)
        ds['sea_water_conservative_temperature'] = gsw.CT_from_t(ds.sea_water_absolute_salinity, ds.sea_water_temperature, ds.sea_water_pressure)
        ds['sea_water_density'] = gsw.density.rho(ds.sea_water_absolute_salinity, ds.sea_water_conservative_temperature, ds.sea_water_pressure) #Overwrites OOI sea_water_density derivation.
    
        return ds
        
        
        
    def get_data(self):
        ds_list = self.open_datasets() # Open all datasets.
        ds_list = self.preprocess(ds_list)
        ds = self.combine_data(ds_list)
        ds = self.process(ds)
        ds = ds[sorted(ds.data_vars)] # Sort variables alphabetically because it is less obtrusive.
        return ds
        

    
    
    
#     def process(self):
#         # COMBINE DATA GOES HERE
#         ds_list = self.open_datasets()
#         ds_list = [ds.drop_vars(DROP_CTD_VARS, errors = 'ignore') for ds in ds_list]
#         if self.process_qartod is True:
#             ds_list = [self.nan_by_qartod(ds, self.nan_flags) for ds in ds_list]
#         if self.drop_qartod is True:
#             with multiprocessing.Pool(len(ds_list)) as pool:
#                 ds_list = pool.map(drop_qartod_test_vars, ds_list)
#         ds = self.combine_data(ds_list)
        
        
        
#         # ADDITIONAL CLEANUP GOES HERE
#         ds['sea_water_practical_salinity'] = xr.where((ds.sea_water_practical_salinity < 2) & (ds.sea_water_practical_salinity > 2) , np.nan, ds.sea_water_practical_salinity)
#         ds['sea_water_temperature'] = xr.where((ds.sea_water_temperature < -5) & (ds.sea_water_temperature > 50) , np.nan, ds.sea_water_temperature)

        
        
#         # ADDITIONAL PROCESSING GOES HERE
#         ds['sea_water_absolute_salinity'] = gsw.SA_from_SP(ds.sea_water_practical_salinity, ds.sea_water_pressure, ds.longitude, ds.latitude)
#         ds['conservative_temperature'] = gsw.CT_from_t(ds.sea_water_absolute_salinity, ds.sea_water_temperature, ds.sea_water_pressure)
#         ds['sea_water_density'] = gsw.density.rho(ds.sea_water_absolute_salinity, ds.conservative_temperature, ds.sea_water_pressure) #Overwrites OOI sea_water_density derivation.
        
#         # Sort variables alphabetically because it is less obtrusive.
#         return ds
        
        
            