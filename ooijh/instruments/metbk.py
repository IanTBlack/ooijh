from datetime import datetime
import gsw
import multiprocessing
import xarray as xr
import os

from ooijh.core import KDATA
from ooijh.drops import DROP_METBK_VARS, drop_qartod_test_vars


class METBK(KDATA):
    def __init__(self, site: str, node: str, instrument: str = 'METBK', stream: str = 'metbk_a_dcl',
                 begin_datetime : datetime = datetime(2014,1,1), end_datetime: datetime = datetime(2040,12,31,23,59,59),
                 process_qartod: bool = True, nan_flags: list = [4,9], drop_qartod: bool = True):

        """
        A class for obtaining preprocessed and processed METBK data.

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
                         begin_datetime, end_datetime, process_qartod, nan_flags, drop_qartod)
    
    
    def preprocess(self, ds_list: list) -> list:
        """
        The preprocess function is for operations or actions that may be dataset dependent.

        :param ds_list: A list of xarray datasets.
        :return: A list of preprocessed xarray datasets.
        """
        ds_list = [ds.drop_vars(DROP_METBK_VARS, errors = 'ignore') for ds in ds_list]
        if self.process_qartod is True:
            ds_list = [self.nan_by_qartod(ds, self.nan_flags) for ds in ds_list]
        if self.drop_qartod is True:
            with multiprocessing.Pool(os.cpu_count()-1, maxtasksperchild = 1) as pool:
                ds_list = pool.map(drop_qartod_test_vars, ds_list)
        return ds_list
    
        
    def process(self, ds: xr.Dataset) -> xr.Dataset:
        """
        The process function is for operations or actions that can be performed after multiple datasets of the same
            type are combined.

        :param ds: The input xarray dataset.
        :return: The output processed xarray dataset.
        """

        __rename = {'met_salsurf': 'sea_water_practical_salinity',
                   'met_current_direction': 'current_direction',
                   'met_current_speed':'current_speed',
                   'eastward_velocity': 'sea_water_eastward_velocity',
                   'northward_velocity': 'sea_water_northward_velocity',
                   'met_barpres': 'barometric_pressure',
                    'met_relwind_speed': 'relative_wind_speed',
                    'met_relwind_direction': 'relative_wind_direction',
                   'sea_surface_temperature': 'sea_water_temperature'}
        ds = ds.rename(__rename)
    
        # Compute additional data products.
        ds['sea_water_pressure'] = (['time'],gsw.p_from_z(ds.ct_depth.values * -1, ds.latitude.values))
        ds = ds.drop_vars(['ct_depth'], errors = 'ignore')
        ds['sea_water_absolute_salinity'] = (['time'],gsw.SA_from_SP(ds.sea_water_practical_salinity.values, 
                                                             ds.sea_water_pressure.values, 
                                                             ds.longitude.values, ds.latitude.values))
        ds['sea_water_conservative_temperature'] = (['time'],gsw.CT_from_t(ds.sea_water_absolute_salinity.values, 
                                                                   ds.sea_water_temperature.values, 
                                                                   ds.sea_water_pressure.values))
        ds['sea_water_density'] = (['time'],gsw.density.rho(ds.sea_water_absolute_salinity.values, 
                                                    ds.sea_water_conservative_temperature.values, 
                                                    ds.sea_water_pressure.values)) #Overwrites OOI sea_water_density derivation.
        
        
        ds['sea_water_absolute_salinity'].attrs['units'] = 'ASU'
        ds['sea_water_absolute_salinity'].attrs['units_tex'] = r'$\frac{g}{kg}$'
        ds['sea_water_absolute_salinity'].attrs['description'] = 'The mass fraction of salt in seawater. Defined under TEOS-10.'

        ds['sea_water_conservative_temperature'].attrs['units'] = 'degrees Celsius'
        ds['sea_water_conservative_temperature'].attrs['units_tex'] = r'$^{\circ}C$'
        ds['sea_water_conservative_temperature'].attrs['description'] = 'Conservative temperature represents the heat content of the ocean conserved during mixing and pressure changes. Further defined under TEOS-10.'
        
        ds['sea_water_density'].attrs['units'] = 'kilograms per cubic meter'
        ds['sea_water_density'].attrs['units_tex'] = r'$\frac{kg}{m^3}$'
        ds['sea_water_density'].attrs['description'] = 'Sea water density calculated under TEOS-10.'
        
        
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
        