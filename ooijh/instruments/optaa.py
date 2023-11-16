from datetime import datetime
import multiprocessing
import numpy as np
import os
import re
from scipy import interpolate
from struct import calcsize
import warnings
import xarray as xr

from ooijh.core import KDATA, _USER_DIR
from ooijh.drops import DROP_OPTAA_VARS, drop_qartod_test_vars


class OPTAA(KDATA):
    def __init__(self, site: str, node: str, instrument: str = 'OPTAA', stream: str = 'optaa',
                 begin_datetime: datetime = datetime(2014,1,1), end_datetime: datetime = datetime(2040,12,31,23,59,59),
                 process_qartod: bool = True, nan_flags: list = [4,9], drop_qartod: bool = True,
                 dev_filepaths: None or list = None, tscor_filepath: None or os.path.abspath = None):
        """
        A class for obtaining preprocessed and processed OPTAA data.

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
        :param dev_filepaths: A list of filepaths, where each filepath matches the index of the files attribute.
        :param tscor_filepath: The filepath to an ACS TS4.cor file.
        """

        super().__init__(site.upper(), node.upper(), instrument.upper(), stream.lower(), 
                         begin_datetime, end_datetime, process_qartod, nan_flags, drop_qartod)
        
        if self.dev_filepaths is None:
            raise FileNotFoundError('No dev files supplied.')
        else:
            self.dev_filepaths = dev_filepaths #Filepath index much match file index.
        if tscor_filepath is None:
            warnings.warn('No TS4.cor file supplied. Defaulting to internal package file.')
            self.tscor_filepath = f"{_USER_DIR}/ooijh/__cals__/TS4.cor"
        else:
            self.tscor_filepath = tscor_filepath
    
    
    def preprocess(self, ds_list: list) -> list:
        """
        The preprocess function is for operations or actions that may be dataset dependent.

        :param ds_list: A list of xarray datasets.
        :return: A list of preprocessed xarray datasets.
        """

        if len(ds_list) != len(self.dev_filepaths):
            raise ValueError('dev_filepaths must be a list of filepaths equivalent to each file in KDATA.files.')
        
        # Preprocess data.
        ds_list = [ds.drop_vars(DROP_OPTAA_VARS, errors = 'ignore') for ds in ds_list]
        if self.process_qartod is True:
            ds_list = [self.nan_by_qartod(ds, self.nan_flags) for ds in ds_list]
        if self.drop_qartod is True:
            with multiprocessing.Pool(len(ds_list)) as pool:
                ds_list = pool.map(drop_qartod_test_vars, ds_list)
        
        # Processes each file based on supplied dev files.
        preprocessed_ds_list = []    
        for i in range(len(ds_list)):
            ds = ds_list[i]
            dev = Dev(self.dev_filepaths[i])
            tscor = TSCor(self.tscor_filepath)
            __rename = {'external_temp_raw': 'raw_external_temperature',
                       'internal_temp_raw': 'raw_internal_temperature',
                       'a_signal_counts': 'a_signal',
                       'c_signal_counts': 'c_signal',
                       'a_reference_counts': 'a_reference',
                       'c_reference_counts': 'c_reference',
                       'elapsed_run_time': 'elapsed_time',
                       'c_reference_dark_counts': 'c_reference_dark',
                       'a_reference_dark_counts': 'a_reference_dark',
                       'c_signal_dark_counts': 'c_signal_dark',
                       'a_signal_dark_counts': 'a_signal_dark'}
            ds = ds.rename(__rename)
            ds = ds.assign_coords({'wavelength_a': np.unique(ds.wavelength_a), 'wavelength_c':np.unique(ds.wavelength_c)})
            for a in ['a_signal','a_reference']:
                ds[a] = (['time','wavelength_a'], ds[a].data)
            for c in ['c_signal','c_reference']:
                ds[c] = (['time','wavelength_c'], ds[c].data)
            ds = ds.drop_dims('wavelength')
            ds['internal_temperature'] = compute_internal_temperature(ds.raw_internal_temperature)
            ds['external_temperature'] = compute_external_temperature(ds.raw_external_temperature)
            ds['a_uncorr'] = compute_uncorrected(ds.a_signal, ds.a_reference, dev)
            ds['c_uncorr'] = compute_uncorrected(ds.c_signal, ds.c_reference, dev)
            ds['a_pg'] = compute_pg(ds.a_uncorr, 'a', ds.internal_temperature, dev)
            ds['c_pg'] = compute_pg(ds.c_uncorr, 'c', ds.internal_temperature, dev)
            ds['a_pg_ts'] = compute_ts(ds.a_pg, 'a', ds.sea_water_temperature, 
                                       ds.sea_water_practical_salinity, dev, tscor)
            ds['c_pg_ts'] = compute_ts(ds.c_pg, 'c', ds.sea_water_temperature, 
                                       ds.sea_water_practical_salinity, dev, tscor)

            # Interpolation to common wavelength bins make scattering correction easier.
            wvls = np.arange(410, 721,1) # Most ACS' measure between 410 and 720 nm.
            ds = ds.interp({'wavelength_a': wvls, 'wavelength_c': wvls}) # Interpolate to 1nm bins.
            ds = ds.reset_index(['wavelength_a',
                                 'wavelength_c'], drop = True).assign_coords(wavelength = wvls).rename({'wavelength_a':'wavelength',
                                                                                                        'wavelength_c':'wavelength'}) # Rename wavelength dimension.
            ds['a_pg_ts_baseline'] = scattering_correction_baseline(ds.a_pg_ts)
            ds['a_pg_ts_fixed'] = scattering_correction_fixed(ds.a_pg_ts, ds.c_pg_ts)
            ds['a_pg_ts_proportional'] = scattering_correction_proportional(ds.a_pg_ts, ds.c_pg_ts)
            preprocessed_ds_list.append(ds)
        return preprocessed_ds_list
    
    def process(self, ds: xr.Dataset) -> xr.Dataset:
        """
        The process function is for operations or actions that can be performed after multiple datasets of the same
            type are combined.

        :param ds: The input xarray dataset.
        :return: The output processed xarray dataset.
        """
        return ds
        
        
    def get_data(self) -> xr.Dataset:
        """
        The get_data function performs the preprocess and process functions and serves up data to the end user.

        :return: An xarray dataset containing the data of interest.
        """
        ds_list = self.open_datasets() # Open all datasets.
        ds_list = self.preprocess(ds_list, self.dev_filepaths, self.tscor_filepath)
        ds = self.combine_data(ds_list)
        ds = self.process(ds)
        ds = ds[sorted(ds.data_vars)] # Sort variables alphabetically because it is less obtrusive.
        return ds
        
        
class Dev():
    """
    A class for parsing ACS calibration (.dev) files. 
    These files are necessary for performing corrections.
  
    This class can be instantiated and data can be accessed as class attributes.
    Alternatively, the to_ds() function will format calibration data so that it can be used with
    operations based in xarray. The to_nc() function will export calibration data as a netcdf.
    """
    def __init__(self, filepath: os.path.abspath) -> None:
        """
        Parse the .dev file.
        
        :param filepath: The location of the dev file.
        """
        self.filepath = os.path.normpath(filepath)
        self.__read_dev()
        self.__parse_metadata()
        self.__parse_tbins()
        self.__parse_offsets()
        self.__check_parse()
        self.__build_packet_header()

    def __read_dev(self) -> None:
        """Import the .dev file as a text file."""
        
        with open(self.filepath, 'r') as _file:
            self._lines = _file.readlines()

    def __parse_metadata(self) -> None:
        """Parse the .dev file for sensor metadata."""
        
        metadata_lines = [line for line in self._lines if 'C and A offset' not in line]
        for line in metadata_lines:
            if 'ACS Meter' in line:
                self.sensor_type = re.findall('(.*?)\n', line)[0]
            elif 'Serial' in line:
                self.sn_hexdec = re.findall('(.*?)\t', line)[0]
                self.sn = 'ACS' + str(int(self.sn_hexdec[-6:], 16)).zfill(3)
            elif 'structure version' in line:
                self.structure_version = int(re.findall('(.*?)\t', line)[0])
            elif 'tcal' in line:
                self.tcal, self.ical = [float(v) for v in re.findall(': [+-]?([0-9]*[.]?[0-9]+) C', line)]
                cal_date_str = re.findall('file on (.*?)[.].*?\n', line)[0].replace(' ', '')
                try:
                    self.cal_date = datetime.strptime(cal_date_str, '%m/%d/%Y').strftime('%Y-%m-%d')
                except:
                    self.cal_date = datetime.strptime(cal_date_str, '%m/%d/%y').strftime('%Y-%m-%d')
            elif 'Depth calibration' in line:
                (self.depth_cal_1,
                 self.depth_cal_2) = [float(v) for v in re.findall('[+-]?([0-9]*[.]?[0-9]+)\t', line)]
            elif 'Baud' in line:
                self.baudrate = int(re.findall('(.*?)\t', line)[0])
            elif 'Path' in line:
                self.path_length = float(re.findall('(.*?)\t', line)[0])
            elif 'wavelengths' in line:
                self.output_wavelengths = int(re.findall('(.*?)\t', line)[0])
            elif 'number of temperature bins' in line:
                self.num_tbins = int(re.findall('(.*?)\t', line)[0])
            elif 'maxANoise' in line:
                (self.max_a_noise, self.max_c_noise, self.max_a_nonconform, self.max_c_nonconform,
                 self.max_a_difference, self.max_c_difference, self.min_a_counts,
                 self.min_c_counts, self.min_r_counts, self.max_tempsdev,
                 self.max_depth_sdev) = [float(v) for v in re.findall('[+-]?([0-9]*[.]?[0-9]+)\t', line)]

                
    def __parse_tbins(self) -> None:
        """Parse the .dev file for temperature bin information."""
        
        line = [line for line in self._lines if '; temperature bins' in line][0]
        tbins = line.split('\t')
        tbins = [v for v in tbins if v]
        tbins = [v for v in tbins if v != '\n']
        tbins = [float(v) for v in tbins if 'temperature bins' not in v]
        self.tbins = np.array(tbins)

    def __parse_offsets(self) -> None:
        """Parse the .dev file for a and c offsets."""
        
        offset_lines = [line for line in self._lines if 'C and A offset' in line]
        c_wavelengths = []
        a_wavelengths = []
        c_offsets = []
        a_offsets = []
        c_deltas = []
        a_deltas = []
        for line in offset_lines:
            offsets, c_delta, a_delta = line.split('\t\t')[:-1]
            wavelength_c, wavelength_a, _, offset_c, offset_a = offsets.split('\t')
            wavelength_c = float(wavelength_c.replace('C', ''))
            wavelength_a = float(wavelength_a.replace('A', ''))
            c_wavelengths.append(wavelength_c)
            a_wavelengths.append(wavelength_a)
            offset_c = float(offset_c)
            offset_a = float(offset_a)
            c_offsets.append(offset_c)
            a_offsets.append(offset_a)
            c_delta = np.array([float(v) for v in c_delta.split('\t')])
            a_delta = np.array([float(v) for v in a_delta.split('\t')])
            c_deltas.append(c_delta)
            a_deltas.append(a_delta)
        self.wavelength_c = np.array(c_wavelengths)
        self.wavelength_a = np.array(a_wavelengths)
        self.offset_c = np.array(c_offsets)
        self.offset_a = np.array(a_offsets)
        self.delta_t_c = np.array(c_deltas)
        self.delta_t_a = np.array(a_deltas)
        self.f_delta_t_c = interpolate.interp1d(self.tbins, self.delta_t_c, axis=1, assume_sorted=True, copy=False,
                                                bounds_error=False,
                                                fill_value=(self.delta_t_c[:, 1], self.delta_t_c[:, -1]))
        self.f_delta_t_a = interpolate.interp1d(self.tbins, self.delta_t_a, axis=1, assume_sorted=True, copy=False,
                                                bounds_error=False,
                                                fill_value=(self.delta_t_a[:, 1], self.delta_t_a[:, -1]))

    def __build_packet_header(self) -> None:
        """
        Build a packet descriptor for parsing binary ACS packets.
        Only used when reading raw binary from a file or over serial.
        """
        
        self.PACKET_REGISTRATION = b'\xff\x00\xff\x00'
        self.LEN_PACKET_REGISTRATION = len(self.PACKET_REGISTRATION)
        self.PACKET_HEADER = '!HBBlHHHHHHHIBB'
        self.LEN_PACKET_HEADER = calcsize(self.PACKET_HEADER)

        self.packet_header = self.PACKET_HEADER
        for i in range(self.output_wavelengths):
            self.packet_header += 'HHHH'
        self.packet_length = self.LEN_PACKET_REGISTRATION + calcsize(self.packet_header)

        
    def __check_parse(self) -> None:
        """Verify that the parse obtained the correct informatoin."""
        
        if len(self.wavelength_c) != len(self.wavelength_a):
            raise ValueError('Mismatch between number of wavelengths extracted for A and C.')
        if self.delta_t_c.shape != (len(self.wavelength_c), self.num_tbins):
            raise ValueError('Mismatch between length of C wavelengths and number of temperature bins.')
        if self.delta_t_a.shape != (len(self.wavelength_a), self.num_tbins):
            raise ValueError('Mismatch between length of A wavelengths and number of temperature bins.')

            
    def to_ds(self) -> xr.Dataset:
        """
        Export class attributes as an xr.Dataset.
        
        :return: An xarray dataset containing calibration information.
        """
        
        ds = xr.Dataset()
        ds = ds.assign_coords({'wavelength_a': self.wavelength_a})
        ds = ds.assign_coords({'wavelength_c': self.wavelength_c})
        ds = ds.assign_coords({'temperature_bins': self.tbins})

        ds['offsets_a'] = (['wavelength_a'], np.array(self.offset_a))
        ds['delta_t_a'] = (['wavelength_a', 'temperature_bins'], np.array(self.delta_t_a))
        
        ds['offsets_c'] = (['wavelength_c'], np.array(self.offset_c))
        ds['delta_t_c'] = (['wavelength_c', 'temperature_bins'], np.array(self.delta_t_c))

        ds.attrs['sensor_type'] = self.sensor_type
        ds.attrs['serial_number'] = self.sn
        ds.attrs['factor_calibration_date'] = self.cal_date
        ds.attrs['output_wavelengths'] = self.output_wavelengths
        ds.attrs['number_temp_bins'] = self.num_tbins
        ds.attrs['path_length'] = self.path_length
        ds.attrs['tcal'] = self.tcal
        ds.attrs['ical'] = self.ical
        ds.attrs['baudrate'] = self.baudrate
        ds.attrs['dev_structure_version'] = self.structure_version
        return ds


    
    def to_nc(self, out_filepath: os.path.abspath) -> None:
        """
        Export .dev data as a netcdf.

        :param out_filepath: 
        """

        split = os.path.splitext(out_filepath)
        if split[-1] != '.nc':
            out_filepath += '.nc'
        ds = self.to_ds()
        ds.to_netcdf(out_filepath, engine = 'netcdf4')
        

class TSCor():
    def __init__(self, filepath: os.path.abspath) -> None:
        """
        Parse the .cor file and assign data as attributes.

        :param filepath: The filepath of the TS4.cor file.
        """
        self.filepath = os.path.normpath(filepath)
        self.__read_cor()
        self.__parse_lines()
    
    
    def __read_cor(self):
        with open(self.filepath, 'r') as _file:
            self._lines = _file.readlines()


    def __parse_lines(self):
        """Parse the lines of the .cor file to get correction information."""

        wavelengths = []
        psi_t = []
        psi_s_c = []
        psi_s_a = []
        for line in self._lines:
            line_data = line.split('\t')
            line_data = [v.replace('\n', '') for v in line_data]
            line_data = [v.replace(' ', '') for v in line_data]
            if line_data == ['']:
                break
            line_data = [float(v) for v in line_data]
            wavelengths.append(line_data[0])
            psi_t.append(line_data[1])
            psi_s_c.append(line_data[2])
            psi_s_a.append(line_data[3])
        if len(wavelengths) != len(psi_t) != len(psi_s_c) != len(psi_s_a):
            raise ValueError('Mismatch in length of TS4cor file.')
        else:
            self.wavelengths = wavelengths
            self.psi_t = psi_t
            self.psi_s_c = psi_s_c
            self.psi_s_a = psi_s_a


    def to_ds(self) -> xr.Dataset:
        """
        Export class attributes to an xarray dataset.

        :return: An xarray dataset containing correction data.
        """
        ds = xr.Dataset()
        ds = ds.assign_coords({'wavelength': self.wavelengths})
        ds['psi_t'] = (['wavelength'], self.psi_t)
        ds['psi_s_c'] = (['wavelength'], self.psi_s_c)
        ds['psi_s_a'] = (['wavelength'], self.psi_s_a)
        return ds


    def to_nc(self, out_filepath: os.path.abspath) -> None:
        """
        Export .cor data as a netcdf.

        :param out_filepath: The save location of netcdf containing .cor information.
        """

        split = os.path.splitext(out_filepath)
        if split[-1] != '.nc': 
            out_filepath += '.nc'
        ds = self.to_ds()
        ds.to_netcdf(out_filepath, engine = 'netcdf4')

        
def compute_internal_temperature(counts: xr.DataArray) -> xr.DataArray:
    """
    Compute internal temperature from the raw thermistor counts.
    
    :param counts: Raw thermistor counts.
    :return: An xr.DataArray of internal temperature.
    """
    counts = counts.astype('int64')
    volts = 5 * counts / 65535
    resistance = 10000 * volts / (4.516 - volts)
    internal_temperature = 1 / (
                0.00093135 + 0.000221631 * np.log(resistance) + 0.000000125741 * np.log(resistance) ** 3) - 273.15
    return internal_temperature


def compute_external_temperature(counts: xr.DataArray) -> xr.DataArray:
    """
    Compute external temperature from the raw thermistor counts.
    
    :param counts: Raw thermistor counts.
    :return: An xr.DataArray of external temperature.
    """
    counts = counts.astype('int64')
    a = -7.1023317e-13
    b = 7.09341920e-08
    c = -3.87065673e-03
    d = 95.8241397
    external_temperature = a * counts ** 3 + b * counts ** 2 + c * counts + d
    return external_temperature


def compute_uncorrected(signal_counts: xr.DataArray, reference_counts: xr.DataArray, dev: object or None) -> xr.DataArray:
    """
    Compute uncorrected a or c values. 
    
    :param signal_counts: Raw signal counts from the ACS.
    :param reference_counts: Reference counts from the ACS.
    :param dev: The corresponding Dev object for the dataset.
    :return: An xarray of uncorrected data (no delta t or water offset applied).
    """
    if dev is None:
        x = 0.25 # Assume a standard ACS path length of 0.25 meters.
    else:
        x = dev.path_length
    uncorr = (1 / x) * np.log(signal_counts / reference_counts)
    return uncorr


def compute_pg(uncorrected: xr.DataArray, channel: str, internal_temperature: xr.DataArray, dev: object) -> xr.DataArray:
    """
    Compute the measured absorption and attenuation. pg refers the effect of particles and gelbstoff.
    In the ACS manual, this value is called a_m/c_m. The absorption or attenuation with the effect of water removed.
    
    :param uncorrected: The uncorrected signal.
    :param channel: 'a' or 'c', used to determing the delta_t function.
    :param internal_temperature: The internal temperature of the sensor.
    :param dev: The corresponding Dev object for the dataset.
    """
    
    if channel.lower() == 'a':
        delta_t = dev.f_delta_t_a(internal_temperature).T
        offsets = dev.offset_a
    elif channel.lower() == 'c':
        delta_t = dev.f_delta_t_c(internal_temperature).T
        offsets = dev.offset_c
    offsets, _ = np.meshgrid(offsets, internal_temperature)
    pg = (offsets - uncorrected) - delta_t
    return pg


def compute_ts(pg: xr.DataArray, channel: str, 
                  temperature: xr.DataArray, salinity: xr.DataArray,
                  dev: object, tscor: object) -> xr.DataArray:
    """
    Compute temperature-salinity corrected data.
    
    :param pg: The measured or particlate/gelbstoff signal.
    :param channel: 'a' or 'c', used for obtaining the correct correction coeffs.
    :param temperature: sea water temperature. conservative temperature is acceptable per the manual
    :param salinity: sea water salinity. absolute salinity is acceptable per the manual.
    :param dev: The corresponding Dev object for the dataset.
    :param tscor: The corresponding TSCor object for the dataset.
    
    :return: An xr.DataArray containing temperature-salinity corrected data.
    """
    tsds = tscor.to_ds()
    if channel.lower() == 'a':
        wvls = list(pg.wavelength_a.values)
        tsds = tsds.sel(wavelength=wvls, method='nearest')
        psi_temp, delta_t = np.meshgrid(tsds.psi_t.values, temperature.values - dev.tcal)
        psi_sal, s = np.meshgrid(tsds.psi_s_a.values, salinity.values)
        pg_ts = pg - ((psi_temp * delta_t) + (psi_sal * s))
    elif channel.lower() == 'c':
        wvls = list(pg.wavelength_c.values)
        tsds = tsds.sel(wavelength=wvls, method='nearest')
        psi_temp, delta_t = np.meshgrid(tsds.psi_t.values, temperature.values - dev.tcal)
        psi_sal, s = np.meshgrid(tsds.psi_s_c.values, salinity.values)
        pg_ts = pg - ((psi_temp * delta_t) + (psi_sal * s))
    return pg_ts


def scattering_correction_baseline(a_pg_ts: xr.DataArray, reference_wavelength: int = 715):
    """
    METHOD 1
    Perform baseline scattering correction on absorption data. 
    Although this function will work with data that is not on a common wavelength bin, 
    it is recommended that you apply this function with interpolated data.
     
    :param a_pg_ts: TS corrected absorption.
    :param reference_wavelength: The reference wavelength. Usually 715.
    :return: Absorption corrected using the baseline method.
    """
    ref = a_pg_ts.sel(wavelength=reference_wavelength, method='nearest')
    scatcorr = a_pg_ts - ref
    return scatcorr


def scattering_correction_fixed(a_pg_ts: xr.DataArray, c_pg_ts: xr.DataArray, F: float = 0.18) -> xr.DataArray:
    """
    METHOD 2
    Perform fixed scattering correction on absorption data.
    Use of this method requires that the absorption and attenuation data be on common wavelength bins.
    
    :param a_pg_ts: TS corrected absorption.
    :param c_pg_ts: TS corrected attenuation.
    :param F: The fixed offset.
    :return: Absorption corrected using the baseline method.
    """
    scatcorr = a_pg_ts - F * (c_pg_ts - a_pg_ts)
    return scatcorr


def scattering_correction_proportional(a_pg_ts: xr.DataArray, c_pg_ts: xr.DataArray,
                                       reference_wavelength: int = 715) -> xr.DataArray:
    """
    METHOD 3
    Perform proportional scattering correction on absorption data.
    Use of this method requires that the absorption and attenuation data be on common wavelength bins.
    
    :param a_pg_ts: TS corrected absorption.
    :param c_pg_ts: TS corrected attenuation.
    :param reference_wavelength: The reference wavelength. Usually 715.
    :return: Absorption corrected using the baseline method.
    """
    ref_a = a_pg_ts.sel(wavelength=reference_wavelength, method='nearest')
    ref_c = c_pg_ts.sel(wavelength=reference_wavelength, method='nearest')
    scatcorr = a_pg_ts - ((ref_a / (ref_c - ref_a)) * (c_pg_ts - a_pg_ts))
    return scatcorr