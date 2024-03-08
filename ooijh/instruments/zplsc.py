"""IN DEVELOPMENT"""

from datetime import datetime
import echopype as ep
import fsspec
import os
import re
import warnings
import xarray as xr

from ooijh.core import _USER_DIR

class ZPLSC():
    
    SITE_PATH = {'CE02SHBP': f"{_USER_DIR}/ooi/rsn_cabled/rsn_data/DVT_Data/mj01c/ZPLSCB101_10.33.13.7",
                 'CE04OSPS': f"{_USER_DIR}/ooi/rsn_cabled/rsn_data/DVT_Data/pc01b/ZPLSCB102_10.33.10.143"}
    
    def __init__(self, site: str, begin_datetime: datetime, end_datetime: datetime):
        self.bdt = begin_datetime
        self.edt = end_datetime
        self._base_dir = self.SITE_PATH[site.upper()]
        self.files = self.find_files()
    
    
    def find_files(self) -> list:
        local = fsspec.filesystem('file')
        files = []
        year_dirs = [yd for yd in local.glob(self._base_dir + '/*') if 'log' not in yd and 'config' not in yd and '_Cal' not in yd and 'Collier' not in yd]
        year_dirs = [yd for yd in year_dirs if self.bdt.year <= int(yd[-4:]) <= self.edt.year]
        for yd in year_dirs:
            month_dirs = [md for md in local.glob(yd + '/*')]
            for md in month_dirs: 
                day_dirs = [dd for dd in local.glob(md + '/*') if self.bdt <= datetime.strptime(dd[-10:],'%Y/%m/%d') <= self.edt]
                for dd in day_dirs:
                    day_files = local.glob(dd + '/*.raw')
                    for day_file in day_files:
                        filename = os.path.basename(day_file)
                        file_dt = datetime.strptime(re.findall('OOI-D(\d{8}-T\d{6}).raw',filename)[0], '%Y%m%d-T%H%M%S')
                        if self.bdt <= file_dt <= self.edt:
                            files.append(day_file)
        files = sorted(files)
        return files
        
        
    def open_file(self, filepath: os.path.abspath) -> object:
        ed = ep.open_raw(raw_file = filepath, sonar_model = 'ek60')
        return ed
    
    
    def combine_data(self, filepaths: list) -> object:
        ed_list = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore") # Catch user warnings that echopype may throw. Usually time related warnings.
            for fp in filepaths:
                _ed = self.open_file(fp)
                if len(_ed.platform.channel) == 3:
                    ed_list.append(_ed)
                else:
                    msg = f"File: {fp} only has {len(_ed.platform.channel)} channels. This file will not be added to the combined dataset."
                    warnings.warn(msg)
            ed = ep.combine_echodata(ed_list)
        return ed
    
    
    def process(self, ed, depth_bin: float, time_bin: str) -> xr.Dataset:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore") # Catch user warnings that echopype may throw. Usually empty slice warnings.
            sv = ep.calibrate.compute_Sv(ed).compute()
            sv = sv.sortby('ping_time')
            mvbs = ep.commongrid.compute_MVBS(sv, range_meter_bin = depth_bin, ping_time_bin = time_bin)
            mvbs = mvbs.assign_coords(depth = ('echo_range', mvbs['echo_range'].values[::-1]))
            mvbs = mvbs.swap_dims({'echo_range':'depth'})
            mvbs = ep.consolidate.swap_dims_channel_frequency(mvbs)
            mvbs = mvbs.drop_vars(['channel'], errors = 'ignore')
            mvbs = mvbs.drop_vars(['echo_range'], errors = 'ignore')
            mvbs = mvbs.rename({'Sv':'sv','ping_time':'time','frequency_nominal':'frequency'})
        return mvbs

    
    def get_data(self, depth_bin: float = 0.2, time_bin: str = '10s') -> xr.Dataset:
        ed = self.combine_data(self.files)
        mvbs = self.process(ed, depth_bin, time_bin)
        return mvbs
        