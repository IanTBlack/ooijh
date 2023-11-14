import cv2
from datetime import datetime
import fsspec
import os
import re
import warnings

from ooijh.core import _USER_DIR

class CAMDS():
    
    SITE_PATH = {'CE04OSBP':f"{_USER_DIR}/ooi/rsn_cabled/rsn_data/DVT_Data/lv01c/CAMDSB106_10.33.9.6",
                 'CE02SHBP': f"{_USER_DIR}/ooi/rsn_cabled/rsn_data/DVT_Data/mj01c/CAMDSB107_10.33.13.8",
                 'RS03INT1': f"{_USER_DIR}/ooi/rsn_cabled/rsn_data/DVT_Data/mj03c/CAMDSB303_10.31.8.5",
                'RS01SUM2': f"{_USER_DIR}/ooi/rsn_cabled/rsn_data/DVT_Data/mj01b/CAMDSB103_10.33.7.5"}
    
    def __init__(self, site: str, begin_datetime: datetime, end_datetime: datetime):
        self.bdt = begin_datetime
        self.edt = end_datetime
        self._base_dir = self.SITE_PATH[site.upper()]
        self.files = self.find_files()
    
    def find_files(self):
        local = fsspec.filesystem('file')
        files = []
        year_dirs = [yd for yd in local.glob(self._base_dir + '/*') if 'log' not in yd and 'db' not in yd and 'DS_Store' not in yd and 'EbcAIM' not in yd]
        year_dirs = [yd for yd in year_dirs if self.bdt.year <= int(yd[-4:]) <= self.edt.year]
        for yd in year_dirs:
            month_dirs = [md for md in local.glob(yd + '/*')]
            for md in month_dirs: 
                day_dirs = [dd for dd in local.glob(md + '/*') if self.bdt <= datetime.strptime(dd[-10:],'%Y/%m/%d') <= self.edt]
                for dd in day_dirs:
                    day_files = local.glob(dd + '/*.jpg')
                    if len(day_files) == 0:
                        day_files = local.glob(dd + '/*.png')
                    for day_file in day_files:
                        filename = os.path.basename(day_file)
                        file_dt = datetime.strptime(re.findall('(\d{8}T\d{6})',filename)[0], '%Y%m%dT%H%M%S')
                        if self.bdt <= file_dt <= self.edt:
                            files.append(day_file)
        files = sorted(files)
        return files
        
        
#     def open_file(self, filepath):
#         ed = ep.open_raw(raw_file = filepath, sonar_model = 'ek60')
#         return ed
    
    
#     def combine_data(self, filepaths):
#         ed_list = []
#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore") # Catch user warnings that echopype may throw. Usually time related warnings.
#             for fp in filepaths:
#                 _ed = self.open_file(fp)
#                 if len(_ed.platform.channel) == 3:
#                     ed_list.append(_ed)
#                 else:
#                     msg = f"File: {fp} only has {len(_ed.platform.channel)} channels. This file will not be added to the combined dataset."
#                     raise warnings.warn(msg)
#             ed = ep.combine_echodata(ed_list)
#         return ed
    
    
#     def process(self, ed, depth_bin, time_bin):
#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore") # Catch user warnings that echopype may throw. Usually empty slice warnings.
#             sv = ep.calibrate.compute_Sv(ed).compute()
#             mvbs = ep.commongrid.compute_MVBS(sv, range_meter_bin = depth_bin, ping_time_bin = time_bin)
#             mvbs = mvbs.assign_coords(depth = ('echo_range', mvbs['echo_range'].values[::-1]))
#             mvbs = mvbs.swap_dims({'echo_range':'depth'})
#             mvbs = ep.consolidate.swap_dims_channel_frequency(mvbs)
#             mvbs = mvbs.drop_vars(['channel'], errors = 'ignore')
#             mvbs = mvbs.drop_dims(['echo_range'], errors = 'ignore')
#             mvbs = mvbs.rename({'Sv':'sv'})
#         return mvbs

    
#     def get_data(self, depth_bin = 0.2, time_bin = '10s'):
#         ed = self.combine_data(self.files)
#         mvbs = self.process(ed, depth_bin, time_bin)
#         return mvbs
        