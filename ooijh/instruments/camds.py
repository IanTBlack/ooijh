"""IN DEVELOPMENT"""

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
        

    def open_image(self, filepath):
        im = cv2.imread(file)
        return im
    
