import xarray as xr


DROP_METHODS= ['Available']


DROP_STREAMS = ['data_header','system_data','config','engineering','metadata','power','blank','dark','diagnostics',
                'dspec','fcoeff','log9','m_wvs','cspp_eng','calibration','coefficients','secondary_node','status',
                'hardware','voltage','station','winch_eng','record_cal','dev1','primary_node','clock_data',
                'identification_string','settings','hyd_o','msg_counts','superv','error','cg_dcl','cg_stc_eng','wfp_eng']


DROP_SYS_VARS = ['driver_timestamp','id','internal_timestamp','port_timestamp','provenance','preferred_timestamp',
                 'ingestion_timestamp', 'suspect_timestamp','dcl_controller_timestamp','profiler_timestamp']


DROP_CTD_VARS = ['ctd_time', 'pressure_temp', 'conductivity','temperature', 'date_time_string','pressure',
                 'oxy_calphase','ctd_tc_oxygen','oxygen_from_calphase','do_stable_sample-corrected_dissolved_oxygen',
                 'corrected_dissolved_oxygen','oxygen','do_stable_sample-oxygen_from_calphase','oxy_temp']


DROP_PHSEN_VARS = ['record_type','voltage_battery','thermistor_start','thermistor_end','light_measurements',
                   'record_time','reference_light_measurments', 'reference_light_measurements_dim_0', 'spectrum',
                   'signal_intensity_578_dim_0', 'ph_light_measurements_dim_0','signal_intensity_434_dim_0',
                   'ph_light_measurements','signal_intensity_578','phsen_battery_volts','signal_intensity_434',
                   'record_length','checksum','reference_light_measurements','unique_id','sea_water_practical_salinity','phsen_abcdef_signal_intensity_578_dim_0', 'phsen_abcdef_signal_intensity_434_dim_0','reference_light_measurements_dim_0', 'phsen_abcdef_signal_intensity_434', 'phsen_abcdef_signal_intensity_578','passed_checksum']


DROP_VEL3D_VARS =['ensemble_counter','correlation_beam_1','correlation_beam_2','correlation_beam_3',
                  'turbulent_velocity_north','turbulent_velocity_east','turbulent_velocity_vertical',
                  'analog_input_1','analog_input_2','amplitude_beam_1','amplitude_beam_2','amplitude_beam_3']


DROP_PCO2_VARS = ['spectrum','record_type','absorbance_ratio_434','absorbance_ratio_620','voltage_battery',
                  'light_measurements','pco2w_b_absorbance_blank_434','record_length','pco2w_b_absorbance_blank_620',
                  'record_time', 'checksum','thermistor_raw', 'unique_id','date_time_string','humidity', 'humidity_temperature','irga_detector_temperature','avg_irga_temperature','current_a2d','zero_a2d','supply_voltage''irga_source_temperature','eastward_velocity','gas_stream_pressure','longwave_irradiance','northward_velocity','precipitation','supply_voltage','shortwave_irradiance']


DROP_ADCP_VARS = ['corrected_echo_intensity_beam4','ensemble_number', 'corrected_echo_intensity_beam3',
                  'corrected_echo_intensity_beam2', 'error_seawater_velocity', 'corrected_echo_intensity_beam1',
                  'pitch', 'velocity_beam4', 'velocity_beam1', 'velocity_beam3','velocity_beam2','percent_good_beam1',
                  'percent_good_beam3','percent_good_beam2','percent_good_beam4', 'heading',
                  'sysconfig_vertical_orientation','roll','ctdbp_no_sample-depth','pressure','temperature',
                  'int_ctd_pressure', 'echo_intensity_beam2','salinity','correlation_magnitude_beam4',
                  'correlation_magnitude_beam2','correlation_magnitude_beam3','correlation_magnitude_beam1',
                  'echo_intensity_beam4','echo_intensity_beam3','echo_intensity_beam1','depth', 'num_cells',
                  'depth_from_pressure','non_zero_pressure','non_zero_depth','bin_1_distance',
                  'percent_transforms_reject','percent_good_4beam','percent_good_3beam',
                  'ctdbp_cdef_dcl_instrument_recovered-depth','error_velocity','ctdbp_cdef_dcl_instrument-depth',
                  'water_velocity_up','water_velocity_east','water_velocity_north','percent_bad_beams']


DROP_OPTAA_VARS = ['beam_attenuation','optical_absorption','meter_type','record_length','num_wavelengths',
                   'packet_type','pressure_counts','checksum','serial_number']


DROP_FLORT_VARS = ['raw_signal_cdom','measurement_wavelength_chl','measurement_wavelength_cdom',
                   'measurement_wavelength_beta','raw_internal_temp','raw_signal_beta','raw_signal_chl',
                   'seawater_scattering_coefficient','total_volume_scattering_coefficient',
                   'sea_water_practical_salinity','sea_water_temperature']


DROP_METBK_VARS = ['met_windavg_mag_corr_east','met_windavg_mag_corr_north','barometric_pressure']


DROP_DOSTA_VARS = ['blue_amplitude','blue_phase','calibrated_phase','dosta_abcdjm_cspp_tc_oxygen','int_ctd_pressure',
                   'product_number','red_amplitude','red_phase','serial_number','temp_compensated_phase',
                   'sea_water_practical_salinity','sea_water_temperature','raw_temperature',
                   'estimated_oxygen_concentration','estimated_oxygen_saturation']


DROP_NUTNR_VARS = ['aux_fitting_1', 'aux_fitting_2', 'checksum','date_of_sample', 'humidity', 'lamp_time', 
                   'nutnr_absorbance_at_254_nm','nutnr_absorbance_at_350_nm','nutnr_bromide_trace',
                   'nutnr_current_main', 'nutnr_dark_value_used_for_fit','nutnr_fit_base_1','nutnr_fit_base_2',
                   'nutnr_integration_time_factor','nutnr_nitrogen_in_nitrate','nutnr_spectrum_average',
                   'nutnr_voltage_int', 'sea_water_practical_salinity', 'sea_water_temperature','serial_number',
                   'spectral_channels','temp_interior','temp_lamp','temp_spectrometer','time_of_sample','voltage_lamp',
                   'voltage_main','wavelength','nitrate_concentration','ctd_time_uint32','day_of_year','year','ctdpf_j_cspp_instrument_recovered-sea_water_practical_salinity','ctdpf_j_cspp_instrument_recovered-sea_water_temperature']


DROP_VELPT_VARS = ['amplitude_beam1','amplitude_beam2','amplitude_beam3','analog1','battery_voltage_dv',
                   'date_time_string', 'error_code', 'heading_decidegree', 'pitch_decidegree', 'roll_decidegree',
                   'sea_water_pressure_mbar', 'sound_speed_dms','status','temperature_centidegree','velocity_beam1',
                   'velocity_beam2', 'velocity_beam3']


DROP_SPKIR_VARS = ['va_sense','vin_sense','frame_counter','instrument_id','passed_checksum','internal_temperature',
                   'sample_delay','serial_number','timer','channel_array']


DROP_VEL3D_VARS = ['amplitude_beam_1','amplitude_beam_2','amplitude_beam_3','analog_input_1','analog_input_2',
                   'correlation_beam_1','correlation_beam_2','correlation_beam_3','ensemble_counter',
                   'sea_water_pressure_mbar','turbulent_velocity_north','turbulent_velocity_east',
                   'turbulent_velocity_vertical']

DROP_FDCHP_VARS = ['v_num_datacollection','time_datacollection','pitch_min','pitch_max','u_corr_std','z_accel_std',
                   'wind_v_max','wind_v_min','w_corr_std','x_ang_rate_std','y_ang_rate_min','v_corr_std','heading_std','status_datacollection','roll_min','roll_max','heading_max','speed_of_sound_max','wind_w_std','z_ang_rate_min','x_ag_rate_min','roll_std','x_accel_std','fdchp_status_1','fdchp_status_2','day','millisecond','month','minute','heading_min','second','pitch_std','year','y_ang_rate_std','z_ang_rate_std','hour','instrument_start_timestamp','speed_of_sound_min','speed_of_sound_std','wind_u_min','wind_u_max','wind_u_std','wind_v_min','wind_v_max','wind_v_std','x_accel_max','x_accel_min','x_ang_rate_max','x_ang_rate_min','y_accel_max','y_accel_min','y_accel_std','z_accel_min','z_accel_max']

DROP_WAVSS_VARS = ['time_string','serial_number','date_string']

DROP_MOPAK_VARS = []

DROP_PARAD_VARS = ['date_string', 'parad_j_par_counts_output','time_string']

def drop_qc_test_vars(ds: xr.Dataset) -> xr.Dataset:
    """
    Drop qc_executed and qc_results variables.
    
    :param ds: The input xr.Dataset.
    :return: A slightly cleaner xr.Dataset.
    """
    
    qc_vars = [v for v in ds.data_vars if 'qc_' in v]
    ds = ds.drop_vars(qc_vars, errors = 'ignore')
    return ds


def drop_qartod_test_vars(ds: xr.Dataset) -> xr.Dataset:
    """
    Drop qartod_executed and qartod_results variables.
    
    :param ds: The input xr.Dataset.
    :return: A slightly cleaner xr.Dataset.
    """
    
    qartod_vars = [v for v in ds.data_vars if 'qartod_' in v]
    ds = ds.drop_vars(qartod_vars, errors = 'ignore')
    return ds
