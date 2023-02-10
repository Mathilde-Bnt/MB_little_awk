# =====================================================================================
# =========================== Imports required ========================================
# =====================================================================================

import ddensity
import snowtemp
import ddensity_ice

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import csv
import math
from scipy.stats import sem


# ===========================================================================================
# =========================== Pre-processing dataset ========================================
# ===========================================================================================

def fill_in_missing_variables(ds, var_to_copy):
    '''
    Function to artificially add a 'mean' variable in a dataset
    Args:
        ds: clean dataset
        var_to_copy: string, name of variable to be copied into a 'mean' variable
    Returns:
    '''
    ds['mean'] = ds[var_to_copy]


def fill_in_missing_times(ds, day_date, output_dir_name, ds_corresponds_to_file=False):
    '''
    Function that reads a dayly netcdf dataset with variable 'mean' (if it exists), and creates and saves a new netcdf 
    data file with regular timestamps, containing the original data, and nan values where there were no data points originally
    Args:
        ds: day-long dataset from netcdf file, containing a 'mean' variable and the time, x and y coordinates
        day_date: string of the day's date, eg. April 30th 2021 would be '2021-04-30'
        output_dir_name: output directory where the new 'filled-in' netcdf file will be saved
        ds_corresponds_to_file: boolean, is True if ds is the dayly netcdf dataset corresponding to the day_date, default False to avoid errors
    Returns:
    '''
    ds = ds.ffill(dim='time')
    
    # Create artificial regular timestamps every 15 minutes
    timestamps = pd.date_range(start=day_date,end=day_date+'T23:45',freq='15min')
    xstamps = ds.x.values
    ystamps = ds.y.values
    
    # Initialize data
    nan_fill_ins = np.zeros((len(timestamps), len(ystamps), len(xstamps)))
    nan_fill_ins[:][:][:] = np.nan

    # If the netcdf file already exists, keep its data if it matches the timestamps
    if ds_corresponds_to_file:
        for timestamp in ds.time.values:
            minutes_since_start_of_day = (float(timestamp)/1000000000 - pd.to_datetime(day_date).timestamp()) / 60   # in min
            # If the timestamp already exists, keep the data
            if minutes_since_start_of_day%15 == 0:
                index = int(minutes_since_start_of_day // 15)    # index of timestamp in timestamps
                nan_fill_ins[index] = np.array(ds.sel(time=timestamp)['mean'])
        
    # Create dataset
    ds_filled = xr.Dataset({
        'mean': xr.DataArray(
                    data   = nan_fill_ins,
                    dims   = ['time', 'y', 'x'],
                    coords = {'time': timestamps, 'x':xstamps, 'y':ystamps}
                    )
                }
        )
    
    # Name formatting
    split_day = day_date.split('-')
    file_date = split_day[0] + split_day[1] + split_day[2]
    
    # Save data in a new netcdf file
    ds_filled.to_netcdf(output_dir_name+file_date+'.nc')
    return()


def median_space_filtering(ds, min_periods_val, x_span=11, y_span=11):
    '''
    Function to apply median filtering in space
    Args:
        ds: clean data dataset
        x_span: x-span on which to apply filtering, unit [index]
        y_span: y-span on which to apply filtering, unit [index]
    Returns:
    '''
    ds['snow_surface'] = ds['mean'].rolling({'x': x_span, 'y': y_span}, min_periods=min_periods_val, center=True).median()
    print(f'---> Median filtering in space with a window [{x_span}, {y_span}]')
    return()


def median_time_filtering(ds, min_periods_val, time_window=11):
    '''
    Function to apply median filtering in time
    Args:
        ds: clean data dataset
        time_window: time windo on which to apply median filter [index]
    Returns:
    '''
    ds['snow_surface'] = ds['snow_surface'].rolling(time=time_window, min_periods=min_periods_val, center=True).median()
    print(f'---> Median filtering in time with a window of {time_window}')
    return()


def define_summer_surface(ds, start, end):
    '''
    Function to define summer surface by taking the median along time axis a series of summer scans
    Args:
        ds: dataset to use
        start: start date
        end: end date

    Returns:
        add summer surface to dataset

    '''
    ds['summer_surface'] = ds['snow_surface'].sel(time = slice(start, end)).median(dim='time')
    ds.summer_surface.attrs = {'units':'m', 'standard_name':'summer_surface',
                               'long_name':'Summer surface'}
    print(f'---> Summer surface defined based on scans from {start} to {end}')
    return()


# ==========================================================================================
# =========================== Detecting snow events ========================================
# ==========================================================================================

def get_snow_events(ds, x_sel, y_sel, time_window_std, std_threshold):
    '''
    Function that computes the dates of start and end times of snow events (snow accumulation and erosion)
    The events are defined as periods during which the snow-depth rolling standard deviation is higher than 
    a given threshold
    We distinguish between accumulation and erosion of snow by looking at the snow depth before and after the event
    Args:
        ds: clean dataset with 'snow_surface' variable
        x_sel: x coordinate of the point of interest (index)
        y_sel: y coordinate of the point of interest (index)
        time_window_std: size of the rolling window to compute standard deviation
        std_threshold: standard deviation threshold above which the curve is considered to have strong variations > snow event
    Returns:
        start_accumulation_indices: list of indices (in ds) of the times corresponding to the start of an accumulation event
        start_erosion_indices: list of indices (in ds) of the times corresponding to the end of an accumulation event
        end_accumulation_indices: list of indices (in ds) of the times corresponding to the start of an erosion event
        end_erosion_indices: list of indices (in ds) of the times corresponding to the end of an erosion event
    '''
    # Compute standard deviation values around each point
    stdev = ds.isel(x=x_sel, y=y_sel).snow_surface.rolling(time=time_window_std, center=True).std(dim='time').values
    
    # Define snow events' timing
    snow_events_occurrences = stdev > std_threshold   # booleans
    snow_events_occurrences = np.diff(snow_events_occurrences.astype(int))   # 1 or 0
    start_time_indices = np.where(snow_events_occurrences==1)[0]
    end_time_indices = np.where(snow_events_occurrences==-1)[0]
    
    # Initialize variables
    start_accumulation_indices, start_erosion_indices, end_accumulation_indices, end_erosion_indices = [], [], [], []
    
    # Identify snow events as accumulations or erosions
    for index in range(len(start_time_indices)):
        
        start_date = start_time_indices[index]
        end_date = end_time_indices[index]
    
        start_snow_height = float(ds.snow_surface.isel(x=x_sel, y=y_sel, time=start_date))
        end_snow_height = float(ds.snow_surface.isel(x=x_sel, y=y_sel, time=end_date))
    
        if start_snow_height < end_snow_height:
            # Accumulation if the snow height rises during the event
            start_accumulation_indices.append(start_time_indices[index])
            end_accumulation_indices.append(end_time_indices[index])
        else:
            # Erosion if the snow height falls during the event
            start_erosion_indices.append(start_time_indices[index])
            end_erosion_indices.append(end_time_indices[index])
        
    return(start_accumulation_indices, start_erosion_indices, end_accumulation_indices, end_erosion_indices)


def get_change_in_snow_depth(ds, start_events, end_events, index_of_event, x_sel, y_sel):
    '''
    Function to compute snow height difference (absolute value) before and after an event
    Args:
        ds: dataset containing the snow-depth data ('snow_surface' variable)
        start_events: list of time indices at which the events of interest (accumulation or erosion) started
        end_events: list of time indices at which the events of interest (accumulation or erosion) ended
        index_of_event: index of the event of interest in the lists of time indices
        x_sel: index of the x-coordinate of the point of interest
        y_sel: index of the y-coordinate of the point of interest
    Returns:
        difference in snow-depth between the start and end of the event
    '''
    # Get timing of event
    start_date = start_events[index_of_event]
    end_date = end_events[index_of_event]
    
    # Get snow heights at start and end of event
    start_snow_height = float(ds.snow_surface.isel(x=x_sel, y=y_sel, time=start_date))
    end_snow_height = float(ds.snow_surface.isel(x=x_sel, y=y_sel, time=end_date))
    difference = abs(end_snow_height - start_snow_height)
    
    return(difference)


# ======================================================================================================
# =========================== Simulating the snowpack evolution ========================================
# ======================================================================================================

def get_met_forcing(simulation_start_date='2021-12-06 00:00:00', file_start_date='2021-12-01 00:00:00',
                        file_name='/home/mabonnet/Desktop/data/202202_finse_livox/met_obs.csv'):
    '''
    Function to get surface temperature and wind speed forcing from meteorological files
    Args:
        simulation_start_date: date of beginning of simulation, format eg. '2021-04-28 16:30:00' for April 28th 2021, 4:30PM
        file_start_date: date from the meteorological file where we wish the forcing data to start
                (first timestamp in the forcing data), format eg. '2021-04-28 16:30:00' for April 28th 2021, 4:30PM
        file_name: string, path to the meteorological file
    Returns:
        met_temp_data: list containing surface temperatures for each meteorological timestamp (degrees Celcius)
        met_time_data: list containing the meteorological timestamps (s since start of simulation)
        met_wind_data: list containing wind speed for each meteorological timestamp (m.s^-1)
    '''
    print('In get_met_forcing() - Warning: check the format of your file corresponds to the indices given in the functions (wind speed 5, surface temperature 8, time 0).')

    # Initialize lists
    time_series = []
    air_temp_series = []
    surf_temp_series = []
    wind_speed_series = []

    with open(file_name) as met_data:
        reader = csv.reader(met_data)
        next(reader)
        start_date_in_sec = pd.to_datetime(simulation_start_date).timestamp()

        for row in reader:
            # Save timestamps
            time_series.append(pd.to_datetime(row[0]).timestamp() - start_date_in_sec)

            # Save air temperature data
            if row[6]!='':
                air_temp_series.append(float(row[6]))
            else:
                air_temp_series.append(None)

            # Save surface temperature data
            if row[8]!='':
                surf_temp_series.append(float(row[8]))
            else:
                surf_temp_series.append(None)

            # Save wind speed data
            if row[5]!='':
                wind_speed_series.append(float(row[5]))
            else:
                wind_speed_series.append(None)

    # Start of forcing index
    index_start = time_series.index(pd.to_datetime(file_start_date).timestamp() - start_date_in_sec)

    # Only keep the data after the desired start of the forcing
    met_time_data = time_series[index_start:]  # times since start of simulation, in s
    air_temp_series = air_temp_series[index_start:]
    met_temp_data = surf_temp_series[index_start:]
    met_wind_data = wind_speed_series[index_start:]

    return(met_time_data, met_temp_data, met_wind_data)


def is_detected_ice(ds, x_sel, y_sel, end_events, start_events, index_of_end_event, index_of_start_event, dt, slope_threshold, min_duration_in_s):
    '''
    Function that determines whether ice formed during a compaction period (between two distinct snow events), based on the flatness
    of the lidar curve (flat = ice)
    Args:
        ds: clean dataset
        x_sel: index of the x-coordinate of the point of interest
        y_sel: index of the y-coordinate of the point of interest
        end_events: list of time indices at which the events of interest (accumulation or erosion,
                corresponding to the event right before compaction time period) ended
        start_events: list of time indices at which the events of interest (accumulation or erosion,
                corresponding to the event right before compaction time period) started
        index_of_end_event: index (in end_events) of the event right before the compaction time period
        index_of_start_event: index (in start_events) of the event right after the compaction time period
        dt: timestep between two iterations of the simulation, in s
        slope_threshold: value below which the avrage slope is considered negligeable (m.s^-1) > ice
        min_duration_in_s: minimum duration of an ice event (s), under which we consider there is no ice formation
    Returns:
        boolean, True if the compaction period is considered flat enough, False otherwise
    '''
    # Get the indices (in ds) of the start and end dates of the compaction period
    start_date_of_ice = end_events[index_of_end_event]
    end_date_of_ice = start_events[index_of_start_event]
    
    # Get the change in snow height during the period
    start_snow_height = float(ds.snow_surface.isel(x=x_sel, y=y_sel, time=start_date_of_ice))
    end_snow_height = float(ds.snow_surface.isel(x=x_sel, y=y_sel, time=end_date_of_ice))
    snow_level_change = abs(end_snow_height - start_snow_height)
    
    duration_of_ice = (float(ds.time.values[end_date_of_ice]) - float(ds.time.values[start_date_of_ice])) / 1000000000  # in s
    
    # If the compaction period was too short, no time for ice to form
    if duration_of_ice < min_duration_in_s:
        return(False)
    
    # If the slope of the lidar curve during the compaction period is small enough, consider there has been no compaction > ice
    else:
        slope_of_curve_during_ice = snow_level_change / duration_of_ice    # absolute value

        if slope_of_curve_during_ice < slope_threshold:
            return(True)
        else:
            return(False)


def simulate_snowpack_evolution(ds, x_sel, y_sel, nb_iterations, end_accumulation_times, end_erosion_times,
                                start_accumulation, end_accumulation, start_erosion, end_erosion,
                                jj, dt, ro_layer, ro_water, ro_ice, t_old, tf, tsfc_default, dy_snow, age_layers, gamma, cp_snow,
                                melt_flag, a1, a2, met_temp_data=[None], met_wind_data=[None], met_time_data=[0],
                                fit_top_of_snowfall_to_curve=False, erode_several_layers=False, detect_ice=False,
                                start_accumulation_times=None, start_erosion_times=None, a1_vector=None, slope_threshold=None,
                                min_duration_in_s=None):
    '''
    Function that simulates the evolution of the snowpack over a certain period of time: at each timestamp, the new density, height and
    temperature of each layer is computed, according to SnowModel (2006)
    Options: use correct meteorological forcing, fit the top of snowfalls to the lidar data curve, erode more than just the top layer,
    detect ice layers based on the flatness of the lidar curve (and stop their compaction)
    Args:
        ds: clean dataset
        x_sel: x-coordinate of the point of interest (index)
        y_sel: y-coordinate of the point of interest (index)
        nb_iterations: number of iterations
        
        end_accumulation_times: list of ending times of accumulations in seconds since data starting date
        end_erosion_times: list of ending times of erosions in seconds since data starting date
        start_accumulation: list of the indices of starting times of accumulations in ds
        end_accumulation: list of the indices of ending times of accumulations in ds
        start_erosion: list of the indices of starting times of erosions in ds
        end_erosion: list of the indices of ending times of erosions in ds
        
        jj: number of layers initially present
        dt: timestep (s)
        ro_layer(1*max_nb_of_layers) array containing density value (kg.m^-3) for each layer
        ro_water: density of water (kg.m^-3)
        ro_ice: density of ice (kg.m^-3)
        t_old: (1*max_nb_of_layers) array containing temperature value (degrees Celcius) for each layer
        tf: ice fusion temperature (degrees Celcius)
        tsfc_default: surface temperature (degrees Celcius)
        dy_snow: (1*max_nb_of_layers) array containing depth value (m) for each layer
        age_layers: (1*max_nb_of_layers) array containing age (s) of each layer
        gamma: (1*max_nb_of_layers) array containing zeros
        cp_snow: thermal capacity of snow (J.kg^-1.K^-1)
        melt_flag: (1*max_nb_of_layers) array containing melt value (1 or 0) for each layer
        a1, a2: exponential parameters, empirically calibrated, a1 in m^-1.s^-1, a2 in m^3.kg^-1

        met_temp_data: array containing surface temperatures for each meteorological timestamp (degrees Celcius),
                default [None], i.e. use tsfc=tsfc_default
        met_wind_data: array containing wind speed for each meteorological timestamp (m.s^-1), default [None], i.e. use new_snow_ro=150
        met_time_data: array containing the meteorological timestamps (s since start of simulation), default [0]

        fit_top_of_snowfall_to_curve: boolean, if True the height of snowfalls will be such that the snow depth
                matches the lidar curve depth at the end of snowfalls, default False
        erode_several_layers: boolean, if True the erosions may affect more than just the top layer
                (remove the measured snow height difference from the simulated layers), default False

        detect_ice: boolean, if True, flat portions of the lidar curve are detected as ice and stop being compacted, default False
        start_accumulation_times:  list of starting times of accumulations in seconds since data starting date,
                useful for the ice detection option, default None
        start_erosion_times:  list of starting times of erosions in seconds since data starting date,
                useful for the ice detection option, default None
        a1_vector: array containing the a1 values (m^-1.s^-1) for each layers (null for ice layers that do not compact),
                useful for the ice detection option, default None
        slope_threshold: value below which the slope of a compaction period is considered negligeable (m.s^-1),
                useful for the ice detection option, default None
        min_duration_in_s: minimum duration of an ice event (s), under which we consider there is no ice formation,
                useful for the ice detection option, default None

    Returns:
        ro_layer_evolution: list of the states of layers' density (kg.m^-3) through time,
                format [[ro_layer_1, ro_layer_2, ro_layer_3, ...]_time_1, [...]_time_2, [...], ...]
        depth_evolution: list of the states of layers' depth (m) through time,
                format [[depth_layer_1, depth_layer_2, depth_layer_3, ...]_time_1, [...]_time_2, [...], ...]
        temperature_evolution: list of the states of layers' temperature (degrees Celcius) through time,
                format [[temp_layer_1, temp_layer_2, temp_layer_3, ...]_time_1, [...]_time_2, [...], ...]
        ice_layers_times: list of the simulation times' indices at which ice was detected
    '''
    # Check all parameters are well defined for ice detection
    if detect_ice:
        if start_accumulation_times==None or start_erosion_times==None or a1_vector==None or slope_threshold==None or min_duration_in_s==None:
            print('Values are missing for ice detection')
            return()
        ice_layers_times = []

    # Initialize arrays to keep track of variables in time
    ro_layer_evolution, depth_evolution, temperature_evolution = [], [], []

    # Initialize surface temperature
    if met_temp_data[0] != None:
        tsfc = float(met_temp_data[0])
    else:
        tsfc = tsfc_default
    
    # Initialize new snow density
    if met_wind_data[0] != None and float(met_wind_data[0]) > 6:
        new_snow_ro = 250
    else:
        new_snow_ro = 150

    # Initialize indices of next accumulation/erosion events coming up (updated when their time is past)
    accumulation_index, erosion_index, temperature_index, wind_index = 0, 0, 0, 0

    # Initialize ice detection marker
    if detect_ice:
        is_iced = False
    
    # Run model
    for i in range(nb_iterations):

        # Update surface temperature if needed
        if temperature_index<len(met_temp_data) and met_time_data[temperature_index]!=None and i*dt>=met_time_data[temperature_index]:
            if met_temp_data[temperature_index] != None:
                tsfc = float(met_temp_data[temperature_index])
            temperature_index += 1
        
        # Update new snow density if needed
        if wind_index<len(met_wind_data) and met_time_data[wind_index]!=None and i*dt>=met_time_data[wind_index]:
            if met_wind_data[wind_index] != None:
                if float(met_wind_data[wind_index]) > 6:
                    new_snow_ro = 250
            else:
                new_snow_ro = 150
            wind_index += 1

        # Detection of accumulations
        if accumulation_index<len(end_accumulation_times) and i*dt>=end_accumulation_times[accumulation_index]:     # if an accumulation was passed

            if detect_ice:
                if accumulation_index+1 >= len(start_accumulation_times) and erosion_index >= len(start_erosion_times):   # got to the end
                    if is_detected_ice(ds, x_sel, y_sel, end_accumulation, range(len(ds.time.values)), accumulation_index, -1, dt, slope_threshold, min_duration_in_s):
                        is_iced=True
                elif accumulation_index+1 < len(start_accumulation_times) and (erosion_index>=len(start_erosion_times)
                            or start_accumulation_times[accumulation_index+1]<start_erosion_times[erosion_index]):  # next event will be an accumulation
                    if is_detected_ice(ds, x_sel, y_sel, end_accumulation, start_accumulation, accumulation_index, accumulation_index+1, dt, slope_threshold, min_duration_in_s):
                        is_iced=True
                elif erosion_index < len(start_erosion_times) and (accumulation_index+1>=len(start_accumulation_times)
                            or start_accumulation_times[accumulation_index+1]>start_erosion_times[erosion_index]):  # next event will be an erosion
                    if is_detected_ice(ds, x_sel, y_sel, end_accumulation, start_erosion, accumulation_index, erosion_index, dt, slope_threshold, min_duration_in_s):
                        is_iced=True
                    
                if is_iced:
                    ice_layers_times.append(i)
                    for i in range(jj+1):
                        a1_vector[i] = 0
                else:
                    a1_vector[jj] = a1
                is_iced = False
            
            if fit_top_of_snowfall_to_curve:
                snow_depth_total = sum(dy_snow[i] for i in range(0,jj))
                ddepth = float(ds.snow_surface.isel(x=x_sel, y=y_sel, time=end_accumulation[accumulation_index])) - snow_depth_total
            else:
                ddepth = get_change_in_snow_depth(ds, start_accumulation, end_accumulation, accumulation_index, x_sel, y_sel)
            ro_layer[jj] = new_snow_ro
            t_old[jj] = tsfc
            dy_snow[jj] = ddepth
            age_layers[jj] = (nb_iterations-i) * dt     # age in seconds at end of simulation
            if t_old[jj] <= 0:
                melt_flag[jj] = 0
            else:
                melt_flag[jj] = 1
            jj += 1
            accumulation_index += 1     # next accumulation that will come up
    
        # Detection of erosions
        if erosion_index<len(end_erosion_times) and i*dt>=end_erosion_times[erosion_index]:     # if an erosion was passed

            if detect_ice:
                if accumulation_index >= len(start_accumulation_times) and erosion_index+1 >= len(start_erosion_times):   # got to the end
                    if is_detected_ice(ds, x_sel, y_sel, end_erosion, range(len(ds.time.values)), erosion_index, -1, dt, slope_threshold, min_duration_in_s):
                        is_iced=True
                elif accumulation_index < len(start_accumulation_times) and (erosion_index+1>=len(start_erosion_times)
                            or start_accumulation_times[accumulation_index]<start_erosion_times[erosion_index+1]):   # next event will be an accumulation
                    if is_detected_ice(ds, x_sel, y_sel, end_erosion, start_accumulation, erosion_index, accumulation_index, dt, slope_threshold, min_duration_in_s):
                        is_iced=True
                elif erosion_index+1 < len(start_erosion_times) and (accumulation_index>=len(start_accumulation_times)
                            or start_accumulation_times[accumulation_index]>start_erosion_times[erosion_index+1]):   # next event will be an erosion
                    if is_detected_ice(ds, x_sel, y_sel, end_erosion, start_erosion, erosion_index, erosion_index+1, dt, slope_threshold, min_duration_in_s):
                        is_iced=True
                    
                if is_iced:
                    ice_layers_times.append(i)
                    for i in range(jj+1):
                        a1_vector[i] = 0
                else:
                    a1_vector[jj] = a1
                is_iced = False

            ddepth = get_change_in_snow_depth(ds, start_erosion, end_erosion, erosion_index, x_sel, y_sel)
            erosion_index += 1     # next erosion that will come up
            if erode_several_layers:
                # Get rid of all the layers that are smaller than the measured loss of snow height
                while jj>0 and dy_snow[jj-1] <= ddepth:
                    ddepth -= dy_snow[jj-1]
                    jj -= 1
                    dy_snow[jj] = 0
                    ro_layer[jj] = 0
                    t_old[jj] = 0
                    age_layers[jj] = 0
                    melt_flag[jj] = 0
            if jj>0:
                if dy_snow[jj-1] > ddepth:
                    dy_snow[jj-1] = dy_snow[jj-1] - ddepth
                else:
                    jj -= 1
                    dy_snow[jj] = 0
                    ro_layer[jj] = 0
                    t_old[jj] = 0
                    age_layers[jj] = 0
                    melt_flag[jj] = 0
    
        # Update layers' parameters (fortran code model)
        if detect_ice:
            ro_layer, dy_snow = ddensity_ice.ddensity_ml(ro_layer, tf, dt, ro_water, ro_ice, t_old, jj, dy_snow, a1_vector, a2)
        else:
            ro_layer, dy_snow = ddensity.ddensity_ml(ro_layer, tf, dt, ro_water, ro_ice, t_old, jj, dy_snow, a1, a2)
        t_old = snowtemp.snowtemp_ml(gamma, t_old, tsfc, jj, dt, ro_layer, cp_snow, tf, dy_snow, melt_flag)
    
        # Keep track of events
        ro_layer_evolution.append(ro_layer)
        depth_evolution.append(dy_snow)
        temperature_evolution.append(t_old)
        
        if not detect_ice:
            ice_layers_times = None
        
    return(ro_layer_evolution, depth_evolution, temperature_evolution, ice_layers_times)


# ====================================================================================================
# =========================== Plotting the snowpack evolution ========================================
# ====================================================================================================

def plot_simul_and_signal(ds, x_sel, y_sel, depth_evolution, nb_layers_to_plot, data_start_date, dt, nb_iterations,
                          start_accumulation, end_accumulation, start_erosion, end_erosion, ice_layers_times=None,
                          my_title='Comparison between lidar-measured and simulated snow depth', save_file=False,
                          my_file_name='my_fig.png', my_figsize=(16, 7)):
    '''
    Function to plot the simulated snowpack layers and lidar signal on the same plot
    Args:
        ds: clean dataset
        x_sel: x-coordinate of the point of interest
        y_sel: y-coordinate of the point of interest
        depth_evolution: list of the states of layers' depth (m) through time, 
                format [[depth_layer_1, depth_layer_2, depth_layer_3, ...]_time_1, [...]_time_2, [...], ...]
        nb_layers_to_plot: number of layers to plot, starting from the bottom
        data_start_date: first date of the dataset, pandas datetime format
        dt: timestep used in snowpack simulation (s)
        nb_iterations: number of iterations used in snowpack simulation
        
        start_accumulation: list of the indices (in ds) of starting times of accumulations
        end_accumulation: list of the indices (in ds) of ending times of accumulations
        start_erosion: list of the indices (in ds) of starting times of erosions
        end_erosion: list of the indices (in ds) of ending times of erosions

        ice_layers_times: list of time indices where an ice layer was detected, default None
        
        my_title: title of the figure, default 'Comparison between lidar-measured and simulated snow depth'
        save_file: boolean, default False
        my_file_name: name to be given to the saved file, default 'my_fig.png'
        my_figsize: figure size, default (16, 7)
    Returns:
    '''
    # Construct the data to plot: height of each layer as a function of time
    layers = np.zeros((nb_layers_to_plot, len(depth_evolution)))
    # First layer
    for i in range(len(depth_evolution)):
        layers[0][i] = depth_evolution[i][0]
    # Next layers
    for layer_index in range(1, nb_layers_to_plot):
        for i in range(len(depth_evolution)):
            layers[layer_index][i] = depth_evolution[i][layer_index] + layers[layer_index-1][i]
    
    # Define figure and timestamps
    fig = plt.figure(figsize=my_figsize)
    times = pd.date_range(start=data_start_date,freq=str(dt)+'S',periods=nb_iterations)
    
    # Plot each layer
    for layer_index in range(nb_layers_to_plot):
        plt.plot(times, layers[layer_index], label='layer '+str(layer_index+1))

    # Plot a star at each ice layer detection
    if ice_layers_times != None:
        for time_index in ice_layers_times:
            plt.plot(times[time_index], layers[-1][time_index]+0.001, c='y', marker='*', markersize=15, label='ice layer detected')
    
    # Plot the lidar signal
    ds.isel(x=x_sel, y=y_sel).snow_surface.plot(c='k', alpha=0.2, label='lidar signal')

    # Plot the start and end of detected snow events on the lidar curve
    ds.isel(x=x_sel, y=y_sel, time=start_accumulation).snow_surface.plot(c='b', marker='^', markersize=6, linestyle='None', label='start accum.')
    ds.isel(x=x_sel, y=y_sel, time=end_accumulation).snow_surface.plot(c='g', marker='^', markersize=6, linestyle='None', label='end accum.')
    ds.isel(x=x_sel, y=y_sel, time=start_erosion).snow_surface.plot(c='m', marker='v', markersize=6, linestyle='None', label='start erosion')
    ds.isel(x=x_sel, y=y_sel, time=end_erosion).snow_surface.plot(c='r', marker='v', markersize=6, linestyle='None', label='end erosion')
    
    plt.ylabel('snow depth (m)')
    plt.xlabel('time (date)')

    plt.legend()
    plt.title(my_title)
    
    if save_file:
        fig.savefig(my_file_name)
    
    plt.show(block=False)

    return()


# ====================================================================================================
# =========================== Measuring closeness of simulations =====================================
# ====================================================================================================

def rmse_measure(simul_total_height_array, lidar_height_array):
    '''
    Function that computes the root mean square error between two arrays
    Args:
        simul_total_height_array: array containing the total height (m) of the simulated snowpack at each timestamp
        lidar_height_array: array containing the height (m) of the snowpack measured by the lidar at each timestamp
    Returns:
        rmse: root mean square error between the two arrays
    '''
    mse = np.square(np.subtract(simul_total_height_array, lidar_height_array)).mean()
    rmse = math.sqrt(mse)
    
    return(rmse)


def stderr_measure(simul_total_height_array, lidar_height_array):
    '''
    Function that computes the standard error of the difference of two arrays (measures to what extent the two series
    are "parallel", the smaller the better)
    Args:
        simul_total_height_array: array containing the total height (m) of the simulated snowpack at each timestamp
        lidar_height_array: array containing the height (m) of the snowpack measured by the lidar at each timestamp
    Returns:
        stderr: standard error of the difference of the two arrays
    '''
    difference = simul_total_height_array - lidar_height_array
    stderr = sem(difference)
    
    return(stderr)


def p_correl_measure(simul_total_height_array, lidar_height_array):
    '''
    Function that computes the Pearson correlation between two arrays
    Args:
        simul_total_height_array: array containing the total height (m) of the simulated snowpack at each timestamp
        lidar_height_array: array containing the height (m) of the snowpack measured by the lidar at each timestamp
    Returns:
        p_correl: Pearson correlation between the two arrays
    '''
    dataframe = {
        "Array_1": simul_total_height_array,
        "Array_2": lidar_height_array
    }

    data = pd.DataFrame(dataframe)
    p_correl = data.corr().iloc[0, 1]
    
    return(p_correl)


def all_measures(a1, a2, simul_total_height_array, lidar_height_array):
    '''
    Function that computes the rmse, standard error of the difference and Pearson correlation between two arrays
    Args:
        a1: a1 value used to obtain simul_total_height_array
        a2: a2 value used to obtain simul_total_height_array
        simul_total_height_array: array containing the total height of the simulated snowpack at each timestamp
        lidar_height_array: array containing the height of the snowpack measured by the lidar at each timestamp
    Returns:
        a1, a2: unchanged values, used for identification of the results
        all_measrs: tuple containing the computed values in the order (rmse, stde, p_corr)
    '''
    all_measrs = rmse_measure(simul_total_height_array, lidar_height_array), stderr_measure(simul_total_height_array, lidar_height_array), p_correl_measure(simul_total_height_array, lidar_height_array)
    
    return(a1, a2, all_measrs)


# ====================================================================================================
# =========================== Computing physical profiles of simulations =============================
# ====================================================================================================

def get_data_from_ref(index_of_ref_layer, ro_layer_array, depth_array):
    '''
    Function that computes some characteristic values (SWE, height, average density) above a "ref" interface
    Args:
        index_of_ref_layer: index of the layer that is directly underneath the ref interface
        ro_layer_array: array of the densities of each layers at the end of the simulation
        depth_array: array of the heights of each layers at the end of the simulation
    Returns:
        swe_from_ref: value of the SWE above the given interface
        height_from_ref: height of the snowpack above the given interface
        ave_density_from_ref: averaged density of the snow above the given interface
    '''
    swe_from_ref = np.dot(np.array(ro_layer_array[index_of_ref_layer+1:]), np.array(depth_array[index_of_ref_layer+1:])) / 1000
    height_from_ref = sum(depth_array[i] for i in range(index_of_ref_layer+1, len(depth_array)))
    ave_density_from_ref = np.dot(np.array(ro_layer_array[index_of_ref_layer+1:]), np.array(depth_array[index_of_ref_layer+1:])) / height_from_ref
    
    return(swe_from_ref, height_from_ref, ave_density_from_ref)


def get_depth_layers_indices(bottom_depth, sampling_length, depth_array):
    '''
    Function that creates an array of the layer index at each (regularly) sampled depth (useful for computing profiles)
    Args:
        bottom_depth: depth of the lowest sample, in meters, counted positively from the surface of the snowpack downwards (included)
        sampling_length: height between two consecutive sampled depths
        depth_array: array of the heights of each layers, in meters, at the end of the simulation
    Returns:
        layers_array: array containing the indices of the layers corresponding to each sampled depth (which layer was sampled)
    '''
    current_layer_index = next(len(depth_array) - i for i, j in enumerate(reversed(depth_array), 1) if j != 0)  # index of top layer
    remaining_depth_in_current_layer = depth_array[current_layer_index]
    current_depth = 0
    layers_array = [current_layer_index]    # top of the snowpack
    leftover_depth = sampling_length
    
    while current_depth < bottom_depth and current_layer_index > -1:
        while remaining_depth_in_current_layer >= sampling_length and current_depth < bottom_depth:   # next sample is in current layer
            layers_array.append(current_layer_index)
            remaining_depth_in_current_layer -= sampling_length
            current_depth += sampling_length
        while remaining_depth_in_current_layer < leftover_depth and current_depth < bottom_depth and current_layer_index > -1:   # next sample is not in current layer
            current_layer_index -= 1   # go to next layer
            leftover_depth -= remaining_depth_in_current_layer   # what is left to "subtract" to this new layer
            remaining_depth_in_current_layer = depth_array[current_layer_index]
            if remaining_depth_in_current_layer >= leftover_depth:   # next sample is in current layer
                layers_array.append(current_layer_index)
                remaining_depth_in_current_layer -= leftover_depth
                current_depth += sampling_length
                leftover_depth = sampling_length
                
    layers_array = list(reversed(layers_array))
    
    return(layers_array)


def density_profile(bottom_depth, sampling_length, depth_array, ro_array):
    '''
    Function that creates an array of the layer density at each (regularly) sampled depth
    Args:
        bottom_depth: depth of the lowest sample, in meters, counted positively from the surface of the snowpack downwards (included)
        sampling_length: height between two consecutive sampled depths
        depth_array: array of the heights of each layers, in meters, at the end of the simulation
        ro_array: array of the densities of each layers, in kg.m^-3, at the end of the simulation
    Returns:
        ro_profile: array containing the simulated densities at each sampled depth
    '''
    indices_array = get_depth_layers_indices(bottom_depth, sampling_length, depth_array)
    ro_profile = [ro_array[i] for i in indices_array]
    
    return(ro_profile)


def temp_profile(bottom_depth, sampling_length, depth_array, temp_array):
    '''
    Function that creates an array of the layer temperature at each (regularly) sampled depth
    Args:
        bottom_depth: depth of the lowest sample, in meters, counted positively from the surface of the snowpack downwards (included)
        sampling_length: height between two consecutive sampled depths
        depth_array: array of the heights of each layers, in meters, at the end of the simulation
        temp_array: array of the temperatures of each layers at the end of the simulation
    Returns:
        temp_profile: array containing the simulated temperatures at each sampled depth
    '''
    indices_array = get_depth_layers_indices(bottom_depth, sampling_length, depth_array)
    temp_profile = [temp_array[i] for i in indices_array]
    
    return(temp_profile)
