# =====================================================================================
# =========================== Imports required ========================================
# =====================================================================================

import ddensity
import snowtemp

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import csv


# ===========================================================================================
# =========================== Pre-processing dataset ========================================
# ===========================================================================================

def fill_in_missing_variables(ds, var_to_copy):
    '''
    Function to artificially add 'snow_surface' variable in dataset
    Args:
        ds: clean data dataset
        var_to_copy: string, name of variable to be copied into a 'mean' variable
    Returns:
    '''
    ds['mean'] = ds[var_to_copy]


def fill_in_missing_times(ds, day_date, output_dir_name, ds_corresponds_to_file=False):
    '''
    Function that reads a dayly netcdf dataset with variable 'mean' (if it exists), and creates and saves a new netcdf data file with regular timestamps,
    containing the original data, and nan values where there were no data points originally
    Args:
        ds: day-long dataset from netcdf file, containing a 'mean' variable and the x and y coordinates in the desired format
        day_date: string of the day's date, eg. April 30th 2021 would be '2021-04-30'
        output_dir_name: output directory where the new 'filled-in' netcdf file will be saved
        ds_corresponds_to_file: boolean, is True of ds is the dayly netcdf dataset corresponding to the day_date, default False to avoid errors
    Returns:
    '''
    ds = ds.ffill(dim='time')
    
    timestamps = pd.date_range(start=day_date,end=day_date+'T23:45',freq='15min')   # TODO modular
    xstamps = ds.x.values
    ystamps = ds.y.values
    
    nan_fill_ins = np.zeros((len(timestamps), 91, 201))   # TODO modular
    nan_fill_ins[:][:][:] = np.nan

    if ds_corresponds_to_file:
        for timestamp in ds.time.values:
            minutes_since_start_of_day = (float(timestamp)/1000000000 - pd.to_datetime(day_date).timestamp()) / 60   # in min
            if minutes_since_start_of_day%15 == 0:
                index = int(minutes_since_start_of_day // 15)    # index of timestamp in timestamps
                nan_fill_ins[index] = np.array(ds.sel(time=timestamp)['mean'])
        
    # create dataset
    ds_filled = xr.Dataset({
        'mean': xr.DataArray(
                    data   = nan_fill_ins,
                    dims   = ['time', 'y', 'x'],
                    coords = {'time': timestamps, 'x':xstamps, 'y':ystamps}
                    )
                }
        )
    
    split_day = day_date.split('-')
    file_date = split_day[0] + split_day[1] + split_day[2]
    
    ds_filled.to_netcdf(output_dir_name+file_date+'.nc')   # TODO use paths instead
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
    
    # note that a skipna argument is available (useful?) in xarray.DataArray.median()

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
    
    # note that a skipna argument is available (useful?) in xarray.DataArray.median()


# ==========================================================================================
# =========================== Detecting snow events ========================================
# ==========================================================================================

def get_snow_events(ds, x_sel, y_sel, time_window_std, std_threshold):
    '''
    Function that computes the dates of start and end times of snow events (snow accumulation and erosion)
    The events are defined as periods during which the snow-depth rolling standard deviation is higher than 
    a given threshold (set to 0.02 here)
    We distinguish between accumulation and erosion of snow by looking at the snow depth before and after the event
    Args:
        ds: clean dataset with 'snow_surface' variable
        x_sel: x coordinate of the point of interest (index)
        y_sel: y coordinate of the point of interest (index)
        time_window_std: size of the rolling window to compute standard deviation
        std_threshold: standard deviation threshold above which the curve is considered to have strong variations > snow event
    Returns:
        start_accumulation_indices: list of time indices corresponding to the start of an accumulation event
        start_erosion_indices: list of time indices corresponding to the end of an accumulation event
        end_accumulation_indices: list of time indices corresponding to the start of an erosion event
        end_erosion_indices: list of time indices corresponding to the end of an erosion event
    '''
    
    stdev = ds.isel(x=x_sel, y=y_sel).snow_surface.rolling(time=time_window_std, center=True).std(dim='time').values
    
    snow_events_occurrences = stdev > std_threshold   # booleans
    snow_events_occurrences = np.diff(snow_events_occurrences.astype(int))   # 1 or 0
    start_time_indices = np.where(snow_events_occurrences==1)[0]
    end_time_indices = np.where(snow_events_occurrences==-1)[0]
    
    start_accumulation_indices = []
    start_erosion_indices = []
    end_accumulation_indices = []
    end_erosion_indices = []
    
    for index in range(len(start_time_indices)):
        
        start_date = start_time_indices[index]
        end_date = end_time_indices[index]
    
        start_snow_height = float(ds.snow_surface.isel(x=x_sel, y=y_sel, time=start_date))
        end_snow_height = float(ds.snow_surface.isel(x=x_sel, y=y_sel, time=end_date))
    
        if start_snow_height < end_snow_height:
            start_accumulation_indices.append(start_time_indices[index])
            end_accumulation_indices.append(end_time_indices[index])
        else:
            start_erosion_indices.append(start_time_indices[index])
            end_erosion_indices.append(end_time_indices[index])
        
    return(start_accumulation_indices, start_erosion_indices, end_accumulation_indices, end_erosion_indices)


def get_change_in_snow_depth(ds, start_events, end_events, index_of_event, x_sel, y_sel):
    '''
    Function to get snow height difference (absolute value) before and after an event
    Args:
        ds: dataset containing the snow-depth data ('snow_surface' variable)
        start_events: list of time indices at which the events of interest (accumulation or erosion) started
        end_events: list of time indices at which the events of interest (accumulation or erosion) ended
        index_of_event: index of event of interest in the lists of time indices
        x_sel: index of the x-coordinate of the point of interest
        y_sel: index of the y-coordinate of the point of interest
    Returns:
        difference in snow-depth between the start and end of the event
    '''
    
    start_date = start_events[index_of_event]
    end_date = end_events[index_of_event]
    
    start_snow_height = float(ds.snow_surface.isel(x=x_sel, y=y_sel, time=start_date))
    end_snow_height = float(ds.snow_surface.isel(x=x_sel, y=y_sel, time=end_date))
    difference = abs(end_snow_height - start_snow_height)
    
    return(difference)


# ======================================================================================================
# =========================== Simulating the snowpack evolution ========================================
# ======================================================================================================

def get_met_forcing(simulation_start_date='2021-12-06 00:00:00', file_start_date='2021-12-01 00:00:00', file_name='/home/mabonnet/Desktop/data/202202_finse_livox/met_obs.csv'):
    '''
    Function to get surface temperature forcing from meteorological files
    Args:
        simulation_start_date: date of beginning of simulation, format eg. '2021-04-28 16:30:00' for April 28th 2021, 4:30PM
        file_start_date: date from the meteorological file where we wish the forcing data to start (first timestamp in the forcing data), format eg. '2021-04-28 16:30:00' for April 28th 2021, 4:30PM
        file_name: path to the meteorological file
    Returns:
        met_temp_data: array containing surface temperatures for each meteorological timestamp (degrees Celcius)
        met_time_data: array containing the meteorological timestamps (since start of simulation, in s)
    '''

    time_series = []
    air_temp_series = []
    surf_temp_series = []

    with open(file_name) as met_data:
        reader = csv.reader(met_data)
        header = next(reader)
        start_date_in_sec = pd.to_datetime(simulation_start_date).timestamp()

        for row in reader:
            time_series.append(pd.to_datetime(row[0]).timestamp() - start_date_in_sec)    # TODO get index of rows automatically?

            if row[6]!='':
                air_temp_series.append(float(row[6]))
            else:
                air_temp_series.append(None)

            if row[8]!='':
                surf_temp_series.append(float(row[8]))
            else:
                surf_temp_series.append(None)

        index_start = time_series.index(pd.to_datetime(file_start_date).timestamp() - start_date_in_sec)

        time_series = time_series[index_start:]
        air_temp_series = air_temp_series[index_start:]
        surf_temp_series = surf_temp_series[index_start:]
    
    met_time_data = time_series  # times since start of simulation, in s
    met_temp_data = surf_temp_series

    return(met_time_data, met_temp_data)


def simulate_snowpack_evolution(ds, x_sel, y_sel, nb_iterations, end_accumulation_times, end_erosion_times,
                                start_accumulation, end_accumulation, start_erosion, end_erosion,
                                jj, dt, ro_layer, ro_water, ro_ice, t_old, tf, tsfc_default, dy_snow, age_layers, gamma, cp_snow, melt_flag, a1, a2,
                                met_temp_data=[None], met_time_data=[0],
                                new_snow_ro=150, fit_top_of_snowfall_to_curve=False):
    '''
    Function that simulates the evolution of the snowpack over a certain period of time
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
        
        jj: number of layers initially present                    # TODO add units and default values
        dt: timestep (s)
        ro_layer(1*max_nb_of_layers) array containing density value (kg per m**3) for each layer
        ro_water: density of water (kg per m**3)
        ro_ice: density of ice (kg per m**3)
        t_old: (1*max_nb_of_layers) array containing temperature value (degrees Celcius) for each layer
        tf: ice fusion temperature (degrees Celcius)
        tsfc_default: surface temperature (degrees Celcius)
        dy_snow: (1*max_nb_of_layers) array containing depth value (m) for each layer
        age_layers: (1*max_nb_of_layers) array containing age (s) of each layer
        gamma: (1*max_nb_of_layers) array containing zeros
        cp_snow: thermal capacity of snow
        melt_flag: (1*max_nb_of_layers) array containing melt value (1 or 0) for each layer
        a1, a2: exponential parameters, empirically calibrated

        met_temp_data: array containing surface temperatures for each meteorological timestamp (degrees Celcius), default [None], i.e. use tsfc=tsfc_default
        met_time_data: array containing the meteorological timestamps (since start of simulation, in s), default [0]
        
        new_snow_ro: density of newly fallen snow in kg per m**3, default 150

        fit_top_of_snowfall_to_curve: boolean, if True the height of snowfalls will be such that the snow depth is the same as the one measured by lidar, default False
    Returns:
        ro_layer_evolution: list of the states of layers' density through time, format [[ro_layer_1, ro_layer_2, ro_layer_3, ...]_time_1, [...]_time_2, [...], ...]
        depth_evolution: list of the states of layers' depth through time, format [[depth_layer_1, depth_layer_2, depth_layer_3, ...]_time_1, [...]_time_2, [...], ...]
        temperature_evolution: list of the states of layers' temperature through time, format [[temp_layer_1, temp_layer_2, temp_layer_3, ...]_time_1, [...]_time_2, [...], ...]
    '''
    # Initialize arrays to keep track of variables in time
    ro_layer_evolution = []
    depth_evolution = []
    temperature_evolution = []

    if met_temp_data[0] != None:
        tsfc = float(met_temp_data[0])
    else:
        tsfc = tsfc_default

    # Initialize indices of next accumulation/erosion events coming up (updated when their time is past)
    accumulation_index = 0
    erosion_index = 0
    temperature_index = 0
    
    for i in range(nb_iterations):

        if temperature_index<len(met_temp_data) and met_time_data[temperature_index]!=None and i*dt>=met_time_data[temperature_index]:
            if met_temp_data[temperature_index] != None:
                tsfc = float(met_temp_data[temperature_index])
            temperature_index += 1

        if accumulation_index<len(end_accumulation_times) and i*dt>=end_accumulation_times[accumulation_index]:
            if fit_top_of_snowfall_to_curve:
                snow_depth_total = sum(dy_snow[i] for i in range(0,jj))
                ddepth = float(ds.snow_surface.isel(x=x_sel, y=y_sel, time=end_accumulation[accumulation_index])) - snow_depth_total
            else:
                ddepth = get_change_in_snow_depth(ds, start_accumulation, end_accumulation, accumulation_index, x_sel, y_sel)
            ro_layer[jj] = new_snow_ro
            t_old[jj] = tsfc
            dy_snow[jj] = ddepth
            age_layers[jj] = (nb_iterations-i) * dt           # age in seconds at end of simulation
            if t_old[jj] <= 0:
                melt_flag[jj] = 0
            else:
                melt_flag[jj] = 1
            jj += 1
            accumulation_index += 1
    
        if erosion_index<len(end_erosion_times) and i*dt>=end_erosion_times[erosion_index]:
            ddepth = get_change_in_snow_depth(ds, start_erosion, end_erosion, erosion_index, x_sel, y_sel)
            erosion_index += 1
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
    
        # Update layers' parameters
        ro_layer, dy_snow = ddensity.ddensity_ml(ro_layer, tf, dt, ro_water, ro_ice, t_old, jj, dy_snow, a1, a2)
        t_old = snowtemp.snowtemp_ml(gamma, t_old, tsfc, jj, dt, ro_layer, cp_snow, tf, dy_snow, melt_flag)
    
        # Keep track of events
        ro_layer_evolution.append(ro_layer)
        depth_evolution.append(dy_snow)
        temperature_evolution.append(t_old)
        
    return(ro_layer_evolution, depth_evolution, temperature_evolution)


# ====================================================================================================
# =========================== Plotting the snowpack evolution ========================================
# ====================================================================================================

def plot_simul_and_signal(ds, x_sel, y_sel, depth_evolution, nb_layers_to_plot, data_start_date, dt, nb_iterations,
                          start_accumulation, end_accumulation, start_erosion, end_erosion, ice_layers_times_indices=None,
                          my_title='Comparison between lidar-measured and simulated snow depth', save_file=False, my_file_name='my_fig.png', my_figsize=(15, 7)):
    '''
    Function to plot the simulated snowpack and lidar signal on the same plot
    Args:
        ds: clean dataset
        x_sel: x-coordinate of the point of interest
        y_sel: y-coordinate of the point of interest
        depth_evolution: list of the states of layers' depth through time, format [[depth_layer_1, depth_layer_2, depth_layer_3, ...]_time_1, [...]_time_2, [...], ...]
        nb_layers_to_plot: number of layers to plot
        data_start_date: first date of the dataset, pandas datetime format
        dt: timestep used in snowpack simulation
        nb_iterations: number of iterations used in snowpack simulation
        
        start_accumulation: list of the indices of starting times of accumulations in ds
        end_accumulation: list of the indices of ending times of accumulations in ds
        start_erosion: list of the indices of starting times of erosions in ds
        end_erosion: list of the indices of ending times of erosions in ds

        ice_layers_times_indices: list of time indices where an ice ayer was detected, default None
        
        my_title: title of the figure, default 'Comparison between lidar-measured and simulated snow depth'
        save_file: boolean, default False
        my_file_name: name to be given to the saved file, default 'my_fig.png'
        my_figsize: figure size, default (15, 7)
    Returns:
    '''
    layers = np.zeros((nb_layers_to_plot, len(depth_evolution)))
    for i in range(len(depth_evolution)):
        layers[0][i] = depth_evolution[i][0]
    for layer_index in range(1, nb_layers_to_plot):
        for i in range(len(depth_evolution)):
            layers[layer_index][i] = depth_evolution[i][layer_index] + layers[layer_index-1][i]
    
    fig = plt.figure(figsize=my_figsize)
    times = pd.date_range(start=data_start_date,freq=str(dt)+'S',periods=nb_iterations)
    
    for layer_index in range(nb_layers_to_plot):
        plt.plot(times, layers[layer_index], label='layer '+str(layer_index+1))

    if ice_layers_times_indices != None:
        for time_index in ice_layers_times_indices:
            plt.plot(times[time_index], layers[-1][time_index]+0.001, c='y', marker='*', markersize=15, label='ice layer detected')
        
    ds.isel(x=x_sel, y=y_sel).snow_surface.plot(c='k', alpha=0.2)

    ds.isel(x=x_sel, y=y_sel, time=start_accumulation).snow_surface.plot(c='b', marker='^', markersize=6, linestyle='None', label='start accum.')
    ds.isel(x=x_sel, y=y_sel, time=end_accumulation).snow_surface.plot(c='g', marker='^', markersize=6, linestyle='None', label='end accum.')
    ds.isel(x=x_sel, y=y_sel, time=start_erosion).snow_surface.plot(c='m', marker='v', markersize=6, linestyle='None', label='start erosion')
    ds.isel(x=x_sel, y=y_sel, time=end_erosion).snow_surface.plot(c='r', marker='v', markersize=6, linestyle='None', label='end erosion')
    
    plt.legend()
    plt.title(my_title)
    
    if save_file:
        fig.savefig(my_file_name)
    
    plt.show(block=False)

    return()

