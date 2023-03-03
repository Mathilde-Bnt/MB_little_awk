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
import math
from scipy.stats import sem


# ===========================================================================================
# =========================== Pre-processing dataset ========================================
# ===========================================================================================

def fill_in_missing_variables(ds, var_to_copy):
    '''
    Function to artificially create a 'mean' variable in a dataset
    Args:
        ds: clean dataset
        var_to_copy: string, name of variable to be copied into a 'mean' variable
    Returns:
    '''
    ds['mean'] = ds[var_to_copy]


def fill_in_missing_times(ds, day_date, output_dir_name, ds_corresponds_to_file):
    '''
    Function that reads a dayly netcdf dataset with variable 'mean', and creates and saves a new netcdf 
    data file with regular timestamps, containing the original data, and nan values where there were no data points originally
    Args:
        ds: day-long dataset from netcdf file, containing a 'mean' variable and the correct x and y coordinates
        day_date: string of the day's date, eg. April 30th 2021 would be '2021-04-30'
        output_dir_name: output directory where the new 'filled-in' netcdf file will be saved
        ds_corresponds_to_file: boolean, is True if ds is the dayly netcdf dataset corresponding to the day_date, is False if 
                the created netcdf should contain only nan
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

def get_snow_events(ds, x_isel, y_isel, time_window_std, std_threshold):
    '''
    Function that computes the dates of start and end times of snow events (snow accumulation and erosion)
    The events are defined as periods during which the snow depth rolling standard deviation is higher than 
    a given threshold
    We distinguish between accumulation and erosion of snow by looking at the snow depth before and after the event
    Args:
        ds: clean dataset with 'snow_depth' variable
        x_isel: x coordinate of the point of interest (index)
        y_isel: y coordinate of the point of interest (index)
        time_window_std: size of the rolling window to compute standard deviation
        std_threshold: standard deviation threshold above which the curve is considered to have strong variations > snow event
    Returns:
        start_accumulation_indices: list of indices (in ds) of the times corresponding to the start of an accumulation event
        start_erosion_indices: list of indices (in ds) of the times corresponding to the end of an accumulation event
        end_accumulation_indices: list of indices (in ds) of the times corresponding to the start of an erosion event
        end_erosion_indices: list of indices (in ds) of the times corresponding to the end of an erosion event
    '''
    # Compute standard deviation values around each point
    stdev = ds.isel(x=x_isel, y=y_isel).snow_depth.rolling(time=time_window_std, center=True).std(dim='time').values
    
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
    
        start_snow_height = float(ds.snow_depth.isel(x=x_isel, y=y_isel, time=start_date))
        end_snow_height = float(ds.snow_depth.isel(x=x_isel, y=y_isel, time=end_date))
    
        if start_snow_height < end_snow_height:
            # Accumulation if the snow height rises during the event
            start_accumulation_indices.append(start_time_indices[index])
            end_accumulation_indices.append(end_time_indices[index])
        else:
            # Erosion if the snow height falls during the event
            start_erosion_indices.append(start_time_indices[index])
            end_erosion_indices.append(end_time_indices[index])
        
    return(start_accumulation_indices, start_erosion_indices, end_accumulation_indices, end_erosion_indices)


def get_change_in_snow_depth(ds, start_events, end_events, index_of_event, x_isel, y_isel):
    '''
    Function to compute snow height difference (absolute value) before and after an event
    Args:
        ds: dataset containing the snow depth data ('snow_depth' variable)
        start_events: list of time indices at which the events of interest (accumulation or erosion) started
        end_events: list of time indices at which the events of interest (accumulation or erosion) ended
        index_of_event: index of the event of interest in the lists of time indices
        x_isel: index of the x-coordinate of the point of interest
        y_isel: index of the y-coordinate of the point of interest
    Returns:
        difference in snow depth between the start and end of the event
    '''
    # Get timing of event
    start_date = start_events[index_of_event]
    end_date = end_events[index_of_event]
    
    # Get snow heights at start and end of event
    start_snow_height = ds.snow_depth.isel(x=x_isel, y=y_isel, time=start_date)
    end_snow_height = ds.snow_depth.isel(x=x_isel, y=y_isel, time=end_date)
    difference = np.abs(end_snow_height - start_snow_height)
    
    return(difference)


# ======================================================================================================
# =========================== Simulating the snowpack evolution ========================================
# ======================================================================================================

def get_met_forcing_df(met_file_name, dt, sim_start_date, sim_end_date, format_time='%d.%m.%Y %H:%M'):
    '''
    Function to get surface temperature, wind speed forcing and new snow density from meteorological files in a dataframe
    Args:
        file_name: string, path to the meteorological csv file, expected to be extracted from seklima.no website, and containing
                Time, Air_temp and Wind_speed data
        dt: time in seconds between two simulation iterations
        simulation_start_date: date of start of simulation in pd.datetime format
    Returns:
        met_df: dataframe containing the interpolated meteorological data to match simuation timestamps (time in pd.datetime format, air temperature(degrees Celcius),
                wind speed (m.s^-1) and new snow density (kg.m^⁻3)), indices are the time in seconds since the start of the simulation
    '''
    print('WARNING: check your meteorological file is in the right format (cf docstring of get_met_forcing())')
    
    met_df = pd.read_csv(met_file_name)
    met_df['Time'] = pd.to_datetime(met_df['Time'], format=format_time)
    met_df.Wind_speed.loc[met_df.Wind_speed == '-'] = 0
    met_df.Wind_speed = met_df.Wind_speed.astype(float)
    met_df.Air_temp = met_df.Air_temp.astype(float)
    
    met_df = met_df.set_index('Time')
    start = pd.to_datetime(sim_start_date).strftime('%Y-%m-%d %H:%M')
    end = pd.to_datetime(sim_end_date).strftime('%Y-%m-%d %H:%M')
    met_df = met_df.resample(str(dt)+'S', origin=start, closed='right').interpolate()[start:]
    
    met_df['Ro_new_snow'] = 150
    met_df.Ro_new_snow.loc[met_df.Wind_speed >= 6] = 250
    
    return(met_df[:end])



def simulate_snowpack_evolution_df(ds, met_df, x_isel, y_isel, nb_iterations, max_nb_of_layers, end_accumulation_times, end_erosion_times,
                                start_accumulation, end_accumulation, start_erosion, end_erosion,
                                dt, ro_water, ro_ice, tf, cp_snow,
                                a1, a2):
    '''
    Function that simulates the evolution of the snowpack over a certain period of time: at each timestamp, the new density, height and
    temperature of each layer is computed, according to SnowModel (2006)
    Options: use correct (rather than constant) meteorological forcing in dataframe form
    Args:
        ds: clean dataset (x, y, time)
        x_isel: x-coordinate of the point of interest (index)
        y_isel: y-coordinate of the point of interest (index)
        nb_iterations: number of iterations
        max_nb_of_layers: int, maximum number of layers that can be detected at the point of interest
        
        end_accumulation_times: list of ending times of accumulations in seconds since data starting date
        end_erosion_times: list of ending times of erosions in seconds since data starting date
        start_accumulation: list of the indices of starting times of accumulations in ds
        end_accumulation: list of the indices of ending times of accumulations in ds
        start_erosion: list of the indices of starting times of erosions in ds
        end_erosion: list of the indices of ending times of erosions in ds
        
        jj: number of layers initially present
        dt: timestep (s)
        
        ro_water: density of water (kg.m^-3)
        ro_ice: density of ice (kg.m^-3)
        
        tf: ice fusion temperature (degrees Celcius)        
        cp_snow: thermal capacity of snow (J.kg^-1.K^-1)
        a1, a2: exponential parameters, empirically calibrated, a1 in m^-1.s^-1, a2 in m^3.kg^-1

        met_data: dataframe as constructed in the get_met_forcing_df() function, default None
        use_true_met: boolean, if True then met_data is assumed not to be None, and is used 
                as forcing in the simulation, default False

    Returns:
        ro_layer_evolution: list of the states of layers' density (kg.m^-3) at each timestamp,
                format [[ro_layer_1, ro_layer_2, ro_layer_3, ...], [...], ...]
        thickness_evolution: list of the states of layers' thickness (m) at each timestamp,
                format [[thickness_layer_1, thickness_layer_2, thickness_layer_3, ...], [...], ...]
        temperature_evolution: list of the states of layers' temperature (degrees Celcius) at each timestamp,
                format [[temp_layer_1, temp_layer_2, temp_layer_3, ...], [...], ...]
        age_layers_evolution: list of the age of each layer (seconds until the end of the simulation) at each timestamp,
                format [[age_layer_1, age_layer_2, age_layer_3, ...], [...], ...]
    '''

    # INITIALIZATION

    # Initialize arrays to keep track of variables in time
    ro_layer_evolution = np.empty((nb_iterations, max_nb_of_layers))
    thickness_evolution = np.empty((nb_iterations, max_nb_of_layers))
    temperature_evolution = np.empty((nb_iterations, max_nb_of_layers))
    age_layers_evolution = np.empty((nb_iterations, max_nb_of_layers))
    melt_flag_evolution = age_layers_evolution.copy()
    
    # ro_layer(1*max_nb_of_layers) array containing density value (kg.m^-3) for each layer
    # t_old: (1*max_nb_of_layers) array containing temperature value (degrees Celcius) for each layer
    # dy_snow: (1*max_nb_of_layers) array containing thickness value (m) for each layer
    # age_layers: (1*max_nb_of_layers) array containing age (s) of each layer
    # gamma: (1*max_nb_of_layers) array containing zeros
    # melt_flag: (1*max_nb_of_layers) array containing melt value (1 or 0) for each layer
    
    # Define structures to store snow parameters

    ro_layer = np.empty((max_nb_of_layers))      # store the density of layers at a given timestamp, starting from the bottom
    t_old = np.empty((max_nb_of_layers))         # store the temperature of layers at a given timestamp, starting from the bottom
    dy_snow = np.empty((max_nb_of_layers))       # store the thickness of layers at a given timestamp, starting from the bottom
    gamma = np.empty((max_nb_of_layers))         # store the thermal conductivity of layers at a given timestamp, starting from the bottom
    melt_flag = np.empty((max_nb_of_layers))     # store the melt flag (boolean, True if there is melt in the layer) of layers at a given timestamp, starting from the bottom
    age_layers = np.empty((max_nb_of_layers))       # store the age (seconds until the end of the simulation) of layers at a given timestamp, starting from the bottom

    # Initialize indices of next accumulation/erosion events coming up (updated when their time is past)
    accumulation_index, erosion_index = 0, 0

    # RUN MODEL
    i=0
    jj = 0                           # (initial number of layers)
    for tstep, row in met_df.iterrows(): 
        # Update surface temperature and snow density if needed
        tsfc = row.Air_temp
        
        # Detection of accumulations
        if accumulation_index<len(end_accumulation_times) and i*dt>=end_accumulation_times[accumulation_index]:     # if an accumulation was passed
            
            ddepth = get_change_in_snow_depth(ds, start_accumulation, end_accumulation, accumulation_index, x_isel, y_isel)
            ro_layer[jj] = row.Ro_new_snow
            t_old[jj] = tsfc
            dy_snow[jj] = ddepth
            age_layers[jj] = (nb_iterations-i) * dt     # age in seconds at end of simulation

            jj += 1
            accumulation_index += 1     # next accumulation that will come up
    
        # Detection of erosions
        if erosion_index<len(end_erosion_times) and i*dt>=end_erosion_times[erosion_index]:     # if an erosion was passed

            ddepth = get_change_in_snow_depth(ds, start_erosion, end_erosion, erosion_index, x_isel, y_isel)
            erosion_index += 1     # next erosion that will come up
            if jj>0:
                if dy_snow[jj-1] > ddepth:
                    dy_snow[jj-1] = dy_snow[jj-1] - ddepth
                else:
                    jj -= 1
                    dy_snow[jj] = np.nan
                    ro_layer[jj] = np.nan
                    t_old[jj] = np.nan
                    age_layers[jj] = np.nan

        melt_flag[t_old>0] = 1
        melt_flag[t_old<=0] = 0          
    
        # Update layers' parameters (fortran code model)
        ro_layer, dy_snow = ddensity.ddensity_ml(ro_layer, tf, dt, ro_water, ro_ice, t_old, jj, dy_snow, a1, a2)
        t_old = snowtemp.snowtemp_ml(gamma, t_old, tsfc, jj, dt, ro_layer, cp_snow, tf, dy_snow, melt_flag)

        # Keep track of events
        ro_layer_evolution[i] = ro_layer
        thickness_evolution[i] = dy_snow
        temperature_evolution[i] = t_old
        age_layers_evolution[i] = age_layers
        melt_flag_evolution[i] = melt_flag
        i += 1
        
    return(ro_layer_evolution, thickness_evolution, temperature_evolution, age_layers_evolution, melt_flag_evolution)


# ====================================================================================================
# =========================== Plotting the snowpack evolution ========================================
# ====================================================================================================

def plot_simul_and_signal(ds, x_isel, y_isel, thickness_evolution, nb_layers_to_plot, data_start_date, dt, nb_iterations,
                          start_accumulation, end_accumulation, start_erosion, end_erosion, ice_layers_times=None,
                          fig_title='Comparison between lidar-measured and simulated snow depth', save_file=False,
                          fig_file_name='my_fig.png', fig_figsize=(16, 7), show_layers_in_legend=False):
    '''
    Function to plot the simulated snowpack layers and lidar signal on the same plot
    Args:
        ds: clean dataset
        x_isel: x-coordinate of the point of interest
        y_isel: y-coordinate of the point of interest
        thickness_evolution: list of the states of layers' thickness (m) at each timestamp, 
                format [[thickness_layer_1, thickness_layer_2, thickness_layer_3, ...], [...], ...]
        nb_layers_to_plot: number of layers to plot, starting from the bottom
        data_start_date: first date of the dataset, pandas datetime format
        dt: timestep used in snowpack simulation (s)
        nb_iterations: number of iterations used in snowpack simulation
        
        start_accumulation: list of the indices (in ds) of starting times of accumulations
        end_accumulation: list of the indices (in ds) of ending times of accumulations
        start_erosion: list of the indices (in ds) of starting times of erosions
        end_erosion: list of the indices (in ds) of ending times of erosions

        ice_layers_times: list of time indices where an ice layer was detected, default None
        
        fig_title: title of the figure, default 'Comparison between lidar-measured and simulated snow depth'
        save_file: boolean, default False
        fig_file_name: name to be given to the saved file, default 'my_fig.png'
        fig_figsize: figure size, default (16, 7)
        show_layers_in_legend: boolean, is True if all layer numbers should appear in the legend (more text)
    Returns:
    '''
    # Construct the data to plot: height of each layer as a function of time
    layers = np.zeros((nb_layers_to_plot, len(thickness_evolution)))
    
    # First layer
    layers[0] = np.array(thickness_evolution)[:,0]
    # Next layers
    for layer_index in range(1, nb_layers_to_plot):
        layers[layer_index] = np.array(thickness_evolution)[:,layer_index] + layers[layer_index-1]

    # Define figure and timestamps
    fig = plt.figure(figsize=fig_figsize)
    times = pd.date_range(start=data_start_date,freq=str(dt)+'S',periods=nb_iterations)
    
    # Plot each layer
    for layer_index in range(nb_layers_to_plot):
        if show_layers_in_legend:
            plt.plot(times, layers[layer_index], label='layer '+str(layer_index+1))
        else:
            plt.plot(times, layers[layer_index])

    # Plot a star at each ice layer detection
    if ice_layers_times != None:
        for time_index in ice_layers_times:
            plt.plot(times[time_index], layers[-1][time_index]+0.001, c='y', marker='*', markersize=15, label='ice layer detected')
    
    # Plot the lidar signal
    ds.isel(x=x_isel, y=y_isel).snow_depth.plot(c='k', alpha=0.2, label='lidar signal')

    # Plot the start and end of detected snow events on the lidar curve
    ds.isel(x=x_isel, y=y_isel, time=start_accumulation).snow_depth.plot(c='b', marker='^', markersize=6, linestyle='None', label='start accum.')
    ds.isel(x=x_isel, y=y_isel, time=end_accumulation).snow_depth.plot(c='g', marker='^', markersize=6, linestyle='None', label='end accum.')
    ds.isel(x=x_isel, y=y_isel, time=start_erosion).snow_depth.plot(c='m', marker='v', markersize=6, linestyle='None', label='start erosion')
    ds.isel(x=x_isel, y=y_isel, time=end_erosion).snow_depth.plot(c='r', marker='v', markersize=6, linestyle='None', label='end erosion')
    
    plt.ylabel('snow depth (m)')
    plt.xlabel('time (date)')

    plt.legend()
    plt.title(fig_title)
    
    if save_file:
        fig.savefig(fig_file_name)
    
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

def get_data_from_ref(index_of_ref_layer, ro_layer_array, thickness_array):
    '''
    Function that computes some characteristic values (SWE, height, average density) above a "ref" interface
    Args:
        index_of_ref_layer: index of the layer that is directly underneath the ref interface
        ro_layer_array: array of the densities of each layer at the end of the simulation
        thickness_array: array of the thickness of each layer at the end of the simulation
    Returns:
        swe_from_ref: value of the SWE above the given interface, in m
        height_from_ref: height of the snowpack above the given interface, in m
        ave_density_from_ref: averaged density of the snow above the given interface, in kg.m^-3
    '''
    swe_from_ref = np.dot(np.array(ro_layer_array[index_of_ref_layer+1:]), np.array(thickness_array[index_of_ref_layer+1:])) / 1000
    height_from_ref = sum(thickness_array[i] for i in range(index_of_ref_layer+1, len(thickness_array)))
    ave_density_from_ref = np.dot(np.array(ro_layer_array[index_of_ref_layer+1:]), np.array(thickness_array[index_of_ref_layer+1:])) / height_from_ref
    
    return(swe_from_ref, height_from_ref, ave_density_from_ref)


def get_depth_layers_indices(bottom_depth, sampling_length, thickness_array, start_from_top=False):
    '''
    Function that creates an array of the layer index at each regularly sampled depth (useful to compare with regularly sampled observation profiles),
    starting from the top or bottom of the snowpack
    Args:
        bottom_depth: depth of the lowest sample, in meters, counted positively from the surface of the snowpack downwards (included) if start_from_top=True,
                      counted positively from the ground upwards if start_from_top=False
        sampling_length: height between two consecutive sampled depths
        thickness_array: array of the thickness of each layer, in meters, at the end of the simulation
        start_from_top: boolean, is True if the sampling starts at the top of the snowpack (included) downwards, and False if it starts
                        at bottom_depth (included) upwards. Default False
    Returns:
        layers_array: array containing the indices of the layers corresponding to each sampled depth (which layer was sampled)
    '''
    if start_from_top:
        current_layer_index = next(len(thickness_array) - i for i, j in enumerate(reversed(thickness_array), 1) if j != 0)  # index of top layer
        remaining_depth_in_current_layer = thickness_array[current_layer_index]
        current_depth = 0      # from the top of the snowpack
        layers_array = [current_layer_index]    # top of the snowpack
        leftover_depth = sampling_length    # depth to go to the next sample
    
        while current_depth < bottom_depth and current_layer_index > -1:
            while remaining_depth_in_current_layer >= leftover_depth and current_depth < bottom_depth and current_layer_index > -1:   # next sample is in current layer
                layers_array.append(current_layer_index)
                remaining_depth_in_current_layer -= leftover_depth
                current_depth += sampling_length
                leftover_depth = sampling_length
            while remaining_depth_in_current_layer < leftover_depth and current_depth < bottom_depth and current_layer_index > -1:   # next sample is not in current layer
                current_layer_index -= 1   # go to next layer
                leftover_depth -= remaining_depth_in_current_layer   # what is left to "subtract" to this new layer
                remaining_depth_in_current_layer = thickness_array[current_layer_index]
                
        layers_array = list(reversed(layers_array))
        
    else: # sampling starts from bottom_depth included
        current_layer_index = 0    # index of bottom layer
        remaining_depth_in_current_layer = thickness_array[current_layer_index]
        layers_array = []
        leftover_depth = bottom_depth    # depth to go to the next sample
        
        # Get all other indices
        while current_layer_index < len(thickness_array):
            while remaining_depth_in_current_layer >= leftover_depth and current_layer_index < len(thickness_array):   # next sample is in current layer
                layers_array.append(current_layer_index)
                remaining_depth_in_current_layer -= leftover_depth
                leftover_depth = sampling_length
            while remaining_depth_in_current_layer < leftover_depth and current_layer_index < len(thickness_array):   # next sample is not in current layer
                current_layer_index += 1   # go to next layer
                leftover_depth -= remaining_depth_in_current_layer   # what is left to "subtract" to this new layer
                if current_layer_index < len(thickness_array):
                    remaining_depth_in_current_layer = thickness_array[current_layer_index]
    
    return(layers_array)


def density_profile(bottom_depth, sampling_length, thickness_array, ro_array, start_from_top=False):
    '''
    Function that creates an array of the layer density at each (regularly) sampled depth
    Args:
        bottom_depth: depth of the lowest sample, in meters, counted positively from the surface of the snowpack downwards (included)
        sampling_length: height between two consecutive sampled depths
        thickness_array: array of the thickness of each layer, in meters, at the end of the simulation
        ro_array: array of the densities of each layer, in kg.m^-3, at the end of the simulation
        start_from_top: boolean, is True if the sampling starts at the top of the snowpack (included) downwards, and False if it starts
                        at bottom_depth (included) upwards. Default False
    Returns:
        ro_profile: array containing the simulated densities at each sampled depth
    '''
    indices_array = get_depth_layers_indices(bottom_depth, sampling_length, thickness_array, start_from_top)
    ro_profile = [ro_array[i] for i in indices_array]
    
    return(ro_profile)


def temp_profile(bottom_depth, sampling_length, thickness_array, temp_array, start_from_top=False):
    '''
    Function that creates an array of the layer temperature at each (regularly) sampled depth
    Args:
        bottom_depth: depth of the lowest sample, in meters, counted positively from the surface of the snowpack downwards (included)
        sampling_length: height between two consecutive sampled depths
        thickness_array: array of the thickness of each layer, in meters, at the end of the simulation
        temp_array: array of the temperatures of each layer at the end of the simulation
        start_from_top: boolean, is True if the sampling starts at the top of the snowpack (included) downwards, and False if it starts
                        at bottom_depth (included) upwards. Default False
    Returns:
        temp_profile: array containing the simulated temperatures at each sampled depth
    '''
    indices_array = get_depth_layers_indices(bottom_depth, sampling_length, thickness_array, start_from_top)
    temp_profile = [temp_array[i] for i in indices_array]
    
    return(temp_profile)
