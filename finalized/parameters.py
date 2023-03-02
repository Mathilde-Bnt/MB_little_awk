import numpy as np
import xarray as xr
import little_awk_functions

# Coordinates of the point of interest
x_isel = 10            # index of x-coordinate of point of interest in dataset
y_isel = 10            # index of y-coordinate of point of interest in dataset

# Definition of the summer surface
start_summer_surface = '2022-10-07'
end_summer_surface =  '2022-10-15'

# Snow events detection parameters
time_window_std = 25             # index (number of timepoints that are taken into account when computing standard deviation)
std_threshold = 0.022            # m (mean height variation above which a point is considered to be part of a snow event)

# Initial state for compaction/temperature model, with 0 layers

# Adaptable parameters
tsfc = -5           # degrees Celcius (default surface temperature)
dt = 100            # s (time between each iteration of the simulation)
a1 = 0.0005         # m^-1.s^-1 (compaction coefficient)
a2 = 0.016          # m^3.kg^-1 (compaction coefficient)

max_nb_of_layers = 35    # (number of layers that can be detected by the model)

# Choose options for simulation
use_true_met = True                                # True if use the meteorological forcing from real data

# 'Fixed' parameters
tf = 0                           # degrees Celcius (fusion temperature of ice)
ro_water = 1000                  # kg.m^-3 (density of water)
ro_ice = 910                     # kg.m^-3 (density of ice)
cp_snow = 2106                   # J.kg^-1.K^-1 (thermal capacity of snow)
jj = 0                           # (initial number of layers)

# Meteorological forcing

if use_true_met:
    met_time, met_temp, met_wind = little_awk_functions.get_met_forcing(simulation_start_date='2022-10-07 00:00:00',
                        file_start_date='2022-10-07 00:00:00',
                        file_name='/home/mabonnet/Desktop/data/Data_2023/finse_meteo_obs.csv')
else:
    met_time, met_temp, met_wind = [0], [None], [None]

# Define structures to store snow parameters

ro_layer = np.zeros((max_nb_of_layers, 1))      # store the density of layers at a given timestamp, starting from the bottom
t_old = np.zeros((max_nb_of_layers, 1))         # store the temperature of layers at a given timestamp, starting from the bottom
dy_snow = np.zeros((max_nb_of_layers, 1))       # store the thickness of layers at a given timestamp, starting from the bottom
gamma = np.zeros((max_nb_of_layers, 1))         # store the thermal conductivity of layers at a given timestamp, starting from the bottom
melt_flag = np.zeros((max_nb_of_layers, 1))     # store the melt flag (boolean, True if there is melt in the layer) of layers at a given timestamp, starting from the bottom
age_layers = np.zeros((max_nb_of_layers))       # store the age (seconds until the end of the simulation) of layers at a given timestamp, starting from the bottom
