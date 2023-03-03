import numpy as np
import xarray as xr
import little_awk_functions


data_dir = '/home/simonfi/github/MB_little_awk/data'

simulation_start = '2022-11-01 00:00:00'
simulation_end = '2023-02-06 12:00:00'

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
dt = 100            # s (time between each iteration of the simulation)
a1 = 0.0005         # m^-1.s^-1 (compaction coefficient)
a2 = 0.016          # m^3.kg^-1 (compaction coefficient)

max_nb_of_layers = 35    # (number of layers that can be detected by the model)

# 'Fixed' parameters
tf = 0                           # degrees Celcius (fusion temperature of ice)
ro_water = 1000                  # kg.m^-3 (density of water)
ro_ice = 910                     # kg.m^-3 (density of ice)
cp_snow = 2106                   # J.kg^-1.K^-1 (thermal capacity of snow)
jj = 0                           # (initial number of layers)




