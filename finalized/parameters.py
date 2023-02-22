import numpy as np
import xarray as xr
import little_awk_functions

# Coordinates of the point of interest
x_sel = 10
y_sel = 10

# Definition of the summer surface
start_summer_surface = '2022-10-07'
end_summer_surface =  '2022-10-15'

# Snow events detection parameters
time_window_std = 25
std_threshold = 0.022

# Initial state for compaction/temperature model, with 0 layers

# Adaptable parameters
tsfc = -5           # degrees Celcius
dt = 100            # s
a1 = 0.0005         # m^-1.s^-1
a2 = 0.016          # m^3.kg^-1

max_nb_of_layers = 25

# Choose options for simulation
use_true_met = True
simul_fit_top_of_snowfall_to_curve = False
simul_erode_several_layers = False
simul_detect_ice = False

# If you wish to detect ice layers
simul_slope_threshold = None            # m.s^-1
simul_min_duration_in_s = None          # s

# 'Fixed' parameters
tf = 0                           # degrees Celcius
ro_water = 1000                  # kg.m^-3
ro_ice = 910                     # kg.m^-3
cp_snow = 2106                   # J.kg^-1.K^-1
jj = 0

# Meteorological forcing

if use_true_met:
    met_time, met_temp, met_wind = little_awk_functions.get_met_forcing(simulation_start_date='2022-10-07 00:00:00',
                        file_start_date='2022-10-07 00:00:00',
                        file_name='/home/mabonnet/Desktop/data/Data_2023/finse_meteo_obs.csv')
else:
    met_time, met_temp, met_wind = [0], [None], [None]

# Define structures to store snow parameters

ro_layer = np.zeros((max_nb_of_layers, 1))
t_old = np.zeros((max_nb_of_layers, 1))
dy_snow = np.zeros((max_nb_of_layers, 1))
gamma = np.zeros((max_nb_of_layers, 1))
melt_flag = np.zeros((max_nb_of_layers, 1))
age_layers = np.zeros((max_nb_of_layers, 1))

if simul_detect_ice:
    simul_a1_vector = np.zeros((max_nb_of_layers, 1))
else:
    simul_a1_vector = None
