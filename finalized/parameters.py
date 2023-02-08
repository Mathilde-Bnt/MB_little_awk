import numpy as np
import xarray as xr
import little_awk_functions

# Coordinates of the point of interest
x_sel=10
y_sel=10

# Summer surface dataset

# If needed, create new dataset        # TODO make this check automatic to avoid issues if files exist already or not
# summer_data_raw = xr.open_mfdataset('/home/mabonnet/Desktop/data/2021_2022_livox_surfaces/202107*.nc')
# little_awk_functions.fill_in_missing_variables(summer_data_raw, 'surface')
# little_awk_functions.make_summer_netcdf(summer_data_raw, '/home/mabonnet/github/MB_little_awk/')

# Load dataset
summer_data = xr.open_dataset('/home/mabonnet/github/MB_little_awk/current_development/summer_surface.nc')

# Snow events detection parameters
time_window_std = 25
std_threshold = 0.015

# Initial state for compaction/temperature model, with 0 layers

# Adaptable parameters
tsfc = -5           # degrees Celcius
dt = 100            # s
a1 = 0.0013         # m^-1.s^-1
a2 = 0.021          # m^3.kg^-1

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
    met_time, met_temp, met_wind = little_awk_functions.get_met_forcing()
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
