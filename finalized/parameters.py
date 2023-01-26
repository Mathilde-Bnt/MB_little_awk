import numpy as np
import little_awk_functions

# Coordinates of the point of interest
x_sel=10
y_sel=10

# Snow events detection parameters
time_window_std = 25
std_threshold = 0.015

# Initial state for compaction/temperature model, with 0 layers

# Adaptable parameters
tsfc = -5
dt = 100
a1 = 0.0013
a2 = 0.021

max_nb_of_layers = 25

# Choose options for simulation
use_true_met = False   # set to True if want to use the correct temperature forcing

simul_fit_top_of_snowfall_to_curve = False

# 'Fixed' parameters
tf = 0
ro_water = 1000
ro_ice = 910
cp_snow = 2106
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