import numpy as np
import little_awk_functions


x_sel=10
y_sel=10

time_window_std = 25
std_threshold = 0.015

# Initial state for compaction/temperature model, with 0 layers

# Adaptable parameters
tsfc = -5
dt = 100
a1 = 0.0013
a2 = 0.021

max_nb_of_layers = 25

use_true_temp = False   # set to True if want to use the correct temperature forcing

simul_new_snow_ro = 180
simul_fit_top_of_snowfall_to_curve = False

# 'Fixed' parameters
tf = 0
ro_water = 1000
ro_ice = 910
cp_snow = 2106
jj = 0

# Meteorological forcing

if use_true_temp:
    met_time, met_temp = little_awk_functions.get_met_forcing()
else:
    met_time, met_temp = [0], [tsfc]

# Define structures to store snow parameters

ro_layer = np.zeros((max_nb_of_layers, 1))
t_old = np.zeros((max_nb_of_layers, 1))
dy_snow = np.zeros((max_nb_of_layers, 1))
gamma = np.zeros((max_nb_of_layers, 1))
melt_flag = np.zeros((max_nb_of_layers, 1))
age_layers = np.zeros((max_nb_of_layers, 1))