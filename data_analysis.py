# -*- coding: utf-8 -*-

"""
Created on Fri Feb 25 10:09:39 2022

@author: Lucas & Gaspard
"""

#%% Cell 0 ##

## module to install
## =================

# pip install xarray
# pip install netCDF4

#%% Cell 1 ##

## Importing modules
## =================

import matplotlib
import matplotlib.colors as clrs
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.cm as cm
import xarray as xr
import random as rd
import pandas as pd
import numpy as np
import imageio
import pickle
import copy
import os

from datetime import datetime
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import find_peaks


# Command to be executed in the console before starting to display the curves in a pop-up window
# %matplotlib qt


#%% Cell 2 ##

## Parameters
## ==========

# Indices of the point studied on the footprint
x_index = 106 # value of x as an index
y_index = 48  # value of y as an index

# Time interval over which the average is calculated to define the altitude of the point studied in summer
dt = 100
#start
#end

# parameters for smoothing
time_parameter    = 11 # the time median is performed on time_parameter points
space_parameter_x = 11 # the spatial median is performed on parametre_space_x points
space_parameter_y = 11 # the spatial median is performed on parametre_space_y points

# parameters for peak / fall detection
valeur_width_pics = 0.1           # parameters of the find_peaks method
value_prominence_falls = 0.01 #0.01
value_width_falls = 0.1

#%% Cell 3 ##

## Loading the useful dataset
# ===========================

# Loading
path = r'C:\Users\lucas\OneDrive\Bureau\ECL\UE\PAr\netcdf\petit\2021*.nc' #\2022*.nc allows to analyse only data from January, which makes calculations easier

ds_initial = xr.open_mfdataset(path) # dataset that will not be modified
ds = ds_initial                      # dataset that we will modify

# Reduce the size of the dataset
ds = ds.drop('max')
ds = ds.drop('min')
ds = ds.drop('stdev')
ds = ds.drop('count')
ds = ds.drop('idw')


#%% Cell 4 ##

## Location on the footprint
# ==========================

def n_d(axe,n):  # number to distance
    """ Convert a rank as an index into a distance (m) along x and y """
    if axe == 'x':
        d=-19.15+n* 0.1
    if axe == 'y':
        d=4.55-n*0.1
    return d

def d_n(axe,d):  # distance to number
    """ Convert a distance (m) along x and y into a rank as an index """
    if axe == 'x':
       n=(d+19.15)/0.1
    if axe == 'y':
       n=(4.55-d)/0.1
    return round(n)

'''
# If you want to give values in metres and not in indices

x_distance = -8.5 # value in metres
y_distance = -0.3     # value in metres

x_distance_in_ds = ds['mean'].sel(x=x_distance,y=y_distance,method="nearest").coords['x'].to_numpy()
y_distance_in_ds = ds['mean'].sel(x=x_distance,y=y_distance,method="nearest").coords['y'] .to_numpy()

x_index = d_n('x',x_distance_in_ds)
y_index = d_n('y',y_distance_in_ds)
'''


#%% Cell 5 ##

## Suppression of outliers
# ========================

# Definition of the summer altitude to define too small outliers
def average_beginning_curve(dt,ds,x_index,y_index) : # dt = time over which we will take an average
    '''the altitude at a point (x,y) is defined as the average height over a period dt during the summer'''
    list_of_beginning_values = ds['mean'].isel(x=x_index,y=y_index).to_numpy()[:dt]
    newlist = [x for x in list_of_beginning_values if np.isnan(x) == False]
    average = np.mean(newlist)
    return(average)

def average_bis(ds,start,end) :
    dd = ds['mean'].sel(time = slice(start,end)).mean(dim = 'time')
    return(dd.compute())

# ne plus passer en numpy
# faire une time range 
    
    
# Removal of outliers
def remove_outliers(ds,start,end) :
    '''we start by removing the extreme and therefore aberrant values from the dataset'''

    # Definition of the minimum snow altitude not to be exceeded
    down_safety_interval  = 0.1  # delete the values that are 10 cm below the average calculated on the first values
    average = average_bis(ds,start,end)
    low_limit =  average - down_safety_interval

    # Definition de l'altitude maximale de neige supposee n'etre pas depassee
    up_safety_interval = 4  # en m
    limit_high = average + up_safety_interval

    # Suppression des valeurs aberrantes
    ds_clean = ds.where(ds['mean'] < limit_high)
    ds_clean = ds.where(ds['mean'] > low_limit)

    return(ds_clean)


#%% Cell 6 ##

## Lissage de la courbe
# =====================

# Filtre median temporel
def median_time_filter(ds, time_parameter) :
    '''on applique une mediane de fenêtre glissante sur un nombre de points time_parameter'''
    ds_median_time_filter = ds.rolling(time=time_parameter, center=True).median(dim = 'time')
    return(ds_median_time_filter)

# Filtre median en espace
def median_space_filter(ds , space_parameter_x , space_parameter_y) :
    '''on applique une mediane de fenêtre glissante sur un nombre de points
    - space_parameter_x selon x
    - space_parameter_y selon y '''
    ds_median_space_filter = ds.rolling({'x':space_parameter_x,'y':space_parameter_y}, center=True).median(dim = ['x','y'])
    return(ds_median_space_filter)

# a rectifier

ds_clean = remove_outliers(ds,start,end)
ds_median_time_filter  = median_time_filter (ds_clean, time_parameter)
ds_median_space_filter = median_space_filter(ds_median_time_filter, space_parameter_x , space_parameter_y)


#%% Cell 7 ##

## Characteristic points recovery
# ===============================

sd = ds_median_space_filter['mean'].isel(x=x_index,y=y_index)
first_derivative_np = np.gradient(sd)
second_derivative_np = np.gradient(first_derivative_np)


# Peak recovery with second derivative
# ====================================

ds_median_space_filter = ds_median_space_filter.assign(second_derivative = second_derivative_np)
ds_median_space_filter = ds_median_space_filter.assign(opposite_second_derivative = - second_derivative_np)  # the derivative is inverted to obtain the minimums

# Recovering peak indices in dataset
peaks_pics, properties = find_peaks(ds_median_space_filter['opposite_second_derivative'], prominence=value_prominence_falls, width=value_width_falls)

# Converting indices into dates
list_time_pics = sd.isel(time=peaks_pics).time

'''
list_dates_axis_abscisses = sd.coords['time'].values
list_time_pics = []
for i in peaks_pics:
    list_time_pics.append(list_dates_axis_abscisses[i])
'''

# Recovering the snow height at the peaks
'''
order_pics = sd.isel(time=peaks_pics)['mean']
order_pics_second_derivative = sd.isel(time=peaks_pics)['mean']
'''
order_pics = []
order_pics_second_derivative = []
for i in peaks_pics:
    order_pics.append(sd.isel(time=i).to_numpy())
    order_pics_second_derivative.append(ds_median_space_filter['second_derivative'][i].to_numpy())


# Recovery of the beginning of the snowfall with second drift
# ===========================================================

ds_median_space_filter = ds_median_space_filter.assign(second_derivative = second_derivative_np) 

# Recovering peak indices in dataset
peaks_falls, properties = find_peaks(ds_median_space_filter['second_derivative'], prominence=value_prominence_falls, width=value_width_falls)

# Converting indices into dates
list_falls_times = []
for i in peaks_falls:
    list_falls_times.append(list_dates_axis_abscisses[i])

# Recovery of snow height at the start of the fall

order_falls = []
order_falls_second_derivative = []
for i in peaks_falls:
    order_falls.append(sd.isel(time=i).to_numpy())
    order_falls_second_derivative.append(ds_median_space_filter['second_derivative'][i].to_numpy())



#%% Cell 8 ## dunes filtering using first derivative

## point qui ne vont pas derivee premiere
# =======================================

valeur_width_pics_first_derivative = 0.5           # parameters of the find_peaks method
value_prominence_falls_first_derivative = 0.022
value_width_falls_first_derivative = 0.5

sd = ds_median_space_filter['mean'].isel(x=x_index,y=y_index)


# recupération de spoints problématiques avec la derivee premiere
# ====================================

ds_median_space_filter = ds_median_space_filter.assign(opposite_first_derivative = - first_derivative_np)  # the derivative is inverted to obtain the minimums
ds_median_space_filter = ds_median_space_filter.assign(first_derivative = first_derivative_np) 

# Recovering peak indices in dataset
peaks_pics_first_derivative, properties_first_derivative = find_peaks(ds_median_space_filter['opposite_first_derivative'], prominence=value_prominence_falls_first_derivative, width=valeur_width_pics_first_derivative)



# Recovery of snow height
order_first = []
for i in peaks_pics_first_derivative:
    order_first.append(ds_median_space_filter['first_derivative'][i].to_numpy())

# suppression des points en trop sur la courbe

peaks_falls = list(peaks_falls)
for i in range(len(peaks_pics_first_derivative)) :
    j=0
    while j < len(peaks_falls) :
        if peaks_falls[j]-1 == peaks_pics_first_derivative[i] :
            peaks_falls.pop(j)
            order_falls.pop(j)
            order_falls_second_derivative.pop(j)
        j+=1


#%% Cell 9 ##

## Removal of excess points to have a fall / peak alternation
# ===========================================================

# lists are created by indicating whether the point is a fall or a peak

list_tuple_falls = []
for fall in peaks_falls :
    list_tuple_falls.append((fall,'fall'))

liste_tuple_pics = []
for pic in peaks_pics :
    liste_tuple_pics.append((pic,'pic'))

liste_tuple_points_caracteristics = list_tuple_falls + liste_tuple_pics

liste_tuple_points_caracteristics.sort() # sorting the list with all characteristic points

# The first characteristic point is forced to be a snowfall

if liste_tuple_points_caracteristics[0][1] == 'pic' :
    list_tuple_falls.append((0,'fall'))

# Removal of excess points by case disjunction

def alternation_characteristic_points(liste_tuple_points_caracteristics) :
    index = 0

    if len(liste_tuple_points_caracteristics)>=3:
        while index <= (len(liste_tuple_points_caracteristics) - 3) :
            if liste_tuple_points_caracteristics[index][1] == 'pic' :
                if liste_tuple_points_caracteristics[index+1][1] == 'pic' :
                    if liste_tuple_points_caracteristics[index+2][1] == 'pic' : # PPP
                        liste_tuple_points_caracteristics.pop(index)
                        liste_tuple_points_caracteristics.pop(index+1)
                        index += 0
                    else : # PPF
                        liste_tuple_points_caracteristics.pop(index+1)
                        index += 1
                else :
                    if liste_tuple_points_caracteristics[index+2][1] == 'pic' : # PFP
                        index += 2
                    else : # PFF
                        liste_tuple_points_caracteristics.pop(index+1)
                        index += 1
            else :
                if liste_tuple_points_caracteristics[index+1][1] == 'pic' :
                    if liste_tuple_points_caracteristics[index+2][1] == 'pic' : # FPP
                        liste_tuple_points_caracteristics.pop(index+1)
                        index += 1
                    else : # FPF
                        index += 2
                else :
                    if liste_tuple_points_caracteristics[index+2][1] == 'pic' : # FFP
                        liste_tuple_points_caracteristics.pop(index+1)
                        index += 1
                    else : # FFF
                        liste_tuple_points_caracteristics.pop(index+1)
                        liste_tuple_points_caracteristics.pop(index+1)
                        index += 0

    elif len(liste_tuple_points_caracteristics) == 2 :
        if liste_tuple_points_caracteristics[1][1] == 'fall' : # If the first 2 are falls
            liste_tuple_points_caracteristics.pop(1)

    return(liste_tuple_points_caracteristics)

liste_tuple_points_caracteristics = alternation_characteristic_points(liste_tuple_points_caracteristics)


# Recovering the snow height at the peaks
order_pics_with_alternation = []
order_pics_second_derivative_with_alternation = []
peaks_pics_ecreme = [liste_tuple_points_caracteristics[i][0] for i in range(len(liste_tuple_points_caracteristics)) if liste_tuple_points_caracteristics[i][1] == 'pic']
for i in peaks_pics_ecreme:
    order_pics_with_alternation.append(sd.isel(time=i).to_numpy())
    order_pics_second_derivative_with_alternation.append(ds_median_space_filter['second_derivative'][i].to_numpy())


# Recovery of snow height at the beginning of the fall
order_falls_with_alternation = []
order_falls_second_derivative_with_alternation = []
peaks_falls_with_alternation = [liste_tuple_points_caracteristics[i][0] for i in range(len(liste_tuple_points_caracteristics)) if liste_tuple_points_caracteristics[i][1] == 'fall']

for i in peaks_falls_with_alternation:
    order_falls_with_alternation.append(sd.isel(time=i).to_numpy())
    order_falls_second_derivative_with_alternation.append(ds_median_space_filter['second_derivative'][i].to_numpy())



#%% Cell 10 ##

## Display
# ========

# Link to the image storage folder

link = r'C:\Users\lucas\OneDrive\Bureau\ECL\UE\PAr\git_manteau_neigeux\Finse\modele_manteau_neigeux\images_2022_03_29'


# layout of the selected layers

def construction_couches(peaks_falls_with_alternation,sd) :
    ''''Creation of a matrix where each row corresponds to a layer. The top row corresponds to the last fallen layer'''
    matrice = np.empty((len(peaks_falls_with_alternation) + 1,len(sd)))  # Matrice de Nan avec une ligne supplementaire en haut etant la courbe en entiere
    matrice[0:len(peaks_falls_with_alternation)] = snow_height # a mettre si on ne veut pas de Nan dans la matrice
    matrice[-1] = snow_height

    for i in range(len(peaks_falls_with_alternation)-2) :
        matrice[i][liste_tuple_points_caracteristics[2*i][0]:liste_tuple_points_caracteristics[2*(i+1)][0]] = matrice[-1][liste_tuple_points_caracteristics[2*i][0]:liste_tuple_points_caracteristics[2*(i+1)][0]]
        for k in range(liste_tuple_points_caracteristics[2*(i+1)][0], len(sd)) :
            matrice[i][k] = matrice[-1][liste_tuple_points_caracteristics[2*(i+1)][0]]
    return(matrice)


def compaction(matrice,snow_density) : # TO DO
    matrice_bis = matrice
    # ????????????????????
    return(matrice_bis)

if 1 :
    snow_height = sd.to_numpy()
    snow_density = 6
    matrice = construction_couches(peaks_falls_with_alternation,sd)
    matrice = compaction(matrice,snow_density)

    # Vector of the legends = dates of appearance of the layers
    vecteur_dates = []

    if liste_tuple_points_caracteristics[0][1] == 'fall' :
        i=0
        while i < len(liste_tuple_points_caracteristics)-1 :
            vecteur_dates.append([str(list_dates_axis_abscisses[liste_tuple_points_caracteristics[i][0]]),str(list_dates_axis_abscisses[liste_tuple_points_caracteristics[i+1][0]])])
            i+=2
    if liste_tuple_points_caracteristics[0][1] == 'pic' :
        i=1
        while i < len(liste_tuple_points_caracteristics)-1 :
            vecteur_dates.append([str(list_dates_axis_abscisses[liste_tuple_points_caracteristics[i][0]]),str(list_dates_axis_abscisses[liste_tuple_points_caracteristics[i+1][0]])])
            i+=2

    # Affichage
    plt.figure()
    plt.plot(list_dates_axis_abscisses,matrice[0])
    plt.plot(list_dates_axis_abscisses,matrice[-1]) # list_dates_axis_abscisses
    plt.title("Time interval for the appearance of snowfalls")
    if len(peaks_falls_with_alternation) != 1 :
        for i in range(1,len(peaks_falls_with_alternation)) :
            plt.plot(list_dates_axis_abscisses,matrice[i],label=vecteur_dates[i-1]) # list_dates_axis_abscisses
        plt.legend()

    plt.savefig(link + 'final_result')
