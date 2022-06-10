import xarray as xr
import numpy as np
import os

def define_summer_surface(ds, start, end):
    '''
    Function to define summer surface by averaging along time axis a series of summer scans
    Args:
        ds: dataset to use
        start: start date
        end: end date

    Returns:
        add summer surface to dataset

    '''
    ds['summer_surf'] = ds['mean'].sel(time = slice(start,end)).mean(dim = 'time')
    # add here code to include metadata to variable

def remove_outliers(ds, lower_thresh=0.1, upper_thresh=4):
    '''
    Function to replace potential abherant outliers relative to the summer surface with NaNs
    Args:
        ds: dataset
        thresholds: lower and upper thresholds, in same unit as ds['mean']

    Returns:

    '''
    ds['mean'] = ds.where(ds['mean'] < ds.summer_surf + upper_thresh)
    ds['mean'] = ds.where(ds['mean'] > ds.summer_surf - lower_thresh)

def median_spacetime_filtering(ds, time_window=11, x_span=11, y_span=11):
    '''
    TO BE CHECKED!!!!
    Function to apply median filtering in time and space
    Args:
        ds: clean data dataset
        time_window: time windo on which to apply median filter [index]
        x_span: x-span on which to apply filtering, unit [index]
        y_span: y-span on which to apply filtering, unit [index]

    Returns:

    '''
    ds['surfs'] = ds['mean'].rolling(time=time_window, center=True).median(dim = 'time')
    ds['surfs'] = ds['surfs'].rolling({'x':x_span,'y':y_span}, center=True).median(dim = ['x','y'])

'''
next steps:
    - compute first and second derivative in respect to time
    - extract peaks and identify snowfall start and end
        - find faster method than scipy.signal.find_peak()
    - remove effect due to dune migration
    - 
'''



if __name__ == '__main__':
    ds = xr.open_mfdataset(path)