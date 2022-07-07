import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import dask
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
    ds['summer_surf'] = ds['mean'].sel(time = slice(start,end)).mean(dim='time')
    # add here code to include metadata to variable

def remove_outliers(ds, lower_thresh=0.1, upper_thresh=4):
    '''
    Function to replace potential abherant outliers relative to the summer surface with NaNs
    Args:
        ds: dataset
        thresholds: lower and upper thresholds, in same unit as ds['mean']

    Returns:

    '''
    ds = ds.where(ds['mean'] < ds.summer_surface + upper_thresh, np.nan)
    ds = ds.where(ds['mean'] > ds.summer_surface - lower_thresh, np.nan)

def median_spacetime_filtering(ds, time_window=11, x_span=11, y_span=11):
    '''
    Function to apply median filtering in time and space
    Args:
        ds: clean data dataset
        time_window: time windo on which to apply median filter [index]
        x_span: x-span on which to apply filtering, unit [index]
        y_span: y-span on which to apply filtering, unit [index]

    Returns:

    '''
    ds['snow_surface'] = ds['mean'].rolling(time=time_window, center=True).median()
    ds['snow_surface'] = ds['snow_surface'].rolling({'x': x_span, 'y': y_span}, center=True).median()

def surface_to_depth(ds, fname):
    '''
    Function to store summer reference and snow surfaces to an independent netcdf for further processing with dask.
    Args:
        ds:
        fname:

    Returns:

    '''
    ds['snow_depth'] = ds.snow_surface - ds.summer_surface
    ds[['snow_depth']].to_netcdf('test.nc', engine='h5netcdf',
                                          encoding={
                                              'snow_depth': {'dtype': 'float32', 'zlib': True}
                                          })
    print('Snow depth maps saved in file {}'.format(fname))




'''
next steps:
    - compute first and second derivative in respect to time
    - extract peaks and identify snowfall start and end
        - find faster method than scipy.signal.find_peak()
        - find local extrema and minima with scipy.signal.argrelextrema()
        - compute z-score of derivative and apply a threshold to find peaks
    - remove effect due to dune migration
    - 
'''



if __name__ == '__main__':
    ds = xr.open_mfdataset(path, engine='h5netcdf', chunks={'x': 10,'y': 10})
    define_summer_surface(ds, '2022-02-03', '2022-02-04')
    remove_outliers(ds)
    median_spacetime_filtering(ds)
    surface_to_netcdf(ds, 'sssrf.nc')
    ds = None
    ds = xr.open_mfdataset('sssrf.nc', engine='h5netcdf', chunks={'x': 10,'y': 10})
    grad1 = dask.array.gradient(ds['snow_surface'], axis=0, edge_order=1)
    grad2 = dask.array.gradient(grad1, axis=0, edge_order=1)

    gradients = ds.drop_vars(list(ds.keys()))
    gradients['first'] = (['time', 'y', 'x'], grad1)
    gradients['second'] = (['time', 'y', 'x'], grad2)


    fig, ax = plt.subplots(3,1,sharex=True)
    ax[0].plot(ds['surfs'].isel(y=35,x=50))
    ax[1].plot(a[:,35,50].compute())
    ax[2].plot(b[:,35,50].compute())
    plt.show()