'''
Set of tools to derive snowpack depth and layers from gridded lidar scans
S. Filhol 2022

'''

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
    ds['summer_surface'] = ds['mean'].sel(time = slice(start, end)).mean(dim='time')
    ds.summer_surface.attrs = {'units':'m', 'standard_name':'summer_surface',
                               'long_name':'Summer surface'}
    print(f'---> Summer surface defined based on scans from {start} to {end}')


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
    print(f'---> Outliers removed')

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
    print(f'---> Median filtering in time with a window of{time_window}')
    ds['snow_surface'] = ds['snow_surface'].rolling({'x': x_span, 'y': y_span}, center=True).median()
    print(f'---> Median filtering in space with a window [{x_span}, {y_span}]')

def to_depth_nc(ds, fname):
    '''
    Function to store summer reference and snow surfaces to an independent netcdf for further processing with dask.
    Args:
        ds:
        fname:

    Returns:

    '''

    # Add logic to split into multiple netcdf file if too big. For instance use timestamp. One file per month or so.
    ds['snow_depth'] = ds.snow_surface - ds.summer_surface
    ds.snow_depth.attrs = {'units':'m', 'standard_name':'snow_depth',
                               'long_name':'Snow depth'}
    ds[['snow_depth', 'summer_surface']].to_netcdf('test.nc', engine='h5netcdf',
                                          encoding={
                                              'snow_depth': {'dtype': 'float32', 'zlib': True},
                                              'summer_surface': {'dtype': 'float32', 'zlib': True}
                                          })

    print(f'---> Snow depth maps saved to file {fname}')


def compute_time_derivatives(ds):
    '''
    Function to compute 1st and 2nd derivative of the snow surface height along the time axis
    Args:
        ds (xarray dataset): dataset containing snow surface data organized in ['time', 'y', 'x']

    Returns:
        A xarray dataset with the first and second order derivative
    '''
    grad1 = dask.array.gradient(ds['snow_surface'], axis=0, edge_order=1)
    grad2 = dask.array.gradient(grad1, axis=0, edge_order=1)

    gradients = ds.drop_vars(list(ds.keys()))
    gradients['first'] = (['time', 'y', 'x'], grad1)
    gradients['second'] = (['time', 'y', 'x'], grad2)
    print('---> Time derivatives computed')
    return gradients


def to_netcdf(ds, path, fname='tmp.nc', n=16):
    ''''
    Function to store dataset to netcdf with compression
    Args:
            ds (dataset): dataset to store as netcdf
            path (str): path where to save netcdf
            fname (str): name of netcdf file to save
            n (int): number of digits to account for during compression. Used in computing scaling factor

    '''
    def compute_scaling_and_offset(da, n=16):
        """
        Compute offset and scale factor for int conversion.

        Args:
            da (dataarray): of a given variable
            n (int): number of digits to account for
        """
        vmin = float(da.min().values)
        vmax = float(da.max().values)

        # stretch/compress data to the available packed range
        scale_factor = (vmax - vmin) / (2 ** n - 1)
        # translate the range to be symmetric about zero
        add_offset = vmin + 2 ** (n - 1) * scale_factor

        return scale_factor, add_offset

    encod_dict = {}
    # loop through each variable of the dataset
    for var in list(ds.keys()):
        scale_factor, add_offset = compute_scaling_and_offset(ds[var], n=10)
        encod_dict.update({var:{"zlib": True,
                               "complevel": 9,
                               'dtype':'int16',
                               'scale_factor':scale_factor,
                               'add_offset':add_offset}})
    ds.to_netcdf(path + fname, encoding=encod_dict)
    print('---> File {} saved'.format(fname))

'''
next steps:
    - compute first and second derivative in respect to time
    - extract peaks and identify snowfall start and end
        - find faster method than scipy.signal.find_peak()
        - find local extrema and minima with scipy.signal.argrelextrema()
        - compute z-score of derivative and apply a threshold to find peaks
    - remove effect due to dune migration
    
Notes:
    Find a simple but robust technique to identify peaks in derivative signal that minimize computation.
'''



if __name__ == '__main__':
    path = '/home/simonfi/github/little_awk/example/data/netcdfs/'
    nc_filepattern = '20*.nc'
    ds = xr.open_mfdataset(path + nc_filepattern, engine='h5netcdf', chunks={'x': 10,'y': 10})
    define_summer_surface(ds, '2021-07-01', '2021-07-02')
    remove_outliers(ds)
    median_spacetime_filtering(ds)
    to_netcdf(ds, '.', 'sssrf.nc')
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