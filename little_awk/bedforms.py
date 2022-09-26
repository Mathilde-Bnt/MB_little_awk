'''
Set of tools and function to extract bedforms

'''

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import dask
import os

import extract_snowpack as es

def plot_surface(ref_surf, ds):