'''
Run example script to test and develop snow layers picking
S. Filhol, July 2022

TODO:
- grab one day of summer data to create summer surface
- explore how to detect peaks in derivatives using dask:
    - extract peaks and identify snowfall start and end
        - find faster method than scipy.signal.find_peak()
        - find local extrema and minima with scipy.signal.argrelextrema()
        - compute z-score of derivative and apply a threshold to find peaks


'''
from little_awk import pcl_to_netcdf as pn
from little_awk import extract_snowpack as es
import matplotlib.pyplot as plt
import os, glob
import dask
import xarray as xr
proj_dir = 'data'
#=========  0. Converting Livox bin to Netcdf =======
pn.convert_bin_to_las(project_dir=proj_dir, deleteBin=True)
pn.rotate_crop_filter_pcl(project_dir=proj_dir)
pn.las_to_tif(project_dir=proj_dir)
pn.tif_to_netcdf(project_dir=proj_dir)

# logic to remove all intermediate files
for folder in ['las_raw', 'las_clean', 'TIFs']:
    flist = glob.glob(proj_dir + os.sep + folder + os.sep + '*')
    for file in flist:
        os.remove(file)


#=========  2. Summer and snow depth map =======
ds = xr.open_mfdataset('data/netcdfs/*', engine='h5netcdf', chunks={'x': 20,'y': 20})
es.define_summer_surface(ds, '2021-07-01', '2021-07-02')
es.remove_outliers(ds)
es.median_spacetime_filtering(ds)
es.to_depth_nc(ds, 'data/sdepths.nc')

#=========  2. Extracting Snow Layers =======
ds = None
ds = xr.open_mfdataset('data/sdepths.nc', engine='h5netcdf', chunks={'x': 20,'y': 20})
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

