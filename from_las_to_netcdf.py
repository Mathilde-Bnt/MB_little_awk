"""
Script to process Livox point clouds
Original work from S. Filhol and L. Girod
March 2021

Gaspard and Lucas
March 2022
"""

import numpy as np
import configparser, logging
#import openpylivox as opl
import glob, os
import pdal, json
import pandas as pd


#This Function rotate the point clouds, removes points with intensity below 20 and crops to the studied area
def ref_and_crop(z_range='[-5:0.5]',crop_corners='([-20, -1], [-5, 5])',path_to_data=r'C:\Users\Gaspard\Documents\Ecole\5_Centrale\Cours\PAr\donnees 2022\gaspard_laz'):


    file_list = glob.glob(path_to_data + '\laz_raw\*.laz')
    file_list = [os.path.normpath(i) for i in file_list]
    file_list.sort()



    i = 8055
    for file in file_list:

        pip_filter_json = json.dumps(
            {
                "pipeline":
                    [
                        file,
                        #rotation matrix
                        {
                            "type":"filters.transformation",
                            "matrix": '-0.34448594  0.93707407  0.05675957  2.51637959 -0.00583132  0.05832322 -0.9982807   0.35913649 -0.93877339 -0.34422466 -0.01462715  9.57211494 0. 0. 0. 1.'
                         },
                         #removes points with intensity below 20
                         {
                             "type": "filters.range",
                             "limits": "Z[-6:1],Intensity![0:20]"
                         },
                         #cropping
                        {
                                "type": "filters.crop",
                                "bounds": crop_corners
                            },
                        #creation of the new file
                        {
                            "type":"writers.las",
                            "filename": path_to_data + "\\laz_referenced_crop_not020\\" + file.split('\\')[-1]
                        }
                    ]
            })
        pipeline = pdal.Pipeline(pip_filter_json)
        pipeline.execute()
        print("fichier : "+ file.split('\\')[-1]+" traité")





'''
Conversion las/laz to .tif
'''

def laz_to_tif(GSD= 0.1, origin_x=-20,origin_y=-4.5,height=10,width=20,method='pdal',path_to_data=r'C:\Users\Gaspard\Documents\Ecole\5_Centrale\Cours\PAr\donnees 2022\gaspard_laz',delete_las=False,tif_to_zip=False):
#path_to_data=r'D:\gaspard_laz'

    file_list = glob.glob(path_to_data + '\\laz_referenced_crop_not020\\*.laz')
    file_list = [os.path.normpath(i) for i in file_list]
    file_list.sort()

    i = 8055

    for file in file_list:

        #attention où on split
        tst_data = pd.to_datetime(file.split('\\')[-1][:-4],format="%Y.%m.%dT%H-%M-%S")
        #if (tst_data.second + 60*tst_data.minute + 3600*tst_data.hour) % sampling_interval == 0:
            # Compute DEM with PDAL
        if method == 'pdal':
            pip_filter_json = json.dumps(
                {
                    "pipeline":
                        [
                            file,
                            {
                                "type": "writers.gdal",
                                "gdaldriver": "GTiff",
                                "output_type": "all",
                                "resolution": str(GSD),
                                "bounds": "([-20,0],[-4.5,4.5])",
                                "filename": path_to_data + "\\TIFs\\" +  file.split('\\')[-1][:-4] + ".tif"
                             }
                        ]
                })
            pipeline = pdal.Pipeline(pip_filter_json)
            pipeline.execute()
            i-=1

            print("reste :"+ str(i) +" fichier : "+ file.split('\\')[-1][:-4] + " traité")





























'''
Conversion tif to  netcdf
'''


# ENV rxr

import glob, argparse, datetime, configparser, logging
import os
import pandas as pd
import xarray as xr
from osgeo import gdal





def tif_to_ds(path_to_data=r'C:\Users\Gaspard\Documents\Ecole\5_Centrale\Cours\PAr\donnees_2022\gaspard_laz',filename_format='%Y%m%d.nc'):

    # list filename
    file_list = glob.glob(path_to_data + '\\TIFs\\*.tif')
    file_list = [os.path.normpath(i) for i in file_list]
    file_list.sort()

    # create dataframe of file metadata
    meta = pd.DataFrame({'fname':file_list})
    #extract timestamp from filename
    meta['tst'] = pd.to_datetime(meta.fname.apply(lambda x: x.split('\\')[-1][:-4]),format = '%Y.%m.%dT%H-%M-%S')

    i = len(meta.tst.dt.strftime('%m-%d').unique())

    # Create on netcdf file per day
    for date in meta.tst.dt.strftime('%m-%d').unique():
        # create time variable
        time_var = xr.Variable('time', meta.tst.loc[meta.tst.dt.strftime('%m-%d')==date])
        # open raster files in datarray
        geotiffs_da = xr.concat([xr.open_rasterio(i, engine = "rasterio") for i in meta.fname.loc[meta.tst.dt.strftime('%m-%d')==date]], dim=time_var)
        # drop all NaN values
        geotiffs_da = geotiffs_da.where(geotiffs_da!=-9999, drop=True)
        # rename variables with raster band names
        #pdb.set_trace()
        var_name = dict(zip(geotiffs_da.band.values, geotiffs_da.descriptions))
        geotiffs_ds = geotiffs_da.to_dataset('band')
        geotiffs_ds = geotiffs_ds.rename(var_name)

        # save to netcdf file
        fname_nc = path_to_data + '\\NC\\' + meta.tst.loc[meta.tst.dt.strftime('%m-%d')==date].iloc[0].strftime(filename_format)

        geotiffs_ds.to_netcdf(fname_nc, engine = "h5netcdf")
        print('Reste : '+ str(i) + 'File saved: ', fname_nc)
        i-=1
        # clear memory cache before next loop
        geotiffs_da = None
        geotiffs_ds = None





if __name__ == '__main__':
        #path to working directory where file will be created
        #path_to_data
        ref_and_crop(path_to_data = path_to_data)
        laz_to_tif()
