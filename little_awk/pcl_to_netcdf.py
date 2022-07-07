"""
Script to process Livox point clouds
Original work from S. Filhol and L. Girod, March 2021
Modified by Gaspard and Lucas, March 2022
Modified S. Filhol, June 2022

TODO:
- add functionality to delete .las and tif file once converted.
- add metadata field in netcdf file
- possibility to speed up this process parallelizing the file processing on multiple core with
    'from multiprocessing.dummy import Pool as ThreadPool'


"""

import glob
import os
import json
import pandas as pd
import pdal
from tqdm import tqdm
import xarray as xr
import openpylivox as opl

def make_directory(project_dir='/home/data'):
    os.mkdir(project_dir + os.sep + 'bin')
    os.mkdir(project_dir + os.sep + 'las_raw')
    os.mkdir(project_dir + os.sep + 'las_clean')
    os.mkdir(project_dir + os.sep + 'TIFs')
    os.mkdir(project_dir + os.sep + 'netcdfs')

def convert_bin_to_las(project_dir='/home/data',
                       bin_folder='bin',
                       file_pattern='*.bin',
                       las_folder='las_raw',
                       output_path='./',
                       deleteBin = False):
    """
    Function to convert a list of Livox binary files from .bin to .las
    Args:
        bin_flist (str list): list of .bin filepath to convert
        output_path (str): path where to export file
        deleteBin (bool): delete bin file, True/False
    """
    file_list = glob.glob(project_dir + os.sep + bin_folder + os.sep + file_pattern)
    file_list.sort()
    if file_list.__len__() == 0:
        print('WARNING: No file found')
        return

    for file in file_list:
        opl.convertBin2LAS(file, deleteBin=deleteBin)
        if os.path.isfile(file + '.las'):
            os.rename(file + '.las', project_dir + os.sep + las_folder + os.sep + file.split('/')[-1][:-4] + '.las')
    return


def rotate_crop_filter_pcl(z_range=[-6, 1],
                           crop_extent='([-20, -1], [-5, 5])',
                           rotation='-0.34448594  0.93707407  0.05675957  2.51637959 -0.00583132  0.05832322 -0.9982807   0.35913649 -0.93877339 -0.34422466 -0.01462715  9.57211494 0. 0. 0. 1.',
                           project_dir='to/my/proj',
                           las_raw_folder='las_raw',
                           las_clean_folder='las_clean',
                           file_pattern='*.las'):
    '''
    Functino to rotate the point clouds, removes points with intensity below 20 and crops to the studied area
    Args:
        z_range: range to crop along the Z-axis [zmin, zmax]
        crop_extent: extent to crop pcl ([xmin, xmax],[ymin, ymax])
        rotation: rotation matrix
        project_dir:
        file_pattern:

    Returns:
        save processed file to folder las_clean/
    '''
    print('---> rotate, crop, and filter intensity of point cloud')
    file_list = glob.glob(project_dir + os.sep + las_raw_folder + os.sep + file_pattern)
    file_list = [os.path.normpath(i) for i in file_list]
    file_list.sort()

    for file in tqdm(file_list):
        pip_filter_json = json.dumps(
            {
                "pipeline":
                    [
                        file,
                        # rotation matrix
                        {
                            "type": "filters.transformation",
                            "matrix": rotation
                         },
                         # removes points with intensity below 20
                         {
                             "type": "filters.range",
                             "limits": "Z[{}:{}],Intensity![0:20]".format(z_range[0], z_range[1])
                         },
                         # cropping in X and Y
                        {
                                "type": "filters.crop",
                                "bounds": crop_extent
                            },
                        # write to new file
                        {
                            "type": "writers.las",
                            "filename": project_dir + os.sep + las_clean_folder + os.sep + file.split('/')[-1]
                        }
                    ]
            }
        )
        pipeline = pdal.Pipeline(pip_filter_json)
        pipeline.execute()
        print("file : " + file.split(os.sep)[-1] + " processed")


def las_to_tif(resolution= 0.1,
               bounds='([-20,0],[-4.5,4.5])',
               project_dir='to/project',
               file_pattern='*.las'):
    '''
    Function to extract geotiff from point clouds

    Args:
        resolution: ground resolution in the same unit as point cloud
        bounds: bounds of raster ong XY axis ([],[})
        project_dir: project directory
        file_pattern: file pattern to convert

    Returns:
        save geotif file in folder '/tifs'
    '''
    print('---> Converting las to tif')
    file_list = glob.glob(project_dir + os.sep + 'las_clean' + os.sep + file_pattern)
    file_list = [os.path.normpath(i) for i in file_list]
    file_list.sort()

    for file in tqdm(file_list):
        # Extract timestamp from filename
        tst_data = pd.to_datetime(file.split(os.sep)[-1][:-4], format="%Y.%m.%dT%H-%M-%S")
        #if (tst_data.second + 60*tst_data.minute + 3600*tst_data.hour) % sampling_interval == 0:

        # Compute DEM with PDAL
        pip_filter_json = json.dumps(
            {
                "pipeline":
                    [
                        file,
                        {
                            "type": "writers.gdal",
                            "gdaldriver": "GTiff",
                            "output_type": "all",
                            "resolution": str(resolution),
                            "bounds": bounds,
                            "filename": project_dir + os.sep + "TIFs" + os.sep + file.split(os.sep)[-1][:-4] + ".tif"
                         }
                    ]
            }
        )
        pipeline = pdal.Pipeline(pip_filter_json)
        pipeline.execute()


def tif_to_netcdf(project_dir='to/my/project',
              file_pattern='*.tif',
              output_format='%Y%m%d.nc',
              ):
    '''
    Function to compile geotiff to daily netcdf file
    Args:
        project_dir: project directory
        file_pattern:
        output_format:

    Returns:
        save netcdf file

    Todo:
    - add metadata field based on YAML file
    - replace deprecated xarray rasterio function by rioxarray

    '''
    print('---> compile geotiffs into daily netcdfs')
    # list filename
    file_list = glob.glob(project_dir + os.sep + 'TIFs' + os.sep + file_pattern)
    file_list = [os.path.normpath(i) for i in file_list]
    file_list.sort()

    # create dataframe of file metadata
    meta = pd.DataFrame({'fname': file_list})
    #extract timestamp from filename
    meta['tst'] = pd.to_datetime(meta.fname.apply(lambda x: x.split(os.sep)[-1][:-4]), format='%Y.%m.%dT%H-%M-%S')

    i = len(meta.tst.dt.strftime('%m-%d').unique())

    # Create on netcdf file per day
    for date in meta.tst.dt.strftime('%m-%d').unique():
        # create time variable
        time_var = xr.Variable('time', meta.tst.loc[meta.tst.dt.strftime('%m-%d') == date])
        # open raster files in datarray
        geotiffs_da = xr.concat([xr.open_rasterio(i, engine="rasterio") for i in meta.fname.loc[meta.tst.dt.strftime('%m-%d') == date]], dim=time_var)
        # drop all NaN values
        geotiffs_da = geotiffs_da.where(geotiffs_da != -9999, drop=True)
        # rename variables with raster band names
        #pdb.set_trace()
        var_name = dict(zip(geotiffs_da.band.values, geotiffs_da.descriptions))
        geotiffs_ds = geotiffs_da.to_dataset('band')
        geotiffs_ds = geotiffs_ds.rename(var_name)

        # save to netcdf file
        fname_nc = project_dir + os.sep + 'netcdfs' + os.sep + meta.tst.loc[meta.tst.dt.strftime('%m-%d') == date].iloc[0].strftime(output_format)

        geotiffs_ds.to_netcdf(fname_nc, engine="h5netcdf")
        print('Remaining : ' + str(i) + 'File saved: ', fname_nc)
        i -= 1
        # clear memory cache before next loop
        geotiffs_da = None
        geotiffs_ds = None


if __name__ == '__main__':
    import yaml
    import argparse
    import logging

    '''
    Organize project with the following folders:
        - myproject/
            - laz_raw/
            - las_clean/
            - TIFs/
            - netcdfs/
    '''


    logging.info('Parse arguments')
    parser = argparse.ArgumentParser()
    parser.add_argument('--project_directory', '-dir', help='Path to project directory', default='./myproject')
    parser.add_argument('--awk_config', '-cf', help='Path to config file of little_awk project', default='config.yml')
    args = parser.parse_args()

    logfile = args.project_dir + 'awk_proj.log'
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s %(levelname)s %(threadName)s %(name)s %(message)s",
                        # datefmt='%m-%d %H:%M',
                        filename=logfile,
                        filemode='w')
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    logging.getLogger('').addHandler(console)
    logging.info(args)

    # Open config file into a dictionnary
    with open(args.project_dir + args.awk_config, 'r') as file:
        conf = yaml.safe_load(file)

    rotate_crop_filter_pcl(project_dir=args.project_dir,
                           z_range=conf['pcl_to_netcdf']['Z_range'],
                           crop_extent=conf['pcl_to_netcdf']['XY_extent'],
                           rotation=conf['pcl_to_netcdf']['rotation'],
                           file_pattern='*.las'
                           )

    las_to_tif(project_dir=args.project_dir,
               resolution= conf['pcl_to_netcdf']['resolution'],
               bounds=conf['pcl_to_netcdf']['XY_extent'],
               file_pattern='*.las')

    tif_to_netcdf(project_dir=args.project_dir,
                  file_pattern='*.tif',
                  output_format='%Y%m%d.nc',)


