# little_awk
Package for processing time-lapse lidar, and etract snowpack stratigraphy

## Installation
Steps to create a Python environemnt valid to work with Little_awk library
'''
mamba create -n "awky" python=3.7 ipython
conda activate awky
conda install -c conda-forge libgdal
conda install -c conda-forge gdal
conda install -c conda-forge pdal python-pdal
pip install xarray matplotlib pandas numpy scipy plyfile geopandas pyyaml rasterio dask
pip install laspy==1.7


git clone https://github.com/ArcticSnow/OpenPyLivox
pip install -e OpenPyLivox
'''
**Valid on Ubuntu 20.04**


**TODO**:
- [ ] Implement new snow detection algorithm
- [ ] Figure out compaction/temperature model
  - use existing model tool such as snowpack.f from snowmodel or a simplified version of Crocus
- [ ] implement method to visualize dunes and sastrugi over the whole season.
  - figure out what could be a good reference surface that would highlight the new dunes/sastrugi
- [ ] bring all configuration into a config file in YAML
- [ ] add compatibility code with config file
- [ ] 