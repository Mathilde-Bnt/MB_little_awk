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
pip install xarray matplotlib pandas numpy scipy plyfile geopandas pyyaml

git clone https://github.com/ArcticSnow/OpenPyLivox
pip install -e OpenPyLivox
'''
**Valid on Ubuntu 20.04**


**TODO**:
- [ ] bring all configuration into a config file in YAML
- [ ] add compatibility code with config file
- [ ] 