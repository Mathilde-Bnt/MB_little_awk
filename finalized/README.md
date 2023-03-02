This repository contains all the necessary tools to analyze the output of a time-lapse lidar.
The signal should already have been converted to netcdf format, using the Little Awk tools.

## Description of files

- ddensity_ml.f95 and associated cpython file (ml stands for multi-layer) contain the compaction subroutine taken from SnowModel (2006), which updates the density of each layer at simulation timestamps.
- snowtemp_ml.f95 and associated cpython file contain the temperature subroutines taken from SnowModel (2006), which updates the temperature of each layer at simulation timestamps.

- little_awk_functions.py contains only functions that are usefull to analyze lidar data and run snowpack simulations. This file calls the ddensity_ml and snowtemp_ml routines.
- parameters.py is a file in which parameters and data storage structures are defined before launching simulations. This file uses functions from little_awk_functions.py.

- all .ipynb files are codes that simulate a snowpack or analyze the lidar data in some way. Their individual functions are described below. They usually start by running little_awk_functions.py and parameters.py. Some have additional functions and parameters of their own.

## Typical pipeline once the netcdf files are created from Lidar data

1. Apply the tools
- Make sure all timestamps contain data > notebook 01
- Create the netcdf files containing the filled data > notebook 00
- Run a first simulation > notebook 02
- Predict 2023 snowpack ineer structure > notebook 03

2. Analyze results
- Identify "marker" layers using the meteorological data > notebook 11
- Adjust parameters of simulation (sensitivity analysis and A1-A2 calibration) > notebooks 12 and 13
- Update parameters to be used > parameters.py
- Run validation scripts to evaluate the correctness of outputs > notebook 14
- Make a precipitation forcing file for Crocus > notebook 15

3. Make figures
- Illustrate the detection of snow events > notebook 21
- Illustrate the general match of simulated layers and lidar signal > notebook 02
- Make a simulated snowpit to compare to experimental images > notebook 23

## Compilation code Fortran
```sh
f2py -c snowtemp_ml.f95 -m snowtemp
```

This creates a module that is readable in Python by running
```sh
import snowtemp
```
