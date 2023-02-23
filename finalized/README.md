This repository contains all the necessary tools to analyze the output of a time-lapse lidar.
The signal should already have been converted to netcdf format, using the Little Awk tools.

Typical pipeline once these files are ready:

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

3. Make figures
- Illustrate the detection of snow events > notebook 21
- Illustrate the general match of simulated layers and lidar signal > notebook 02
- Make a simulated snowpit to compare to experimental images > notebook 23

## conpilation code fortran
```sh
f2py -c snowtemp_ml.f95 -m snowtemp_test
```

