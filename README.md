"main_script.m" contains the bulk of the code (both correlation length analysis and float finding algorithm)\
  - Since this is the most well-commented version of the code, it is recommended that readers start here to understand the analysis.
  - "gft.m" is a helper function that can be downloaded from https://github.com/mmazloff/gft

"DCall_1deg_update.mat" contains the correlation lengths with 1 deg lat/lon resolution that were reported in the paper\
  - This file has matrix `DCL` with dimensions (longitude x latitude x 6 output variables).
  - The 6 output variables are semimajor length, semiminor length, angle of rotation, $r^2$, prior zonal length, and prior meridional length.
  - `DCL` is called in main_script.m (line 228-230) to plot the paper figures, but can be replaced with `QQ` if you have run the iterative loop (lines 83-205).

"find_floats_USE.m" is the last working version of the function to find floats iteratively using the float finding algorithm\
"float_boat.m" is the version of the float finding algorithm where the available float locations are limited to the ship track\
"plot_float_coverage.m" is the function used to plot the float coverage
