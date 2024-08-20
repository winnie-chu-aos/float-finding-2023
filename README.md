"main_script.m" contains the bulk of the code (both correlation length analysis and float finding algorithm)
    -> since this is the most well-commented version of the code, it is recommended that readers start here to understand the analysis
"DCall_1deg_update.mat" contains the correlation lengths with 1 deg lat/lon resolution that were reported in the paper
"find_floats_USE.m" is the last working version of the function to find floats iteratively using the float finding algorithm
"float_boat.m" is the version of the float finding algorithm where the available float locations are limited to the ship track
"plot_float_coverage.m" is the function used to plot the float coverage
