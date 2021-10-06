# CommonBasisFunction
Python module for calculating Common Basis Function from NetCDF files.
This is based on the outline of the CBF method from Lee et al. 2018.

The module requires the following packages:  iris, numpy, scipy, eofs

The module's primary function has the following usage:

    cbf,cbf_perc,eof,eof_perc = cbfs.cbfs(model_file_path, model_variable_name, observations_file_path, observations_variable_name, nmode)
  
  Inputs:
  
    model_file_path             = Path to NetCDF file containing variable of interest
    model_variable_name         = Name of variable of interest in the model file
    observations_file_path      = Path to NetCDF file containing observations of variable of interest
    observations_variable_name  = Name of variable in the observation file
    nmode (optional)            = number of the EOF mode to use as basis for CBF, must be 10 or less.

  Outputs:
  
    cbf       = common base function pattern for the chosen mode of the EOF
    cbf_perc  = percentage of variability explained by mode in the CBF
    eof       = iris cube of pattern of eof for chosen mode
    eof_perc  = percentages of variability explained by chosen mode in EOF


The module is contained in cbfs.py, this is the only file in this repository required to use this module. 

A test script and example data files are provided for checking the CBFS module is working correctly. 
The module can be tested by running 

    test_cbfs.py

This test script requires the following packages: iris, numpy, scipy, eofs, matplotlib, cartopy

The data file in the data directory are present only for the test_cbfs.py script and are no necessary for using the module. 

Files in Data directory include mean sea level pressure for the period of 1900-2010 over the domain of 20N-90N and 0-360E.

The files are in NetCDF format and are for 

    - ERA 20th Century reanalysis
    - NorESM2-LM 
    - NorESM2-MM
    - NorCPM1

The NorESM/NorCPM files are taken from the first ensemble members of the CMIP6 Historical runs for the respective version of the Norwegian Earth System Model. 
They are provided here purely for example purposes.


