'''
cbfs.py
Functions to calculate the Common Base Function for a given field in a netCDF file. 
Calculation functions include:
    cbf_calc = calculates the cbf given a netCDF containing the single variable from a model and another netCDF containing the observations (reanalysis).

Given NetCDF files should contain data on the same lon-lat grid for the models and the observations.
Code uses the full grid and time period given, so to analyse data for a season and/or sub-domain (e.g. DJF NAO over North Atlantic), ensure given NetCDF files contain data on the desired sub-domain and for the given season only. These can be extracted using CDO or NCO before calling this script. 

The CBF analysis methodology is described in the appendix of Lee et al. 2018.
This code utilises the Iris package for handling structured data.

Written by Stephen Outten 2021
'''


import numpy as np
import iris
from eofs.iris import Eof
from scipy import stats
import sys


def load_data(fn,vn):
    '''
    Loads the data from a NetCDF with a given filename, fn, for a given variable name, vn.
    Returns error message and exits if data does not load correctly.
    Confirms resulting iris cube contains 3 dimensions with two of which being latitude and longitude.
    Identifies the time, lat, and lon coordinates and returns their names. 
    usage:  cube, coords = load_data(fn,vn)
        fn = file name with path for NetCDF file
        vn = variable name within NetCDF file
        cube = data loaded as cube
        coords = list containing names of the time, latitude and longitude coordinates in that order.
    '''
    # Load data, return error and exit if fail
    try: 
        cube = iris.load_cube(fn, vn)
    except Exception as e:
        print('The following error occurred while trying to load {}:'.format(fn))
        print(e)
        sys.exit()
    
    # Check data dimensions and extract coord names and dims
    coords = [coord.name() for coord in cube.coords()]
    if len(coords) != 3:
        print('Data field has more than three dimensions. Should have latitude, longitude, and time.')
        sys.exit()
    tcoord = [x for x in coords if x not in ['lat','latitude','lon','longitude']][0]
    ltcoord = [x for x in coords if x in ['lat','latitude']][0]
    lncoord = [x for x in coords if x in ['lon','longitude']][0]

    # Check latitude orientation
    if cube.coord(ltcoord).points[0] > cube.coord(ltcoord).points[1]:
        cube.data = np.flip(cube.data, cube.coord_dims(ltcoord)[0])
        cube.coord(ltcoord).points = np.flip(cube.coord(ltcoord).points, 0)
   
    # Guess latitude and longitude boundaries
    cube.coord(ltcoord).guess_bounds()
    cube.coord(lncoord).guess_bounds()

    return cube, [tcoord, ltcoord, lncoord]



def calc_grid_area(field, coords):
    '''
    Calculates the grid area and reduces it to 2D by removing the time coordinate.
    Usage: grid_areas = calc_grid_areas(cube, coords)
        cube = iris cube of data containing time, latitude, and longitude dimensions
        coords = list of names of coordinate in order of [time, lat, lon]
        grid_areas = 2D array of grid areas
    '''
    if field.coord_dims(coords[0])[0] == 1:
        grid_areas = iris.analysis.cartography.area_weights(field)[:,0,:]
    elif field.coord_dims(coords[0])[0] == 2:
        grid_areas = iris.analysis.cartography.area_weights(field)[:,:,0]
    else:
        grid_areas = iris.analysis.cartography.area_weights(field)[0,:,:]
    return grid_areas



def calc_eof(cube, coords, nmode=1):
    '''
    Calculates the EOF on a given cube for a give mode.
    Usage: eof, perc = calc_eof(cube, coords, nmode)
        cube = iris cube of data containing time, latitude, and longitude dimensions 
        coords = list of names of coordinate in order of [time, lat, lon]
        nmode = number of the mode to calculate (max=10, default=1)
        eof = iris cube of pattern of eof for chosen mode
        perc = percentages of variability explained by chosen mode
    '''
    if nmode > 10:
        print('Function only calculates the normalised EOF for a given mode up to 10. Requested EOF mode is {}'.format(nmode))
        print('Please select a lower EOF mode.')
        sys.exit()
    grid_areas = calc_grid_area(cube,coords)
    cube_anom = cube - cube.collapsed(coords[0], iris.analysis.MEAN)
    eof_solver = Eof(cube_anom, weights='area')
    eof = eof_solver.eofsAsCovariance(neofs=nmode)     
    perc = eof_solver.varianceFraction(nmode).data[nmode-1]*100
    return eof[nmode-1,:,:], perc



def calc_cbf(cube, eof, coords):
    '''
    Calculates the CBF on a given cube using a given EOF.
    Usage: cbf, perc = calc_cbf(cube, eof, coords)
        cube = iris cube of data containing time, latitude, and longitude dimensions
        eof = 2D iris cube of normalised EOF pattern calculated from observations
        coords = list of names of coordinate in order of [time, lat, lon]
        cbf = numpy array containing pattern of CBF
        perc = percentages of variability explained by CBF 
    '''
    a,b,c = cube.shape
    grid_areas = calc_grid_area(cube,coords)
    cube_anom = cube - cube.collapsed(coords[0], iris.analysis.MEAN)
    field = np.asarray(cube_anom.data)
    eofn = np.asarray(((eof - eof.collapsed([coords[2], coords[1]], iris.analysis.MEAN, weights=grid_areas)) / eof.collapsed([coords[1], coords[2]],iris.analysis.STD_DEV)).data).flatten()
    pc = np.nan * np.ones(a)
    for i in range(a):
        pc[i] = np.dot(field[i,:,:].flatten(),eofn)
    slope = np.nan * np.ones((b,c))
    intercept = np.nan * np.ones((b,c))
    for i in range(b):
        for j in range(c):
            slope[i,j], intercept[i,j], _, _, _ = stats.linregress(pc, field[:,i,j])
    weights = grid_areas/np.sum(grid_areas)
    cbf_var = np.sum(np.var(np.moveaxis(slope[..., None] * pc, -1, 0), 0) * weights)
    cbf_perc = 100 * cbf_var/(np.sum(np.var(cube.data, 0) * weights))
    cbf = slope * np.std(pc)

    return cbf, cbf_perc



def cbfs(model_fn, model_vn, obs_fn, obs_vn, nmode=1):
    '''
    Calculates the Common Base Function as based on the work of Lee et al. 2018.
    Note: Variables in model NetCDF and observations NetCDF must be the same size. 
    Usage: cbf,cbf_perc,eof,eof_perc = cbfs(model_file_path, model_variable_name, observations_file_path, observations_variable_name, nmode)
        model_file_path = Path to NetCDF file containing variable of interest
        model_variable_name = Name of variable of interest in the model file
        observations_file_path = Path to NetCDF file containing observations of variable of interest
        observations_variable_name = Name of variable in the observation file
        nmode = number of the EOF mode to use as basis for CBF, must be 10 or less.
        cbf = common base function pattern for the chosen mode of the EOF
        cbf_perc = percentage of variability explained by mode in the CBF
        eof = iris cube of pattern of eof for chosen mode
        eof_perc = percentages of variability explained by chosen mode in EOF
    '''
    # Check request mode is 10  or less.
    if nmode > 10:
        print('Module only calculates CBFs based on EOF modes up to 10. Requested EOF mode is {}'.format(nmode))
        print('Please select a lower EOF mode to use as the basis of the CBF.')
        sys.exit()
        
    # Load observations and calculate EOF for requested mode
    odat, ocoord = load_data(obs_fn, obs_vn)
    eof,eof_perc = calc_eof(odat, ocoord, nmode)

    # Check how many models are requested and load and calculate CBF for each in turn
    if type(model_fn) == str:
        mdat, mcoord = load_data(model_fn, model_vn)
        if (mdat.coord(mcoord[1]).shape[0], mdat.coord(mcoord[2]).shape[0]) != (odat.coord(ocoord[1]).shape[0], odat.coord(ocoord[2]).shape[0]):
            print('Observations and model have different latitude and longitudes. They must be on the same spatial grid. Recommend re-gridding with CDO.')
            sys.exit()
        cbf,cbf_perc = calc_cbf(mdat, eof, mcoord)
    elif type(model_fn) == list:
        _,b,c = odat.shape
        a = len(model_fn)
        cbf = np.nan * np.ones((a,b,c))
        cbf_perc = np.nan * np.ones((a))
        for ii in range(a):
            mdat, mcoord = load_data(model_fn[ii], model_vn)
            if (mdat.coord(mcoord[1]).shape[0], mdat.coord(mcoord[2]).shape[0]) != (odat.coord(ocoord[1]).shape[0], odat.coord(ocoord[2]).shape[0]):
                print('Observations and model have different latitude and longitudes for model_fn {}. Models must use same spatial grid as observations. Recommend re-gridding model data with CDO.'.format(model_fn[ii]))
                print('Skipping this model.')
                continue
            cbf_temp,perc_temp = calc_cbf(mdat, eof, mcoord)
            cbf[ii,:,:] = cbf_temp.copy()
            cbf_perc[ii] = perc_temp
    else:
        print('model_fn must either be a string containing a single filename including path or a list containing strings of filenames including paths.')
        print('Given model_fn is of type {}, so this ends here.'.format(type(model_fn)))
        sys.exit()
    
    return cbf, cbf_perc, eof, eof_perc
