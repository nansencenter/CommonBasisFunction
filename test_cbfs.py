'''
test_cbfs.py
Calls the cbfs module with the approriate inputs.
Used for testing the functions in the module.
This code requires iris, scipy, and numpy.

Module assumes the NetCDF files given are the correct season and region.
Example file are for ERA20th Century reanalysis and the CMIP6 Historical runs for three variants of the NorESM model.
The field in the example files is mean sea level pressure and the domain covers 20N-90N for all longitudes.
Given these files, the modules calculates the NAM (AO) pattern.

There is a flag for single model runs.
If loading a single model, this test all the sub-functions in the module directly.
If loading multiple models, this uses the primary cbfs function.
In both cases, this code plots the EOF from ERA20C and the CBF for NorESM2-LM.

There is a flag toturn off the plotting example.
Plotting requires matplotlib and cartopy.


Written by Stephen Outten October 2021.
'''

# Single model mode flag
single_model = True
iplot = True  

# ************************* Do Not Modify Below This Line **********************

import numpy as np
import sys
import warnings
import cbfs 
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib
import iris.plot as iplt

warnings.filterwarnings('ignore')

FolderName = 'Data/'
seasons = ['djf', 'mam', 'jja', 'son']
iseas = 1
obsFileName = 'era20c_psl_1900_2010_nh_'
obs_vn = 'sp'
model_vn = 'psl'
imod = 2    
modelFileName = ['psl_Amon_MRI-ESM2-0_historical_r1i1p1f1_gn_190001-201012_',
                'psl_Amon_NorCPM1_historical_r1i1p1f1_gn_190001-201012_',
                'psl_Amon_NorESM2-LM_historical_r1i1p1f1_gn_190001-201012_',
                'psl_Amon_NorESM2-MM_historical_r1i1p1f1_gn_190001-201012_']



# Set up a filename for observations and a model.
obs_fn = FolderName + obsFileName + seasons[iseas] + '.nc'
if single_model:
    model_fn = FolderName + modelFileName[imod] + seasons[iseas] + '.nc'
else:
    model_fn = [FolderName + modelFileName[imod-1] + seasons[iseas] + '.nc',
                FolderName + modelFileName[imod] + seasons[iseas] + '.nc',
                FolderName + modelFileName[imod+1] + seasons[iseas] + '.nc']


# Load data and perform EOF and CBF
if single_model:
    odat,ocoord = cbfs.load_data(obs_fn,obs_vn)
    mdat,mcoord = cbfs.load_data(model_fn,model_vn)
    eof,eof_perc = cbfs.calc_eof(odat,ocoord,nmode=1)
    cbf,cbf_perc = cbfs.calc_cbf(mdat,eof,mcoord)
else:
    cbf,cbf_perc,eof,eof_perc = cbfs.cbfs(model_fn, model_vn, obs_fn, obs_vn, 1)
    odat,ocoord = cbfs.load_data(obs_fn,obs_vn)


# Show percentage variance explained
if single_model:
    print('Percentage variance explained for EOF is {0:5.3}%, and for CBF is {1:5.3}%.'.format(eof_perc, cbf_perc))
else:
    print('Percentage variance explained for EOF is {0:5.3}%, and for CBF is {1}%.'.format(eof_perc, [int(i*100)/100 for i in cbf_perc]))


# Create plot
if iplot:
    if single_model:
        cbf_plot = cbf.copy()
    else:
        cbf_plot = cbf[1,:,:].copy()
    eof = np.asarray(eof.data)
    clevs = np.linspace(-6, 6, 21)
    extents = [-179.9, 180, 20, 90]
    lon = odat.coord('longitude').points
    lat = odat.coord('latitude').points
    proj = ccrs.PlateCarree(central_longitude=0.0)

    ax1 = plt.subplot(2, 1, 1, projection=proj)
    ax1.coastlines()
    ax1.set_global()
    ax1.contourf(lon, lat, eof[:,:]/100, levels=clevs, cmap=plt.cm.RdBu_r, transform=proj)
    ax1.axis(extents)
    ax1.set_title('ERA 20th Century')

    ax2 = plt.subplot(2, 1, 2, projection=proj)
    ax2.coastlines()
    ax2.set_global()
    cs = ax2.contourf(lon, lat, cbf_plot[:,:]/100, levels=clevs, cmap=plt.cm.RdBu_r, transform=proj)
    ax2.axis(extents)
    ax2.set_title('NorESM2-LM Historical')

    f1 = plt.gcf()
    f1.subplots_adjust(right=0.85)
    cbar_ax = f1.add_axes([0.88, 0.18, 0.02, 0.65])
    f1.colorbar(cs, cax=cbar_ax)

    plt.show()

