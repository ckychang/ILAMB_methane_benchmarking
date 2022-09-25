#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 09:24:54 2019

@author: ckychang
"""

"""
Converts csv files downloaded from the Ameriflux archive into a CF-compliant netCDF4 file.
"""
import glob
import numpy as np
from netCDF4 import Dataset
# declare parameters
mon = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
mon_doy = np.cumsum(mon)
# calculate the start and end DOY for each month
mon_start_doy = np.zeros(12)
mon_end_doy = np.zeros(12)
for m in range(0,12):
    if m==0:
        mon_start_doy[m] = 1
        mon_end_doy[m] = mon_doy[m]
    else:
        mon_start_doy[m] = mon_doy[m-1]+1
        mon_end_doy[m] = mon_doy[m]

# list all BU files
files_fch4_BU = sorted(glob.glob("../first_gridcell_outputs/regrid_180x360/regrid_*ilamb.nc"))
for file in range(0, len(files_fch4_BU)):
    fch4_sat_no_mask = np.zeros((120, 180, 360)) # masked array of zeros
    fch4_sat_w_mask = np.ma.masked_array(np.zeros((120, 180, 360)),mask=True) # masked array of zeros
    dataset = Dataset(files_fch4_BU[file])
    fch4 = np.array(dataset.variables['FCH4'])
    time = np.array(dataset.variables['time'])
    time_bounds = np.array(dataset.variables['time_bounds'])
    LATs = np.array(dataset.variables['lat'])
    LONs = np.array(dataset.variables['lon'])
    filename_no_mask = files_fch4_BU[file].replace('../first_gridcell_outputs/regrid_180x360','./no_mask').replace('regrid_','wetland_')
    filename_w_mask = files_fch4_BU[file].replace('../first_gridcell_outputs/regrid_180x360','./w_mask').replace('regrid_','wetland_')
    fch4[fch4<=0] = np.nan
    if (files_fch4_BU[file]=='../first_gridcell_outputs/regrid_180x360/regrid_TRIPLEX-GHG_CRUJRA_WAD2M_SWAMPS_GCP2019_MASK11_ilamb.nc'):
        fch4 = fch4/12.0*16.0
    print("processing %s" % (files_fch4_BU[file]))
    for t in range(0, 120):
        fch4_sat_no_mask[t,:,:] = np.nan_to_num(fch4[t,:,:], nan=0.0)
        fch4_sat_w_mask[t,:,:] = np.ma.masked_invalid(fch4[t,:,:])
 	# Open a dataset for writing
    dset = Dataset(filename_no_mask,mode="w")
 	# Create netCDF dimensions
    dset.createDimension("time",size=time.size)
    dset.createDimension("lon",size=LONs.size)
    dset.createDimension("lat",size=LATs.size)
    dset.createDimension("nb"  ,size=2)
 	# Create netCDF variables
    T  = dset.createVariable("time"       ,time.dtype   ,("time"       ))
    TB = dset.createVariable("time_bounds",time_bounds.dtype   ,("time","nb"  ))
    X  = dset.createVariable("lat"        ,LATs.dtype ,("lat"       ))
    Y  = dset.createVariable("lon"        ,LONs.dtype ,("lon"       ))
    FCH4_out  = dset.createVariable("FCH4",fch4_sat_no_mask.dtype,("time","lat","lon"))
 	# Load data and encode attributes
    T [...]    = time
    T.units    = "days since 1850-01-01"
    T.calendar = "noleap"
    T.bounds   = "time_bounds"
    TB[...]    = time_bounds
    X [...]    = LATs
    X.units    = "degrees_north"
    Y [...]    = LONs
    Y.units    = "degrees_east"
    FCH4_out[...]     = fch4_sat_no_mask
    FCH4_out.units    = "g C m-2 d-1"
    dset.close()
    # Open a dataset for writing
    dset = Dataset(filename_w_mask,mode="w")
 	# Create netCDF dimensions
    dset.createDimension("time",size=time.size)
    dset.createDimension("lon",size=LONs.size)
    dset.createDimension("lat",size=LATs.size)
    dset.createDimension("nb"  ,size=2)
 	# Create netCDF variables
    T  = dset.createVariable("time"       ,time.dtype   ,("time"       ))
    TB = dset.createVariable("time_bounds",time_bounds.dtype   ,("time","nb"  ))
    X  = dset.createVariable("lat"        ,LATs.dtype ,("lat"       ))
    Y  = dset.createVariable("lon"        ,LONs.dtype ,("lon"       ))
    FCH4_out  = dset.createVariable("FCH4",fch4_sat_w_mask.dtype,("time","lat","lon"))
 	# Load data and encode attributes
    T [...]    = time
    T.units    = "days since 1850-01-01"
    T.calendar = "noleap"
    T.bounds   = "time_bounds"
    TB[...]    = time_bounds
    X [...]    = LATs
    X.units    = "degrees_north"
    Y [...]    = LONs
    Y.units    = "degrees_east"
    FCH4_out[...]     = fch4_sat_w_mask
    FCH4_out.units    = "g C m-2 d-1"
    dset.close()
    print("finish encoding %s" % (files_fch4_BU[file]))
# list all TD files
files_fch4_TD = sorted(glob.glob("../first_gridcell_outputs/regrid_180x360/regrid_*wetl.nc"))
for file in range(0, len(files_fch4_TD)):
    fch4_sat_no_mask = np.zeros((120, 180, 360)) # masked array of zeros
    fch4_sat_w_mask = np.ma.masked_array(np.zeros((120, 180, 360)),mask=True) # masked array of zeros
    dataset = Dataset(files_fch4_TD[file])
    fch4 = np.array(dataset.variables['FCH4'])
    time = np.array(dataset.variables['time'])
    time_bounds = np.array(dataset.variables['time_bounds'])
    LATs = np.array(dataset.variables['lat'])
    LONs = np.array(dataset.variables['lon'])
    filename_no_mask = files_fch4_TD[file].replace('../first_gridcell_outputs/regrid_180x360','./no_mask').replace('regrid_','wetland_')
    filename_w_mask = files_fch4_TD[file].replace('../first_gridcell_outputs/regrid_180x360','./w_mask').replace('regrid_','wetland_')
    fch4[fch4<=0.0] = np.nan
    print("processing %s" % (files_fch4_TD[file]))
    for t in range(0, 120):
        fch4_sat_no_mask[t,:,:] = np.nan_to_num(fch4[t,:,:], nan=0.0)
        fch4_sat_w_mask[t,:,:] = np.ma.masked_invalid(fch4[t,:,:])
	# Open a dataset for writing
    dset = Dataset(filename_no_mask,mode="w")
	# Create netCDF dimensions
    dset.createDimension("time",size=time.size)
    dset.createDimension("lon",size=LONs.size)
    dset.createDimension("lat",size=LATs.size)
    dset.createDimension("nb"  ,size=2)
	# Create netCDF variables
    T  = dset.createVariable("time"       ,time.dtype   ,("time"       ))
    TB = dset.createVariable("time_bounds",time_bounds.dtype   ,("time","nb"  ))
    X  = dset.createVariable("lat"        ,LATs.dtype ,("lat"       ))
    Y  = dset.createVariable("lon"        ,LONs.dtype ,("lon"       ))
    FCH4_out  = dset.createVariable("FCH4",fch4_sat_no_mask.dtype,("time","lat","lon"))
	# Load data and encode attributes
    T [...]    = time
    T.units    = "days since 1850-01-01"
    T.calendar = "noleap"
    T.bounds   = "time_bounds"
    TB[...]    = time_bounds
    X [...]    = LATs
    X.units    = "degrees_north"
    Y [...]    = LONs
    Y.units    = "degrees_east"
    FCH4_out[...]     = fch4_sat_no_mask
    FCH4_out.units    = "g C m-2 d-1"
    dset.close()
    # Open a dataset for writing
    dset = Dataset(filename_w_mask,mode="w")
	# Create netCDF dimensions
    dset.createDimension("time",size=time.size)
    dset.createDimension("lon",size=LONs.size)
    dset.createDimension("lat",size=LATs.size)
    dset.createDimension("nb"  ,size=2)
	# Create netCDF variables
    T  = dset.createVariable("time"       ,time.dtype   ,("time"       ))
    TB = dset.createVariable("time_bounds",time_bounds.dtype   ,("time","nb"  ))
    X  = dset.createVariable("lat"        ,LATs.dtype ,("lat"       ))
    Y  = dset.createVariable("lon"        ,LONs.dtype ,("lon"       ))
    FCH4_out  = dset.createVariable("FCH4",fch4_sat_w_mask.dtype,("time","lat","lon"))
	# Load data and encode attributes
    T [...]    = time
    T.units    = "days since 1850-01-01"
    T.calendar = "noleap"
    T.bounds   = "time_bounds"
    TB[...]    = time_bounds
    X [...]    = LATs
    X.units    = "degrees_north"
    Y [...]    = LONs
    Y.units    = "degrees_east"
    FCH4_out[...]     = fch4_sat_w_mask
    FCH4_out.units    = "g C m-2 d-1"
    dset.close()
    print("finish encoding %s" % (files_fch4_TD[file]))
print("done")
