"""
Converts csv files downloaded from the Ameriflux archive into a CF-compliant netCDF4 file.
"""
import glob
import pandas as pd
import numpy as np
from netCDF4 import Dataset
# from cf_units import Unit

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

# list of sites
csvs = glob.glob("*.csv")
# read in site info
data = pd.read_excel("Powell_site_LATLON_V3.xlsx")
IDs = data['Site Code']
sites = ' '.join(IDs).replace('-','').split()
LATs = np.array(data['LAT'])
LONs = np.array(data['LON'])
# # build time array that covers 2006-2018
nyears    = 2018-2006+1
month_bnd = np.asarray([0,31,59,90,120,151,181,212,243,273,304,334,365],dtype=float)
tbnd  = np.asarray([((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[:-1]).flatten(),
                    ((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[+1:]).flatten()]).T
tbnd += (2006-1850)*365
t     = tbnd.mean(axis=1)
t_label_agg  = np.asarray([((np.arange(nyears)*365)[:,np.newaxis]+month_bnd[:-1]).flatten()]).T
t_label_agg += (2006-1850)*365

# FCH4_agg = np.ma.masked_array(np.zeros((len(time), len(sites))),mask=True) # masked array of zeros
FCH4_F_ANN_agg = np.ma.masked_array(np.zeros(((t_label_agg.size), \
                                              len(sites))),mask=True) # masked array of zeros

lat_agg = np.zeros(len(sites))
lon_agg = np.zeros(len(sites))
count = 0

# loop through sites
wetland = 0
for site in sites:
    print("processing %s" % (site))
    filename = ' '.join([site.upper(),'_V2.csv']).replace(' ','')
    if filename in csvs: 
        print("reading %s" % (filename))
        site_data = pd.read_csv(filename)
        FCH4_F_ANN = np.ma.masked_values(np.array(site_data['FCH4_F_ANN_mean']), \
                                     -9999)
            # convert nmol/m2/s to gC/m2/d
        FCH4_F_ANN = np.ma.masked_where(FCH4_F_ANN == -9999, FCH4_F_ANN)*12*(10**-9)*86400
        # calculate obs year and DOY
        years = np.array(site_data['Year'])
        years_label = np.unique(np.array(site_data['Year']))
        DOYs = np.array(site_data['DOY'])
        t_label = np.array([])
        FCH4_F_ANN_label = np.array([])
        for yr in years_label:
            # monthly means
            for m in range(0, 12):
                # locate DOY in each month
                idx = np.array(np.where((DOYs<=mon_end_doy[m])&(DOYs>=mon_start_doy[m])&(years == yr)))
                idx = idx.reshape([idx.size,1])
                # if the monthly data is available
                if idx.size==mon[m]:
                    FCH4_F_ANN_label = np.append(FCH4_F_ANN_label, np.nanmean(FCH4_F_ANN[idx]))
                    # days since 1850
                    t_label = np.append(t_label, DOYs[idx[0]]+(yr-1850)*365-1)
            # end monthly loop in each site
        idx_start = np.array(np.where(t_label_agg == t_label[0]))
        idx_end = np.array(np.where(t_label_agg == t_label[-1]))
        if (idx_start.size>0 and idx_end.size>0):
            idx = np.arange(idx_start[0], idx_end[0]+1, dtype=np.int)
            FCH4_F_ANN_agg[idx,count] = FCH4_F_ANN_label
            wetland +=1
    count +=1
    print("finish encoding %s" % (site))

# Open a dataset for writing
dset = Dataset('FCH4_F_ANN_monthly_wetland_tier1.nc',mode="w")
# Create netCDF dimensions
dset.createDimension("time",size=t.size)
dset.createDimension("data",size=LATs.size)
dset.createDimension("nb"  ,size=2)
# Create netCDF variables
T  = dset.createVariable("time"       ,t.dtype   ,("time"       ))
TB = dset.createVariable("time_bounds",tbnd.dtype   ,("time","nb"  ))
X  = dset.createVariable("lat"        ,LATs.dtype ,("data"       ))
Y  = dset.createVariable("lon"        ,LONs.dtype ,("data"       ))
FCH4_F_ANN_out  = dset.createVariable("FCH4",FCH4_F_ANN_agg.dtype,("time","data"))
# Load data and encode attributes
T [...]    = t
T.units    = "days since 1850-01-01"
T.calendar = "noleap"
T.bounds   = "time_bounds"
TB[...]    = tbnd
X [...]    = LATs
X.units    = "degrees_north"
Y [...]    = LONs
Y.units    = "degrees_east"
FCH4_F_ANN_out[...]     = FCH4_F_ANN_agg
FCH4_F_ANN_out.units    = "g C m-2 d-1"
dset.close()
print("finish encoding FCH4_F_ANN_monthly_wetland_tier1.nc")

