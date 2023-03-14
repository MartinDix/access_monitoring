# Time series of monthly min and max SSS from 2D monthly means
# saved with the COSIMA output.
# There are zero values over land, so need to use a mask field.

# Arguments are archive directory and run id name.
# Output file is created in current working directory

import netCDF4, sys, os, glob, datetime, cf_units
from pathlib import Path
import numpy as np, xarray as xr

archivedir = sys.argv[1]
runid = sys.argv[2]

prefix = 'ocean-2d-surface_salt-1-monthly-mean-ym_'
# Check for any new complete years
# Assume running in 6 month chunks here.
flist = sorted(list(glob.glob(os.path.join(archivedir,prefix+'[0-9]*_07.nc'))))
if not flist:
    print("Nothing to process")
    sys.exit(0)

lastfile = flist[-1]
lastfile_year = int(lastfile[-10:-6])

def createfile(filename,lastfile):
    # Create a new file for the annual mean with time dimension
    dnew = netCDF4.Dataset(filename, 'w')
    # Any file will do for getting dimensions etc
    d = netCDF4.Dataset(lastfile)

    dnew.createDimension('nv', 2)
    time = d.variables['time']
    time_bnds = d.variables['time_bnds']
    dnew.createDimension('time', None)
    dnew.createVariable('time', np.float64, ('time',))
    tnew = dnew.variables['time']
    for attr in time.ncattrs():
        # Missing value in the bounds shouldn't occur
        # Raw MOM files incorrectly set calendar as gregorian rather
        # than proleptic, but fixed in post-processing
        if attr not in ("missing_value", "_FillValue", "calendar", "calendar_type"):
            setattr(tnew, attr, getattr(time,attr))
    tnew.calendar = "proleptic_gregorian"
    dnew.createVariable('time_bnds', np.float64, ('time', 'nv'))
    tnew = dnew.variables['time_bnds']
    for attr in time_bnds.ncattrs():
        # Missing value in the bounds shouldn't occur
        if attr not in ("missing_value", "_FillValue"):
            setattr(tnew, attr, getattr(time_bnds,attr))

    # Variables
    sss_in = d.variables['surface_salt']
    dnew.createVariable('surface_salt_min', np.float32, ('time', ), fill_value=sss_in._FillValue)
    sss_min = dnew.variables['surface_salt_min']
    for attr in sss_in.ncattrs():
        if attr != "_FillValue":
            setattr(sss_min, attr, getattr(sss_in,attr))
    sss_min.cell_methods = "area: min"
    dnew.createVariable('surface_salt_max', np.float32, ('time', ), fill_value=sss_in._FillValue)
    sss_max = dnew.variables['surface_salt_max']
    for attr in sss_in.ncattrs():
        if attr != "_FillValue":
            setattr(sss_max, attr, getattr(sss_in,attr))
    sss_max.cell_methods = "area: max"
    d.close()
    return dnew

fname = 'ocean_surface_salt_minmax_%s.nc' % runid
if Path(fname).exists():
    dout = netCDF4.Dataset(fname, 'r+')
else:
    dout = createfile(fname, lastfile)

time = dout.variables['time']
time_bnds = dout.variables['time_bnds']
taxis = cf_units.Unit(time.units, time.calendar)

sss_min = dout.variables['surface_salt_min']
sss_max = dout.variables['surface_salt_max']
nt = len(time)

if nt > 0:
    lastdate = netCDF4.num2date(time[-1], time.units, time.calendar)
    lastyear = lastdate.year
else:
    # Set lastyear to year of first file -1 so loop starts correctly
    firstfile = flist[0]
    lastyear = int(firstfile[-10:-6]) - 1

print("LASTYEAR", lastyear, lastfile_year)

# Ocean grid fraction
mask = xr.open_dataset('/g/data/ik11/inputs/access-om2/input_20201102/mom_025deg/ocean_mask.nc').mask

if lastfile_year > lastyear:
    print('Data to process', runid, lastyear+1, lastfile_year)
    # Loop is over expected years, so missing files will cause an error.
    # Note that raw MOM data wrongly has Gregorian calendar but this is
    # fixed in post-processing. Still do an explicit decode to get
    # the correct units for the bounds
    for year in range(lastyear+1, lastfile_year+1):
        d = xr.open_mfdataset(os.path.join(archivedir,'%s%4.4d_*.nc' % (prefix,year)),
                              combine='by_coords', decode_times=False)
        d.time.attrs["calendar"] = "proleptic_gregorian"
        # MOM files set this w/o a base time which confuses xarray
        d.time_bnds.attrs['units'] = d.time.units
        d = xr.decode_cf(d, use_cftime=True)

        offset = len(time)
        # Check whether the dates match
        if offset:
            lastdate = netCDF4.num2date(time_bnds[-1,1], time.units, time.calendar)
            newdate = d.time.values[0]
            # Should be ~ 16 days between end of year and middle of Jan
            if not 15 <= (newdate-lastdate).days <= 17:
                    print("Date mismatch", lastdate, newdate, newdate-lastdate, file=sys.stderr)
                    raise Exception('Date mismatch')
        for t in range(12):
            sss_min[offset+t] = d['surface_salt'][t].values[mask > 0].min()
            sss_max[offset+t] = d['surface_salt'][t].values[mask > 0].max()

            time_bnds[offset+t] =  taxis.date2num(d.time_bnds[t])
            time[offset+t] = time_bnds[offset+t].mean()

        d.close()
        dout.sync()  # Sync to disk once per year

dout.close()
