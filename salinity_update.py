# Annual mean SSS from 2D monthly means saved with the COSIMA output

# Arguments are archive directory and run id name.
# Output file is created in current working directory

import netCDF4, sys, os, glob, datetime
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

    # Dimensions and coordinate variables
    for dim in ['xt_ocean', 'yt_ocean']:
        din = d.variables[dim]
        dnew.createDimension(dim, len(din))
        dnew.createVariable(dim, np.float64, (dim,))
        var = dnew.variables[dim]
        for attr in din.ncattrs():
            setattr(var, attr, getattr(din,attr))
        var[:] = din[:]

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
    dnew.createVariable('surface_salt', np.float32, ('time', 'yt_ocean', 'xt_ocean'), fill_value=sss_in._FillValue)
    sss = dnew.variables['surface_salt']
    for attr in sss_in.ncattrs():
        if attr != "_FillValue":
            setattr(sss, attr, getattr(sss_in,attr))
    d.close()
    return dnew

fname = 'ocean_surface_salt_%s.nc' % runid
if Path(fname).exists():
    dout = netCDF4.Dataset(fname, 'r+')
else:
    dout = createfile(fname, lastfile)

time = dout.variables['time']
time_bnds = dout.variables['time_bnds']

sss = dout.variables['surface_salt']
nt = len(time)

if nt > 0:
    lastdate = netCDF4.num2date(time[-1], time.units, time.calendar)
    lastyear = lastdate.year
else:
    # Set lastyear to year of first file -1 so loop starts correctly
    firstfile = flist[0]
    lastyear = int(firstfile[-10:-6]) - 1

print("LASTYEAR", lastyear, lastfile_year)

if lastfile_year > lastyear:
    print('Data to process', runid, lastyear+1, lastfile_year)
    # Loop is over expected years, so missing files will cause an error.
    for year in range(lastyear+1, lastfile_year+1):
        d = xr.open_mfdataset(os.path.join(archivedir,'%s%4.4d_*.nc' % (prefix,year)),
                              combine='by_coords', use_cftime=True)

        offset = len(time)
        # Check whether the dates match
        if offset:
            lastdate = netCDF4.num2date(time_bnds[-1,1], time.units, 'proleptic_gregorian')
            newdate = d.time.values[0]
            # Should be ~ 16 days between end of year and middle of Jan
            if not 15 <= (newdate-lastdate).days <= 17:
                    print("Date mismatch", lastdate, newdate, newdate-lastdate, file=sys.stderr)
                    raise Exception('Date mismatch')

        for vname in dout.variables:
            v = dout.variables[vname]
            if not (vname.startswith('time') or vname in ('xt_ocean', 'yt_ocean')):
                v[offset] = d[vname].mean(dim='time').values
        date0 = datetime.datetime(year, 1, 1)
        date1 = datetime.datetime(year+1, 1, 1)
        time_bnds[offset,0] = netCDF4.date2num(date0, units=time.units, calendar=time.calendar)
        time_bnds[offset,1] = netCDF4.date2num(date1, units=time.units, calendar=time.calendar)
        time[offset] = time_bnds[offset].mean()

        d.close()
        dout.sync()  # Sync to disk once per year

dout.close()
