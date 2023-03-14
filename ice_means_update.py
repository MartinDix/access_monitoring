# Calculate hemispheric means of the ice area and volume from CICE
# monthly history files

# Arguments are archive directory and run id name

import netCDF4, sys, os, glob
from pathlib import Path
import numpy as np

archivedir = sys.argv[1]
runid = sys.argv[2]

# Check for any new complete years
flist = sorted(list(glob.glob(os.path.join(archivedir,'iceh_m.[0-9]*-12.nc'))))
if not flist:
    print("Nothing to process")
    sys.exit(0)

lastfile = flist[-1]
lastfile_year = int(lastfile[-10:-6])

def createfile(filename,lastfile):
    # Create a new file for the means
    dnew = netCDF4.Dataset(filename, 'w')
    # Any file will do for getting dimensions etc
    d = netCDF4.Dataset(lastfile)

    # Dimensions and coordinate variables
    time = d.variables['time']
    time_bounds = d.variables['time_bounds']
    dnew.createDimension('time', None)
    dnew.createVariable('time', np.float64, ('time',))
    tnew = dnew.variables['time']
    for attr in time.ncattrs():
        # Missing value in the bounds shouldn't occur
        # calendar icorrectly set as gregorian
        if attr not in ("missing_value", "_FillValue", "calendar"):
            setattr(tnew, attr, getattr(time,attr))
    tnew.calendar = "proleptic_gregorian"
    dnew.createDimension('nv', 2)
    dnew.createVariable('time_bounds', np.float64, ('time','nv'))
    tnew = dnew.variables['time_bounds']
    for attr in time_bounds.ncattrs():
        if attr not in ("missing_value", "_FillValue", "calendar"):
            setattr(tnew, attr, getattr(time_bounds,attr))
    tnew.calendar = "proleptic_gregorian"
    # Variables
    for vname in ('ice_area_nh', 'ice_area_sh'):
        dnew.createVariable(vname, np.float32, ('time',))
        var = dnew.variables[vname]
        var.units = 'm2'
        var.standard_name = "sea_ice_area"
    for vname in ('ice_vol_nh', 'ice_vol_sh'):
        dnew.createVariable(vname, np.float32, ('time',))
        var = dnew.variables[vname]
        var.units = 'm3'
        var.standard_name = "sea_ice_volume"
    for vname in ('ice_ext_nh', 'ice_ext_sh'):
        dnew.createVariable(vname, np.float32, ('time',))
        var = dnew.variables[vname]
        var.units = 'm3'
        var.standard_name = "sea_ice_extent"

    d.close()
    return dnew

fname = 'ice_mean_%s.nc' % runid
if Path(fname).exists():
    dout = netCDF4.Dataset(fname, 'r+')
else:
    dout = createfile(fname, lastfile)

time = dout.variables['time']
nt = len(time)

if nt%12 != 0:
    raise Exception("Unexpected state: In file %s. nt=%d is not a multiple of 12" % (fname, nt))

time_bounds = dout.variables['time_bounds']
ice_area_nh = dout.variables['ice_area_nh']
ice_area_sh = dout.variables['ice_area_sh']
ice_vol_nh = dout.variables['ice_vol_nh']
ice_vol_sh = dout.variables['ice_vol_sh']
ice_ext_nh = dout.variables['ice_ext_nh']
ice_ext_sh = dout.variables['ice_ext_sh']

if nt > 0:
    lastdate = netCDF4.num2date(time[-1], time.units,time.calendar)
    lastyear = lastdate.year
else:
    # Set lastyear to year of first file -1 so loop starts correctly
    firstfile = flist[0]
    lastyear = int(firstfile[-10:-6]) - 1

if lastfile_year > lastyear:
    print('Ice data to process', runid, lastyear+1, lastfile_year)
    # Loop is over expected years, so missing files will cause an
    # error.
    for year in range(lastyear+1, lastfile_year+1):
        # Processing glitches sometimes mean there's a month missing.
        # Check this first to avoid partial updates
        flist_year = glob.glob(os.path.join(archivedir,'iceh_m.%4.4d-[0-9]*.nc'%year))
        if len(flist_year) != 12:
            print("Missing files for year %d:" % year, flist_year, file=sys.stderr)
            raise Exception("Missing files")
        for mon in range(1,13):
            fname = os.path.join(archivedir,'iceh_m.%4.4d-%2.2d.nc' % (year, mon))
            print(fname)
            d = netCDF4.Dataset(fname)
            time_in = d.variables['time']
            time_bounds_in = d.variables['time_bounds']
            aice = d.variables['aice']
            hi = d.variables['hi'] # Grid box mean thickness
            area = d.variables['tarea']
            nlat = area.shape[0]

            offset = len(time)
            # Ice file dates are end of month, so instead check the bounds
            # Check whether the dates match
            if offset:
                lastdate = netCDF4.num2date(time_bounds[-1][1], time.units, time.calendar)
                newdate = netCDF4.num2date(time_bounds_in[0][0], time_in.units, 'proleptic_gregorian')
                if newdate != lastdate:
                    print("Date mismatch", lastdate, newdate, newdate-lastdate, file=sys.stderr)
                    raise Exception('Date mismatch')

            for t in range(len(time_in)):
                ice_area_sh[offset+t] = np.ma.sum(aice[t,:nlat//2,:] * area[:nlat//2,:])
                ice_area_nh[offset+t] = np.ma.sum(aice[t,nlat//2:,:] * area[nlat//2:,:])
                ice_vol_sh[offset+t] = np.ma.sum(hi[t,:nlat//2,:] * area[:nlat//2,:])
                ice_vol_nh[offset+t] = np.ma.sum(hi[t,nlat//2:,:] * area[nlat//2:,:])
                extent = np.ma.where(aice[t] > 0.15, 1., 0.)
                ice_ext_sh[offset+t] = np.ma.sum(extent[:nlat//2,:] * area[:nlat//2,:])
                ice_ext_nh[offset+t] = np.ma.sum(extent[nlat//2:,:] * area[nlat//2:,:])
                # Handle possible calendar changes
                date = netCDF4.num2date(time_bounds_in[t], time_in.units, 'proleptic_gregorian')
                time_bounds[offset+t] = netCDF4.date2num(date, time.units, time.calendar)
                # Input file time is end of month, so time here
                # as the middle of the bounds
                date = netCDF4.num2date(time_bounds_in[t].mean(), time_in.units, 'proleptic_gregorian')
                time[offset+t] = netCDF4.date2num(date, time.units, time.calendar)
            d.close()
        dout.sync()  # Sync to disk once per year

dout.close()
