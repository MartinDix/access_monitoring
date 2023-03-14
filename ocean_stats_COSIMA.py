# Global mean statistics KE from ocean-scalar files saved with the COSIMA output
# File names of form ocean-scalar-1-daily-ym_0001_01.nc
# Only create annual variables

# Arguments are archive directory and run id name.
# Output file is created in current working directory

import netCDF4, sys, os, glob, datetime
from pathlib import Path
import numpy as np, xarray as xr

archivedir = sys.argv[1]
runid = sys.argv[2]

prefix = 'ocean-scalar-1-daily-ym_'
# Check for any new complete years
# Assume running in 6 month chunks here.
flist = sorted(list(glob.glob(os.path.join(archivedir,prefix+'[0-9]*_07.nc'))))
if not flist:
    print("Nothing to process")
    sys.exit(0)

lastfile = flist[-1]
lastfile_year = int(lastfile[-10:-6])

# Map from variable names used here to those in the MOM file
name_map = {'ke_tot':'ke_tot', 'sst':'temp_surface_ave', 'sea_level':'eta_global',
           'sss':'salt_surface_ave', 'ohc':'total_ocean_heat', 'temp':'temp_global_ave',
           'year':'time'}

def createfile(filename,lastfile):
    # Create a new file for the global mean with time dimension
    dnew = netCDF4.Dataset(filename, 'w')
    # Any file will do for getting dimensions etc
    d = netCDF4.Dataset(lastfile)

    # Dimensions and coordinate variables
    dnew.createDimension('nv', 2)
    time = d.variables['time']
    # time_bounds = d.variables['time_bounds']
    # for newvar in ('time', 'year'):
    for newvar in ('year',):
        dnew.createDimension(newvar, None)
        dnew.createVariable(newvar, np.float64, (newvar,))
        tnew = dnew.variables[newvar]
        for attr in time.ncattrs():
            # Missing value in the bounds shouldn't occur
            # MOM files incorrectly set calendar as gregorian rather
            # than proleptic
            if attr not in ("missing_value", "_FillValue", "calendar", "calendar_type"):
                setattr(tnew, attr, getattr(time,attr))
        tnew.calendar = "proleptic_gregorian"
        # This case isn't handled correctly by the previous iteration
        if newvar == 'year':
            tnew.bounds = 'year_bounds'
        dnew.createVariable('%s_bounds' % newvar, np.float64, (newvar,'nv'))
        tnew = dnew.variables['%s_bounds' % newvar]
        # for attr in time_bounds.ncattrs():
        #     # Missing value in the bounds shouldn't occur
        #     if attr not in ("missing_value", "_FillValue"):
        #         setattr(tnew, attr, getattr(time_bounds,attr))

    # Variables
    dnew.createVariable('ke_tot_ann', np.float32, ('year',))
    ke = dnew.variables['ke_tot_ann']
    ke.units = '1e15 J'
    ke.long_name = "Globally integrated ocean kinetic energy"
    dnew.createVariable('sst_ann', np.float32, ('year'))
    sst = dnew.variables['sst_ann']
    sst.units = 'K'
    sst.standard_name = "sea_surface_temperature"
    dnew.createVariable('sea_level_ann', np.float32, ('year'))
    sea_level = dnew.variables['sea_level_ann']
    sea_level.units = 'm'
    sea_level.standard_name = "sea_surface_height_above_geoid"
    dnew.createVariable('sss_ann', np.float32, ('year'))
    sss = dnew.variables['sss_ann']
    sss.units = 'psu'
    sss.standard_name = "sea_surface_salinity"
    dnew.createVariable('ohc_ann', np.float32, ('year'))
    ohc = dnew.variables['ohc_ann']
    ohc.units = '1e25 J'
    ohc.long_name = "ocean_heat_content"
    dnew.createVariable('temp_ann', np.float32, ('year'))
    temp = dnew.variables['temp_ann']
    temp.units = 'deg_C'
    temp.long_name =  "Global mean temp in liquid seawater"
    temp:standard_name = "sea_water_potential_temperature"
    # dnew.createVariable('temp_zmean', np.float32, ('year'))
    # tmean = dnew.variables['temp_zmean']
    # tmean.units = 'K'
    # tmean.standard_name = "sea_water_potential_temperature"
    # dnew.createVariable('salt_zmean', np.float32, ('year'))
    # smean = dnew.variables['salt_zmean']
    # smean.units = 'psu'
    # smean.standard_name = "sea_water_salinity"

    # for vname in ('ke_tot','sst', 'sea_level', 'sss', 'ohc', 'temp_zmean', 'salt_zmean'):
    #     # Create annual mean
    #     var = dnew.variables[vname]
    #     annname = '%s_ann' % vname
    #     dnew.createVariable(annname, np.float32, ('year',))
    #     annvar = dnew.variables[annname]
    #     for attr in var.ncattrs():
    #         setattr(annvar, attr, getattr(var,attr))
    d.close()
    return dnew

fname = 'ocean_stats_%s.nc' % runid
if Path(fname).exists():
    dout = netCDF4.Dataset(fname, 'r+')
else:
    dout = createfile(fname, lastfile)

# time = dout.variables['time']
# nt = len(time)

# if nt%12 != 0:
#     raise Exception("Unexpected state: In file %s. nt=%d is not a multiple of 12" % (fname, nt))

# time_bounds = dout.variables['time_bounds']
yearvar = dout.variables['year']
year_bounds = dout.variables['year_bounds']
ke = dout.variables['ke_tot_ann']
sst = dout.variables['sst_ann']
sea_level = dout.variables['sea_level_ann']
sss = dout.variables['sss_ann']
ohc = dout.variables['ohc_ann']
# tmean = dout.variables['temp_zmean']
# smean = dout.variables['salt_zmean']
nt = len(yearvar)

if nt > 0:
    lastdate = netCDF4.num2date(yearvar[-1], yearvar.units, yearvar.calendar)
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
        # Note that raw MOM data wrongly has Gregorian calendar but this is
        # fixed in post-processing
        d = xr.open_mfdataset(os.path.join(archivedir,'%s%4.4d_*.nc' % (prefix,year)),
                              combine='by_coords', use_cftime=True)
        offset = len(yearvar)
        # Check whether the dates match
        if offset:
            lastdate = netCDF4.num2date(yearvar[-1], yearvar.units, yearvar.calendar)
            newdate = d.time.values[0]
            if not newdate.year == lastdate.year+1:
                print("Date mismatch", lastdate, newdate, file=sys.stderr)
                raise Exception('Date mismatch')

        for vname in dout.variables:
            v = dout.variables[vname]
            if not vname.startswith('year'):
                # Strip the _ann
                v[offset] = d[name_map[vname[:-4]]].mean(dim='time').values
        date0 = datetime.datetime(year, 1, 1)
        date1 = datetime.datetime(year+1, 1, 1)
        year_bounds[offset,0] = netCDF4.date2num(date0, units=yearvar.units, calendar=yearvar.calendar)
        year_bounds[offset,1] = netCDF4.date2num(date1, units=yearvar.units, calendar=yearvar.calendar)
        yearvar[offset] = year_bounds[offset].mean()

        d.close()
        dout.sync()  # Sync to disk once per year

dout.close()
