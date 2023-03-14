# Global mean KE from ocean_scalar files

# Arguments are archive directory and run id name

import netCDF4, sys, os, glob
from pathlib import Path
import numpy as np

archivedir = sys.argv[1]
runid = sys.argv[2]

# Check for any new complete years
flist = sorted(list(glob.glob(os.path.join(archivedir,'ocean_scalar.nc-[0-9]*1231'))))
if not flist:
    print("Nothing to process")
    sys.exit(0)

lastfile = flist[-1]
lastfile_year = int(lastfile[-8:-4])

def isleap(year):
    # Proleptic gregreogian
    if year % 100 == 0:
        return year%400==0
    else:
        return year%4==0

def createfile(filename,lastfile):
    # Create a new file for the global mean with time dimension
    dnew = netCDF4.Dataset(filename, 'w')
    # Any file will do for getting dimensions etc
    d = netCDF4.Dataset(lastfile)

    # Dimensions and coordinate variables
    dnew.createDimension('nv', 2)
    time = d.variables['time']
    # This is inconsistent across model runs for some reason
    tbounds_name = getattr(time, "bounds")
    assert tbounds_name[:5] == 'time_'
    bounds_name = tbounds_name[5:]
    time_bounds = d.variables[tbounds_name]
    for newvar in ('time', 'year'):
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
            tnew.bounds = f'year_{bounds_name}'
        dnew.createVariable(f'{newvar}_{bounds_name}', np.float64, (newvar,'nv'))
        tnew = dnew.variables[f'{newvar}_{bounds_name}']
        for attr in time_bounds.ncattrs():
            # Missing value in the bounds shouldn't occur
            if attr not in ("missing_value", "_FillValue"):
                setattr(tnew, attr, getattr(time_bounds,attr))

    # Variables
    dnew.createVariable('ke_tot', np.float32, ('time',))
    ke = dnew.variables['ke_tot']
    ke.units = '1e15 J'
    ke.long_name = "Globally integrated ocean kinetic energy"
    for vname in ('ke_tot',):
        # Create annual mean
        var = dnew.variables[vname]
        annname = '%s_ann' % vname
        dnew.createVariable(annname, np.float32, ('year',))
        annvar = dnew.variables[annname]
        for attr in var.ncattrs():
            setattr(annvar, attr, getattr(var,attr))
    d.close()
    return dnew

fname = 'ocean_ke_%s.nc' % runid
if Path(fname).exists():
    dout = netCDF4.Dataset(fname, 'r+')
else:
    dout = createfile(fname, lastfile)

time = dout.variables['time']
nt = len(time)

if nt%12 != 0:
    raise Exception("Unexpected state: In file %s. nt=%d is not a multiple of 12" % (fname, nt))

try:
    time_bounds = dout.variables['time_bounds']
    bounds_name = 'bounds'
except KeyError:
    time_bounds = dout.variables['time_bnds']
    bounds_name = 'bnds'
ke = dout.variables['ke_tot']
ke_ann = dout.variables['ke_tot_ann']
yearvar = dout.variables['year']
year_bounds = dout.variables[f'year_{bounds_name}']

if nt > 0:
    lastdate = netCDF4.num2date(time[-1], time.units, time.calendar)
    lastyear = lastdate.year
else:
    # Set lastyear to year of first file -1 so loop starts correctly
    firstfile = flist[0]
    lastyear = int(firstfile[-8:-4]) - 1

if lastfile_year > lastyear:
    print('Data to process', runid, lastyear+1, lastfile_year)
    # Loop is over expected years, so missing files will cause an
    # error.
    for year in range(lastyear+1, lastfile_year+1):
        flist_year = glob.glob(os.path.join(archivedir,'ocean_scalar.nc-%4.4d[0-9]*' % (year)))
        if isleap(year):
            mwts = np.array([31,29,31,30,31,30,31,31,30,31,30,31])/366.
        else:
            mwts = np.array([31,28,31,30,31,30,31,31,30,31,30,31])/365.
        # Initial check that we have 12 months. May get a failure if
        # mppcombine is still running.
        nm = 0
        for f in sorted(flist_year):
            d = netCDF4.Dataset(f)
            time_in = d.variables['time']
            nm += len(time_in)
        if nm != 12:
            print("Missing files for year %d, nm=%d:" % (year, nm), flist_year, file=sys.stderr)
            raise Exception("Missing files")
        for f in sorted(flist_year):
            print(f)
            d = netCDF4.Dataset(f)
            time_in = d.variables['time']
            time_bounds_in = d.variables[f'time_{bounds_name}']
            ke_in = d.variables['ke_tot']

            offset = len(time)
            # Check whether the dates match
            if offset:
                lastdate = netCDF4.num2date(time[-1], time.units, time.calendar)
                newdate = netCDF4.num2date(time_in[0], time_in.units, 'proleptic_gregorian')
                if not 25 <= (newdate-lastdate).days <= 35:
                    print("Date mismatch", lastdate, newdate, newdate-lastdate, file=sys.stderr)
                    raise Exception('Date mismatch')

            # Global area fractions as a function of level
            for t in range(len(time_in)):
                ke[offset+t] = ke_in[t]
                # Handle possible changes in the base date
                # Ocean model files incorrectly have calendar attribute
                # set as gregorian, but really use proleptic_gregorian
                date = netCDF4.num2date(time_in[t], time_in.units, 'proleptic_gregorian')
                print("DATE", date)
                time[offset+t] = netCDF4.date2num(date, time.units, time.calendar)
                mon = date.month - 1  # Convert to an index
                date = netCDF4.num2date(time_bounds_in[t], time_in.units, 'proleptic_gregorian')
                time_bounds[offset+t] = netCDF4.date2num(date, time.units, time.calendar)
                # Update the annual means
                annt = (offset+t)//12
                for vname in dout.variables:
                    v = dout.variables[vname]
                    if 'time' in v.dimensions:
                        if vname.startswith('time'):
                            annvar = dout.variables[vname.replace('time','year')]
                        else:
                            annvar = dout.variables['%s_ann' % vname]
                        if vname == f'time_{bounds_name}':
                            if mon==0:
                                annvar[annt] = v[offset+t]
                            elif mon==11:
                                # End bounds
                                annvar[annt,1] = v[offset+t,1]
                        else:
                            if mon==0:
                                annvar[annt] = mwts[mon]*v[offset+t]
                            else:
                                annvar[annt] += mwts[mon]*v[offset+t]
            d.close()

        dout.sync()  # Sync to disk once per year

dout.close()
