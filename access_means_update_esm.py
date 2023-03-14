# Update files of monthly and annual means
# This version works with ESM files of form SSP-126-ext-05.pa-2108011001.nc

import netCDF4, sys, os, glob, datetime
from pathlib import Path
import numpy as np

allow_missing = True

if len(sys.argv) > 1:
    runid = sys.argv[1]
else:
    # Get it from the directory name
    runid = os.path.basename(os.getcwd())

# Check for any new complete years
flist = sorted(list(glob.glob('%s.pa*012001.nc' % runid)))
if not flist:
    print("Nothing to process")
    sys.exit(0)

lastfile = flist[-1]
lastfile_year = int(lastfile[-13:-9])

def createfile(filename,lastfile):
    d = netCDF4.Dataset(filename, 'w')
    # Any file will do for the variable definition
    dnew = netCDF4.Dataset(lastfile)
    for dname, dim in dnew.dimensions.items():
        if dim.isunlimited():
            d.createDimension(dname,None)
        else:
            d.createDimension(dname,dim.size)
    for vname, var in dnew.variables.items():
        if hasattr(var,'_FillValue'):
            newv = d.createVariable(vname, var.datatype, var.dimensions,
                                    zlib=True, fill_value=getattr(var,'_FillValue'))
        else:
            newv = d.createVariable(vname, var.datatype, var.dimensions,
                                    zlib=True)
        for attr in var.ncattrs():
            if attr != '_FillValue':
                newv.setncattr(attr, getattr(var,attr))
        # If this variable doesn't have a time dimension, copy
        # from input
        if 'time' not in var.dimensions:
            d.variables[vname][:] = dnew.variables[vname][:]
    return d

output_file = "%s.nc" % runid
if Path(output_file).exists():
    d = netCDF4.Dataset(output_file, 'r+')
else:
    d = createfile(output_file, lastfile)

time = d.variables['time']
nt = len(time)

if nt%12 != 0:
    raise Exception("Unexpected state: nt=%d is not a multiple of 12 in run %s" % (nt, runid))

if nt > 0:
    lastdate = netCDF4.num2date(time[-1], time.units,time.calendar)
    lastyear = lastdate.year
else:
    # Set lastyear to year of first file -1 so loop starts correctly
    firstfile = flist[0]
    lastyear = int(firstfile[-13:-9]) - 1

if lastfile_year > lastyear:
    print('Data to process', lastyear+1, lastfile_year)
    # Loop is over expected years, so missing files will cause an
    # error.
    for year in range(lastyear+1, lastfile_year+1):
        # Processing glitches sometimes mean there's a month missing.
        # Check this first to avoid partial updates
        flist_year = glob.glob('%s.pa-%4.4d???001.nc' % (runid, year))
        if len(flist_year) != 12:
            print("Missing files for year %d:" % year, flist_year, file=sys.stderr)
            raise Exception("Missing files in %s" % runid)
        for mon in range(12):
            fname = '%s.pa-%4.4d%3.3d001.nc' % (runid, year, mon+1)
            print(fname)
            dnew = netCDF4.Dataset(fname)
            # Allow for some variables to be missing
            for vname in d.variables:
                try:
                    v = d.variables[vname]
                    if 'time' in v.dimensions:
                        # Get the same variable from the new file
                        vnew = dnew.variables[vname]
                        # Append the data
                        v[nt] = vnew[:]
                except KeyError:
                    if allow_missing:
                        print("Warning missing %s in %s" % (vname, fname))
                    else:
                        raise Exception("Missing %s in %s" % (vname, fname))
            nt += 1
            dnew.close()
d.close()

# Update the annual means
ann_file = "%s_ann.nc" % runid
if Path(ann_file).exists():
    d = netCDF4.Dataset(ann_file, 'r+')
else:
    d = createfile(ann_file, lastfile)

time = d.variables['time']
nt = len(time)

if nt > 0:
    lastdate = netCDF4.num2date(time[-1], time.units,time.calendar)
    lastyear = lastdate.year
else:
    # Set lastyear to year of first file -1 so loop starts correctly
    print("FLIST", flist)
    firstfile = flist[0]
    lastyear = int(firstfile[-13:-9]) - 1

def isleap(year):
    # Proleptic gregreogian
    if year % 100 == 0:
        return year%400==0
    else:
        return year%4==0

if lastfile_year > lastyear:
    print('Data to process for annual means', lastyear+1, lastfile_year)
    for year in range(lastyear+1, lastfile_year+1):
        if isleap(year):
            mwts = np.array([31,29,31,30,31,30,31,31,30,31,30,31])/366.
        else:
            mwts = np.array([31,28,31,30,31,30,31,31,30,31,30,31])/365.
        for mon in range(12):
            fname = '%s.pa-%4.4d%3.3d001.nc' % (runid, year, mon+1)
            print(fname)
            dnew = netCDF4.Dataset(fname)
            for vname in d.variables:
                v = d.variables[vname]
                if 'time' in v.dimensions:
                    # Get the same variable from the new file
                    try:
                        vnew = dnew.variables[vname]
                        # Append the data with monthly weights to
                        # calculate the annual means
                        if vname == 'time_bnds':
                            if mon==0:
                                v[nt] = vnew[:]
                            elif mon==11:
                                # End bounds
                                v[nt,1] = vnew[:,1]
                        else:
                            if mon==0:
                                v[nt] = mwts[mon]*vnew[:]
                            else:
                                v[nt] += mwts[mon]*vnew[0]
                    except KeyError:
                        if allow_missing:
                            print("Warning missing %s in %s" % (vname, fname))
                        else:
                            raise Exception("Missing %s in %s" % (vname, fname))
            dnew.close()
        nt += 1

d.close()
