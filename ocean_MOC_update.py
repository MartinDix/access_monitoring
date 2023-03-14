# Calculate diagnostics of the ocean model meridional overturning
# Original ferret calculation doesn't handle staggered grids properly, so just use
# simple numpy arrays here

# Basin file is

# Arguments are archive directory and run id name

import netCDF4, sys, os, glob
from pathlib import Path
import numpy as np

archivedir = sys.argv[1]
runid = sys.argv[2]

# Check for any new complete years
flist = sorted(list(glob.glob(os.path.join(archivedir,'ocean_month.nc-[0-9]*1231'))))
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

def process(d):
    # These are on the yu_ocean grid
    ty_trans = d.variables['ty_trans'][:]
    ty_trans_gm = d.variables['ty_trans_gm'][:]

    # To match ferret behaviour, choose the region using the coordinate bounds
    # depth[k] has bounds depth_edges[k] and depth_edges[k+1]
    depth = d.variables['st_ocean'][:]
    depth_edges = d.variables['st_edges_ocean'][:]
    # Level that includes depth 3000 is first for which bottom edges is > 3000
    ilev_3000 = np.argmax(depth_edges > 3000) - 1
    ilev_5000 = np.argmax(depth_edges > 5000) - 1
    # print("Levels", ilev_3000, ilev_5000, depth[ilev_3000], depth[ilev_5000])

    lat = d.variables['grid_yu_ocean'][:]
    # lat[k] has bounds lat_edges[k] and lat_edges[k+1], except at N extreme which isn't needed here
    lat_edges = d.variables['grid_yt_ocean'][:]
    ilat_m60 = np.argmax(lat_edges > -60) - 1
    ilat_30 = np.argmax(lat_edges > 30) - 1
    ilat_60 = np.argmax(lat_edges > 60) - 1
    ilat_26 = np.argmin(abs(lat-26)) - 1
    # print("Lats", ilat_m60, ilat_30, lat[ilat_m60], lat[ilat_30])

    dmask = netCDF4.Dataset("/g/data/p66/ars599/plot_paper_data/lsmask_20110618.nc")
    # Only want the top level of the mask
    mask = dmask.variables['mask_tucell'][0]

    moc_glb = ty_trans.sum(axis=3).cumsum(axis=1) + ty_trans_gm.sum(axis=3)
    # aabwf = MOC_gbl[y=-80:-60@min,z=1:3000@min]
    aabwf = moc_glb[:,:ilev_3000+1,:ilat_m60+1].min(axis=(1,2))
    # sodc = MOC_gbl[y=-60:30@min, z=3000:5000@min]
    sodc = moc_glb[:,ilev_3000:ilev_5000+1,ilat_m60:ilat_30+1].min(axis=(1,2))

    # Atlantic defined by mask=2, 4. Want to mask out non-Atlantic points so invert this
    mask_atl = np.logical_not(np.logical_or(mask==2, mask==4))
    # Automatic broadcasting doesn't work
    # ty_trans_atl = np.ma.masked_array(ty_trans, mask_atl[np.newaxis, np.newaxis,:,:])
    mask_atl = np.broadcast_to(mask_atl, ty_trans.shape)
    ty_trans_atl = np.ma.masked_array(ty_trans, mask_atl)
    ty_trans_gm_atl = np.ma.masked_array(ty_trans_gm, mask_atl)
    moc_atl = ty_trans_atl.sum(axis=3).cumsum(axis=1) + ty_trans_gm_atl.sum(axis=3)
    # nadwf = MOC_atl[y=30:60@max, z=1:3000@max]
    nadwf = moc_atl[:,:ilev_3000+1,ilat_30:ilat_60+1].max(axis=(1,2))
    amoc26n = moc_atl[:,:ilev_3000+1,ilat_26].max(axis=1)

    return aabwf, sodc, nadwf, amoc26n

def createfile(filename,lastfile):
    # Create a new file for the global mean with time and level dimensions
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
    dnew.createVariable('nadwf', np.float32, ('time',))
    nadwf = dnew.variables['nadwf']
    nadwf.units = '1e9 kg/s'
    nadwf.standard_name = "ocean_y_overturning_mass_streamfunction"
    nadwf.long_name = "North Atlantic deep water formation"
    dnew.createVariable('amoc26n', np.float32, ('time',))
    amoc26n = dnew.variables['amoc26n']
    amoc26n.units = '1e9 kg/s'
    amoc26n.standard_name = "ocean_y_overturning_mass_streamfunction"
    amoc26n.long_name = "Maximum AMOC at 26N"
    dnew.createVariable('aabwf', np.float32, ('time',))
    aabwf = dnew.variables['aabwf']
    aabwf.units = '1e9 kg/s'
    aabwf.standard_name = "ocean_y_overturning_mass_streamfunction"
    aabwf.long_name = "Antarctic bottom water formation"
    dnew.createVariable('sodc', np.float32, ('time',))
    sodc = dnew.variables['sodc']
    sodc.units = '1e9 kg/s'
    sodc.standard_name = "ocean_y_overturning_mass_streamfunction"
    sodc.long_name = "Deep cell originating from Southern Ocean"
    for vname in ('nadwf', 'amoc26n', 'aabwf', 'sodc'):
        # Create annual mean
        var = dnew.variables[vname]
        annname = '%s_ann' % vname
        dnew.createVariable(annname, np.float32, ('year',))
        annvar = dnew.variables[annname]
        for attr in var.ncattrs():
            setattr(annvar, attr, getattr(var,attr))
    d.close()
    return dnew

fname = 'ocean_MOC_%s.nc' % runid
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
nadwf = dout.variables['nadwf']
amoc26n = dout.variables['amoc26n']
aabwf = dout.variables['aabwf']
sodc = dout.variables['sodc']

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
    # Loop is over expected years, so missing files will cause an error.
    for year in range(lastyear+1, lastfile_year+1):
        flist_year = glob.glob(os.path.join(archivedir,'ocean_month.nc-%4.4d[0-9]*' % (year)))
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
            d = netCDF4.Dataset(f)
            time_in = d.variables['time']
            time_bounds_in = d.variables[f'time_{bounds_name}']
            offset = len(time)
            # Check whether the dates match
            if offset:
                lastdate = netCDF4.num2date(time[-1], time.units, time.calendar)
                newdate = netCDF4.num2date(time_in[0], time_in.units, 'proleptic_gregorian')
                if not 25 <= (newdate-lastdate).days <= 35:
                    print("Date mismatch", lastdate, newdate, newdate-lastdate, file=sys.stderr)
                    raise Exception('Date mismatch')

            result = process(d)

            for t in range(len(time_in)):

                aabwf[offset+t] = result[0][t]*1e-9
                sodc[offset+t] = result[1][t]*1e-9
                nadwf[offset+t] = result[2][t]*1e-9
                amoc26n[offset+t] = result[3][t]*1e-9

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
