# Calculate global means of the ocean potential temperature with
# volume weighting.
# Iris seems to get confused by the 2D lat/lon coordinates so do a
# simple numpy sum

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

def createfile(filename,lastfile):
    # Create a new file for the global mean with time and level dimensions
    dnew = netCDF4.Dataset(filename, 'w')
    # Any file will do for getting dimensions etc
    d = netCDF4.Dataset(lastfile)

    # Dimensions and coordinate variables
    lev = d.variables['st_ocean']
    dnew.createDimension('nv', 2)
    dnew.createDimension('lev', len(lev))
    dnew.createVariable('lev', np.float64, ('lev',))
    lnew = dnew.variables['lev']
    lnew[:] = lev[:]
    lnew.units = 'm'
    lnew.positive = 'down'
    lnew.long_name = 'tcell zstar depth'
    lnew.cartesian_axis = 'Z'
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
    dnew.createVariable('temp', np.float32, ('time', 'lev'))
    dnew.createVariable('temp_zmean', np.float32, ('time'))
    tmean = dnew.variables['temp']
    tmean.units = 'K'
    tmean.standard_name = "sea_water_conservative_temperature"
    t_zmean = dnew.variables['temp_zmean']
    t_zmean.units = 'K'
    t_zmean.standard_name = "sea_water_conservative_temperature"
    dnew.createVariable('salt', np.float32, ('time', 'lev'))
    dnew.createVariable('salt_zmean', np.float32, ('time'))
    smean = dnew.variables['salt']
    smean.units = 'psu'
    smean.standard_name = "sea_water_salinity"
    s_zmean = dnew.variables['salt_zmean']
    s_zmean.units = 'psu'
    s_zmean.standard_name = "sea_water_salinity"
    dnew.createVariable('sea_level', np.float32, ('time'))
    sea_level = dnew.variables['sea_level']
    sea_level.units = 'm'
    sea_level.standard_name = "sea_surface_height_above_geoid"
    dnew.createVariable('sst_nh', np.float32, ('time'))
    sst_nh = dnew.variables['sst_nh']
    sst_nh.units = 'K'
    sst_nh.standard_name = "sea_surface_temperature"
    dnew.createVariable('sst_sh', np.float32, ('time'))
    sst_sh = dnew.variables['sst_sh']
    sst_sh.units = 'K'
    sst_sh.standard_name = "sea_surface_temperature"
    dnew.createVariable('sss_nh', np.float32, ('time'))
    sss_nh = dnew.variables['sss_nh']
    sss_nh.units = 'psu'
    sss_nh.standard_name = "sea_surface_salinity"
    dnew.createVariable('sss_sh', np.float32, ('time'))
    sss_sh = dnew.variables['sss_sh']
    sss_sh.units = 'psu'
    sss_sh.standard_name = "sea_surface_salinity"
    dnew.createVariable('ohc', np.float32, ('time'))
    ohc = dnew.variables['ohc']
    ohc.units = '1e24 J'
    ohc.long_name = "ocean_heat_content"
    dnew.createVariable('acc_drake', np.float32, ('time'))
    acc_drake = dnew.variables['acc_drake']
    acc_drake.units = 'Sv'
    acc_drake.long_name = "ACC in Drake passage"
    for vname in ('temp', 'temp_zmean', 'salt', 'salt_zmean', 'sea_level',
                  'sst_nh', 'sst_sh', 'sss_nh', 'sss_sh', 'ohc', 'acc_drake'):
        # Create annual mean
        var = dnew.variables[vname]
        annname = '%s_ann' % vname
        if len(var.shape) == 1:
            dnew.createVariable(annname, np.float32, ('year',))
        else:
            dnew.createVariable(annname, np.float32, ('year', 'lev'))
        annvar = dnew.variables[annname]
        for attr in var.ncattrs():
            setattr(annvar, attr, getattr(var,attr))
    d.close()
    return dnew

fname = 'ocean_mean_%s.nc' % runid
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
temp = dout.variables['temp']
temp_zmean = dout.variables['temp_zmean']
salt = dout.variables['salt']
salt_zmean = dout.variables['salt_zmean']
sea_level = dout.variables['sea_level']
sst_nh = dout.variables['sst_nh']
sst_sh = dout.variables['sst_sh']
sss_nh = dout.variables['sss_nh']
sss_sh = dout.variables['sss_sh']
ohc = dout.variables['ohc']
acc_drake = dout.variables['acc_drake']

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
            print(f)
            d = netCDF4.Dataset(f)
            time_in = d.variables['time']
            time_bounds_in = d.variables[f'time_{bounds_name}']
            temp_in = d.variables['temp']
            salt_in = d.variables['salt']
            sea_level_in = d.variables['sea_level']
            trans_in = d.variables['tx_trans']
            dht = d.variables['dht']
            area = d.variables['area_t']
            lat = d.variables['grid_yt_ocean']
            # Find limit of SH by counting latitudes < 0
            sh_limit = (lat[:] < 0).sum()
            area_sh = np.array(area[:])
            area_sh[sh_limit:] = 0.
            area_nh = np.array(area[:])
            area_nh[:sh_limit] = 0.
            nlev = temp.shape[1]

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
                weights3 = np.zeros(temp_in[0].shape)
                for k in range(nlev):
                    weights3[k] = dht[t,k] * area
                # Create a temporary copy so that it's only uncompressed once
                tempx = np.ma.array(temp_in[t])
                for k in range(nlev):
                    temp[offset+t,k] = np.ma.average(tempx[k], weights=weights3[k])
                temp_zmean[offset+t] = np.ma.average(tempx, weights=weights3)
                saltx = np.ma.array(salt_in[t])
                for k in range(nlev):
                    salt[offset+t,k] = np.ma.average(saltx[k], weights=weights3[k])
                salt_zmean[offset+t] = np.ma.average(saltx, weights=weights3)
                sea_level[offset+t] = np.ma.average(sea_level_in[t], weights=area[:])
                sst_nh[offset+t] = np.ma.average(temp_in[t,0], weights=area_nh)
                sst_sh[offset+t] = np.ma.average(temp_in[t,0], weights=area_sh)
                sss_nh[offset+t] = np.ma.average(salt_in[t,0], weights=area_nh)
                sss_sh[offset+t] = np.ma.average(salt_in[t,0], weights=area_sh)
                cp = 3992   # The sea water specific heat capacity (J(kg*K))
                rho = 1035  # Sea water density is roughly (kg/m^3)
                ohc[offset+t] = 1e-24 * cp * rho * np.ma.sum(tempx*weights3)
                # Dave's ferret calculation
                # acc_drake = TX_TRANS[D=1,I=213,J=33:50@SUM,K=@SUM]/1e9
                acc_drake[offset+t] = 1e-9 * np.ma.sum(trans_in[t,:,32:50,212])

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
