# Annual mean overturning circulation
# Optional argument nocosima means to use old 1 degree style output
# Arguments are archive directory and run id name.
# Output file is created in current working directory

import netCDF4, sys, os, glob, datetime, argparse
from pathlib import Path
import numpy as np, xarray as xr

parser = argparse.ArgumentParser(description="Calculate annual mean overtruning circulation")
parser.add_argument('--nocosima', action='store_true',
                    default=False, help="Original CM2 style multi-variable files")
parser.add_argument('archivedir', help='Archive directory')
parser.add_argument('runid', help='Run ID')
args = parser.parse_args()

# Check for any new complete years
# Assume running in 6 month chunks here.
if args.nocosima:
    # For 1 degree model
    prefix = 'ocean_month.nc-'
    flist = sorted(list(glob.glob(os.path.join(args.archivedir,prefix+'[0-9]*1231'))))
else:
    prefix = 'ocean-3d-ty_trans_rho-1-monthly-mean-ym_'
    prefix_gm = 'ocean-3d-ty_trans_rho_gm-1-monthly-mean-ym_'
    flist = sorted(list(glob.glob(os.path.join(args.archivedir,prefix+'[0-9]*_07.nc'))))
if not flist:
    print("Nothing to process")
    sys.exit(0)

lastfile = flist[-1]
if args.nocosima:
    # 1 degree
    lastfile_year = int(lastfile[-8:-4])
else:
    lastfile_year = int(lastfile[-10:-6])


def createfile(filename,lastfile):
    # Create a new file for the annual mean with time dimension
    dnew = netCDF4.Dataset(filename, 'w')
    # Any file will do for getting dimensions etc
    d = netCDF4.Dataset(lastfile)

    dnew.createDimension('nv', 2)
    time = d.variables['time']
    try:
        time_bnds = d.variables['time_bnds']
    except KeyError:
        time_bnds = d.variables['time_bounds']
    dnew.createDimension('time', None)
    dnew.createVariable('time', np.float64, ('time',))
    tnew = dnew.variables['time']
    for attr in time.ncattrs():
        # Missing value in the bounds shouldn't occur
        # Raw MOM files incorrectly set calendar as gregorian rather
        # than proleptic, but fixed in post-processing
        if attr not in ("missing_value", "_FillValue", "calendar", "calendar_type", "bounds"):
            setattr(tnew, attr, getattr(time,attr))
    tnew.bounds = "time_bnds"
    tnew.calendar = "proleptic_gregorian"
    dnew.createVariable('time_bnds', np.float64, ('time', 'nv'))
    tnew = dnew.variables['time_bnds']
    for attr in time_bnds.ncattrs():
        # Missing value in the bounds shouldn't occur
        if attr not in ("missing_value", "_FillValue"):
            setattr(tnew, attr, getattr(time_bnds,attr))

    # Variables
    dnew.createVariable('amoc', np.float32, ('time',))
    vnew = dnew.variables['amoc']
    vnew.long_name = "AMOC at 26N in density space"
    vnew.units = "Sv"
    d.close()
    return dnew

def get_AMOC(ty_trans_rho, ty_trans_rho_gm):
    # From https://github.com/wghuneke/ACCESS-CM2-025_Analysis/blob/main/notebooks/Playing_Around_ACCESS_CM2.ipynb

    latitude = 26
    xt_slice = slice(-103, -5)

    psi_amoc = ty_trans_rho.sel(grid_xt_ocean=xt_slice).sum('grid_xt_ocean').sel(method='Nearest',grid_yu_ocean=latitude)*1.0e-9
    psi_amoc_sum = psi_amoc.cumsum('potrho') - psi_amoc.sum('potrho')

    psi_amoc_gm = ty_trans_rho_gm.sel(grid_xt_ocean=xt_slice).sum('grid_xt_ocean').sel(method='Nearest',grid_yu_ocean=latitude)*1.0e-9
    psi_amoc_sum = psi_amoc_sum + psi_amoc_gm

    psi_amoc_sum = psi_amoc_sum.sel(potrho=slice(1035.5,None)).max('potrho')
    psi_amoc_sum = psi_amoc_sum.resample(time='Y').mean().compute()

    return psi_amoc_sum

fname = 'amoc_%s.nc' % args.runid
if Path(fname).exists():
    dout = netCDF4.Dataset(fname, 'r+')
else:
    dout = createfile(fname, lastfile)

time = dout.variables['time']
time_bnds = dout.variables['time_bnds']

amoc = dout.variables['amoc']
nt = len(time)

if nt > 0:
    lastdate = netCDF4.num2date(time[-1], time.units, time.calendar)
    lastyear = lastdate.year
else:
    # Set lastyear to year of first file -1 so loop starts correctly
    firstfile = flist[0]
    if args.nocosima:
        # 1 degree
        lastyear = int(firstfile[-8:-4]) - 1
    else:
        lastyear = int(firstfile[-10:-6]) - 1

print("LASTYEAR", lastyear, lastfile_year)

if lastfile_year > lastyear:
    print('Data to process', args.runid, lastyear+1, lastfile_year)
    # Loop is over expected years, so missing files will cause an error.
    for year in range(lastyear+1, lastfile_year+1):
        print(year)
        d = xr.open_mfdataset(os.path.join(args.archivedir,'%s%4.4d*' % (prefix,year)),
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

        if args.nocosima:
            amoc[offset] = get_AMOC(d.ty_trans_rho, d.ty_trans_rho_gm).values[0]
        else:
            d_gm = xr.open_mfdataset(os.path.join(args.archivedir,'%s%4.4d*' % (prefix_gm,year)),
                    combine='by_coords', use_cftime=True)
            amoc[offset] = get_AMOC(d.ty_trans_rho, d_gm.ty_trans_rho_gm).values[0]
        date0 = datetime.datetime(year, 1, 1)
        date1 = datetime.datetime(year+1, 1, 1)
        time_bnds[offset,0] = netCDF4.date2num(date0, units=time.units, calendar=time.calendar)
        time_bnds[offset,1] = netCDF4.date2num(date1, units=time.units, calendar=time.calendar)
        time[offset] = time_bnds[offset].mean()

        d.close()
        dout.sync()  # Sync to disk once per year

dout.close()
