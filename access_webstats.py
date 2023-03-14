# Monitoring graphs from ACCESS-CM2 simulations

import iris
import matplotlib.pyplot as plt
import iris.plot as iplt
import warnings, datetime, cf_units, datetime, cftime, os, argparse, numpy as np, nc_time_axis
import iris_tools as itools
import json

parser = argparse.ArgumentParser(description="Plot CMIP6 runs")
parser.add_argument('--minyr', dest='minyr', type=int,
                    default=1850, help="First year to plot")
parser.add_argument('--maxyr', dest='maxyr', type=int,
                    default=2000, help="Last year to plot")
parser.add_argument('--modelyears', dest='modelyears', action='store_true',
                    default=False, help="Use model years with no offsets (for control runs)")
parser.add_argument('--savefig', dest='savefig', action='store_true',
                    default=False, help="Save plot to file")
parser.add_argument('--savedir', dest='savedir',
                    default='/g/data/p66/accessdev-web/mrd599/access-cm2/control',
                    help='Directory to save plot')
parser.add_argument('--fade', dest='fade', nargs='+', default=[], help='Fade these runs in plot')
parser.add_argument('runs', nargs='+', help='Run names for plotting')
args = parser.parse_args()

plt.style.use('mrd_presentation')

dpi=150
lw=1.5
# Work around to allow comments. Use json5 instead?
with open('/g/data/p66/accessdev-web/mrd599/access-cm2/exptIDs.json') as jsonfile:
    jsondata = ''.join(line for line in jsonfile if not line.strip().startswith('//'))
    label = json.loads(jsondata)
color = {}
k=0
for run in args.runs:
    if run in args.fade:
        color[run] = 'grey'
    else:
        color[run] = f'C{k%10}'
        k += 1

def shift_time(cube):
    if args.modelyears:
        return
    pi_start = 950 # Start year of official piControl
    time = cube.coord('time')
    if time.units.num2date(time.points[0]).year < 1850:
        # Needs to be shifted
        tpi = cftime.DatetimeProlepticGregorian(pi_start,1,1)
        this = cftime.DatetimeProlepticGregorian(1850,1,1)
        delta = this-tpi
        time.points = time.points + delta.days
        time.bounds = time.bounds + delta.days

def set_plot_props(axes):
    axes.set_xticks(ticks)
    axes.xaxis.set_major_formatter(year_fmt)
    axes.set_xlim(t0,t1)
    legend_info = axes.get_legend_handles_labels()
    if len(legend_info[1]) > 6:
        axes.legend(fontsize=8,ncol=2)
    else:
        axes.legend()
    axes.grid(True)

def savefig(filename):
    if args.savefig:
        plt.savefig(os.path.join(args.savedir,filename), dpi=150, bbox_inches='tight',
            facecolor='white',
            metadata={"History": "%s created by access_webstats.py at %s" %
            (filename, datetime.datetime.today().strftime("%Y-%m-%d %H:%M"))})
    else:
        plt.show()
    plt.close(plt.gcf())

def setfade(l, run):
    # iplt.plt returns a list
    if run in args.fade:
        l[0].set_linewidth(1)

# Plot time axis is relative to 01-jan-2000
taxis = cf_units.Unit("days since 2000-01-01 00:00", calendar="proleptic_gregorian")
t0 = taxis.date2num(datetime.datetime(args.minyr,1,1,0,0,0))
t1 = taxis.date2num(datetime.datetime(args.maxyr,1,1,0,0,0))
year_fmt = nc_time_axis.CFTimeFormatter("%Y", calendar='proleptic_gregorian')

t_start = cftime.DatetimeProlepticGregorian(args.minyr, 1, 1)
t_end = cftime.DatetimeProlepticGregorian(args.maxyr, 12, 31)
time_selection = iris.Constraint(time=lambda c: t_start <= c.point <= t_end)

if args.maxyr - args.minyr > 1000:
    tickfreq = 200
elif args.maxyr - args.minyr > 600:
    tickfreq = 100
elif args.maxyr - args.minyr > 200:
    tickfreq = 50
else:
    tickfreq = 20
ticks = np.array([taxis.date2num(cftime.DatetimeProlepticGregorian(y,1,1,0,0,0))
                for y in range(max(tickfreq,args.minyr),args.maxyr+1,tickfreq)])

# Atmospheric model fields

cubes = {}
for run in args.runs:
    with warnings.catch_warnings():
        # Warnings about masked coordinates
        warnings.simplefilter("ignore")
        clist = iris.load(f'/g/data/p66/mrd599/access_stats/{run}/{run}_ann.nc',
                          ("tas", "rlut", "rsdt", "rsut") )
    cd = {}
    for c in clist:
        shift_time(c)
        # Subset the data to the given year range. If simply use set_xlim then
        # matplotlib doesn't set the y range correctly.
        c = c.extract(time_selection)
        # A single matching time gets removed as a dimension so need the
        # extra check on the length of the time coordinate.
        if c and len(c.coord('time').points) > 1:
            cd[c.var_name] = c
    cubes[run] = cd

# Global mean surface air temperature

fig, axes = plt.subplots()
var = 'tas'
for run in args.runs:
    if var in cubes[run]: # May have a missing cube because of time constraint.
        l = iplt.plot(itools.global_mean(cubes[run][var]), linewidth=lw, color=color[run], label=label[run])
        setfade(l, run)

axes.set_title('Global mean surface air temperature')
axes.set_xlabel('Year')
axes.set_ylabel('K')
set_plot_props(axes)
savefig('tas.png')

# Hemispheric mean surface temperature

fig, axes = plt.subplots()
var = 'tas'
for run in args.runs:
    if var in cubes[run]:
        l = iplt.plot(itools.NH_mean(cubes[run][var]), linewidth=lw, label=label[run], color=color[run])
        setfade(l, run)
        l = iplt.plot(itools.SH_mean(cubes[run][var]), linewidth=lw, dashes=(5,1), color=color[run])
        setfade(l, run)

axes.set_title('Hemispheric mean surface air temperature')
axes.set_xlabel('Year')
axes.set_ylabel('K')
set_plot_props(axes)
savefig('ts_hem.png')

# TOA net flux

fig, axes = plt.subplots()
for run in args.runs:
    c = cubes[run]
    if 'rsdt' not in c:
        # May be missing because of time subsetting
        continue
    net = c['rsdt'] - c['rsut'] - c['rlut']
    netg = itools.global_mean(net)
    if args.maxyr - args.minyr > 300:
        iplt.plot(netg, linewidth=0.5,alpha=0.5, color=color[run])
    else:
        iplt.plot(netg, linewidth=1,alpha=0.75, color=color[run])
    if len(netg.data) >= 10:
        netg_mean = netg.rolling_window('time', iris.analysis.MEAN, 10)
        l = iplt.plot(netg_mean, linewidth=lw, label=label[run], color=color[run])
        setfade(l, run)
    else:
        # Temporarily plot the annual values
        l = iplt.plot(netg, linewidth=lw, label=label[run], color=color[run])
        setfade(l, run)

axes.set_title('TOA net flux')
axes.set_xlabel('Year')
axes.set_ylabel('W/m$^2$')
set_plot_props(axes)
savefig('net.png')

fig, axes = plt.subplots()
var = 'rlut'
for run in args.runs:
    if var in cubes[run]:
        l = iplt.plot(itools.global_mean(cubes[run][var]), linewidth=lw, label=label[run], color=color[run])
        setfade(l, run)

axes.set_title('Global mean RLUT')
axes.set_xlabel('Year')
axes.set_ylabel('W/m$^2$')
set_plot_props(axes)
savefig('rlut.png')

fig, axes = plt.subplots()
for run in args.runs:
    if var in cubes[run]:
        l = iplt.plot(itools.global_mean(cubes[run]['rsdt']-cubes[run]['rsut']), linewidth=lw, color=color[run], label=label[run])
        setfade(l, run)

axes.set_title('Global mean TOA net SW')
axes.set_xlabel('Year')
axes.set_ylabel('W/m$^2$')
set_plot_props(axes)
savefig('toa_net_SW.png')

# Ocean model
ocubes = {}
for run in args.runs:
    with warnings.catch_warnings():
        # Warnings about masked coordinates
        warnings.simplefilter("ignore")
        clist = iris.load(f'/g/data/p66/mrd599/access_stats/{run}/ocean_mean_{run}.nc',
           ('temp_ann', 'salt_ann', 'temp_zmean_ann', 'sea_level_ann', 'acc_drake_ann'))
    cd = {}
    for c in clist:
        shift_time(c)
        c = c.extract(time_selection)
        if c and len(c.coord('time').points) > 1:
            cd[c.var_name] = c
    ocubes[run] = cd

# SST is top level of the temperature. Note that the older data is in degrees C rather than K
fig, axes = plt.subplots()
tzero = 273.15
for run in args.runs:
    if 'temp_ann' in ocubes[run]:
        if ocubes[run]['temp_ann'][0,0].data > tzero:
            l = iplt.plot(ocubes[run]['temp_ann'][:,0]-tzero, linewidth=lw, color=color[run], label=label[run])
        else:
            l = iplt.plot(ocubes[run]['temp_ann'][:,0], linewidth=lw, color=color[run], label=label[run])
        setfade(l, run)

axes.set_title('Global mean SST')
axes.set_xlabel('Year')
axes.set_ylabel('degrees C')
set_plot_props(axes)
savefig('sst.png')

# Sea surface salinity
fig, axes = plt.subplots()
for run in args.runs:
    if 'salt_ann' in ocubes[run]:
        l = iplt.plot(ocubes[run]['salt_ann'][:,0], linewidth=lw, color=color[run], label=label[run])
        setfade(l, run)

axes.set_title('Global mean sea surface salinity')
axes.set_xlabel('Year')
axes.set_ylabel('PSU')
set_plot_props(axes)
savefig('sss.png')

# 3D mean ocean temperature
fig, axes = plt.subplots()
for run in args.runs:
    if 'temp_zmean_ann' in ocubes[run]:
        if ocubes[run]['temp_zmean_ann'][0].data > tzero:
            l = iplt.plot(ocubes[run]['temp_zmean_ann']-tzero, linewidth=lw, color=color[run], label=label[run])
        else:
            l = iplt.plot(ocubes[run]['temp_zmean_ann'], linewidth=lw, color=color[run], label=label[run])
        setfade(l, run)

axes.set_title('3D ocean mean temperature')
axes.set_xlabel('Year')
axes.set_ylabel('degrees C')
set_plot_props(axes)
savefig('ocean_3d_temp.png')

# Sea level
fig, axes = plt.subplots()
for run in args.runs:
    if 'sea_level_ann' in ocubes[run]:
        l = iplt.plot(ocubes[run]['sea_level_ann'], linewidth=lw, color=color[run], label=label[run])
        setfade(l, run)

axes.set_title('Sea level')
axes.set_xlabel('Year')
axes.set_ylabel('m')
set_plot_props(axes)
savefig('sea_level.png')

fig, axes = plt.subplots()
for run in args.runs:
    if 'acc_drake_ann' in ocubes[run]:
        l = iplt.plot(ocubes[run]['acc_drake_ann'], linewidth=lw, color=color[run], label=label[run])
        setfade(l, run)

axes.set_title('Drake passage transport')
axes.set_xlabel('Year')
axes.set_ylabel('Sv')
set_plot_props(axes)
savefig('drake_passage.png')

# Meriodional overturning
mcubes = {}
for run in args.runs:
    with warnings.catch_warnings():
        # Warnings about masked coordinates
        warnings.simplefilter("ignore")
        clist = iris.load(f'/g/data/p66/mrd599/access_stats/{run}/ocean_MOC_{run}.nc')
    cd = {}
    for c in clist:
        shift_time(c)
        c = c.extract(time_selection)
        if c and len(c.coord('time').points) > 1:
            cd[c.var_name] = c
    mcubes[run] = cd

for var in ('nadwf_ann', 'amoc26n_ann', 'aabwf_ann', 'sodc_ann'):
    fig, axes = plt.subplots()
    for run in args.runs:
        if var in mcubes[run]:
            l = iplt.plot(mcubes[run][var], linewidth=lw, color=color[run], label=label[run])
            setfade(l, run)

    axes.set_title(mcubes[run][var].long_name)
    axes.set_xlabel('Year')
    axes.set_ylabel('Sv')
    set_plot_props(axes)
    # Remove the _ann from var for the filename
    savefig(f'{var[:-4]}.png')

icubes = {}
for run in args.runs:
    with warnings.catch_warnings():
        # Warnings about masked coordinates
        warnings.simplefilter("ignore")
        clist = iris.load(f'/g/data/p66/mrd599/access_stats/{run}/ice_mean_{run}.nc')
    cd = {}
    for c in clist:
        shift_time(c)
        c = c.extract(time_selection)
        if c and len(c.coord('time').points) > 1:
            cd[c.var_name] = c
    icubes[run] = cd

for var in ('area', 'vol'):
    for month in [2,8]:
        if month==2:
            month_name = 'March'
        else:
            month_name = 'September'
        for hem in ('NH', 'SH'):
            fig, axes = plt.subplots()
            for run in args.runs:
                vname = f'ice_{var}_{hem.lower()}'
                if vname in icubes[run]:
                    area = icubes[run][vname][month::12] * 1e-12
                    l = iplt.plot(area, linewidth=lw, color=color[run], label=label[run])
                    setfade(l, run)

            axes.set_xlabel('Year')
            if var == 'area':
                axes.set_title(f'Sea ice area {month_name} {hem}')
                axes.set_ylabel('Million square km')
            else:
                axes.set_title(f'Sea ice volume {month_name} {hem}')
                axes.set_ylabel('1000 km$^3$')
            set_plot_props(axes)
            savefig(f"ice_{var}_{month_name}_{hem}.png")

# Ocean KE
fig, axes = plt.subplots()
for run in args.runs:
    with warnings.catch_warnings():
        # Warnings about masked coordinates
        warnings.simplefilter("ignore")
        c = iris.load_cube('/g/data/p66/mrd599/access_stats/{0}/ocean_ke_{0}.nc'.format(run),
        'ke_tot_ann')
    shift_time(c)
    c = c.extract(time_selection)
    if c and len(c.coord('time').points) > 1:
        l = iplt.plot(c, linewidth=lw, color=color[run], label=label[run])
        setfade(l, run)

axes.set_title('Global mean ocean KE')
axes.set_xlabel('Year')
axes.set_ylabel('1e15 J')
set_plot_props(axes)
savefig('ocean_ke.png')

# Set up index.html from template.html with current datetime
# Note that cylc in UTC mode resets the time zone so be explicit
# about using UTC
now = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
f = open(os.path.join(args.savedir,"index.html"), 'w')
for l in open(os.path.join(args.savedir,"template.html")).readlines():
    if l.startswith("DATETIME"):
        f.write(f"<p>Generated at {now}</p>\n")
    else:
        f.write(l)
f.close()
