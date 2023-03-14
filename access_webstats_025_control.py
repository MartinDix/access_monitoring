# ACCESS-CM2 PD control with 1 and 0.25 degree ocean

import iris
import matplotlib.pyplot as plt
import iris.plot as iplt
import warnings, datetime, cf_units, datetime, cftime, os, argparse, nc_time_axis, numpy as np
import iris_tools as itools

parser = argparse.ArgumentParser(description="Plot CMIP6 runs")
parser.add_argument('--minyr', dest='minyr', type=int,
                    default=1, help="First year to plot")
parser.add_argument('--maxyr', dest='maxyr', type=int,
                    default=1400, help="Last year to plot")
parser.add_argument('--savefig', dest='savefig', action='store_true',
                    default=False, help="Save plot to file")
parser.add_argument('--savedir', dest='savedir',
                    default='/g/data/p66/accessdev-web/mrd599/access-cm2/control025',
                    help='Directory to save plot')
args = parser.parse_args()

plt.style.use('mrd_presentation')

lw=1.5
# runs = ['bz687', 'ch495', 'cj877', 'cq880']
# Plot in this order so that the 2 0.25 degree runs are more clearly separated by color
runs = ['cj877', 'cq880', 'bz687']
# label = {'bz687':'1$\degree$', 'ch495':'0.25$\degree$', 'cj877':'0.25$\degree$ TOPO5'}
label = {'bz687':'1$\degree$', 'ch495':'0.25$\degree$ old topo', 'cj877':'0.25$\degree$',
         'cq880':'0.25$\degree$ new GM'}
# Put 1 degree at the back to make the 0.25 results clearer
zorder = {'cj877':1, 'cq880':2, 'bz687':0}

def savefig(filename):
    if args.savefig:
        plt.savefig(os.path.join(args.savedir,filename), dpi=150, bbox_inches='tight',
                    facecolor='white',
                    metadata={"History": "%s created by access_webstats_025_control.py at %s" %
                    (filename, datetime.datetime.today().strftime("%Y-%m-%d %H:%M"))})
    else:
        plt.show()
    plt.close(plt.gcf())

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

cubes = {}
for run in runs:
    with warnings.catch_warnings():
        # Warnings about masked coordinates
        warnings.simplefilter("ignore")
        clist = iris.load(f'/g/data/p66/mrd599/access_stats/{run}/{run}_ann.nc')
    cd = {c.var_name:c for c in clist}
    cubes[run] = cd

# Plot time axis is relative to 01-jan-2000
taxis = cf_units.Unit("days since 2000-01-01 00:00", calendar="proleptic_gregorian")
t0 = taxis.date2num(datetime.datetime(args.minyr,1,1,0,0,0))
t1 = taxis.date2num(datetime.datetime(args.maxyr,1,1,0,0,0))

t_start = cftime.DatetimeProlepticGregorian(args.minyr, 1, 1)
t_end = cftime.DatetimeProlepticGregorian(args.maxyr, 12, 31)
time_selection = iris.Constraint(time=lambda c: t_start <= c.point <= t_end)

ny = args.maxyr - args.minyr + 1  # Number of years on graph
if ny >= 1000:
    tickfreq = 200
elif ny >= 600:
    tickfreq = 100
elif ny >= 200:
    tickfreq = 50
elif ny >= 100:
    tickfreq = 20
else:
    tickfreq = 10
ticks = np.array([taxis.date2num(cftime.DatetimeProlepticGregorian(y,1,1,0,0,0))
                for y in range(max(tickfreq,args.minyr),args.maxyr+1,tickfreq)])
# With nc_time_axis >= 1.4
year_fmt = nc_time_axis.CFTimeFormatter("%Y", calendar='proleptic_gregorian')


# Subset the data to the given year range. If simply use set_xlim then
# matplotlib doesn't set the y range correctly.
for run in runs:
    for var in ["ts", "tas", "ts_sea", "rlut", "rsdt", "rsut"]:
        cubes[run][var] = cubes[run][var].extract(time_selection)
        # A single matching time gets removed as a dimension so need the
        # extra check on the length of the time coordinate.
        if cubes[run][var] and len(cubes[run][var].coord('time').points) == 1:
            cubes[run][var] = None


# Global mean surface air temperature

fig, axes = plt.subplots()

var = 'tas'
for run in runs:
    # May have an empty cube because of time constraint.
    if cubes[run][var]:
        iplt.plot(itools.global_mean(cubes[run][var]), label=label[run], zorder=zorder[run])

axes.set_title('Global mean surface air temperature')
axes.set_xlabel('Year')
axes.set_ylabel('K')
set_plot_props(axes)
savefig('tas.png')

# Hemispheric mean surface temperature

fig, axes = plt.subplots()

var = 'tas'
for k, run in enumerate(runs):
    color = f'C{k%10}'
    if cubes[run][var]: # May have an empty cube because of time constraint
        iplt.plot(itools.NH_mean(cubes[run][var]),linewidth=lw, label=label[run], zorder=zorder[run], color=color)
for k, run in enumerate(runs):
    color = f'C{k%10}'
    if cubes[run][var]: # May have an empty cube because of time constraint
        iplt.plot(itools.SH_mean(cubes[run][var]),linewidth=lw, dashes=(5,1), color=color)

axes.set_title('Hemispheric mean surface air temperature')
axes.set_xlabel('Year')
axes.set_ylabel('K')
set_plot_props(axes)
savefig('ts_hem.png')


# TOA net flux

fig, axes = plt.subplots()

for k, run in enumerate(runs):
    c = cubes[run]
    if not c['rsdt']:
        # May be empty because of time subsetting
        continue
    net = c['rsdt'] - c['rsut'] - c['rlut']
    netg = itools.global_mean(net)
    color = f'C{k%10}'
    iplt.plot(netg, linewidth=0.5,alpha=0.5, color=color)
    # if len(netg.data) >= 10:
    #     netg_mean = netg.rolling_window('time', iris.analysis.MEAN, 10)
    #     iplt.plot(netg_mean, linewidth=lw, label=label[run], zorder=zorder[run], color=color)
    # else:
    # Temporarily plot the annual values
    iplt.plot(netg, linewidth=lw, label=label[run], zorder=zorder[run], color=color)

axes.set_title('TOA net flux')
axes.set_xlabel('Year')
axes.set_ylabel('W/m$^2$')
set_plot_props(axes)
axes.grid(True)
savefig('net.png')

fig, axes = plt.subplots()

var = 'rlut'
for run in runs:
    if cubes[run][var]: # May have an empty cube because of time constraint
        iplt.plot(itools.global_mean(cubes[run][var]), linewidth=lw, label=label[run], zorder=zorder[run])

axes.set_title('Global mean RLUT')
axes.set_xlabel('Year')
axes.set_ylabel('W/m$^2$')
set_plot_props(axes)
savefig('rlut.png')

fig, axes = plt.subplots()

for run in runs:
    if cubes[run][var]: # May have an empty cube because of time constraint
        iplt.plot(itools.global_mean(cubes[run]['rsdt']-cubes[run]['rsut']), linewidth=lw, label=label[run], zorder=zorder[run])

axes.set_title('Global mean TOA net SW')
axes.set_xlabel('Year')
axes.set_ylabel('W/m$^2$')
set_plot_props(axes)
savefig('toa_net_SW.png')

# Ocean model
ocubes = {}
for run in runs:
    with warnings.catch_warnings():
        # Warnings about masked coordinates
        warnings.simplefilter("ignore")
        if run == 'bz687':
            clist = iris.load(f'/g/data/p66/mrd599/access_stats/{run}/ocean_mean_{run}.nc')
        else:
            # 0.25 degree run has a different format
            clist = iris.load(f'/g/data/p66/mrd599/access_stats/{run}/ocean_stats_{run}.nc')
    # Apply time selection
    cd = {c.var_name:c.extract(time_selection) for c in clist}
    ocubes[run] = cd

fig, axes = plt.subplots()
tzero = 273.15
for run in runs:
    try:
        # 0.25 in degrees C and only single level
        iplt.plot(ocubes[run]['sst_ann'], linewidth=lw, label=label[run], zorder=zorder[run])
    except KeyError:
        # SST is top level of the temperature.
        iplt.plot(ocubes[run]['temp_ann'][:,0]-tzero, linewidth=lw, label=label[run], zorder=zorder[run])

axes.set_title('Global mean SST')
axes.set_xlabel('Year')
axes.set_ylabel('degrees C')
set_plot_props(axes)
savefig('sst.png')

# Sea surface salinity
fig, axes = plt.subplots()
for run in runs:
    try:
        iplt.plot(ocubes[run]['sss_ann'], linewidth=lw, label=label[run], zorder=zorder[run])
    except KeyError:
        iplt.plot(ocubes[run]['salt_ann'][:,0], linewidth=lw, label=label[run], zorder=zorder[run])

axes.set_title('Global mean sea surface salinity')
axes.set_xlabel('Year')
axes.set_ylabel('PSU')
set_plot_props(axes)
savefig('sss.png')

# 3D mean ocean temperature
fig, axes = plt.subplots()
for run in runs:
    try:
        iplt.plot(ocubes[run]['temp_zmean_ann']-tzero, linewidth=lw, label=label[run], zorder=zorder[run])
    except KeyError:
        iplt.plot(ocubes[run]['temp_ann'], linewidth=lw, label=label[run], zorder=zorder[run])

axes.set_title('3D ocean mean temperature')
axes.set_xlabel('Year')
axes.set_ylabel('degrees C')
set_plot_props(axes)
savefig('ocean_3d_temp.png')

# Sea level
fig, axes = plt.subplots()
for run in runs:
    iplt.plot(ocubes[run]['sea_level_ann'], linewidth=lw, label=label[run], zorder=zorder[run])

axes.set_title('Sea level')
axes.set_xlabel('Year')
axes.set_ylabel('m')
set_plot_props(axes)
savefig('sea_level.png')

# Drake passage
fig, axes = plt.subplots()
for run in runs:
    drake = iris.load_cube(f'/g/data/p66/mrd599/access_stats/{run}/drake_passage_{run}.nc')
    drake = drake.extract(time_selection)
    iplt.plot(drake, linewidth=lw, label=label[run], zorder=zorder[run])

axes.set_title('Drake passage transport')
axes.set_xlabel('Year')
axes.set_ylabel('Sv')
set_plot_props(axes)
savefig('drake_passage.png')

# AMOC
fig, axes = plt.subplots()
for run in runs:
    amoc = iris.load_cube(f'/g/data/p66/mrd599/access_stats/{run}/amoc_{run}.nc')
    amoc = amoc.extract(time_selection)
    iplt.plot(amoc, linewidth=lw, label=label[run], zorder=zorder[run])

axes.set_title('AMOC at 26N')
axes.set_xlabel('Year')
axes.set_ylabel('Sv')
set_plot_props(axes)
savefig('amoc.png')



icubes = {}
for run in runs:
    with warnings.catch_warnings():
        # Warnings about masked coordinates
        warnings.simplefilter("ignore")
        clist = iris.load(f'/g/data/p66/mrd599/access_stats/{run}/ice_mean_{run}.nc')
    # Apply time selection
    cd = {c.var_name:c.extract(time_selection) for c in clist}
    icubes[run] = cd

for var in ('area', 'vol'):
    for month in [2,8]:
        if month==2:
            month_name = 'March'
        else:
            month_name = 'September'
        for hem in ('NH', 'SH'):
            fig, axes = plt.subplots()
            for run in runs:
                vname = f'ice_{var}_{hem.lower()}'
                if vname in icubes[run]:
                    area = icubes[run][vname][month::12] * 1e-12
                    l = iplt.plot(area, linewidth=lw, label=label[run], zorder=zorder[run])

            axes.set_xlabel('Year')
            if var == 'area':
                axes.set_title(f'Sea ice area {month_name} {hem}')
                axes.set_ylabel('Million square km')
            else:
                axes.set_title(f'Sea ice volume {month_name} {hem}')
                axes.set_ylabel('1000 km$^3$')
            set_plot_props(axes)
            savefig(f"ice_{var}_{month_name}_{hem}.png")

# ke_runs = ['PIspinup','bz687']
# # Ocean KE
# cubes = {}
# for run in ke_runs:
#     with warnings.catch_warnings():
#         # Warnings about masked coordinates
#         warnings.simplefilter("ignore")
#         clist = iris.load('/g/data/p66/mrd599/access_stats/{0}/ocean_ke_{0}.nc'.format(run))
#     cd = {c.var_name:c for c in clist}
#     cubes[run] = cd


# fig, axes = plt.subplots()

# for run in ke_runs:
#     iplt.plot(cubes[run]['ke_tot_ann'], linewidth=lw, label=label[run], zorder=zorder[run])

# axes.set_title('Global mean ocean KE')
# axes.set_xlabel('Year')
# axes.set_ylabel('1e15 J')
# set_plot_props(axes)
# savefig('ocean_ke.png')


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
