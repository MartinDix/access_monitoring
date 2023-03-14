# ACCESS Climate model monitoring scripts

These scripts implement a simple monitoring scheme for ACCESS-CM2 simulations on gadi

It's driven by a cylc suite u-be329 that runs once per day (included here for reference).

CM2 runs save data after each cycle to an archive directory in /scratch or /g/data. This has a structure like (note u- part of suite name is dropped)
```
archive/cu929
   history
      atm
      cpl
      ice
      ocan
  restart
      atm
      cpl
      ice
      ocn
```

Model runs to be processed are set in the ACTIVE_RUNS list in the suite rose-suite.conf file. For each of these runs the suite scripts `um2netcdf.sh` and `coupled_update.sh` are run.

For newer CM2 runs that include post-processing to netCDF as part of the suite `um2netcdf.sh` simply calls `access_means_update_new.py` with a list of variable IDs. For older runs it converts selected monthly mean fields to netCDF and then calls `access_means_update.py`. The update script checks for new model output and appends to files of monthly means and also calculates annual means.

There's a similar process for ocean model output. However standard CM2 runs put all monthly means in a single file while COSIMA style output used in the 0.25 degree runs uses separate files for each variable. This requires different post-processing scripts. There's also a slightly different collection of diagnostics in each case.

There are two plotting scripts, `access_webstats.py` and `access_webstats_025_control.py`. The way these are run is controlled by `suite.rc` (e.g. to set different time ranges for different classes of simulations). These scripts calculate global means and create the graphs in `/g/data/p66/accessdev-web/` which is accessible via accessdev. This uses a pre-existing template for the `index.html` file.

## ISSUES
The accumulated means files were originally designed for a slightly different purpose and include variables that aren't plotted. The processing should calculate the global means rather than leaving it to the plotting script. This could then simply plot everything it finds.
