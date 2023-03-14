set -e
RUN=$1

OUTPUT_DIR=/g/data/p66/mrd599/access_stats/${RUN}

mkdir -p $OUTPUT_DIR

cd $OUTPUT_DIR

# On gadi, ARCHIVEDIR may be either in /g/data or /scratch
ARCHIVEDIR=/g/data/$PROJECT/$USER/archive
if [ ! -d $ARCHIVEDIR/$RUN ]; then
    ARCHIVEDIR=/scratch/$PROJECT/$USER/archive
fi

# For the 0.25 degree model skip the ocean updates
# These runs use COSIMA diag_table with files ocean-scalar-1-daily*
if compgen -G $ARCHIVEDIR/$RUN/history/ocn/ocean-scalar-1-daily* > /dev/null; then
    python ~/src/python/access/ocean_stats_COSIMA.py $ARCHIVEDIR/$RUN/history/ocn $RUN
    python ~/src/python/access/ocean_transport_update.py $ARCHIVEDIR/$RUN/history/ocn $RUN
    python ~/src/python/access/ocean_overturning_update.py $ARCHIVEDIR/$RUN/history/ocn $RUN
    python ~/src/python/access/salinity_update.py $ARCHIVEDIR/$RUN/history/ocn $RUN
    python ~/src/python/access/salinity_minmax.py $ARCHIVEDIR/$RUN/history/ocn $RUN
else
    python ~/src/python/access/ocean_means_update.py $ARCHIVEDIR/$RUN/history/ocn $RUN
    python ~/src/python/access/ocean_MOC_update.py $ARCHIVEDIR/$RUN/history/ocn $RUN
    python ~/src/python/access/ocean_ke_update.py $ARCHIVEDIR/$RUN/history/ocn $RUN
fi

python ~/src/python/access/ice_means_update.py $ARCHIVEDIR/$RUN/history/ice $RUN
