RUN=$1

OUTPUT_DIR=/g/data/p66/mrd599/access_stats/${RUN}

mkdir -p $OUTPUT_DIR
# On gadi, ARCHIVEDIR may be either in /g/data or /scratch
ARCHIVEDIR=/g/data/$PROJECT/$USER/archive
if [ ! -d $ARCHIVEDIR/$RUN ]; then
    ARCHIVEDIR=/scratch/$PROJECT/$USER/archive
fi

cd ${ARCHIVEDIR}/${RUN}/history/atm

export PYTHONPATH=$PYTHONPATH:~/src/python/umfile

# If there are file names of form ${RUN}a.pm???????.nc use new style
if compgen -G "${RUN}a.pm???????.nc" > /dev/null; then

    python ~/src/python/access/access_means_update_new.py -r ${RUN} -v 23 24 31 32 507 1201 1207 1208 1209 1210 1211 1235 2201 2204 2205 2206 2207 2208 3217 3234 3235 3236 3332 3258 5216 8223 8225 3223 3232 3353 8234 8235 8245 26004 3202 3314 3298  -o $OUTPUT_DIR -n
    # Don't fail on missing data
    if [[ $? != 0 ]]; then
        echo "Problem from access_means_update", ${RUN}
    fi

else

    find . -name "${RUN}a.pm???????" | parallel "[[ -f $OUTPUT_DIR/{}.nc ]] || echo {}"
    find . -name "${RUN}a.pm???????" | parallel "[[ -f $OUTPUT_DIR/{}.nc ]] || ( python ~mrd599/src/python/umfile/um_fields_subset.py -v  23,24,31,32,507,1201,1207,1208,1209,1210,1211,1235,2201,2204,2205,2206,2207,2208,3217,3234,3235,3236,3332,3258,5216,8223,8225,3223,3232,3353,8234,8235,8245,26004,3202,3314,3298 -i {} -o $PBS_JOBFS/{} &&  python ~mrd599/src/python/umfile/um2netcdf_iris_mon.py $PBS_JOBFS/{} $OUTPUT_DIR/{}.nc && rm $PBS_JOBFS/{} ) "

    cd $OUTPUT_DIR

    python ~mrd599/src/python/access/access_means_update.py

fi
