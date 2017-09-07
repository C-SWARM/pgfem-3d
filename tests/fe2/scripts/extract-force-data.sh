#!/bin/bash

# This shell script extracts the force data from the console output of
# PGFem3D. The console ouput must of course be piped into a file for
# processing, i.e. `PGFem3D ${opts} > ${logfile}`. Only the forces
# associated with the first surface specified in the `traction.in`
# file are extracted from the output.
#
# extract-force-data.sh ${logfile} ${outfile}
#
# Return non-zero status on error

if [ $# -lt 2 ]; then
    echo "ERROR: must specify input and output filenames!"
    exit 1;
fi

scriptname=$0
logfile=$1
outfile=$2


# extract the data from the logfile
grep -A 1 '^Forces on' $logfile | sed -e '/^Forces/d' -e '/^--/d' > $outfile



