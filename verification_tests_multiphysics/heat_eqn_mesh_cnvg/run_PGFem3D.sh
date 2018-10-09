#!/bin/bash

# read working directory path
output_dir=out

# define constant for running PGFem3D

if [ $# -lt 4 ]; then
  echo "usage: ./run_PGFem3D.sh [file_base] [NP] [task] [network]"
  echo "       [file_base] : filebase name"
  echo "       [NP]        : number of process"
  echo "       [task]      : task name"
  echo "       [network] : -isir or -pwc"
  exit 1
fi

# read options
filebase=$1
NP=$2
task=$3
network=$4

# run simulations
#branch=master

TEST_DIR=$PWD
IN_DIR=${TEST_DIR}/${filebase}_${NP}CPU_${task}
#exe=$PGFEM3D_INSTALL/${branch}/bin/PGFem3D

input=${IN_DIR}/${filebase}_
OUTDIR=$output_dir/mesh_cnvg/${task}/${filebase}_${NP}CPU
output=${OUTDIR}/${filebase}

opts="-SS -cm 1 -noLS -no-compute-macro -noCCE -maxit 3000 -kdim 1000 -V ${override} ${network}"

time_start=`date`
echo "============================================================="
echo "start run: mpirun -np ${NP} ${exe} ${opts} ${input} ${output}"
echo "============================================================="
echo "Start time: ${time_start}"
echo "============================================================="
echo "job ID: $JOB_ID"
echo "host name: $HOSTNAME"
echo "============================================================="

mpirun -np $NP $exe $opts $input $output
time_end='date'

echo "============================================================="
echo "End time: ${time_end}"
echo "============================================================="
