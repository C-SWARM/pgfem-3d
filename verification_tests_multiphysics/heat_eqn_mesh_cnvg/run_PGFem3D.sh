#!/bin/bash

# read working directory path
work_dir=out

# define constant for running PGFem3D

if [ $# -lt 3 ]; then
  echo "usage: ./run_PGFem3D.sh [file_base] [NP] [task]"
  echo "       [file_base] : filebase name"
  echo "       [NP]        : number of process"
  echo "       [task]      : task name"
  exit 1
fi

# read options
filebase=$1
NP=$2
task=$3

# run simulations
branch=master

name=${filebase}
TEST_DIR=$PWD
IN_DIR=${TEST_DIR}/${name}_${NP}CPU_${task}
exe=$PGFEM3D_INSTALL/${branch}/bin/PGFem3D

input=${IN_DIR}/${name}_
OUTDIR=$work_dir/mesh_cnvg/${task}/${name}_${NP}CPU
output=${OUTDIR}/${name}

opts="-SS -cm 1 -noLS -no-compute-macro -noCCE -maxit 3000 -kdim 1000 -V ${override}"

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

