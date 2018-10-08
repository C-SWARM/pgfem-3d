#!/bin/bash

# number of precesses is fixed to 8 for this example
NP=8
filebase=box

if [ $# -lt 1 ]; then
  echo "usage: ./run_test.sh [test_name] [gen_mesh] [transport_type]"
  echo "       [test_name]      : test_name = tension or compression"
  echo "       [gen_mesh]       : if gen_mesh = 1: cleanup everything and re-generate mesh and run the test"
  echo "                           otherwise: run the test"
  echo "       [transport_type] : mpi or photon communication backend"
  exit 1
fi

##########################################
# read options
##########################################
if [ $# -lt 2 ]; then
  test_name=$1
  gen_mesh=1
else
  test_name=$1
  gen_mesh=$2
  transport_type=$3
fi

##########################################
# if gen_mesh = 1, generate new inputs
##########################################
if [ $gen_mesh = 1 ]; then
  rm -rf ${filebase}_${NP}CPU_${test_name}
  cp ${filebase}.out.header_${test_name} ${filebase}.out.header
  ../local_makeset.pl -np $NP -f ${filebase} clean -d 0.36
  mv ${filebase}_${NP}CPU ${filebase}_${NP}CPU_${test_name}
fi

##########################################
# run simulations
##########################################

#branch=master

TEST_DIR=$PWD
exe=../../build/src/PGFem3D

# inputs
IN_DIR=${TEST_DIR}/${filebase}_${NP}CPU_${test_name}
input=${IN_DIR}/${filebase}_

# outputs
OUTDIR=${test_name}/${filebase}_${NP}CPU
output=${OUTDIR}/${filebase}

override="-override-solver-file ${filebase}_0.in.st_${test_name}"
opts="-SS ${transport_type} -cm 1 -noLS -noCCE -no-compute-macro -maxit 3000 -kdim 1000 -V ${override}"
mpirun -np $NP $exe $opts $input $output
