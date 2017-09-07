#!/bin/bash

nproc=$1
BRANCH=master
pgfem_path=$PGFEM3D_INSTALL/$BRANCH/bin

exe="$pgfem_path/PGFem3D -SS"
input=./pack_${nproc}CPU/pack_
output=./out/${nproc}CPU/pack
opts='-ms -disp -V -maxit 3000 -kdim 1000 -override-pre-disp pack.load'

mpirun -np $nproc $exe $opts $input $output
#ddt -np $nproc $exe $opts $input $output
