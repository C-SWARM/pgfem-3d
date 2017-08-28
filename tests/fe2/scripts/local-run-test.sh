#!/bin/bash

branch=master
exe=$PGFEM3D_INSTALL/$branch/bin/PGFem3D

macro_nproc=1
micro_nproc=1
total_nproc=3

macro_in_dir=../in/macro/macro_${macro_nproc}CPU
micro_in_dir=../in/micro/pack_${micro_nproc}CPU
macro_opath=../cur/macro
micro_opath=../cur/micro

# macro-scale options
macro_opts="-disp -coh"
macro_opt_block="-macro-start $macro_opts -ipath $macro_in_dir -opath $macro_opath macro_ macro -macro-end"

# micro-scale options
micro_opts="-ms -disp"
micro_opt_block="-micro-start $micro_opts -ipath $micro_in_dir -opath $micro_opath pack_ pack -micro-end"

run="mpirun -np $total_nproc $exe -MS -macro-np $macro_nproc -micro-group-size $micro_nproc $macro_opt_block $micro_opt_block"

echo $run
$run |& tee test.log

./extract-force-data.sh test.log ../cur/force.dat

rm test.log