#!/bin/bash
#$ -pe mpi-16 16
#$ -q *@@cswarm
#$ -o out-server.log

module load mvapich2/2.3.1/gcc/8.3.0 

exe=/afs/crc.nd.edu/user/k/ksaha/pgfem3d/pgfem_3d/src/PGFem3D

macro_nproc=3
micro_nproc=3
total_nproc=6

macro_in_dir=../in/macro/macro_${macro_nproc}CPU
micro_in_dir=../in/micro/pack_${micro_nproc}CPU
macro_opath=../cur/macro
micro_opath=../cur/micro

# macro-scale options
# macro_opts="-disp -coh -V -isir"
macro_opts="-disp -coh -isir"
macro_opt_block="-macro-start $macro_opts -ipath $macro_in_dir -opath $macro_opath macro_ macro -macro-end"

# micro-scale options
# micro_opts="-ms -disp -V -isir"
micro_opts="-ms -disp -isir"
micro_opt_block="-micro-start $micro_opts -ipath $micro_in_dir -opath $micro_opath pack_ pack -micro-end"

run="mpirun -np $total_nproc $exe -MS -macro-np $macro_nproc -micro-group-size $micro_nproc $macro_opt_block $micro_opt_block"

echo $run
$run |& tee out-server.log

./extract-force-data.sh out-server.log ./force.dat
