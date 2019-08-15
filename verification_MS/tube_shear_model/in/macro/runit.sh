#!/bin/bash

nproc=$1
#BRANCH=master
BRANCH=develop
#pgfem_path=$PGFEM3D_INSTALL/$BRANCH/bin

#module load mvapich2/2.2-gcc-7.1.0-mlx
module load mvapich2/2.2-gcc-8.1.0-mlx
#pgfem_path=/scratch365/cswarm/skim/pgfem_3d/build_develop/bin
pgfem_path=/scratch365/cswarm/tphan2/pgfem_3d/src/

exe="$pgfem_path/PGFem3D -SS"


input=./macro_${nproc}CPU/macro_
output=./out/${nproc}CPU/macro
opts="-cm 1 -coh -V"

mpirun -np $nproc $exe $opts $input $output
#ddt -np $nproc $exe $opts $input $output
