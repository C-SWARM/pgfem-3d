#!/bin/bash
###$ -N hyp 
#$ -pe mpi-16 16
#$ -q *@@cswarm
#$ -o out.log


module load mvapich2/2.1-intel-15.0-mlx

NP=16
echo $NP


TEST_DIR=$PWD
name=box
#opts="-SS -disp -V -kdim 1000 -maxit 1000"
#opts="-SS -disp -V -noCCE -noLS -no-compute-macro"
opts="-SS -disp -V "
input=${name}_${NP}CPU/${name}_
output=${TEST_DIR}/out/${name}_${NP}CPU/${name}

run="mpirun -np $NP /afs/crc.nd.edu/user/i/iviti/pgBuild/master/bin/PGFem3D $opts $input $output"

echo $run
$run
