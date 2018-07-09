#!/bin/csh
###$ -N hyp 
#$ -pe mpi-16 256 
#$ -q *@@cswarm
#$ -o out.log 


# define constant for reading options

if ($#argv < 4) then
  echo "usage: ./run.sh [NP] [el_size] [el_order] [method|#]"
  echo "       [NP]     : number of process"
  echo "       [el_size]  : element size";
  echo "       [el_order] : element order linear=1, quadratic=2";
  echo "       [method|#] : disp, cm 1, cm3f 1"
  exit 0;
endif

setenv NP $1
setenv el_size $2
setenv el_order $3
setenv method $4
setenv numb $5

module load CRC_default/1.0
module load pgfem3d/mvapich2-2.2-gcc-7.1.0-mlx-develop

setenv TEST_DIR $PWD
setenv name beam
setenv exe /afs/crc.nd.edu/user/s/skim43/PGFem3D_tutorial/pgfem_3d/build_develop/bin/PGFem3D
setenv opts "-SS -${method} ${numb} -restart -1 -noLS -no-compute-macro -maxit 5000 -kdim 3000 -V"

setenv input ${method}_${name}_${NP}CPU/${name}_
setenv output ${TEST_DIR}/out_${method}S${el_size}O${el_order}/${name}_

setenv run "mpirun -np $NP $exe  $opts $input $output"

echo $run
$run



