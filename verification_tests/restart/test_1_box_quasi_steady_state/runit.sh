#!/bin/bash

NP=4

branch=adaptive_time_stepping
name=box
TEST_DIR=$PWD
IN_DIR=${TEST_DIR}/${name}_${NP}CPU
exe=$PGFEM3D_INSTALL/${branch}/bin/PGFem3D
#override="-override-pre-disp initial_disp.in"

input=${IN_DIR}/${name}_
OUTDIR=${TEST_DIR}/out/${name}_${NP}CPU
output=${OUTDIR}/${name}

# get restart
RN=-1
opts="-SS -cm 1 -restart $RN -noLS -no-compute-reactions -no-compute-macro -kdim 1000 -maxit 2000 -V ${override}"
mpirun -np $NP $exe $opts $input $output >& out_run_all.log

RN=3
opts="-SS -cm 1 -restart $RN -noLS -no-compute-reactions -no-compute-macro -kdim 1000 -maxit 2000 -V ${override}"
mpirun -np $NP $exe $opts $input $output >& out_restart.log

grep -m 5 "Damage thresh alpha:" out_run_all.log | tail -n 1 >& Damage_all.txt
grep -m 1 "Damage thresh alpha:" out_restart.log | tail -n 1 >& Damage_restart.txt

/cswarm/tools/bin/numdiff Damage_all.txt Damage_restart.txt
rm -rf CRYSTAL_ORIENTATION out *.log *.txt
