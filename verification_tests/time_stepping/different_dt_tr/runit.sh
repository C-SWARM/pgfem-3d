#!/bin/bash

NP=4

branch=adaptive_time_stepping
name=box
TEST_DIR=$PWD
IN_DIR=${TEST_DIR}/../inputs/${name}_${NP}CPU_tr
exe=$PGFEM3D_INSTALL/${branch}/bin/PGFem3D
#override_solver="-override-solver-file box_0.in.st_3steps"
override_disp="-override-pre-disp initial_disp.in"
override="$override_solver $override_disp"

input=${IN_DIR}/${name}_
opts="-SS -cm 1 -noLS -no-compute-reactions -no-compute-macro -kdim 1000 -maxit 2000 -V ${override}"

OUTDIR_3steps=${TEST_DIR}/out_3steps/${name}_${NP}CPU
output=${OUTDIR_3steps}/${name}
mpirun -np $NP $exe $opts -override-solver-file box_0.in.st_3steps $input $output >& out_run_3steps.log

OUTDIR_5steps=${TEST_DIR}/out_5steps/${name}_${NP}CPU
output=${OUTDIR_5steps}/${name}
mpirun -np $NP $exe $opts -override-solver-file box_0.in.st_5steps $input $output >& out_run_5steps.log


# NP should not be greater than 10

for ((n=0; n<$NP; n++)); do
  node=`grep -m 1 -o "[0-9]*" ${input}$n.in | head -n 1`
  head ${OUTDIR_3steps}/restart/STEP_00002/${name}_${n}_2.res -n $node >& disp_3steps_${n}.txt
  head ${OUTDIR_5steps}/restart/STEP_00004/${name}_${n}_4.res -n $node >& disp_5steps_${n}.txt
  /cswarm/tools/bin/numdiff -a 1.0e-5 disp_3steps_${n}.txt disp_5steps_${n}.txt
done

rm -rf *.txt *.log
rm -rf out_3steps out_5steps CRYSTAL_ORIENTATION
