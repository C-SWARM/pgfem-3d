#!/bin/bash

if [ -r /opt/crc/Modules/current/init/bash ]; then
  source /opt/crc/Modules/current/init/bash
fi

module load mvapich2/2.2-gcc-7.1.0-mlx intel

exe=../../src/PGFem3D
map=../gen_restart_from_NP2NP

Froms=(4 12)
Tos=(12 4)

testno=${#Froms[@]}

for (( ia=0; ia<${testno}; ia++ )); do
  # run PGFem3D for 4 cpu case
  opt=`./gen.PGFem3D.options.sh ${Froms[$ia]} -1 out.ref`;
  mpirun -np ${Froms[$ia]} ${exe} ${opt}
done


# test mapping Froms -> Tos
for (( ia=0; ia<${testno}; ia++ )); do
  from=${Froms[$ia]}
  to=${Tos[$ia]}
  cp -r out.ref/test_${from}CPU out

  opt=`./gen.PGFem3D.options.sh ${from} 4 out`;

  ${map} map.info.${from}to${to} ${opt}

  opt=`./gen.PGFem3D.options.sh ${to} 4 out`;
  mpirun -np ${to} ${exe} ${opt}
  
  # verification
  for (( ib=0; ib<${to}; ib++ )); do
    ref=out.ref/test_${to}CPU/restart/Mechanical/STEP_000009/test_${ib}_9.res
    test=out/test_${to}CPU/restart/Mechanical/STEP_000009/test_${ib}_9.res
    numdiff -a 1.0e-6 -r 1.0e-6 ${ref} ${test}
  done
  rm -rf out
done

rm -rf out.ref
