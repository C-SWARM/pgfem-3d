#!/bin/bash

NP=1
echo ${HOME}
exe=PGFem3D

W_DIR=${TEST_DIR}/verification_tests/box_patch_LT
I_FORMAT=inputs/box_LT_${NP}CPU/box_LT_
O_FORMAT=out/box_LT_${NP}CPU/box

mpirun -np ${NP} ${exe} -SS -disp -V ${W_DIR}/shear/${I_FORMAT} ${W_DIR}/shear/${O_FORMAT}
mpirun -np ${NP} ${exe} -SS -disp -V ${W_DIR}/tension_disp/${I_FORMAT} ${W_DIR}/tension_disp/${O_FORMAT}
mpirun -np ${NP} ${exe} -SS -disp -V ${W_DIR}/tension_pressure/${I_FORMAT} ${W_DIR}/tension_pressure/${O_FORMAT}

mpirun -np ${NP} src/shear -SS -disp -V ${W_DIR}/shear/${I_FORMAT} ${W_DIR}/shear/${O_FORMAT}
mpirun -np ${NP} src/tension_disp -SS -disp -V ${W_DIR}/tension_disp/${I_FORMAT} ${W_DIR}/tension_disp/${O_FORMAT}
mpirun -np ${NP} src/tension_pressure -SS -disp -V ${W_DIR}/tension_pressure/${I_FORMAT} ${W_DIR}/tension_pressure/${O_FORMAT}

