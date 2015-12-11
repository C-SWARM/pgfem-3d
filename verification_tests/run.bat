#!/bin/bash

BRANCH_NAME=$(git branch | grep '*' | sed -e 's/* //')

NP=1
exe=$PGFEM3D_INSTALL/$BRANCH_NAME/bin/PGFem3D 

W_DIR=${TEST_DIR}/verification_tests/box_patch_LT
I_FORMAT=inputs/box_LT_${NP}CPU/box_LT_
O_FORMAT=out/box_LT_${NP}CPU/box

#mpirun -np ${NP} ${exe} -SS -disp -V ${W_DIR}/shear/${I_FORMAT} ${W_DIR}/shear/${O_FORMAT}
#mpirun -np ${NP} ${exe} -SS -disp -V ${W_DIR}/tension_disp/${I_FORMAT} ${W_DIR}/tension_disp/${O_FORMAT}
#mpirun -np ${NP} ${exe} -SS -disp -V ${W_DIR}/tension_pressure/${I_FORMAT} ${W_DIR}/tension_pressure/${O_FORMAT}

#mpirun -np ${NP} src/shear -SS -disp -V ${W_DIR}/shear/${I_FORMAT} ${W_DIR}/shear/${O_FORMAT}
#mpirun -np ${NP} src/tension_disp -SS -disp -V ${W_DIR}/tension_disp/${I_FORMAT} ${W_DIR}/tension_disp/${O_FORMAT}
#mpirun -np ${NP} src/tension_pressure -SS -disp -V ${W_DIR}/tension_pressure/${I_FORMAT} ${W_DIR}/tension_pressure/${O_FORMAT}

BOX_COMP=${TEST_DIR}/verification_tests/box_compression
mpirun -np ${NP} ${exe} -SS -cm 1 -V ${BOX_COMP}/input/box_LT_ ${BOX_COMP}/out/box_LT
#mpirun -np ${NP} src/crystal_plasticity -SS -cm 1 ${BOX_COMP}/input/box_LT_ ${BOX_COMP}/out/box_LT
