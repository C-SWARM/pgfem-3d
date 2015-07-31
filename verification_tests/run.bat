#!/bin/bash

NP=1
TEST_EXAMPLE_SHEAR=${TEST_DIR}/verification_tests/box_patch_LT/shear
TEST_EXAMPLE_TDISP=${TEST_DIR}/verification_tests/box_patch_LT/tension_disp
TEST_EXAMPLE_TPRES=${TEST_DIR}/verification_tests/box_patch_LT/tension_pressure

mpirun -np ${NP} PGFem3D -SS -disp -V ${TEST_EXAMPLE_SHEAR}/inputs/box_LT_${NP}CPU/box_LT_ ${TEST_EXAMPLE_SHEAR}/out/box_LT_${NP}CPU/box
mpirun -np ${NP} PGFem3D -SS -disp -V ${TEST_EXAMPLE_TDISP}/inputs/box_LT_${NP}CPU/box_LT_ ${TEST_EXAMPLE_TDISP}/out/box_LT_${NP}CPU/box
mpirun -np ${NP} PGFem3D -SS -disp -V ${TEST_EXAMPLE_TPRES}/inputs/box_LT_${NP}CPU/box_LT_ ${TEST_EXAMPLE_TPRES}/out/box_LT_${NP}CPU/box

#mpirun -np ${NP} src/shear -SS -disp -V ${TEST_EXAMPLE_SHEAR}/inputs/box_LT_${NP}CPU/box_LT_ ${TEST_EXAMPLE_SHEAR}/out/box_LT_${NP}CPU/box
#mpirun -np ${NP} src/tension_disp -SS -disp -V ${TEST_EXAMPLE_TDISP}/inputs/box_LT_${NP}CPU/box_LT_ ${TEST_EXAMPLE_TDISP}/out/box_LT_${NP}CPU/box
#mpirun -np ${NP} src/tension_pressure -SS -disp -V ${TEST_EXAMPLE_TPRES}/inputs/box_LT_${NP}CPU/box_LT_ ${TEST_EXAMPLE_TPRES}/out/box_LT_${NP}CPU/box

