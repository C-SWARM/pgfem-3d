#!/bin/bash

TEST_DIR=$PWD
echo ${TEST_DIR}

mpirun -np 1 PGFem3D -SS -disp ${TEST_DIR}/inputs/box_LT_ ${TEST_DIR}/out/box_LT
