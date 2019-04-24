#!/bin/bash

if [ $# -lt 3 ]; then
  echo "usage: ./run.test.sh [NP] [restart #]"
  echo "       [NP]         : number of process"
  echo "       [restart #]  : restart number"
  echo "       [out dir}    : output directory name"
  exit 1
fi

NP=$1
restart=$2
out=$3

filebase=test
run="-SS -cm 1 -restart ${restart} -no-compute-macro -no-compute-reactions -V -kdim 3000 -maxit 3000  -isir -noCCE -noLS ${filebase}_${NP}CPU/${filebase}_ ${out}/${filebase}_${NP}CPU/${filebase}"

echo $run
