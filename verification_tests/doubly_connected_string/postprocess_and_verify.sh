#!/bin/bash

NP=4

data_dir=out/box_${NP}CPU/restart/Mechanical

cp postprocess/read_out.sh $data_dir/.
cp postprocess/bind_files.m $data_dir/. 
cp postprocess/analyze_and_compare.m $data_dir/.
cp frequency.ref $data_dir/.

cd $data_dir 
./read_out.sh
octave bind_files.m
octave analyze_and_compare.m 
