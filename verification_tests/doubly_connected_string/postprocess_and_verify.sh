#!/bin/bash

NP=4

cp postprocess/read_out.sh out/box_${NP}CPU/restart/
cp postprocess/bind_files.m out/box_${NP}CPU/restart/
cp postprocess/analyze_and_compare.m out/box_${NP}CPU/restart/
cp frequency.ref out/box_${NP}CPU/restart/

cd out/box_${NP}CPU/restart/
./read_out.sh
octave bind_files.m
octave analyze_and_compare.m 
