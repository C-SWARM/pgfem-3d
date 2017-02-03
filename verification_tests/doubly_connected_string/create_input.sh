#!/bin/bash

####set -x

filebase=box
NP=4

icpc src/create_dot_st.cc -o create_dot_st.o
./create_dot_st.o > $filebase\_0.in.st
./local_makeset.pl clean
./local_makeset.pl -np ${NP}
\rm -r partitions.${NP}
cp traction.in $filebase\_${NP}\CPU

icpc src/create_dot_initial.cc -o create_dot_initial.o
cp create_dot_initial.o $filebase\_${NP}\CPU
cd $filebase\_${NP}\CPU

for i in `seq 0 ${NP}`;
do
./create_dot_initial.o $filebase\_$i.in > $filebase\_$i.initial
done

for j in `seq 1 ${NP}`;
do
perl -pi -e '$_ = "" if ( $. == 1 );' $filebase\_$j.initial
done


