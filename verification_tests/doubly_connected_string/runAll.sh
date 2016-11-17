#!/bin/bash

#set -x


filebase=box
np=16

mv out.log out.old
g++ makeST.cpp -o makeST.o
./makeST.o > $filebase\_0.in.st
./local_makeset.pl clean
./local_makeset.pl -np $np
cp traction.in $filebase\_$np\CPU

g++ createInit.cpp -o createInit.o
cp createInit.o $filebase\_$np\CPU
cd $filebase\_$np\CPU

for i in `seq 0 $np`;
do

./createInit.o $filebase\_$i.in > $filebase\_$i.initial

done

for j in `seq 1 $np`;
do
perl -pi -e '$_ = "" if ( $. == 1 );' $filebase\_$j.initial

done


cd ..
qsub script.sh





