#!/bin/csh

module purge 
module load opt_local/1.0 
module load mpich/3.1.2-gcc
module load hypre/2.4.0b-mpich-3.1.2-gcc-4.9.2
module load atlas/3.8.3
module load gcc
 

#cd /afs/crc.nd.edu/user/i/iviti/Generalizsed_constitutive_model_gcc 
#git pull
#echo "++++++++++++ finished pulling Generalized_constitutive_model +++++++++++"
#setenv MKLROOT /afs/crc.nd.edu/user/i/iviti/Generalizsed_constitutive_model_gcc/mkl_include
#make CC=mpicc CXX=mpicxx CFLAGS="-std=c99 -O3 -g" CXXFLAGS="-std=c99 -O3 -g"
#make clean
#make
#echo "++++++++++++ finished compiling Generalized_constitutive_model +++++++++++"
#echo "    "

cp ./build/convert2cc/share/config.site-gcc ./build/convert2cc/share/config.site

#cd /afs/crc.nd.edu/user/k/ksaha/Ivan/pgfem_3d 
#git pull
make distclean
echo "++++++++++++ finished pulling pgfem_3d +++++++++++"
setenv PGFEM3D_INSTALL $PWD/build
./reconf_git_branch.sh
#make clean
make 
echo "++++++++++++ finished compiling pgfem_3d +++++++++++"
#make install 
echo "++++++++++++ finished building pgfem_3d +++++++++++"
