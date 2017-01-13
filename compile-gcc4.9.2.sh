#!/bin/csh

module purge 
module load opt_local/1.0 
module load mpich/3.1.2-gcc
module load gcc/4.9.2
 

cd /afs/crc.nd.edu/user/k/ksaha/Ivan/Generalizsed_constitutive_model_gcc4.9.2
#git pull
setenv MKLROOT /afs/crc.nd.edu/user/k/ksaha/Ivan/Generalizsed_constitutive_model_gcc4.9.2/mkl_include 
make clean
make CC=mpicc CXX=mpicxx CFLAGS="-std=c99 -O3 -g" CXXFLAGS="-std=c99 -O3 -g"
#echo "++++++++++++ finished compiling Generalized_constitutive_model +++++++++++"
#echo "    "


cd /afs/crc.nd.edu/user/k/ksaha/Ivan/pgfem_3d 
cp config.site-gcc4.9.2 ./build/convert2cc/share/config.site
#git pull
make distclean
#make clean
echo "++++++++++++ finished pulling pgfem_3d +++++++++++"
setenv PGFEM3D_INSTALL $PWD/build
./reconf_git_branch.sh
make 
echo "++++++++++++ finished compiling pgfem_3d +++++++++++"
#make install 
echo "++++++++++++ FInished building pgfem_3d +++++++++++"
