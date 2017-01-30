#!/bin/csh

module purge 
module load opt_local/1.0 
module load mvapich2/2.2-gcc-6.2.0-mlx 
module load intel/16.0


cd /afs/crc.nd.edu/user/k/ksaha/Ivan/Generalizsed_constitutive_model
#git pull
#echo "++++++++++++ finished pulling Generalized_constitutive_model +++++++++++"
####setenv MKLROOT /afs/crc.nd.edu/user/i/iviti/Generalizsed_constitutive_model/mkl_include
#make clean
#make
#echo "++++++++++++ finished compiling Generalized_constitutive_model +++++++++++"
#echo "    "


cd /afs/crc.nd.edu/user/k/ksaha/Ivan/pgfem_3d 
cp config.site-gcc6.2.0 ./build/convert2cc/share/config.site
#git pull
make distclean
#make clean
echo "++++++++++++ finished pulling pgfem_3d +++++++++++"
setenv PGFEM3D_INSTALL $PWD/build
./reconf_git_branch.sh
make -j 8
echo "++++++++++++ finished compiling pgfem_3d +++++++++++"
make install 
echo "++++++++++++ finished building pgfem_3d +++++++++++"
