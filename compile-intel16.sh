#!/bin/csh

module purge 
module load opt_local/1.0 
module load mvapich2/2.2-intel-16.0-mlx
module load gcc/6.2.0    # to get support for c++14


cd /afs/crc.nd.edu/user/k/ksaha/Ivan/Generalizsed_constitutive_model 
#git pull
make clean
make
#echo "++++++++++++ finished compiling Generalized_constitutive_model +++++++++++"
#echo "    "


cd /afs/crc.nd.edu/user/k/ksaha/Ivan/pgfem_3d 
cp config.site-intel16 ./build/convert2cc/share/config.site
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
