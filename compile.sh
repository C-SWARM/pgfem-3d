#!/bin/csh

module purge 
module load opt_local/1.0 
module load pgfem3d/mvapich2-2.1-intel-15.0-mlx
module load gcc/6.2.0    # to get support for c++14


#cd /afs/crc.nd.edu/group/cswarm/Generalizsed_constitutive_model
#git pull
#echo "++++++++++++ finished pulling Generalized_constitutive_model +++++++++++"
#make clean
#make
#echo "++++++++++++ finished compiling Generalized_constitutive_model +++++++++++"
#echo "    "

cp config.site-intel ./build/convert2cc/share/config.site

#cd /afs/crc.nd.edu/group/cswarm/pgfem_3d
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
