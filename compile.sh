#!/bin/csh

module purge 
module load opt_local/1.0 
module load pgfem3d/mvapich2-2.1-intel-15.0-mlx


#cd /afs/crc.nd.edu/group/cswarm/Generalizsed_constitutive_model
#git pull
#echo "++++++++++++ finished pulling Generalized_constitutive_model +++++++++++"
#make clean
#make
#echo "++++++++++++ finished compiling Generalized_constitutive_model +++++++++++"
#echo "    "

cp /afs/crc.nd.edu/user/i/iviti/pgconvert/pgfem_3d/build/convert2cc/share/config.site-intel /afs/crc.nd.edu/user/i/iviti/pgconvert/pgfem_3d/build/convert2cc/share/config.site

#cd /afs/crc.nd.edu/group/cswarm/pgfem_3d
#git pull
make distclean
echo "++++++++++++ finished pulling pgfem_3d +++++++++++"
setenv PGFEM3D_INSTALL /afs/crc.nd.edu/user/i/iviti/pgconvert/pgfem_3d/build
./reconf_git_branch.sh
#make clean
make 
echo "++++++++++++ finished compiling pgfem_3d +++++++++++"
#make install 
echo "++++++++++++ finished building pgfem_3d +++++++++++"
