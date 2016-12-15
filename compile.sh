#!/bin/csh

setenv CPLUS_INCLUDE_PATH /afs/crc.nd.edu/group/cswarm/SuiteSparse/src/SuiteSparse-2.1.1/UMFPACK/Include:/afs/crc.nd.edu/group/cswarm/SuiteSparse/src/SuiteSparse-2.1.1/UFconfig:/afs/crc.nd.edu/group/cswarm/SuiteSparse/src/SuiteSparse-2.1.1/AMD/Include:/afs/crc.nd.edu/group/cswarm/hypre/2.4.0b/intel/15.0/mvapich2/2.1/include:/afs/crc.nd.edu/user/k/ksaha/Ivan/pgfem_3d/lib/VTK_IO/include



#/afs/crc.nd.edu/group/cswarm/Generalizsed_constitutive_model/utils/include:/afs/crc.nd.edu/group/cswarm/Generalizsed_constitutive_model/material/include:/afs/crc.nd.edu/group/cswarm/Generalizsed_constitutive_model/elasticity/include:/afs/crc.nd.edu/group/cswarm/Generalizsed_constitutive_model/crystal_plasticity/include

setenv LD_LIBRARY_PATH /afs/crc.nd.edu/group/cswarm/Generalizsed_constitutive_model/lib

module purge 
module load opt_local/1.0 
module load pgfem3d/mvapich2-2.1-intel-15.0-mlx
module load gcc/6.2.0    # to get support for c++14

autoreconf -if

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
#make clean
echo "++++++++++++ finished pulling pgfem_3d +++++++++++"
setenv PGFEM3D_INSTALL $PWD/build
./reconf_git_branch.sh
make -j 8
echo "++++++++++++ finished compiling pgfem_3d +++++++++++"
#make install 
echo "++++++++++++ finished building pgfem_3d +++++++++++"
