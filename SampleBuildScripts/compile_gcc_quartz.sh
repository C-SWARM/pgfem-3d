#!/bin/csh

### Load the modules
module switch intel/16.0.3 gcc/6.1.0
module load mvapich2/2.2
module load mkl/11.3.3
setenv CPLUS_INCLUDE_PATH /usr/gapps/notredame/external_lib_quartz/ttl/install_ttl_02162017/include  #Path the ttl header files


### User needs to set any dependent libraries' path and a directory where the pgfem3d to be built
### You may already have hypre, suitespare and vtk built in your systerm, if so just set the path in
setenv path_to_hypre /usr/gapps/notredame/external_lib_quartz/hypre/2.4.0b/gcc/6.1.0/mvapich2/2.2
setenv path_to_suitesparse /usr/gapps/notredame/external_lib_quartz/suitesparse/src/SuiteSparse-2.1.1_gcc6.1.0
setenv path_to_gcm /g/g90/saha4/Generalizsed_constitutive_model_cab
setenv path_to_pgfem3d /g/g90/saha4/pgfem_3d
setenv path_to_pgfem3d_build /usr/gapps/notredame/pgfem_3d_install_quartz_gcc
#setenv path_to_vtk /opt/crc/vtk/5.10.1/gcc


### This takes to the GCM directory, pull the latest code and compile it
cd $path_to_gcm    
git pull
make clean
make CXX=mpicxx CXXFLAGS="-Wall -std=c++14 -Ofast -fpermissive -march=core-avx2" 
echo "++++++++++++ Finished Compiling GCM +++++++++++"


### This takes to the pgfem3d directory and pull the latest version
cd $path_to_pgfem3d
git pull

### Add the branch name to the pgfem3d build and compile
setenv branch_name `git branch | grep '*' | awk '{ print $2; }'`     #To get the name of the current branch
setenv PGFEM3D_INSTALL $path_to_pgfem3d_build/$branch_name           #Give the name of the build directory 

make distclean
autoreconf -if                                                       #To generate configure file


./configure --prefix=$PGFEM3D_INSTALL                               \
--with-mpi=yes                                                      \
CXX=mpicxx                                                          \
CXXFLAGS="-Wall -std=c++14 -Ofast -fpermissive -march=core-avx2"    \
--with-hypre-dir=$path_to_hypre                                     \
--with-suitesparse-dir=$path_to_suitesparse                         \
--with-cnstvm-dir=$path_to_gcm                                      \
--enable-vtk=no 
#--with-vtk-include="-I$path_to_vtk/include/vtk-6.0"  \
#--with-vtk-libs="-Wl,-rpath=$path_to_vtk/lib -L$path_to_vtk/lib -lvtkIOXML-6.0 -lvtkIOXMLParser-6.0 -lvtkIOCore-6.0 -lvtkzlib-6.0 -lvtkCommonExecutionModel-6.0 -lvtkCommonDataModel-6.0 -lvtkCommonMisc-6.0 -lvtkCommonSystem-6.0 -lvtkCommonTransforms-6.0 -lvtkCommonMath-6.0 -lvtkIOGeometry-6.0 -lvtkCommonCore-6.0 -lvtksys-6.0 -lvtkexpat-6.0"

make -j 8
make install

echo "++++++++++++ Finished Compiling pgfem3d +++++++++++"
