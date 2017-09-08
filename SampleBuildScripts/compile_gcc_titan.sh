#!/bin/csh

### Load the modules
module load autoconf/2.69
module unload PrgEnv-pgi/5.2.82
module load PrgEnv-gnu/5.2.82
module switch gcc/4.9.3 gcc/6.2.0
module load intel/16.0.3.210   # for mkl support
setenv CPLUS_INCLUDE_PATH /ccs/home/kv8/ttl/install_ttl_02032017/include     #Path the ttl header files


### User needs to set any dependent libraries' path and a directory where the pgfem3d to be built
### You may already have hypre, suitespare and vtk built in your systerm, if so just set the path in
setenv path_to_hypre /ccs/home/kv8/hypre/2.4.0b/gcc/6.2.0/cray-mpich/7.5.2
setenv path_to_suitesparse /ccs/home/kv8/SuiteSparse/src/SuiteSparse-2.1.1_gcc-metis_unsigned_long_int
setenv path_to_gcm /ccs/home/kv8/Generalizsed_constitutive_model_gcc
setenv path_to_pgfem3d /ccs/home/kv8/pgfem_3d
setenv path_to_pgfem3d_build /ccs/proj/csc188/nd/pgfem3d_build/gcc_6.2.0_cray-mpich_7.5.2
#setenv path_to_vtk /opt/crc/vtk/5.10.1/gcc


### This takes to the GCM directory, pull the latest code and compile it
cd $path_to_gcm    
git pull
make clean
make CXX=CC CXXFLAGS="-Wall -std=c++14 -fpermissive -O3" 
echo "++++++++++++ Finished Compiling GCM +++++++++++"


### This takes to the pgfem3d directory and pull the latest version
cd $path_to_pgfem3d
git pull

### Add the branch name to the pgfem3d build and compile
setenv branch_name `git branch | grep '*' | awk '{ print $2; }'`     #To get the name of the current branch
setenv PGFEM3D_INSTALL $path_to_pgfem3d_build/$branch_name           #Give the name of the build directory 

make distclean
autoreconf -if                                                       #To generate configure file


./configure --prefix=$PGFEM3D_INSTALL              \
--with-mpi=yes                                     \
CXX=CC                                             \
CXXFLAGS="-Wall -std=c++14 -fpermissive -O3"       \
--with-hypre-dir=$path_to_hypre                    \
--with-suitesparse-dir=$path_to_suitesparse        \
--with-cnstvm-dir=$path_to_gcm                     \
--enable-vtk=no 
#--with-vtk-include="-I$path_to_vtk/include/vtk-6.0"  \
#--with-vtk-libs="-Wl,-rpath=$path_to_vtk/lib -L$path_to_vtk/lib -lvtkIOXML-6.0 -lvtkIOXMLParser-6.0 -lvtkIOCore-6.0 -lvtkzlib-6.0 -lvtkCommonExecutionModel-6.0 -lvtkCommonDataModel-6.0 -lvtkCommonMisc-6.0 -lvtkCommonSystem-6.0 -lvtkCommonTransforms-6.0 -lvtkCommonMath-6.0 -lvtkIOGeometry-6.0 -lvtkCommonCore-6.0 -lvtksys-6.0 -lvtkexpat-6.0"

make -j 8
make install

echo "++++++++++++ Finished Compiling pgfem3d +++++++++++"
