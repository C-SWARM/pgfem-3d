#!/bin/csh

### Load the modules
module purge 
module load opt_local/1.0 ttl mvapich2/2.2-gcc-6.2.0-mlx 
module load intel/16.0 # for mkl support


### User needs to set any dependent libraries' path and a directory where the pgfem3d to be built
### You may already have hypre, suitespare and vtk built in your systerm, if so just set the path in
setenv path_to_hypre /afs/crc.nd.edu/group/cswarm/hypre/2.4.0b/gcc/6.2.0/mvapich2/2.2                 
setenv path_to_suitesparse /afs/crc.nd.edu/group/cswarm/SuiteSparse/src/SuiteSparse-2.1.1_gcc6.2.0
setenv path_to_gcm /afs/crc.nd.edu/group/cswarm/Generalizsed_constitutive_model                       
setenv path_to_pgfem3d /afs/crc.nd.edu/group/cswarm/pgfem_3d
setenv path_to_pgfem3d_build /afs/crc.nd.edu/group/cswarm/pgfem_3d/build_pgfem3d_mvapich2-2.2-gcc-6.2.0-mlx
#setenv path_to_vtk /opt/crc/vtk/5.10.1/gcc


### This takes to the GCM directory, pull the latest code and compile it
cd $path_to_gcm    
git pull
make clean
make CC=mpicxx CXXFLAGS="-Wall -std=c++14 -fpermissive -O3"
echo "++++++++++++ Finished Compiling GCM +++++++++++"


### This takes to pgfem3d directory and pullt the latest code 
cd $path_to_pgfem3d
git pull

### Add the branch name to the pgfem3d build and compile
setenv branch_name `git branch | grep '*' | awk '{ print $2; }'`     #To get the name of the current branch
setenv PGFEM3D_INSTALL $path_to_pgfem3d_build/$branch_name           #Give the name of the build directory 

make distclean
autoreconf -if                                                       #To generate configure file


./configure --prefix=$PGFEM3D_INSTALL              \
--with-mpi=yes                                     \
CXX=mpicxx                                         \
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
