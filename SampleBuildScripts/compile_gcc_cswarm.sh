#!/bin/csh

### Load the modules
module purge 
module load opt_local/1.0 ttl mpich/3.2-gcc-6.2.0 
module load intel/16.0 # for mkl support


### Compile pgfem_3d 
git pull


setenv branch_name `git branch | grep '*' | awk '{ print $2; }'`     #To get the name of the current branch
setenv PGFEM3D_INSTALL $PWD/build_pgfem3d/$branch_name               #Give the name of the build directory 

make distclean
autoreconf -if                                                       #To generate configure file

setenv vtk_dir /opt/crc/vtk/5.10.1/gcc

./configure --prefix=$PGFEM3D_INSTALL \
--with-mpi=yes \
CXX=mpicxx     \
CXXFLAGS="-Wall -std=c++14 -fpermissive -Ofast -g" \
--with-hypre-dir=/afs/crc.nd.edu/group/cswarm/hypre/2.4.0b/gcc/6.2.0/mpich/3.2 \
--with-suitesparse-dir=/afs/crc.nd.edu/group/cswarm/SuiteSparse/src/SuiteSparse-2.1.1_gcc6.2.0 \
--with-cnstvm-dir=/cswarm/Jenkins_lib/Generalizsed_constitutive_model  \
--enable-vtk=yes \
--with-vtk-include="-I$vtk_dir/include/vtk-6.0"  \
--with-vtk-libs="-Wl,-rpath=$vtk_dir/lib -L$vtk_dir/lib -lvtkIOXML-6.0 -lvtkIOXMLParser-6.0 -lvtkIOCore-6.0 -lvtkzlib-6.0 -lvtkCommonExecutionModel-6.0 -lvtkCommonDataModel-6.0 -lvtkCommonMisc-6.0 -lvtkCommonSystem-6.0 -lvtkCommonTransforms-6.0 -lvtkCommonMath-6.0 -lvtkIOGeometry-6.0 -lvtkCommonCore-6.0 -lvtksys-6.0 -lvtkexpat-6.0"

make -j 8
make install
