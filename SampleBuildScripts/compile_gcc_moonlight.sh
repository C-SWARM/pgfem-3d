#!/bin/csh

### Load the modules
module purge
#module load mvapich2/2.2   #compilation failed
module load openmpi/1.10.5 
module load mkl/11.3.3
module load gcc/6.2.0
setenv CPLUS_INCLUDE_PATH /turquoise/usr/projects/notredame/external_lib/install_ttl_02032017/include


### User needs to set any dependent libraries' path and a directory where the pgfem3d to be built
### You may already have hypre, suitespare and vtk built in your systerm, if so just set the path in
setenv path_to_hypre /turquoise/usr/projects/notredame/external_lib/hypre/2.4.0b/gcc/6.2.0/openmpi/1.10.5
setenv path_to_suitesparse /turquoise/usr/projects/notredame/external_lib/SuiteSparse/src/SuiteSparse-2.1.1_gcc6.2.0
setenv path_to_gcm /turquoise/usr/projects/notredame/Generalizsed_constitutive_model 
setenv path_to_pgfem3d /turquoise/usr/projects/notredame/pgfem_3d 
setenv path_to_pgfem3d_build /turquoise/usr/projects/notredame/pgfem_3d_install_moonlight_gcc 
#setenv path_to_vtk /opt/crc/vtk/5.10.1/gcc


### This takes to the GCM directory, pull the latest code and compile it
cd $path_to_gcm
#git pull
make clean
make CXX=mpicxx CXXFLAGS="-Wall -std=c++14 -fpermissive -Ofast -g"                  
echo "++++++++++++ Finished Compiling GCM +++++++++++"


### This takes to the pgfem3d directory and pull the latest version
cd $path_to_pgfem3d
#git pull

### Add the branch name to the pgfem3d build and compile
setenv branch_name `git branch | grep '*' | awk '{ print $2; }'`     #To get the name of the current branch
setenv PGFEM3D_INSTALL $path_to_pgfem3d_build/$branch_name           #Give the name of the build directory 

make distclean
autoreconf -if                                                       #To generate configure file


./configure --prefix=$PGFEM3D_INSTALL                               \
--with-mpi=yes                                                      \
CXX=mpicxx                                                          \
CXXFLAGS="-Wall -std=c++14 -fpermissive -Ofast -g"                  \
--with-hypre-dir=$path_to_hypre                                     \
--with-suitesparse-dir=$path_to_suitesparse                         \
--with-cnstvm-dir=$path_to_gcm                                      \
--enable-vtk=no
#--with-vtk-include="-I$path_to_vtk/include/vtk-6.0"  \
#--with-vtk-libs="-Wl,-rpath=$path_to_vtk/lib -L$path_to_vtk/lib -lvtkIOXML-6.0 -lvtkIOXMLParser-6.0 -lvtkIOCore-6.0 -lvtkzlib-6.0 -lvtkCommonExecutionModel-6.0 -lvtkCommonDataModel-6.0 -lvtkCommonMisc-6.0 -lvtkCommonSystem-6.0 -lvtkCommonTransforms-6.0 -lvtkCommonMath-6.0 -lvtkIOGeometry-6.0 -lvtkCommonCore-6.0 -lvtksys-6.0 -lvtkexpat-6.0"

make -j 8
make install

echo "++++++++++++ Finished Compiling pgfem3d +++++++++++"
