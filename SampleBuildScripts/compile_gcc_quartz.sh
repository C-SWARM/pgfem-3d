#!/bin/csh

### Load the modules
setenv CPLUS_INCLUDE_PATH /usr/gapps/notredame/external_lib_quartz/ttl/install_ttl_02162017/include
module switch intel/16.0.3 gcc/6.1.0 
module load mvapich2/2.2   
module load mkl/11.3.3


### Build GCM library
cd /g/g90/saha4/Generalizsed_constitutive_model_cab                  #PATH_to_GCM_direcotry
git pull
make clean
make CXX=mpicxx CXXFLAGS="-std=c++14 -Ofast -fpermissive -g -march=core-avx2"
echo "... Finished compiling GCM ... " 
echo ""


### Compile pgfem_3d 
cd /g/g90/saha4/pgfem_3d                                             #PATH_to_pgfem_3d_directory
git pull

setenv PGFEM3D_INSTALL $PWD/build_pgfem3d                            #Give the name of the build directory 
setenv confg_file_name SampleBuildScripts/config.site_gcc_quartz     #Give the name of the config.site file

setenv branch_name `git branch | grep '*' | awk '{ print $2; }'`     #To get the name of the current branch
mkdir -pv $PGFEM3D_INSTALL/$branch_name/share
cp $confg_file_name $PGFEM3D_INSTALL/$branch_name/share/config.site

#make clean
make distclean
sh reconf_git_branch.sh
make -j 12
make install
 
echo "... Finished compiling pgfem_3d ... " 
