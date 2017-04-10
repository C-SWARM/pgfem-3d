#!/bin/csh

### Load the modules
module purge 
module load opt_local/1.0 mvapich2/2.2-gcc-6.2.0-mlx 
module load intel/16.0 # for mkl support
module load ttl
#module load pmtl/4.0
#module load boost/1.58-mvapich2-2.2-gcc-6.2.0-mlx


### Build GCM library
cd /afs/crc.nd.edu/group/cswarm/Generalizsed_constitutive_model      #PATH_to_GCM_direcotry
git pull
make clean
make
echo "... Finished compiling GCM ... " 
echo ""


### Compile pgfem_3d 
cd /scratch365/cswarm/ksaha/pgfem_3d                                 #PATH_to_pgfem_3d_directory
git pull

setenv PGFEM3D_INSTALL $PWD/build_pgfem3d                            #Give the name of the build directory 
setenv confg_file_name SampleBuildScripts/config.site_gcc_cswarm     #Give the name of the config.site file

setenv branch_name `git branch | grep '*' | awk '{ print $2; }'`     #To get the name of the current branch
mkdir -pv $PGFEM3D_INSTALL/$branch_name/share
cp $confg_file_name $PGFEM3D_INSTALL/$branch_name/share/config.site

make clean
#make distclean
sh reconf_git_branch.sh
make -j 12
make install
 
echo "... Finished compiling pgfem_3d ... " 
