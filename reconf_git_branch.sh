#!/bin/bash

# This is a script to aid in configuring/installing for each git branch.

# get the name of the current git branch
branch_name=$(git branch | grep '*' | sed -e 's/* //')

# ensure that PGFEM3D_INSTALL is set
if [ -z $PGFEM3D_INSTALL ]; then
    echo "please set PGFEM3D_INSTALL and try again"
    exit
fi

# ensure that the install directory exists
git_install_prefix=$PGFEM3D_INSTALL/$branch_name
mkdir -pv $git_install_prefix/share

# notify user if prefix/share/config.site doesn't exist
if [ ! -e $git_install_prefix/share/config.site ]; then
    echo "WARNING: $git_install_prefix/share/config.site does not exist!"
fi

echo "Ensuring that the correct scripts are used"
autoreconf -if

echo "Configuring with ./configure --prefix=$git_install_prefix"
./configure --prefix=$git_install_prefix
