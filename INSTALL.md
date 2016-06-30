This is the current version of PGFem3D, a parallel finite element
solver developed under the advisory of Prof. Karel Matous at the
University of Notre Dame.

PGFem3D relies on several external libraries. They are:

VTK-6.0	      	       : http://vtk.org/VTK/resources/software.html
hypre-2.40b   	       : http://computation.llnl.gov/casc/hypre/software.html
suitesparse-2.1.1      : http://www.cise.ufl.edu/research/sparse/SuiteSparse/
           	         (only the AMD and UMFPACK routines are utilized)
Intel MKL	       : included w/ the Intel compiler suite

While many of these libraries have newer versions available, these
versions will not break our build tree and are known to work properly
for our cases.

PGFem3D also depends on the internal library
`Generalizsed_constitutive_model` available at
git@gitlab-cswarm.crc.nd.edu:Constitutive_Model/Generalizsed_constitutive_model.git

We are using the GNU autotools to build this project. If a configure
script does not exist, recreate it using 'autoreconf -if'. We
currently requrie the following options to be specified in the
configuration step:

  --with-hypre-dir=<path>
  --with-suitesparse-dir=<path>
  --with-cnstvm-dir=<path>

Use of the VTK libraries is off by default but may be enabled with the
following options:

  --enable-vtk
  --with-vtk-include=<include line>
  --with-vtk-libs=<link line>

See `./configure -h` for more information. After configuration, simply
`make && make install`.

### Tips for Developers ###

It is suggested that you set a config.site file for each Git branch
you are actively developing such as

`.../pgfem3d_install/branch_name/share/config.site`

then options will automatically be set when you envoke

`./configure --prefix=.../pgfem3d_install/branch_name`

A sample `config.site` might look like the following:

# required options
with_hypre_dir=/path/to/hypre-2.4.0b
with_suitesparse_dir=/path/to/suitesparse
with_cnstvm_dir=/path/to/Generalizsed_constitutive_model

# enable/disable vtk
enable_vtk=yes

# required extra opts for enable_vtk=yes
vtk_dir=/path/to/vtk
with_vtk_include="-I$vtk_dir/include/vtk-6.0"
with_vtk_libs="-Wl,-rpath=$vtk_dir/lib -L$vtk_dir/lib \
-lvtkIOXML-6.0 -lvtkIOXMLParser-6.0 -lvtkIOCore-6.0 -lvtkzlib-6.0 \
-lvtkCommonExecutionModel-6.0 -lvtkCommonDataModel-6.0 \
-lvtkCommonMisc-6.0 -lvtkCommonSystem-6.0 -lvtkCommonTransforms-6.0 \
-lvtkCommonMath-6.0 -lvtkIOGeometry-6.0 -lvtkCommonCore-6.0 \
-lvtksys-6.0 -lvtkexpat-6.0"