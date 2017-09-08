This is the current version of PGFem3D, a parallel finite element
solver developed under the advisory of Prof. Karel Matous at the
University of Notre Dame.

PGFem3D relies on several external libraries. They are:

- hypre (v2.40b or above)            : http://computation.llnl.gov/casc/hypre/software.html
- suitesparse (v2.1.1 or above)      : http://faculty.cse.tamu.edu/davis/suitesparse.html
                                     (only the AMD and UMFPACK routines are utilized)
- VTK (v5.10.1 or above) - optional  : http://www.vtk.org/
- Intel MKL                          : included w/ the Intel compiler suite, or can downloaded
- Tensor Template Library (TTL)      : included with the code

Currently we are using hypre (v2.40b), suitesparse (v2.1.1), VTK (v6.0.0) and they are
known to work properly for our cases. However, if you use a newer version of these software
some of them may break our build tree and this can be easily fixed by adjusting library
paths into configure.ac file of PGFem3D.


PGFem3D also depends on the internal library
`Generalizsed_constitutive_model` available at
git@gitlab-cswarm.crc.nd.edu:Constitutive_Model/Generalizsed_constitutive_model.git

We are using the GNU autotools to build this project. If a configure
script does not exist, create/recreate it using 'autoreconf -if'. We
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

### Quick Compilation Guide ###

Sample compilation scripts along with a README file are available at 
/pgfem_3d/SampleBuildScripts. You can simply copy a script to the pgfem_3d
root directory and follow the README file to build PGFem3D. 
