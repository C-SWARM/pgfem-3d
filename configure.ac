#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
m4_define([PGFEM3D_MAJOR], [0])
m4_define([PGFEM3D_MINOR], [1])
AC_INIT([PGFem3D], [PGFEM3D_MAJOR.PGFEM3D_MINOR], [bugs@cswarm.nd.edu])

# WE OVERWRITE "-O2 -g" default flags, because there is no world where
# ./configure CFLAGS="-O2 -g" is better than ./configure CXXFLAGS="-O3 -g" as a
# default build.
# - uncomment if you need C
# : ${CFLAGS="-O3 -g"}
: ${CXXFLAGS="-O3 -g"}

AC_LANG([C++])

AC_CONFIG_SRCDIR([src/super_main.cc])
AC_CONFIG_HEADERS([include/config.h])
AC_CONFIG_MACRO_DIRS([m4])

AC_PROG_RANLIB
# AX_PROG_CC_MPI([test "x$with_mpi" != xno], [have_mpicc=yes])
# AC_PROG_CC_STDC
# AX_PROG_FC_MPI([test "x$with_mpi" != xno], [have_mpifc=yes])
AX_PROG_CXX_MPI([test "x$with_mpi" != xno], [have_mpicxx=yes])
AX_CXX_COMPILE_STDCXX_14([ext],[mandatory])
AC_C_RESTRICT                                   # Works for C++ too

# check for the ability to build tests
AC_CHECK_PROG([found_numdiff], [numdiff], [yes])

# Initialize autotools and libtool
AM_INIT_AUTOMAKE([1.11 nostdinc dist-bzip2 subdir-objects foreign tar-ustar -Wall
                  -Werror parallel-tests color-tests])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

# Allow the user to disable MPI.
AC_ARG_WITH(mpi,
  [AS_HELP_STRING([--with-mpi{=yes,no,auto}],
                  [Compile with MPI (parallelization) support @<:@default=yes@:>@])],
  [], [with_mpi=yes])

# Consolidated MPI checks
# AS_IF([test "x$have_mpicxx" != xyes -o "x$have_mpicc" != xyes],
# AS_IF([test "x$have_mpicxx" != xyes -o "x$have_mpicc" != xyes -o "x$have_mpifc" != xyes],
AS_IF([test "x$have_mpicxx" != xyes],
  [AS_CASE($with_mpi,
   [yes], [AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])],
   [auto], [AC_MSG_WARN([No MPI compiler, disabling MPI])],
   [])])

AC_ARG_WITH(special-link,
  [AS_HELP_STRING([--with-special-link=<linker>],
		         [Link the executable with the specified linking tool @<:@default=$(CXX)@:>@])],
  [special_link=yes],
  [special_link=no])

AS_IF([test "x$special_link" = xyes],
  [AC_MSG_NOTICE([Linking PGFem3D with "$with_special_link"])
   AC_SUBST([PGFem3D_linker],["$with_special_link"])],
  [AC_SUBST([PGFem3D_linker],['$(CXX)'])])

# Checks for libraries.
# HYPRE
AC_ARG_WITH(hypre-dir,
  [AS_HELP_STRING([--with-hypre-dir=<path>],
			      [Provide location of root directory for standard HYPRE installation.])],
  [AC_MSG_CHECKING([location of HYPRE])
   AC_MSG_RESULT(["$with_hypre_dir"])],
  [with_hypre_dir=no])

AS_IF([test "x$with_hypre_dir" != xno],
  [AC_SUBST([hypre_include],["-I$with_hypre_dir/include"])
   AC_SUBST([hypre_lib],["-L$with_hypre_dir/lib -lHYPRE"])],
  [AC_SUBST([hypre_lib],["-lHYPRE"])])

# Constitutive_model
AC_ARG_WITH(cnstvm-dir,
  [AS_HELP_STRING([--with-cnstvm-dir=<path>],
                  [Provide location of root directory for standard Constitutive Model installation.])],
  [], [AC_MSG_FAILURE([--with-cnstvm-dir required!])])

AC_MSG_CHECKING([location of Constitutive Model])
AC_MSG_RESULT(["$with_cnstvm_dir"])
AC_SUBST([cnstvm_include],["-I$with_cnstvm_dir/utils/include -I$with_cnstvm_dir/material/include -I$with_cnstvm_dir/elasticity/include -I$with_cnstvm_dir/crystal_plasticity/include"])
AC_SUBST([cnstvm_lib],["-L$with_cnstvm_dir/lib -lConstitutiveModel"])

# MKL
AC_ARG_WITH(mkl,
  [AS_HELP_STRING([--without-mkl],
			      [Do not use Intel MKL libraries])],
  [], [with_mkl=yes])

# CBLAS / LAPACK
AC_ARG_WITH(mkl-override,
  [AS_HELP_STRING([--with-mkl-override=<include line>],
			      [Specify include line for MKL override headers (requires --without-mkl)])])
AC_ARG_WITH(cblas-include,
  [AS_HELP_STRING([--with-cblas-include=<include line>],
			      [Specify include line for access to CBLAS (requires --without-mkl)])])
                  
AC_ARG_WITH(cblas-lib,
  [AS_HELP_STRING([--with-cblas-lib=<link line>],
			      [Specify link line for access to CBLAS (requires --without-mkl)])])
                  
AC_ARG_WITH(lapack-include,
  [AS_HELP_STRING([--with-lapack-include=<include line>],
			      [Specify include line for access to LAPACK (requires --without-mkl)])])
                  
AC_ARG_WITH(lapack-lib,
  [AS_HELP_STRING([--with-lapack-lib=<link line>],
			      [Specify link line for access to LAPACK (requires --without-mkl)])])

AS_IF([test "x$with_mkl" = xyes],
  [AC_MSG_NOTICE([Using Intel MKL libraries])
   AC_SUBST([mkl_include],['-m64 -I$(MKLROOT)/include'])
   AC_SUBST([mkl_lib],['-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl'])],
  [AC_MSG_NOTICE([Using specified CBLAS/LAPACK])
   AC_SUBST([mkl_include],["$with_mkl_override $with_cblas_include $with_lapack_include"])
   AC_SUBST([mkl_lib],["$with_cblas_lib $with_lapack_lib"])])

# suitesparse
AC_ARG_WITH(suitesparse-dir,
  [AS_HELP_STRING([--with-suitesparse-dir=<path>],
			      [Specify path to the root directory of a standard suitesparse installation.])],
  [AC_MSG_CHECKING([location of SuiteSparse])
   AC_MSG_RESULT(["$with_suitesparse_dir"])],
  [with_suitesparse_dir=no])

# special case "." argument for source tarball paths
AS_IF([test "x$with_suitesparse_dir" == x.],
  [AC_SUBST([suitesparse_include],["-I$with_suitesparse_dir/UMFPACK/Include -I$with_suitesparse_dir/AMD/Include -I$with_suitesparse_dir/UFconfig"])
   AC_SUBST([suitesparse_lib],["-L$with_suitesparse_dir/UMFPACK/Lib -L$with_suitesparse_dir/AMD/Lib -lumfpack -lamd"])],
  # else we assume path points to an install target
  [AS_IF([test "x$with_suitesparse_dir" != xno],
    [AC_SUBST([suitesparse_include],["-I$with_suitesparse_dir/include"])
     AC_SUBST([suitesparse_lib],["-L$with_suitesparse_dir/lib -lumfpack -lamd"])],
    # else fallback to user environment
    [AC_SUBST([suitesparse_lib],["-lumfpack -lamd"])])])

# TTL
AC_ARG_WITH(ttl-dir,
  [AS_HELP_STRING([--with-ttl-dir=<path>],
			      [Specify path to TTL installation.])],
  [AC_MSG_CHECKING([location of TTL])
   AC_MSG_RESULT(["$with_ttl_dir"])],
  [with_ttl_dir=no])

AS_IF([test "x$with_ttl_dir" != xno],
  [AC_SUBST([ttl_include],["-I$with_ttl_dir/include"])],
  [])

# This macro will provide options for the user to enable and specify the path
# and version of VTK, and will do some configuration-time searching to see if
# these paths work and to set VTK_CPPFLAGS, VTK_LDFLAGS, and VTK_LIBS.
PGFEM3D_VTK([5.10.1], [have_vtk=yes])

AS_IF([test "x$have_vtk" != xyes],
  [AC_DEFINE(NO_VTK_LIB,1,[Do not use VTK libraries])],
  [AC_SUBST([VTK_CPPFLAGS], ["$VTK_CPPFLAGS"])
   AC_SUBST([VTK_LDFLAGS], ["$VTK_LDFLAGS"])
   AC_SUBST([VTK_LIBS], ["$VTK_LIBS"])])

AM_CONDITIONAL([HAVE_VTK], [test "x$have_vtk" == xyes])

#tests 
AC_ARG_ENABLE(tests,
  [AS_HELP_STRING([--enable-tests],
		          [Enable tests (requires numdiff) @<:@default=yes@:>@])],
  [], [enable_tests=yes])

AC_ARG_WITH(tests-nprocs,
  [AS_HELP_STRING([--with-tests-nprocs=<n>],
                  [The maximum number of processes to use during testing @<:@default=1@:>@])],
  [], [with_tests_nprocs=1])

AC_SUBST([TESTS_NPROCS], [$with_tests_nprocs])

# Checks for header files.
AC_CHECK_HEADERS([limits.h stdlib.h string.h sys/time.h unistd.h])

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset mkdir pow sqrt strchr strerror strrchr])

# Deal with tests
AS_IF([test "x$enable_tests" == xyes -a "x$found_numdiff" != xyes],
  [AC_MSG_WARN([Missing numdiff, disabling tests])],
  [AC_CONFIG_FILES([tests/box_compression/run.sh
                    tests/box_shear/run.sh
                    tests/box_tension_disp/run.sh],
                   [chmod +x tests/box_compression/run.sh
                    chmod +x tests/box_shear/run.sh
                    chmod +x tests/box_tension_disp/run.sh])])

AM_CONDITIONAL([BUILD_TESTS], [test "x$enable_tests" == xyes -a "x$found_numdiff" == xyes])

AC_CONFIG_FILES([Makefile
                 lib/VTK_IO/src/Makefile
                 energy_equation/src/Makefile
                 src/Makefile
                 verification_MMS/Makefile
                 verification_tests/Makefile
                 verification_tests/src/Makefile
                 tests/Makefile
                 tests/box_compression/Makefile
                 tests/box_shear/Makefile
                 tests/box_tension_disp/Makefile])

AC_OUTPUT
