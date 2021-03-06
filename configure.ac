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
: ${CXXFLAGS="-O3"}

AC_LANG([C++])

AC_CONFIG_SRCDIR([src/super_main.cc])
AC_CONFIG_HEADERS([include/config.h])
AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_MACRO_DIR([m4])

AX_PROG_CXX_MPI([test "x$with_mpi" != xno], [have_mpicxx=yes])
AX_CXX_COMPILE_STDCXX_14([ext],[mandatory])

# Initialize autotools and libtool
AM_INIT_AUTOMAKE([1.11 nostdinc dist-bzip2 subdir-objects foreign tar-ustar -Wall
                  -Werror parallel-tests color-tests])
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
LT_PREREQ([2.4.2])
LT_INIT([disable-shared])

# Allow parallel configuration. We do this early so that it shows up before the
# rest of the 'enable' options. 
AC_ARG_ENABLE([parallel-config],
 [AS_HELP_STRING([--enable-parallel-config],
                 [Enable parallel configuration @<:@default=no@:>@])],
 [], [enable_parallel_config=no])

CONFIG_LAPACK([have_lapack=yes])
CONFIG_TTL([have_ttl=yes])
CONFIG_MTL([have_mtl=yes])
CONFIG_HYPRE([have_hypre=yes])
CONFIG_TRILINOS([have_trilinos=yes])
CONFIG_SUITESPARSE([have_suitesparse=yes])
CONFIG_VTK([5.10.1])
CONFIG_GCM([have_gcm=yes])
CONFIG_MSNET([contrib/msnet],[msnet])
CONFIG_TESTS()

AC_ARG_WITH(special-link,
  [AS_HELP_STRING([--with-special-link=<linker>],
		         [Link the executable with the specified linking tool @<:@default=$(CXX)@:>@])],
  [special_link=yes],
  [special_link=no])

AS_IF([test "x$special_link" = xyes],
  [AC_MSG_NOTICE([Linking PGFem3D with "$with_special_link"])
   AC_SUBST([PGFem3D_linker],["$with_special_link"])],
  [AC_SUBST([PGFem3D_linker],['$(CXX)'])])

AS_IF([test "x$have_mpicxx" = xyes],
  [AC_DEFINE([HAVE_PGFEM_MPI], [1], [MPI support available])]
)

# set automake conditionals that our network build requires
AM_CONDITIONAL([HAVE_PGFEM_MPI], [test "x$have_mpicxx" == xyes])
AM_CONDITIONAL([HAVE_MSNET], [test "x$have_msnet" == xyes])
AM_CONDITIONAL([BUILD_MSNET], [test "x$build_msnet" == xyes])

AS_IF([test "x$have_trilinos" = xyes],
	    [AC_DEFINE([HAVE_TRILINOS], [1], [Trilinos solvers avaiable])])

AM_CONDITIONAL([HAVE_TRILINOS], [test "x$have_trilinos" == xyes])

# https://gitlab-cswarm.crc.nd.edu/Constitutive_Model/Generalizsed_constitutive_model/issues/9
dnl LDFLAGS="$LDFLAGS -Wl,--allow-multiple-definition"

AC_CONFIG_FILES([Makefile
		 contrib/Makefile
                 lib/VTK_IO/src/Makefile
                 energy_equation/src/Makefile
                 FEMLIB/src/Makefile
                 src/Makefile
                 src/solvers/Makefile
                 src/solvers/hypre/Makefile
		 src/solvers/trilinos/Makefile
		 src/multiscale/Makefile
		 src/communication/Makefile
                 utils/Makefile
                 verification_MMS/Makefile
                 verification_tests/Makefile
                 verification_tests/src/Makefile
                 tests/Makefile
                 tests/box/Makefile
                 tests/box/compression/Makefile
                 tests/box/shear/Makefile
                 tests/box/tension_disp/Makefile
                 tests/box/tension_pressure/Makefile
                 tests/cohesive/Makefile
                 tests/constitutive_model/Makefile
                 tests/crystal_plasticity/Makefile
                 tests/multiphysics/Makefile
                 tests/time_stepping/Makefile
                 tests/fe2/Makefile])

# Configure the test scripts. We do this even if running the scripts has been
# disabled. 
AC_CONFIG_FILES([tests/box/compression/run.sh
                 tests/box/common.sh
                 tests/cohesive/common.sh
                 tests/time_stepping/4.cpu.restart.sh
                 tests/time_stepping/4.cpu.different_dt_tr.sh
                 tests/crystal_plasticity/common.sh
                 tests/multiphysics/const_T_with_constrained_BC.sh
                 tests/multiphysics/heat_transfer_on_deformed_config.sh
                 tests/fe2/run.sh],
                [chmod +x tests/box/compression/run.sh
                 chmod +x tests/box/common.sh
                 chmod +x tests/cohesive/common.sh
                 chmod +x tests/time_stepping/4.cpu.restart.sh
                 chmod +x tests/time_stepping/4.cpu.different_dt_tr.sh
                 chmod +x tests/crystal_plasticity/common.sh
                 chmod +x tests/multiphysics/const_T_with_constrained_BC.sh
                 chmod +x tests/multiphysics/heat_transfer_on_deformed_config.sh
                 chmod +x tests/fe2/run.sh])

# For the crystal plasticity test case we need to move the files that are in
# tests/crystal_plasticity/CO into the build tree. They are needed by the test
# and I don't know how to specify their path to PGFem3D (it expects them
# relative to the execution directory). I don't know how to copy the entire
# directory so we just copy each file manually.
AC_CONFIG_FILES([tests/crystal_plasticity/CO/co_0.in:verification_tests/crystal_plasticity/CO/co_0.in
                 tests/crystal_plasticity/CO/co_1.in:verification_tests/crystal_plasticity/CO/co_1.in
                 tests/crystal_plasticity/CO/co_2.in:verification_tests/crystal_plasticity/CO/co_2.in
                 tests/crystal_plasticity/CO/co_3.in:verification_tests/crystal_plasticity/CO/co_3.in
                 tests/crystal_plasticity/CO/co_4.in:verification_tests/crystal_plasticity/CO/co_4.in
                 tests/crystal_plasticity/CO/co_5.in:verification_tests/crystal_plasticity/CO/co_5.in
                 tests/crystal_plasticity/CO/co_6.in:verification_tests/crystal_plasticity/CO/co_6.in
                 tests/crystal_plasticity/CO/co_7.in:verification_tests/crystal_plasticity/CO/co_7.in])

# The verification_MMS build wants to use objects from the /src directory, but
# with different flags (-D). The easiest way I can think of to do this is to
# build links to the files in the verification_MMS directory and depend on those
# different sources.
AC_CONFIG_LINKS([verification_MMS/dynamics.cc:src/dynamics.cc
                 verification_MMS/super_main.cc:src/super_main.cc])

AC_OUTPUT
