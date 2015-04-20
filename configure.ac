#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([PGFem3D], [X.X], [bugs@cswarm.nd.edu])
AM_INIT_AUTOMAKE([-Wall foreign])

AC_CONFIG_SRCDIR([src/super_main.c])
AC_CONFIG_HEADERS([include/config.h])

# Checks for programs.
# If --with-mpi=auto is used, try to find MPI, but use standard C compiler if it is not found.
# If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
# If --with-mpi=no is used, use a standard C compiler instead.
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],
		 [compile with MPI (parallelization) support. If none is found,
		 MPI is not used. Default: yes])
		 ],,[with_mpi=yes])
#
AX_PROG_CXX_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
		use_mpi=no
 		if test x"$with_mpi" = xyes; then
		  AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
 		else
		  AC_MSG_WARN([No MPI compiler found, won't use MPI.])
 		fi
])

AX_PROG_CC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
		use_mpi=no
 		if test x"$with_mpi" = xyes; then
		  AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
 		else
		  AC_MSG_WARN([No MPI compiler found, won't use MPI.])
 		fi
])
AC_PROG_CC_STDC
AC_PROG_RANLIB

# Checks for libraries.

AC_ARG_WITH(parmetis-dir, [AS_HELP_STRING([--with-parmetis-dir=<path>],
			      [Provide location of root directory for standard ParMETIS installation.])
			      ],
			      [AC_MSG_CHECKING([location of ParMETIS])
			      AC_MSG_RESULT(["$with_parmetis_dir"])],
			      [AC_MSG_CHECKING([location of ParMETIS])
			      AC_MSG_FAILURE([ParMETIS: Must provide package location!])
			      ])

AC_SUBST([parmetis_include],[-I"$with_parmetis_dir"])
AC_SUBST([parmetis_lib],["-L$with_parmetis_dir -lparmetis -lmetis"])

# HYPRE
AC_ARG_WITH(hypre-dir, [AS_HELP_STRING([--with-hypre-dir=<path>],
			      [Provide location of root directory for standard HYPRE installation.])
			      ],
			      [AC_MSG_CHECKING([location of HYPRE])
			      AC_MSG_RESULT(["$with_hypre_dir"])],
			      [AC_MSG_CHECKING([location of HYPRE])
			      AC_MSG_FAILURE([HYPRE: must provide package location!])
			      ])
AC_SUBST([hypre_include],["-I$with_hypre_dir/include"])
AC_SUBST([hypre_lib],["-L$with_hypre_dir/lib -lHYPRE"])

# MKL
AC_ARG_WITH(mkl, [AS_HELP_STRING([--without-mkl],
			      [Do not use Intel MKL libraries])
			      ],,[with_mkl=yes])
AS_IF([test "x$with_mkl" = xyes],
	    [AC_MSG_NOTICE([Using Intel MKL libraries])],
	    [AC_MSG_FAILURE([Currently must use Intel MKL libraries])])
AC_SUBST([mkl_include],['-I$(MKLROOT)/include'])
AC_SUBST([mkl_lib],['-Wl,--start-group \
$(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a \
$(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
$(MKLROOT)/lib/intel64/libmkl_sequential.a \
$(MKLROOT)/lib/intel64/libmkl_core.a \
$(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
-Wl,--end-group -lpthread'])

# suitesparse
AC_ARG_WITH(suitesparse-dir, [AS_HELP_STRING([--with-suitesparse-dir=<path>],
			      [Specify path to the root directory of a standard suitesparse installation.])
			      ],
			      [AC_MSG_CHECKING([location of SuiteSparse])
			      AC_MSG_RESULT(["$with_suitesparse_dir"])],
			      [AC_MSG_CHECKING([location of SuiteSparse])
			      AC_MSG_FAILURE([Must specify location of SuiteSparse])
			      ])
AC_SUBST([suitesparse_include],["-I$with_suitesparse_dir/UMFPACK/Include \
-I$with_suitesparse_dir/AMD/Include \
-I$with_suitesparse_dir/UFconfig"])
AC_SUBST([suitesparse_lib],["-L$with_suitesparse_dir/UMFPACK/Lib -lumfpack \
-L$with_suitesparse_dir/AMD/Lib -lamd"])

# vtk
AC_ARG_ENABLE(vtk, [AS_HELP_STRING([--enable-vtk],
		 [Output in VTK binary output. Requires additional options
		 --with-vtk-include and --with-vtk-libs. Default:no]
		 )],,[enable_vtk=no]
	    )
AC_ARG_WITH(vtk-include, [AS_HELP_STRING([--with-vtk-include=<include line>],
		 [Provide the include line for VTK.]
		 )],[have_vtk_include=yes],[have_vtk_include=no]
	    )
AC_ARG_WITH(vtk-libs, [AS_HELP_STRING([--with-vtk-libs=<link line>],
		 [Provide the link line for VTK.]
		 )],[have_vtk_libs=yes],[have_vtk_libs=no]
	    )
AS_IF([test "x$enable_vtk" = xyes && test "x$have_vtk_include" = xyes && test "x$have_vtk_libs" = xyes],
dnl using VTK and everything is set
[
AC_MSG_NOTICE([Using VTK Libraries])
AC_SUBST([vtk_include],["$with_vtk_include"])
AC_SUBST([vtk_lib],["$with_vtk_libs"])
],
dnl do not use VTK
[
AC_MSG_NOTICE([No VTK support])
AC_SUBST([vtk_include],[])
AC_SUBST([vtk_lib],[])
AC_DEFINE(NO_VTK_LIB,1,[Do not use VTK libraries])
])

# Checks for header files.
AC_CHECK_HEADERS([limits.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_C_INLINE
#AC_TYPE_MODE_T
#AC_C_RESTRICT
#AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset mkdir pow sqrt strchr strerror strrchr])

AC_CONFIG_FILES([Makefile
                 lib/Renumbering/src/Makefile
                 lib/VTK_IO/src/Makefile
                 src/Makefile])
AC_OUTPUT