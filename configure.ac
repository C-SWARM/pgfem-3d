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
		  AC_MSG_WARN([Will not use MPI.])
 		fi
])

AX_PROG_CC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
		use_mpi=no
 		if test x"$with_mpi" = xyes; then
		  AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
 		else
		  AC_MSG_WARN([Will not use MPI.])
 		fi
])

AX_PROG_FC_MPI([test x"$with_mpi" != xno],[use_mpi=yes],[
		use_mpi=no
		if test x"$with_mpi" = xyes; then
		  AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
		else
		  AC_MSG_WARN([Will not use MPI.])
		fi
])

AC_PROG_CC_STDC
AC_PROG_RANLIB

AC_ARG_WITH(special-link, [AS_HELP_STRING([--with-special-link=<linker>],
		 [Link the executable with the specified linking tool. Default: "$(CXX)"])
		 ],[special_link=yes],[special_link=no])

AS_IF([test "x$special_link" = xyes],
	    [
	    AC_MSG_NOTICE([Linking PGFem3D with "$with_special_link"])
	    AC_SUBST([PGFem3D_linker],["$with_special_link"])
	    ],[
	    AC_SUBST([PGFem3D_linker],['$(CXX)'])
	    ])

# Checks for libraries.
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

# Constitutive_model
AC_ARG_WITH(cnstvm-dir, [AS_HELP_STRING([--with-cnstvm-dir=<path>],
                              [Provide location of root directory for standard Constitutive Model installation.])
                              ],
                              [AC_MSG_CHECKING([location of Constitutive Model])
                              AC_MSG_RESULT(["$with_cnstvm_dir"])],
                              [AC_MSG_CHECKING([location of Constitutive Model])
                              AC_MSG_FAILURE([Constitutive Model: must provide package location!])
                              ])
AC_SUBST([cnstvm_include],["-I$with_cnstvm_dir/include"])
AC_SUBST([cnstvm_lib],["-L$with_cnstvm_dir/src -lCM "])
#AC_SUBST([cnstvm_lib],["-L$with_cnstvm_dir/src -lCM -L/opt/crc/scilib/lapack/gnu -llapack -llapack_pic"])


# MKL
AC_ARG_WITH(mkl, [AS_HELP_STRING([--without-mkl],
			      [Do not use Intel MKL libraries])
			      ],,[with_mkl=yes])

# CBLAS / LAPACK
AC_ARG_WITH(mkl-override, [AS_HELP_STRING([--with-mkl-override=<include line>],
			      [Specify include line for MKL override headers (requires --without-mkl)])
			      ],,[with_mkl_override=''])
AC_ARG_WITH(cblas-include, [AS_HELP_STRING([--with-cblas-include=<include line>],
			      [Specify include line for access to CBLAS (requires --without-mkl)])
			      ],,[with_cblas_include=''])
AC_ARG_WITH(cblas-lib, [AS_HELP_STRING([--with-cblas-lib=<link line>],
			      [Specify link line for access to CBLAS (requires --without-mkl)])
			      ],,[with_cblas_lib=''])
AC_ARG_WITH(lapack-include, [AS_HELP_STRING([--with-lapack-include=<include line>],
			      [Specify include line for access to LAPACK (requires --without-mkl)])
			      ],,[with_lapack_include=''])
AC_ARG_WITH(lapack-lib, [AS_HELP_STRING([--with-lapack-lib=<link line>],
			      [Specify link line for access to LAPACK (requires --without-mkl)])
			      ],,[with_lapack_lib=''])

AS_IF([test "x$with_mkl" = xyes],
	    [
	    AC_MSG_NOTICE([Using Intel MKL libraries])
	    AC_SUBST([mkl_include],['-I$(MKLROOT)/include'])
	    AC_SUBST([mkl_lib],['-Wl,--start-group \
	    $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
	    $(MKLROOT)/lib/intel64/libmkl_sequential.a \
	    $(MKLROOT)/lib/intel64/libmkl_core.a \
	    -Wl,--end-group -lpthread'])
	    ],
	    [
	    AC_MSG_NOTICE([Using specified CBLAS/LAPACK])
	    AC_SUBST([mkl_include],["$with_mkl_override $with_cblas_include $with_lapack_include"])
	    AC_SUBST([mkl_lib],["$with_cblas_lib $with_lapack_lib"])
	    ])

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
                 lib/VTK_IO/src/Makefile
                 src/Makefile])
AC_OUTPUT
