AC_DEFUN([CONFIG_SUITESPARSE], [
 AC_ARG_WITH(suitesparse,
  [AS_HELP_STRING([--with-suitesparse=<path>],
			      [Specify path to the root directory of a standard suitesparse installation.])],
  [], [with_suitesparse=no])

 # special case "." argument for source tarball paths
 AS_IF([test "x$with_suitesparse" == x.],
  [suitesparse_include="-I$with_suitesparse/UMFPACK/Include -I$with_suitesparse/AMD/Include -I$with_suitesparse/UFconfig"
   suitesparse_ldflags="-L$with_suitesparse/UMFPACK/Lib -Wl,-rpath,$with_suitesparse/UMFPACK/Lib -L$with_suitesparse/AMD/Lib -Wl,-rpath,$with_suitesparse/AMD/Lib"])

 # else we assume path points to an install target
 AS_IF([test "x$with_suitesparse" != xno],
    [suitesparse_include="-I$with_suitesparse/include"
     suitesparse_ldflags="-L$with_suitesparse/lib -Wl,-rpath,$with_suitesparse/lib"])

 old_CPPFLAGS=$CPPFLAGS
 CPPFLAGS="$CPPFLAGS $suitesparse_include"
 AC_CHECK_HEADER([umfpack.h], [], [AC_MSG_ERROR(Can't find SuiteSparse headers)])
 CPPFLAGS=$old_CPPFLAGS  

 LDFLAGS="$LDFLAGS $suitesparse_ldflags"
 AC_SEARCH_LIBS([umfpack_dl_solve], [umfpack], [],
   [AC_MSG_ERROR(Can't link SuiteSparse as $suitesparse_ldflags -lumfpack)])

 AC_SUBST([SUITESPARSE_INCLUDE], [$suitesparse_include])
 $1
])
