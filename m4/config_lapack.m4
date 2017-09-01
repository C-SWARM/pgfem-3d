AC_DEFUN([CONFIG_LAPACK], [
 AC_ARG_ENABLE(mkl,
  [AS_HELP_STRING([--enable-mkl], [Use Intel MKL @<:@default=yes@:>@])],
  [], [enable_mkl=yes])

 AC_ARG_WITH(mkl,
  [AS_HELP_STRING([--with-mkl], [Path to Intel MKL @<:@default=$MKLROOT@:>@])],
  [], [with_mkl=$MKLROOT])

 mkl_include="-m64 -I$with_mkl/include $with_mkl_override"
 mkl_ldflags="-L$with_mkl/lib/intel64 -Wl,-rpath,$with_mkl/lib/intel64"
 mkl_libs="-Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"

 AC_ARG_WITH(mkl-override,
  [AS_HELP_STRING([--with-mkl-override=<include line>],
			      [Specify include line for MKL override headers])])

 AC_ARG_WITH(cblas-include,
  [AS_HELP_STRING([--with-cblas-include=<include line>],
			      [Specify include line for access to CBLAS])])
                  
 AC_ARG_WITH(cblas-lib,
  [AS_HELP_STRING([--with-cblas-lib=<link line>],
			      [Specify link line for access to CBLAS])])
                  
 AC_ARG_WITH(lapack-include,
  [AS_HELP_STRING([--with-lapack-include=<include line>],
			      [Specify include line for access to LAPACK])])
                  
 AC_ARG_WITH(lapack-lib,
  [AS_HELP_STRING([--with-lapack-lib=<link line>],
			      [Specify link line for access to LAPACK])])

 AS_IF([test "x$enable_mkl" == xyes],
  [AC_MSG_NOTICE([Using Intel MKL libraries])
   cblas_h="mkl_cblas.h"
   lapack_h="mkl_lapack.h"
   lapack_include="$mkl_include"
   lapack_ldflags="$mkl_ldflags"
   lapack_libs="$mkl_libs"],
  [AC_MSG_NOTICE([Using specified CBLAS/LAPACK])
   cblas_h="cblas.h"
   lapack_h="lapack.h"  
   lapack_include="$with_mkl_override $with_cblas_include $with_lapack_include"
   lapack_libs="$with_cblas_lib $with_lapack_lib"])

 old_CPPFLAGS=$CPPFLAGS
 CPPFLAGS="$CPPFLAGS $lapack_include"
 AC_CHECK_HEADERS([$cblas_h $lapack_h], [], [AC_MSG_ERROR(Could not find cblas headers)])
 CPPFLAGS=$old_CPPFLAGS
 
 dnl old_LDFLAGS=$LDFLAGS
 dnl LDFALGS="$LDFLAGS $lapack_ldflags"
 dnl AC_SEARCH_LIBS(cblas_dgemm, [$lapack_libs], [], [AC_MSG_ERROR(Failed to link cblas)])
 dnl LDFLAGS=$old_LDFLAGS

 AC_SUBST([CBLAS_H], ["$cblas_h"])
 AC_SUBST([LAPACK_H], ["$lapack_h"])
 AC_SUBST([LAPACK_INCLUDE], ["$lapack_include"])
 AC_SUBST([LAPACK_LDFLAGS], ["$lapack_ldflags"])
 AC_SUBST([LAPACK_LIBS], ["$lapack_libs"])
 $1
])
