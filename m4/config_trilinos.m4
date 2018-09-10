AC_DEFUN([CONFIG_TRILINOS], [
 # allow the user to enable or disable trilinos
 AC_ARG_ENABLE([trilinos],
  [AS_HELP_STRING([--enable-trilinos],
                  [Enable the use of trilinos @<:@default=yes@:>@])],
  [], [enable_trilinos=yes])

 # allow the user to override where trilinos is
 AC_ARG_WITH(trilinos,
  [AS_HELP_STRING([--with-trilinos=<path>],
                  [Provide location of root directory for standard TRILINOS installation.])])

 # turn the with_trilinos string into useful path variables.
 AS_IF([test "x$with_trilinos" != x],
  [trilinos_include="-I$with_trilinos/include"
   trilinos_ldflags="-L$with_trilinos/lib -Wl,-rpath,$with_trilinos/lib"])

 # Trilinos export
 #AC_ARG_WITH(trilinos,
 # [AS_HELP_STRING([--trilinos-inc=<path>],
#		  [Provide location of Trilinos inc])])

 #AS_IF([test "x$trilinos_inc" != x],
  #[trilinos_include="$trilinos_include -I$trilinos_inc"])

 # output trilinos-related build variables
 AS_IF([test "x$enable_trilinos" == xyes],
  [old_CPPFLAGS=$CPPFLAGS
   CPPFLAGS="$CPPFLAGS $trilinos_include"
   AC_CHECK_HEADER([Trilinos_version.h], [],
     [AC_MSG_ERROR(Can't find TRILINOS headers $trilinos_include)])
   CPPFLAGS=$old_CPPFLAGS

   #LDFLAGS="$LDFLAGS $trilinos_ldflags"
   #AC_SEARCH_LIBS([_GLOBAL__sub_I__ZN5Belos13Belos_VersionB5cxx11Ev], [], [],
   #  [AC_MSG_ERROR(Can't link TRILINOS with $trilinos_ldflags -lbelos)])
   #$1])

 AC_SUBST([TRILINOS_INCLUDE], ["$trilinos_include"])
])
