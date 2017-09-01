AC_DEFUN([CONFIG_HYPRE], [
 # allow the user to enable or disable hypre
 AC_ARG_ENABLE([hypre],
  [AS_HELP_STRING([--enable-hypre],
                  [Enable the use of hypre @<:@default=yes@:>@])],
  [], [enable_hypre=yes])

 # allow the user to override where hypre is
 AC_ARG_WITH(hypre,
  [AS_HELP_STRING([--with-hypre=<path>],
                  [Provide location of root directory for standard HYPRE installation.])])

 # turn the with_hypre string into useful path variables.
 AS_IF([test "x$with_hypre" != x],
  [hypre_include="-I$with_hypre/include"
   hypre_ldflags="-L$with_hypre/lib -Wl,-rpath,$with_hypre/lib"])

 # output hypre-related build variables
 AS_IF([test "x$enable_hypre" == xyes],
  [AC_SUBST([HYPRE_INCLUDE], ["$hypre_include"])
   AC_SUBST([HYPRE_LDFLAGS], ["$hypre_ldflags"])
   AC_SUBST([HYPRE_LIBS], ["-lHYPRE"])
   $1],
  [$2])
])
