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
  [old_CPPFLAGS=$CPPFLAGS
   CPPFLAGS="$CPPFLAGS $hypre_include"
   AC_CHECK_HEADER([HYPRE.h], [],
     [AC_MSG_ERROR(Can't find HYPRE headers $hypre_include)])
   CPPFLAGS=$old_CPPFLAGS

   LDFLAGS="$LDFLAGS $hypre_ldflags"
   AC_SEARCH_LIBS([HYPRE_IJMatrixCreate], [HYPRE], [],
     [AC_MSG_ERROR(Can't link HYPRE with $hypre_ldflags -lHYPRE)])
   $1])

 AC_SUBST([HYPRE_INCLUDE], ["$hypre_include"])
])
