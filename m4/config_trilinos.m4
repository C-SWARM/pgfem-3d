AC_DEFUN([CONFIG_TRILINOS], [
 # allow the user to enable or disable trilinos
 AC_ARG_ENABLE([trilinos],
  [AS_HELP_STRING([--enable-trilinos],
                  [Enable the use of trilinos @<:@default=no@:>@])],
		  [enable_trilinos=yes], [enable_trilinos=no])

 

 # allow the user to override where trilinos is
 AC_ARG_WITH(trilinos,
  [AS_HELP_STRING([--with-trilinos=<path>],
                  [Provide location of root directory for standard TRILINOS installation.])])

 # turn the with_trilinos string into useful path variables.
 AS_IF([test "x$with_trilinos" != x],
  [trilinos_include="-I$with_trilinos/include"
   trilinos_ldflags="-L$with_trilinos/lib -Wl,-rpath,$with_trilinos/lib"])

 # output trilinos-related build variables
 AS_IF([test "x$enable_trilinos" == xyes],
  [old_CPPFLAGS=$CPPFLAGS
   CPPFLAGS="$CPPFLAGS $trilinos_include"
   AC_CHECK_HEADER([Trilinos_version.h], [],
     [AC_MSG_ERROR(Can't find TRILINOS headers $trilinos_include)])
   CPPFLAGS=$old_CPPFLAGS


  TRILINOS_LIBADD="-lmuelu-adapters -lmuelu-interface -lmuelu -lintrepid2 -lstratimikos -lstratimikosbelos -lstratimikosamesos2 -lifpack2-adapters -lifpack2 -lamesos2 -lshylu_nodetacho -lshylu_nodehts -lshylu_nodetacho -lshylu_nodehts -lbelosxpetra -lbelostpetra -lbelos -lml -lzoltan2 -lpamgen_extras -lpamgen -lgaleri-xpetra -lxpetra-sup -lxpetra -lthyratpetra  -lthyracore -lthyratpetra -lthyracore -ltrilinosss -ltpetraext -ltpetrainout -ltpetra -lkokkostsqr -ltpetraclassiclinalg -ltpetraclassicnodeapi -ltpetraclassic -ltpetraext -ltpetrainout -ltpetra -lkokkostsqr -ltpetraclassiclinalg -ltpetraclassicnodeapi -ltpetraclassic -lshards -lzoltan -lsacado -lrtop -lkokkoskernels -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchosparser -lteuchoscore -lteuchoskokkoscomm -lteuchoskokkoscompat -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchosparser -lteuchoscore -lkokkosalgorithms -lkokkoscontainers -lkokkoscore -lkokkosalgorithms -lkokkoscontainers -lkokkoscore -lgtest"
   
 LDFLAGS="$LDFLAGS $trilinos_ldflags"
   
 $1])

 AC_SUBST([TRILINOS_LIBADD], ["$TRILINOS_LIBADD"])
 AC_SUBST([TRILINOS_INCLUDE], ["$trilinos_include"])
])
