# Constitutive_model
AC_DEFUN([CONFIG_GCM], [
 AC_ARG_WITH(cnstvm,
   [AS_HELP_STRING([--with-cnstvm=<path>],
                   [Directory for standard Constitutive Model.])],
   [], [with_cnstvm=Generalizsed_constitutive_model])

 gcm_include="-I$with_cnstvm/utils/include -I$with_cnstvm/material/include -I$with_cnstvm/elasticity/include -I$with_cnstvm/crystal_plasticity/include"
 gcm_ldflags="-L$with_cnstvm/lib"
 gcm_libs="-lConstitutiveModel"

 old_CPPFLAGS=$CPPFLAGS
 CPPFLAGS="$CPPFLAGS $gcm_include"
 AC_CHECK_HEADERS([math_help.h], [],
   [AC_MSG_ERROR(Failed to find the Constitutive Model at $with_cnstvm)])
 CPPFLAGS=$old_CPPFLAGS
 
 old_LDFLAGS=$LDFLAGS
 LDFLAGS="$LDFLAGS $gcm_ldflags"
 AC_CHECK_LIB(ConstitutiveModel, construct_slip_system, [],
   [AC_MSG_ERROR(Failed to link the Constitutive Model at $with_cnstvm)])
 LDFLAGS=$old_LDFLAGS

 AC_SUBST([GCM_INCLUDE], ["$gcm_include"])
 AC_SUBST([GCM_LDFLAGS], ["$gcm_ldflags"])
 AC_SUBST([GCM_LIBS], ["$gcm_libs"])
 $1
])
