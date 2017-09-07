# Constitutive_model
AC_DEFUN([CONFIG_GCM], [
 AC_ARG_WITH(cnstvm,
   [AS_HELP_STRING([--with-cnstvm=<path>],
                   [Directory for standard Constitutive Model.])],
   [], [with_cnstvm=Generalizsed_constitutive_model])

 gcm_include="-I$with_cnstvm/utils/include -I$with_cnstvm/material/include -I$with_cnstvm/elasticity/include -I$with_cnstvm/crystal_plasticity/include -I$with_cnstvm/data_structure/include"
 gcm_ldflags="-L$with_cnstvm/lib"

 old_CPPFLAGS=$CPPFLAGS
 CPPFLAGS="$CPPFLAGS $gcm_include"
 AC_CHECK_HEADERS([math_help.h], [],
   [AC_MSG_ERROR(Failed to find the Constitutive Model at $with_cnstvm)])
 CPPFLAGS=$old_CPPFLAGS
 
 LDFLAGS="$LDFLAGS $gcm_ldflags"
 AC_CHECK_LIB(ConstitutiveModel, construct_slip_system, [],
   [AC_MSG_ERROR(Failed to link the Constitutive Model at $with_cnstvm)])

 AC_SUBST([GCM_INCLUDE], ["$gcm_include"])
 $1
])
