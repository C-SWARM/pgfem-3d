# Constitutive_model
AC_DEFUN([CONFIG_GCM], [
 AC_ARG_WITH(cnstvm,
   [AS_HELP_STRING([--with-cnstvm=<path>],
                   [Directory for standard Constitutive Model.])],
 [], [with_cnstvm=Generalizsed_constitutive_model])

 LIBGCM_CPPFLAGS="-I$with_cnstvm/constitutive_model_handle/include -I$with_cnstvm/utils/include -I$with_cnstvm/material/include -I$with_cnstvm/elasticity/include -I$with_cnstvm/crystal_plasticity/include -I$with_cnstvm/poro_viscoplasticity/include -I$with_cnstvm/damage/include -I$with_cnstvm/data_structure/include"
 LIBGCM_LIBADD="$with_cnstvm/lib/libConstitutiveModel.la"
  
 AC_SUBST([LIBGCM_CPPFLAGS], ["$LIBGCM_CPPFLAGS"])
 AC_SUBST([LIBGCM_LIBADD], ["$LIBGCM_LIBADD"]) 
 $1
])
