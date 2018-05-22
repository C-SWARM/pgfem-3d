AC_DEFUN([_HAVE_PHOTON], [
  AC_DEFINE([HAVE_PHOTON], [1], [Photon support available])
  AC_DEFINE([HAVE_AS_GLOBAL], [1], [We have global memory])
  AC_DEFINE([HAVE_AS_REGISTERED], [1], [We have registered memory])
  have_photon=yes
])

AC_DEFUN([_PKG_PHOTON], [
 pkg=$1
 
 # search for a photon pkg-config package
 PKG_CHECK_MODULES([PHOTON], [$pkg],
   [_HAVE_PHOTON
    CFLAGS="$CFLAGS $PHOTON_CFLAGS"
    CPPFLAGS="$CPPFLAGS $PHOTON_CFLAGS"
    CXXFLAGS="$CXXFLAGS $PHOTON_CFLAGS"
    LDFLAGS="$LDFLAGS $PHOTON_LDFLAGS"
    LIBS="$LIBS $PHOTON_LIBS"])
])

AC_DEFUN([_LIB_PHOTON], [
 # look for photon in the path
 AC_CHECK_HEADER([photon.h],
   [AC_CHECK_LIB([photon], [photon_init],
     [_HAVE_PHOTON
      LIBS="$LIBS -lphoton"])])
])

AC_DEFUN([_WITH_PHOTON], [
 pkg=$1
 
 # handle the with_photon option, if enable_photon is selected
 AS_CASE($with_photon,
   [no], [AC_MSG_ERROR([--enable-photon=$enable_photon excludes --without-photon])],
   
   # system means that we look for a library in the system path, or a
   # default-named pkg-config package
   [system], [_LIB_PHOTON
              AS_IF([test "x$with_photon" != xyes], [_PKG_PHOTON($pkg)])],

   # any other string is interpreted as a custom pkg-config package
   [_PKG_PHOTON($with_photon)])
])

AC_DEFUN([CONFIG_PHOTON], [
 pkg=$1
 
 # Allow the user to override the way we try and find photon.
 AC_ARG_ENABLE([photon],
   [AS_HELP_STRING([--enable-photon],
                   [Enable the photon network @<:@default=no@:>@])],
   [], [enable_photon=no])

 AC_ARG_WITH(photon,
   [AS_HELP_STRING([--with-photon{=system,PKG}],
                   [How we find photon @<:@default=system@:>@])],
   [], [with_photon=system])

 AS_IF([test "x$enable_photon" != xno], [_WITH_PHOTON($pkg)])
])
