AC_ARG_VAR([MSNET_CARGS], [Additional arguments passed to Multiscale Network contrib])

AC_DEFUN([_HAVE_MSNET], [
  AC_DEFINE([HAVE_MSNET], [1], [Msnet support available])
  AC_DEFINE([HAVE_AS_GLOBAL], [1], [We have global memory])
  AC_DEFINE([HAVE_AS_REGISTERED], [1], [We have registered memory])
  have_msnet=yes
])

AC_DEFUN([_CONFIG_MSNET], [
 contrib=$1

 # add some default CARGS
 MSNET_CARGS="$MSNET_CARGS"
 ACX_CONFIGURE_DIR([$contrib], [$contrib], ["$MSNET_CARGS"])
 _HAVE_MSNET
 MSNET_INCLUDE="$MSNET_CPPFLAGS -I\$(top_srcdir)/$1/include"
 MSNET_LIBADD="$MSNET_LIBADD \$(top_builddir)/$1/src/libmsnet.la"
])

AC_DEFUN([_PKG_MSNET], [
 pkg=$1
 
 # search for a msnet pkg-config package
 PKG_CHECK_MODULES([MSNET], [$pkg],
   [_HAVE_MSNET
    CFLAGS="$CFLAGS $MSNET_CFLAGS"
    CPPFLAGS="$CPPFLAGS $MSNET_CFLAGS"
    CXXFLAGS="$CXXFLAGS $MSNET_CFLAGS"
    LDFLAGS="$LDFLAGS $MSNET_LDFLAGS"
    LIBS="$LIBS $MSNET_LIBS"]
    )
])

AC_DEFUN([_LIB_MSNET], [
 # look for msnet in the path
 build_msnet=no
 AC_CHECK_HEADER([msnet/Network.hpp],
  [_HAVE_MSNET
   LIBS="$LIBS -lmsnet"])])
])

AC_DEFUN([_WITH_MSNET], [
 pkg=$1
 
 # handle the with_msnet option, if enable_msnet is selected
 AS_CASE($with_msnet,
   [no], [AC_MSG_ERROR([--enable-msnet=$enable_msnet excludes --without-msnet])],

   # contrib means we should just go ahead and build the library
   [contrib], [build_msnet=yes],
   [yes], [build_msnet=yes],

   # system means that we look for a library in the system path, or a
   # default-named pkg-config package
   [system], [_LIB_MSNET
              AS_IF([test "x$have_msnet" != xyes], [_PKG_MSNET($pkg)])],

   # any other string is interpreted as a custom pkg-config package
   [_PKG_MSNET($with_msnet)])
])

AC_DEFUN([CONFIG_MSNET], [
 contrib=$1
 pkg=$2
 build_msnet=no

 # Allow the user to override the way we try and find msnet.
 AC_ARG_ENABLE([msnet],
   [AS_HELP_STRING([--enable-msnet],
                   [Enable the msnet network @<:@default=yes@:>@])],
   [], [enable_msnet=yes])

 AC_ARG_WITH(msnet,
   [AS_HELP_STRING([--with-msnet{=contrib,system,PKG}],
                   [How we find msnet @<:@default=contrib@:>@])],
   [], [with_msnet=contrib])

 AS_IF([test "x$enable_msnet" != xno], [_WITH_MSNET($pkg)])
 AS_IF([test "x$build_msnet" != xno], [_CONFIG_MSNET($contrib)])

 AC_SUBST(MSNET_INCLUDE)
 AC_SUBST(MSNET_LIBADD)
])
