# ======================================================================================
# Author: Francesco Montorsi
# RCS-ID: $Id: vtk.m4,v 1.1 2005/11/20 14:47:40 frm Exp $
#
# Implements the AM_OPTIONS_VTK, to add the --with-vtk=path option, and the
# AM_PATH_VTK macro used to detect VTK presence, location and version.
# ======================================================================================
#
# 18.8.2017 - Substantially modified by Luke D'Alessandro for PGFem3D framework.

#
#  AM_PATH_VTK([minimum-version], [action-if-found], [action-if-not-found])
#  ------------------------------------------------------------------------
#
#  NOTE: [minimum-version] must be in the form [X.Y.Z]
#
AC_DEFUN([AM_PATH_VTK], [
 # Use the VTK paths to build some variables that we will need.
 vtk_cppflags="-I$with_vtk/include/vtk$with_vtk_version"
 vtk_ldflags="-L$with_vtk/lib -Wl,-rpath,$with_vtk/lib"
 vtk_libs="-lvtkIOXML$with_vtk_version -lvtkIOXMLParser$with_vtk_version -lvtkIOCore$with_vtk_version -lvtkCommonExecutionModel$with_vtk_version -lvtkCommonDataModel$with_vtk_version -lvtkCommonMisc$with_vtk_version -lvtkCommonSystem$with_vtk_version -lvtkCommonTransforms$with_vtk_version -lvtkCommonMath$with_vtk_version -lvtkIOGeometry$with_vtk_version -lvtkCommonCore$with_vtk_version -lvtksys$with_vtk_version"

 # Make sure that the vtk headers are found.
 AC_CHECK_FILE([$with_vtk/include/vtk$with_vtk_version/vtkVersionMacros.h], [found_vtk_headers=yes])

 AC_MSG_CHECKING([if VTK is installed in $with_vtk])

 AS_IF([test "x$found_vtk_headers" != xyes],
   [AC_MSG_RESULT([no])
    $3],
   [AC_MSG_RESULT([yes])
    # now, eventually check version
    AS_IF([test -n "$1"],
      [# A version was specified... parse the version string in $1
       # The version of VTK that we need:
       maj=`echo $1 | sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
       min=`echo $1 | sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
       rel=`echo $1 | sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`

       AC_MSG_CHECKING([if VTK version is at least $maj.$min.$rel])
     
       # Compare required version of VTK against installed version:
       #
       # Note that in order to be able to compile the following test program,
       # we need to add to the current flags, the VTK settings...
       old_CPPFLAGS=$CPPFLAGS
       old_LDFLAGS=$LDFLAGS
       CPPFLAGS="$CPPFLAGS $vtk_cppflags"
       LDFLAGS="$LDFLAGS $vtk_ldflags $vtk_libs"
       #
       # check if the installed VTK is greater or not
       AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
           [#include <vtkVersionMacros.h>
            #include <stdio.h>
           ],
           [printf("VTK version is: %d.%d.%d", VTK_MAJOR_VERSION, VTK_MINOR_VERSION, VTK_BUILD_VERSION);
            #if (VTK_MAJOR_VERSION < $maj) || \
                (VTK_MAJOR_VERSION == $maj && VTK_MINOR_VERSION < $min) || \
                (VTK_MAJOR_VERSION == $maj && VTK_MINOR_VERSION == $min && VTK_BUILD_VERSION < $rel)
            #error Installed VTK is too old !
            #endif
           ])],
         [vtkVersion="OK"])
       #
       # restore all flags without VTK values
       CPPFLAGS=$old_CPPFLAGS
       LDFLAGS=$old_LDFLAGS

       # Execute $2 if version is ok, otherwise execute $3
       AS_IF([test "$vtkVersion" = "OK"],
         [AC_MSG_RESULT([yes])
          $2],
         [AC_MSG_RESULT([no])
          $3])],
      [# A target version number was not provided... execute $2 unconditionally
       #
       # if we don't have to check for minimum version (because the user did not
       # set that option), then we can execute here the block action-if-found
       $2])])

  # a hack to figure out if the local system has
  # -lvtkzlib$with_vtk_version
  # -lvtkexpat$with_vtk_version
  old_LDFLAGS=$LDFLAGS
  LDFLAGS="$LDFLAGS $vtk_ldflags"
  AC_CHECK_LIB([vtkzlib$with_vtk_version], [vtk_zlib_zlibVersion], [vtk_libs="$vtk_libs -lvtkzlib$with_vtk_version"])
  AC_CHECK_LIB([vtkexpat$with_vtk_version], [vtk_expat_XML_ExpatVersion], [vtk_libs="$vtk_libs -lvtkexpat$with_vtk_version"])
  LDFLAGS=$old_LDFLAGS
])# AM_PATH_VTK

AC_DEFUN([CONFIG_VTK], [
 AC_ARG_ENABLE(vtk,
   [AS_HELP_STRING([--enable-vtk],
 		          [Output in VTK binary output @<:@default=yes@:>@])],
   [], [enable_vtk=yes])

 AC_ARG_WITH([vtk],
   [AC_HELP_STRING([--with-vtk],
                   [The prefix where VTK is installed @<:@default=/usr@:>@])],
   [], [with_vtk="/usr"])

 AC_ARG_WITH([vtk-version],
   [AC_HELP_STRING([--with-vtk-version],
                  [VTK's include directory name is vtk-suffix, e.g. vtk-5.10/. What's the suffix? @<:@default=-5.10.1@:>@])],
   [], [with_vtk_version="-5.10.1"])

 AS_IF([test "x$enable_vtk" == xyes],
   [AM_PATH_VTK([$1],
                [have_vtk=yes],
                [AC_ERROR(Could not find VTK $with_vtk_version)])])
                
 AS_IF([test "x$have_vtk" != xyes],
   [AC_DEFINE(NO_VTK_LIB, 1, [Do not use VTK libraries])])

 AC_SUBST([VTK_CPPFLAGS], ["$vtk_cppflags"])
 AC_SUBST([VTK_LDFLAGS], ["$vtk_ldflags"])
 AC_SUBST([VTK_LIBS], ["$vtk_libs"])
])
