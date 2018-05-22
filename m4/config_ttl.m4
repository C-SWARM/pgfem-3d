AC_DEFUN([CONFIG_TTL], [
 AC_ARG_WITH(ttl,
   [AS_HELP_STRING([--with-ttl=<path>], [Specify path to TTL installation.])])

 AS_IF([test "x$with_ttl" != x],
   [ttl_include="-I$with_ttl/include"])

 old_CPPFLAGS=$CPPFLAGS
 CPPFLAGS="$CPPFLAGS $ttl_include"
 CFLAGS="$CFLAGS $ttl_include"
 AC_CHECK_HEADER(ttl/ttl.h, [], [AC_MSG_ERROR(Failed to find ttl.h)])
 CPPFLAGS=$old_CPPFLAGS

 AC_SUBST([TTL_INCLUDE], ["$ttl_include"])
 $1
])
