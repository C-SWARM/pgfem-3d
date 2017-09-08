AC_DEFUN([CONFIG_TESTS], [
 AC_ARG_ENABLE(tests,
   [AS_HELP_STRING([--enable-tests], [Enable tests @<:@default=yes@:>@])],
   [], [enable_tests=yes])

 AC_ARG_WITH(tests-nprocs,
   [AS_HELP_STRING([--with-tests-nprocs=<n>],
                   [The maximum number of processes to use during testing @<:@default=1@:>@])],
   [], [with_tests_nprocs=1])

 # check for the ability to build tests
 AC_CHECK_PROG([found_numdiff], [numdiff], [yes])
 AC_CHECK_PROG([found_matlab], [matlab], [yes])

 AS_IF([test "x$enable_tests" == xyes -a "x$found_numdiff" != xyes],
   [AC_MSG_WARN([Missing numdiff, disabling tests])])
 AS_IF([test "x$enable_tests" == xyes -a "x$found_matlab" != xyes],
   [AC_MSG_WARN([Missing matlab, disabling some tests])])
 AS_IF([test "x$enable_tests" == xyes -a "x$have_vtk" != xyes],
   [AC_MSG_WARN([Missing VTK, disabling some tests])])

 AC_SUBST([TESTS_NPROCS], [$with_tests_nprocs])

 AM_CONDITIONAL([BUILD_TESTS], [test "x$enable_tests" == xyes])
 AM_CONDITIONAL([HAVE_VTK], [test "x$have_vtk" == xyes])
 AM_CONDITIONAL([HAVE_NUMDIFF], [test "x$found_numdiff" == xyes])
 AM_CONDITIONAL([HAVE_MATLAB], [test "x$have_matlab" == xyes])
])
