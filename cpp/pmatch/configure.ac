AC_PREREQ([2.68])
AC_INIT([libsrvf-pmatch], [0.1], [drobinso@math.fsu.edu])
AC_CONFIG_SRCDIR([src/paretoset.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.
AC_PROG_CXX
LT_INIT

# Checks for header files.
AC_CHECK_HEADERS([unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])


# If unit tests enabled, enable recursion into tests subdirectory, 
# and also check for Boost.
enable_unit_tests=no
AC_ARG_ENABLE([unit-tests],
  [AS_HELP_STRING([--enable-unit-tests],
    [Enable unit tests])],
  [enable_unit_tests="yes"], [])

AM_CONDITIONAL([COND_UNIT_TESTS],
  [test "$enable_unit_tests" = "yes"])

if test "$enable_pmatch_tests" = "yes"; then
  AC_LANG(C++)
  saved_LDFLAGS="$LDFLAGS"
  LDFLAGS="$LDFLAGS -lboost_unit_test_framework"
  AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([#include <boost/test/unit_test.hpp>],[])],dnl
    [LIBS="-lboost_unit_test_framework $LIBS"]dnl
    AC_DEFINE([HAVE_LIBBOOST_UNIT_TEST_FRAMEWORK],[],[Boost unit test]),dnl
    [AC_MSG_WARN([Boost unit testing libraries not found])])dnl
  LDFLAGS="$saved_LDFLAGS"
fi

AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
