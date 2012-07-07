AC_PREREQ([2.68])
AC_INIT([libsrvf-tests], [0.1])
AC_CONFIG_SRCDIR([src/plotwin_tests.cc])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([.])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.
AC_PROG_CXX

# Check for Boost
AC_LANG(C++)
saved_LDFLAGS="$LDFLAGS"
LDFLAGS="$LDFLAGS -lboost_unit_test_framework"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <boost/test/unit_test.hpp>],[])],dnl
  [LIBS="-lboost_unit_test_framework $LIBS"]dnl
  AC_DEFINE([HAVE_LIBBOOST_UNIT_TEST_FRAMEWORK],[],[Boost unit test]),dnl
  [AC_MSG_WARN([Boost unit testing libraries not found])])dnl
LDFLAGS="$saved_LDFLAGS"

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
