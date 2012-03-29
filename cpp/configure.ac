AC_INIT([libsrvf], [0.0.0+], [drobinso@math.fsu.edu])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_CONFIG_SRCDIR([src/matrix.h])
AC_CONFIG_HEADERS([config.h])

AC_PROG_CXX
AC_PROG_RANLIB

# Check for GSL
AC_CHECK_LIB([gsl],[gsl_matrix_alloc])

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


AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
