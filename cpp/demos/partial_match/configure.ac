AC_PREREQ([2.68])
AC_INIT([srvf_partialmatch], [0.1], [drobinso@math.fsu.edu])
AC_CONFIG_SRCDIR([src/paretoset.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([.])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.
AC_PROG_CXX

# Checks for libraries.

# Check for FLTK
AC_LANG(C++)
saved_LDFLAGS="$LDFLAGS"
LDFLAGS="$LDFLAGS -lfltk2"
AC_LINK_IFELSE(
  [AC_LANG_PROGRAM([#include <fltk/Widget.h>],[fltk::Widget w(0,0,0,0,"")])],dnl
  [LIBS="-lfltk2 -lfltk2_gl $LIBS"]dnl
  AC_DEFINE([HAVE_FLTK],[],[FLTK]),dnl
  [AC_MSG_WARN([FLTK not found])])dnl
LDFLAGS="$saved_LDFLAGS"

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


# Checks for header files.
AC_CHECK_HEADERS([unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 tests/Makefile])
AC_OUTPUT
