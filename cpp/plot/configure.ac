AC_PREREQ([2.68])
AC_INIT([libsrvf-plot], [0.1], [drobinso@math.fsu.edu])
AC_CONFIG_SRCDIR([src/render.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.
AC_PROG_CXX
AC_PROG_LIBTOOL
LT_INIT

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

# Check for OpenGL
AC_CHECK_LIB([GL],[glBegin])
AC_CHECK_LIB([GLU],[gluLookAt])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
