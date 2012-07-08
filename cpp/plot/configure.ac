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

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Check for FLTK
#AC_LANG(C++)
#saved_LDFLAGS="$LDFLAGS"
#LDFLAGS="$LDFLAGS -lfltk2"
#AC_LINK_IFELSE(
#  [AC_LANG_PROGRAM([#include <fltk/Widget.h>],[fltk::Widget w(0,0,0,0,"")])],dnl
#  [LIBS="-lfltk2 -lfltk2_gl $LIBS"]dnl
#  AC_DEFINE([HAVE_FLTK],[],[FLTK]),dnl
#  [AC_MSG_WARN([FLTK not found])])dnl
#LDFLAGS="$saved_LDFLAGS"

fltk_config_prog=fltk-config
AC_ARG_WITH([fltk-config],
  [AS_HELP_STRING([--with-fltk-config=PATH],
    [Path to the fltk-config program (optional)])],
  [fltk_config_prog="$withval"],
  [])

plot_enabled="no"
AC_ARG_ENABLE([plot],
  [AS_HELP_STRING([--enable-plot],
    [Enable plotting (requires FLTK and OpenGL)])],
  [if test "$enableval" = "yes"; then plot_enabled="yes"; fi],
  [])

if test "$plot_enabled"="yes"; then
  AC_CHECK_PROG([fltk_config_prog_found],"$fltk_config_prog",[yes],[no])
  if test "$fltk_config_prog_found" = "yes"; then
    FLTK_CFLAGS="`$fltk_config_prog --use-gl --cflags`"
    FLTK_LDFLAGS="`$fltk_config_prog --use-gl --ldflags`"
    AC_SUBST(FLTK_CFLAGS)
    AC_SUBST(FLTK_LDFLAGS)
  else
    AC_MSG_ERROR([cannot find $fltk_config_prog])
  fi
else
  AC_MSG_ERROR([pass --enable-plot to configure to enable plotting])
fi

# Check for OpenGL
AC_CHECK_LIB([GL],[glBegin])
AC_CHECK_LIB([GLU],[gluLookAt])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
