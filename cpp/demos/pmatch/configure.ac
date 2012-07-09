AC_PREREQ([2.68])
AC_INIT([pmatch_demo], [0.1], [drobinso@math.fsu.edu])
AC_CONFIG_SRCDIR([src/ui.h])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Check for GSL
AC_CHECK_LIB([gsl],[gsl_matrix_alloc],[],[AC_MSG_ERROR([gsl not found])])

# Checks for header files.
AC_CHECK_HEADERS([unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
fltk_config_prog=fltk-config
AC_ARG_WITH([fltk-config],
  [AS_HELP_STRING([--with-fltk-config=PATH],
    [Path to the fltk-config program (optional)])],
  [fltk_config_prog="$withval"],
  [])
AC_CHECK_PROG([fltk_config_prog_found],"$fltk_config_prog",[yes],[no])
if test "$fltk_config_prog_found" = "yes"; then
  FLTK_CFLAGS="`$fltk_config_prog --use-gl --cflags`"
  FLTK_LDFLAGS="`$fltk_config_prog --use-gl --ldflags`"
  AC_SUBST(FLTK_CFLAGS)
  AC_SUBST(FLTK_LDFLAGS)
else
  AC_MSG_ERROR([cannot find $fltk_config_prog])
fi

# The SRVF sources and libraries in this directory tree
SRVF_CPPFLAGS="-I../../../srvf/src"
SRVF_CPPFLAGS+=" -I../../../plot/src"
SRVF_CPPFLAGS+=" -I../../../pmatch/src"
SRVF_LDFLAGS="-L../../../srvf/src/.libs"
SRVF_LDFLAGS+=" -L../../../plot/src/.libs"
SRVF_LDFLAGS+=" -L../../../pmatch/src/.libs"
SRVF_LIBS=" -lsrvf-pmatch -lsrvf-plot -lsrvf"
AC_SUBST(SRVF_CPPFLAGS)
AC_SUBST(SRVF_LDFLAGS)
AC_SUBST(SRVF_LIBS)

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
