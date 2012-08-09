AC_INIT([libsrvf], [0.1], [drobinso@math.fsu.edu])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([srvf/matrix.h])
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])


# Checks for common programs
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_INSTALL
LT_INIT


# Checks for standard types
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])


# Check for GSL
AC_CHECK_LIB([gsl],[gsl_matrix_alloc],[LIBS="-lgsl -lgslcblas $LIBS"],[AC_MSG_ERROR([gsl not found])],[-lgslcblas])


# Check for OpenMP
AC_LANG(C++)
AC_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"


# Path to fltk-config
fltk_config_prog=fltk-config
AC_ARG_WITH([fltk-config],
  [AS_HELP_STRING([--with-fltk-config=PATH],
    [Path to the fltk-config program (optional)])],
  [fltk_config_prog="$withval"],
  [])


# Optional plotting support
plot_enabled="no"
AC_ARG_ENABLE(plot,
  [AS_HELP_STRING([--enable-plot],
    [Enable plotting support (requires FLTK 1 and OpenGL)])],
  [if test "$enableval" = yes; then 
    plot_enabled="yes"
  fi],
  [])
AM_CONDITIONAL([COND_PLOT],[test "$plot_enabled" = "yes"])

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

  # Check for OpenGL
  AC_CHECK_LIB([GL],[glBegin])
  AC_CHECK_LIB([GLU],[gluLookAt])
fi


AC_CONFIG_FILES([Makefile
                 srvf/Makefile
                 src/Makefile
                 tests/Makefile
                 tests/data/Makefile
                 plot/Makefile
                 demos/Makefile
                 demos/pmatch/Makefile
                 demos/pmatch/src/Makefile
                 demos/pmatch-gui/Makefile
                 demos/pmatch-gui/src/Makefile
                 demos/shapedist/Makefile
                 demos/shapedist/src/Makefile])
AC_OUTPUT
