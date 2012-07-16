AC_PREREQ([2.68])
AC_INIT([libsrvf], [0.1], [drobinso@math.fsu.edu])
AC_CONFIG_SRCDIR([src/matrix.cc])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])


# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
LT_INIT


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])


# Checks for library functions.
AC_CHECK_FUNCS([sqrt])


# Check for GSL
AC_CHECK_LIB([gsl],[gsl_matrix_alloc],[],
  [AC_MSG_ERROR([gsl not found])],[-lgslcblas])


# Check for OpenMP
AC_LANG(C++)
AC_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"


# Optional FLTK support.  Only used by the unit tests.
fltk_config_prog="fltk-config"
AC_ARG_WITH([fltk-config],
  [AS_HELP_STRING([--with-fltk-config=PATH],
    [Path to the fltk-config program])],
  [fltk_config_prog="$withval"],
  [])

plot_enabled="no"
AC_ARG_ENABLE([plot],
  [AS_HELP_STRING([--enable-plot],
    [Enable plotting (requires FLTK and OpenGL)])],
  [if test "$enableval" = "yes"; then plot_enabled="yes"; fi],
  [])

if test "$plot_enabled" = yes; then :
  AC_CHECK_PROG([fltk_config_prog_found],"$fltk_config_prog",[yes],[no])
  if test "$fltk_config_prog_found" = "yes"; then
    FLTK_CFLAGS="`$fltk_config_prog --use-gl --cflags`"
    FLTK_LDFLAGS="`$fltk_config_prog --use-gl --ldflags`"
    AC_SUBST(FLTK_CFLAGS)
    AC_SUBST(FLTK_LDFLAGS)
  else
    AC_MSG_ERROR([cannot find $fltk_config_prog])
  fi
fi


AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
