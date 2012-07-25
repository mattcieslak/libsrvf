AC_PREREQ([2.68])
AC_INIT([libsrvf], [0.1], [drobinso@math.fsu.edu])
AC_CONFIG_SRCDIR([src/main.cc])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])


# Checks for programs.
AC_PROG_CXX
AC_PROG_CC


# Check for OpenMP
AC_LANG(C++)
AC_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"


# Checks for header files.
AC_CHECK_HEADERS([unistd.h])


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T


# Check for GSL
AC_CHECK_LIB([gsl],[gsl_matrix_alloc],[LIBS="-lgsl -lgslcblas $LIBS"],[AC_MSG_ERROR([gsl not found])],[-lgslcblas])


# The SRVF sources and libraries in this directory tree
SRVF_CPPFLAGS="-I../../../srvf/src"
SRVF_CPPFLAGS+=" -I../../../pmatch/src"
SRVF_LDFLAGS="-L../../../srvf/src/.libs"
SRVF_LDFLAGS+=" -L../../../pmatch/src/.libs"
SRVF_LIBS=" -lsrvf-pmatch -lsrvf"
AC_SUBST(SRVF_CPPFLAGS)
AC_SUBST(SRVF_LDFLAGS)
AC_SUBST(SRVF_LIBS)

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
