AC_INIT([libsrvf], [0.1], [drobinso@math.fsu.edu])
AC_CONFIG_SRCDIR([src/matrix.cc])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([src/matrix.h])
AC_CONFIG_AUX_DIR([.])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_PROG_CXX
LT_INIT

# Check for GSL
AC_CHECK_LIB([gsl],[gsl_matrix_alloc])

# Check for OpenMP
AC_LANG(C++)
AC_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
