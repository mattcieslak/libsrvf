AC_PREREQ([2.68])
AC_INIT([libsrvf-pmatch], [0.1], [drobinso@math.fsu.edu])
AC_CONFIG_SRCDIR([src/paretoset.h])
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


# Checks for header files.
AC_CHECK_HEADERS([unistd.h])


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_CHECK_TYPES([ptrdiff_t])


# Checks for library functions.
AC_CHECK_FUNCS([sqrt])


# OpenMP
AC_LANG(C++)
AC_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"


AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
