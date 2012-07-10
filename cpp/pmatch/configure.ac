AC_PREREQ([2.68])
AC_INIT([libsrvf-pmatch], [0.1], [drobinso@math.fsu.edu])
AC_CONFIG_SRCDIR([src/paretoset.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_AUX_DIR([.])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programs.
AC_PROG_CXX
LT_INIT

# Checks for header files.
AC_CHECK_HEADERS([unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])


AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile])
AC_OUTPUT
